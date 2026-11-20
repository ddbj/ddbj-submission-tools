import json

import pytest

from ddbj_auth import cli, oauth, sts
from ddbj_auth.cli import ISSUERS, main
from ddbj_auth.oauth import AuthError
from ddbj_auth.store import list_apis, load, save
from tests.fake_idp import FakeIdp, make_jwt

STG = ISSUERS["staging"]
TOKENS = {"access_token": make_jwt(aud="kura", preferred_username="u1"), "refresh_token": make_jwt(exp=2000000000, n=1)}
ROTATED = {"access_token": make_jwt(aud=["kura", "account"], preferred_username="u1", n=2),
           "refresh_token": make_jwt(exp=2000000100, n=2)}


@pytest.fixture
def fake(monkeypatch, tmp_path):
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
    monkeypatch.delenv("DDBJ_TOKEN", raising=False)
    idp = FakeIdp()
    idp.device["interval"] = 0
    monkeypatch.setattr(oauth, "post_form", idp.post_form)
    monkeypatch.setattr(cli.webbrowser, "open", lambda *a, **k: False)
    return idp


def _login(fake):
    fake.token_queue.append(TOKENS)
    assert main(["login", "--for", "kura", "--env", "staging"]) == 0


def test_login_then_token(fake, capsys):
    _login(fake)
    assert "as u1" in capsys.readouterr().err
    rec = load("staging", "kura")
    assert rec["issuer"] == STG and rec["client_id"] == "kura" and rec["refresh_expires_at"] == 2000000000
    fake.token_queue.append(ROTATED)
    assert main(["token", "--for", "kura", "--env", "staging"]) == 0
    out = capsys.readouterr()
    assert out.out.strip() == ROTATED["access_token"] and ROTATED["access_token"] not in out.err
    assert load("staging", "kura")["refresh_token"] == ROTATED["refresh_token"]


def test_token_not_logged_in(fake, capsys):
    assert main(["token", "--for", "kura"]) == 1
    assert "ddbj-auth login --for kura" in capsys.readouterr().err


def test_token_after_refresh_expired_keeps_file(fake, capsys):
    _login(fake)
    before = load("staging", "kura")
    fake.token_queue.append({"error": "invalid_grant"})
    assert main(["token", "--for", "kura", "--env", "staging"]) == 1
    assert load("staging", "kura") == before
    assert "login" in capsys.readouterr().err


def test_token_wrong_audience(fake, capsys):
    _login(fake)
    fake.token_queue.append({"access_token": make_jwt(aud=["account"]), "refresh_token": make_jwt(n=3)})
    assert main(["token", "--for", "kura", "--env", "staging"]) == 1
    assert "aud" in capsys.readouterr().err


def test_ddbj_token_env_checked(fake, monkeypatch, capsys):
    good = make_jwt(iss=STG, aud="kura")
    monkeypatch.setenv("DDBJ_TOKEN", good)
    assert main(["token", "--for", "kura", "--env", "staging"]) == 0
    assert capsys.readouterr().out.strip() == good
    monkeypatch.setenv("DDBJ_TOKEN", make_jwt(iss=ISSUERS["production"], aud="kura"))
    assert main(["token", "--for", "kura", "--env", "staging"]) == 1
    assert "iss" in capsys.readouterr().err
    monkeypatch.setenv("DDBJ_TOKEN", make_jwt(iss=STG, aud="other"))
    assert main(["token", "--for", "kura", "--env", "staging"]) == 1
    assert "aud" in capsys.readouterr().err


def test_status_never_prints_tokens(fake, capsys):
    _login(fake)
    capsys.readouterr()
    assert main(["status", "--env", "staging"]) == 0
    out = capsys.readouterr().out
    assert "kura" in out and "u1" in out and TOKENS["refresh_token"] not in out


def test_status_not_logged_in(fake, capsys):
    assert main(["status", "--env", "staging"]) == 0
    assert "not logged in" in capsys.readouterr().out


def test_logout_revokes_and_deletes(fake):
    _login(fake)
    assert main(["logout", "--env", "staging"]) == 0
    assert fake.urls[-1].endswith("/revoke") and fake.requests[-1]["token"] == TOKENS["refresh_token"]
    assert list_apis("staging") == []


def test_logout_deletes_even_if_revoke_fails(fake, capsys):
    _login(fake)
    fake.fail_revoke = AuthError("cannot reach the IdP: timed out")
    assert main(["logout", "--for", "kura", "--env", "staging"]) == 0
    assert list_apis("staging") == [] and "could not revoke" in capsys.readouterr().err


def test_unknown_api_is_usage_error(fake, capsys):
    assert main(["login", "--for", "nope"]) == 2
    assert "kura" in capsys.readouterr().err


def test_unknown_api_with_client_and_audience(fake):
    fake.token_queue.append({"access_token": make_jwt(aud="https://x/mcp", preferred_username="u1"),
                             "refresh_token": make_jwt(n=4)})
    assert main(["login", "--for", "mcp", "--env", "staging", "--client-id", "mcp-c", "--audience", "https://x/mcp"]) == 0
    assert load("staging", "mcp")["client_id"] == "mcp-c"


def test_idp_down_is_exit_1_without_traceback(fake, monkeypatch, capsys):
    def down(url, data):
        raise AuthError("cannot reach the IdP: timed out")
    monkeypatch.setattr(oauth, "post_form", down)
    assert main(["login", "--for", "kura", "--env", "staging"]) == 1
    err = capsys.readouterr().err
    assert "cannot reach the IdP" in err and "Traceback" not in err


def test_offline_adds_scope(fake):
    fake.token_queue.append(TOKENS)
    assert main(["login", "--for", "kura", "--env", "staging", "--offline"]) == 0
    assert fake.requests[0]["scope"] == "openid offline_access"


def test_saved_record_kept_when_save_unrelated(fake):
    save("staging", "other", {"refresh_token": "x"})
    _login(fake)
    assert list_apis("staging") == ["kura", "other"]


@pytest.mark.parametrize(("display", "platform", "opened"), [(None, "linux", False), (":0", "linux", True), (None, "darwin", True)])
def test_browser_opened_only_with_a_gui(fake, monkeypatch, display, platform, opened):
    calls = []
    monkeypatch.setattr(cli.webbrowser, "open", lambda *a, **k: calls.append(a) or True)
    monkeypatch.delenv("WAYLAND_DISPLAY", raising=False)
    if display:
        monkeypatch.setenv("DISPLAY", display)
    else:
        monkeypatch.delenv("DISPLAY", raising=False)
    monkeypatch.setattr(cli.sys, "platform", platform)
    _login(fake)
    assert bool(calls) is opened


@pytest.mark.parametrize("api", ["../x", "a/b", ""])
def test_bad_api_name_is_usage_error(fake, api, capsys):
    assert main(["login", "--for", api, "--env", "staging", "--client-id", "c", "--audience", "a"]) == 2


def test_parallel_token_calls_share_rotation(fake, monkeypatch, capsys):
    import threading
    import time

    _login(fake)
    current = {"rt": TOKENS["refresh_token"], "n": 10}
    guard = threading.Lock()

    def rotating(url, data):
        time.sleep(0.05)                     # widen the race window
        with guard:
            if data.get("refresh_token") != current["rt"]:
                return {"error": "invalid_grant"}     # reuse of a rotated refresh token
            current["n"] += 1
            current["rt"] = make_jwt(n=current["n"])
            return {"access_token": make_jwt(aud="kura", n=current["n"]), "refresh_token": current["rt"]}
    monkeypatch.setattr(oauth, "post_form", rotating)
    codes = []
    threads = [threading.Thread(target=lambda: codes.append(main(["token", "--for", "kura", "--env", "staging"]))) for _ in range(4)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    assert codes == [0, 0, 0, 0]


def test_login_without_refresh_token_is_clear_error(fake, capsys):
    fake.token_queue.append({"access_token": make_jwt(aud="kura")})
    assert main(["login", "--for", "kura", "--env", "staging"]) == 1
    err = capsys.readouterr().err
    assert "refresh token" in err and "Traceback" not in err


def test_ctrl_c_exits_130(fake, monkeypatch, capsys):
    def interrupted(*a, **k):
        raise KeyboardInterrupt
    monkeypatch.setattr(oauth, "device_login", interrupted)
    assert main(["login", "--for", "kura", "--env", "staging"]) == 130


STS_OK = """<AssumeRoleWithWebIdentityResponse xmlns="https://sts.amazonaws.com/doc/2011-06-15/">
<AssumeRoleWithWebIdentityResult><Credentials>
<AccessKeyId>AKID</AccessKeyId><SecretAccessKey>SECRET</SecretAccessKey>
<SessionToken>STOKEN</SessionToken><Expiration>2026-11-20T12:00:00Z</Expiration>
</Credentials></AssumeRoleWithWebIdentityResult></AssumeRoleWithWebIdentityResponse>"""
STS_ERR = """<ErrorResponse><Error><Code>InvalidParameterValue</Code>
<Message>token audience matches none of the configured client IDs</Message></Error></ErrorResponse>"""


@pytest.fixture
def fake_sts(monkeypatch):
    class Calls(list):
        reply = (200, STS_OK)
    calls = Calls()

    def post(url, data):
        calls.append((url, dict(data)))
        return calls.reply
    monkeypatch.setattr(sts, "post", post)
    return calls


def test_aws_credentials_prints_credential_process_json(fake, fake_sts, capsys):
    _login(fake)
    capsys.readouterr()
    fake.token_queue.append(ROTATED)
    assert main(["aws-credentials", "--for", "kura", "--env", "staging"]) == 0
    out = capsys.readouterr()
    assert json.loads(out.out) == {"Version": 1, "AccessKeyId": "AKID", "SecretAccessKey": "SECRET",
                                   "SessionToken": "STOKEN", "Expiration": "2026-11-20T12:00:00Z"}
    assert ROTATED["access_token"] not in out.out + out.err
    url, form = fake_sts[0]
    assert url == "https://kura-s3-staging.ddbj.nig.ac.jp"
    assert form["Action"] == "AssumeRoleWithWebIdentity" and form["RoleArn"] == "arn:aws:iam::role/KuraUserRole"
    assert form["WebIdentityToken"] == ROTATED["access_token"]


def test_aws_credentials_sts_refusal_is_exit_1_without_token(fake, fake_sts, capsys):
    _login(fake)
    capsys.readouterr()
    fake.token_queue.append(ROTATED)
    fake_sts.reply = (400, STS_ERR)
    assert main(["aws-credentials", "--for", "kura", "--env", "staging"]) == 1
    out = capsys.readouterr()
    assert "InvalidParameterValue" in out.err and out.out == ""
    assert ROTATED["access_token"] not in out.err


def test_aws_credentials_uses_ddbj_token(fake, fake_sts, monkeypatch, capsys):
    tok = make_jwt(iss=ISSUERS["production"], aud="kura")
    monkeypatch.setenv("DDBJ_TOKEN", tok)
    assert main(["aws-credentials", "--for", "kura"]) == 0
    assert fake_sts[0][0] == "https://kura-s3.ddbj.nig.ac.jp" and fake_sts[0][1]["WebIdentityToken"] == tok


def test_aws_credentials_only_for_s3_apis(fake, fake_sts, capsys):
    assert main(["aws-credentials", "--for", "other"]) == 2
    assert "kura" in capsys.readouterr().err and fake_sts == []


def test_status_shows_bucket(fake, capsys):
    fake.token_queue.append({"access_token": make_jwt(aud="kura", preferred_username="u1", sub="0f3c-uuid"),
                             "refresh_token": make_jwt(exp=2000000000)})
    assert main(["login", "--for", "kura", "--env", "staging"]) == 0
    capsys.readouterr()
    main(["status", "--env", "staging"])
    assert "0f3c-uuid" in capsys.readouterr().out
