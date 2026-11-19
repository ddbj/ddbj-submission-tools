import pytest

from ddbj_auth import cli, oauth
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
