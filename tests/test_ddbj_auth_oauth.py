import base64
import hashlib
import urllib.error

import pytest
from tests.fake_idp import ISS, FakeIdp, make_jwt

from ddbj_auth import oauth
from ddbj_auth.oauth import AuthError, device_login, has_audience, pkce_pair, refresh, revoke

TOKENS = {"access_token": make_jwt(aud="kura", preferred_username="u1"), "refresh_token": make_jwt(exp=2000000000),
          "expires_in": 300}


@pytest.fixture
def fake(monkeypatch):
    idp = FakeIdp()
    monkeypatch.setattr(oauth, "post_form", idp.post_form)
    return idp


class Clock:
    def __init__(self):
        self.t = 1000.0
        self.slept: list[float] = []

    def now(self):
        return self.t

    def sleep(self, s):
        self.slept.append(s)
        self.t += s


def _login(fake, clock=None, scope="openid"):
    clock = clock or Clock()
    shown = []
    return device_login(ISS, "kura", scope, show=lambda *a: shown.append(a), sleep=clock.sleep, now=clock.now), shown


def test_pkce_pair_s256():
    v, c = pkce_pair()
    assert 43 <= len(v) <= 128 and "=" not in c
    assert c == base64.urlsafe_b64encode(hashlib.sha256(v.encode()).digest()).decode().rstrip("=")


def test_device_login_polls_until_tokens(fake):
    fake.token_queue = [{"error": "authorization_pending"}, {"error": "slow_down"}, TOKENS]
    clock = Clock()
    result, shown = _login(fake, clock)
    assert result == TOKENS
    assert clock.slept == [5, 5, 10]          # interval 5; after slow_down, 10
    assert shown == [(f"{ISS}/device", "ABCD-EFGH", f"{ISS}/device?user_code=ABCD-EFGH")]


def test_device_login_sends_pkce_and_scope(fake):
    fake.token_queue = [TOKENS]
    _login(fake, scope="openid offline_access")
    dev, tok = fake.requests[0], fake.requests[1]
    assert fake.urls[0] == f"{ISS}/protocol/openid-connect/auth/device"
    assert dev["code_challenge_method"] == "S256" and dev["scope"] == "openid offline_access" and dev["client_id"] == "kura"
    assert fake.urls[1] == f"{ISS}/protocol/openid-connect/token"
    assert tok["grant_type"] == "urn:ietf:params:oauth:grant-type:device_code" and tok["device_code"] == "dc"
    assert tok["code_verifier"]


@pytest.mark.parametrize("err", ["access_denied", "expired_token"])
def test_device_login_stops_on(fake, err):
    fake.token_queue = [{"error": err}]
    with pytest.raises(AuthError, match=err):
        _login(fake)


def test_device_login_times_out(fake):
    fake.device["expires_in"] = 10
    fake.token_queue = [{"error": "authorization_pending"}] * 10
    with pytest.raises(AuthError, match="timed out"):
        _login(fake)


def test_refresh_invalid_grant(fake):
    fake.token_queue = [{"error": "invalid_grant", "error_description": "Token is not active"}]
    with pytest.raises(AuthError, match="login expired"):
        refresh(ISS, "kura", "r")


def test_refresh_returns_tokens(fake):
    fake.token_queue = [TOKENS]
    assert refresh(ISS, "kura", "r") == TOKENS
    assert fake.requests[-1] == {"grant_type": "refresh_token", "client_id": "kura", "refresh_token": "r"}


def test_network_error_is_auth_error(monkeypatch):
    def boom(*a, **k):
        raise urllib.error.URLError("timed out")
    monkeypatch.setattr(oauth.urllib.request, "urlopen", boom)
    with pytest.raises(AuthError, match="cannot reach the IdP") as e:
        oauth.post_form(f"{ISS}/protocol/openid-connect/token", {"refresh_token": "eyJsecret"})
    assert "eyJ" not in str(e.value)


@pytest.mark.parametrize(("aud", "ok"), [("kura", True), (["kura", "account"], True), (["account"], False), (None, False)])
def test_has_audience(aud, ok):
    assert has_audience({"aud": aud} if aud is not None else {}, "kura") is ok


def test_revoke_posts_token_and_hint(fake):
    revoke(ISS, "kura", "r")
    assert fake.urls[-1] == f"{ISS}/protocol/openid-connect/revoke"
    assert fake.requests[-1] == {"client_id": "kura", "token": "r", "token_type_hint": "refresh_token"}
