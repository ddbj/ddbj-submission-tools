"""Keycloak endpoints: device flow with PKCE (RFC 8628 / 7636), refresh and revoke.

Every HTTP call goes through post_form(), which tests replace with a fake IdP.
Error messages never contain a token.
"""
from __future__ import annotations

import base64
import hashlib
import json
import secrets
import time
import urllib.error
import urllib.parse
import urllib.request
from typing import Callable

DEVICE_GRANT = "urn:ietf:params:oauth:grant-type:device_code"


class AuthError(Exception):
    """Login is needed, the IdP refused, or the IdP cannot be reached."""


def _endpoint(issuer: str, name: str) -> str:
    return f"{issuer.rstrip('/')}/protocol/openid-connect/{name}"


def post_form(url: str, data: dict[str, str]) -> dict:
    req = urllib.request.Request(url, data=urllib.parse.urlencode(data).encode(),
                                 headers={"Content-Type": "application/x-www-form-urlencoded",
                                          "Accept": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            body = resp.read()
    except urllib.error.HTTPError as exc:
        body = exc.read()
        if exc.code >= 500:
            raise AuthError(f"the IdP returned HTTP {exc.code}") from None
    except (urllib.error.URLError, OSError) as exc:
        reason = getattr(exc, "reason", exc)
        raise AuthError(f"cannot reach the IdP: {reason}") from None
    try:
        out = json.loads(body or b"{}")
    except ValueError:
        raise AuthError("the IdP returned a response that is not JSON") from None
    return out if isinstance(out, dict) else {}


def pkce_pair() -> tuple[str, str]:
    verifier = secrets.token_urlsafe(64)[:96]
    challenge = base64.urlsafe_b64encode(hashlib.sha256(verifier.encode()).digest()).decode().rstrip("=")
    return verifier, challenge


def _error(resp: dict) -> str:
    desc = resp.get("error_description")
    return f"{resp['error']}" + (f" ({desc})" if desc else "")


def device_login(issuer: str, client_id: str, scope: str,
                 show: Callable[[str, str, str | None], None],
                 sleep: Callable[[float], None] = time.sleep,
                 now: Callable[[], float] = time.time) -> dict:
    verifier, challenge = pkce_pair()
    dev = post_form(_endpoint(issuer, "auth/device"),
                    {"client_id": client_id, "scope": scope,
                     "code_challenge": challenge, "code_challenge_method": "S256"})
    if "error" in dev or "device_code" not in dev:
        raise AuthError(f"the IdP refused the login: {_error(dev) if 'error' in dev else 'no device code'}")
    show(dev["verification_uri"], dev["user_code"], dev.get("verification_uri_complete"))
    interval = float(dev.get("interval", 5))
    deadline = now() + float(dev.get("expires_in", 600))
    while True:
        sleep(interval)
        if now() > deadline:
            raise AuthError("login timed out; run login again")
        resp = post_form(_endpoint(issuer, "token"),
                         {"grant_type": DEVICE_GRANT, "device_code": dev["device_code"],
                          "client_id": client_id, "code_verifier": verifier})
        err = resp.get("error")
        if err is None:
            return resp
        if err == "authorization_pending":
            continue
        if err == "slow_down":
            interval += 5
            continue
        raise AuthError(f"login failed: {_error(resp)}")


def refresh(issuer: str, client_id: str, refresh_token: str) -> dict:
    resp = post_form(_endpoint(issuer, "token"),
                     {"grant_type": "refresh_token", "client_id": client_id, "refresh_token": refresh_token})
    if resp.get("error") == "invalid_grant":
        raise AuthError("login expired or was revoked")
    if "error" in resp:
        raise AuthError(f"the IdP refused the refresh: {_error(resp)}")
    return resp


def revoke(issuer: str, client_id: str, refresh_token: str) -> None:
    resp = post_form(_endpoint(issuer, "revoke"),
                     {"client_id": client_id, "token": refresh_token, "token_type_hint": "refresh_token"})
    if "error" in resp:
        raise AuthError(f"the IdP refused the revocation: {_error(resp)}")


def has_audience(token_claims: dict, audience: str) -> bool:
    aud = token_claims.get("aud")
    return aud == audience if isinstance(aud, str) else isinstance(aud, list) and audience in aud
