"""S3 STS: exchange an access token for temporary S3 credentials (AssumeRoleWithWebIdentity).

The call needs no AWS credentials, so it is a plain form POST. post() is the only HTTP call;
tests replace it. Error messages never contain the token.
"""
from __future__ import annotations

import http.client
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

from .oauth import AuthError

FIELDS = ("AccessKeyId", "SecretAccessKey", "SessionToken", "Expiration")


def post(url: str, data: dict[str, str]) -> tuple[int, str]:
    req = urllib.request.Request(url, data=urllib.parse.urlencode(data).encode(),
                                 headers={"Content-Type": "application/x-www-form-urlencoded"})
    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            return resp.status, resp.read().decode("utf-8", "replace")
    except urllib.error.HTTPError as exc:
        return exc.code, exc.read().decode("utf-8", "replace")
    except (urllib.error.URLError, OSError, http.client.HTTPException) as exc:
        raise AuthError(f"cannot reach the STS: {getattr(exc, 'reason', exc)}") from None


def _find(root: ET.Element, name: str) -> str | None:
    for el in root.iter():
        if el.tag.rsplit("}", 1)[-1] == name:
            return (el.text or "").strip()
    return None


def assume_role(endpoint: str, role_arn: str, token: str, duration: int = 3600) -> dict[str, str]:
    """Temporary credentials {AccessKeyId, SecretAccessKey, SessionToken, Expiration}."""
    status, body = post(endpoint, {
        "Action": "AssumeRoleWithWebIdentity", "Version": "2011-06-15", "RoleArn": role_arn,
        "RoleSessionName": "ddbj-auth", "WebIdentityToken": token, "DurationSeconds": str(duration),
    })
    try:
        root = ET.fromstring(body)
    except ET.ParseError:
        raise AuthError(f"the STS returned HTTP {status} and a response that is not XML") from None
    if _find(root, "Code") is not None:
        raise AuthError(f"the STS refused the token: {_find(root, 'Code')} ({_find(root, 'Message')})")
    creds = {k: _find(root, k) for k in FIELDS}
    if status != 200 or not all(creds.values()):
        raise AuthError(f"the STS returned HTTP {status} without credentials")
    return creds
