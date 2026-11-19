"""Refresh tokens on disk, one file per environment and API.

${XDG_CONFIG_HOME:-~/.config}/ddbj/<env>/<api>.json, directory 0700, file 0600.
Writes take a lock and replace the file atomically, so two processes refreshing at
the same time (refresh-token rotation) never leave a half-written file.
Access tokens are never stored.
"""
from __future__ import annotations

import base64
import binascii
import contextlib
import fcntl
import json
import os
import tempfile
from pathlib import Path


def config_dir(env: str) -> Path:
    base = os.environ.get("XDG_CONFIG_HOME") or os.path.join(os.path.expanduser("~"), ".config")
    return Path(base) / "ddbj" / env


def _path(env: str, api: str) -> Path:
    return config_dir(env) / f"{api}.json"


@contextlib.contextmanager
def _locked(env: str, api: str):
    d = config_dir(env)
    d.mkdir(parents=True, exist_ok=True)
    os.chmod(d, 0o700)
    with open(d / f"{api}.json.lock", "a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        try:
            yield d
        finally:
            fcntl.flock(lock, fcntl.LOCK_UN)


def load(env: str, api: str) -> dict | None:
    try:
        return json.loads(_path(env, api).read_text(encoding="utf-8"))
    except FileNotFoundError:
        return None


def save(env: str, api: str, record: dict) -> None:
    with _locked(env, api) as d:
        fd, tmp = tempfile.mkstemp(dir=d, prefix=f".{api}.json.tmp")
        try:
            os.chmod(tmp, 0o600)
            with os.fdopen(fd, "w", encoding="utf-8") as fh:
                json.dump(record, fh)
            os.replace(tmp, _path(env, api))
        except BaseException:
            with contextlib.suppress(FileNotFoundError):
                os.unlink(tmp)
            raise


def delete(env: str, api: str) -> bool:
    with _locked(env, api):
        try:
            _path(env, api).unlink()
            return True
        except FileNotFoundError:
            return False


def list_apis(env: str) -> list[str]:
    d = config_dir(env)
    return sorted(p.stem for p in d.glob("*.json")) if d.is_dir() else []


def claims(jwt_token: str) -> dict:
    """Payload of a JWT, NOT verified (the API that receives the token verifies it)."""
    parts = jwt_token.split(".")
    if len(parts) != 3:
        return {}
    try:
        payload = base64.urlsafe_b64decode(parts[1] + "=" * (-len(parts[1]) % 4))
        out = json.loads(payload)
    except (binascii.Error, ValueError):
        return {}
    return out if isinstance(out, dict) else {}
