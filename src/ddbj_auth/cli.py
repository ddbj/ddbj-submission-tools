"""ddbj-auth: log in to the DDBJ account from the command line (device flow), one login per API.

  ddbj-auth login  --for kura [--env staging|production] [--offline]
  ddbj-auth token  --for kura          # prints an access token for kura (stdout only)
  ddbj-auth aws-credentials --for kura # temporary S3 credentials, for credential_process
  ddbj-auth status [--env ...]
  ddbj-auth logout [--for kura] [--env ...]

Use the token inside a command, never print it into a conversation or a log:
  curl -H "Authorization: Bearer $(ddbj-auth token --for kura)" ...

For aws cli, boto3 or rclone, put this in ~/.aws/config (bucket = the sub that status shows):
  [profile kura]
  region = us-east-1
  endpoint_url = https://kura-s3.ddbj.nig.ac.jp
  credential_process = ddbj-auth aws-credentials --for kura
"""
from __future__ import annotations

import argparse
import datetime
import json
import os
import re
import sys
import webbrowser

from . import oauth, store, sts
from .oauth import AuthError, has_audience

ISSUERS = {
    "staging": "https://idp-staging.ddbj.nig.ac.jp/realms/master",
    "production": "https://idp.ddbj.nig.ac.jp/realms/master",
}
APIS = {"kura": {"client_id": "kura", "audience": "kura",
                 # S3 API and STS share one endpoint (kura docs/sp-integration.md)
                 "sts": {"staging": "https://kura-s3-staging.ddbj.nig.ac.jp",
                         "production": "https://kura-s3.ddbj.nig.ac.jp"},
                 "role_arn": "arn:aws:iam::role/KuraUserRole"}}


class UsageError(Exception):
    pass


_API_NAME = re.compile(r"[A-Za-z0-9_-]+")


def _api(args) -> dict:
    if not _API_NAME.fullmatch(args.api or ""):
        raise UsageError(f"bad API name {args.api!r}: use letters, digits, '-' and '_'")
    known = APIS.get(args.api, {})
    client_id = getattr(args, "client_id", None) or known.get("client_id")
    audience = getattr(args, "audience", None) or known.get("audience")
    if not client_id or not audience:
        raise UsageError(f"unknown API {args.api!r}: choose from {', '.join(sorted(APIS))}, "
                         "or give --client-id and --audience")
    return {"client_id": client_id, "audience": audience}


def _show(uri: str, code: str, complete: str | None) -> None:
    print(f"Open {complete or uri} and approve the login (code {code}).", file=sys.stderr)
    print("Approve only a code you started yourself.", file=sys.stderr)
    # Only where a graphical browser is likely: on a headless server Python would start a
    # text browser (lynx, w3m) in the terminal, hiding the code and blocking the polling
    gui = sys.platform in ("darwin", "win32") or os.environ.get("DISPLAY") or os.environ.get("WAYLAND_DISPLAY")
    if complete and gui:
        try:
            webbrowser.open(complete)
        except Exception:  # noqa: BLE001 - no browser on a headless server: the URL is enough
            pass


def _login(args) -> int:
    api = _api(args)
    issuer = ISSUERS[args.env]
    scope = "openid offline_access" if args.offline else "openid"
    tokens = oauth.device_login(issuer, api["client_id"], scope, show=_show)
    if not tokens.get("refresh_token"):
        raise AuthError("the IdP returned no refresh token; the client must allow refresh tokens")
    access = store.claims(tokens.get("access_token", ""))
    if not has_audience(access, api["audience"]):
        raise AuthError(f"the access token's aud does not include {api['audience']!r}; check the client's audience mapper")
    store.save(args.env, args.api, {
        "refresh_token": tokens["refresh_token"], "client_id": api["client_id"], "issuer": issuer,
        "audience": api["audience"], "username": access.get("preferred_username"), "sub": access.get("sub"),
        "refresh_expires_at": store.claims(tokens["refresh_token"]).get("exp"),
    })
    print(f"logged in to {args.api} ({args.env}) as {access.get('preferred_username')}", file=sys.stderr)
    return 0


def _check_env_token(token: str, issuer: str, audience: str) -> None:
    c = store.claims(token)
    if c.get("iss") != issuer:
        raise AuthError(f"DDBJ_TOKEN: iss is {c.get('iss')!r}, not {issuer!r} (--env)")
    if not has_audience(c, audience):
        raise AuthError(f"DDBJ_TOKEN: aud does not include {audience!r}")


def _access_token(args) -> str:
    issuer = ISSUERS[args.env]
    env_token = os.environ.get("DDBJ_TOKEN")
    if env_token:
        _check_env_token(env_token, issuer, _api(args)["audience"])
        return env_token
    if not _API_NAME.fullmatch(args.api or ""):
        raise UsageError(f"bad API name {args.api!r}")
    # Hold the lock across load, refresh and save: with refresh-token rotation, a second caller
    # must use the token the first one received (the IdP may reject a reused one)
    with store.locked(args.env, args.api):
        rec = store.load(args.env, args.api)
        if rec is None:
            raise AuthError(f"not logged in: run ddbj-auth login --for {args.api} --env {args.env}")
        try:
            tokens = oauth.refresh(rec["issuer"], rec["client_id"], rec["refresh_token"])
        except AuthError as exc:
            raise AuthError(f"{exc}: run ddbj-auth login --for {args.api} --env {args.env}") from None
        if tokens.get("refresh_token"):
            store.write(args.env, args.api, rec | {
                "refresh_token": tokens["refresh_token"],
                "refresh_expires_at": store.claims(tokens["refresh_token"]).get("exp"),
            })
    if not has_audience(store.claims(tokens.get("access_token", "")), rec["audience"]):
        raise AuthError(f"the access token's aud does not include {rec['audience']!r}")
    return tokens["access_token"]


def _token(args) -> int:
    print(_access_token(args))
    return 0


def _aws_credentials(args) -> int:
    known = APIS.get(args.api, {})
    if "sts" not in known:
        s3_apis = ", ".join(sorted(k for k, v in APIS.items() if "sts" in v))
        raise UsageError(f"{args.api!r} has no S3 API: choose from {s3_apis}")
    creds = sts.assume_role(known["sts"][args.env], known["role_arn"], _access_token(args))
    # The format of credential_process; aws cli and boto3 call again before Expiration
    print(json.dumps({"Version": 1} | creds))
    return 0


def _when(exp) -> str:
    if not exp:
        return "the IdP session"
    return datetime.datetime.fromtimestamp(int(exp), datetime.timezone.utc).isoformat(timespec="minutes")


def _status(args) -> int:
    apis = store.list_apis(args.env)
    if not apis:
        print(f"not logged in ({args.env})")
    for api in apis:
        rec = store.load(args.env, api) or {}
        print(f"{api}\t{rec.get('username')}\tsub {rec.get('sub') or '-'}\trefresh until {_when(rec.get('refresh_expires_at'))}")
    return 0


def _logout(args) -> int:
    for api in [args.api] if args.api else store.list_apis(args.env):
        rec = store.load(args.env, api)
        if rec is None:
            print(f"{api}: not logged in", file=sys.stderr)
            continue
        try:
            oauth.revoke(rec["issuer"], rec["client_id"], rec["refresh_token"])
        except AuthError as exc:
            print(f"{api}: could not revoke at the IdP ({exc}); removed the local token", file=sys.stderr)
        store.delete(args.env, api)
        print(f"logged out of {api} ({args.env})", file=sys.stderr)
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="ddbj-auth", description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__.split("\n", 1)[1])
    sub = ap.add_subparsers(dest="cmd", required=True)

    def common(p, api_required: bool):
        p.add_argument("--env", choices=sorted(ISSUERS), default="production")
        p.add_argument("--for", dest="api", required=api_required, help=f"API ({', '.join(sorted(APIS))}, or any with --client-id/--audience)")

    p = sub.add_parser("login", help="log in with the device flow")
    common(p, True)
    p.add_argument("--offline", action="store_true", help="ask for an offline token (survives the IdP session)")
    p.add_argument("--client-id")
    p.add_argument("--audience")
    p = sub.add_parser("token", help="print an access token for the API")
    common(p, True)
    p.add_argument("--audience")
    common(sub.add_parser("aws-credentials", help="print temporary S3 credentials (credential_process format)"), True)
    common(sub.add_parser("status", help="show logins (never the tokens)"), False)
    common(sub.add_parser("logout", help="revoke and remove the stored login"), False)
    args = ap.parse_args(argv)
    try:
        return {"login": _login, "token": _token, "aws-credentials": _aws_credentials,
                "status": _status, "logout": _logout}[args.cmd](args)
    except UsageError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    except AuthError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print("interrupted", file=sys.stderr)
        return 130


if __name__ == "__main__":
    sys.exit(main())
