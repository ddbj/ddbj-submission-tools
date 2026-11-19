"""ddbj-auth: log in to the DDBJ account from the command line (device flow), one login per API.

  ddbj-auth login  --for kura [--env staging|production] [--offline]
  ddbj-auth token  --for kura          # prints an access token for kura (stdout only)
  ddbj-auth status [--env ...]
  ddbj-auth logout [--for kura] [--env ...]

Use the token inside a command, never print it into a conversation or a log:
  curl -H "Authorization: Bearer $(ddbj-auth token --for kura)" ...
"""
from __future__ import annotations

import argparse
import datetime
import os
import sys
import webbrowser

from . import oauth, store
from .oauth import AuthError, has_audience

ISSUERS = {
    "staging": "https://idp-staging.ddbj.nig.ac.jp/realms/master",
    "production": "https://idp.ddbj.nig.ac.jp/realms/master",
}
APIS = {"kura": {"client_id": "kura", "audience": "kura"}}


class UsageError(Exception):
    pass


def _api(args) -> dict:
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
    if complete:
        try:
            webbrowser.open(complete)
        except Exception:  # noqa: BLE001 - no browser on a headless server: the URL is enough
            pass


def _login(args) -> int:
    api = _api(args)
    issuer = ISSUERS[args.env]
    scope = "openid offline_access" if args.offline else "openid"
    tokens = oauth.device_login(issuer, api["client_id"], scope, show=_show)
    access = store.claims(tokens.get("access_token", ""))
    if not has_audience(access, api["audience"]):
        raise AuthError(f"the access token's aud does not include {api['audience']!r}; check the client's audience mapper")
    store.save(args.env, args.api, {
        "refresh_token": tokens["refresh_token"], "client_id": api["client_id"], "issuer": issuer,
        "audience": api["audience"], "username": access.get("preferred_username"),
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


def _token(args) -> int:
    issuer = ISSUERS[args.env]
    env_token = os.environ.get("DDBJ_TOKEN")
    if env_token:
        _check_env_token(env_token, issuer, _api(args)["audience"])
        print(env_token)
        return 0
    rec = store.load(args.env, args.api)
    if rec is None:
        raise AuthError(f"not logged in: run ddbj-auth login --for {args.api} --env {args.env}")
    try:
        tokens = oauth.refresh(rec["issuer"], rec["client_id"], rec["refresh_token"])
    except AuthError as exc:
        raise AuthError(f"{exc}: run ddbj-auth login --for {args.api} --env {args.env}") from None
    if tokens.get("refresh_token"):
        store.save(args.env, args.api, rec | {
            "refresh_token": tokens["refresh_token"],
            "refresh_expires_at": store.claims(tokens["refresh_token"]).get("exp"),
        })
    if not has_audience(store.claims(tokens.get("access_token", "")), rec["audience"]):
        raise AuthError(f"the access token's aud does not include {rec['audience']!r}")
    print(tokens["access_token"])
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
        print(f"{api}\t{rec.get('username')}\trefresh until {_when(rec.get('refresh_expires_at'))}")
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
    common(sub.add_parser("status", help="show logins (never the tokens)"), False)
    common(sub.add_parser("logout", help="revoke and remove the stored login"), False)
    args = ap.parse_args(argv)
    try:
        return {"login": _login, "token": _token, "status": _status, "logout": _logout}[args.cmd](args)
    except UsageError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    except AuthError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
