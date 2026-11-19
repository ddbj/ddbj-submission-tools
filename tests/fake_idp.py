"""A scripted stand-in for the Keycloak endpoints that ddbj_auth.oauth.post_form talks to."""
import base64
import json

ISS = "https://idp.example.org/realms/master"


def make_jwt(**claims) -> str:
    def b64(obj):
        return base64.urlsafe_b64encode(json.dumps(obj).encode()).decode().rstrip("=")
    return f"{b64({'alg': 'none'})}.{b64(claims)}.sig"


class FakeIdp:
    def __init__(self):
        self.urls: list[str] = []
        self.requests: list[dict] = []
        self.device = {"device_code": "dc", "user_code": "ABCD-EFGH", "verification_uri": f"{ISS}/device",
                       "verification_uri_complete": f"{ISS}/device?user_code=ABCD-EFGH", "expires_in": 600, "interval": 5}
        self.token_queue: list[dict] = []
        self.fail_revoke: Exception | None = None

    def post_form(self, url: str, data: dict) -> dict:
        self.urls.append(url)
        self.requests.append(dict(data))
        if url.endswith("/auth/device"):
            return dict(self.device)
        if url.endswith("/revoke"):
            if self.fail_revoke:
                raise self.fail_revoke
            return {}
        return self.token_queue.pop(0)
