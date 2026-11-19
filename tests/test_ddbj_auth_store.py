import base64
import json
import os
import stat
import threading

from ddbj_auth.store import claims, config_dir, delete, list_apis, load, save


def _b64(obj) -> str:
    return base64.urlsafe_b64encode(json.dumps(obj).encode()).decode().rstrip("=")


def test_save_load_roundtrip_and_modes(tmp_path, monkeypatch):
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
    save("staging", "kura", {"refresh_token": "r"})
    assert load("staging", "kura") == {"refresh_token": "r"}
    assert stat.S_IMODE(os.stat(config_dir("staging")).st_mode) == 0o700
    assert stat.S_IMODE(os.stat(config_dir("staging") / "kura.json").st_mode) == 0o600


def test_default_dir_under_home(tmp_path, monkeypatch):
    monkeypatch.delenv("XDG_CONFIG_HOME", raising=False)
    monkeypatch.setenv("HOME", str(tmp_path))
    assert config_dir("production") == tmp_path / ".config" / "ddbj" / "production"


def test_load_missing_returns_none(tmp_path, monkeypatch):
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
    assert load("staging", "nope") is None


def test_delete_and_list(tmp_path, monkeypatch):
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
    save("staging", "b", {"refresh_token": "1"})
    save("staging", "a", {"refresh_token": "2"})
    assert list_apis("staging") == ["a", "b"]
    assert delete("staging", "a") and not delete("staging", "a")
    assert list_apis("staging") == ["b"]


def test_concurrent_saves_leave_valid_json(tmp_path, monkeypatch):
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
    threads = [threading.Thread(target=save, args=("staging", "kura", {"refresh_token": str(i)})) for i in range(20)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    assert load("staging", "kura")["refresh_token"] in {str(i) for i in range(20)}
    assert not [p for p in config_dir("staging").iterdir() if "tmp" in p.name]


def test_claims_decodes_payload_without_padding():
    assert claims("x." + _b64({"aud": ["kura"], "exp": 5}) + ".y") == {"aud": ["kura"], "exp": 5}
    assert claims("garbage") == {}
    assert claims("a.!!!.c") == {}
