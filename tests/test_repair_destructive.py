from ddbj_gff.repair import REGISTRY, Operation, driver


def _dummy(name="dummy-destructive"):
    return Operation(name, "test-only", requires_sequence=False,
                     detect=lambda doc, ctx: [], apply=lambda doc, ctx, sel: [],
                     destructive=True)


def test_existing_operations_are_not_destructive():
    for name in ("internal-stop-to-misc", "utr-absent-to-partial-mrna",
                 "missing-start-stop-to-partial-cds"):
        assert REGISTRY[name].destructive is False


def test_apply_all_names_excludes_destructive(monkeypatch):
    monkeypatch.setitem(REGISTRY, "dummy-destructive", _dummy())
    monkeypatch.setattr(driver, "DEFAULT_ORDER",
                        ["dummy-destructive", "internal-stop-to-misc"])
    assert driver.apply_all_names() == ["internal-stop-to-misc"]


def test_apply_all_names_skips_unregistered(monkeypatch):
    monkeypatch.setattr(driver, "DEFAULT_ORDER",
                        ["not-registered", "internal-stop-to-misc"])
    assert driver.apply_all_names() == ["internal-stop-to-misc"]


def test_list_marks_destructive(monkeypatch, capsys):
    from ddbj_gff.repair import cli
    monkeypatch.setitem(REGISTRY, "dummy-destructive", _dummy())
    assert cli.main(["--list"]) == 0
    out = capsys.readouterr().out
    assert "dummy-destructive" in out
    assert "[destructive]" in out
    assert "internal-stop-to-misc" in out
