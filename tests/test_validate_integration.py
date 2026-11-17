from pathlib import Path
from ddbj_gff import parse
from ddbj_gff.errors import Severity
from ddbj_gff.validate import validate

FIX = Path(__file__).parent / "validate_fixtures"


def test_valid_insdc_fixture_has_no_errors():
    diags = validate(parse((FIX / "valid_insdc.gff3").read_text()))
    errors = [d for d in diags if d.severity == Severity.ERROR]
    assert errors == [], f"unexpected errors: {[(d.code, d.message) for d in errors]}"



