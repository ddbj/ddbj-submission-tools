"""Package layout: everything lives under ddbj_submission_tools; command names never change."""
import importlib.util
import subprocess
import sys
import zipfile
from importlib.metadata import entry_points
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# The portal (nuc) calls these by name; never rename them.
COMMANDS = {"egapx2mss", "mss_builder", "mss2ff", "batch_wgs_builder", "tsa2mss", "ann_meta",
            "gff2mss", "ann_fix", "gff_from_tool", "gff_normalize", "gff_repair", "gff_validate",
            "flatfile2gff", "ddbj-auth"}


def test_commands_unchanged_and_loadable():
    eps = {ep.name: ep for ep in entry_points(group="console_scripts")
           if ep.dist and ep.dist.name == "ddbj-submission-tools"}
    assert set(eps) == COMMANDS
    for ep in eps.values():
        assert ep.value.startswith("ddbj_submission_tools."), ep.value
        assert callable(ep.load())


def test_no_top_level_packages_left():
    for old in ["common", "ann_fix", "ann_meta", "ddbj_gff", "ddbj_auth", "gff2mss", "mss_builder",
                "mss2ff", "egapx2mss", "tsa2mss", "batch_wgs_builder"]:
        assert importlib.util.find_spec(old) is None, old


def test_egapx_bin_dir_is_repo_bin():
    from ddbj_submission_tools.nuc.egapx2mss.asn_tools import DEFAULT_BIN_DIR
    assert DEFAULT_BIN_DIR == ROOT / "bin"


def test_wheel_contains_package_data(tmp_path):
    subprocess.run([sys.executable, "-m", "pip", "wheel", str(ROOT), "--no-deps", "--no-build-isolation",
                    "-q", "-w", str(tmp_path)], check=True)
    [wheel] = tmp_path.glob("*.whl")
    names = zipfile.ZipFile(wheel).namelist()
    assert "ddbj_submission_tools/nuc/common/data/source_qualifiers.json" in names
    assert "ddbj_submission_tools/nuc/common/data/submission_categories.json" in names
    assert any(n.startswith("ddbj_submission_tools/nuc/ddbj_gff/validate/data/") for n in names)
    assert not any(n.startswith("common/") for n in names)
