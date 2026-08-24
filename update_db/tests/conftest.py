"""Shared machinery for the `update_db` guard suite.

Run from the repository root:

    python -m pytest update_db/tests

Two mechanics that every test here inherits, both measured while prototyping this suite
(issue #351) rather than assumed:

1. **`update_db/scripts/` is not a package.** `utils_update.py` is imported by its two siblings
   as a flat `from utils_update import ...`, which only resolves with `scripts/` on `sys.path`.
   We load it by file path instead, so nothing here depends on the working directory or on
   `scripts/` becoming a package, and no `sys.path` mutation leaks into another test session.

2. **No interpreter on a development machine can import `utils_update.py`.** It imports `bs4` at
   module scope, and `bs4` is declared as a dependency nowhere in this repository -- the
   `update_db` dependencies exist only inside the Docker image the tooling runs in. So `bs4` and
   `requests` are stubbed in `sys.modules` before the load.

   **The stubs raise on use.** That is the load-bearing half. A guard that quietly exercised a
   fake network would be worse than no guard, so if a code path under test ever reaches
   `BeautifulSoup` or `requests.get`, this suite fails loudly instead of measuring a mock.
   `pandas` and `numpy` are the real thing; every runner that can host this suite has them.

Between them, these two keep the whole suite inside `pip install pytest pandas` -- the exact
dependency set the repository's other CI guards already install -- with no network, no Docker
and no emulated image.
"""

from __future__ import annotations

import gzip
import importlib.util
import sys
import types
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
UTILS_UPDATE = REPO_ROOT / "update_db" / "scripts" / "utils_update.py"


def _stub_absent_imports() -> None:
    """Stand in for the module-scope imports no development machine has, raising if used."""
    if "bs4" not in sys.modules:
        bs4 = types.ModuleType("bs4")

        def _BeautifulSoup(*args, **kwargs):  # noqa: N802 - mirrors the real name
            raise AssertionError(
                "the update_db guard suite stubs bs4, and a code path under test reached "
                "BeautifulSoup. Either the test is exercising the wrong function or a "
                "primitive has started scraping; do not replace this with a working mock."
            )

        bs4.BeautifulSoup = _BeautifulSoup
        sys.modules["bs4"] = bs4

    if "requests" not in sys.modules:
        requests = types.ModuleType("requests")

        def _get(*args, **kwargs):
            raise AssertionError(
                "the update_db guard suite stubs requests, and a code path under test "
                "attempted a download. This suite must never touch the network."
            )

        requests.get = _get
        sys.modules["requests"] = requests


def load_utils_update():
    """Load `update_db/scripts/utils_update.py` from its path, under a private module name."""
    _stub_absent_imports()
    spec = importlib.util.spec_from_file_location("utils_update_under_test", UTILS_UPDATE)
    module = importlib.util.module_from_spec(spec)
    sys.modules["utils_update_under_test"] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="session")
def utils():
    """The module under test."""
    return load_utils_update()


# ---------------------------------------------------------------------------------------------
# Staged databases that pass every check.
# ---------------------------------------------------------------------------------------------
#
# Each builder writes what its update script writes on a good day, read off the script rather
# than off the validator -- otherwise the fixture would agree with the code by construction and
# the checks would be testing themselves. The tests then break exactly one thing at a time.
#
#   update_clinvar_funcotator.sh  ->  clinvar/<build>/clinvar_<YYYYMMDD>.vcf + .idx + config
#   get_new_hgnc.sh               ->  hgnc/<build>/hgnc_<%b%d%Y>.tsv + hgnc.config, no index
#   update_dbsnp.sh               ->  dbsnp/<build>/<GCF accession>.gz + .tbi + dbSNP.config

GOOD_VERSION = {
    "clinvar": "20250812",
    "hgnc": "Aug202026",
    "dbsnp": "b156",
    "acmg_rec": "v3.2",
}

_DBSNP_DATA = "GCF_000001405.40.gz"


def _config(**fields) -> str:
    """A Funcotator config: the fields that matter, in the real files' `key = value` shape."""
    body = "".join(f"{key} = {value}\n" for key, value in fields.items())
    return body + "type = vcf\nncbi_build_version = \n"


def build_stage(root: Path, name: str, builds=("hg38",)) -> Path:
    """Write a staged `name` database under `root` that satisfies all five checks."""
    for build in builds:
        build_dir = root / name / build
        build_dir.mkdir(parents=True)

        if name == "clinvar":
            version = GOOD_VERSION["clinvar"]
            data = build_dir / f"clinvar_{version}.vcf"
            data.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\nchr1\t100\n")
            (build_dir / f"{data.name}.idx").write_bytes(b"\x00index\n")
            (build_dir / "clinvar_vcf.config").write_text(
                _config(name="ClinVar_VCF", version=version, src_file=data.name)
            )

        elif name == "hgnc":
            version = GOOD_VERSION["hgnc"]
            data = build_dir / f"hgnc_{version}.tsv"
            data.write_text("HGNC ID\tApproved Symbol\nHGNC:5\tA1BG\n")
            (build_dir / "hgnc.config").write_text(
                _config(name="HGNC", version=version, src_file=data.name)
            )

        elif name == "acmg_rec":
            version = GOOD_VERSION["acmg_rec"]
            data = build_dir / f"acmg_{version}_Aug202026_test_cleaned.txt"
            data.write_text("Disease_Name\tgene\nHereditary breast cancer\tBRCA1\n")
            (build_dir / "acmg_rec.config").write_text(
                _config(name="ACMG_recommendation", version=version, src_file=data.name)
            )

        elif name == "dbsnp":
            version = GOOD_VERSION["dbsnp"]
            data = build_dir / _DBSNP_DATA
            with gzip.open(data, "wt") as handle:
                handle.write("##fileformat=VCFv4.2\n#CHROM\tPOS\nchr1\t100\n")
            (build_dir / f"{data.name}.tbi").write_bytes(b"\x00tabix\n")
            # The trailing space is the real script's: it writes `version = b$version `.
            (build_dir / "dbSNP.config").write_text(
                _config(name="dbSNP", version=f"{version} ", src_file=data.name)
            )

        else:
            raise AssertionError(
                f"no fixture builder for {name!r}. A new STAGE_SPEC entry needs one here, or "
                "its checks are covered by nothing."
            )
    return root / name


def build_live(root: Path, name: str, builds=("hg38",), version="20250101") -> Path:
    """A live (already installed) database, for the germline drift check to compare."""
    config = {
        "clinvar": "clinvar_vcf.config",
        "hgnc": "hgnc.config",
        "dbsnp": "dbSNP.config",
        "acmg_rec": "acmg_rec.config",
    }
    for build in builds:
        build_dir = root / name / build
        build_dir.mkdir(parents=True)
        (build_dir / config[name]).write_text(_config(name=name, version=version))
    return root / name
