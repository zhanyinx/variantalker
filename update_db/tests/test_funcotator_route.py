"""The Funcotator route, measured against a POPULATED database with no network and no Docker.

Issue #353. Every test here replaces the update script with a stub that reproduces one measured
failure shape and then asserts what happened to the live databases -- because the way this bug
survived was being exercised against empty throwaway directories, where a run that destroys
everything and a run that works are indistinguishable: both exit 0 and both leave the log looking
plausible.

The stubs are not inventions. Each one reproduces a shape that was measured in
`yinxiu/variantalker_db:latest` with `--network none`, or on macOS with a script that exits
non-zero:

* `dash` -- the container's `/bin/sh`. `update_clinvar_funcotator.sh` has no `function` keyword,
  so dash does not fail to parse it; it fails at `[[ ... ]]` at *runtime*, with every conditional
  false. The script exits **0** having created nothing but an empty `clinvar/`, and that empty
  directory was copied over both live ClinVar databases while the log said
  `SUCCESS: clinvar updated successfully.`
* `blank version` -- what #339's `sh`->`bash` fix produces on its own. The download fails, `$date`
  is never set, and the stage carries a 0-byte `clinvar_.vcf` beside a config whose `version =` is
  empty. Funcotator accepts that as a valid VCF data source with zero records, so this is the
  SUBTLER outcome, not the fixed one -- which is why the interpreter change and the validation had
  to land in the same commit.
* `no index` -- `gatk` missing from PATH.
* `single build` -- a correct `--build hg38` run. The script honours the flag and stages one build;
  it was the Python around it that replaced the whole directory.

Two independent facts these tests pin, both of which the route got wrong in a way no log revealed:
a failed step must leave the live database **byte-identical**, and a single-build run must leave the
build it was not asked for **byte-identical**.
"""

from __future__ import annotations

import ast
import os
import shutil
import stat
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))

from conftest import GOOD_VERSION, load_utils_update  # noqa: E402

UTILS = load_utils_update()

BUILDS = ("hg38", "hg19")


# ---------------------------------------------------------------------------------------------
# A populated installation, and stub update scripts that reproduce the measured shapes.
# ---------------------------------------------------------------------------------------------


OLD_VERSION = {"clinvar": "20250101", "hgnc": "Jan012025", "dbsnp": "b155", "acmg_rec": "v3.1"}

# What the INSTALLED database's data file is called. It has to be named the way each step looks for
# it -- `update_hgnc` globs `hgnc_*.tsv` against the live tree -- or a test would measure a missing
# database rather than the behaviour it meant to.
OLD_DATA = {
    "clinvar": f"clinvar_{OLD_VERSION['clinvar']}.vcf",
    "hgnc": f"hgnc_{OLD_VERSION['hgnc']}.tsv",
    "dbsnp": "GCF_000001405.39.gz",
    "acmg_rec": f"acmg_{OLD_VERSION['acmg_rec']}_Jan012025_test_cleaned.txt",
}


def _write_live(root: Path, name: str, builds=BUILDS, version=None):
    """A live database with real-shaped files in it, so a destructive run has something to lose."""
    version = version or OLD_VERSION[name]
    config_name = UTILS.STAGE_SPEC[name]["config"]
    data_name = OLD_DATA[name] if version == OLD_VERSION[name] else f"{name}_{version}.data"
    for build in builds:
        build_dir = root / name / build
        build_dir.mkdir(parents=True, exist_ok=True)
        if name == "hgnc":
            # Deliberately different content from the staged copy, so `assert_frame_equal` sees a
            # change: a live copy identical to the stage is the CURRENT case, tested separately.
            (build_dir / data_name).write_text("HGNC ID\tApproved Symbol\nHGNC:5\tA1BG_OLD\n")
        else:
            (build_dir / data_name).write_text(f"{name} {build} {version} payload\n")
        (build_dir / config_name).write_text(
            f"name = {name}\nversion = {version}\nsrc_file = {data_name}\n"
        )


def _snapshot(root: Path) -> dict:
    """Every file under `root`, by relative path, with its bytes -- so 'untouched' can be asserted
    rather than inferred from a directory still existing. An empty directory where a database used
    to be is the exact shape this map exists for, and `os.path.isdir` cannot tell the difference."""
    out = {}
    for path in sorted(root.rglob("*")):
        if path.is_file():
            out[str(path.relative_to(root))] = path.read_bytes()
    return out


SCRIPT_PATH = {
    "clinvar": "update_clinvar/update_clinvar_funcotator.sh",
    "hgnc": "update_hgnc/get_new_hgnc.sh",
    "dbsnp": "update_dbsnp/update_dbsnp.sh",
}

DATA_NAME = {
    "clinvar": f"clinvar_{GOOD_VERSION['clinvar']}.vcf",
    "hgnc": f"hgnc_{GOOD_VERSION['hgnc']}.tsv",
    "dbsnp": "GCF_000001405.40.gz",
}

# How each script writes its data file. Per database rather than one shape for all three, because
# each is a different kind of file and the checks read them: the dbSNP VCF is bgzipped (a valid
# gzip stream), HGNC is a tab table whose first line is its header, ClinVar is a plain VCF. A rig
# that wrote VCF bytes into `hgnc_*.tsv` would fail validation for the wrong reason and the test
# would still look green.
DATA_COMMAND = {
    "clinvar": "printf '##fileformat=VCFv4.2\\n#CHROM\\tPOS\\nchr1\\t100\\n' > '{path}'",
    "hgnc": "printf 'HGNC ID\\tApproved Symbol\\nHGNC:5\\tA1BG\\n' > '{path}'",
    "dbsnp": "printf '##fileformat=VCFv4.2\\n#CHROM\\tPOS\\nchr1\\t100\\n' | gzip > '{path}'",
}


def _stub_body(name: str, shape: str, builds) -> str:
    """The shell body of a stub update script, reproducing one measured failure shape."""
    if shape == "hard_fail":
        # Nothing staged and a non-zero exit -- what `set -eo pipefail` now produces on a failed
        # download, and what a stub scriptdir reproduces on macOS with no Docker at all.
        return "echo 'download failed' >&2\nexit 4\n"

    lines = [f"mkdir -p {name}"]
    if shape == "dash":
        # Exits 0 having created nothing but the top-level stage directory. Measured under dash:
        # `[[ ... ]]` is a runtime command-not-found, so every build conditional was false.
        return "\n".join(lines) + "\n"

    config = UTILS.STAGE_SPEC[name]["config"]
    index = UTILS.STAGE_SPEC[name]["index"]
    for build in builds:
        d = f"{name}/{build}"
        lines.append(f"mkdir -p {d}")
        if shape == "blank_version":
            # The download failed, `$date` never resolved: a 0-byte `clinvar_.vcf` beside a config
            # whose `version =` is empty. #339's interpreter fix alone produces exactly this.
            lines.append(f": > '{d}/clinvar_.vcf'")
            lines.append(
                f"printf 'name = {name}\\nversion = \\nsrc_file = clinvar_.vcf\\n' > '{d}/{config}'"
            )
            continue

        data = DATA_NAME[name]
        lines.append(DATA_COMMAND[name].format(path=f"{d}/{data}"))
        # `no_index` is gatk missing from PATH: data and config fine, nothing indexed them.
        if index and shape != "no_index":
            lines.append(f"printf 'index' > '{d}/{data}{index}'")
        lines.append(
            f"printf 'name = {name}\\nversion = {GOOD_VERSION[name]}\\nsrc_file = {data}\\n' "
            f"> '{d}/{config}'"
        )
    return "\n".join(lines) + "\n"


def _make_scriptdir(root: Path, name: str, shape: str, builds=BUILDS, extra: str = "") -> str:
    """A stub `scriptdir` whose script reproduces one measured shape for one database."""
    script = root / "scriptdir" / SCRIPT_PATH[name]
    script.parent.mkdir(parents=True, exist_ok=True)
    body = _stub_body(name, shape, builds)
    if shape == "hard_fail":
        script.write_text("#!/bin/bash\n" + body)
    else:
        script.write_text("#!/bin/bash\nset -eo pipefail\n" + body + extra + "\nexit 0\n")
    script.chmod(script.stat().st_mode | stat.S_IEXEC)
    return str(root / "scriptdir")


@pytest.fixture
def installation(tmp_path, monkeypatch):
    """A populated somatic + germline installation, with the working directory set to the scratch
    space the update scripts stage into -- which is how the real tool is invoked."""
    somatic = tmp_path / "somatic"
    germline = tmp_path / "germline"
    backup = tmp_path / "backup" / "20260820"
    scratch = tmp_path / "scratch"
    for d in (somatic, germline, backup, scratch):
        d.mkdir(parents=True)
    monkeypatch.chdir(scratch)
    return {
        "somatic": somatic,
        "germline": germline,
        "backup": backup,
        "root": tmp_path,
        "log": tmp_path / "update.log",
    }


def _run_log(installation):
    return UTILS.RunLog(str(installation["log"]))


def _verdicts(installation):
    return installation["log"].read_text() if installation["log"].exists() else ""


# ---------------------------------------------------------------------------------------------
# A failed step must leave the live databases byte-identical.
# ---------------------------------------------------------------------------------------------


@pytest.mark.parametrize("shape", ["dash", "blank_version", "no_index", "hard_fail"])
def test_a_failed_clinvar_update_leaves_both_live_databases_byte_identical(
    shape, installation
):
    """The measurement that matters, in the form it was originally made: a populated fixture with
    real-shaped ClinVar files in BOTH the somatic and germline trees.

    Before this ticket, all four of these shapes replaced both databases with an empty directory
    -- germline with no backup ever taken, because `update_clinvar` copied the stage to germline
    and backed up only somatic -- and wrote `SUCCESS: clinvar updated successfully.` at exit 0."""
    _write_live(installation["somatic"], "clinvar")
    _write_live(installation["germline"], "clinvar")
    before_somatic = _snapshot(installation["somatic"])
    before_germline = _snapshot(installation["germline"])

    scriptdir = _make_scriptdir(installation["root"], "clinvar", shape)
    run_log = _run_log(installation)
    exit_code = UTILS.run_steps(
        [
            (
                "clinvar",
                lambda: UTILS.update_clinvar(
                    db_dir=str(installation["somatic"]),
                    backup_dir=str(installation["backup"]),
                    run_log=run_log,
                    scriptdir=scriptdir,
                    db_germline_dir=str(installation["germline"]),
                    build="both",
                ),
            )
        ],
        run_log,
        list(BUILDS),
    )

    assert _snapshot(installation["somatic"]) == before_somatic
    assert _snapshot(installation["germline"]) == before_germline
    assert exit_code == 1
    log = _verdicts(installation)
    assert "FAILED:" in log
    assert "UPDATED:" not in log
    assert "SUCCESS" not in log


def test_the_failure_says_which_check_refused_the_stage(installation):
    """A `FAILED` line that does not say what was checked is the failure mode `StepFailed` exists to
    replace -- `update_hgnc`'s message-less `assert`, and the bare `except:` that reported a curl
    download problem when no download had been reached."""
    _write_live(installation["somatic"], "clinvar")
    _write_live(installation["germline"], "clinvar")
    scriptdir = _make_scriptdir(installation["root"], "clinvar", "blank_version")
    run_log = _run_log(installation)
    UTILS.run_steps(
        [
            (
                "clinvar",
                lambda: UTILS.update_clinvar(
                    db_dir=str(installation["somatic"]),
                    backup_dir=str(installation["backup"]),
                    run_log=run_log,
                    scriptdir=scriptdir,
                    db_germline_dir=str(installation["germline"]),
                    build="both",
                ),
            )
        ],
        run_log,
        list(BUILDS),
    )
    log = _verdicts(installation)
    assert "clinvar_.vcf" in log, "the failure does not name the file it refused"
    assert "live database is untouched" in log


def test_a_good_update_replaces_both_builds_and_backs_up_the_old_one(installation):
    """The other half of the same property: when the stage IS a database, it is installed -- so
    these tests cannot pass by refusing everything."""
    _write_live(installation["somatic"], "clinvar")
    _write_live(installation["germline"], "clinvar")
    scriptdir = _make_scriptdir(installation["root"], "clinvar", "good")
    run_log = _run_log(installation)
    exit_code = UTILS.run_steps(
        [
            (
                "clinvar",
                lambda: UTILS.update_clinvar(
                    db_dir=str(installation["somatic"]),
                    backup_dir=str(installation["backup"]),
                    run_log=run_log,
                    scriptdir=scriptdir,
                    db_germline_dir=str(installation["germline"]),
                    build="both",
                ),
            )
        ],
        run_log,
        list(BUILDS),
    )

    assert exit_code == 0
    new = GOOD_VERSION["clinvar"]
    for tree in ("somatic", "germline"):
        for build in BUILDS:
            installed = installation[tree] / "clinvar" / build
            assert (installed / f"clinvar_{new}.vcf").is_file()
            assert UTILS.read_config_field(str(installed / "clinvar_vcf.config"), "version") == new
    # The old version is recoverable, for both trees, because germline was byte-identical to it.
    assert (installation["backup"] / "clinvar" / "hg38" / OLD_DATA["clinvar"]).is_file()
    log = _verdicts(installation)
    assert f"UPDATED: clinvar hg38 {OLD_VERSION['clinvar']} -> {new}" in log
    assert "RESULT:" in log


# ---------------------------------------------------------------------------------------------
# A single-build run must leave the build it was not asked for byte-identical.
# ---------------------------------------------------------------------------------------------


@pytest.mark.parametrize("requested,untouched", [("hg38", "hg19"), ("hg19", "hg38")])
def test_a_single_build_clinvar_update_leaves_the_other_build_alone(
    requested, untouched, installation
):
    """`--build hg38` means "update hg38 and leave hg19 alone", settled with the person who added
    the flag. The whole-directory `rm -rf` + `cp -r` was simply the wrong granularity for it: a
    successful single-build run deleted the other build from BOTH trees, silently, and no later run
    noticed -- the "no existing database" checks all looked at the build that WAS requested."""
    _write_live(installation["somatic"], "clinvar")
    _write_live(installation["germline"], "clinvar")
    keep_somatic = _snapshot(installation["somatic"] / "clinvar" / untouched)
    keep_germline = _snapshot(installation["germline"] / "clinvar" / untouched)

    scriptdir = _make_scriptdir(installation["root"], "clinvar", "good", builds=(requested,))
    run_log = _run_log(installation)
    exit_code = UTILS.run_steps(
        [
            (
                "clinvar",
                lambda: UTILS.update_clinvar(
                    db_dir=str(installation["somatic"]),
                    backup_dir=str(installation["backup"]),
                    run_log=run_log,
                    scriptdir=scriptdir,
                    db_germline_dir=str(installation["germline"]),
                    build=requested,
                ),
            )
        ],
        run_log,
        [requested],
    )

    assert exit_code == 0
    assert _snapshot(installation["somatic"] / "clinvar" / untouched) == keep_somatic
    assert _snapshot(installation["germline"] / "clinvar" / untouched) == keep_germline
    assert (
        installation["somatic"] / "clinvar" / requested / f"clinvar_{GOOD_VERSION['clinvar']}.vcf"
    ).is_file()
    # Ruling 3: only the build being replaced is backed up, so a single-build run does not copy a
    # whole-genome VCF for a build that is not moving.
    assert not (installation["backup"] / "clinvar" / untouched).exists()
    # Ruling 5 as amended by the build-scope decision: the unrequested build gets no log line. Its
    # absence from the verdicts is the point -- five "not requested" skips per run would be noise,
    # and would put five entries in the RESULT line's skip list.
    assert untouched not in _verdicts(installation)


def test_a_single_build_dbsnp_update_leaves_the_other_build_alone(installation):
    """The same defect in the one step that already inspected something. Its `files_created` gate
    is an `or` -- it asks whether ANYTHING was staged, not whether every live build was accounted
    for -- so `--build hg38` passed the gate and replaced the whole `dbsnp` directory in both
    trees."""
    _write_live(installation["somatic"], "dbsnp")
    _write_live(installation["germline"], "dbsnp")
    keep = _snapshot(installation["somatic"] / "dbsnp" / "hg19")

    # The real script prints SKIPPED only from inside its own build conditional, so a --build hg38
    # run never mentions hg19 at all: the step cannot tell "hg19 was not asked for" from "hg19 is
    # fine". Under the per-build swap it does not need to.
    scriptdir = _make_scriptdir(installation["root"], "dbsnp", "good", builds=("hg38",))
    run_log = _run_log(installation)
    exit_code = UTILS.run_steps(
        [
            (
                "dbsnp",
                lambda: UTILS.update_dbsnp(
                    run_log=run_log,
                    db_dir=str(installation["somatic"]),
                    backup_dir=str(installation["backup"]),
                    scriptdir=scriptdir,
                    db_germline_dir=str(installation["germline"]),
                    build="hg38",
                ),
            )
        ],
        run_log,
        ["hg38"],
    )

    assert exit_code == 0
    assert _snapshot(installation["somatic"] / "dbsnp" / "hg19") == keep
    assert f"UPDATED: dbsnp hg38 {OLD_VERSION['dbsnp']} -> {GOOD_VERSION['dbsnp']}" in _verdicts(
        installation
    )


def test_dbsnp_reports_a_script_side_skip_per_build_without_failing_the_run(installation):
    """`update_dbsnp` was the only step that read a machine-readable signal out of its child, and
    that survives -- as the house pattern rather than as a special case. A per-build skip is benign:
    a host whose install is deliberately partial must not exit 1 every night."""
    _write_live(installation["somatic"], "dbsnp")
    _write_live(installation["germline"], "dbsnp")
    scriptdir = _make_scriptdir(
        installation["root"],
        "dbsnp",
        "good",
        builds=("hg38",),
        extra="echo 'SKIPPED: hg19 dbSNP update - no existing database found in db/hg19/'",
    )
    run_log = _run_log(installation)
    exit_code = UTILS.run_steps(
        [
            (
                "dbsnp",
                lambda: UTILS.update_dbsnp(
                    run_log=run_log,
                    db_dir=str(installation["somatic"]),
                    backup_dir=str(installation["backup"]),
                    scriptdir=scriptdir,
                    db_germline_dir=str(installation["germline"]),
                    build="both",
                ),
            )
        ],
        run_log,
        list(BUILDS),
    )

    log = _verdicts(installation)
    assert "SKIPPED: dbsnp hg19" in log
    assert "UPDATED: dbsnp hg38" in log
    assert exit_code == 0, "a per-build skip must not affect the exit code"
    assert "(skipped: dbsnp hg19" in log, "RESULT must name the skips"


def test_dbsnp_skipping_every_requested_build_is_a_failure(installation):
    """The escalation this generalises from `update_dbsnp`'s own rule: one build missing is a
    warning, every requested build missing means the step did nothing at all."""
    _write_live(installation["somatic"], "dbsnp")
    scriptdir = _make_scriptdir(
        installation["root"],
        "dbsnp",
        "dash",
        extra=(
            "echo 'SKIPPED: hg38 dbSNP update - no existing database found'\n"
            "echo 'SKIPPED: hg19 dbSNP update - no existing database found'"
        ),
    )
    run_log = _run_log(installation)
    exit_code = UTILS.run_steps(
        [
            (
                "dbsnp",
                lambda: UTILS.update_dbsnp(
                    run_log=run_log,
                    db_dir=str(installation["somatic"]),
                    backup_dir=str(installation["backup"]),
                    scriptdir=scriptdir,
                    build="both",
                ),
            )
        ],
        run_log,
        list(BUILDS),
    )
    assert exit_code == 1
    assert "no requested build was updated" in _verdicts(installation)


# ---------------------------------------------------------------------------------------------
# HGNC: the build that had never been able to work.
# ---------------------------------------------------------------------------------------------


@pytest.mark.parametrize("build", ["hg19", "hg38"])
def test_hgnc_can_update_either_build(build, installation):
    """`update_hgnc --build hg19` had never worked: both globs named `hg38`, and the STAGED one
    raised a bare `IndexError` because the script stages only the build it was asked for. Reported
    as "There was a problem with the curl download", which was wrong twice over.

    Fixing only the message would have turned a crash into a permanent, well-worded FAILED on every
    `--build hg19` run -- documenting the bug instead of fixing it. So this asserts the update
    actually happens."""
    _write_live(installation["somatic"], "hgnc")
    other = "hg38" if build == "hg19" else "hg19"
    keep = _snapshot(installation["somatic"] / "hgnc" / other)

    scriptdir = _make_scriptdir(installation["root"], "hgnc", "good", builds=(build,))
    run_log = _run_log(installation)
    exit_code = UTILS.run_steps(
        [
            (
                "hgnc",
                lambda: UTILS.update_hgnc(
                    db_dir=str(installation["somatic"]),
                    backup_dir=str(installation["backup"]),
                    run_log=run_log,
                    scriptdir=scriptdir,
                    build=build,
                ),
            )
        ],
        run_log,
        [build],
    )

    assert exit_code == 0
    assert f"UPDATED: hgnc {build} {OLD_VERSION['hgnc']} -> {GOOD_VERSION['hgnc']}" in _verdicts(
        installation
    )
    assert (
        installation["somatic"] / "hgnc" / build / f"hgnc_{GOOD_VERSION['hgnc']}.tsv"
    ).is_file()
    assert _snapshot(installation["somatic"] / "hgnc" / other) == keep


def test_hgnc_with_no_live_copy_for_a_build_skips_that_build(installation):
    """"No existing database found" is a skip everywhere else in the tool, so it is a skip here --
    and because it is the only build requested, `_settle_step` escalates it to a step-level
    failure, which is the outcome the contract's ruling 7 asked for."""
    _write_live(installation["somatic"], "hgnc", builds=("hg38",))
    scriptdir = _make_scriptdir(installation["root"], "hgnc", "good", builds=("hg19",))
    run_log = _run_log(installation)
    exit_code = UTILS.run_steps(
        [
            (
                "hgnc",
                lambda: UTILS.update_hgnc(
                    db_dir=str(installation["somatic"]),
                    backup_dir=str(installation["backup"]),
                    run_log=run_log,
                    scriptdir=scriptdir,
                    build="hg19",
                ),
            )
        ],
        run_log,
        ["hg19"],
    )
    log = _verdicts(installation)
    assert "SKIPPED: hgnc hg19" in log
    assert "install it manually first" in log
    assert exit_code == 1


# ---------------------------------------------------------------------------------------------
# Germline is derivable from somatic, and the step proves it before overwriting.
# ---------------------------------------------------------------------------------------------


def test_germline_drift_stops_the_update_before_anything_is_destroyed(installation):
    """Germline receives the same staged bytes somatic does and gets no backup of its own, so
    `backup_dir` is a germline restore path only if the two were in sync. Nothing guaranteed that.
    Two config reads turn the unenforced assumption into a checked precondition."""
    _write_live(installation["somatic"], "clinvar", version="20250101")
    _write_live(installation["germline"], "clinvar", version="20240101")
    before = _snapshot(installation["germline"])

    scriptdir = _make_scriptdir(installation["root"], "clinvar", "good")
    run_log = _run_log(installation)
    exit_code = UTILS.run_steps(
        [
            (
                "clinvar",
                lambda: UTILS.update_clinvar(
                    db_dir=str(installation["somatic"]),
                    backup_dir=str(installation["backup"]),
                    run_log=run_log,
                    scriptdir=scriptdir,
                    db_germline_dir=str(installation["germline"]),
                    build="both",
                ),
            )
        ],
        run_log,
        list(BUILDS),
    )

    assert exit_code == 1
    assert _snapshot(installation["germline"]) == before
    log = _verdicts(installation)
    assert "not derivable from the somatic one" in log
    assert "somatic is 20250101 but germline is 20240101" in log


# ---------------------------------------------------------------------------------------------
# The vocabulary, at the level of a whole run.
# ---------------------------------------------------------------------------------------------


def test_a_leftover_stage_from_a_crashed_run_does_not_become_the_new_database(installation):
    """A stale staged file beside a fresh one is the one way a validated stage could still be
    wrong. The stage is cleared before the script runs, and two files matching the same pattern are
    refused as ambiguous rather than resolved by guessing."""
    _write_live(installation["somatic"], "clinvar")
    _write_live(installation["germline"], "clinvar")
    leftover = Path("clinvar/hg38")
    leftover.mkdir(parents=True)
    (leftover / "clinvar_19990101.vcf").write_text("stale\n")

    scriptdir = _make_scriptdir(installation["root"], "clinvar", "good")
    run_log = _run_log(installation)
    exit_code = UTILS.run_steps(
        [
            (
                "clinvar",
                lambda: UTILS.update_clinvar(
                    db_dir=str(installation["somatic"]),
                    backup_dir=str(installation["backup"]),
                    run_log=run_log,
                    scriptdir=scriptdir,
                    db_germline_dir=str(installation["germline"]),
                    build="both",
                ),
            )
        ],
        run_log,
        list(BUILDS),
    )
    assert exit_code == 0
    installed = installation["somatic"] / "clinvar" / "hg38"
    assert not (installed / "clinvar_19990101.vcf").exists()
    assert (installed / f"clinvar_{GOOD_VERSION['clinvar']}.vcf").is_file()


def test_no_funcotator_step_can_write_the_word_success():
    """`SUCCESS` is retired on this route. It was doing four different jobs -- "updated", "already
    the latest", "no update needed", and at `update_clinvar` "did not happen at all" -- and since
    the destination is *no step reports SUCCESS for work that did not happen*, the word itself is
    where the ambiguity lived.

    Scoped to the functions this ticket converted. The Annovar route still has its own, which is
    #354's work, so a tree-wide assertion would ship red -- and the moment it can pass, it should
    be made, by whichever of the two lands second."""
    converted = {
        "check_oncotator",
        "check_dna_repair_genes",
        "check_oreganno",
        "update_cosmic",
        "update_dbsnp",
        "update_gencode",
        "update_clinvar",
        "update_hgnc",
        "update_acmg_rec",
    }
    tree = ast.parse(Path(UTILS.__file__).read_text())
    seen, offending = set(), []
    for node in ast.walk(tree):
        if not isinstance(node, ast.FunctionDef) or node.name not in converted:
            continue
        seen.add(node.name)
        docstring = ast.get_docstring(node, clean=False)
        for child in ast.walk(node):
            # String CONSTANTS only, and never the docstring -- a comment or a docstring that
            # describes the retired vocabulary is documentation, not a verdict being written. An
            # earlier draft of this test scanned the text and failed on its own explanation.
            if isinstance(child, ast.Constant) and isinstance(child.value, str):
                if child.value == docstring or "SUCCESS" not in child.value:
                    continue
                offending.append(f"{node.name}: {child.value[:60]!r}")

    missing = converted - seen
    assert not missing, f"these converted steps have been renamed, so nothing covers them: {missing}"
    assert not offending, "converted steps must not write SUCCESS: " + "; ".join(offending)


def test_the_entry_point_exits_non_zero_when_a_step_failed():
    """The machine-readable half of the contract. Nothing in this repository parses the log -- the
    only references to it are the two lines that write it -- so the log is a human-legibility
    contract and the exit code is the signal. `update_funcotator.py` had none: it returned None
    and exited 0 whatever happened."""
    entry = Path(UTILS.__file__).parent / "update_funcotator.py"
    text = entry.read_text()
    assert "sys.exit(main())" in text, "the entry point must carry the run's exit code"
    assert "return run_steps(" in text, "main must return the driver's exit code"

    # The bare `except:` is gone, checked against the parse tree rather than the text -- the text
    # still explains what used to be there, and should.
    tree = ast.parse(text)
    bare = [
        node.lineno
        for node in ast.walk(tree)
        if isinstance(node, ast.ExceptHandler) and node.type is None
    ]
    assert not bare, (
        f"a bare `except:` is back at line(s) {bare}; it is what reported a curl download problem "
        "for every possible HGNC failure, including one where no download was attempted"
    )


def test_the_stub_scripts_are_a_faithful_rig():
    """The anti-vacuity control for this file. If the stub scriptdir did not actually run -- a
    permissions slip, a bad shebang, a path typo -- every test above would still pass, because a
    script that never runs stages nothing and "nothing was staged" is what most of them assert.

    So: prove the rig can produce a database that PASSES validation, and prove the stub really is
    executed by bash."""
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        cwd = os.getcwd()
        try:
            os.chdir(root)
            scriptdir = _make_scriptdir(root, "clinvar", "good")
            result = UTILS.run_update_script(
                f"{scriptdir}/{SCRIPT_PATH['clinvar']}", ["--build", "both"]
            )
            assert result.returncode == 0, result.stderr
            versions = UTILS._require_valid_stage("clinvar", list(BUILDS))
            assert versions == {b: GOOD_VERSION["clinvar"] for b in BUILDS}
        finally:
            os.chdir(cwd)

    # And a failing stub really does fail, rather than being silently unrunnable.
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        scriptdir = _make_scriptdir(root, "clinvar", "hard_fail")
        result = UTILS.run_update_script(f"{scriptdir}/{SCRIPT_PATH['clinvar']}")
        assert result.returncode == 4, "the stub scriptdir is not being executed"


def test_the_route_runs_its_scripts_under_bash():
    """Restates the interpreter fix as a property of this route specifically, and keeps working
    after the tree-wide xfail guard is retired: whatever `run_update_script` does, it must not hand
    a bash script to a POSIX shell."""
    assert shutil.which("bash"), "this suite needs bash, which every runner that hosts it has"
    probe = Path(__file__).parent / "_probe_not_written"
    assert not probe.exists()
    with subprocess.Popen(
        ["bash", "-c", "echo ${BASH_VERSION}"], stdout=subprocess.PIPE, text=True
    ) as proc:
        out, _ = proc.communicate()
    assert out.strip(), "bash did not identify itself"

    import inspect

    source = inspect.getsource(UTILS.run_update_script)
    assert '"bash"' in source
    assert '"sh"' not in source
