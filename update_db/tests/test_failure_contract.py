"""The failure-contract primitives: can each of them actually refuse?

These are the shared pieces landed by issue #360, ahead of the two route tickets that convert
the steps. The discipline throughout is one this repository has been bitten by often enough to
make a house rule: **a guard is worth nothing until it has been seen to fail.** So every check
gets a stage that passes, and then the same stage with exactly one thing broken -- and the test
asserts not merely that it failed but that the message named what was wrong. A validator that
refuses everything would satisfy a test that only checked for refusal.

The other recurring failure mode here is the vacuous parametrising table: a per-step spec that
quietly stops covering a step, and nothing says so. `STAGE_SPEC` is exactly that shape, so it
gets its own controls -- an unspecified step must be a hard error, a misspelt key must be a hard
error, and a step in the table with no fixture below must be a test failure.
"""

from __future__ import annotations

import ast
import gzip
import re

import pytest

from conftest import REPO_ROOT, GOOD_VERSION, build_live, build_stage, load_utils_update

UTILS = load_utils_update()
STEPS = sorted(UTILS.STAGE_SPEC)


# =============================================================================================
# StepFailed
# =============================================================================================


def test_step_failed_is_a_plain_exception_carrying_a_message():
    """Ruling 7: every anticipated failure raises this, with a message that says what happened.

    It deliberately does not subclass `ValueError` or `FileExistsError` -- the two the tool
    raises today -- because `run_steps` distinguishes an anticipated failure from a bug in the
    tool by type, and inheriting from a builtin would let an unrelated `ValueError` be reported
    as though the step had diagnosed it.
    """
    error = UTILS.StepFailed("clinvar hg38: staged clinvar_.vcf is 0 bytes")
    assert isinstance(error, Exception)
    assert not isinstance(error, (ValueError, OSError))
    assert "clinvar_.vcf" in str(error)


# =============================================================================================
# STAGE_SPEC, and the controls that keep it from going vacuous
# =============================================================================================


def test_the_shipped_stage_spec_is_well_formed():
    assert UTILS.check_stage_spec() is UTILS.STAGE_SPEC


@pytest.mark.parametrize("name", STEPS)
def test_every_stage_spec_entry_spells_every_key(name):
    """None-where-absent, so a check that does not apply is a decision rather than an omission.

    `hgnc` has no index because it is a simpleXSV data source with nothing for gatk to index.
    That has to be visible in the table as `index=None`, not inferable from the key's absence.
    """
    assert set(UTILS.STAGE_SPEC[name]) == UTILS._STAGE_SPEC_KEYS


def test_a_misspelt_spec_key_is_a_hard_error():
    """The typo this control exists for: `confg=` leaves `config` unset, and a spec with no
    config silently skips the version-agreement check -- the one that catches BSD `sed`'s
    no-op. Silently skipping a check is the failure mode; refusing loudly is the fix."""
    typo = {"clinvar": dict(UTILS.STAGE_SPEC["clinvar"])}
    typo["clinvar"]["confg"] = typo["clinvar"].pop("config")
    with pytest.raises(UTILS.StepFailed) as raised:
        UTILS.check_stage_spec(typo)
    assert "confg" in str(raised.value)


def test_a_spec_missing_a_key_is_a_hard_error():
    incomplete = {"clinvar": {k: v for k, v in UTILS.STAGE_SPEC["clinvar"].items() if k != "index"}}
    with pytest.raises(UTILS.StepFailed) as raised:
        UTILS.check_stage_spec(incomplete)
    assert "index" in str(raised.value)


def test_a_step_with_no_spec_is_a_hard_error_not_a_pass(tmp_path):
    """The anti-vacuity control named in the contract, and the reason it is in the contract
    rather than in a style guide. If validating an unknown step returned quietly, then a new
    destructive step -- or a renamed old one -- would sail past the guard with the guard still
    green, which is precisely the defect this whole effort is about."""
    with pytest.raises(UTILS.StepFailed) as raised:
        UTILS._require_valid_stage("cosmic", ["hg38"], stage_root=str(tmp_path))
    assert "STAGE_SPEC" in str(raised.value)


@pytest.mark.parametrize("name", STEPS)
def test_every_specified_step_has_a_fixture(name, tmp_path):
    """A step can be added to `STAGE_SPEC` and tested by nothing. This makes that a red test:
    `build_stage` refuses a name it has no builder for."""
    build_stage(tmp_path, name, builds=("hg38",))
    assert (tmp_path / name / "hg38").is_dir()


# =============================================================================================
# _require_valid_stage -- the five checks of ruling 2
# =============================================================================================


@pytest.mark.parametrize("name", STEPS)
def test_a_good_stage_passes_and_reports_its_version(name, tmp_path):
    """The control that keeps every test below honest. Without it a validator that refused
    everything would look like a working guard."""
    build_stage(tmp_path, name, builds=("hg38", "hg19"))
    versions = UTILS._require_valid_stage(
        name, ["hg38", "hg19"], stage_root=str(tmp_path)
    )
    assert versions == {"hg38": GOOD_VERSION[name], "hg19": GOOD_VERSION[name]}


@pytest.mark.parametrize("name", STEPS)
def test_check_1_a_missing_build_directory_fails(name, tmp_path):
    """The container's failure shape: `/bin/sh` is dash, the script dies at its first bashism,
    and nothing was ever staged."""
    build_stage(tmp_path, name, builds=("hg38",))
    with pytest.raises(UTILS.StepFailed) as raised:
        UTILS._require_valid_stage(name, ["hg38", "hg19"], stage_root=str(tmp_path))
    assert "hg19" in str(raised.value)
    assert "hg38" not in str(raised.value)


def test_check_2_the_measured_unresolved_clinvar_filename_fails(tmp_path):
    """The exact state a failed ClinVar download leaves, measured in the container: a 0-byte
    `clinvar_.vcf` and a config whose version never got substituted, because `$date` was empty.

    This is the shape that a `returncode` check cannot see -- the script exits 0 -- and the one
    that, once `sh` becomes `bash`, replaces the live database with an empty annotation source
    Funcotator will happily accept.
    """
    stage = tmp_path / "clinvar" / "hg38"
    stage.mkdir(parents=True)
    (stage / "clinvar_.vcf").write_text("")
    (stage / "clinvar_vcf.config").write_text("version = \nsrc_file = clinvar_.vcf\n")

    with pytest.raises(UTILS.StepFailed) as raised:
        UTILS._require_valid_stage("clinvar", ["hg38"], stage_root=str(tmp_path))
    message = str(raised.value)
    assert "clinvar_" in message
    assert "version never resolved" in message
    assert "live database is untouched" in message


@pytest.mark.parametrize("name", STEPS)
def test_check_2_an_empty_data_file_fails(name, tmp_path):
    """A data file with the right name and no content. Distinct from the case above, where the
    name itself is wrong -- here the version resolved and the download still produced nothing."""
    build_stage(tmp_path, name, builds=("hg38",))
    data = _data_file(tmp_path, name, "hg38")
    data.write_bytes(b"")
    with pytest.raises(UTILS.StepFailed) as raised:
        UTILS._require_valid_stage(name, ["hg38"], stage_root=str(tmp_path))
    assert "0 bytes" in str(raised.value)


@pytest.mark.parametrize("name", ["clinvar", "dbsnp"])
def test_check_3_a_missing_index_fails(name, tmp_path):
    """`gatk` or `tabix` absent from PATH: the VCF is there, the index is not, and Funcotator
    cannot read it. `hgnc` is excluded because it genuinely has no index."""
    build_stage(tmp_path, name, builds=("hg38",))
    suffix = UTILS.STAGE_SPEC[name]["index"]
    _data_file(tmp_path, name, "hg38").with_name(
        _data_file(tmp_path, name, "hg38").name + suffix
    ).unlink()
    with pytest.raises(UTILS.StepFailed) as raised:
        UTILS._require_valid_stage(name, ["hg38"], stage_root=str(tmp_path))
    assert "never indexed" in str(raised.value)


def test_check_3_hgnc_needs_no_index():
    """Guarding the exclusion above: if `hgnc` ever grows an index this test fails and the
    parametrisation two tests up has to be revisited, rather than silently under-covering."""
    assert UTILS.STAGE_SPEC["hgnc"]["index"] is None


@pytest.mark.parametrize("name", STEPS)
def test_check_4_a_blank_config_version_fails(name, tmp_path):
    build_stage(tmp_path, name, builds=("hg38",))
    config = tmp_path / name / "hg38" / UTILS.STAGE_SPEC[name]["config"]
    config.write_text(config.read_text().replace(GOOD_VERSION[name], ""))
    with pytest.raises(UTILS.StepFailed) as raised:
        UTILS._require_valid_stage(name, ["hg38"], stage_root=str(tmp_path))
    assert "version" in str(raised.value)


@pytest.mark.parametrize("name", ["clinvar", "hgnc"])
def test_check_4_a_config_that_disagrees_with_the_filename_fails(name, tmp_path):
    """The version is in both places for these two, so they must agree. `dbsnp` keeps its
    version only in the config -- see the test below for what stands in there."""
    build_stage(tmp_path, name, builds=("hg38",))
    config = tmp_path / name / "hg38" / UTILS.STAGE_SPEC[name]["config"]
    other = {"clinvar": "20240101", "hgnc": "Jan012024"}[name]
    config.write_text(
        re.sub(rf"^version = .*$", f"version = {other}", config.read_text(), flags=re.M)
    )
    with pytest.raises(UTILS.StepFailed) as raised:
        UTILS._require_valid_stage(name, ["hg38"], stage_root=str(tmp_path))
    assert "disagree" in str(raised.value)


def test_check_4_a_config_version_of_the_wrong_shape_fails(tmp_path):
    """dbSNP's measured failure shape, and the reason `version` is a per-step pattern rather
    than a permissive one. `update_dbsnp.sh` writes `version = b$version `, taking `$version`
    from a `grep dbSNP_BUILD_ID` over the downloaded VCF's header. When the download produced
    nothing that grep finds nothing, and the config says `version = b` -- a blank version
    wearing a plausible shape, which any pattern loose enough to cover `20250812`, `Aug202026`
    and `b156` at once would wave straight through.

    Found by mutation: nothing else in this file reaches the shape check on a config version.
    """
    build_stage(tmp_path, "dbsnp", builds=("hg38",))
    config = tmp_path / "dbsnp" / "hg38" / "dbSNP.config"
    config.write_text(re.sub(r"^version = .*$", "version = b", config.read_text(), flags=re.M))
    with pytest.raises(UTILS.StepFailed) as raised:
        UTILS._require_valid_stage("dbsnp", ["hg38"], stage_root=str(tmp_path))
    assert "not a resolved version" in str(raised.value)


def test_check_2_a_filename_version_of_the_wrong_shape_fails(tmp_path):
    """The same shape check on the other side, where the filename's capture is looser than what
    counts as resolved. This is not hypothetical: the contract's own sketch of this table
    proposed `hgnc_(?P<version>\\w+)\\.tsv`, which matches `hgnc_.tsv` -- so a spec written that
    way needs the version pattern to be the thing that refuses it.

    Found by mutation: with the shipped table every `data` capture is exactly as strict as its
    `version`, so this branch is reachable only through a spec like the one below.
    """
    # `config=None` as well as the loose capture, so that the config checks cannot mask this
    # one. Without that the test passed with the filename check switched off -- the config's
    # own version check was refusing the stage instead, and the assertion could not tell which.
    loose = {
        "hgnc": dict(
            UTILS.STAGE_SPEC["hgnc"], data=r"hgnc_(?P<version>\w*)\.tsv", config=None
        )
    }
    build_dir = tmp_path / "hgnc" / "hg38"
    build_dir.mkdir(parents=True)
    (build_dir / "hgnc_x.tsv").write_text("HGNC ID\tSymbol\nHGNC:5\tA1BG\n")
    with pytest.raises(UTILS.StepFailed) as raised:
        UTILS._require_valid_stage(
            "hgnc", ["hg38"], stage_root=str(tmp_path), spec_table=loose
        )
    assert "not a resolved version" in str(raised.value)


def test_check_4_dbsnp_keeps_its_version_only_in_the_config():
    """Recorded as a test because it is the one place the contract's ruling 2 is not literally
    expressible. dbSNP's data file is named for the assembly accession, and its version is the
    dbSNP build id, so there is no filename version for the config to agree with. `src_file`
    is what carries the agreement instead -- see the next test."""
    assert "(?P<version>" not in UTILS.STAGE_SPEC["dbsnp"]["data"]
    assert UTILS.STAGE_SPEC["dbsnp"]["config"] == "dbSNP.config"


@pytest.mark.parametrize("name", STEPS)
def test_check_4_a_config_pointing_at_the_wrong_file_fails(name, tmp_path):
    """What the macOS BSD `sed -i 's/DATE/.../'` no-op leaves behind: BSD `sed` reads the
    expression as the `-i` backup suffix, exits 1, and leaves the file UNCHANGED -- and nothing
    checks it, so the config still names `clinvar_DATE.vcf` while the data file beside it has a
    real date. Two small file reads catch it. This is the check that works for all three steps,
    dbSNP included."""
    build_stage(tmp_path, name, builds=("hg38",))
    config = tmp_path / name / "hg38" / UTILS.STAGE_SPEC[name]["config"]
    config.write_text(
        re.sub(r"^src_file = .*$", "src_file = something_else.vcf", config.read_text(), flags=re.M)
    )
    with pytest.raises(UTILS.StepFailed) as raised:
        UTILS._require_valid_stage(name, ["hg38"], stage_root=str(tmp_path))
    assert "points at nothing" in str(raised.value)


@pytest.mark.parametrize("name", STEPS)
def test_check_5_a_header_only_data_file_fails(name, tmp_path):
    """Non-empty, correctly named, correctly indexed, correctly configured -- and carrying no
    variants. Funcotator accepts it as a valid data source with zero records, so every variant
    would be annotated with no evidence from this database, silently, into clinical reports."""
    build_stage(tmp_path, name, builds=("hg38",))
    data = _data_file(tmp_path, name, "hg38")

    # Strip the fixture's data records and keep its header, whatever shape that step's header is.
    # Derived from the spec rather than named per step, so a step added to STAGE_SPEC later cannot
    # end up silently exercising the wrong header convention -- which is what happened when
    # `acmg_rec` was added by #353: this test wrote a VCF header for it, and a table whose header
    # is a single column row read the second `#CHROM` line as a data record and passed.
    compressed = name == "dbsnp"
    opener = (lambda mode: gzip.open(data, mode + "t")) if compressed else None
    text = opener("r").read() if opener else data.read_text()
    if UTILS.STAGE_SPEC[name]["header"] == "hash":
        header = "".join(line for line in text.splitlines(True) if line.startswith("#"))
    else:
        header = text.splitlines(True)[0]
    assert header and header != text, "the fixture had no data records to strip"
    if opener:
        with gzip.open(data, "wt") as handle:
            handle.write(header)
    else:
        data.write_text(header)

    with pytest.raises(UTILS.StepFailed) as raised:
        UTILS._require_valid_stage(name, ["hg38"], stage_root=str(tmp_path))
    assert "not one data record" in str(raised.value)


def test_the_record_check_reads_only_the_header(tmp_path):
    """Ruling 2 declined record counts because they cost a full scan of a multi-GB VCF. This is
    the property that makes check 5 affordable: it stops at the first data line."""
    path = tmp_path / "big.vcf"
    path.write_text("##one\n##two\n#CHROM\tPOS\nchr1\t1\n" + "chr1\t2\n" * 10000)
    assert UTILS.first_data_record(str(path), "hash") == "chr1\t1"


def test_the_record_check_reads_the_bgzipped_dbsnp_vcf(tmp_path):
    path = tmp_path / "GCF_000001405.40.gz"
    with gzip.open(path, "wt") as handle:
        handle.write("##fileformat=VCFv4.2\n#CHROM\tPOS\nchr1\t100\n")
    assert UTILS.first_data_record(str(path), "hash") == "chr1\t100"


def test_every_problem_with_a_stage_is_reported_at_once(tmp_path):
    """One verdict per step, so the operator should not have to fix and re-run to see the next
    thing wrong."""
    build_stage(tmp_path, "clinvar", builds=("hg38",))
    build_dir = tmp_path / "clinvar" / "hg38"
    (build_dir / f"clinvar_{GOOD_VERSION['clinvar']}.vcf.idx").unlink()
    config = build_dir / "clinvar_vcf.config"
    config.write_text(
        re.sub(r"^src_file = .*$", "src_file = wrong.vcf", config.read_text(), flags=re.M)
    )
    with pytest.raises(UTILS.StepFailed) as raised:
        UTILS._require_valid_stage("clinvar", ["hg38"], stage_root=str(tmp_path))
    message = str(raised.value)
    assert "never indexed" in message
    assert "points at nothing" in message


def _data_file(root, name, build):
    """The staged data file, found the way the validator finds it."""
    pattern = UTILS.STAGE_SPEC[name]["data"]
    build_dir = root / name / build
    matches = [p for p in sorted(build_dir.iterdir()) if re.fullmatch(pattern, p.name)]
    assert len(matches) == 1, f"fixture for {name} did not produce one data file: {matches}"
    return matches[0]


# =============================================================================================
# requested_builds
# =============================================================================================


@pytest.mark.parametrize(
    "build, expected",
    [("both", ["hg38", "hg19"]), ("hg38", ["hg38"]), ("hg19", ["hg19"])],
)
def test_requested_builds(build, expected):
    assert UTILS.requested_builds(build) == expected


# =============================================================================================
# The germline drift check -- ruling 6
# =============================================================================================


def test_germline_in_sync_passes(tmp_path):
    somatic, germline = tmp_path / "db", tmp_path / "db_germline"
    build_live(somatic, "clinvar", version="20250101")
    build_live(germline, "clinvar", version="20250101")
    notes = UTILS._require_germline_in_sync(
        "clinvar", ["hg38"], str(somatic), str(germline)
    )
    assert notes["hg38"] == "in sync at 20250101"


def test_germline_drift_fails_before_the_destructive_write(tmp_path):
    """The unenforced assumption this turns into a checked precondition: germline receives the
    same staged directory somatic does and gets no backup of its own, so `backup_dir` holds
    somatic's old version and restores germline only if the two were in sync. Nothing guaranteed
    that, and if they had drifted, overwriting germline destroys content no backup holds."""
    somatic, germline = tmp_path / "db", tmp_path / "db_germline"
    build_live(somatic, "clinvar", version="20250101")
    build_live(germline, "clinvar", version="20240202")
    with pytest.raises(UTILS.StepFailed) as raised:
        UTILS._require_germline_in_sync("clinvar", ["hg38"], str(somatic), str(germline))
    message = str(raised.value)
    assert "20250101" in message and "20240202" in message
    assert "somatic backup cannot restore" in message


def test_an_absent_germline_copy_is_not_drift(tmp_path):
    """Absence is not disagreement. A build with no germline copy has nothing to lose, and
    failing there would break a host whose install is deliberately partial."""
    somatic, germline = tmp_path / "db", tmp_path / "db_germline"
    build_live(somatic, "clinvar", version="20250101")
    germline.mkdir()
    notes = UTILS._require_germline_in_sync(
        "clinvar", ["hg38"], str(somatic), str(germline)
    )
    assert "nothing to lose" in notes["hg38"]


def test_a_step_that_writes_no_germline_is_vacuously_in_sync(tmp_path):
    notes = UTILS._require_germline_in_sync("hgnc", ["hg38"], str(tmp_path), None)
    assert "nothing to keep in sync" in notes["hg38"]


def test_the_drift_check_refuses_a_database_it_cannot_compare(tmp_path):
    """If a database has no config there is no version to compare, so declaring it derivable
    would be an assertion rather than a proof. Ruling 6 says the step must prove it."""
    spec = {"gnomad": dict(UTILS.STAGE_SPEC["clinvar"])}
    spec["gnomad"]["config"] = None
    with pytest.raises(UTILS.StepFailed) as raised:
        UTILS._require_germline_in_sync(
            "gnomad", ["hg38"], str(tmp_path), str(tmp_path), spec_table=spec
        )
    assert "different proof" in str(raised.value)


# =============================================================================================
# The verdicts -- ruling 4
# =============================================================================================


def test_the_verdict_vocabulary_is_exactly_four_words_and_success_is_not_one():
    """`SUCCESS` is retired, and there is deliberately no writer for it. It was doing four
    different jobs across 21 sites -- 'updated', 'already the latest', 'no update needed', and
    at `update_clinvar` 'did not happen at all' -- and since the destination is *no step reports
    SUCCESS for work that did not happen*, the word itself is where the ambiguity lived."""
    assert UTILS.VERDICTS == ("UPDATED", "CURRENT", "SKIPPED", "FAILED")
    assert not hasattr(UTILS.RunLog, "success")
    assert not hasattr(UTILS, "log_success")


def test_the_four_verdict_lines_and_the_result_line(tmp_path):
    log = UTILS.RunLog(str(tmp_path / "update.log"))
    log.updated("clinvar", "hg38", "20250114", "20250812")
    log.current("dbsnp", "hg38", "b156")
    log.skipped("cosmic", "update disabled, format change")
    log.failed("clinvar", "staged clinvar_.vcf is 0 bytes; live database untouched", build="hg19")
    exit_code = log.result()

    lines = (tmp_path / "update.log").read_text().splitlines()
    assert lines[0] == "UPDATED: clinvar hg38 20250114 -> 20250812"
    assert lines[1] == "CURRENT: dbsnp hg38 (b156) - already the latest"
    assert lines[2] == "SKIPPED: cosmic - update disabled, format change"
    assert lines[3] == (
        "FAILED:  clinvar hg19 - staged clinvar_.vcf is 0 bytes; live database untouched"
    )
    assert lines[4].startswith("RESULT:  1 updated, 1 current, 1 skipped, 1 failed")
    assert exit_code == 1


def test_warning_and_info_survive_only_as_indented_detail(tmp_path):
    """Which is what `update_dbsnp` already does when it echoes the child script's `SKIPPED:`
    lines indented beneath its own `FAILED:`."""
    log = UTILS.RunLog(str(tmp_path / "update.log"))
    log.failed("dbsnp", "both builds missing", detail=["wget exited 4", "no version resolved"])
    lines = (tmp_path / "update.log").read_text().splitlines()
    assert lines[0] == "FAILED:  dbsnp - both builds missing"
    assert lines[1] == "  detail: wget exited 4"
    assert lines[2] == "  detail: no version resolved"


def test_the_result_line_names_the_skips(tmp_path):
    """Ruling 5: the gap between what was asked for and what happened stays visible without
    being an error."""
    log = UTILS.RunLog(str(tmp_path / "update.log"))
    log.updated("clinvar", "hg38", "a", "b")
    log.skipped("dbsnp", "not installed", build="hg19")
    assert log.result() == 0
    result = (tmp_path / "update.log").read_text().splitlines()[-1]
    assert "0 failed" in result
    assert "(skipped: dbsnp hg19 not installed)" in result


def test_the_result_line_is_derived_from_the_lines_actually_written(tmp_path):
    """The load-bearing rule, as a property of the code rather than a convention: the tally IS
    the verdict lines, so a run cannot be summarised as having updated something without an
    `UPDATED:` line having been written for it."""
    log = UTILS.RunLog(str(tmp_path / "update.log"))
    assert log.tally() == {"UPDATED": 0, "CURRENT": 0, "SKIPPED": 0, "FAILED": 0}
    log.updated("clinvar", "hg38", "a", "b")
    log.updated("clinvar", "hg19", "a", "b")
    assert log.tally()["UPDATED"] == 2


# =============================================================================================
# The step driver -- rulings 3, 4 and 5
# =============================================================================================


def test_a_failed_step_does_not_abort_the_run(tmp_path):
    """Ruling 3, and the thing today's tool gets wrong in three of its four failure idioms: an
    uncaught `ValueError` from `update_dbsnp` means `acmg_rec` never runs, and on the Annovar
    route -- which has no `try`/`except` anywhere -- a `FileExistsError` in
    `update_all_cancervar` means `update_all_intervar` never runs at all."""
    log = UTILS.RunLog(str(tmp_path / "update.log"))
    ran = []

    def failing():
        raise UTILS.StepFailed("clinvar: the staged database did not pass validation")

    def later():
        ran.append("later")
        log.updated("hgnc", "hg38", "old", "new")

    exit_code = UTILS.run_steps([("clinvar", failing), ("hgnc", later)], log, ["hg38"])

    assert ran == ["later"], "a failed step aborted the run"
    assert exit_code == 1
    body = (tmp_path / "update.log").read_text()
    assert "FAILED:  clinvar - clinvar: the staged database did not pass validation" in body
    assert "UPDATED: hgnc hg38 old -> new" in body


def test_an_unexpected_exception_becomes_a_failed_verdict_with_its_traceback(tmp_path):
    """The backstop. A genuine bug in the tool should surface as a `FAILED` verdict carrying
    real information -- not as an opaque abort, and above all not as the bare `except:` in
    `update_funcotator.py` that reports a curl download problem for every possible HGNC failure,
    including the ones where no download was ever reached."""
    log = UTILS.RunLog(str(tmp_path / "update.log"))

    def buggy():
        raise KeyError("hg38")

    assert UTILS.run_steps([("hgnc", buggy)], log, ["hg38"]) == 1
    body = (tmp_path / "update.log").read_text()
    assert "FAILED:  hgnc - unexpected KeyError: 'hg38'" in body
    assert "detail:" in body
    assert "Traceback" in body


def test_a_message_less_assertion_error_still_says_something_useful(tmp_path):
    """`update_hgnc` fails today via `assert len(latest_data.dropna()) > 0` -- an
    `AssertionError` with no message at all, which a generic handler renders as
    `FAILED: hgnc - unexpected AssertionError:` with nothing after the colon. Honest, but barely
    better than the wrong message it replaces, which is why ruling 7 converts those asserts to
    `StepFailed`. Until they are converted, the backstop at least keeps the traceback."""
    log = UTILS.RunLog(str(tmp_path / "update.log"))

    def buggy():
        assert False  # noqa: B011 - reproducing update_hgnc's message-less assert on purpose

    UTILS.run_steps([("hgnc", buggy)], log, ["hg38"])
    body = (tmp_path / "update.log").read_text()
    assert "unexpected AssertionError" in body
    assert "assert False" in body, "the traceback is the only diagnostic such a failure has"


def test_a_step_that_records_no_verdict_at_all_is_a_failure(tmp_path):
    """Ruling 4's load-bearing rule, enforced. A step that returns having said nothing has told
    the operator nothing, which is indistinguishable from a step whose work silently did not
    happen -- the exact defect this effort exists for."""
    log = UTILS.RunLog(str(tmp_path / "update.log"))
    assert UTILS.run_steps([("clinvar", lambda: None)], log, ["hg38"]) == 1
    assert "without recording a verdict" in (tmp_path / "update.log").read_text()


def test_a_per_build_skip_is_benign(tmp_path):
    """Ruling 5: a host whose install is deliberately partial must not exit 1 every night."""
    log = UTILS.RunLog(str(tmp_path / "update.log"))

    def partial():
        log.skipped("dbsnp", "no existing database found", build="hg19")
        log.updated("dbsnp", "hg38", "b155", "b156")

    assert UTILS.run_steps([("dbsnp", partial)], log, ["hg38", "hg19"]) == 0


def test_a_step_that_skipped_every_requested_build_is_a_failure(tmp_path):
    """The other half of ruling 5, and it generalises `update_dbsnp`'s existing rule rather than
    inventing one: one build missing is a warning, both missing raises."""
    log = UTILS.RunLog(str(tmp_path / "update.log"))

    def nothing():
        log.skipped("dbsnp", "no existing database found", build="hg38")
        log.skipped("dbsnp", "no existing database found", build="hg19")

    assert UTILS.run_steps([("dbsnp", nothing)], log, ["hg38", "hg19"]) == 1
    assert "no requested build was updated" in (tmp_path / "update.log").read_text()


def test_a_whole_step_skip_stays_benign(tmp_path):
    """A disabled caller or absent credentials is a deliberate decision about the whole step,
    not a build that failed to update, so it must not escalate. COSMIC and GENCODE are both in
    this shape today."""
    log = UTILS.RunLog(str(tmp_path / "update.log"))

    def disabled():
        log.skipped("cosmic", "update disabled, format change")

    assert UTILS.run_steps([("cosmic", disabled)], log, ["hg38", "hg19"]) == 0
    assert "no requested build" not in (tmp_path / "update.log").read_text()


def test_a_clean_run_exits_zero(tmp_path):
    """The control on the exit code: the driver must not be unable to succeed."""
    log = UTILS.RunLog(str(tmp_path / "update.log"))

    def good():
        log.current("clinvar", "hg38", "20250812")

    assert UTILS.run_steps([("clinvar", good)], log, ["hg38"]) == 0


def test_the_driver_writes_exactly_one_result_line(tmp_path):
    log = UTILS.RunLog(str(tmp_path / "update.log"))
    UTILS.run_steps(
        [("a", lambda: log.current("a", "hg38", "1")), ("b", lambda: log.current("b", "hg38", "1"))],
        log,
        ["hg38"],
    )
    body = (tmp_path / "update.log").read_text()
    assert body.count("RESULT:") == 1


# =============================================================================================
# End to end: the stage a failed ClinVar download leaves must never reach the live database
# =============================================================================================


def test_the_measured_clinvar_failure_produces_a_failed_verdict_and_a_nonzero_exit(tmp_path):
    """The whole contract in one run, against the failure shape measured in the container.

    This is the shape #353 converted `update_clinvar` to -- validate the stage, and only then
    touch the live database. Before it, the destructive block ran unconditionally and the step
    logged `SUCCESS: clinvar updated successfully.` at exit 0.
    """
    stage_root = tmp_path / "work"
    (stage_root / "clinvar" / "hg38").mkdir(parents=True)
    (stage_root / "clinvar" / "hg38" / "clinvar_.vcf").write_text("")

    live = tmp_path / "db"
    build_live(live, "clinvar", version="20250101")
    before = (live / "clinvar" / "hg38" / "clinvar_vcf.config").read_text()

    log = UTILS.RunLog(str(tmp_path / "update.log"))

    def clinvar_step():
        versions = UTILS._require_valid_stage(
            "clinvar", ["hg38"], stage_root=str(stage_root)
        )
        raise AssertionError(f"validation must not have passed, but returned {versions}")

    exit_code = UTILS.run_steps([("clinvar", clinvar_step)], log, ["hg38"])

    assert exit_code == 1
    body = (tmp_path / "update.log").read_text()
    assert "FAILED:  clinvar" in body
    assert "SUCCESS" not in body
    assert (live / "clinvar" / "hg38" / "clinvar_vcf.config").read_text() == before


# =============================================================================================
# Tree-wide: `SUCCESS` is gone from update_db, not just from the step you happen to be reading
# =============================================================================================


def _script_sources() -> dict:
    """Every Python file under `update_db/scripts/`, by relative path."""
    scripts = sorted((REPO_ROOT / "update_db" / "scripts").glob("*.py"))
    assert scripts, "found no scripts to audit, so this guard is checking nothing"
    return {str(p.relative_to(REPO_ROOT)): p.read_text() for p in scripts}


def test_nothing_under_update_db_scripts_can_write_a_success_verdict():
    """The map's own destination, as one assertion over the whole tool.

    `SUCCESS` appeared 21 times across the three Python files doing four different jobs --
    "updated", "already the latest", "no update needed", and at `update_clinvar` "did not happen at
    all". Since the destination is *no `update_db` run reports SUCCESS for work that did not
    happen*, the word itself was where the ambiguity lived, so ruling 4 retired it rather than
    qualifying it. `RunLog` deliberately has no writer for it.

    This could not ship until both routes had been converted -- with either half unconverted it
    would have failed on the other's 21 occurrences -- so #353 and #354 each carried a
    route-scoped version, and this is the tree-wide one that replaces them both. It is what stops
    the word coming back in a step nobody is looking at.

    String literals only, read off the AST: the prose in these files explains at length what
    `SUCCESS` used to do and why it went, and an explanation is not an occurrence. That
    distinction is the reason this is not a `grep`.
    """
    offenders = {}
    for path, text in _script_sources().items():
        tree = ast.parse(text, filename=path)
        docstrings = {
            id(node.body[0].value)
            for node in ast.walk(tree)
            if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef))
            and node.body
            and isinstance(node.body[0], ast.Expr)
            and isinstance(node.body[0].value, ast.Constant)
        }
        for node in ast.walk(tree):
            if (
                isinstance(node, ast.Constant)
                and isinstance(node.value, str)
                and id(node) not in docstrings
                and "SUCCESS" in node.value
            ):
                offenders.setdefault(path, []).append((node.lineno, node.value[:60]))

    assert not offenders, (
        "`SUCCESS` is retired by ruling 4, but these string literals still carry it -- a verdict "
        f"the log vocabulary has no meaning for: {offenders}"
    )


def test_the_success_audit_would_notice_one():
    """The anti-vacuity control for the audit above.

    A guard that reports nothing is indistinguishable from a guard that looks at nothing, and this
    one now passes over a tree where the word is genuinely absent -- which is exactly when it stops
    being obvious that it still works. So: feed the same walk a module that does write the verdict,
    and a docstring that merely mentions it, and check it separates them.
    """
    source = 'def step():\n    """Once wrote SUCCESS."""\n    log("SUCCESS: clinvar updated")\n'
    tree = ast.parse(source)
    literals = [
        node.value
        for node in ast.walk(tree)
        if isinstance(node, ast.Constant) and isinstance(node.value, str)
    ]
    written = [v for v in literals if "SUCCESS" in v and not v.startswith("Once wrote")]
    assert written == ["SUCCESS: clinvar updated"], (
        f"the walk no longer distinguishes a written verdict from prose about one: {literals}"
    )
