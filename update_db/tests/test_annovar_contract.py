"""The failure contract as it applies to the Annovar route (issue #354).

The Funcotator route can destroy an annotation database: its scripts stage into a directory and
Python copies that over the live tree, so a stage that is empty replaces a real database with
nothing. The Annovar route cannot do that -- the build is a filename prefix and the version is
part of the name, so every write lands on a new file and the previous database stays where it is.

Which means the failure this route actually has is a different one, and these tests are about
that one: **a database is only in use because two files NAME it**, the annotation script's
hardcoded `-protocol` list and the config's `database_names` line. Downloading a table changes
nothing until those names move. So the properties worth guarding are

* a name never moves to a table that was not written, and
* a name always moves when the log says it did, and
* the two names never move independently of each other.

Everything here runs offline with stub downloaders. `annotate_variation.pl` and `vt` are absent
from the container as well as from every development machine, so driving the individual functions
against fixtures is not a compromise for the sake of the test -- it is the only way this route can
be exercised at all, and the properties above do not need either tool to be observable.

    python -m pytest update_db/tests/test_annovar_contract.py
"""

from __future__ import annotations

import ast
import os
import shutil
import stat
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
from conftest import REPO_ROOT  # noqa: E402


# ---------------------------------------------------------------------------------------------
# Fixtures: what annovar's humandb actually holds.
# ---------------------------------------------------------------------------------------------
#
# Read off the producers rather than off the validator, so the fixtures cannot agree with the
# code by construction: `update_clinvar_annovar.sh` writes `<build>_<name>.txt` through
# `gatk VariantsToTable` and indexes it with `utils/index_annovar.pl`, whose last lines print
# `#BIN\t<bin size>\t<size of the indexed file>`; `annotate_variation.pl -downdb` writes the same
# pair, and `annotate_variation.pl:2328` reads that header back.

CLINVAR_HEADER = "#Chr\tStart\tEnd\tRef\tAlt\tCLNALLELEID\tCLNDN\tCLNDISDB\tCLNREVSTAT\tCLNSIG\n"
CLINVAR_RECORD = "1\t100\t100\tA\tG\t12345\tsome_disease\tMONDO:1\treviewed\tPathogenic\n"


def write_annovar_db(humandb: Path, build: str, database: str, *, index=True, body=None) -> Path:
    """A database table annovar can annotate from, with a matching index."""
    table = humandb / f"{build}_{database}.txt"
    table.write_text(CLINVAR_HEADER + CLINVAR_RECORD if body is None else body)
    if index:
        write_annovar_index(table)
    return table


def write_annovar_index(table: Path, declared_size=None) -> Path:
    """The index `index_annovar.pl` writes: a `#BIN` header naming the table's size."""
    size = table.stat().st_size if declared_size is None else declared_size
    index = Path(str(table) + ".idx")
    index.write_text(f"#BIN\t1000\t{size}\n1\t0\t0\t{size}\n")
    return index


def write_stub_downloader(path: Path, *, builds=("hg19", "hg38"), returncode=0) -> Path:
    """Stand in for `annotate_variation.pl -downdb`, which no machine here has.

    Called exactly as the real one is -- `-buildver <build> -downdb -webfrom annovar <db> <dir>`
    -- and writes the pair of files the real one leaves behind, but only for the builds it was
    told to succeed for. That is what makes the partial-download case reachable offline.
    """
    succeed = " ".join(builds)
    path.write_text(
        "#!/bin/bash\n"
        f'for ok in {succeed or "__none__"}; do\n'
        '  if [ "$2" = "$ok" ]; then\n'
        '    printf \'%b\' "'
        + CLINVAR_HEADER.replace("\t", "\\t").replace("\n", "\\n")
        + CLINVAR_RECORD.replace("\t", "\\t").replace("\n", "\\n")
        + '" > "$7/$2_$6.txt"\n'
        '    size=$(wc -c < "$7/$2_$6.txt" | tr -d " ")\n'
        '    printf \'#BIN\\t1000\\t%s\\n\' "$size" > "$7/$2_$6.txt.idx"\n'
        "  fi\n"
        "done\n"
        f"exit {returncode}\n"
    )
    path.chmod(path.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    return path


def annovar_table(rows) -> pd.DataFrame:
    """The frame `get_annovar_databases` scrapes, without the scrape."""
    return pd.DataFrame(
        [(build, name, "", "2026-01-01") for name, build in rows],
        columns=["build", "name", "explanation", "date"],
    )


@pytest.fixture
def humandb(tmp_path: Path) -> Path:
    path = tmp_path / "humandb"
    path.mkdir()
    return path


@pytest.fixture
def annotation_script(tmp_path: Path) -> Path:
    """A stand-in for CancerVar.py: the versioned names appear more than once, as they really do.

    In the real file each database version is named in several `-protocol` strings and again as a
    key in the `Funcanno_flgs` dictionary, which is why the substitution has always been global.
    """
    script = tmp_path / "CancerVar.py"
    script.write_text(
        "-protocol refGene,ensGene,knownGene,esp6500siv2_all,exac03,avsnp151,dbnsfp47a,"
        "dbscsnv11,dbnsfp47a_interpro,clinvar_20260102,cosmic98,icgc28,gnomad_genome\n"
        'flags = {"avsnp151": 0, "cosmic98": 0, "dbnsfp47a": 0, "icgc28": 0}\n'
        "-protocol refGene,avsnp151,dbnsfp47a,dbnsfp47a_interpro,clinvar_20260102,cosmic98,icgc28\n"
    )
    return script


@pytest.fixture
def config_file(tmp_path: Path) -> Path:
    """The CancerVar config's real `database_names` line, order and all."""
    path = tmp_path / "config.init.CancerVar"
    path.write_text(
        "[other]\n"
        "database_names = refGene ensGene knownGene esp6500siv2_all 1000g2015aug exac03 "
        "dbscsnv11 dbnsfp47a_interpro gnomad_genome avsnp151 dbnsfp47a icgc28 cosmic98 "
        "clinvar_20260102\n"
        "otherkey = value\n"
    )
    return path


NOT_TO_UPDATE = [
    "refGene",
    "ensGene",
    "knownGene",
    "esp6500siv2_all",
    "1000g2015aug",
    "exac03",
    "dbscsnv11",
    "gnomad_genome",
    "rmsk",
]


def verdicts(run_log):
    """`[(verdict, db, detail-free reason)]`, for asserting on what a run reported."""
    return [(entry[0], entry[1]) for entry in run_log.entries]


def database_names(config_file: Path) -> list:
    for line in config_file.read_text().splitlines():
        if line.startswith("database_names"):
            return line.split("=", 1)[1].split()
    raise AssertionError(f"no database_names line in {config_file}")


# ---------------------------------------------------------------------------------------------
# The index is an optimisation, and that is a measured claim about the consumer.
# ---------------------------------------------------------------------------------------------


def test_a_good_index_has_nothing_to_report(utils, humandb):
    table = write_annovar_db(humandb, "hg38", "clinvar_20260820")
    assert utils.describe_annovar_index(str(table)) is None


@pytest.mark.parametrize(
    "break_it, expected",
    [
        ("absent", "scan the whole table"),
        ("empty", "does not start with an index header"),
        ("truncated", "does not start with an index header"),
        ("stale", "out of date"),
    ],
)
def test_an_unusable_index_is_described_rather_than_raised(utils, humandb, break_it, expected):
    """Every one of these is a footnote, not a failure, and the reason is annovar's own code.

    `annotate_variation.pl:2324-2364` requires the index's first line to match `BIN\\t<int>\\t<int>`
    and the size it declares to equal the actual size of the table, and when either is wrong it
    prints *"ANNOVAR can still generate correct results without index file"* and scans the table
    instead. A missing index takes the same path. So refusing an update over its index would
    refuse an update annovar handles correctly -- which is why this route's validation is about
    the table, and the index only ever produces a detail line.

    The `truncated` case is the one that actually happens: the indexer's output is redirected into
    the `.idx`, so the shell creates the file whether or not the indexer ran.
    """
    table = write_annovar_db(humandb, "hg38", "clinvar_20260820", index=break_it != "absent")
    index = Path(str(table) + ".idx")
    if break_it == "empty":
        index.write_text("")
    elif break_it == "truncated":
        index.write_text("#BIN\t1000")
    elif break_it == "stale":
        write_annovar_index(table, declared_size=table.stat().st_size + 512)

    note = utils.describe_annovar_index(str(table))
    assert note is not None, f"a {break_it} index should be reported"
    assert expected in note


def test_a_missing_index_does_not_fail_the_database(utils, humandb):
    """The check that keeps the ruling above honest, at the level that decides the verdict."""
    write_annovar_db(humandb, "hg38", "avsnp152", index=False)
    note = utils.require_valid_annovar_db(str(humandb), "hg38", "avsnp152")
    assert note and "results are unaffected" in note


# ---------------------------------------------------------------------------------------------
# The table is the correctness artifact.
# ---------------------------------------------------------------------------------------------


def test_a_good_database_passes(utils, humandb):
    write_annovar_db(humandb, "hg19", "avsnp152")
    assert utils.require_valid_annovar_db(str(humandb), "hg19", "avsnp152") is None


@pytest.mark.parametrize(
    "shape, body, expected",
    [
        ("missing", None, "was not created"),
        ("empty", "", "0 bytes"),
        ("header only", CLINVAR_HEADER, "not one data record"),
    ],
)
def test_a_table_that_carries_no_data_is_a_failure(utils, humandb, shape, body, expected):
    """The three shapes a failed download leaves behind, and the third is the dangerous one.

    A header-only table is a database annovar reads without complaint and annotates nothing from,
    which is the same shape as the empty ClinVar VCF this whole map exists for -- the difference
    being that here it is caught before any name starts pointing at it.
    """
    if body is not None:
        write_annovar_db(humandb, "hg19", "avsnp152", body=body)
    with pytest.raises(utils.StepFailed) as raised:
        utils.require_valid_annovar_db(str(humandb), "hg19", "avsnp152")
    assert expected in str(raised.value)


# ---------------------------------------------------------------------------------------------
# Adoption: moving the name, and knowing whether it moved.
# ---------------------------------------------------------------------------------------------


def test_adoption_rewrites_every_occurrence(utils, annotation_script):
    count = utils.adopt_annovar_database(str(annotation_script), "avsnp151", "avsnp152")
    text = annotation_script.read_text()
    assert count == 3
    assert "avsnp151" not in text
    assert text.count("avsnp152") == 3
    # Untouched neighbours: a global replace of one database must not disturb another.
    assert text.count("clinvar_20260102") == 2


def test_adopting_what_is_already_adopted_is_not_a_failure(utils, annotation_script):
    utils.adopt_annovar_database(str(annotation_script), "avsnp151", "avsnp152")
    assert utils.adopt_annovar_database(str(annotation_script), "avsnp151", "avsnp152") == 0


def test_adoption_refuses_a_script_that_names_neither_version(utils, tmp_path):
    """A database downloaded but never named is a download nothing will ever read."""
    script = tmp_path / "CancerVar.py"
    script.write_text("-protocol refGene,knownGene\n")
    with pytest.raises(utils.StepFailed) as raised:
        utils.adopt_annovar_database(str(script), "avsnp151", "avsnp152")
    assert "names neither" in str(raised.value)


def test_adoption_refuses_a_rename_whose_old_name_is_inside_the_new_one(utils, tmp_path):
    """`dbnsfp47a` -> `dbnsfp47a_interpro` is this shape, and a global replace cannot express it:
    the result names both databases at once and a second run extends the name again. The two are
    kept apart today only because `get_dbname_version` splits them into separate entries, which
    is a property of the caller rather than of the substitution."""
    script = tmp_path / "CancerVar.py"
    original = "-protocol refGene,clinvar_2026\n"
    script.write_text(original)
    with pytest.raises(utils.StepFailed) as raised:
        utils.adopt_annovar_database(str(script), "clinvar_2026", "clinvar_20260102")
    assert "inside the new one" in str(raised.value)
    assert script.read_text() == original, "a refusal must not have edited the file first"


def test_adoption_refuses_a_versionless_name_rather_than_corrupting_the_script(utils, tmp_path):
    """`get_dbname_version` yields an empty version for a name with no digits, which makes the
    old name a bare prefix -- and a global replace of `dbnsfp` by `dbnsfp47a` also rewrites
    `dbnsfp31a_interpro`, and rewrites the replacement's own prefix on a second pass. No live
    configuration reaches this, and it is now a refusal rather than a silent corruption."""
    script = tmp_path / "CancerVar.py"
    script.write_text("-protocol dbnsfp,dbnsfp31a_interpro\n")
    with pytest.raises(utils.StepFailed) as raised:
        utils.adopt_annovar_database(str(script), "dbnsfp", "dbnsfp47a")
    assert "carries no version" in str(raised.value)
    assert script.read_text() == "-protocol dbnsfp,dbnsfp31a_interpro\n"


@pytest.mark.skipif(
    shutil.which("sed") is None, reason="no sed on this host, so there is nothing to compare with"
)
def test_the_shell_out_this_replaced_really_does_nothing_on_this_host(tmp_path):
    """The anti-vacuity control for the whole `sed` half of this ticket.

    `subprocess.run(f"sed -i 's/old/new/g' {file}", shell=True)` was the previous implementation,
    and it is GNU-only: BSD `sed` reads the expression as `-i`'s backup suffix and the path as the
    script, so it exits non-zero having changed nothing. Nothing inspected that, so on macOS the
    annotation script kept naming the previous database while the log reported an update.

    This test asserts what `sed -i` does *here*, whichever host that is, so the claim is measured
    on the machine reading it rather than taken on trust. On GNU it documents that the old code
    worked; on BSD it demonstrates the defect the Python substitution removes.
    """
    target = tmp_path / "script.py"
    original = "-protocol avsnp151\n"
    target.write_text(original)
    completed = subprocess.run(
        f"sed -i 's/avsnp151/avsnp152/g' {target}",
        shell=True,
        capture_output=True,
        text=True,
    )
    if completed.returncode == 0:
        assert target.read_text() == "-protocol avsnp152\n", "GNU sed should have substituted"
    else:
        assert target.read_text() == original, (
            "BSD sed is expected to leave the file unchanged when it fails; if it now edits the "
            "file this test's premise needs revisiting"
        )
        assert not list(tmp_path.glob("script.py*[!y]")), "no backup file should appear"


# ---------------------------------------------------------------------------------------------
# The config line: order is load-bearing, and only adopted names may appear on it.
# ---------------------------------------------------------------------------------------------


def test_the_config_line_keeps_its_order(utils, config_file):
    """CancerVar's `database_names` lines up with the `-operation` codes in the annotation script,
    so its order is not cosmetic.

    The previous implementation rebuilt the line by partitioning into "not updated" and "updated"
    and concatenating, which moved `dbnsfp47a_interpro` past `gnomad_genome`; the manual swap of
    exactly those two entries that followed it existed to undo that. Substituting in place cannot
    reorder anything, so this asserts the property the swap was compensating for.
    """
    configs = config_file.read_text().splitlines(keepends=True)
    before = database_names(config_file)
    out = utils.rewrite_database_names(configs, {"avsnp151": "avsnp152"})
    config_file.write_text("".join(out))
    after = database_names(config_file)

    assert after == [{"avsnp151": "avsnp152"}.get(n, n) for n in before]
    assert after.index("dbnsfp47a_interpro") < after.index("gnomad_genome")
    assert len(after) == len(before)


def test_a_database_that_did_not_update_keeps_its_old_name(utils, config_file):
    configs = config_file.read_text().splitlines(keepends=True)
    out = utils.rewrite_database_names(configs, {})
    assert "".join(out) == config_file.read_text()


def test_rewriting_a_config_with_no_database_names_line_is_a_failure(utils):
    with pytest.raises(utils.StepFailed) as raised:
        utils.rewrite_database_names(["[section]\n"], {"a1": "a2"})
    assert "no `database_names` line" in str(raised.value)


# ---------------------------------------------------------------------------------------------
# multi_update: the name moves only when every build it needs is there.
# ---------------------------------------------------------------------------------------------


def _multi_update(utils, humandb, annotation_script, tmp_path, rows, name_version, **stub):
    downloader = write_stub_downloader(tmp_path / "annotate_variation.pl", **stub)
    run_log = utils.RunLog(str(tmp_path / "update.log"))
    renames = utils.multi_update(
        annovar_databases=annovar_table(rows),
        annovar_db_path=str(humandb),
        annovar_download_script=str(downloader),
        annotation_script=str(annotation_script),
        run_log=run_log,
        name_version=name_version,
    )
    return renames, run_log


def test_a_database_that_downloads_for_every_build_is_adopted(
    utils, humandb, annotation_script, tmp_path
):
    renames, run_log = _multi_update(
        utils,
        humandb,
        annotation_script,
        tmp_path,
        [("avsnp152", "hg19"), ("avsnp152", "hg38")],
        {"avsnp": "151"},
    )
    assert renames == {"avsnp151": "avsnp152"}
    assert verdicts(run_log) == [("UPDATED", "avsnp")]
    assert "avsnp151" not in annotation_script.read_text()
    for build in ("hg19", "hg38"):
        assert (humandb / f"{build}_avsnp152.txt").is_file()


def test_a_database_already_at_the_latest_version_is_current_not_updated(
    utils, humandb, annotation_script, tmp_path
):
    renames, run_log = _multi_update(
        utils,
        humandb,
        annotation_script,
        tmp_path,
        [("avsnp151", "hg19"), ("avsnp151", "hg38")],
        {"avsnp": "151"},
    )
    assert renames == {}
    assert verdicts(run_log) == [("CURRENT", "avsnp")]
    assert "avsnp151" in annotation_script.read_text()


def test_a_download_that_writes_nothing_leaves_every_name_alone(
    utils, humandb, annotation_script, tmp_path
):
    """The shape the whole contract exists for, in this route's own terms: the update did not
    happen, so nothing may claim it did and nothing may start pointing at it."""
    before = annotation_script.read_text()
    renames, run_log = _multi_update(
        utils,
        humandb,
        annotation_script,
        tmp_path,
        [("avsnp152", "hg19"), ("avsnp152", "hg38")],
        {"avsnp": "151"},
        builds=(),
        returncode=1,
    )
    assert renames == {}
    assert verdicts(run_log) == [("FAILED", "avsnp")]
    assert annotation_script.read_text() == before
    reason = run_log.entries[0][3]
    assert "was not adopted" in reason and "avsnp151" in reason


def test_a_database_that_downloads_for_only_one_build_is_not_adopted(
    utils, humandb, annotation_script, tmp_path
):
    """The check that the annotation script's build-blindness is respected.

    `-protocol` carries no build: CancerVar runs against one `--buildver` at a time off the same
    list. So advancing the name after only `hg38` landed would leave every `hg19` run pointing at
    a table that was never written -- which is the failure this route is capable of, and it is
    invisible until somebody runs the other build.
    """
    before = annotation_script.read_text()
    renames, run_log = _multi_update(
        utils,
        humandb,
        annotation_script,
        tmp_path,
        [("avsnp152", "hg19"), ("avsnp152", "hg38")],
        {"avsnp": "151"},
        builds=("hg38",),
    )
    assert renames == {}
    assert verdicts(run_log) == [("FAILED", "avsnp")]
    assert annotation_script.read_text() == before
    assert "hg19_avsnp152.txt was not created" in run_log.entries[0][3]
    # The half that did land stays on disk: it is additive, so it costs nothing and the next run
    # finds it already present.
    assert (humandb / "hg38_avsnp152.txt").is_file()


def test_a_failing_database_does_not_stop_the_next_one(
    utils, humandb, annotation_script, tmp_path
):
    """Ruling 3, and it is safe here for the same reason as everywhere else: a failed database is
    a no-op, so carrying on cannot compound anything."""
    renames, run_log = _multi_update(
        utils,
        humandb,
        annotation_script,
        tmp_path,
        [("avsnp152", "hg19"), ("avsnp152", "hg38"), ("cosmic99", "hg38")],
        {"avsnp": "151", "cosmic": "98"},
        builds=("hg38",),
    )
    assert verdicts(run_log) == [("FAILED", "avsnp"), ("UPDATED", "cosmic")]
    assert renames == {"cosmic98": "cosmic99"}


# ---------------------------------------------------------------------------------------------
# ClinVar: both builds, and the failure the previous check could not see.
# ---------------------------------------------------------------------------------------------


def write_stub_clinvar_script(path: Path, *, builds=("hg19", "hg38"), returncode=0) -> Path:
    """Stand in for `update_clinvar_annovar.sh`, whose real run needs `vt` and `gatk`."""
    body = ["#!/bin/bash", "output=''", "name=''"]
    body += [
        'while [ $# -gt 0 ]; do',
        '  case "$1" in',
        '    --output) output="$2"; shift 2;;',
        '    --name) name="$2"; shift 2;;',
        '    *) shift;;',
        "  esac",
        "done",
    ]
    for build in builds:
        body += [
            f'printf \'%b\' "'
            + CLINVAR_HEADER.replace("\t", "\\t").replace("\n", "\\n")
            + CLINVAR_RECORD.replace("\t", "\\t").replace("\n", "\\n")
            + f'" > "$output/{build}_$name.txt"',
            f'size=$(wc -c < "$output/{build}_$name.txt" | tr -d " ")',
            f'printf \'#BIN\\t1000\\t%s\\n\' "$size" > "$output/{build}_$name.txt.idx"',
        ]
    body.append(f"exit {returncode}")
    path.write_text("\n".join(body) + "\n")
    path.chmod(path.stat().st_mode | stat.S_IEXEC)
    return path


def _clinvar(utils, humandb, annotation_script, tmp_path, current_databases, **stub):
    scriptdir = tmp_path / "scriptdir"
    scriptdir.mkdir(exist_ok=True)
    write_stub_clinvar_script(scriptdir / "update_clinvar_annovar.sh", **stub)
    vt = tmp_path / "vt"
    vt.mkdir(exist_ok=True)
    run_log = utils.RunLog(str(tmp_path / "update.log"))
    renames = utils.update_clinvar_annovar(
        annovar_db_path=str(humandb),
        annotation_script=str(annotation_script),
        current_databases=current_databases,
        run_log=run_log,
        scriptdir=str(scriptdir),
        today="20260820",
        vt=str(vt),
    )
    return renames, run_log


def test_clinvar_is_adopted_when_both_builds_are_written(
    utils, humandb, annotation_script, tmp_path
):
    renames, run_log = _clinvar(
        utils, humandb, annotation_script, tmp_path, ["avsnp151", "clinvar_20260102"]
    )
    assert renames == {"clinvar_20260102": "clinvar_20260820"}
    assert verdicts(run_log) == [("UPDATED", "clinvar")]
    assert "clinvar_20260102" not in annotation_script.read_text()


def test_clinvar_is_not_adopted_when_only_hg19_is_written(
    utils, humandb, annotation_script, tmp_path
):
    """The regression test for the check this ticket replaced.

    The previous code ran the rebuild and then asked one question --
    `os.path.exists(f"{humandb}/hg19_clinvar_{today}.txt")` -- so an `hg38` half that never
    happened was invisible, although the script writes both and although CancerVar is run against
    both builds. Every name moved anyway, and the log read `SUCCESS: clinvar updated to <today>`.
    """
    before = annotation_script.read_text()
    with pytest.raises(utils.StepFailed) as raised:
        _clinvar(
            utils,
            humandb,
            annotation_script,
            tmp_path,
            ["clinvar_20260102"],
            builds=("hg19",),
        )
    assert "hg38_clinvar_20260820.txt was not created" in str(raised.value)
    assert annotation_script.read_text() == before, "no name may move on a failed rebuild"


def test_clinvar_already_built_today_is_current(utils, humandb, annotation_script, tmp_path):
    """A same-day re-run neither rebuilds nor claims an update, and says which."""
    for build in ("hg19", "hg38"):
        write_annovar_db(humandb, build, "clinvar_20260820")
    renames, run_log = _clinvar(
        utils, humandb, annotation_script, tmp_path, ["clinvar_20260820"], builds=()
    )
    assert renames == {}
    assert verdicts(run_log) == [("CURRENT", "clinvar")]


def test_a_bad_leftover_from_today_does_not_block_the_rebuild(
    utils, humandb, annotation_script, tmp_path
):
    """The same-day short-circuit tests whether today's tables PASS, not whether they exist.

    Otherwise a run that left a header-only table behind would make the rebuild look unnecessary
    for the rest of the day, and every later run would fail on the leftover without ever trying
    again.
    """
    for build in ("hg19", "hg38"):
        write_annovar_db(humandb, build, "clinvar_20260820", body=CLINVAR_HEADER)
    renames, run_log = _clinvar(
        utils, humandb, annotation_script, tmp_path, ["clinvar_20260102"]
    )
    assert renames == {"clinvar_20260102": "clinvar_20260820"}
    assert verdicts(run_log) == [("UPDATED", "clinvar")]


def test_a_config_naming_no_clinvar_is_a_failure(utils, humandb, annotation_script, tmp_path):
    with pytest.raises(utils.StepFailed) as raised:
        _clinvar(utils, humandb, annotation_script, tmp_path, ["avsnp151"])
    assert "names no clinvar database" in str(raised.value)


# ---------------------------------------------------------------------------------------------
# The route end to end: the script and the config always agree.
# ---------------------------------------------------------------------------------------------


def _route(utils, humandb, annotation_script, config_file, tmp_path, rows, **stub):
    scriptdir = tmp_path / "scriptdir"
    scriptdir.mkdir(exist_ok=True)
    write_stub_clinvar_script(scriptdir / "update_clinvar_annovar.sh", **stub)
    vt = tmp_path / "vt"
    vt.mkdir(exist_ok=True)
    downloader = write_stub_downloader(tmp_path / "annotate_variation.pl")
    run_log = utils.RunLog(str(tmp_path / "update.log"))
    exit_code = utils.run_steps(
        [
            (
                "cancervar",
                lambda: utils.update_all_cancervar(
                    annovar_databases=annovar_table(rows),
                    annovar_db_path=str(humandb),
                    annovar_download_script=str(downloader),
                    cancervar_config_file=str(config_file),
                    cancervar_script=str(annotation_script),
                    email=None,
                    run_log=run_log,
                    not_to_update=NOT_TO_UPDATE,
                    password=None,
                    scriptdir=str(scriptdir),
                    today="20260820",
                    vt=str(vt),
                ),
            )
        ],
        run_log,
        ["hg19", "hg38"],
    )
    return exit_code, run_log


ROWS_ALL_NEW = [
    ("avsnp152", "hg19"),
    ("avsnp152", "hg38"),
    ("dbnsfp48a", "hg19"),
    ("dbnsfp48a", "hg38"),
    ("dbnsfp48a_interpro", "hg19"),
    ("dbnsfp48a_interpro", "hg38"),
    ("icgc29", "hg19"),
    ("icgc29", "hg38"),
]


def test_a_whole_route_adopts_everything_it_validated(
    utils, humandb, annotation_script, config_file, tmp_path
):
    exit_code, run_log = _route(
        utils, humandb, annotation_script, config_file, tmp_path, ROWS_ALL_NEW
    )
    names = database_names(config_file)
    assert exit_code == 0
    assert "clinvar_20260820" in names
    assert "avsnp152" in names and "avsnp151" not in names
    # COSMIC is a deliberate step-level skip and must not colour the exit code (ruling 5).
    assert ("SKIPPED", "cosmic") in verdicts(run_log)
    assert "cosmic98" in names


def test_a_clinvar_failure_still_leaves_the_script_and_the_config_agreeing(
    utils, humandb, annotation_script, config_file, tmp_path
):
    """The defect the ticket names, and the reason the config write is not left to the driver.

    The annotation script is rewritten per database as the loop goes, while the config is written
    once at the end. Previously a failed ClinVar raised out of the middle of that -- and nothing
    in `update_annovar.py` caught anything -- so the script came out naming the new `avsnp` while
    the config still named the old one, and Intervar was never updated at all. Now the failure is
    a verdict, the config write still happens, and the only thing ClinVar's failure costs is
    ClinVar.
    """
    exit_code, run_log = _route(
        utils,
        humandb,
        annotation_script,
        config_file,
        tmp_path,
        ROWS_ALL_NEW,
        builds=("hg19",),
    )
    names = database_names(config_file)
    script = annotation_script.read_text()

    assert exit_code == 1, "a failed step must make the run exit non-zero"
    assert ("FAILED", "clinvar") in verdicts(run_log)
    # ClinVar did not move, anywhere.
    assert "clinvar_20260102" in names
    assert "clinvar_20260820" not in names
    assert "clinvar_20260820" not in script
    # Everything that did validate moved in both places, and they agree.
    assert "avsnp152" in names and "avsnp152" in script
    assert "avsnp151" not in names and "avsnp151" not in script


def test_one_database_that_cannot_be_adopted_does_not_cost_the_others(
    utils, humandb, annotation_script, config_file, tmp_path
):
    """Adoption is the one per-database step that writes to a file every other database shares.

    So a failure there is the one that can leave the annotation script rewritten for the databases
    that came earlier and the config rewritten for none of them -- the disagreement this ticket
    exists to remove, arrived at from inside the loop rather than from ClinVar. This drives it by
    removing one database's name from the annotation script, so its adoption has nothing to
    rewrite while its neighbours' adoptions succeed.

    The mutation harness is what asked for this test: with the per-database `try` removed, every
    other assertion in this file still passed.
    """
    annotation_script.write_text(annotation_script.read_text().replace("icgc28", "icgc99"))
    exit_code, run_log = _route(
        utils, humandb, annotation_script, config_file, tmp_path, ROWS_ALL_NEW
    )
    names = database_names(config_file)
    script = annotation_script.read_text()

    assert exit_code == 1
    assert ("FAILED", "icgc") in verdicts(run_log)
    # The one that could not be adopted kept its name on the config line...
    assert "icgc28" in names and "icgc29" not in names
    # ...and everything else still moved, in both files, which means the config was written.
    assert "avsnp152" in names and "avsnp152" in script
    assert "clinvar_20260820" in names and "clinvar_20260820" in script


def test_every_name_the_config_carries_is_named_by_the_annotation_script(
    utils, humandb, annotation_script, config_file, tmp_path
):
    """The property that makes the two files' agreement checkable rather than asserted per case.

    Only the versioned names are compared: the annotation script's `-protocol` list is a subset of
    `database_names` -- `check_downdb` pre-fetches from the config, the annotation runs off the
    protocol -- so the direction that matters is that the script never names a version the config
    does not.
    """
    _route(utils, humandb, annotation_script, config_file, tmp_path, ROWS_ALL_NEW, builds=("hg19",))
    names = set(database_names(config_file))
    script = annotation_script.read_text()
    for versioned in ("avsnp", "clinvar_", "cosmic"):
        named_by_script = {
            token
            for token in script.replace(",", " ").replace('"', " ").split()
            if token.startswith(versioned)
        }
        assert named_by_script <= names, (
            f"{named_by_script - names} is named by the annotation script but not by the config, "
            "so annotation would use a database nothing pre-fetches"
        )


# ---------------------------------------------------------------------------------------------
# `SUCCESS` is retired -- checked here per route, and tree-wide in `test_failure_contract.py`. The
# tree-wide one could not ship until both routes were converted, since either half on its own left
# 21 occurrences in the other; #354 being the second of the two, it lands there. These stay because
# they name the route, so a failure says which half regressed.
# ---------------------------------------------------------------------------------------------

ANNOVAR_FUNCTIONS = (
    "get_annovar_databases",
    "get_dbname_version",
    "multi_update",
    "update_cosmic_annovar",
    "read_database_names",
    "update_clinvar_annovar",
    "_update_annovar_route",
    "update_all_cancervar",
    "update_all_intervar",
)


def _annovar_sources():
    """The source of each Annovar-route function, by name, from the module's own AST."""
    path = REPO_ROOT / "update_db" / "scripts" / "utils_update.py"
    text = path.read_text()
    tree = ast.parse(text)
    found = {
        node.name: ast.get_source_segment(text, node)
        for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name in ANNOVAR_FUNCTIONS
    }
    missing = set(ANNOVAR_FUNCTIONS) - set(found)
    assert not missing, (
        f"{sorted(missing)} no longer exist in utils_update.py, so this check has stopped "
        "covering the route it names. Rename them here or drop them deliberately."
    )
    return found


def test_no_annovar_function_writes_a_success_verdict():
    """`SUCCESS` was doing four different jobs across this tool -- "updated", "already the
    latest", "no update needed" and, at the ClinVar step, "did not happen at all". Since the
    destination is that no step reports success for work that did not happen, the word is where
    the ambiguity lives and it is retired rather than qualified (ruling 4).

    Read off the string literals rather than the raw text, for the same reason as the `sed` check
    below: prose explaining that `SUCCESS` was retired is not an occurrence of it.
    """
    offenders = {}
    for name, source in _annovar_sources().items():
        tree = ast.parse(source)
        docstrings = {
            id(node.body[0].value)
            for node in ast.walk(tree)
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef, ast.Module))
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
                offenders.setdefault(name, []).append(node.value[:60])
    assert not offenders, f"these Annovar functions still write SUCCESS: {offenders}"


def test_the_annovar_entry_point_reports_an_exit_code():
    """`update_annovar.py` had no `try`/`except` anywhere and never set an exit code, so a caller
    checking `$?` saw 0 whether or not half the databases updated."""
    source = (REPO_ROOT / "update_db" / "scripts" / "update_annovar.py").read_text()
    assert "sys.exit(main())" in source
    assert "run_steps(" in source


def test_every_annovar_shell_out_names_bash():
    """The route-scoped half of the interpreter guard next door, and it can be green today.

    `test_shell_interpreters.py` asserts the tree-wide property and is a strict xfail until both
    route tickets have landed, so while it is marked it cannot notice a regression on this route.
    This one can: every literal interpreter in an Annovar function must be `bash`. It matters most
    for `update_clinvar_annovar.sh`, which opens with `function usage {` -- a parse error in dash,
    the container's `/bin/sh` -- so under the previous invocation the whole script exited 2 having
    run nothing, and its own `exit 1` after every stage never got the chance to fire.
    """
    posix = {"sh", "/bin/sh"}
    offenders = {}
    for name, source in _annovar_sources().items():
        for node in ast.walk(ast.parse(source)):
            if not isinstance(node, ast.Call) or not node.args:
                continue
            if getattr(node.func, "attr", getattr(node.func, "id", None)) != "run":
                continue
            argv = node.args[0]
            head = None
            if isinstance(argv, (ast.List, ast.Tuple)) and argv.elts:
                first = argv.elts[0]
                if isinstance(first, ast.Constant):
                    head = first.value
            elif isinstance(argv, ast.Constant) and isinstance(argv.value, str):
                head = argv.value.split(None, 1)[0]
            elif isinstance(argv, ast.JoinedStr) and argv.values:
                first = argv.values[0]
                if isinstance(first, ast.Constant) and isinstance(first.value, str):
                    head = first.value.split(None, 1)[0] if first.value.split() else None
            if head in posix:
                offenders.setdefault(name, []).append(head)
    assert not offenders, (
        f"these Annovar shell-outs hand a bash script to a POSIX shell: {offenders}"
    )


def test_no_annovar_function_shells_out_to_sed():
    """The four GNU-only `sed -i` calls are gone, and none may come back: they are the silent
    no-op on macOS, and a substitution nothing inspects cannot support a verdict.

    Read off the AST rather than the text, because the text is not evidence: the comments that
    explain why these calls went name the command they replaced, and a substring scan cannot tell
    an explanation from a call. That distinction is the whole reason the interpreter guard next
    door walks the AST too.
    """
    offenders = {}
    for name, source in _annovar_sources().items():
        for node in ast.walk(ast.parse(source)):
            if not isinstance(node, ast.Call) or not node.args:
                continue
            func = node.func
            if getattr(func, "attr", getattr(func, "id", None)) != "run":
                continue
            command = ast.unparse(node.args[0])
            if "sed" in command:
                offenders[name] = command
    assert not offenders, f"these Annovar functions shell out to sed again: {offenders}"


# ---------------------------------------------------------------------------------------------
# The shell script's own half of the contract.
# ---------------------------------------------------------------------------------------------

CLINVAR_SCRIPT = REPO_ROOT / "update_db" / "annovar4cancervar_intervar" / "update_clinvar_annovar.sh"


def test_the_clinvar_rebuild_script_stops_at_the_first_failure():
    """`pipefail` is the load-bearing half and `set -e` alone would not do: both `date=`
    assignments take their value from a pipeline, and a pipeline's status is its last command's,
    so a failed `zcat` is invisible without it and the script goes on to build `clinvar_hg38_.vcf`.

    Asserted as a statement rather than as a substring, because `"set -eo pipefail" in text` is
    satisfied by a commented-out `# set -eo pipefail` -- which is how the first version of this
    test read, and the mutation harness is what found it.
    """
    lines = [line.strip() for line in CLINVAR_SCRIPT.read_text().splitlines()]
    assert "set -eo pipefail" in lines, (
        "the script must set it as a statement; a commented-out line satisfies a substring check "
        "while changing nothing"
    )
    assert any(line.startswith("date=$(zcat") for line in lines), (
        "the pipeline this setting exists for has moved; re-derive whether pipefail still covers it"
    )


def test_the_clinvar_rebuild_script_never_refuses_with_a_success_status():
    """Three refusals -- too few arguments, no `vt` directory, no output directory -- were a bare
    `exit`, which carries the status of the last command run. That was an `echo`, so the script
    reported success for having declined to do anything."""
    lines = CLINVAR_SCRIPT.read_text().splitlines()
    bare = [
        number
        for number, line in enumerate(lines, start=1)
        if line.strip() in ("exit", "exit;")
    ]
    # `help` legitimately exits 0: it was asked to print the usage and it did. Every other bare
    # `exit` in this script was a refusal reporting success.
    assert len(bare) == 1, (
        "unexpected bare `exit` statements, which carry the previous command's status, on lines "
        f"{bare}"
    )
    inside_help = [line.strip() for line in lines[: bare[0]]]
    assert "function help {" in inside_help, "the only bare exit should be the one inside `help`"
