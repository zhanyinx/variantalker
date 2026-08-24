# The step functions near the top of this module are annotated with `RunLog`, which the failure
# contract defines at the BOTTOM -- deliberately, so that the line numbers the tickets working on
# this file measured stay usable. Without this, those annotations would be evaluated at definition
# time and the module would not import.
from __future__ import annotations

import argparse
from ftplib import FTP
import base64
import glob
import gzip
import datetime
import io
import logging
import os
import re
import subprocess
import tempfile
import traceback
from typing import Union


# import gsutil
from bs4 import BeautifulSoup
import numpy as np
import pandas as pd
import requests

from pandas.testing import assert_frame_equal


logging.basicConfig(
    format="%(asctime)s | %(levelname)s: %(message)s", level=logging.NOTSET
)


class FolderType:
    """Custom type supporting folders."""

    def __init__(self):
        pass

    def __call__(self, value):  # noqa: D102
        if not os.path.isdir(value):
            raise argparse.ArgumentTypeError(
                f"Input value must be folder and must exist. '{value}' is not."
            )
        return value


class FileType:
    """Custom type supporting folders."""

    def __init__(self):
        pass

    def __call__(self, value):  # noqa: D102
        if not os.path.isfile(value):
            raise argparse.ArgumentTypeError(
                f"Input value must be folder and must exist. '{value}' is not."
            )
        return value


def get_version(file: str):
    """
    file: configuration file used to extract database version.
    """
    with open(file) as f:
        for line in f:
            if line.startswith("version = "):
                v = line.split("=")[1]
                return v.strip()


ONCOTATOR_CURRENT = "oncotator_v1_ds_April052016.tar.gz"
DNA_REPAIR_GENES_URL = "https://www.mdanderson.org/research/departments-labs-institutes/labs/wood-laboratory/resources.html"
DNA_REPAIR_GENES_CURRENT = (
    "Table information was last updated by R. Wood and B. Dennehey, November 2025."
)
OREGANNO_CURRENT = "ORegAnno_Combined_2016.01.19.tsv"


def check_oncotator(run_log: RunLog):
    """Check whether the manually-maintained oncotator bundle has moved on. Hard coded.

    A check, not an update: it touches no database, so its verdict is `CURRENT` when the upstream
    bundle still matches the installed one and `SKIPPED` when a human has to act. Neither affects
    the exit code, because neither is a failure of this run.

    The two `os._exit(1)` calls this used to make are retired (#352 ruling 3). They bypassed the
    step loop and the `RESULT:` summary, and -- because `os._exit` skips stdio flushing and they
    sat inside the `with open(...)` block -- they could lose the very `WARNING` line they had just
    written, killing the process at exit 1 with a log that never said why.
    """
    data = []
    with FTP("ftp.broadinstitute.org", user="gsapubftp-anonymous") as ftp:
        ftp.cwd("bundle/oncotator")
        ftp.dir(data.append)
    name = "achilles, cancer-gene-census, familial & simple-uniprot"

    if not len(data):
        raise StepFailed(
            f"oncotator: the Broad bundle listing came back empty, so nothing was compared "
            f"against the installed {ONCOTATOR_CURRENT}"
        )
    if len(data) > 1:
        run_log.skipped(
            "oncotator",
            f"{len(data)} entries in the Broad bundle listing where one was expected - a new "
            f"database may be available and needs a manual update. Installed: "
            f"{ONCOTATOR_CURRENT}",
            detail=[f"{name}"] + [line for line in data],
        )
        return

    search = re.search(r"^.+\s+(\d+)\s+(oncotator_.+)$", data[0])
    if search is None:
        raise StepFailed(
            f"oncotator: could not read a filename out of the Broad bundle listing entry "
            f"{data[0]!r}, so the installed {ONCOTATOR_CURRENT} was not compared against anything"
        )
    size, latest = search.groups()
    if latest != ONCOTATOR_CURRENT:
        run_log.skipped(
            "oncotator",
            f"a newer database is available and needs a manual update: {latest!r} ({size} bytes). "
            f"Installed: {ONCOTATOR_CURRENT}",
            detail=[name],
        )
    else:
        run_log.current("oncotator", None, ONCOTATOR_CURRENT, reason="no manual update needed")


def check_dna_repair_genes(run_log: RunLog):
    """Check whether the DNA repair gene table has been revised upstream. (hard coded)

    The first live step of a Funcotator run, and the reason a TOTAL network outage never used to
    reach the destructive ClinVar step: this raised an uncaught `ConnectionError` and took the
    process down first. Under the contract it becomes a `FAILED` verdict and the run continues --
    which is safe precisely because every later step now validates its stage before touching
    anything.
    """
    r = requests.get(DNA_REPAIR_GENES_URL)
    soup = BeautifulSoup(r.text, features="lxml")
    stamp = None
    for line in soup.get_text().splitlines():
        if "Table information was last updated" in line:
            stamp = line.strip()

    # Was an unbound local: when the page stopped carrying that line -- a redesign, a captive
    # portal, an error page served with HTTP 200 -- this raised `UnboundLocalError` from the
    # comparison below, which said nothing about the page.
    if stamp is None:
        raise StepFailed(
            "dna_repair_genes: the Wood laboratory resources page carries no 'Table information "
            f"was last updated' line, so nothing was compared against the installed table "
            f"({DNA_REPAIR_GENES_URL})"
        )

    if stamp == DNA_REPAIR_GENES_CURRENT:
        run_log.current(
            "dna_repair_genes",
            None,
            "November 2025",
            reason="the source table has not been revised since the installed copy",
        )
    else:
        run_log.skipped(
            "dna_repair_genes",
            f"the source table has been revised and needs a manual update, at "
            f"{DNA_REPAIR_GENES_URL}",
            detail=[f"upstream now says: {stamp}"],
        )


def check_oreganno(run_log: RunLog):
    """Check whether a newer ORegAnno dump has been published. (hard coded)"""
    df = (
        pd.read_html("http://www.oreganno.org/dump/")[0][
            ["Name", "Last modified", "Size"]
        ]
        .dropna()
        .reset_index(drop=True)
    )
    if not len(df):
        raise StepFailed(
            "oreganno: the dump listing at http://www.oreganno.org/dump/ held no named files, so "
            f"nothing was compared against the installed {OREGANNO_CURRENT}"
        )

    latest = sorted(df["Name"], reverse=True)[0]
    if latest == OREGANNO_CURRENT:
        run_log.current("oreganno", None, OREGANNO_CURRENT, reason="no manual update needed")
    else:
        run_log.skipped(
            "oreganno",
            f"a newer dump is available and needs a manual update: {latest!r}. Installed: "
            f"{OREGANNO_CURRENT}",
        )


def update_cosmic(
    db_dir: str,
    curr_version: str,
    backup_dir: str,
    run_log: RunLog,
    scriptdir: str,
    email: str,
    password: str,
    build: str = "both",
):
    """Update cosmic, cosmic_tissue, and cosmic_fusion if a newer version is available on the Sanger Institute website.

    NOTE: This function is currently DISABLED due to format changes in COSMIC database.
    COSMIC has changed its download format and authentication method, requiring script updates.

    Re-enabling it is a separate effort, and the code below the early return has NOT been brought
    under the failure contract: its destructive block still runs unconditionally, on a `cosmic*`
    glob rather than a per-build directory. Whoever re-enables it must convert it first -- there is
    no `STAGE_SPEC` entry for `cosmic`, so `_require_valid_stage` would refuse it by design.

    Args:
        db_dir (str): The directory path where the Cosmic database is stored.
        curr_version (str): The current version of the Cosmic database.
        backup_dir (str): The directory path where the backup of the current Cosmic database will be stored.
        run_log (RunLog): The log every verdict is written to.
        scriptdir (str): The directory path where the Cosmic update script is stored.
        email (str): email to log in into cancer.sanger.ac.uk
        password (str): password to login into cancer.sanger.ac.uk
        build (str): Genome build to update: 'hg19', 'hg38', or 'both' (default: 'both').
    """
    # DISABLED: COSMIC update temporarily disabled due to format changes
    run_log.skipped(
        "cosmic",
        f"update disabled - COSMIC changed its download format and authentication method, so a "
        f"manual update is required. Installed: {curr_version}",
    )
    return

    # Original code commented out - update when COSMIC format is fixed
    # # Check if credentials are provided
    # if not email or not password:
    #     with open(file, "a") as f:
    #         f.write(
    #             f"SKIPPED: Cosmic update skipped - no credentials provided. "
    #             f"Please provide --cosmic_email and --cosmic_password to update COSMIC database. "
    #             f"Current version: {curr_version}\n"
    #         )
    #     return
    
    passcode = base64.b64encode(f"{email}:{password}".encode("ascii")).decode()

    # Find latest cosmic version
    r = requests.get("https://cancer.sanger.ac.uk/cosmic/")
    soup = BeautifulSoup(r.text)
    string = soup.find("section", attrs={"id": "index-intro"}).find("h1").text.strip()
    version, date = re.search(r"^COSMIC\s(v\d+),\sreleased\s([\w-]+)$", string).groups()

    if curr_version < version:
        # TODO change path
        update_script = f"{scriptdir}/update_cosmic/updateCosmicDataBase.sh"

        # Run update script. `bash`, not a POSIX shell, like every other call site on this route --
        # `updateCosmicDataBase.sh` is bash-shebanged. Unreachable today, and fixed in passing so
        # that re-enabling COSMIC does not also re-enable the interpreter bug.
        subprocess.run(
            f"bash {update_script} {passcode} {version}",
            shell=True,
            stdout=subprocess.DEVNULL,
        )
        subprocess.run(
            f"cp -r {db_dir}/cosmic* {backup_dir}",
            shell=True,
            stdout=subprocess.DEVNULL,
        )
        subprocess.run(
            f"rm -rf {db_dir}/cosmic*",
            shell=True,
            stdout=subprocess.DEVNULL,
        )
        subprocess.run(f"cp -r cosmic* {db_dir}", shell=True, stdout=subprocess.DEVNULL)
        subprocess.run(f"rm -rf cosmic*", shell=True, stdout=subprocess.DEVNULL)
        run_log.updated("cosmic", None, curr_version, version)
    else:
        run_log.current("cosmic", None, curr_version, reason=f"{version} is the latest published")


def update_dbsnp(
    run_log: RunLog,
    db_dir: str,
    backup_dir: str,
    scriptdir: str,
    db_germline_dir: str = None,
    build: str = "both",
):
    """Update the dbSNP database used in the Funcotator annotation, per requested build.

    This was the one step in the tool that inspected anything at all, and it is the model the rest
    of the route now follows: it runs its script under `bash`, compares md5sums rather than
    assuming, reads the child's machine-readable `SKIPPED:` lines, and gates its destructive block
    on a check. What it did NOT do was distinguish a build it was not asked for from a build that
    is fine -- its `files_created` gate is an `or`, asking whether ANYTHING was staged rather than
    whether every live build was accounted for, so `--build hg38` replaced the whole `dbsnp`
    directory in both trees and the live hg19 copy was gone, logged `SUCCESS: dbSNP updated`.

    Three verdicts now come out of three different inspections:

    * the script printed `SKIPPED: <build> dbSNP` -- it refuses to initialise a database that is
      not installed. Per-build, benign, and the reason is carried through.
    * the script staged nothing for a requested build and exited 0 -- its md5sum comparison found
      the live copy already current. `CURRENT`, because a comparison genuinely ran.
    * the script staged nothing and exited non-zero -- `FAILED`, carrying the exit status. This is
      the case `set -eo pipefail` in the script exists to make visible: without `pipefail` a failed
      download still produced a config saying `version = b`, a version-shaped string with no
      version in it.

    Args:
        run_log (RunLog): The log every verdict is written to.
        db_dir (str): Path to the directory containing funcotator databases.
        backup_dir (str): Path to the directory where a backup copy of the current dbsnp files will be stored.
        scriptdir (str): Path to the directory where the external shell script is located.
        db_germline_dir (str, optional): Path to the directory containing the germline dbsnp files. Default is None.
        build (str): Genome build to update: 'hg19', 'hg38', or 'both' (default: 'both').
    """
    builds = requested_builds(build)

    # Convert relative paths to absolute to ensure correct file access
    db_dir = os.path.abspath(db_dir)
    backup_dir = os.path.abspath(backup_dir)
    scriptdir = os.path.abspath(scriptdir)
    if db_germline_dir:
        db_germline_dir = os.path.abspath(db_germline_dir)

    update_script = f"{scriptdir}/update_dbsnp/update_dbsnp.sh"

    clear_stage("dbsnp")
    result = run_update_script(update_script, ["-d", f"{db_dir}/dbsnp/", "--build", build])
    output = result.stdout + result.stderr
    detail = [
        line for line in output.splitlines() if "SKIPPED:" in line or "To initialize" in line
    ]

    for one in builds:
        if f"SKIPPED: {one} dbSNP" in output:
            run_log.skipped(
                "dbsnp",
                "no existing database found; initialise it manually before it can be updated",
                build=one,
                detail=[line for line in detail if one in line],
            )

    attempted = [b for b in builds if f"SKIPPED: {b} dbSNP" not in output]
    if not attempted:
        # Every requested build was deliberately not attempted. `_settle_step` turns that into a
        # step-level FAILED, which is exactly what the old `raise ValueError` achieved for the
        # `--build both` case -- without aborting the run, so `acmg_rec` still gets to run.
        return

    versions = _require_valid_stage("dbsnp", staged_builds("dbsnp", attempted))

    for one in attempted:
        installed = live_version("dbsnp", one, db_dir)
        if one not in versions:
            if result.returncode != 0:
                raise StepFailed(
                    f"dbsnp {one}: the update script exited {result.returncode} without staging a "
                    "database, so the live one is untouched"
                )
            run_log.current(
                "dbsnp",
                one,
                installed or "version not recorded",
                reason="the published md5sum matches the installed database",
            )
            continue

        _require_germline_in_sync("dbsnp", [one], db_dir, db_germline_dir)
        swap_build("dbsnp", one, db_dir, backup_dir, db_germline_dir=db_germline_dir)
        run_log.updated("dbsnp", one, installed or "not installed", versions[one])

    clear_stage("dbsnp")


def update_gencode(
    db_dir: str,
    curr_version: str,
    backup_dir: str,
    run_log: RunLog,
    scriptdir: str,
    db_germline_dir: str = None,
    build: str = "both",
):
    """
    Updates the gencode database used in the somatic Funcotator annotation.

    NOTE: this step is DISABLED at the call site -- the rebuilt database looks correct but
    Funcotator does not use it properly. It is brought under the failure contract here anyway,
    because the alternative is leaving a step that replaces the live GENCODE database with whatever
    the script happened to leave behind, waiting to be re-enabled.

    There is deliberately no `STAGE_SPEC` entry for `gencode`: its stage is written by the vendored
    `getGencode.sh`, whose output shape has not been measured, and a spec written from guesswork
    would be a check that looks like it is working. By contract an unspecified step is a hard
    error, so re-enabling this step fails loudly with a message saying what to add, rather than
    quietly installing an unvalidated database.

    Args:
        db_dir (str): Path to the directory containing funcotator databases.
        curr_version (str): Version number of the current gencode release.
        backup_dir (str): Path to the directory where a backup copy of the current gencode files will be stored.
        run_log (RunLog): The log every verdict is written to.
        scriptdir (str): Path to the directory where the external shell script is located.
        db_germline_dir (str, optional): Path to the directory containing the germline gencode files. Default is None.
        build (str): Genome build to update: 'hg19', 'hg38', or 'both' (default: 'both').
    """
    builds = requested_builds(build)

    # Skip update if curr_version is None (database doesn't exist for requested build)
    if curr_version is None:
        run_log.skipped("gencode", f"no existing database found for build {build}")
        return

    data = []
    with FTP("ftp.ebi.ac.uk") as ftp:
        ftp.login()
        ftp.cwd("/pub/databases/gencode/Gencode_human/latest_release")
        ftp.dir(data.append)

    # Find latest version
    annotations = []
    for line in data:
        search = re.search(r"gencode\.v(\d+)\.", line)
        if search is not None:
            version = search.groups()[0]
            break

    if curr_version < version:
        # TODO change path, copy also the python script
        update_script = f"{scriptdir}/update_gencode/getGencode.sh"
        # Read current getgencode and replace
        with open(update_script, "r") as f:
            data = f.readlines()

        for idx, line in enumerate(data):
            if line.startswith("LATEST_RELEASE"):
                data[idx] = f"LATEST_RELEASE={version}\n"
        with open(update_script, "w") as f:
            f.writelines(data)

        # `getGencode.sh` is the one vendored script here, and the only one that already set
        # `set -e` -- it is still bash-shebanged, so a POSIX shell was the wrong interpreter for it
        # too. Good error handling behind the wrong interpreter is worth nothing.
        versions, result = stage_and_validate(
            "gencode", update_script, ["--build", build], builds, run_log
        )

        missing = [b for b in builds if b not in versions]
        if missing:
            raise StepFailed(
                f"gencode: the update script staged nothing for {', '.join(missing)} (it exited "
                f"{result.returncode}), so the live database is untouched"
            )

        _require_germline_in_sync("gencode", builds, db_dir, db_germline_dir)
        for one in builds:
            swap_build("gencode", one, db_dir, backup_dir, db_germline_dir=db_germline_dir)
            run_log.updated("gencode", one, curr_version, version)
        clear_stage("gencode")
    else:
        run_log.current(
            "gencode", None, curr_version, reason=f"v{version} is the latest published release"
        )


def update_clinvar(
    db_dir: str,
    backup_dir: str,
    run_log: RunLog,
    scriptdir: str,
    db_germline_dir: str = None,
    build: str = "both",
):
    """Update ClinVar, or leave both live databases exactly as they were.

    This is the step the whole failure contract was written for. As it stood, it had one code path
    and one log line, so it was structurally incapable of reporting anything but success -- and
    what it actually did when the download failed was replace the somatic AND germline ClinVar
    databases with an empty directory and then write `SUCCESS: clinvar updated successfully.` at
    exit 0. Measured three ways (as shipped, with the interpreter fixed, and with the interpreter
    fixed plus `set -e`), all three destroyed both databases, because the destructive sequence ran
    unconditionally before anything looked at what the update script had produced.

    Now: stage, inspect, and only then swap -- per requested build, so the build that was not asked
    for is never named. A failure at any point leaves every live database untouched and says which
    check refused it.

    Args:
        db_dir (str): Directory path for the somatic database.
        backup_dir (str): Directory path for the backup.
        run_log (RunLog): The log every verdict is written to.
        scriptdir (str): Directory path to the script.
        db_germline_dir (str, optional): Directory path for the germline database. Defaults to None.
        build (str): Genome build to update: 'hg19', 'hg38', or 'both' (default: 'both').
    """
    builds = requested_builds(build)

    # TODO remove hard code script
    update_script = f"{scriptdir}/update_clinvar/update_clinvar_funcotator.sh"

    # The interpreter is fixed for the whole route in `run_update_script`, and it is necessary but
    # nowhere near sufficient: with it fixed and the download still failing, the stage carries a
    # 0-byte `clinvar_.vcf` and a blank-version config, which Funcotator loads as a valid data
    # source with zero records. Measured. The validation below is what makes the difference, which
    # is why the two changes had to land together rather than as two commits.
    versions, result = stage_and_validate(
        "clinvar", update_script, ["--build", build], builds, run_log
    )

    missing = [b for b in builds if b not in versions]
    if missing:
        raise StepFailed(
            f"clinvar: the update script staged nothing for {', '.join(missing)} (it exited "
            f"{result.returncode}), so there is no new database to install and the live one is "
            "untouched"
        )

    # Germline receives the same staged directory somatic does, so it is derivable from somatic
    # rather than independent -- and the only backup taken below is somatic's. Prove they were in
    # sync before overwriting, or the somatic backup is not a germline restore path (#352 ruling 6).
    _require_germline_in_sync("clinvar", builds, db_dir, db_germline_dir)

    for one in builds:
        installed = live_version("clinvar", one, db_dir)
        if installed == versions[one]:
            run_log.current("clinvar", one, versions[one])
            continue
        swap_build("clinvar", one, db_dir, backup_dir, db_germline_dir=db_germline_dir)
        run_log.updated("clinvar", one, installed or "not installed", versions[one])

    clear_stage("clinvar")


def update_hgnc(
    db_dir: str, backup_dir: str, run_log: RunLog, scriptdir: str, build: str = "both"
):
    """Update HGNC, per requested build, comparing each build against its own live copy.

    Two defects met in this step, and fixing either one alone would have entrenched the other.

    The **hardcoded `hg38`** in both globs meant `--build hg19` had never been able to work: the
    update script stages `hgnc/hg19` only, so the staged glob returned an empty list and the step
    died on a bare `IndexError`. Accidentally safe -- nothing was destroyed -- but permanently
    broken, and reported as *"There was a problem with the curl download"* when no download had
    failed. Turning those `IndexError`s into a well-worded failure WITHOUT parameterising the build
    would have documented the bug rather than fixed it.

    The **whole-directory replace** meant a successful `--build hg38` run deleted the live hg19
    copy. So the comparison moves inside the per-build loop: with the two builds now allowed to sit
    at different versions (they are byte-identical by construction today, but `--build` is
    per-build and coordinate-free data honours it uniformly), a single hg38 comparison cannot
    decide hg19. A current hg38 beside a stale hg19 now produces `CURRENT` for one and `UPDATED`
    for the other.

    Args:
        db_dir (str): Directory path for the database.
        backup_dir (str): Directory path for the backup.
        run_log (RunLog): The log every verdict is written to.
        scriptdir (str): Directory path to the script.
        build (str): Genome build to update: 'hg19', 'hg38', or 'both' (default: 'both').
    """
    builds = requested_builds(build)

    # TODO update path
    update_script = f"{scriptdir}/update_hgnc/get_new_hgnc.sh"

    versions, result = stage_and_validate(
        "hgnc", update_script, ["--build", build], builds, run_log
    )

    for one in builds:
        if one not in versions:
            raise StepFailed(
                f"hgnc: the update script staged nothing for {one} (it exited "
                f"{result.returncode}), so the live database is untouched"
            )

        staged = glob.glob(os.path.join("hgnc", one, "hgnc_*.tsv"))[0]
        latest_data = pd.read_csv(staged, sep="\t")

        # If a whole column is NA, some column retrieval in the update script went wrong and the
        # curl link needs updating. Was a message-less `assert`, which the generic handler could
        # only report as `unexpected AssertionError:` with nothing after the colon -- honest, but
        # barely more useful than the wrong message it replaced.
        if not len(latest_data.dropna()):
            raise StepFailed(
                f"hgnc {one}: staged {staged} has {len(latest_data)} rows but none complete, so a "
                "whole column came back empty - check the curl link in "
                f"{update_script}"
            )

        current = glob.glob(os.path.join(db_dir, "hgnc", one, "hgnc_*.tsv"))
        if not current:
            # The tool's rule everywhere else: "no existing database found" is a skip, not a
            # failure. Per-build it is benign, so a host with only hg38 installed can still run
            # `--build both` without exiting 1 every night -- but if EVERY requested build is
            # skipped the step compared nothing at all, and `_settle_step` escalates that to
            # FAILED without this function having to know how many builds were asked for.
            run_log.skipped(
                "hgnc",
                f"no existing database at {db_dir}/hgnc/{one}/hgnc_*.tsv to compare against; "
                "install it manually first",
                build=one,
            )
            continue

        current_data = pd.read_csv(current[0], sep="\t")
        try:
            assert_frame_equal(current_data, latest_data, check_names=False)
            same = True
        except AssertionError:
            same = False

        installed = live_version("hgnc", one, db_dir)
        if same:
            run_log.current(
                "hgnc",
                one,
                installed or versions[one],
                reason="the source table did not change since the installed copy",
            )
            continue

        swap_build("hgnc", one, db_dir, backup_dir)
        run_log.updated("hgnc", one, installed or "not installed", versions[one])

    clear_stage("hgnc")


def update_acmg_rec(
    run_log: RunLog,
    db_germline_dir: str,
    backup_dir: str,
    current_version: str,
    build: str = "both",
):
    """Update the acmg_rec Funcotator database, per requested build.

    The one step that stages without a shell script -- it scrapes the ClinVar ACMG page and writes
    the table itself -- and the only one that is germline-only. That makes it the one database
    germline keeps a backup of, because it is derivable from nothing: for `clinvar`, `dbsnp` and
    `gencode` the germline copy receives the same staged bytes somatic does, so somatic's backup
    restores both, but nothing else in the tree can reconstruct this one.

    Coordinate-free, like HGNC: the data is a bare `Disease_Name\\tgene` table and its `hg19` and
    `hg38` copies are byte-identical by construction, existing only because Funcotator requires a
    per-build directory layout. It still honours `--build` per build rather than refreshing both
    copies, deliberately: the alternative needs a maintained list of which databases are
    coordinate-free, and a per-step classification table is exactly the shape that goes vacuous by
    falling behind the code. The accepted cost is that the two copies may sit at different
    versions, which is correct behaviour for a flag that means "leave the other build alone".

    Args:
        run_log (RunLog): The log every verdict is written to.
        db_germline_dir (str): Directory path for the germline database.
        backup_dir (str): Directory path for the backup.
        current_version (str): Current version of the database.
        build (str): Genome build to update: 'hg19', 'hg38', or 'both' (default: 'both').
    """
    builds = requested_builds(build)

    # Skip update if current_version is None (database doesn't exist for requested build)
    if current_version is None:
        run_log.skipped("acmg_rec", f"no existing database found for build {build}")
        return

    today = datetime.date.today().strftime("%b%d%Y")

    r = requests.get("https://www.ncbi.nlm.nih.gov/clinvar/docs/acmg/")
    soup = BeautifulSoup(r.text)

    # get version
    found = re.search(r"ACMG SF (v\d\.\d)", str(soup.findAll("p")))
    if found is None:
        raise StepFailed(
            "acmg_rec: the ClinVar ACMG page carries no 'ACMG SF vN.N' version string, so there "
            f"was nothing to compare the installed {current_version} against"
        )
    version = found[1]

    if version > current_version:
        # get the new data
        table = soup.find("tbody")
        header, *table_rows = table.findChildren("tr")

        disease_names = []
        genes = []

        for row in table_rows:
            columns = row.findChildren("td")
            if len(columns):
                if (
                    columns[0].text.strip() != "MedGen"
                    and columns[1].find("a").text != "ClinVar"
                ):
                    disease_name = columns[0].text.strip()

                try:
                    gene = columns[2].find("a").text
                except IndexError:
                    gene = columns[0].find("a").text

                if gene == "ClinVar":
                    gene = columns[1].find("a").text

                disease_names.append(disease_name)
                genes.append(gene)

        # write to file
        with open(f"acmg_{version}_{today}_test_cleaned.txt", "w") as f:
            f.write("Disease_Name\tgene\n")

            for row in zip(disease_names, genes):
                f.write("\t".join(row) + "\n")

        # write config file, taking the template from a build that is actually installed rather
        # than from a hardcoded hg38 -- on a `--build hg19` host there may be no hg38 copy to read.
        template = None
        for candidate in builds + ["hg38", "hg19"]:
            path = f"{db_germline_dir}/acmg_rec/{candidate}/acmg_rec.config"
            if os.path.isfile(path):
                template = path
                break
        if template is None:
            raise StepFailed(
                f"acmg_rec: no installed acmg_rec.config under {db_germline_dir}/acmg_rec/ to take "
                "the config template from, so the staged database would have no config"
            )

        with open(template, "r") as f:
            data = f.readlines()

        for idx, line in enumerate(data):
            if line.startswith("version"):
                data[idx] = f"version = {version}\n"

            if line.startswith("src_file"):
                data[idx] = f"src_file = acmg_{version}_{today}_test_cleaned.txt\n"
            if line.startswith("preprocessing_script"):
                data[idx] = f"preprocessing_script = update_funcotator.py\n"

        with open("acmg_rec.config", "w") as f:
            f.writelines(data)

        # Stage every requested build, and only the requested ones. The two copies are identical
        # by construction, so the second is a copy of the first rather than a second scrape.
        table = f"acmg_{version}_{today}_test_cleaned.txt"
        clear_stage("acmg_rec")
        for one in builds:
            _run_checked(
                ["mkdir", "-p", os.path.join("acmg_rec", one)],
                f"creating the staged acmg_rec/{one}",
            )
            for name in (table, "acmg_rec.config"):
                _run_checked(
                    ["cp", name, os.path.join("acmg_rec", one, name)],
                    f"staging {name} for {one}",
                )
        for name in (table, "acmg_rec.config"):
            _run_checked(["rm", "-f", name], f"cleaning up {name}")

        # Inspect the stage before replacing anything, exactly as the shell-script steps do -- this
        # step writes its own stage, but a scrape that returned an empty table is no more
        # installable than a download that returned an empty VCF.
        versions = _require_valid_stage("acmg_rec", builds)

        for one in builds:
            installed = live_version("acmg_rec", one, db_germline_dir)
            swap_build("acmg_rec", one, db_germline_dir, backup_dir)
            run_log.updated("acmg_rec", one, installed or current_version, versions[one])

        clear_stage("acmg_rec")
    else:
        run_log.current(
            "acmg_rec", None, current_version, reason=f"{version} is the latest published"
        )


# Cancervar Intervar and annovar


def get_markdown_table(header, data):
    rows = []
    start_idx = data.index(header) + 3  # Title, Header, ---
    for line in data[start_idx:]:
        if not line.startswith("|"):
            break
        rows.append([c.strip() for c in line.split("|") if c.strip()])
    return rows


def get_annovar_databases(url: str):
    """Get annovar databases."""
    # Get md data from annovar docs GitHub
    r = requests.get(url)
    data = [l for l in r.text.split("\n") if l]

    # Parse data from markdown
    rows = get_markdown_table("### - For gene-based annotation", data)
    rows.extend(get_markdown_table("### - For filter-based annotation", data))

    df = pd.DataFrame(rows, columns=["build", "name", "explanation", "date"])
    return df


def get_dbname_version(current_dbs: list):
    """Get current db version and name."""
    result = dict()
    for db in current_dbs:
        if re.search("^(.*?)(\\d.*)", db):
            search = re.search("^(.*?)(\\d.*)", db)
            if "interpro" in search[2]:
                result[search[1] + "_interpro"] = search[2]
            else:
                result[search[1]] = search[2]
        else:
            result[db] = ""
    return result


def multi_update(
    annovar_databases: pd.DataFrame,
    annovar_db_path: FolderType(),
    annovar_download_script: FileType(),
    annotation_script: FileType(),
    run_log: "RunLog",
    name_version: dict,
) -> dict:
    """Update the annovar databases that can be fetched from annovar's own download server.

    One verdict per database, and it is about the name the annotation actually uses rather than
    about a file appearing on disk -- because on this route those are different things. Nothing
    is adopted until every build annovar publishes for the new version has a table that passed
    `require_valid_annovar_db`: the annotation script's `-protocol` list carries no build, so
    advancing the name after only one of two builds landed would point an `hg19` run at a file
    that was never written.

    A database whose update fails is left exactly as it was -- old name in the script, old name
    in the config, old file still on disk and still what gets annotated with -- and the run
    carries on to the next one (ruling 3).

    Args:
        annovar_databases: A pandas DataFrame with columns name, build, explanation and date,
            containing the list of available databases.
        annovar_db_path: Path to the folder where annovar databases are stored.
        annovar_download_script: Path to the script used to download annovar databases.
        annotation_script: Path to the annotation script that needs to be updated.
        run_log: The run's log; every outcome is reported through it.
        name_version: A dictionary mapping database names to their versions.

    Returns:
        dict: current versioned name -> newly adopted name, for the databases that updated. A
            database that is already current, or that failed, is absent -- which is what keeps
            `rewrite_database_names` from ever naming a database no file stands behind.
    """
    renames = {}
    # iterate over databases
    for label, version in name_version.items():
        name = label
        # Capture interpro flag
        interpro = False
        if "_interpro" in name:
            name = name.replace("_interpro", "")
            interpro = True

        select = annovar_databases[
            [name in x for x in annovar_databases["name"].values]
        ].sort_values("name", ascending=False)

        if not interpro:
            select = select[~select["name"].str.contains("interpro")]

        # get most recent database on annovar webpage
        most_recent = select["name"].values[0]
        current = "".join([name, version])

        # update if most recent is more recent
        if most_recent <= current:
            run_log.current(label, None, current)
            continue

        select = select[select.name == most_recent]

        # iterate on build
        notes = []
        problems = []
        for index, row in select.iterrows():
            build = row["build"]
            database = row["name"]
            logging.info(f"Updating {name} to {most_recent} for build {build}")
            outfile, _ = annovar_db_paths(annovar_db_path, build, database)

            # run database creation only if database is missing
            downloaded = False
            if not os.path.exists(outfile):
                completed = subprocess.run(
                    [
                        annovar_download_script,
                        "-buildver",
                        build,
                        "-downdb",
                        "-webfrom",
                        "annovar",
                        database,
                        annovar_db_path,
                    ],
                    capture_output=True,
                    text=True,
                )
                downloaded = True
                if completed.returncode:
                    # Kept as detail rather than as the verdict: the downloader's exit status is
                    # the least trustworthy thing about it -- nothing in this file checked it
                    # before, and the inventory found no step that could -- so what decides the
                    # outcome is the table it was supposed to leave behind.
                    notes.append(
                        f"{build}: the downloader exited {completed.returncode}"
                    )
                    notes.extend((completed.stderr or "").strip().splitlines()[-3:])

            try:
                index_note = require_valid_annovar_db(annovar_db_path, build, database)
            except StepFailed as exc:
                problems.append(str(exc))
                continue

            notes.append(
                f"{build}: {os.path.basename(outfile)} "
                + ("downloaded" if downloaded else "was already present")
            )
            if index_note:
                notes.append(f"{build}: {index_note}")

        if problems:
            run_log.failed(
                label,
                f"{most_recent} was not adopted, so annotation continues with {current}. "
                + "; ".join(problems),
                detail=notes,
            )
            continue

        # Both names move together, and only now: the annotation script here, and the config's
        # `database_names` from the returned mapping once every database has been through.
        #
        # Caught per database rather than left to escape. Adoption is the one step here that
        # writes to a file shared by every database in the loop, so letting a failure out would
        # abandon the run with the annotation script already rewritten for the databases that
        # came earlier and the config not yet rewritten for any of them -- which is precisely the
        # disagreement between the two files that this ticket exists to remove.
        try:
            adopt_annovar_database(annotation_script, current, most_recent)
        except StepFailed as exc:
            run_log.failed(label, str(exc), detail=notes)
            continue

        run_log.updated(label, None, current, most_recent, detail=notes)
        renames[current] = most_recent

    return renames


def update_cosmic_annovar(
    annovar_db_path: FolderType(),
    annotation_script: FileType(),
    curr_version: str,
    email: str,
    run_log: "RunLog",
    password: str,
    scriptdir: FolderType(),
) -> str:
    """Update cosmic, cosmic_tissue, and cosmic_fusion if a newer version is available on the Sanger Institute website.

    NOTE: This function is currently DISABLED due to format changes in COSMIC database.
    COSMIC has changed its download format and authentication method, requiring script updates.

    Re-enabling it is out of scope for this map, but it was moved onto the failure contract along
    with everything else on this route -- verdicts instead of `SUCCESS`, `bash` instead of `sh`,
    and a checked substitution instead of a GNU-only `sed -i` -- so that switching it back on does
    not also switch the old defects back on. Its live caller reports the skip itself; what is
    below the return is the disabled body.

    Args:
        annovar_db_path: Output folder where to save databases files. In general humandb.
        curr_version: The current version of the Cosmic database.
        email: email to log in into cancer.sanger.ac.uk
        run_log: The run's log; every outcome is reported through it.
        password: password to login into cancer.sanger.ac.uk
        scriptdir: The directory path where the Cosmic update script is stored.
    """
    # DISABLED: COSMIC update temporarily disabled due to format changes
    run_log.skipped(
        "cosmic",
        "update disabled - the COSMIC download format and authentication changed, so this "
        "needs a manual update",
        detail=[f"still using cosmic{curr_version}"],
    )
    return curr_version

    # Original code commented out - update when COSMIC format is fixed
    # # Check if credentials are provided
    # if not email or not password:
    #     run_log.skipped(
    #         "cosmic",
    #         "no credentials provided - pass --cosmic_email and --cosmic_password to update it",
    #         detail=[f"still using cosmic{curr_version}"],
    #     )
    #     return curr_version

    passcode = base64.b64encode(f"{email}:{password}".encode("ascii")).decode()

    # Find latest cosmic version
    r = requests.get("https://cancer.sanger.ac.uk/cosmic/")
    soup = BeautifulSoup(r.text)
    string = soup.find("section", attrs={"id": "index-intro"}).find("h1").text.strip()
    version, date = re.search(r"^COSMIC\sv(\d+),\sreleased\s([\w-]+)$", string).groups()
    if curr_version < version and not os.path.exists(
        f"{annovar_db_path}/hg38_cosmic{version}.txt"
    ):
        logging.info("Updating cosmic to {version}")
        # TODO change path
        update_script = f"{scriptdir}/update_cosmic.sh"

        # Run update script. `bash`, not `sh`: every script under update_db/ is a bash script,
        # and the container's /bin/sh is dash, where this one dies on a syntax error.
        subprocess.run(
            ["bash", update_script, version, passcode, annovar_db_path],
            stdout=subprocess.DEVNULL,
        )
        require_valid_annovar_db(annovar_db_path, "hg38", f"cosmic{version}")
        # Update annotation script. Was a GNU-only `sed -i`, a silent no-op on macOS; see
        # `adopt_annovar_database`. Converted along with the rest of them even though this sits
        # behind the disabled return above, so that re-enabling COSMIC does not also re-introduce
        # the defect -- and moved below the validation, so the name cannot advance to a database
        # that was not written.
        adopt_annovar_database(
            annotation_script, f"cosmic{curr_version}", f"cosmic{version}"
        )
        run_log.updated("cosmic", None, f"cosmic{curr_version}", f"cosmic{version}")
    else:
        run_log.current("cosmic", None, f"cosmic{version}")
    return version


def read_database_names(config_file: str):
    """The databases a CancerVar/Intervar config currently names, in the order it names them.

    Read before anything is touched, so that a config this tool cannot understand stops the step
    while the tree is still consistent rather than after half of it has moved.

    Raises:
        StepFailed: if there is no `database_names` line to read.
    """
    with open(config_file, errors="replace") as f:
        configs = f.readlines()
    for line in configs:
        if line.startswith("database_names"):
            return configs, line.split("=", 1)[1].split()
    raise StepFailed(
        f"{config_file} has no `database_names` line, so there is nothing to read the current "
        "database versions from"
    )


def update_clinvar_annovar(
    annovar_db_path: FolderType(),
    annotation_script: FileType(),
    current_databases: list,
    run_log: "RunLog",
    scriptdir: FolderType(),
    today: str,
    vt: FolderType(),
) -> dict:
    """Rebuild ClinVar for annovar from the current NCBI release, and adopt it if it is good.

    Shared by both routes, which ran a verbatim copy of this each. Three changes from that copy,
    all of them the contract:

    * the script runs under **`bash`**. It opens with `function usage {`, which is a parse error
      in dash -- the container's `/bin/sh` -- so under the previous invocation the whole thing
      exited 2 having run nothing, and its own careful `exit 1` after every stage never ran;
    * **both builds are checked**, and checked for content rather than for existence. The
      previous check looked only for `hg19_clinvar_<today>.txt`, so an `hg38` half of the build
      that never happened was invisible even though the script writes both;
    * failure raises `StepFailed` **before** either name moves, so a failed rebuild leaves the
      previous ClinVar in use rather than pointing the annotation at a file that is not there.

    Returns:
        dict: `{current name: adopted name}`, empty if today's build was already in use.
    """
    target = f"clinvar_{today}"
    current = next((x for x in current_databases if x.startswith("clinvar")), None)
    if current is None:
        raise StepFailed(
            "the config names no clinvar database, so there is no current version to update "
            "from and nothing that would start using a new one"
        )

    builds = ("hg38", "hg19")

    def staged_is_good():
        """Whether today's tables are already there and already pass, without rebuilding."""
        try:
            for build in builds:
                require_valid_annovar_db(annovar_db_path, build, target)
        except StepFailed:
            return False
        return True

    # Validate before deciding to rebuild, not merely test for the files' presence: a same-day
    # re-run should skip the rebuild when today's tables are good, but a previous run that left
    # something behind that does NOT pass must not be able to wedge the step for the rest of the
    # day by making the rebuild look unnecessary.
    completed = None
    if not staged_is_good():
        logging.info(f"Updating clinvar to {today}")
        completed = subprocess.run(
            [
                "bash",
                f"{scriptdir}/update_clinvar_annovar.sh",
                "-v",
                vt,
                "--name",
                target,
                "--output",
                annovar_db_path,
            ],
            capture_output=True,
            text=True,
        )

    notes = []
    if completed is not None and completed.returncode:
        notes.append(f"the rebuild script exited {completed.returncode}")
        notes.extend((completed.stderr or "").strip().splitlines()[-3:])

    problems = []
    for build in builds:
        try:
            index_note = require_valid_annovar_db(annovar_db_path, build, target)
        except StepFailed as exc:
            problems.append(str(exc))
            continue
        if index_note:
            notes.append(f"{build}: {index_note}")
    if problems:
        raise StepFailed(
            f"{target} was not adopted, so annotation continues with {current}. "
            + "; ".join(problems + notes)
        )

    if current == target:
        run_log.current("clinvar", None, target, reason="already built today")
        return {}

    adopt_annovar_database(annotation_script, current, target)
    run_log.updated("clinvar", None, current, target, detail=notes)
    return {current: target}


def _update_annovar_route(
    annovar_databases: pd.DataFrame,
    annovar_db_path: FolderType(),
    annovar_download_script: FileType(),
    config_file: str,
    annotation_script: FileType(),
    run_log: "RunLog",
    not_to_update: list,
    scriptdir: FolderType(),
    today: str,
    vt: FolderType(),
    cosmic: bool,
):
    """One annovar route -- CancerVar or Intervar -- from its config to its config.

    The two were near-identical copies, differing only in whether they carry a COSMIC database.

    The ordering here is the point of the ticket, so it is worth stating: the config's
    `database_names` line is rewritten **once, at the end, from the adoptions that actually
    happened**, and a ClinVar failure is caught rather than allowed to escape. Previously the
    annotation script was rewritten per database as the loop went while the config was written
    only at the very end, so any failure in between -- and a failed ClinVar raised -- left the
    script naming one version and the config another, with `update_all_intervar` never running at
    all because nothing in `update_annovar.py` caught it.
    """
    configs, current_databases = read_database_names(config_file)

    # the ones that can be updated simply by downloading them from annovar's server
    to_update = [x for x in current_databases if x not in not_to_update]
    to_update = [x for x in to_update if "clinvar" not in x and "cosmic" not in x]

    # update avsnp, dbnsfp and icgc
    renames = multi_update(
        annovar_databases=annovar_databases,
        annovar_db_path=annovar_db_path,
        annovar_download_script=annovar_download_script,
        annotation_script=annotation_script,
        run_log=run_log,
        name_version=get_dbname_version(to_update),
    )

    # COSMIC update - DISABLED due to format changes in the download and authentication. A
    # step-level skip: deliberate, benign, and it must not affect the exit code (ruling 5).
    # Re-enabling it is a separate effort; when it happens it calls `update_cosmic_annovar` and
    # adds its result to `renames` like any other database.
    if cosmic:
        present = [x for x in current_databases if x.startswith("cosmic")]
        run_log.skipped(
            "cosmic",
            "update disabled - the COSMIC download format and authentication changed",
            detail=[f"still using {present[0]}" if present else "none configured"],
        )

    # update clinvar. Caught here rather than left to the driver: the adoptions above have
    # already been written into the annotation script, so the config write below is what makes
    # the two agree, and it has to happen whether or not ClinVar worked.
    try:
        renames.update(
            update_clinvar_annovar(
                annovar_db_path=annovar_db_path,
                annotation_script=annotation_script,
                current_databases=current_databases,
                run_log=run_log,
                scriptdir=scriptdir,
                today=today,
                vt=vt,
            )
        )
    except StepFailed as exc:
        run_log.failed("clinvar", str(exc))

    # Record the adoptions, in place, so the order CancerVar depends on cannot move.
    with open(config_file, "w") as f:
        f.writelines(rewrite_database_names(configs, renames))


def update_all_cancervar(
    annovar_databases: pd.DataFrame,
    annovar_db_path: FolderType(),
    annovar_download_script: FileType(),
    cancervar_config_file: str,
    cancervar_script: FileType(),
    email: str,
    run_log: "RunLog",
    not_to_update: list,
    password: str,
    scriptdir: FolderType(),
    today: str,
    vt: FolderType(),
):
    """Update all CancerVar-related annovar databases.

    `email` and `password` are the COSMIC credentials, unused while that update is disabled and
    kept so that re-enabling it is a change to one function rather than to the entry point too.
    """
    _update_annovar_route(
        annovar_databases=annovar_databases,
        annovar_db_path=annovar_db_path,
        annovar_download_script=annovar_download_script,
        config_file=cancervar_config_file,
        annotation_script=cancervar_script,
        run_log=run_log,
        not_to_update=not_to_update,
        scriptdir=scriptdir,
        today=today,
        vt=vt,
        cosmic=True,
    )


def update_all_intervar(
    annovar_databases: pd.DataFrame,
    annovar_db_path: FolderType(),
    annovar_download_script: FileType(),
    intervar_config_file: FileType(),
    intervar_script: FileType(),
    run_log: "RunLog",
    not_to_update: list,
    scriptdir: FolderType(),
    today: str,
    vt: FolderType(),
):
    """Update all Intervar-related annovar databases. Intervar carries no COSMIC database."""
    _update_annovar_route(
        annovar_databases=annovar_databases,
        annovar_db_path=annovar_db_path,
        annovar_download_script=annovar_download_script,
        config_file=intervar_config_file,
        annotation_script=intervar_script,
        run_log=run_log,
        not_to_update=not_to_update,
        scriptdir=scriptdir,
        today=today,
        vt=vt,
        cosmic=False,
    )


# =============================================================================================
# The failure contract
# =============================================================================================
#
# The shared machinery every update step uses to decide -- and then to say -- whether it changed
# the live annotation database. Landed by issue #360 ahead of the two route tickets (#353 for the
# Funcotator route, #354 for the Annovar route) so that neither of them invents it and the other
# inherits a fait accompli. The rulings it implements were decided in #352.
#
# Nothing below is wired into a step yet, deliberately: this block may not change a single
# verdict a real run writes. Converting the steps is #353 and #354.
#
# It is appended at the end of the module rather than inserted at the top so that the line
# numbers those two tickets measured against this file stay very nearly valid -- only the two
# added imports shift them, by two lines.
#
# The load-bearing rule the whole block exists to serve: **every verdict must be derived from an
# inspection, not from reaching a line of code.** `update_clinvar` today has one code path and
# one log line, so it is structurally incapable of reporting anything but success. Nothing here
# lets a step say what happened without having looked.


class StepFailed(Exception):
    """One update step failed in a way that was anticipated, and said why.

    Raised by the checks below and caught by `run_steps`, which turns it into a `FAILED:`
    verdict and carries on to the next database. The message is the operator's only account of
    what went wrong, so it must name both what was checked and what was found -- not just that
    something was wrong. A `StepFailed` with a vague message is the failure mode this exception
    exists to replace: `update_hgnc`'s message-less `assert` and `update_funcotator`'s bare
    `except:` (which reports a curl problem when no download was ever reached) are the two
    in-tree examples.
    """


# ---------------------------------------------------------------------------------------------
# The staged-database specification.
# ---------------------------------------------------------------------------------------------
#
# Every destructive step already builds its new database into a CWD-relative directory named
# after the database and then copies that over the live tree. So the missing piece is not a
# staging area -- it is the inspection between the stage and the swap. This table says what a
# staged database has to look like for each step.
#
# Keys, all of them derived by reading the shell script that writes the stage:
#
#   data      regex matched IN FULL against the data file's name inside `<stage>/<build>/`.
#             A `version` named group is optional: two of the three steps encode the version in
#             the filename, dbSNP does not.
#   version   regex a RESOLVED version must match in full, wherever it is read from. This is
#             what makes check 2 mean something: `version = b` (dbSNP's build id never
#             resolved) and `version =` (ClinVar's date never resolved) are both blank versions
#             wearing a plausible shape, and a permissive pattern would wave the first through.
#   unresolved  regex naming the MEASURED failure shape -- the filename the stage carries when
#             the version never resolved. Optional, and only for the diagnostic message; the
#             check itself is `data` not matching.
#   index     suffix appended to the data filename to name its index, or None for a step that
#             has no index.
#   config    Funcotator config filename inside the build directory, or None.
#   header    how to recognise a header line when looking for a real data record: "hash" for
#             VCF-style `#` comment lines, "first_line" for a single column-header row.
#
# ANTI-VACUITY CONTROL, and it is part of the contract rather than a nicety: a step with no
# entry here is a HARD ERROR, not a silent skip. A per-step table is exactly the shape that goes
# vacuous by falling behind the code -- it stops applying to the step it was written for and
# nothing says so. `_require_valid_stage` therefore refuses to validate a step it has no spec
# for, and `check_stage_spec` refuses a spec with an unrecognised key, so a typo cannot silently
# switch a check off.

_STAGE_SPEC_KEYS = {"data", "version", "unresolved", "index", "config", "header"}

STAGE_SPEC = {
    # update_clinvar_funcotator.sh: clinvar/<build>/clinvar_<YYYYMMDD>.vcf, its gatk .idx, and
    # clinvar_vcf.config with DATE substituted. When the download fails, `$date` is empty and
    # the stage carries a 0-byte `clinvar_.vcf` -- the shape this whole map exists for.
    "clinvar": dict(
        data=r"clinvar_(?P<version>\d{8})\.vcf",
        version=r"\d{8}",
        unresolved=r"clinvar_\.vcf",
        index=".idx",
        config="clinvar_vcf.config",
        header="hash",
    ),
    # get_new_hgnc.sh: hgnc/<build>/hgnc_<%b%d%Y>.tsv plus hgnc.config. No index -- it is a
    # simpleXSV data source, not a VCF, so there is nothing for gatk to index.
    "hgnc": dict(
        data=r"hgnc_(?P<version>[A-Za-z]{3}\d{2}\d{4})\.tsv",
        version=r"[A-Za-z]{3}\d{2}\d{4}",
        unresolved=None,
        index=None,
        config="hgnc.config",
        header="first_line",
    ),
    # update_dbsnp.sh: dbsnp/<build>/<GCF accession>.gz, its tabix .tbi, and dbSNP.config. The
    # version is the dbSNP build id (b156) and lives ONLY in the config -- the filename carries
    # the assembly accession instead -- so there is no `version` group to agree with. What
    # replaces that agreement here is `src_file`; see `_validate_build_stage`.
    "dbsnp": dict(
        data=r"GCF_[0-9.]+\.gz",
        version=r"b\d+",
        unresolved=r"\.gz",
        index=".tbi",
        config="dbSNP.config",
        header="hash",
    ),
    # update_acmg_rec, which is the one step that stages without a shell script: it scrapes the
    # ClinVar ACMG page and writes acmg_rec/<build>/acmg_<vN.N>_<%b%d%Y>_test_cleaned.txt beside a
    # config it copies from the live database and rewrites. Added by #353 -- the step validates
    # like every other, and by contract a step with no entry here is a hard error rather than a
    # silent skip, so being unspecified was not an option once it was converted.
    "acmg_rec": dict(
        data=r"acmg_(?P<version>v\d+\.\d+)_[A-Za-z]{3}\d{2}\d{4}_test_cleaned\.txt",
        version=r"v\d+\.\d+",
        unresolved=None,
        index=None,
        config="acmg_rec.config",
        header="first_line",
    ),
}


def check_stage_spec(spec_table: dict = None):
    """Refuse a `STAGE_SPEC` that has drifted out of the shape the checks understand.

    The point is a typo: `confg="clinvar_vcf.config"` would leave `config` unset, and a spec
    with no config silently skips the version-agreement check -- the one check that catches the
    macOS `sed -i` no-op. Read once by `_require_valid_stage`, and asserted directly by the
    test suite.

    Raises:
        StepFailed: if an entry has an unrecognised key, or is missing a required one.
    """
    table = STAGE_SPEC if spec_table is None else spec_table
    for name, spec in table.items():
        unknown = set(spec) - _STAGE_SPEC_KEYS
        if unknown:
            raise StepFailed(
                f"STAGE_SPEC['{name}'] has unrecognised key(s) {sorted(unknown)}; a misspelt "
                f"key silently switches its check off. Known keys: {sorted(_STAGE_SPEC_KEYS)}"
            )
        missing = _STAGE_SPEC_KEYS - set(spec)
        if missing:
            raise StepFailed(
                f"STAGE_SPEC['{name}'] is missing key(s) {sorted(missing)}; spell them "
                "explicitly, with None where the step genuinely has no index or config, so "
                "that an absent check is a decision on the record rather than an omission"
            )
        if spec["header"] not in ("hash", "first_line"):
            raise StepFailed(
                f"STAGE_SPEC['{name}']['header'] is {spec['header']!r}; expected 'hash' or "
                "'first_line'"
            )
    return table


def requested_builds(build: str):
    """The genome builds a `--build` value asks for, in the order the scripts write them.

    Args:
        build (str): 'hg19', 'hg38' or 'both'.

    Returns:
        list: the build names to validate. Note this says nothing about what happens to a build
            that was NOT requested -- that question is open, and is issue #357.
    """
    if build == "both":
        return ["hg38", "hg19"]
    return [build]


def read_config_field(path: str, field: str):
    """One `field = value` line out of a Funcotator config, or None if it is absent or blank.

    `get_version` above does this for `version` only, and stops at the first match. This is the
    same read generalised, because the checks need `src_file` as well: it is what ties a config
    to the data file beside it, and BSD `sed`'s silent no-op breaks exactly that tie.
    """
    if not os.path.isfile(path):
        return None
    with open(path, errors="replace") as f:
        for line in f:
            if line.startswith(f"{field} ="):
                value = line.split("=", 1)[1].strip()
                return value or None
    return None


def first_data_record(path: str, header: str = "hash"):
    """The first line of real data in a staged data file, or None if there is none.

    Reads only as far as the first non-header line, so the cost is the header rather than the
    file -- which is why counting records was ruled out and this was not. Transparently reads
    the bgzipped dbSNP VCF, which is a valid gzip stream.
    """
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", errors="replace") as handle:
        if header == "first_line":
            handle.readline()
        for line in handle:
            if header == "hash" and line.startswith("#"):
                continue
            if line.strip():
                return line.rstrip("\n")
    return None


def _validate_build_stage(name: str, build: str, build_dir: str, spec: dict):
    """The five structural checks of ruling 2, against one build of one staged database.

    Returns:
        tuple: `(version, problems)`. `version` is the resolved version if every check passed,
            otherwise None; `problems` is a list of human-readable strings, empty on success.
            Problems are collected rather than raised one at a time so the operator sees
            everything wrong with the stage in a single verdict.
    """
    # 1. the build subdirectory exists.
    if not os.path.isdir(build_dir):
        return None, [f"no {build} subdirectory in the staged {name} database"]

    present = sorted(os.listdir(build_dir))
    problems = []

    # 2. the data file exists, is non-empty, and its name carries a resolved version.
    matches = [f for f in present if re.fullmatch(spec["data"], f)]
    if not matches:
        detail = f"expected a file matching {spec['data']!r}, found {present or 'nothing'}"
        unresolved = [f for f in present if spec["unresolved"] and re.fullmatch(spec["unresolved"], f)]
        if unresolved:
            # Name the file the operator will find on disk, not the pattern that matched it. An
            # earlier form of this message printed the regex -- `'clinvar_\\.vcf'` -- which is the
            # sort of diagnostic that sends someone looking for a file with a backslash in it.
            detail = (
                f"the staged file is named {unresolved[0]!r}, which is what this step leaves "
                "behind when its download failed and the version never resolved"
            )
        return None, [f"{name} {build}: no usable data file ({detail})"]
    if len(matches) > 1:
        return None, [
            f"{name} {build}: {len(matches)} files match {spec['data']!r} ({matches}), so "
            "which one is the database is ambiguous"
        ]

    data_name = matches[0]
    data_path = os.path.join(build_dir, data_name)
    if os.path.getsize(data_path) == 0:
        problems.append(f"{name} {build}: staged {data_name} is 0 bytes")

    filename_version = None
    captured = re.fullmatch(spec["data"], data_name).groupdict()
    if "version" in captured:
        filename_version = captured["version"]
        if not re.fullmatch(spec["version"], filename_version or ""):
            problems.append(
                f"{name} {build}: {data_name} carries version {filename_version!r}, which is "
                f"not a resolved version (expected {spec['version']!r})"
            )

    # 3. the index for it is present.
    if spec["index"]:
        index_path = data_path + spec["index"]
        if not os.path.isfile(index_path):
            problems.append(
                f"{name} {build}: {data_name}{spec['index']} is missing, so the data file was "
                "never indexed"
            )
        elif os.path.getsize(index_path) == 0:
            problems.append(f"{name} {build}: {data_name}{spec['index']} is 0 bytes")

    # 4. the config's version is resolved, and the config agrees with the data file beside it.
    config_version = None
    if spec["config"]:
        config_path = os.path.join(build_dir, spec["config"])
        if not os.path.isfile(config_path):
            problems.append(f"{name} {build}: {spec['config']} is missing from the stage")
        else:
            config_version = read_config_field(config_path, "version")
            if config_version is None:
                problems.append(
                    f"{name} {build}: {spec['config']} has no version, so the update script's "
                    "date substitution never ran"
                )
            elif not re.fullmatch(spec["version"], config_version):
                problems.append(
                    f"{name} {build}: {spec['config']} says version {config_version!r}, which "
                    f"is not a resolved version (expected {spec['version']!r})"
                )
            elif filename_version is not None and config_version != filename_version:
                problems.append(
                    f"{name} {build}: {spec['config']} says version {config_version!r} but the "
                    f"data file is {data_name}; the config and the database disagree"
                )

            # The same agreement in the spelling that also works for dbSNP, whose version is
            # not in the filename. This is what a BSD `sed -i` no-op looks like: exit 1, file
            # unchanged, nothing checked, so the config still names `clinvar_DATE.vcf`.
            src_file = read_config_field(config_path, "src_file")
            if src_file is not None and src_file != data_name:
                problems.append(
                    f"{name} {build}: {spec['config']} names src_file {src_file!r} but the "
                    f"staged data file is {data_name}; the config points at nothing"
                )

    # 5. at least one non-header data record.
    if os.path.getsize(data_path) > 0:
        try:
            record = first_data_record(data_path, spec["header"])
        except OSError as exc:
            problems.append(f"{name} {build}: {data_name} could not be read ({exc})")
            record = None
        if record is None:
            problems.append(
                f"{name} {build}: {data_name} has a header but not one data record, so it is a "
                "database Funcotator will accept and annotate nothing from"
            )

    if problems:
        return None, problems
    return (filename_version or config_version), []


def _require_valid_stage(name: str, builds, stage_root: str = ".", spec_table: dict = None):
    """Refuse to let a staged database reach the live tree unless it really is a database.

    This is the inspection that ruling 1 puts between the stage and the swap. Call it after the
    update script has run and BEFORE anything destructive; if it returns, the stage is good for
    every build that was requested.

    Args:
        name (str): the database, and also the staged directory's name -- 'clinvar', 'hgnc',
            'dbsnp'. Must have a `STAGE_SPEC` entry.
        builds (list): the builds that were requested, from `requested_builds`.
        stage_root (str): where the update script staged its output. Defaults to the current
            directory, which is where every script in this repo writes it.
        spec_table (dict): override for the spec table, for tests.

    Returns:
        dict: build -> the resolved version found for it.

    Raises:
        StepFailed: if `name` has no spec, or any requested build fails any check. The message
            lists every problem found, not just the first.
    """
    table = check_stage_spec(spec_table)
    if name not in table:
        raise StepFailed(
            f"{name}: no STAGE_SPEC entry, so there is nothing to validate this step's staged "
            "database against. Add an entry rather than skipping the call -- an unspecified "
            "step is a hard error by contract, because a table that quietly stops covering a "
            "step is how a guard goes vacuous"
        )
    spec = table[name]

    versions = {}
    problems = []
    for build in builds:
        version, found = _validate_build_stage(
            name, build, os.path.join(stage_root, name, build), spec
        )
        problems.extend(found)
        if version is not None:
            versions[build] = version

    if problems:
        raise StepFailed(
            f"{name}: the staged database did not pass validation, so the live database is "
            "untouched. " + "; ".join(problems)
        )
    return versions


def _require_germline_in_sync(
    name: str, builds, db_dir: str, db_germline_dir: str, spec_table: dict = None
):
    """Prove the germline copy really is derivable from the somatic one before overwriting it.

    Ruling 6. `clinvar`, `dbsnp` and `gencode` each delete the germline database and copy in
    the same staged directory they just gave somatic -- so germline receives byte-identical
    content and gets no backup of its own. That makes `backup_dir/<db>`, which holds somatic's
    old version, a faithful germline backup if and only if the two were in sync beforehand.
    Nothing guaranteed that. This checks it, by reading two config files.

    Absence is not drift: if either side has no config for a build there is nothing to lose and
    nothing to compare, so it passes and says so. Only two versions that both exist and
    disagree are drift.

    Args:
        name (str): the database.
        builds (list): the builds that were requested.
        db_dir (str): the live somatic database root.
        db_germline_dir (str): the live germline database root. A falsy value means this step
            does not write germline at all, and the check is vacuous by construction.
        spec_table (dict): override for the spec table, for tests.

    Returns:
        dict: build -> a short account of what was compared.

    Raises:
        StepFailed: if any requested build's somatic and germline versions disagree.
    """
    table = check_stage_spec(spec_table)
    if not db_germline_dir:
        return {build: "no germline directory; nothing to keep in sync" for build in builds}
    if name not in table:
        raise StepFailed(
            f"{name}: no STAGE_SPEC entry, so the germline drift check does not know which "
            "config file to read"
        )
    config = table[name]["config"]
    if not config:
        raise StepFailed(
            f"{name}: STAGE_SPEC says this database has no config, so its germline copy cannot "
            "be shown to be in sync by comparing versions. It needs a different proof before "
            "its germline directory may be overwritten"
        )

    notes = {}
    drift = []
    for build in builds:
        somatic = read_config_field(os.path.join(db_dir, name, build, config), "version")
        germline = read_config_field(
            os.path.join(db_germline_dir, name, build, config), "version"
        )
        if somatic is None:
            notes[build] = "no live somatic version to compare"
        elif germline is None:
            notes[build] = f"no live germline copy; nothing to lose (somatic {somatic})"
        elif somatic != germline:
            drift.append(
                f"{build}: somatic is {somatic} but germline is {germline}"
            )
        else:
            notes[build] = f"in sync at {somatic}"

    if drift:
        raise StepFailed(
            f"{name}: the live germline database is not derivable from the somatic one, so "
            "overwriting it would destroy content that the somatic backup cannot restore. "
            + "; ".join(drift)
        )
    return notes


# ---------------------------------------------------------------------------------------------
# The verdicts, and the log they are written to.
# ---------------------------------------------------------------------------------------------
#
# Four verdicts, defined by what happened to the LIVE database, plus exactly one RESULT line.
# `SUCCESS` is retired: across the three Python files it appeared 21 times doing four different
# jobs -- "updated", "already the latest", "no update needed", and at `update_clinvar` "did not
# happen at all". Since the destination is *no step reports SUCCESS for work that did not
# happen*, the word itself is where the ambiguity lives. There is deliberately no writer for it
# here, and there must not be one.
#
# `WARNING` and `INFO` survive only as indented detail lines beneath a verdict, which is what
# `update_dbsnp` already does when it echoes the child script's `SKIPPED:` lines under its own
# `FAILED:`.

UPDATED = "UPDATED"
CURRENT = "CURRENT"
SKIPPED = "SKIPPED"
FAILED = "FAILED"
VERDICTS = (UPDATED, CURRENT, SKIPPED, FAILED)

# Wide enough to align every verdict's text, RESULT included.
_VERDICT_WIDTH = max(len(v) for v in VERDICTS) + 1


def _qualified(db: str, build: str = None):
    """`'clinvar hg38'`, or just `'clinvar'` for a verdict about the whole step."""
    return f"{db} {build}" if build else db


class RunLog:
    """The update log, and the tally that the `RESULT:` line is derived from.

    Every verdict goes through here, which is what makes the summary honest: a run cannot be
    counted as having updated something without an `UPDATED:` line having been written, because
    the count IS the lines written. The alternative -- a step reporting its own outcome to a
    summariser -- is the shape the tool has today, where reaching a line of code is taken as
    evidence that the work behind it happened.
    """

    def __init__(self, path: str):
        self.path = path
        self.entries = []

    def _write(self, verdict: str, text: str, detail=()):
        with open(self.path, "a") as f:
            f.write(f"{verdict + ':':<{_VERDICT_WIDTH}} {text}\n")
            for line in detail:
                if str(line).strip():
                    f.write(f"  detail: {str(line).rstrip()}\n")

    def updated(self, db: str, build: str, old: str, new: str, detail=()):
        """The live database changed, after its stage passed validation."""
        self.entries.append((UPDATED, db, build, f"{old} -> {new}"))
        self._write(UPDATED, f"{_qualified(db, build)} {old} -> {new}", detail)

    def current(self, db: str, build: str, version: str, reason="already the latest", detail=()):
        """A version comparison ran, found nothing newer, and left the live database alone."""
        self.entries.append((CURRENT, db, build, reason))
        self._write(CURRENT, f"{_qualified(db, build)} ({version}) - {reason}", detail)

    def skipped(self, db: str, reason: str, build: str = None, detail=()):
        """Deliberately not attempted. A per-build skip is benign; see `settle`."""
        self.entries.append((SKIPPED, db, build, reason))
        self._write(SKIPPED, f"{_qualified(db, build)} - {reason}", detail)

    def failed(self, db: str, reason: str, build: str = None, detail=()):
        """Attempted and did not succeed. The live database is untouched."""
        self.entries.append((FAILED, db, build, reason))
        self._write(FAILED, f"{_qualified(db, build)} - {reason}", detail)

    def tally(self):
        """verdict -> how many lines of it were written."""
        return {v: sum(1 for e in self.entries if e[0] == v) for v in VERDICTS}

    def result(self):
        """Write the single `RESULT:` line and return the exit code the run should carry.

        Returns:
            int: 1 if any step failed, else 0. The caller spells `sys.exit(...)`, so that the
                driver stays callable in a test without taking the process down.
        """
        counts = self.tally()
        summary = ", ".join(f"{counts[v]} {v.lower()}" for v in VERDICTS)
        skips = [
            f"{_qualified(db, build)} {reason}"
            for verdict, db, build, reason in self.entries
            if verdict == SKIPPED
        ]
        if skips:
            summary += f" (skipped: {'; '.join(skips)})"
        self._write("RESULT", summary)
        return 1 if counts[FAILED] else 0


def run_steps(steps, run_log: RunLog, builds):
    """Run every update step, let each one fail on its own, and summarise once at the end.

    Ruling 3: a failed step does not abort the run. Ruling 1 is what makes that safe -- with the
    stage validated, a failed step is a no-op on the live database, so carrying on costs nothing
    and the operator gets every verdict instead of one plus an abort. Today's tool does the
    opposite in three of its four failure idioms: an uncaught `ValueError` from `update_dbsnp`
    means `acmg_rec` never runs, a `FileExistsError` on the Annovar route means
    `update_all_intervar` never runs (`update_annovar.py` has no `try`/`except` at all), and
    `check_oncotator`'s `os._exit(1)` dies on the spot.

    Args:
        steps (list): `(name, callable)` pairs. The callable takes no arguments -- bind the
            step's arguments at the call site -- and reports its own verdicts to `run_log`.
        run_log (RunLog): the log every verdict is written to.
        builds (list): the builds this run requested, from `requested_builds`.

    Returns:
        int: the exit code, from `RunLog.result`.
    """
    for name, call in steps:
        before = len(run_log.entries)
        try:
            call()
        except StepFailed as exc:
            run_log.failed(name, str(exc))
            continue
        except Exception as exc:  # noqa: BLE001 - the backstop; see below
            # A genuine bug in the tool surfaces as a FAILED verdict carrying the type, the
            # message and the traceback, rather than as an opaque abort or -- worse -- a
            # swallowed reason. This is the replacement for `update_funcotator.py`'s bare
            # `except:`, which reported a curl download problem for every possible HGNC failure.
            run_log.failed(
                name,
                f"unexpected {type(exc).__name__}: {exc}",
                detail=traceback.format_exc().splitlines(),
            )
            continue
        _settle_step(name, run_log, before, builds)
    return run_log.result()


def _settle_step(name: str, run_log: RunLog, before: int, builds):
    """Hold a step that returned quietly to the contract, once it is too late to hide.

    Two rules that only a look at what the step actually wrote can enforce:

    * **Ruling 4.** A step that returned without recording any verdict has told the operator
      nothing, which is indistinguishable from a step whose work silently did not happen. That
      is the defect this map exists for, so it is a failure rather than a pass.
    * **Ruling 5.** A per-build `SKIPPED` is benign and must not affect the exit code -- a host
      whose install is deliberately partial should not exit 1 every night. But a step that
      skipped EVERY build it was asked for did nothing at all, and that is `FAILED`. This
      generalises `update_dbsnp`'s existing rule (one build missing is a warning, both missing
      raises) rather than inventing one. A step-level skip, with no build, is a deliberate
      whole-step skip -- a disabled caller, absent credentials -- and stays benign.
    """
    written = run_log.entries[before:]
    if not written:
        run_log.failed(
            name,
            "the step returned without recording a verdict, so nothing inspected whether the "
            "live database changed",
        )
        return

    per_build_skips = {e[2] for e in written if e[0] == SKIPPED and e[2]}
    only_skips = all(e[0] == SKIPPED for e in written)
    if only_skips and per_build_skips and set(builds) <= per_build_skips:
        run_log.failed(
            name,
            f"no requested build was updated - every one of {sorted(set(builds))} was skipped",
        )


# =============================================================================================
# The Funcotator route: staging and the per-build swap
# =============================================================================================
#
# Applying the contract above to the Funcotator steps (issue #353). Two rulings meet here:
#
# * **#352 ruling 1** -- nothing destructive runs until the staged database has been inspected.
#   Every step already stages into a directory named after the database in the working directory
#   and then copies that over the live tree, so what was missing was never the staging area. It
#   was the inspection between the stage and the swap.
#
# * **#357 ruling 2** -- the swap is per REQUESTED build, and the unrequested build's directory is
#   never named. `--build hg38` means "update hg38 and leave hg19 alone", settled with the person
#   who added the flag; the whole-directory `rm -rf` these steps used was simply the wrong
#   granularity for it. Because the unrequested build is never touched, nothing has to be
#   restored, and correctness does not depend on `backup_dir` existing or being current.
#
# This section is Funcotator-only by measurement, not by convenience: the Annovar route has no
# build axis to get wrong. In `humandb` the build is a filename PREFIX rather than a directory
# (`{annovar_db_path}/{build}_{database}.txt`) and every file is written individually behind its
# own existence check, so nothing there is replaced wholesale and nothing can be lost.
#
# What is deliberately NOT here: atomicity. The `rm -rf` of one live build followed by a `cp -r`
# of its replacement is still interruptible, so a crash or a full disk mid-copy can leave a
# half-written build. #352 ruling 1 accepted that for the whole directory and #357 ruling 2
# narrowed it to a single build; removing it needs a same-filesystem check plus a copy fallback,
# and is a separate ticket if it is ever wanted.


def _run_checked(argv, what: str):
    """Run a command and refuse to carry on if it failed.

    The plainest expression of what this map is about. There are 39 `subprocess.run` calls in this
    module and, before this ticket, not one of them looked at `returncode` -- so a failed `cp` of a
    new database over a live one was indistinguishable from a successful one, and the log said
    `SUCCESS` either way.

    Takes an argv list rather than a command string on purpose: no quoting, no shell, and nothing
    for the interpreter guard to have an opinion about.

    Raises:
        StepFailed: naming the command, its exit status and its stderr.
    """
    result = subprocess.run(argv, capture_output=True, text=True)
    if result.returncode != 0:
        detail = (result.stderr or result.stdout or "").strip().splitlines()
        raise StepFailed(
            f"{what} failed: `{' '.join(argv)}` exited {result.returncode}"
            + (f" ({detail[0]})" if detail else "")
        )
    return result


def clear_stage(name: str, stage_root: str = "."):
    """Remove any leftover staged database before an update script writes a new one.

    A stage left behind by a crashed run is the one way a validated stage could still be wrong:
    a stale `clinvar_20240101.vcf` sitting beside a fresh `clinvar_20260820.vcf` gives two files
    matching the same pattern, and the validator refuses that as an ambiguous database rather than
    guessing. Clearing first means a run's stage contains only what this run put there.
    """
    stage = os.path.join(stage_root, name)
    if os.path.isdir(stage):
        _run_checked(["rm", "-rf", stage], f"clearing the leftover {name} stage")


def staged_builds(name: str, builds, stage_root: str = "."):
    """Which of the requested builds the update script actually staged something for.

    The distinction this draws is the one `update_dbsnp` needs and the rest of the route inherits:
    a build with no staged directory was either deliberately not attempted or failed, and those
    are different verdicts. Deciding which is the caller's job, because only the caller knows
    whether the script said it was skipping.
    """
    return [b for b in builds if os.path.isdir(os.path.join(stage_root, name, b))]


def live_version(name: str, build: str, db_dir: str, spec_table: dict = None):
    """The version the live database reports for one build, or None if it is not installed.

    Read from the same config field the staged copy is checked against, so that `UPDATED: <old> ->
    <new>` names two versions that were both actually read off disk -- rather than one read and one
    assumed, which is how `SUCCESS: hgnc updated to <today>` came to be printed on runs that
    installed nothing.
    """
    table = check_stage_spec(spec_table)
    config = table[name]["config"]
    if not config:
        return None
    return read_config_field(os.path.join(db_dir, name, build, config), "version")


def swap_build(
    name: str,
    build: str,
    db_dir: str,
    backup_dir: str,
    stage_root: str = ".",
    db_germline_dir: str = None,
):
    """Replace ONE build of the live database with the staged one, backing it up first.

    Call only after `_require_valid_stage` has passed for this build -- that ordering IS the
    contract, and the guard below is a backstop for a future caller who forgets it rather than a
    substitute for it.

    The unrequested build is never named here, which is the whole of #357 ruling 2: there is no
    step in this function that could touch it, so no restore-from-backup is needed to put it back.
    Ruling 3 follows from the same shape -- only the build being replaced is backed up, so a
    single-build run does not copy a whole-genome VCF for a build that is not moving.

    Germline is written from the SAME staged directory somatic just received, which is what makes
    it derivable rather than independent (#352 ruling 6) -- and is why the caller must have proved
    the two were in sync first, because the only backup taken here is somatic's.
    """
    staged = os.path.join(stage_root, name, build)
    if not os.path.isdir(staged):
        raise StepFailed(
            f"{name} {build}: asked to install a staged build that does not exist at {staged}; "
            "refusing, because the live database would be replaced by nothing"
        )

    live = os.path.join(db_dir, name, build)
    if os.path.isdir(live):
        backup = os.path.join(backup_dir, name)
        _run_checked(["mkdir", "-p", backup], f"creating {backup}")
        _run_checked(["cp", "-r", live, backup + os.sep], f"backing up {name} {build}")

    _run_checked(["rm", "-rf", live], f"removing the old {name} {build}")
    _run_checked(["mkdir", "-p", os.path.dirname(live)], f"creating {os.path.dirname(live)}")
    _run_checked(["cp", "-r", staged, live], f"installing {name} {build}")

    if db_germline_dir:
        germline = os.path.join(db_germline_dir, name, build)
        _run_checked(["rm", "-rf", germline], f"removing the old germline {name} {build}")
        _run_checked(
            ["mkdir", "-p", os.path.dirname(germline)],
            f"creating {os.path.dirname(germline)}",
        )
        _run_checked(["cp", "-r", staged, germline], f"installing germline {name} {build}")


def run_update_script(script: str, args=()):
    """Run one of the `update_db` shell scripts under an interpreter that can execute it.

    The single place the Funcotator route names an interpreter, which is the point. All nine shell
    scripts under `update_db/` are bash-shebanged and none is POSIX, and in the container
    `/bin/sh` is dash -- so handing one to `sh` broke it in one of two ways, both silent:

    * `function name {` is a **parse** error, so the script exited 2 having done nothing. That is
      what kept `get_new_hgnc.sh` accidentally safe.
    * `[[ ... ]]` is a **runtime** command-not-found, so the script CARRIED ON with every
      conditional false. `update_clinvar_funcotator.sh` is the only live callee with no `function`
      keyword, so it took this path: it exited 0 having created nothing but an empty `clinvar/`,
      and that empty directory was copied over both live ClinVar databases.

    Keeping the interpreter in one function means a new step cannot get it wrong without adding a
    shell-out that the interpreter guard in `update_db/tests/` will see as a new call site.
    """
    return subprocess.run(["bash", script, *args], capture_output=True, text=True)


def stage_and_validate(
    name: str,
    script: str,
    script_args,
    builds,
    run_log: RunLog,
    stage_root: str = ".",
    spec_table: dict = None,
):
    """Run an update script into a clean stage and inspect what it produced.

    Returns:
        tuple: `(versions, result)` -- the resolved version per build that the script staged, and
            the `CompletedProcess`, whose exit status becomes the detail line on a failure.

    Raises:
        StepFailed: if the stage does not pass validation, with the script's exit status appended.
            The live database has not been touched at this point and will not be.
    """
    clear_stage(name, stage_root)
    result = run_update_script(script, script_args)
    present = staged_builds(name, builds, stage_root)
    try:
        versions = _require_valid_stage(name, present, stage_root, spec_table)
    except StepFailed as exc:
        raise StepFailed(f"{exc} The update script exited {result.returncode}.") from exc
    return versions, result


# =============================================================================================
# The failure contract, Annovar side
# =============================================================================================
#
# Issue #354. The shared primitives above are shaped for the Funcotator route, whose scripts
# stage a database into `<stage>/<db>/<build>/` and let Python copy it over the live tree. The
# Annovar route's layout and its risks are both different, and the difference is worth stating
# once here rather than rediscovering it at each call site.
#
# **The build is a filename prefix, not a directory.** Annovar keeps everything in one flat
# `humandb`, as `<build>_<database>.txt` -- `hg19_avsnp151.txt`, `hg38_clinvar_20260102.txt`.
# Every write lands on a NEW filename, because the version is part of the name, so an update
# never overwrites or removes the database it replaces. That is why this route cannot destroy an
# annotation database the way the Funcotator route can, and why there is no germline tree, no
# backup directory and no `rm -rf` anywhere in it (#357).
#
# **So the failure this route actually has is adoption, not destruction.** A database is "in use"
# because two files NAME it: the annotation script's hardcoded `-protocol` list (`CancerVar.py`,
# `Intervar.py`) and the config's `database_names` line. Downloading a file changes nothing until
# those names move. The three ways that goes wrong, all of them live before this ticket:
#
#   1. the names move to a database whose file was never written, so annovar is pointed at
#      nothing;
#   2. the names DO NOT move although the log says the update succeeded, so annotation quietly
#      carries on with the previous database -- which still exists, because this route is
#      additive. This is what the four GNU-only `sed -i 's/.../.../g'` calls did on macOS: BSD
#      `sed` reads the expression as `-i`'s backup suffix, exits 1, and leaves the file
#      UNCHANGED, and nothing checked the result;
#   3. the two names move independently of each other -- the annotation script was rewritten
#      per-database as the loop went while the config was written once at the end, so anything
#      raising in between left them disagreeing.
#
# What that disagreement costs, established by reading the consumers rather than assumed:
# `database_names` is read in exactly one place, `check_downdb` (`CancerVar.py:636`,
# `Intervar.py:509`), which pre-downloads any named database whose file is missing. The
# annotation itself runs off the hardcoded `-protocol` string. So a config that names MORE than
# the protocol costs a wasted download, while a protocol that names a database the config does
# not is the harmful direction -- nothing pre-fetches it. And case 2, the silent no-op, is the
# worst of the three precisely because nothing breaks: every column is populated, from the
# previous release, under a log line claiming the database was updated.
#
# Hence the shape of everything below: validate the file, then move BOTH names together, and
# never report an update the names did not actually reach.

# `index_annovar.pl` writes `#BIN\t<bin size>\t<size of the indexed file>` as its first line, and
# the databases annovar distributes carry the same header -- the script says so, and
# `annotate_variation.pl` reads it back that way.
_ANNOVAR_INDEX_HEADER = re.compile(r"BIN\t(\d+)\t(\d+)")


def annovar_db_paths(annovar_db_path: str, build: str, database: str):
    """The two files annovar reads for one database of one build: the table and its index."""
    data = os.path.join(annovar_db_path, f"{build}_{database}.txt")
    return data, data + ".idx"


def describe_annovar_index(data_path: str):
    """Why annovar will not use this database's index, or None if it will.

    MEASURED against the consumer, `annotate_variation.pl:2324-2364`, and the answer decides
    whether a bad index is a failure or a footnote: annovar requires the index's first line to
    match `BIN\\t<int>\\t<int>` and the size it declares to equal the actual size of the data
    file, and when either is wrong it prints *"ANNOVAR can still generate correct results without
    index file"* and falls back to scanning the table. A missing index takes the same path.

    So the index is an optimisation, not part of the database's correctness, and this returns
    prose for a detail line rather than raising. The inventory's open question about
    `annotate_variation.pl -downdb` leaving a half-updated database is answered here: the half
    that matters is the `.txt`, and it is checked by `require_valid_annovar_db`.

    Returns:
        str: what is wrong with the index, phrased for an operator, or None if it is usable.
    """
    index_path = data_path + ".idx"
    if not os.path.isfile(index_path):
        return (
            f"{os.path.basename(index_path)} is absent, so annovar will scan the whole table "
            "instead of seeking into it; results are unaffected"
        )
    try:
        with open(index_path, errors="replace") as handle:
            first = handle.readline()
    except OSError as exc:
        return f"{os.path.basename(index_path)} could not be read ({exc})"

    match = _ANNOVAR_INDEX_HEADER.search(first)
    if not match:
        return (
            f"{os.path.basename(index_path)} does not start with an index header, so annovar "
            "will report it malformed and ignore it (a truncated or empty index is what a "
            "failed indexing run leaves behind, because the redirection creates the file "
            "whether or not the indexer succeeded)"
        )
    declared = int(match.group(2))
    actual = os.path.getsize(data_path)
    if declared != actual:
        return (
            f"{os.path.basename(index_path)} was built against a {declared}-byte table but "
            f"{os.path.basename(data_path)} is {actual} bytes, so annovar will call the index "
            "out of date and ignore it"
        )
    return None


def require_valid_annovar_db(annovar_db_path: str, build: str, database: str):
    """Refuse to adopt an annovar database file that is not one.

    The Annovar counterpart of `_require_valid_stage`, and deliberately NOT driven by a per-step
    table. `STAGE_SPEC` works for the Funcotator route because it has three hand-built stages
    whose shapes are known when the code is written. The databases here are discovered at run
    time from annovar's own published table, so a per-database table could not be complete even
    in principle -- it would go stale on annovar's release schedule rather than on ours, which is
    the vacuity failure mode `check_stage_spec` exists to prevent, in a form no assertion could
    catch. Every annovar database has the same shape instead, so one rule covers all of them and
    there is no table to fall behind.

    Args:
        annovar_db_path (str): the humandb directory the databases live in.
        build (str): 'hg19' or 'hg38'.
        database (str): the versioned database name, e.g. 'avsnp151', 'clinvar_20260102'.

    Returns:
        str: a note about the index if annovar will not be able to use it, else None. It is a
            detail line, not a failure; see `describe_annovar_index`.

    Raises:
        StepFailed: if the table itself is missing, empty, or carries no data.
    """
    data_path, _ = annovar_db_paths(annovar_db_path, build, database)
    name = os.path.basename(data_path)

    if not os.path.isfile(data_path):
        raise StepFailed(
            f"{database}: {data_path} was not created, so there is no {build} database to "
            "annotate with and the names in the annotation script and the config are unchanged"
        )
    if os.path.getsize(data_path) == 0:
        raise StepFailed(
            f"{database}: {name} is 0 bytes, so it is a database annovar will read and "
            "annotate nothing from"
        )
    try:
        record = first_data_record(data_path, "hash")
    except OSError as exc:
        raise StepFailed(f"{database}: {name} could not be read ({exc})")
    if record is None:
        raise StepFailed(
            f"{database}: {name} has a header but not one data record, so every variant would "
            "be annotated with no evidence from it"
        )
    return describe_annovar_index(data_path)


def adopt_annovar_database(annotation_script: str, current: str, new: str):
    """Point the annotation script at `new` wherever it currently names `current`.

    This replaces `sed -i 's/<current>/<new>/g'`, which was GNU-only: on macOS BSD `sed` takes
    the expression as `-i`'s backup suffix, exits 1 and leaves the file untouched, and no call
    site inspected the result. Doing the substitution in Python removes the host split entirely
    and -- the point of the exercise -- makes the outcome checkable, so a substitution that did
    not happen is a `StepFailed` rather than a `SUCCESS` line over an unchanged file.

    The replacement is literal, where `sed`'s was a regular expression. That is strictly safer
    for database names and it is why no escaping is needed here.

    Args:
        annotation_script (str): `CancerVar.py` or `Intervar.py`.
        current (str): the versioned name the script names today, e.g. 'avsnp151'.
        new (str): the versioned name to adopt, e.g. 'avsnp152'.

    Returns:
        int: how many occurrences were rewritten. Zero when there was nothing to do because the
            script already named `new`.

    Raises:
        StepFailed: if `current` is not a versioned name, if one name contains the other, if the
            script names neither version, or if the rewrite did not take.
    """
    if current == new:
        return 0
    if not re.search(r"\d", current):
        # `get_dbname_version` returns an empty version for a database whose name carries no
        # digits, which makes `current` a bare prefix -- and a global replace of a prefix
        # corrupts every longer name built on it, including the one being written. No live
        # configuration reaches this, and it is a refusal rather than a silent corruption.
        # Checked before the containment rule below because it is the more actionable of the two
        # messages, and a versionless name is almost always a prefix of its replacement as well.
        raise StepFailed(
            f"{current!r} carries no version, so rewriting every occurrence of it to {new!r} "
            "would also rewrite the names it is a prefix of. Add it to --not_to_update, or give "
            "the update a version to match on"
        )
    if current in new:
        # A global replace of a name by something that still contains it is not well defined: the
        # result names the new database and the old one at the same time, and running it twice
        # extends the name twice. `dbnsfp47a` -> `dbnsfp47a_interpro` is the shape in this tree,
        # kept apart today only by `get_dbname_version` splitting them into separate entries.
        raise StepFailed(
            f"renaming {current} to {new} would leave the old name inside the new one, so "
            "rewriting every occurrence of it is not a well-defined substitution"
        )

    with open(annotation_script, errors="replace") as f:
        before = f.read()

    if current not in before:
        if new in before:
            return 0
        raise StepFailed(
            f"{os.path.basename(annotation_script)} names neither {current} nor {new}, so this "
            "update has nothing to adopt and the database would be downloaded but never used"
        )

    after = before.replace(current, new)
    with open(annotation_script, "w") as f:
        f.write(after)

    # The check the shell-out could not do: `sed` reported nothing and was asked nothing, so the
    # previous implementation could not tell a substitution from a no-op. Belt and braces rather
    # than a live check -- with a literal in-process replace and the two refusals above, there is
    # no input that reaches here and fails, which is why the mutation harness lists it as
    # unobservable rather than uncovered. It is the assertion that the file ON DISK names the new
    # database, so it keeps its value if this is ever implemented by shelling out again.
    with open(annotation_script, errors="replace") as f:
        written = f.read()
    if current in written or new not in written:
        raise StepFailed(
            f"rewriting {os.path.basename(annotation_script)} from {current} to {new} did not "
            "take: the file still names the old database. The annotation would keep using it "
            "while the log claimed an update"
        )
    return before.count(current)


def rewrite_database_names(configs: list, renames: dict):
    """Substitute the adopted names into `database_names`, keeping the line's order.

    CancerVar and Intervar are order-sensitive in `database_names`: it lines up with the
    `-operation` codes in the annotation script. The previous implementation rebuilt the line by
    partitioning the databases into "not updated" and "updated" and concatenating, which reorders
    it -- for the CancerVar config it swapped `dbnsfp47a_interpro` past `gnomad_genome`, and the
    manual swap of exactly those two entries that followed existed to put them back. Substituting
    each name where it already stands cannot reorder anything, so that repair has nothing left to
    repair and is gone.

    It also makes the ticket's requirement structural rather than something to remember: a
    database whose update failed simply is not in `renames`, so its OLD name stays on the line,
    and no name can appear that a validated file does not stand behind.

    Args:
        configs (list): the config file's lines, as read.
        renames (dict): old versioned name -> newly adopted name, for adopted updates only.

    Returns:
        list: the lines with the `database_names` line rewritten in place.

    Raises:
        StepFailed: if the config has no `database_names` line to rewrite.
    """
    out = list(configs)
    for idx, line in enumerate(out):
        if line.startswith("database_names"):
            names = line.split("=", 1)[1].split()
            out[idx] = "database_names = " + " ".join(renames.get(n, n) for n in names) + "\n"
            return out
    raise StepFailed(
        "the config file has no `database_names` line, so there is nothing to record this "
        "update in and the databases would be downloaded but never used"
    )
