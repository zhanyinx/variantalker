#!/usr/bin/env python
import argparse
import datetime
import os
import sys

from utils_update import (
    FolderType,
    FileType,
    RunLog,
    check_oncotator,
    check_dna_repair_genes,
    check_oreganno,
    get_version,
    requested_builds,
    run_steps,
    update_gencode,
    update_cosmic,
    update_clinvar,
    update_dbsnp,
    update_acmg_rec,
    update_hgnc,
)


def _parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-sd",
        "--somatic_database",
        type=FolderType(),
        default=None,
        required=True,
        help="Path to the directory containing funcotator somatic databases",
    )
    parser.add_argument(
        "-ce",
        "--cosmic_email",
        type=str,
        required=False,
        default=None,
        help="Email of your cosmic account (cancer.sanger.ac.uk). If not provided, COSMIC update will be skipped.",
    )
    parser.add_argument(
        "-cp",
        "--cosmic_password",
        type=str,
        required=False,
        default=None,
        help="Password of your cosmic account (cancer.sanger.ac.uk). If not provided, COSMIC update will be skipped.",
    )
    parser.add_argument(
        "-gd",
        "--germline_database",
        type=FolderType(),
        default=None,
        required=True,
        help="Path to the directory containing funcotator germline databases.",
    )
    parser.add_argument(
        "-b",
        "--backup",
        type=str,
        default="/hpcnfs/scratch/DIMA/db_versions",
        help="Folder to save former version of databases.",
    )
    parser.add_argument(
        "-s",
        "--scriptdir",
        required=False,
        type=FolderType(),
        default=None,
        help="Directory containing all the scripts for updating databases. If None (default) it uses the one in the repo",
    )
    parser.add_argument(
        "--build",
        type=str,
        choices=["hg19", "hg38", "both"],
        default="both",
        help="Genome build to update: hg19, hg38, or both (default: both)",
    )
    args = parser.parse_args()
    return args


def main():
    args = _parse_args()

    outfile = f"{datetime.date.today().strftime('%Y%m%d')}_funcotator_update.log"

    if not os.path.exists(args.backup):
        os.mkdir(args.backup)
    backup_dir = f"{args.backup}/{datetime.date.today().strftime('%Y%m%d')}"
    if not os.path.exists(backup_dir):
        os.mkdir(backup_dir)

    if not args.scriptdir:
        scriptdir = f"{os.path.dirname(__file__)}/../funcotator"
    else:
        scriptdir = args.scriptdir

    with open(outfile, "w") as f:
        f.write("Starting funcotator db update. \n")

    run_log = RunLog(outfile)
    builds = requested_builds(args.build)

    # The version comparisons the two disabled steps report, read from a build that is actually
    # installed. `build_to_check` used to be "the build that was requested", which is also why no
    # run ever noticed a missing hg19: once a `--build hg38` run had deleted it, every later
    # `--build hg38` run looked only at hg38 and found everything in order.
    def installed_version(root, database, config):
        for build in builds + ["hg38", "hg19"]:
            path = f"{root}/{database}/{build}/{config}"
            if os.path.exists(path):
                return get_version(path)
        return None

    cosmic_version = installed_version(args.somatic_database, "cosmic", "cosmic.config")
    gencode_version = installed_version(args.somatic_database, "gencode", "gencode.config")
    acmg_version = installed_version(args.germline_database, "acmg_rec", "acmg_rec.config")

    # Every step runs, every step reports its own verdict, and a failure carries on to the next
    # database instead of aborting the run. That is only safe because a failed step is now a no-op
    # on the live database: before the staged-database validation landed, carrying on past a
    # failure risked compounding the damage. The exit code is the machine-readable signal -- the
    # log itself is for a human, and nothing in this repository parses it.
    steps = [
        # Achilles, cancer-gene-census, familial & simple-uniprot (same source dataset, check only)
        # ("oncotator", lambda: check_oncotator(run_log)),
        ("dna_repair_genes", lambda: check_dna_repair_genes(run_log)),
        # ("oreganno", lambda: check_oreganno(run_log)),
        #
        # COSMIC update - DISABLED due to format changes in COSMIC database. The COSMIC database
        # has changed its download format and authentication method, so a manual update is
        # required and the installed version is logged.
        (
            "cosmic",
            lambda: run_log.skipped(
                "cosmic",
                "update disabled - COSMIC changed its download format and authentication method, "
                f"so a manual update is required. Installed: {cosmic_version or 'unknown'}",
            ),
        ),
        # Commented out - uncomment when COSMIC update script is fixed for new format
        # (
        #     "cosmic",
        #     lambda: update_cosmic(
        #         curr_version=cosmic_version,
        #         db_dir=args.somatic_database,
        #         backup_dir=backup_dir,
        #         run_log=run_log,
        #         scriptdir=scriptdir,
        #         email=args.cosmic_email,
        #         password=args.cosmic_password,
        #         build=args.build,
        #     ),
        # ),
        #
        # GENCODE update - DISABLED (not working correctly). The database created seems correct
        # but funcotator does not use it properly.
        (
            "gencode",
            lambda: run_log.skipped(
                "gencode",
                "update disabled - funcotator does not use the rebuilt database correctly, so a "
                f"manual update is required. Installed: {gencode_version or 'unknown'}",
            ),
        ),
        # (
        #     "gencode",
        #     lambda: update_gencode(
        #         curr_version=gencode_version,
        #         db_dir=args.somatic_database,
        #         backup_dir=backup_dir,
        #         run_log=run_log,
        #         scriptdir=scriptdir,
        #         db_germline_dir=args.germline_database,
        #         build=args.build,
        #     ),
        # ),
        (
            "clinvar",
            lambda: update_clinvar(
                db_dir=args.somatic_database,
                backup_dir=backup_dir,
                run_log=run_log,
                scriptdir=scriptdir,
                db_germline_dir=args.germline_database,
                build=args.build,
            ),
        ),
        # No bare `except:` here any more. It reported "There was a problem with the curl download"
        # for every possible HGNC failure -- including the `IndexError` from a glob that could
        # never match, where no download had been attempted. The step now raises failures that say
        # what was checked, and `run_steps` keeps a backstop that logs type, message and traceback.
        (
            "hgnc",
            lambda: update_hgnc(
                db_dir=args.somatic_database,
                backup_dir=backup_dir,
                run_log=run_log,
                scriptdir=scriptdir,
                build=args.build,
            ),
        ),
        (
            "dbsnp",
            lambda: update_dbsnp(
                run_log=run_log,
                db_dir=args.somatic_database,
                backup_dir=backup_dir,
                scriptdir=scriptdir,
                db_germline_dir=args.germline_database,
                build=args.build,
            ),
        ),
        (
            "acmg_rec",
            lambda: update_acmg_rec(
                run_log=run_log,
                db_germline_dir=args.germline_database,
                backup_dir=backup_dir,
                current_version=acmg_version,
                build=args.build,
            ),
        ),
    ]

    return run_steps(steps, run_log, builds)


if __name__ == "__main__":
    sys.exit(main())
