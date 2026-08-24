#!/usr/bin/env python
import argparse
import datetime
import os
import sys

from utils_update import (
    FolderType,
    FileType,
    RunLog,
    get_annovar_databases,
    run_steps,
    update_all_cancervar,
    update_all_intervar,
)


def _parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-a",
        "--annovar_db_path",
        type=FolderType(),
        default=None,
        required=True,
        help="Directory where to save the annovar database. Tipycally path/humandb",
    )
    parser.add_argument(
        "-cc",
        "--config_cancervar",
        type=FileType(),
        default=None,
        required=False,
        help="Current cancervar config file with line database_names = comma_separated_databases. If not defined, it takes the one from repository if cloned",
    )
    parser.add_argument(
        "-cs",
        "--cancervar_script",
        type=FileType(),
        default=None,
        required=False,
        help="CancerVar.py. If not defined, it takes the one from repository if cloned",
    )
    parser.add_argument(
        "-is",
        "--intervar_script",
        type=FileType(),
        default=None,
        required=False,
        help="Intervar.py. If not defined, it takes the one from repository if cloned",
    )
    parser.add_argument(
        "-ci",
        "--config_intervar",
        type=FileType(),
        default=None,
        required=False,
        help="Current inter config file with line database_names = comma_separated_databases. If not defined, it takes the one from repository if cloned",
    )
    parser.add_argument(
        "-as",
        "--annovar_download_script",
        type=FileType(),
        default=None,
        required=True,
        help="Path to annovar the script annotate_variation.pl",
    )
    parser.add_argument(
        "-s",
        "--scriptdir",
        type=FolderType(),
        required=False,
        default=None,
        help="Directory containing all the scripts for updating databases. If None (default) it uses the one in the repo",
    )
    parser.add_argument(
        "-ad",
        "--annovar_docs",
        type=str,
        default="https://raw.githubusercontent.com/WGLab/doc-ANNOVAR/master/docs/user-guide/download.md",
        required=False,
        help="Link to annovar download docs on github. Default https://raw.githubusercontent.com/WGLab/doc-ANNOVAR/master/docs/user-guide/download.md",
    )
    parser.add_argument(
        "-v",
        "--vt",
        type=FolderType(),
        default=None,
        required=True,
        help="Folder containing vt software (variant harmonisation).",
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
    args = parser.parse_args()
    return args


def main():
    """Update the annovar databases CancerVar and Intervar annotate with.

    Returns:
        int: the exit code -- 1 if any step failed, else 0. Previously this function had no
            `try`/`except` anywhere and set no exit code at all, so a `FileExistsError` raised
            while updating CancerVar's databases meant Intervar was never updated, and a caller
            checking `$?` saw the same 0 either way.
    """
    args = _parse_args()

    if not args.scriptdir:
        scriptdir = f"{os.path.dirname(__file__)}/../annovar4cancervar_intervar"
    else:
        scriptdir = args.scriptdir

    today = datetime.date.today().strftime("%Y%m%d")
    logfile = f"{today}_cancervar_intervar_update.log"

    # get databases from annovar
    annovar_databases = get_annovar_databases(args.annovar_docs)

    # list of databases not to update
    not_to_update = [
        "refGene",
        "ensGene",
        "knownGene",
        "esp6500siv2_all",
        "1000g2015aug",
        "exac03",
        "dbscsnv11",
        # "dbnsfp31a_interpro",
        "gnomad_genome",
        "rmsk",
    ]

    if not args.config_cancervar:
        cancervar_config_file = (
            f"{os.path.dirname(__file__)}/../../resources/configs/config.init.CancerVar"
        )
    else:
        cancervar_config_file = args.config_cancervar

    if not args.cancervar_script:
        cancervar_script = (
            f"{os.path.dirname(__file__)}/../../resources/CancerVar/CancerVar.py"
        )
    else:
        cancervar_script = args.cancervar_script

    if not args.config_intervar:
        intervar_config_file = (
            f"{os.path.dirname(__file__)}/../../resources/configs/config.init.intervar"
        )
    else:
        intervar_config_file = args.config_intervar

    if not args.intervar_script:
        intervar_script = (
            f"{os.path.dirname(__file__)}/../../resources/InterVar/Intervar.py"
        )
    else:
        intervar_script = args.intervar_script

    with open(logfile, "w") as f:
        f.write("Starting annovar db update for cancervar and intervar. \n")

    run_log = RunLog(logfile)

    # One step per route. Each reports its own verdicts and, if it fails, does so without taking
    # the other down with it -- which is the whole of ruling 3 as it applies here.
    steps = [
        (
            "cancervar",
            lambda: update_all_cancervar(
                annovar_databases=annovar_databases,
                annovar_db_path=args.annovar_db_path,
                annovar_download_script=args.annovar_download_script,
                cancervar_config_file=cancervar_config_file,
                cancervar_script=cancervar_script,
                email=args.cosmic_email,
                run_log=run_log,
                not_to_update=not_to_update,
                password=args.cosmic_password,
                scriptdir=scriptdir,
                today=today,
                vt=args.vt,
            ),
        ),
        (
            "intervar",
            lambda: update_all_intervar(
                annovar_databases=annovar_databases,
                annovar_db_path=args.annovar_db_path,
                annovar_download_script=args.annovar_download_script,
                intervar_config_file=intervar_config_file,
                intervar_script=intervar_script,
                run_log=run_log,
                not_to_update=not_to_update,
                scriptdir=scriptdir,
                today=today,
                vt=args.vt,
            ),
        ),
    ]

    # This route has no `--build` flag and never asks for one: it discovers each database's builds
    # from annovar's published table, and the build is a filename prefix rather than a directory,
    # so nothing here can lose the build it was not asked for (#357). Both builds are named
    # because both are what a run covers.
    #
    # Note what that means for the two rules in `_settle_step`, so nobody reads more into this
    # than is there: the Annovar verdicts are per DATABASE, not per build -- the annotation
    # script's `-protocol` list carries no build, so a database's name advances once for both or
    # not at all. So ruling 5's "skipped every build it was asked for" rule cannot fire on this
    # route, and the rule doing the work here is the other one: a step that returns without
    # recording a verdict is a failure.
    return run_steps(steps, run_log, ["hg19", "hg38"])


if __name__ == "__main__":
    sys.exit(main())
