#!/usr/bin/env python
import argparse
import datetime
import os

from utils_update import (
    FolderType,
    FileType,
    get_annovar_databases,
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
        required=True,
        help="Email of your cosmic account (cancer.sanger.ac.uk).",
    )
    parser.add_argument(
        "-cp",
        "--cosmic_password",
        type=str,
        required=True,
        help="Password of your cosmic account (cancer.sanger.ac.uk).",
    )
    args = parser.parse_args()
    return args


def main():
    """bla"""
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
        "dbnsfp31a_interpro",
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

    with open(logfile, "w") as f:
        f.write("Starting annovar db update for cancervar. \n")
    update_all_cancervar(
        annovar_databases=annovar_databases,
        annovar_db_path=args.annovar_db_path,
        annovar_download_script=args.annovar_download_script,
        cancervar_config_file=cancervar_config_file,
        cancervar_script=cancervar_script,
        email=args.cosmic_email,
        logfile=logfile,
        not_to_update=not_to_update,
        password=args.cosmic_password,
        scriptdir=scriptdir,
        today=today,
        vt=args.vt,
    )

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

    with open(logfile, "a") as f:
        f.write("Starting annovar db update for intervar. \n")

    update_all_intervar(
        annovar_databases=annovar_databases,
        annovar_db_path=args.annovar_db_path,
        annovar_download_script=args.annovar_download_script,
        intervar_config_file=intervar_config_file,
        intervar_script=intervar_script,
        logfile=logfile,
        not_to_update=not_to_update,
        scriptdir=scriptdir,
        today=today,
        vt=args.vt,
    )


if __name__ == "__main__":
    main()
