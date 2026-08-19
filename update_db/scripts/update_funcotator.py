#!/usr/bin/env python
import argparse
import datetime
import os

from utils_update import (
    FolderType,
    FileType,
    check_oncotator,
    check_dna_repair_genes,
    check_oreganno,
    get_version,
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
        help="Path to the directory containing funcotator germline databases. If not provided, germline databases will not be updated.",
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
    # Achilles, cancer-gene-census, familial & simple-uniprot (same source dataset, check only)
    # check_oncotator(outfile)
    check_dna_repair_genes(outfile)
    #check_oreganno(outfile)

    # COSMIC update - DISABLED due to format changes in COSMIC database
    # The COSMIC database has changed its download format and authentication method.
    # Manual update required. Current version will be logged.
    build_to_check = "hg38" if args.build in ["hg38", "both"] else "hg19"
    if os.path.exists(f"{args.somatic_database}/cosmic/{build_to_check}/cosmic.config"):
        curr_version = get_version(f"{args.somatic_database}/cosmic/{build_to_check}/cosmic.config")
    else:
        curr_version = "unknown"
    with open(outfile, "a") as f:
        f.write(
            f"SKIPPED: COSMIC update disabled due to format changes in COSMIC database download. "
            f"Manual update required. Current version: {curr_version}\n"
        )
    # Commented out - uncomment when COSMIC update script is fixed for new format
    # update_cosmic(
    #     curr_version=curr_version,
    #     db_dir=args.somatic_database,
    #     backup_dir=backup_dir,
    #     file=outfile,
    #     scriptdir=scriptdir,
    #     email=args.cosmic_email,
    #     password=args.cosmic_password,
    #     build=args.build,
    # )

    # Update gencode - DISABLED (not working correctly). The database created seems correct but funcotator does not use it properly.
    build_to_check = "hg38" if args.build in ["hg38", "both"] else "hg19"
    if os.path.exists(f"{args.somatic_database}/gencode/{build_to_check}/gencode.config"):
        curr_version = get_version(f"{args.somatic_database}/gencode/{build_to_check}/gencode.config")
    else:
        curr_version = None
    with open(outfile, "a") as f:
        f.write(
            f"SKIPPED: GENCODE update disabled due to issues with funcotator using the updated database. "
            f"Manual update required. Current version: {curr_version}\n"
        )
    # update_gencode(
    #     curr_version=curr_version,
    #     db_dir=args.somatic_database,
    #     backup_dir=backup_dir,
    #     file=outfile,
    #     scriptdir=scriptdir,
    #     db_germline_dir=args.germline_database,
    #     build=args.build,
    # )
    
    # Update clinvar
    update_clinvar(
        db_dir=args.somatic_database,
        backup_dir=backup_dir,
        file=outfile,
        scriptdir=scriptdir,
        db_germline_dir=args.germline_database,
        build=args.build,
    )

    # Update hgnc
    try:
        update_hgnc(
            db_dir=args.somatic_database,
            backup_dir=backup_dir,
            file=outfile,
            scriptdir=scriptdir,
            build=args.build,
        )
    except:
        with open(outfile, "a") as f:
            f.write(
                f"FAILED: hgnc updated failed -- There was a problem with the curl download. \n"
            )

    # Update dbsnp
    update_dbsnp(
        file=outfile,
        db_dir=args.somatic_database,
        backup_dir=backup_dir,
        scriptdir=scriptdir,
        db_germline_dir=args.germline_database,
        build=args.build,
    )
  

    # Update acmg_rec
    build_to_check = "hg38" if args.build in ["hg38", "both"] else "hg19"
    if os.path.exists(f"{args.germline_database}/acmg_rec/{build_to_check}/acmg_rec.config"):
        curr_version = get_version(f"{args.germline_database}/acmg_rec/{build_to_check}/acmg_rec.config")
    else:
        curr_version = None
    update_acmg_rec(
        file=outfile,
        db_germline_dir=args.germline_database,
        backup_dir=backup_dir,
        current_version=curr_version,
        build=args.build,
    )


if __name__ == "__main__":
    main()
