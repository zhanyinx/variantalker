import argparse
from ftplib import FTP
import base64
import glob
import datetime
import io
import logging
import os
import re
import subprocess
import tempfile
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


def check_oncotator(file: str):
    """Check oncotator database and throw a warning if update is needed. Hard coded."""
    # Get list of files
    data = []

    with FTP("ftp.broadinstitute.org", user="gsapubftp-anonymous") as ftp:
        ftp.cwd("bundle/oncotator")
        ftp.dir(data.append)
    name = "achilles, cancer-gene-census, familial & simple-uniprot"

    # Check for any updates / changes
    if not len(data):
        with open(file, "a") as f:
            f.write(f"WARNING: {name} - no databases found! \n")
            os._exit(1)
    if len(data) > 1:
        with open(file, "a") as f:
            f.write(
                f"WARNING: {name} - NEW database found, please update manually! Current database is oncotator_v1_ds_April052016.tar.gz \n"
            )
            os._exit(1)

    search = re.search("^.+\s+(\d+)\s+(oncotator_.+)$", data[0]).groups()
    if search[1] != "oncotator_v1_ds_April052016.tar.gz":
        with open(file, "a") as f:
            f.write(
                f'WARNING: {name} - latest database from {search[0]} "{search[1]}". Current is oncotator_v1_ds_April052016.tar.gz \n'
            )
    else:
        with open(file, "a") as f:
            f.write(f"SUCCESS: no manual update needed - {name}.\n")


def check_dna_repair_genes(file: str):
    """Check dna repair genes. (hard coded)"""
    # Get website data
    r = requests.get(
        "https://www.mdanderson.org/research/departments-labs-institutes/labs/wood-laboratory/resources.html"
    )
    soup = BeautifulSoup(r.text, features="lxml")
    for line in soup.get_text().splitlines():
        if "Table information was last updated" in line:
            string=line.strip()

    # Compare with current state
    previous_string = "Table information was last updated by R. Wood and B. Dennehey, November 2025."
    if string == previous_string:
        with open(file, "a") as f:
            f.write("SUCCESS: no update needed - dna_repair_genes. \n")
    else:
        with open(file, "a") as f:
            f.write(
                f"WARNING: manual update required for dna_repair_genes, found at https://www.mdanderson.org/research/departments-labs-institutes/labs/wood-laboratory/resources.html \n"
            )


def check_oreganno(file: str):
    """Check oreganno database. (hard coded)"""
    # Parse list of files
    df = (
        pd.read_html("http://www.oreganno.org/dump/")[0][
            ["Name", "Last modified", "Size"]
        ]
        .dropna()
        .reset_index(drop=True)
    )
    if sorted(df["Name"], reverse=True)[0] == "ORegAnno_Combined_2016.01.19.tsv":
        with open(file, "a") as f:
            f.write("SUCCESS: no manual update needed - ORegAnno. \n")
    else:
        with open(file, "a") as f:
            f.write(
                f"WARNING: manual update required, found at http://www.oreganno.org/dump/ \n"
            )


def update_cosmic(
    db_dir: str,
    curr_version: str,
    backup_dir: str,
    file: str,
    scriptdir: str,
    email: str,
    password: str,
    build: str = "both",
):
    """Update cosmic, cosmic_tissue, and cosmic_fusion if a newer version is available on the Sanger Institute website.
    
    NOTE: This function is currently DISABLED due to format changes in COSMIC database.
    COSMIC has changed its download format and authentication method, requiring script updates.

    Args:
        db_dir (str): The directory path where the Cosmic database is stored.
        curr_version (str): The current version of the Cosmic database.
        backup_dir (str): The directory path where the backup of the current Cosmic database will be stored.
        file (str): The file path of the log where the outcome of the update process will be recorded.
        scriptdir (str): The directory path where the Cosmic update script is stored.
        email (str): email to log in into cancer.sanger.ac.uk
        password (str): password to login into cancer.sanger.ac.uk
        build (str): Genome build to update: 'hg19', 'hg38', or 'both' (default: 'both').
    """
    # DISABLED: COSMIC update temporarily disabled due to format changes
    with open(file, "a") as f:
        f.write(
            f"SKIPPED: COSMIC update disabled due to format changes in COSMIC database download. "
            f"Manual update required. Current version: {curr_version}\n"
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

        # Run update script
        subprocess.run(
            f"sh {update_script} {passcode} {version}",
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
        with open(file, "a") as f:
            f.write(f"SUCCESS: Cosmic updated from {curr_version} to {version}. \n")
    else:
        with open(file, "a") as f:
            f.write(
                f"SUCCESS: Cosmic, no need for update. Current version {curr_version} is the latest version {version} \n"
            )


def update_dbsnp(
    file: str,
    db_dir: str,
    backup_dir: str,
    scriptdir: str,
    db_germline_dir: str = None,
    build: str = "both",
):
    """
    Updates the dbsnp database used in the Funcotator annotation.

    Args:
        file (str): Path to the file where success status will be written.
        db_dir (str): Path to the directory containing funcotator databases.
        backup_dir (str): Path to the directory where a backup copy of the current dbsnp files will be stored.
        scriptdir (str): Path to the directory where the external shell script is located.
        db_germline_dir (str, optional): Path to the directory containing the germline dbsnp files. Default is None.
        build (str): Genome build to update: 'hg19', 'hg38', or 'both' (default: 'both').

    Raises:
        ValueError: If the dbsnp update fails during somatic Funcotator update.

    """

    # Convert relative paths to absolute to ensure correct file access
    file = os.path.abspath(file)
    db_dir = os.path.abspath(db_dir)
    backup_dir = os.path.abspath(backup_dir)
    scriptdir = os.path.abspath(scriptdir)
    if db_germline_dir:
        db_germline_dir = os.path.abspath(db_germline_dir)
    
    update_script = f"{scriptdir}/update_dbsnp/update_dbsnp.sh"

    # Run update script with absolute paths and build parameter
    results = subprocess.run(
        ["bash", update_script, "-d", f"{db_dir}/dbsnp/", "--build", build],
        capture_output=True,
        text=True,
    )

    # Check if script output indicates skipped updates
    script_output = results.stdout + results.stderr
    
    # Count how many builds were skipped
    hg38_skipped = "SKIPPED: hg38 dbSNP" in script_output
    hg19_skipped = "SKIPPED: hg19 dbSNP" in script_output
    
    # If both were skipped, raise error
    if hg38_skipped and hg19_skipped:
        with open(file, "a") as f:
            f.write(f"FAILED: dbSNP update failed - both hg38 and hg19 databases are missing.\n")
            for line in script_output.split('\n'):
                if "SKIPPED:" in line or "To initialize" in line:
                    f.write(f"  {line}\n")
        raise ValueError("dbSNP update failed - both hg38 and hg19 databases are missing")
    
    # Check if new database files were created (indicates update occurred)
    files_created = os.path.exists("dbsnp/hg38/dbSNP.config") or os.path.exists("dbsnp/hg19/dbSNP.config")
    
    # Handle partial skips and log appropriately
    if hg38_skipped and not hg19_skipped:
        # hg19 was processed (may or may not have been updated)
        version = None
        if os.path.exists("dbsnp/hg19/dbSNP.config"):
            version = get_version("dbsnp/hg19/dbSNP.config")
        elif os.path.exists(f"{db_dir}/dbsnp/hg19/dbSNP.config"):
            version = get_version(f"{db_dir}/dbsnp/hg19/dbSNP.config")
        
        with open(file, "a") as f:
            f.write(f"WARNING: hg38 dbSNP update skipped - no existing database found.\n")
            if version:
                status = "updated" if os.path.exists("dbsnp/hg19/dbSNP.config") else "already up-to-date"
                f.write(f"INFO: hg19 dbSNP {status} (version: {version}).\n")
            else:
                f.write(f"WARNING: hg19 dbSNP version could not be determined.\n")
                
    elif hg19_skipped and not hg38_skipped:
        # hg38 was processed (may or may not have been updated)
        version = None
        if os.path.exists("dbsnp/hg38/dbSNP.config"):
            version = get_version("dbsnp/hg38/dbSNP.config")
        elif os.path.exists(f"{db_dir}/dbsnp/hg38/dbSNP.config"):
            version = get_version(f"{db_dir}/dbsnp/hg38/dbSNP.config")
        
        with open(file, "a") as f:
            f.write(f"WARNING: hg19 dbSNP update skipped - no existing database found.\n")
            if version:
                status = "updated" if os.path.exists("dbsnp/hg38/dbSNP.config") else "already up-to-date"
                f.write(f"INFO: hg38 dbSNP {status} (version: {version}).\n")
            else:
                f.write(f"WARNING: hg38 dbSNP version could not be determined.\n")

    # If new files were created, update the databases
    if files_created:
        # Get version from whichever build was updated
        if os.path.exists("dbsnp/hg38/dbSNP.config"):
            version = get_version("dbsnp/hg38/dbSNP.config")
        else:
            version = get_version("dbsnp/hg19/dbSNP.config")
        
        # Backup old database
        subprocess.run(
            ["bash", "-c", f"cp -r {db_dir}/dbsnp {backup_dir}; rm -rf {db_dir}/dbsnp"],
            stdout=subprocess.DEVNULL,
        )

        # Copy new database to db_dir
        subprocess.run(
            ["bash", "-c", f"cp -r dbsnp {db_dir}/"],
            stdout=subprocess.DEVNULL,
        )

        # Copy to germline dir if specified
        if db_germline_dir:
            subprocess.run(
                ["bash", "-c", f"rm -rf {db_germline_dir}/dbsnp; cp -r dbsnp {db_germline_dir}"],
                stdout=subprocess.DEVNULL,
            )

        # Clean up local dbsnp directory
        subprocess.run(
            ["bash", "-c", "rm -rf dbsnp"],
            stdout=subprocess.DEVNULL,
        )

        with open(file, "a") as f:
            if not hg38_skipped and not hg19_skipped:
                f.write(f'SUCCESS: dbSNP updated to {version}\n')
            elif hg38_skipped or hg19_skipped:
                # One was skipped but the other was updated
                f.write(f'SUCCESS: dbSNP updated to {version} (one build skipped)\n')
    elif not hg38_skipped and not hg19_skipped:
        # No new files created AND neither was skipped = databases are already up-to-date
        # Get version from existing database
        version = None
        if os.path.exists(f"{db_dir}/dbsnp/hg38/dbSNP.config"):
            version = get_version(f"{db_dir}/dbsnp/hg38/dbSNP.config")
        elif os.path.exists(f"{db_dir}/dbsnp/hg19/dbSNP.config"):
            version = get_version(f"{db_dir}/dbsnp/hg19/dbSNP.config")
        
        with open(file, "a") as f:
            if version:
                f.write(f"SUCCESS: dbSNP is already up-to-date (version: {version}).\n")
            else:
                f.write(f"SUCCESS: dbSNP is already up-to-date.\n")
    # If one was skipped and no update occurred, the message was already logged above in partial skip handling


def update_gencode(
    db_dir: str,
    curr_version: str,
    backup_dir: str,
    file: str,
    scriptdir: str,
    db_germline_dir: str = None,
    build: str = "both",
):
    """
    Updates the gencode database used in the somatic Funcotator annotation.

    Args:
        db_dir (str): Path to the directory containing funcotator databases.
        curr_version (str): Version number of the current gencode release.
        backup_dir (str): Path to the directory where a backup copy of the current gencode files will be stored.
        file (str): Path to the file where success status will be written.
        scriptdir (str): Path to the directory where the external shell script is located.
        db_germline_dir (str, optional): Path to the directory containing the germline gencode files. Default is None.
        build (str): Genome build to update: 'hg19', 'hg38', or 'both' (default: 'both').

    Raises:
        ValueError: If the getGencode update fails during somatic Funcotator update.

    """
    
    # Skip update if curr_version is None (database doesn't exist for requested build)
    if curr_version is None:
        with open(file, "a") as f:
            f.write(f"SKIPPED: gencode update - no existing database found for build {build}.\n")
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

        # Run update script with build parameter
        subprocess.run(["sh", update_script, "--build", build], capture_output=True)
        
        subprocess.run(
            f"cp -r {db_dir}/gencode/ {backup_dir}; rm -rf {db_dir}/gencode/",
            shell=True,
            stdout=subprocess.DEVNULL,
        )
        subprocess.run(
            f"cp -r gencode/ {db_dir}/gencode",
            shell=True,
            stdout=subprocess.DEVNULL,
        )

        if db_germline_dir:
            subprocess.run(
                f"rm -r {db_germline_dir}/gencode; cp -r gencode {db_germline_dir}",
                shell=True,
                stdout=subprocess.DEVNULL,
            )

        subprocess.run(f"rm -r gencode", shell=True, stdout=subprocess.DEVNULL)
        with open(file, "a") as f:
            f.write(f"SUCCESS: gencode updated from {curr_version} to {version}. \n")
    else:
        with open(file, "a") as f:
            f.write(
                f"SUCCESS: gencode, no need for update. Current version {curr_version} is the latest version {version} \n"
            )


def update_clinvar(
    db_dir: str,
    backup_dir: str,
    file: str,
    scriptdir: str,
    db_germline_dir: str = None,
    build: str = "both",
):
    """Updates ClinVar.

    Args:
        db_dir (str): Directory path for the database.
        backup_dir (str): Directory path for the backup.
        file (str): File path to log updates.
        scriptdir (str): Directory path to the script.
        db_germline_dir (str, optional): Directory path for the germline database. Defaults to None.
        build (str): Genome build to update: 'hg19', 'hg38', or 'both' (default: 'both').
    """

    # TODO remove hard code script
    update_script = f"{scriptdir}/update_clinvar/update_clinvar_funcotator.sh"

    # Run update script with build parameter
    results = subprocess.run(["sh", update_script, "--build", build], capture_output=True)

    subprocess.run(
        f"cp -r {db_dir}/clinvar {backup_dir}; rm -rf {db_dir}/clinvar",
        shell=True,
        stdout=subprocess.DEVNULL,
    )

    subprocess.run(f"cp -r clinvar {db_dir}", shell=True, stdout=subprocess.DEVNULL)

    if db_germline_dir:
        subprocess.run(
            f"rm -r {db_germline_dir}/clinvar; cp -r clinvar {db_germline_dir}",
            shell=True,
            stdout=subprocess.DEVNULL,
        )
    subprocess.run(f"rm -rf clinvar", shell=True, stdout=subprocess.DEVNULL)
    with open(file, "a") as f:
        f.write(f"SUCCESS: clinvar updated successfully. \n")


def update_hgnc(db_dir: str, backup_dir: str, file: str, scriptdir: str, build: str = "both"):
    """Updates HGNC.

    Args:
        db_dir (str): Directory path for the database.
        backup_dir (str): Directory path for the backup.
        file (str): File path to log updates.
        scriptdir (str): Directory path to the script.
        build (str): Genome build to update: 'hg19', 'hg38', or 'both' (default: 'both').
    """

    today = datetime.date.today().strftime("%b%d%Y")
    # TODO update path
    update_script = f"{scriptdir}/update_hgnc/get_new_hgnc.sh"

    subprocess.run(["sh", update_script, "--build", build], stdout=subprocess.DEVNULL)

    # get current data
    current_file = glob.glob(f"{db_dir}/hgnc/hg38/hgnc_*.tsv")[0]
    current_data = pd.read_csv(current_file, sep="\t")

    # get new data
    latest_file = glob.glob(f"hgnc/hg38/hgnc_*.tsv")[0]
    latest_data = pd.read_csv(latest_file, sep="\t")

    # if a whole column is NA, it means that some column retrival is wrong from update.
    # if this is the case, update the curl link in the update_script
    assert len(latest_data.dropna()) > 0

    # if current and new data differs, save the new and backup the old
    try:
        assert_frame_equal(current_data, latest_data, check_names=False)
        same = True
    except:
        same = False

    if not same:
        subprocess.run(
            f"cp -r {db_dir}/hgnc {backup_dir}; rm -r {db_dir}/hgnc",
            shell=True,
            stdout=subprocess.DEVNULL,
        )
        subprocess.run(
            f"cp -r hgnc {db_dir}",
            shell=True,
            stdout=subprocess.DEVNULL,
        )
        subprocess.run(
            f"rm -r hgnc",
            shell=True,
            stdout=subprocess.DEVNULL,
        )

        with open(file, "a") as f:
            f.write(f"SUCCESS: hgnc updated to {today}. \n")
    else:
        subprocess.run(
            f"rm -r hgnc",
            shell=True,
            stdout=subprocess.DEVNULL,
        )
        with open(file, "a") as f:
            f.write(
                f"SUCCESS: hgnc, no need for update. Database did not change from last time \n"
            )


def update_acmg_rec(
    file: str, db_germline_dir: str, backup_dir: str, current_version: str, build: str = "both"
):
    """Update acmg_rec Funcotator database.
    
    Args:
        file (str): File path to log updates.
        db_germline_dir (str): Directory path for the germline database.
        backup_dir (str): Directory path for the backup.
        current_version (str): Current version of the database.
        build (str): Genome build to update: 'hg19', 'hg38', or 'both' (default: 'both').
    """
    
    # Skip update if current_version is None (database doesn't exist for requested build)
    if current_version is None:
        with open(file, "a") as f:
            f.write(f"SKIPPED: acmg_rec update - no existing database found for build {build}.\n")
        return

    today = datetime.date.today().strftime("%b%d%Y")

    r = requests.get("https://www.ncbi.nlm.nih.gov/clinvar/docs/acmg/")
    soup = BeautifulSoup(r.text)

    # get version
    version = re.search("ACMG SF (v\d\.\d)", str(soup.findAll("p")))[1]

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

        # write config file
        build_to_use = "hg38" if build in ["hg38", "both"] else "hg19"
        with open(f"{db_germline_dir}/acmg_rec/{build_to_use}/acmg_rec.config", "r") as f:
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

        # Create directories based on build parameter
        if build == "both":
            subprocess.run(
                f"mkdir -p acmg_rec/hg38; mv acmg_{version}_{today}_test_cleaned.txt acmg_rec.config acmg_rec/hg38; cp -r acmg_rec/hg38 acmg_rec/hg19",
                shell=True,
                stdout=subprocess.DEVNULL,
            )
        elif build == "hg38":
            subprocess.run(
                f"mkdir -p acmg_rec/hg38; mv acmg_{version}_{today}_test_cleaned.txt acmg_rec.config acmg_rec/hg38",
                shell=True,
                stdout=subprocess.DEVNULL,
            )
        else:  # hg19
            subprocess.run(
                f"mkdir -p acmg_rec/hg19; mv acmg_{version}_{today}_test_cleaned.txt acmg_rec.config acmg_rec/hg19",
                shell=True,
                stdout=subprocess.DEVNULL,
            )

        subprocess.run(
            f"cp -r {db_germline_dir}/acmg_rec {backup_dir}; rm -r {db_germline_dir}/acmg_rec",
            shell=True,
            stdout=subprocess.DEVNULL,
        )
        subprocess.run(
            f"cp -r acmg_rec {db_germline_dir}",
            shell=True,
            stdout=subprocess.DEVNULL,
        )
        subprocess.run(
            f"rm -r acmg_rec",
            shell=True,
            stdout=subprocess.DEVNULL,
        )

        with open(file, "a") as f:
            f.write(f"SUCCESS: acmg updated to {today}. \n")
    else:
        with open(file, "a") as f:
            f.write(f"SUCCESS: acmg is already up to date. \n")


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
    logfile: FileType(),
    name_version: dict,
) -> list:
    """Update annovar databases.

    Args:
        annovar_databases: A pandas DataFrame with columns name,
            build, version, and date, containing the list of available databases.
        annovar_db_path: Path to the folder where annovar databases are stored.
        annovar_download_script: Path to the script used to download annovar databases.
        annotation_script: Path to the annotation script that needs to be updated.
        logfile: Path to the logfile where updates are logged.
        name_version: A dictionary mapping database names to their versions.

    Returns:
        A list of updated databases.

    Raises:
        FileExistsError: If a downloaded database is not found in the expected location.
    """
    updated_db = []
    # iterate over databases
    for name, version in name_version.items():
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

        # update if most recent is more recent
        if most_recent > "".join([name, version]):
            select = select[select.name == most_recent]

            # iterate on build
            for index, row in select.iterrows():
                logging.info(
                    f"Updating {name} to {most_recent} for build {row['build']}"
                )
                build = row["build"]
                database = row["name"]
                outfile = f"{annovar_db_path}/{build}_{database}.txt"

                # run database creation only if database is missing
                if not os.path.exists(outfile):
                    subprocess.run(
                        f"{annovar_download_script} -buildver {build} -downdb -webfrom annovar {database} {annovar_db_path}",
                        shell=True,
                        stdout=subprocess.DEVNULL,
                    )

                    # Raise if database creation failed
                    if not os.path.exists(outfile):
                        raise FileExistsError(
                            f"{outfile} was not created! Something wrong happened!"
                        )

                    # write log
                    with open(logfile, "a") as f:
                        f.write(
                            f"SUCCESS: {name} updated to {most_recent} for build {build}\n"
                        )
                else:
                    # write log
                    with open(logfile, "a") as f:
                        f.write(
                            f"SUCCESS: The most recent version {most_recent} for build {build} already exists.\n"
                        )

                # Update annotation script
                subprocess.run(
                    f"sed -i 's/{''.join([name, version])}/{database}/g' {annotation_script}",
                    shell=True,
                    stdout=subprocess.DEVNULL,
                )
            updated_db.append(most_recent)
        else:
            with open(logfile, "a") as f:
                f.write(
                    f"SUCCESS: no update needed. The most recent {most_recent} is the current one.\n"
                )
            updated_db.append("".join([name, version]))
    return updated_db


def update_cosmic_annovar(
    annovar_db_path: FolderType(),
    annotation_script: FileType(),
    curr_version: str,
    email: str,
    logfile: FileType(),
    password: str,
    scriptdir: FolderType(),
) -> str:
    """Update cosmic, cosmic_tissue, and cosmic_fusion if a newer version is available on the Sanger Institute website.
    
    NOTE: This function is currently DISABLED due to format changes in COSMIC database.
    COSMIC has changed its download format and authentication method, requiring script updates.

    Args:
        annovar_db_path: Output folder where to save databases files. In general humandb.
        curr_version: The current version of the Cosmic database.
        email: email to log in into cancer.sanger.ac.uk
        logfile: The file path of the log where the outcome of the update process will be recorded.
        password: password to login into cancer.sanger.ac.uk
        scriptdir: The directory path where the Cosmic update script is stored.
    """
    # DISABLED: COSMIC update temporarily disabled due to format changes
    with open(logfile, "a") as f:
        f.write(
            f"SKIPPED: COSMIC (Annovar) update disabled due to format changes in COSMIC database download. "
            f"Manual update required. Current version: {curr_version}\n"
        )
    return curr_version
    
    # Original code commented out - update when COSMIC format is fixed
    # # Check if credentials are provided
    # if not email or not password:
    #     with open(logfile, "a") as f:
    #         f.write(
    #             f"SKIPPED: Cosmic (Annovar) update skipped - no credentials provided. "
    #             f"Please provide --cosmic_email and --cosmic_password to update COSMIC database. "
    #             f"Current version: {curr_version}\n"
    #         )
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

        # Run update script
        subprocess.run(
            f"sh {update_script} {version} {passcode} {annovar_db_path}",
            shell=True,
            stdout=subprocess.DEVNULL,
        )
        with open(logfile, "a") as f:
            f.write(f"SUCCESS: Cosmic updated from {curr_version} to {version}. \n")
    else:
        with open(logfile, "a") as f:
            f.write(
                f"SUCCESS: Cosmic, no need for update. Latest version {version} is already present \n"
            )
    # Update annotation script
    subprocess.run(
        f"sed -i 's/cosmic{curr_version}/cosmic{version}/g' {annotation_script}",
        shell=True,
        stdout=subprocess.DEVNULL,
    )
    return version


def update_all_cancervar(
    annovar_databases: pd.DataFrame,
    annovar_db_path: FolderType(),
    annovar_download_script: FileType(),
    cancervar_config_file: str,
    cancervar_script: FileType(),
    email: str,
    logfile: FileType(),
    not_to_update: list,
    password: str,
    scriptdir: FolderType(),
    today: str,
    vt: FolderType(),
):
    """Update all CancerVar-related annovar databases."""

    # read config file
    with open(cancervar_config_file) as f:
        configs = f.readlines()

    # get current databases
    for idx, line in enumerate(configs):
        if line.startswith("database_names = "):
            current_databases = line.split(" ")[2:]
            break

    # get the once that we can automatically update easily by downloading from annovar website
    to_update = [
        x.rsplit("\n")[0]
        for x in current_databases
        if x.rsplit("\n")[0] not in not_to_update
    ]
    updated = [
        x.rsplit("\n")[0]
        for x in current_databases
        if x.rsplit("\n")[0] in not_to_update
    ]

    to_update = [x for x in to_update if "clinvar" not in x]
    to_update = [x for x in to_update if "cosmic" not in x]
    to_update_name_version = get_dbname_version(to_update)

    # update avsnp, dbnsfp and icgc
    updated.extend(
        multi_update(
            annovar_databases=annovar_databases,
            annovar_db_path=annovar_db_path,
            annovar_download_script=annovar_download_script,
            annotation_script=cancervar_script,
            logfile=logfile,
            name_version=to_update_name_version,
        )
    )

    # COSMIC update - DISABLED due to format changes
    # Keep the existing COSMIC version in the database list
    cosmic = list(filter(lambda x: x.startswith("cosmic"), current_databases))[0]
    curr_cosmic_version = re.search("^(.*?)(\\d.*)", cosmic)[2]
    with open(logfile, "a") as f:
        f.write(
            f"SKIPPED: COSMIC (Annovar) update disabled due to format changes. "
            f"Current version: {curr_cosmic_version}\n"
        )
    # Commented out - uncomment when COSMIC update script is fixed
    # version = update_cosmic_annovar(
    #     annovar_db_path=annovar_db_path,
    #     annotation_script=cancervar_script,
    #     curr_version=curr_cosmic_version,
    #     email=email,
    #     logfile=logfile,
    #     password=password,
    #     scriptdir=scriptdir,
    # )
    # Keep existing version
    updated.append(cosmic)

    # update clinvar
    logging.info(f"Updating clinvar to {today}")

    subprocess.run(
        f"sh {scriptdir}/update_clinvar_annovar.sh -v {vt} --name clinvar_{today} --output {annovar_db_path}",
        shell=True,
        capture_output=True,
    )

    # check clinvar database has been generated
    if os.path.exists(f"{annovar_db_path}/hg19_clinvar_{today}.txt"):
        with open(logfile, "a") as f:
            f.write(f"SUCCESS: clinvar updated to {today}\n")
        updated.append(f"clinvar_{today}")
    else:
        raise FileExistsError(
            f"{annovar_db_path}/hg19_clinvar_{today}.txt has not been created!"
        )

    # Update annotation script
    # remove \n from curr_version
    curr_version = list(filter(lambda x: x.startswith("clinvar"), current_databases))[0].strip()
    subprocess.run(
        f"sed -i 's/{curr_version}/clinvar_{today}/g' {cancervar_script}",
        shell=True,
        stdout=subprocess.DEVNULL,
    )
    # re-sorting as cancervar is order sensitive to databases  

    # Manual resorting by swapping gnomad_genome and *_interpro in updated list
    if "gnomad_genome" in updated and "interpro" in "".join(updated):
        gnomad_idx = updated.index("gnomad_genome")
        interpro_idx = next(
            i for i, db in enumerate(updated) if db.endswith("interpro")
        )
        # Swap positions
        updated[gnomad_idx], updated[interpro_idx] = (
            updated[interpro_idx],
            updated[gnomad_idx],
        )
    
    for idx, line in enumerate(configs):
        if line.startswith("database_names"):
            configs[idx] = f"database_names = {' '.join(updated)}\n"
    with open(cancervar_config_file, "w") as f:
        f.writelines(configs)


def update_all_intervar(
    annovar_databases: pd.DataFrame,
    annovar_db_path: FolderType(),
    annovar_download_script: FileType(),
    intervar_config_file: FileType(),
    intervar_script: FileType(),
    logfile: FileType(),
    not_to_update: list,
    scriptdir: FolderType(),
    today: str,
    vt: FolderType(),
):
    """Update all CancerVar-related annovar databases."""

    # read config file
    with open(intervar_config_file) as f:
        configs = f.readlines()

    # get current databases
    for idx, line in enumerate(configs):
        if line.startswith("database_names = "):
            current_databases = line.split(" ")[2:]
            break

    # get the once that we can automatically update easily by downloading from annovar website
    to_update = [
        x.rsplit("\n")[0]
        for x in current_databases
        if x.rsplit("\n")[0] not in not_to_update
    ]
    updated = [
        x.rsplit("\n")[0]
        for x in current_databases
        if x.rsplit("\n")[0] in not_to_update
    ]

    to_update = [x for x in to_update if "clinvar" not in x]
    to_update = [x for x in to_update if "cosmic" not in x]
    to_update_name_version = get_dbname_version(to_update)

    # update avsnp, dbnsfp and icgc
    updated.extend(
        multi_update(
            annovar_databases=annovar_databases,
            annovar_db_path=annovar_db_path,
            annovar_download_script=annovar_download_script,
            annotation_script=intervar_script,
            logfile=logfile,
            name_version=to_update_name_version,
        )
    )

    # update clinvar
    logging.info(f"Updating clinvar to {today}")
    subprocess.run(
        f"sh {scriptdir}/update_clinvar_annovar.sh -v {vt} --name clinvar_{today} --output {annovar_db_path}",
        shell=True,
        capture_output=True,
    )

    # check clinvar database has been generated
    if os.path.exists(f"{annovar_db_path}/hg19_clinvar_{today}.txt"):
        with open(logfile, "a") as f:
            f.write(f"SUCCESS: clinvar updated to {today}\n")
        updated.append(f"clinvar_{today}")
    else:
        raise FileExistsError(
            f"{annovar_db_path}/hg19_clinvar_{today}.txt has not been created!"
        )

    # Update annotation script
    # Fix intervar clinvar replacement
    curr_version = list(filter(lambda x: x.startswith("clinvar"), current_databases))[0].strip()
    subprocess.run(
        f"sed -i 's/{curr_version}/clinvar_{today}/g' {intervar_script}",
        shell=True,
        stdout=subprocess.DEVNULL,
    )

    # re-sorting as cancervar is order sensitive to databases
    for idx, line in enumerate(configs):
        if line.startswith("database_names"):
            configs[idx] = f"database_names = {' '.join(updated)}\n"
    with open(intervar_config_file, "w") as f:
        f.writelines(configs)
