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
        "https://www.mdanderson.org/documents/Labs/Wood-Laboratory/human-dna-repair-genes.html"
    )
    soup = BeautifulSoup(r.text, features="lxml")
    string = soup.find("p", attrs={"align": "right"}).text.strip()

    # Compare with current state
    previous_string = "This table was last modified by R. Wood and M. Lowery on Wednesday 10th June 2020"
    if string == previous_string:
        with open(file, "a") as f:
            f.write("SUCCESS: no manual update needed - dna_repair_genes. \n")
    else:
        with open(file, "a") as f:
            f.write(
                f"WARNING: manual update required, found at https://www.mdanderson.org/documents/Labs/Wood-Laboratory/human-dna-repair-genes.html \n"
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


def check_gnomad(file: str):
    """Check gnomad database. (hard coded)"""
    output = subprocess.check_output(
        ["gsutil", "ls", "gs://broad-public-datasets/funcotator/"]
    )
    files = sorted(
        [
            i.removeprefix("gs://broad-public-datasets/funcotator/").rstrip("/")
            for i in output.decode().split()
        ]
    )

    releases = []
    for file in files:
        search = re.search(r"gnomAD_([\d\.]+)_VCF", file)
        if search is not None:
            releases.append(search.groups()[0])

    if len(releases) == 1 and releases[0] == "2.1":
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
):
    """Update cosmic, cosmic_tissue, and cosmic_fusion if a newer version is available on the Sanger Institute website.

    Args:
        db_dir (str): The directory path where the Cosmic database is stored.
        curr_version (str): The current version of the Cosmic database.
        backup_dir (str): The directory path where the backup of the current Cosmic database will be stored.
        file (str): The file path of the log where the outcome of the update process will be recorded.
        scriptdir (str): The directory path where the Cosmic update script is stored.
        email (str): email to log in into cancer.sanger.ac.uk
        password (str): password to login into cancer.sanger.ac.uk
    """
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
):
    """
    Updates the dbsnp database used in the Funcotator annotation.

    Args:
        file (str): Path to the file where success status will be written.
        db_dir (str): Path to the directory containing funcotator databases.
        backup_dir (str): Path to the directory where a backup copy of the current dbsnp files will be stored.
        scriptdir (str): Path to the directory where the external shell script is located.
        db_germline_dir (str, optional): Path to the directory containing the germline dbsnp files. Default is None.

    Raises:
        ValueError: If the dbsnp update fails during somatic Funcotator update.

    """

    # TODO path
    update_script = f"{scriptdir}/update_dbsnp/update_dbsnp.sh"

    # Run update script
    results = subprocess.run(
        f"sh {update_script} -d {db_dir}/dbsnp/",
        shell=True,
        capture_output=True,
    )

    if os.path.exists("dbsnp/hg38/dbSNP.config"):
        subprocess.run(
            f"cp -r {db_dir}/dbsnp {backup_dir}; rm -rf {db_dir}/dbsnp",
            shell=True,
            stdout=subprocess.DEVNULL,
        )

        subprocess.run(
            f"cp -r dbsnp {db_dir}/",
            shell=True,
            stdout=subprocess.DEVNULL,
        )

        if db_germline_dir:
            subprocess.run(
                f"rm -rf {db_germline_dir}/dbsnp; cp -r dbsnp {db_germline_dir}",
                shell=True,
                stdout=subprocess.DEVNULL,
            )

        subprocess.run(
            f"rm -rf dbsnp",
            shell=True,
            stdout=subprocess.DEVNULL,
        )
        version = get_version("dbsnp/hg38/dbSNP.config")
        with open(file, "a") as f:
            f.write(f'SUCCESS: dbsnp updated to {version} "\n')
    else:
        with open(file, "a") as f:
            f.write(f"SUCCESS: no manual update needed - dbsnp. \n")


def update_gencode(
    db_dir: str,
    curr_version: str,
    backup_dir: str,
    file: str,
    scriptdir: str,
    db_germline_dir: str = None,
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

    Raises:
        ValueError: If the getGencode update fails during somatic Funcotator update.

    """

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
            if line.startswith("LATEST_HG38_RELEASE"):
                data[idx] = f"LATEST_HG38_RELEASE={version}\n"
        with open(update_script, "w") as f:
            f.writelines(data)

        # Run update script
        results = subprocess.run(["sh", update_script], capture_output=True)
        if results.stderr:
            raise ValueError("getGencode has failed in somatic functorator update!")
        subprocess.run(
            f"cp -r {db_dir}/gencode/hg38 {backup_dir}; rm -rf {db_dir}/gencode/hg38",
            shell=True,
            stdout=subprocess.DEVNULL,
        )
        subprocess.run(
            f"cp -r gencode/hg38 {db_dir}/gencode",
            shell=True,
            stdout=subprocess.DEVNULL,
        )

        if db_germline_dir:
            subprocess.run(
                f"rm -r {db_germline_dir}/gencode/hg38; cp -r gencode/hg38 {db_germline_dir}/gencode",
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
):
    """Updates ClinVar.

    Args:
        db_dir (str): Directory path for the database.
        backup_dir (str): Directory path for the backup.
        file (str): File path to log updates.
        scriptdir (str): Directory path to the script.
        db_germline_dir (str, optional): Directory path for the germline database. Defaults to None.
    """

    # TODO remove hard code script
    update_script = f"{scriptdir}/update_clinvar/update_clinvar_funcotator.sh"

    # Run update script
    results = subprocess.run(["sh", update_script], capture_output=True)

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


def update_hgnc(db_dir: str, backup_dir: str, file: str, scriptdir: str):
    """Updates HGNC.

    Args:
        db_dir (str): Directory path for the database.
        backup_dir (str): Directory path for the backup.
        file (str): File path to log updates.
        scriptdir (str): Directory path to the script.
    """

    today = datetime.date.today().strftime("%b%d%Y")
    # TODO update path
    update_script = f"{scriptdir}/update_hgnc/get_new_hgnc.sh"

    subprocess.run(["sh", update_script], stdout=subprocess.DEVNULL)

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
    file: str, db_germline_dir: str, backup_dir: str, current_version: str
):
    """Update acmg_rec Funcotator database."""

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
        with open(f"{db_germline_dir}/acmg_rec/hg38/acmg_rec.config", "r") as f:
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

        subprocess.run(
            f"mkdir -p acmg_rec/hg38; mv acmg_{version}_{today}_test_cleaned.txt acmg_rec.config acmg_rec/hg38; cp -r acmg_rec/hg38 acmg_rec/hg19",
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
        select = annovar_databases[
            [name in x for x in annovar_databases["name"].values]
        ].sort_values("name", ascending=False)

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

    Args:
        annovar_db_path: Output folder where to save databases files. In general humandb.
        curr_version: The current version of the Cosmic database.
        email: email to log in into cancer.sanger.ac.uk
        logfile: The file path of the log where the outcome of the update process will be recorded.
        password: password to login into cancer.sanger.ac.uk
        scriptdir: The directory path where the Cosmic update script is stored.
    """
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

    # update cosmic
    cosmic = list(filter(lambda x: x.startswith("cosmic"), current_databases))[0]
    version = update_cosmic_annovar(
        annovar_db_path=annovar_db_path,
        annotation_script=cancervar_script,
        curr_version=re.search("^(.*?)(\\d.*)", cosmic)[2],
        email=email,
        logfile=logfile,
        password=password,
        scriptdir=scriptdir,
    )
    updated.append(f"cosmic{version}")

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
    curr_version = list(filter(lambda x: x.startswith("clinvar"), current_databases))[0]
    subprocess.run(
        f"sed -i 's/{curr_version}/clinvar_{today}/g' {cancervar_script}",
        shell=True,
        stdout=subprocess.DEVNULL,
    )
    # re-sorting as cancervar is order sensitive to databases
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
    curr_version = list(filter(lambda x: x.startswith("clinvar"), current_databases))[0]
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
