[![GitHub code licence is MIT](https://img.shields.io/badge/license-MIT-brightgreen.svg)]

# Update scripts to get the most recent version of annovar and funcotator databases

## Contents
- [Contents](#contents)
- [Overview](#overview)
- [Installation](#installation)
- [Documentation](#documentation)
- [Usage](#usage)

## Overview
Regular updates are made to the resources used for annotating variants, and as a result, some variants previously classified as Variants of Uncertain Significance (VUS) may be re-annotated as Pathogenic. Therefore, it is crucial to use up-to-date databases for functional annotation. These utilities offer tools that enable automatic updates of funcotator and annovar databases.

## Installation

### Option 1: Using Docker/Singularity (Recommended)

All required software is pre-installed in the Docker image. See [README_DOCKER.md](README_DOCKER.md) for detailed instructions.

```bash
# Using Docker
docker pull yinxiu/variantalker_db:latest

# Using Singularity (HPC)
singularity pull variantalker_db.sif docker://yinxiu/variantalker_db:latest
```

### Option 2: Manual Installation

#### Required Software

1. **Python 3.9+** with the following libraries:
```bash
conda create -n update python=3.9 -y
conda activate update
pip install beautifulsoup4 numpy pandas requests lxml psycopg[binary]
```

2. **vt** - Variant normalization and decomposition tool
```bash
git clone --recursive https://github.com/atks/vt.git
cd vt
make
# Add vt binary to your PATH
```

3. **GATK 4.5.0.0** - Genome Analysis Toolkit
```bash
wget https://github.com/broadinstitute/gatk/releases/download/4.5.0.0/gatk-4.5.0.0.zip
unzip gatk-4.5.0.0.zip
# Add gatk-4.5.0.0/gatk to your PATH
# Requires Java 17
```

4. **htslib** - Provides tabix and bgzip
```bash
wget https://github.com/samtools/htslib/releases/download/1.19/htslib-1.19.tar.bz2
tar -xjf htslib-1.19.tar.bz2
cd htslib-1.19
./configure --prefix=/usr/local
make
make install
```

5. **Perl** - Required for Annovar indexing scripts (usually pre-installed)

6. **Java 17** - Required for GATK
```bash
# Ubuntu/Debian
apt-get install openjdk-17-jdk

# Conda
conda install -c conda-forge openjdk=17
```

7. **System utilities** - wget, curl, git (usually pre-installed)

#### Clone Repository
```bash
git clone git@github.com:zhanyinx/variantalker.git
cd variantalker/update_db
```


## Documentation

These utilities determine if updates are necessary for specific databases. For some databases, a warning will be generated if updates are required without actually updating them, while others will be automatically updated, replacing the current database version.

List of Funcotator databases:

- Achilles, cancer-gene-census, familial & simple-uniprot (from Oncotator; check only)
- dna_repair_genes (check only)
- oreganno (check only, UPDATE November 2023, oreganno.org is still down, skipped)
- **cosmic (DISABLED - update script not compatible with new COSMIC download format)**
- **gencode (DISABLED - the updated database is built correctly but Funcotator does not use it properly)**
- clinvar (update if a newer version is available)
- hgnc (update if a newer version is available)
- dbsnp (update if a newer version is available)
- acmg_rec (update if a newer version is available)

List of Annovar databases:

- refGene (ignored)
- ensGene (ignored)
- knownGene (ignored)
- esp6500siv2_all (ignored)
- 1000g2015aug_all (ignored)
- exac03 (ignored)
- dbscsnv11 (ignored)
- dbnsfp31a_interpro (ignored)
- gnomad_genome (ignored)
- rmsk (ignored)
- avsnp (update if a newer version is available from the Annovar website)
- dbnsfp (update if a newer version is available from the Annovar website)
- clinvar (update if a newer version is available)
- **cosmic (DISABLED - update script not compatible with new COSMIC download format)**
- icgc (update if a newer version is available from the Annovar website)

IMPORTANT for annovar databases: since these databases are used by CancerVar and InterVar, we need to update also the corresponding CancerVar and InterVar. If CancerVar.py and InterVar.py and the corresponding configs are not provided (RECOMMENDED CASE), the tool will automatically update the CancerVar.py and InterVar.py and the relative configs in the [resources](https://github.com/zhanyinx/variantalker/tree/main/resources)

## Usage

**Note:** COSMIC update is currently disabled in both Funcotator and Annovar update scripts due to changes in the COSMIC database download format and authentication method. The scripts will skip COSMIC updates and log the current version. Manual COSMIC updates are required until the update scripts are updated to support the new format. The Funcotator GENCODE update is likewise disabled (the built database is not consumed correctly by Funcotator); the current version is retained and logged.

### Using Docker/Singularity

See [README_DOCKER.md](README_DOCKER.md) for complete Docker and Singularity usage instructions.

### Manual Execution

Funcotator update

```bash
# Update all databases for both genome builds (default)
python3 update_db/scripts/update_funcotator.py \
  -sd /path/to/funcotator_dataSources.v1.7.20200521s \
  -gd /path/to/funcotator_dataSources.v1.7.20200521g \
  -b /path/to/backup

# Update only hg38 databases
python3 update_db/scripts/update_funcotator.py \
  -sd /path/to/funcotator_dataSources.v1.7.20200521s \
  -gd /path/to/funcotator_dataSources.v1.7.20200521g \
  -b /path/to/backup \
  --build hg38

# Update only hg19 databases
python3 update_db/scripts/update_funcotator.py \
  -sd /path/to/funcotator_dataSources.v1.7.20200521s \
  -gd /path/to/funcotator_dataSources.v1.7.20200521g \
  -b /path/to/backup \
  --build hg19
```

Note: `--cosmic_email` and `--cosmic_password` arguments are optional but currently not used as COSMIC update is disabled.

Annovar update

```bash
python3 update_db/scripts/update_annovar.py \
  --annovar_db_path /path/to/humandb \
  --annovar_download_script /path/to/annovar/annotate_variation.pl \
  -v /path/to/vt
```

Note: `--cosmic_email` and `--cosmic_password` arguments are optional but currently not used as COSMIC update is disabled.

**Required Tools:**
- vt software can be found here: [vt](https://github.com/atks/vt)
- GATK 4.5.0.0: [GATK releases](https://github.com/broadinstitute/gatk/releases/tag/4.5.0.0)
- Annovar: [ANNOVAR](https://annovar.openbioinformatics.org/)
