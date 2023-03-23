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
The following libraries are required (installation using pip install libraryname)

```bash
conda create -n update python=3.9 -y
conda activate update
pip install bs4 numpy pandas requests lxml

git clone git@github.com:zhanyinx/variantalker.git
```


## Documentation

These utilities determine if updates are necessary for specific databases. For some databases, a warning will be generated if updates are required without actually updating them, while others will be automatically updated, replacing the current database version.

List of Funcotator databases:

- Achilles, cancer-gene-census, familial & simple-uniprot (from Oncotator; check only)
- dna_repair_genes (check only)
- oreganno (check only)
- cosmic (update if a newer version is available)
- gencode (update if a newer version is available, currently not functioning correctly)
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
- cosmic (update if a newer version is available)
- icgc (update if a newer version is available from the Annovar website)

IMPORTANT for annovar databases: since these databases are used by CancerVar and InterVar, we need to update also the corresponding CancerVar and InterVar. If CancerVar.py and InterVar.py and the corresponding configs are not provided (RECOMMENDED CASE), the tool will automatically update the CancerVar.py and InterVar.py and the relative configs in the [resources](https://github.com/zhanyinx/variantalker/tree/main/resources)

## Usage

Funcotator update

```bash
path_to/update_db/scripts/update_funcotator.py -sd pathto/funcotator_dataSources.v1.7.20200521s -ce your_cosmic_email -cp your_cosmic_password -gd pathto/funcotator_dataSources.v1.7.20200521g -b backup
```

Annovar update

```bash
pathto/update_db/scripts/update_annovar.py --annovar_db_path pathto/humandb --annovar_download_script pathto/annotate_variation.pl -v pathto/vt/ -ce your_cosmic_email         -cp your_cosmic_password
```

where vt software can be found here [vt](https://github.com/atks/vt)
