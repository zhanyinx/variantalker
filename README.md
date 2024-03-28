[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![Active Development](https://img.shields.io/badge/Maintenance%20Level-Actively%20Developed-brightgreen.svg)](https://gist.github.com/cheerfulstoic/d107229326a01ff0f333a1d3476e068d)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

# Variant annotation and prioritization pipeline

## Contents
- [Contents](#contents)
- [Overview](#overview)
- [Installation](#installation)
- [Documentation](#documentation)
- [Usage](#usage)
- [Input](#input)
- [Output](#output)
- [Liability](#liability)

## Overview

Variant annotation in cancer genomics involves identifying and characterizing the genetic changes (variants) that contribute to cancer development and progression. The challenge is that there are many different types of variants that can occur in the genome, and not all of them are relevant to cancer. Therefore, accurate annotation is critical for identifying the key driver mutations and designing targeted therapies. However, this process is complicated by the large number of potential variants, the need to integrate data from multiple sources, and the ongoing discovery of new cancer-associated variants.

We have developed a Nextflow pipeline called variantalker that enables users to annotate variants from VCF files. Our pipeline supports VCF files generated from dragen, nf-sarek, and ION-torrent platforms.

BETA version: we have implemented the possibility to extract biomarkers such as TMB, mutational signatures (apobec, uv and tabacco), clonal TMB (if bam/cram files and sex are provided), expression of specific genes (if RNA-seq data are provided), gene cnv, etc. For more information, look at [here](https://github.com/zhanyinx/variantalker/tree/main/docs/biomarkers/)

## Installation
Clone the repo

```bash
git clone https://github.com/zhanyinx/variantalker.git
```

variantalker relies on [Annovar](https://annovar.openbioinformatics.org/en/latest/) software and [Funcotator](https://gatk.broadinstitute.org/hc/en-us/articles/360035889931-Funcotator-Information-and-Tutorial) databases.

Download the updated databases. Separate repositories for hg19 and hg38 are available.

```bash
wget -r -N --no-parent -nH --cut-dirs=3 -P public_databases/hg38 https://bioserver.ieo.it/repo/dima/hg38 
wget -r -N --no-parent -nH --cut-dirs=3 -P public_databases/hg19 https://bioserver.ieo.it/repo/dima/hg19
```

## Documentation

The pipeline employs several tools to annotate and prioritize variants: 

- [Funcotator](https://gatk.broadinstitute.org/hc/en-us/articles/360035889931-Funcotator-Information-and-Tutorial) for variant annotation
- [CancerVar](https://github.com/WGLab/CancerVar) for somatic variants prioritization
- [InterVar](https://github.com/WGLab/InterVar) for germline variants annotation
- [Annovar](https://annovar.openbioinformatics.org/en/latest/): cancervar and intervar reply on Annovar. 
- [CIViC](https://civicdb.org/): somatic variant classification using CIViC evidence level.
- [AlphaMissense](https://www.science.org/doi/10.1126/science.adg7492): somatic and germline variant prioritization.

To ensure the accuracy of the pipeline, the databases for Funcotator and Annovar must be regularly updated using the provided tools found here: [update utilities](https://github.com/zhanyinx/variantalker/tree/main/update_db).


## Usage

If you are using for the first time, please consider updating the databases following the [instructions](https://github.com/zhanyinx/variantalker/tree/main/update_db). 

Modify the configuration file (nextflow.config) by setting the following parameters:

- funcotator_germline_db: e.g. path2/public_databases/funcotator_dataSources.v1.7.20200521g

- funcotator_somatic_db: e.g. path2/public_databases/funcotator_dataSources.v1.7.20200521s

- annovar_db: e.g. path2/public_databases/humandb

- annovar_software_folder: e.g. path2/annovar

- alpha_mis_genome_basedir: e.g. path2/public_databases

- fasta: path to fasta file used to generate the vcf

- target: path to the target bed file

The main command line for the annotation is the following

```bash
nextflow run path_to/main.nf -c yourconfig -profile singularity --input samplesheet.csv --outdir outdir
```

```bash
nextflow run path_to/main.nf --help --show_hidden_params
```

## Input

variantalker takes as input a csv samplesheet with 4 columns



__IMPORTANT: HEADER is required__ 

| patient        | tumor_tissue   | sample_file       | sample_type  |
| -------------- | -------------- | ----------------- | -------------|
| patient1       | Lung           | path/tumor.vcf.gz | somatic      |
| .....          | .....          | .....             | .....        |

Sample_file must be provided with full path, __not__ relative path

Available sample_type are: somatic, germline, cnv. 

- somatic sample type: it can be tumoronly (single sample) or tumor_normal (multi sample) vcf.gz file. Requires tumor_tissue to be specified

- germline: single sample vcf.gz file. It does not require tumor_tissue

- cnv: for nfcore/sarek, CNVKit output is supported (cnr file). For dragen, vcf.gz file required. It does not require tumor_tissue 

Available tumor_tissue are: Adrenal_Gland Bile_Duct Bladder Blood Bone Bone_Marrow Brain Breast Cancer_all Cervix Colorectal Esophagus Eye Head_and_Neck Inflammatory Intrahepatic Kidney Liver Lung Lymph_Nodes Nervous_System Other Ovary Pancreas Pleura Prostate Skin Soft_Tissue Stomach Testis Thymus Thyroid Uterus

## Output

Output structure:

```
params.outdir
|-- date
|   `-- annotation
|       |-- germline
|       |   `-- patient
|       |       |-- filtered.patient.maf.pass.tsv
|       |       |-- filtered.patient.maf.nopass.tsv
|       |       |-- patient.vcf
|       |       `-- patient.maf
|       `-- somatic
|           `-- patient
|       |       |-- filtered.patient.maf.pass.tsv
|       |       |-- filtered.patient.maf.nopass.tsv
|       |       |-- patient.vcf
|               `-- patient.maf
|       `-- cnv
|           `-- patient
|       |       |-- patient.cnv.annotated.tsv
```

variantalker outputs for each sample multiple files

1) maf file with all the annotations
2) vcf file with the PASS variants 
3) filtered pass file with variants passing the filters (see below).
4) filtered nopass file with variants not passing the filters (see below).
5) cnv annotated file (if cnv samples provided)

Default filters applied:

- "Silent",  "IGR",  "RNA" variant types are filtered out (unless it's pathogenic or likely pathogenic for clinvar/cancervar/intervar)

- minimum coverage 50 (unless it's pathogenic or likely pathogenic for clinvar/cancervar/intervar)

- minimum somatic VAF: 0.01

- minimum germline VAF: 0.2

- InterVar classes to be kept: Pathogenic,Likely pathogenic (logic OR)

- CancerVar classes to be kept: Tier_II_potential,Tier_I_strong (logic OR)

- ReNOVo class to be kept: LP Pathogenic,IP Pathogenic,HP Pathogenic (logic OR)

- CIViC evidence levels to be kept: A,B,C (logic OR)

- no filters on genes (somatic or germline)

Logic OR filters: a variant is kept if at least one of the OR filters is true


# Liability

Variantalker assumes no responsibility for any injury to person or damage to persons or property arising out of, or related to any use of Variantalker, or for any errors or omissions. The user recognizes they are using Liability at their own risk.
