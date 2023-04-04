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

## Overview

Variant annotation in cancer genomics involves identifying and characterizing the genetic changes (variants) that contribute to cancer development and progression. The challenge is that there are many different types of variants that can occur in the genome, and not all of them are relevant to cancer. Therefore, accurate annotation is critical for identifying the key driver mutations and designing targeted therapies. However, this process is complicated by the large number of potential variants, the need to integrate data from multiple sources, and the ongoing discovery of new cancer-associated variants.

We have developed a Nextflow pipeline called variantalker that enables users to annotate variants from VCF files. Our pipeline supports VCF files generated from dragen, nf-sarek, and ION-torrent platforms.

BETA version: we have implemented the possibility to extract biomarkers such as TMB, mutational signatures (apobec, uv and tabacco), clonal TMB and expression of specific genes if RNA-seq is given

## Installation
Clone the repo

```bash
git clone git@github.com:zhanyinx/variantalker.git
```

variantalker relies on [Annovar](https://annovar.openbioinformatics.org/en/latest/) software and [Funcotator](https://gatk.broadinstitute.org/hc/en-us/articles/360035889931-Funcotator-Information-and-Tutorial) databases.

1- get the [Funcotator](https://gatk.broadinstitute.org/hc/en-us/articles/360035889931-Funcotator-Information-and-Tutorial) databases and [Annovar](https://annovar.openbioinformatics.org/en/latest/)

2- Update the databases following the [instructions](https://github.com/zhanyinx/variantalker/tree/main/update_db). 

Work in progress: download updated databases from [missing link]()

## Documentation

The pipeline employs two tools to annotate and prioritize variants: [Funcotator](https://gatk.broadinstitute.org/hc/en-us/articles/360035889931-Funcotator-Information-and-Tutorial) and [CancerVar](https://github.com/WGLab/CancerVar) for somatic variants and InterVar for germline variants, both of which rely on [Annovar](https://annovar.openbioinformatics.org/en/latest/). To ensure the accuracy of the pipeline, the databases for Funcotator and Annovar must be regularly updated using the provided tools found here: [update utilities](https://github.com/zhanyinx/variantalker/tree/main/update_db).


## Usage

If you are using for the first time, update the databases following the [instructions](https://github.com/zhanyinx/variantalker/tree/main/update_db). 

Update in the configuration file (nextflow.config) the path to the databases:

- funcotator germline: e.g. path2/funcotator_dataSources.v1.7.20200521g

- funcotator somatic: e.g. path2/funcotator_dataSources.v1.7.20200521s

- annovar databases: e.g. path2/humandb

- annovar software folder: e.g. path2/annovar


```bash
nextflow run path_to/main.nf -c yourconfig -profile singularity --input samplesheet.csv --output outdir
```

## Input

variantalker takes as input a csv samplesheet with 3 columns



__IMPORTANT: DO NOT INCLUDE THE HEADER__ 

| tumor_tissue   | sample_file       | sample_type  |
| -------------- | ----------------- | -------------|
| Lung           | path/tumor.vcf.gz | somatic      |
| .....          | .....             | .....        |

Available sample_type are: somatic, germline, cnv. 

- somatic sample type: it can be tumor_only (single sample) or tumor_normal (multi sample) vcf.gz file. Requires tumor_tissue to be specified

- germline: single sample vcf.gz file. It does not require tumor_tissue

- cnv: for nfcore/sarek, CNVKit output is supported (cnr file). For dragen, vcf.gz file required. It does not require tumor_tissue 

Available tumor_tissue are: Adrenal_Gland Bile_Duct Bladder Blood Bone Bone_Marrow Brain Breast Cancer_all Cervix Colorectal Esophagus Eye Head_and_Neck Inflammatory Intrahepatic Kidney Liver Lung Lymph_Nodes Nervous_System Other Ovary Pancreas Pleura Prostate Skin Soft_Tissue Stomach Testis Thymus Thyroid Uterus

## Output

Output structure:

```
params.output
|-- date
|   `-- annotation
|       |-- germline
|       |   `-- sampleid
|       |       |-- filtered.sampleid.small_mutations.intervar.escat.renovo.maf.tsv
|       |       |-- sampleid.cnv.annotated.tsv
|       |       `-- sampleid.small_mutations.intervar.escat.renovo.maf
|       `-- somatic
|           `-- sampleid
|               |-- filtered.sampleid.small_mutations.cancervar.escat.maf.tsv
|       |       |-- sampleid.cnv.annotated.tsv
|               `-- sampleid.small_mutations.cancervar.escat.maf
```

variantalker outputs for each sample two files

1) *maf file with all the annotations
2) filtered*tsv file with filtered variants.

Filters applied:

- "Silent", "Intron", "3'UTR", "5'UTR", "IGR", "5'Flank", "3'Flank", "RNA" variant types are filtered out

- only variants ESCAT tier I and II or pathogenic, likely pathogenic, uncertain in any of InterVar, CancerVar or clinvar are kept

- variants with variant allele frequency smaller than 0.05 (somatic) and 0.2 (germline) are filtered out
