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

Download the updated databases 

```bash
wget -r -N --no-parent -nH --cut-dirs=2 -P public_databases https://bioserver.ieo.it/repo/dima/ 
```

## Documentation

The pipeline employs two tools to annotate and prioritize variants: [Funcotator](https://gatk.broadinstitute.org/hc/en-us/articles/360035889931-Funcotator-Information-and-Tutorial) and [CancerVar](https://github.com/WGLab/CancerVar) for somatic variants and InterVar for germline variants, both of which rely on [Annovar](https://annovar.openbioinformatics.org/en/latest/). To ensure the accuracy of the pipeline, the databases for Funcotator and Annovar must be regularly updated using the provided tools found here: [update utilities](https://github.com/zhanyinx/variantalker/tree/main/update_db).


## Usage

If you are using for the first time, update the databases following the [instructions](https://github.com/zhanyinx/variantalker/tree/main/update_db). 

Update in the configuration file (nextflow.config) by setting the path to the databases:

- funcotator_germline_db: e.g. path2/public_databases/funcotator_dataSources.v1.7.20200521g

- funcotator_somatic_db: e.g. path2/public_databases/funcotator_dataSources.v1.7.20200521s

- annovar_db: e.g. path2/public_databases/humandb

- annovar_software_folder: e.g. path2/annovar


```bash
nextflow run path_to/main.nf -c yourconfig -profile singularity --input samplesheet.csv --output outdir
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

- somatic sample type: it can be tumor_only (single sample) or tumor_normal (multi sample) vcf.gz file. Requires tumor_tissue to be specified

- germline: single sample vcf.gz file. It does not require tumor_tissue

- cnv: for nfcore/sarek, CNVKit output is supported (cnr file). For dragen, vcf.gz file required. It does not require tumor_tissue 

Available tumor_tissue are: Adrenal_Gland Bile_Duct Bladder Blood Bone Bone_Marrow Brain Breast Cancer_all Cervix Colorectal Esophagus Eye Head_and_Neck Inflammatory Intrahepatic Kidney Liver Lung Lymph_Nodes Nervous_System Other Ovary Pancreas Pleura Prostate Skin Soft_Tissue Stomach Testis Thymus Thyroid Uterus

__BETA version biomarkers__

To extract biomarkers, you can use the same type of input sheet. There are two types of sample files available: the maf format, which is the output of the annotation analysis, and the .sf file from the Illumina Dragen RNA pipeline. These sample files are categorized based on their sample type: dna for maf files and rna for .sf files.

__BETA version clonal analysis__

Clonal analysis works only with conda enviroment for now. To extract clonal TMB, we utilize the [nfcore/sarek](https://nf-co.re/sarek)'s ascat tool  and [pyclone-vi](https://github.com/Roth-Lab/pyclone-vi). The input file format that is accepted is the same as in nfcore/sarek, but it includes twi additional columns: 

1) cellularity  

2) annotated maf file from the tumor sample for which you want to calculate the clonal tmb.

To run clonal tmb, add 

```
process {
   withName: 'clonal_tmb' {
      conda = 'PATH2/variantalker/resources/envs/clonal_tmb.yaml'
   }
}
```

to the configuration file

Example code run:

```bash
nextflow run path_to/main.nf run -with-tower -c nextflow.config  -profile conda --input sample.csv --output variantalker_output/ --analysis clonal_tmb
```


## Output

Output structure:

```
params.output
|-- date
|   `-- annotation
|       |-- germline
|       |   `-- patient
|       |       |-- filtered.patient.small_mutations.intervar.escat.renovo.maf.tsv
|       |       |-- patient.cnv.annotated.tsv
|       |       `-- patient.small_mutations.intervar.escat.renovo.maf
|       `-- somatic
|           `-- patient
|               |-- filtered.patient.small_mutations.cancervar.escat.maf.tsv
|       |       |-- patient.cnv.annotated.tsv
|               `-- patient.small_mutations.cancervar.escat.maf
|   `-- biomarkers
|       |-- patient
|       |       |-- patient.rna.tpm.csv
|       |       |-- tmb_signatures.patient.txt
```

variantalker outputs for each sample two files

1) *maf file with all the annotations
2) filtered*tsv file with filtered variants.

Default filters applied:

- "Silent", "Intron", "3'UTR", "5'UTR", "IGR", "5'Flank", "3'Flank", "RNA" variant types are filtered out

-  pathogenic, likely pathogenic, uncertain in any of InterVar, CancerVar or clinvar are kept. For somatic samples also variants with ESCAT tier I and II are kept.

- variants with variant allele frequency smaller than 0.01 (somatic) and 0.2 (germline) are filtered out
