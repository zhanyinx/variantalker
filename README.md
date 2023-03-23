[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![Active Development](https://img.shields.io/badge/Maintenance%20Level-Actively%20Developed-brightgreen.svg)](https://gist.github.com/cheerfulstoic/d107229326a01ff0f333a1d3476e068d)

# Variant annotation and prioritization pipeline

## Contents
- [Contents](#contents)
- [Overview](#overview)
- [Installation](#installation)
- [Documentation](#documentation)
- [Usage](#usage)
- [Input](#input)

## Overview
Variant annotation in cancer genomics involves identifying and characterizing the genetic changes (variants) that contribute to cancer development and progression. The challenge is that there are many different types of variants that can occur in the genome, and not all of them are relevant to cancer. Therefore, accurate annotation is critical for identifying the key driver mutations and designing targeted therapies. However, this process is complicated by the large number of potential variants, the need to integrate data from multiple sources, and the ongoing discovery of new cancer-associated variants.

We have developed a Nextflow pipeline called variantalker that enables users to annotate variants from VCF files. Our pipeline supports VCF files generated from dragen, nf-sarek, and ION-torrent platforms.

## Installation
Clone the repo

```bash
git clone git@github.com:zhanyinx/variantalker.git
```

Update the databases following the [instructions](https://github.com/zhanyinx/variantalker/tree/main/update_db)

## Documentation
The pipeline employs two tools to annotate and prioritize variants: [Funcotator](https://gatk.broadinstitute.org/hc/en-us/articles/360035889931-Funcotator-Information-and-Tutorial) and [CancerVar](https://github.com/WGLab/CancerVar) for somatic variants and InterVar for germline variants, both of which rely on [Annovar](https://annovar.openbioinformatics.org/en/latest/). To ensure the accuracy of the pipeline, the databases for Funcotator and Annovar must be regularly updated using the provided tools found here: [update utilities](https://github.com/zhanyinx/variantalker/tree/main/update_db).


## Usage

First of all, set the configuration file (nextflow.config and the config/configs files) accordingly to the databases you have.

```bash
nextflow run path_to/main.nf -c yourconfig
```

## Input

In the config file, depending on which source of data you want to analyse, you have to provide different type of input data.

### Somatic small variant (snp, indel)

```bash
nextflow run path_to/main.nf -c yourconfig --somatic.input_snp_indel samples.tsv
```

Two columns: first column is cancertype, second column is the vcf.gz file.

Available cancertype: Adrenal_Gland Bile_Duct Bladder Blood Bone Bone_Marrow Brain Breast Cancer_all Cervix Colorectal Esophagus Eye Head_and_Neck Inflammatory Intrahepatic Kidney Liver Lung Lymph_Nodes Nervous_System Other Ovary Pancreas Pleura Prostate Skin Soft_Tissue Stomach Testis Thymus Thyroid Uterus

__IMPORTANT: DO NOT INCLUDE THE HEADER__ 

| tumor type     | path to vcf.gz    |
| -------------- | ----------------- |
| Lung           | path/tumor.vcf.gz |
| .....          | .....             | 

### Germline small variant

```bash
nextflow run path_to/main.nf -c yourconfig --germline.input_snp_indel path/to/*/multiple/*vcf.gz
```

