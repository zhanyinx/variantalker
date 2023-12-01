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

The biomarker module of variantalker consists in extracting clinically relevant biomarkers starting from WES and RNAseq. The pipeline, at the moment, supports only Illumina Dragen outputs. This module relies on variantalker annotation output. We assume you have already run successfully the annotation.


## Usage

Modify the configuration file (nextflow.config) by setting the following parameters:

- ascat_genome_basedir: e.g. path2/public_databases

To perform biomarker analysis:

```bash
nextflow run path_to/main.nf -c yourconfig -profile singularity --input samplesheet.csv --outdir outdir --analysis biomarkers
```

Add --clonal_tmb_input samplesheet.clonaltmb.csv (see [format](https://github.com/zhanyinx/clonal_evolution#input)) to perform clonal tmb analysis

To show the whole list of parameters:

```bash
nextflow run path_to/main.nf --help --show_hidden_params
```

## Input

variantalker takes as input a csv samplesheet with 4 columns


__IMPORTANT: HEADER is required__ 

| patient        | sample_file       | sample_type  |
| -------------- | ----------------- | -------------|
| patient1       | path/tumor.vcf.gz | somatic      |
| .....          | .....             | .....        |

The input file is a csv with 3 columns: patient,samplefile,sample_type. Header is required!

The available sample_type are:

- variant_germline (maf file from variantalker)

- variant_somatic (maf file from variantalker)

- rna (dragen illumina rna pipeline, sf output file)

- msi (dragen msi output)

- tmb (dragen tmb output)

- cnv (output from variantalker)

If clonal tmb biomarker calculation is also required, the --clonal_tmb_input parameter must be also specified.
The format of the clonal_tmb_input file can be found [here](https://github.com/zhanyinx/clonal_evolution#input)

```bash
nextflow run path_to/main.nf run -with-tower -c nextflow.config  -profile conda --input sample.csv --outdir variantalker_output/ --analysis biomarkers --clonal_tmb_input sample_clonal_tmb.csv
```

## Output
Output structure:

```
params.outdir
|-- date
|   `-- biomarkers
|       |-- patient
|       |       |-- patient.report.html
|       |       |-- raw files with data for report
```

# Liability

Variantalker assumes no responsibility for any injury to person or damage to persons or property arising out of, or related to any use of Variantalker, or for any errors or omissions. The user recognizes they are using Liability at their own risk.
