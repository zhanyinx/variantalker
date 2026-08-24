[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![Active Development](https://img.shields.io/badge/Maintenance%20Level-Actively%20Developed-brightgreen.svg)](https://gist.github.com/cheerfulstoic/d107229326a01ff0f333a1d3476e068d)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

# The variantalker biomarker module

## Contents
- [Overview](#overview)
- [Usage](#usage)
- [Input](#input)
- [Output](#output)
- [Liability](#liability)

## Overview

The biomarker module within variantalker serves as a component dedicated to the extraction of clinically relevant biomarkers derived from Whole Exome Sequencing (WES) and RNA sequencing (RNAseq) data. It seamlessly integrates with the variantalker annotation output, assuming a prior successful execution of the annotation process.

This module is a **beta version** — that is what the pipeline's own `--analysis` help calls
it, and a biomarker run prints `Biomarkers BETA version!` as it starts.

**What needs Dragen, and what does not.** Two of the biomarkers are produced only when
`params.pipeline` is `dragen` (the default; `sarek` and `iontorrent` are the alternatives):
gene expression, and clonal TMB. TMB, the mutational signatures and the patient report
itself are computed whatever `params.pipeline` says. What the module cannot do without is
the *annotation* output, since the MAFs it reads come from a previous variantalker run.

**Mutational signatures** are COSMIC single-base-substitution activities, estimated with
SigProfilerAssignment against COSMIC `v3.4` (`--cosmic_version`), then collapsed into named
groups by `--cosmic_group`. The default group file, `resources/cosmic_sbs_group.csv`, defines
**14** groups — among them APOBEC, tobacco and UV, but also MMR, POL, HR and BER deficiency,
chemotherapy, treatment and artifact signatures. Point `--cosmic_group` at your own two-column
file (group name, SBS number) to group them differently.

## Usage

The ASCAT reference files used by clonal TMB are found under `params.database_dir`, at
`<database_dir>/<build>/ascat_wes_files/<build>/`. Leave the four `ascat_*` parameters
`null` — the default — and that is where the pipeline looks; a biomarker run logs the
four paths it chose as it starts, so check those lines rather than guessing. Override an
individual file with `--ascat_loci`, `--ascat_alleles`, `--ascat_loci_rt` or
`--ascat_loci_gc` if your layout differs.

To perform biomarker analysis:

> **Note:** on recent Nextflow releases the new strict language parser is the
> default and is not compatible with this pipeline. Prefix the commands below
> with `NXF_SYNTAX_PARSER=v1` (or `export NXF_SYNTAX_PARSER=v1` once). See the
> [main README](../../README.md#quick-start) for details.

```bash
NXF_SYNTAX_PARSER=v1 nextflow run path_to/main.nf -c yourconfig -profile singularity --input samplesheet.csv --outdir outdir --analysis biomarkers
```

Add --clonal_tmb_input samplesheet.clonaltmb.csv (see [format](https://github.com/zhanyinx/clonal_evolution#input)) to perform clonal tmb analysis

To show the whole list of parameters:

```bash
NXF_SYNTAX_PARSER=v1 nextflow run path_to/main.nf --help --show_hidden_params
```

## Input

variantalker takes as input a csv samplesheet with 3 columns


__IMPORTANT: HEADER is required__ 

| patient        | sample_file       | sample_type     |
| -------------- | ----------------- | ----------------|
| patient1       | path/tumor.maf    | variant_somatic |
| .....          | .....             | .....           |

The header must read exactly `patient,sample_file,sample_type` — all three names, spelled
that way. The pipeline checks the header before it does anything else and stops with
`Header missing or CSV file does not contain all of the required columns in the header:
[patient, sample_file, sample_type]`, so a missing underscore in `sample_file` is a
first-second-of-the-run failure rather than a silent one.

One patient has one row per sample_type. The available sample_type are:

- variant_germline (maf file from variantalker)

- variant_somatic (maf file from variantalker)

- variant_somatic_vcf (the somatic vcf the maf was annotated from)

- rna (dragen illumina rna pipeline, sf output file)

- msi (dragen msi output)

- tmb (dragen tmb output)

- hrd (homologous-recombination deficiency output)

- cnv (output from variantalker)

- coverage (per-sample coverage metrics; may repeat for one patient)

**`variant_somatic_vcf` is what TMB and the mutational signatures are computed from**, paired
with that patient's `variant_somatic` row. Omit it and nothing fails and nothing warns: the
report is still built, just with no TMB and no signatures in it. That is the one omission
here that is quiet.

Every sample_type is optional, and a report is produced from whichever rows are present.

If clonal tmb biomarker calculation is also required, the --clonal_tmb_input parameter must be also specified.
The format of the clonal_tmb_input file can be found [here](https://github.com/zhanyinx/clonal_evolution#input)

```bash
NXF_SYNTAX_PARSER=v1 nextflow run path_to/main.nf -c nextflow.config -profile singularity --input sample.csv --outdir variantalker_output/ --analysis biomarkers --clonal_tmb_input sample_clonal_tmb.csv
```

The only container profiles this pipeline defines are `singularity` and `docker` (plus
`debug`). Any other name stops the run immediately with `Unknown configuration profile`,
before a single process is submitted.

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

## Liability

Variantalker assumes no responsibility for any injury to person or damage to persons or property arising out of, or related to any use of Variantalker, or for any errors or omissions. The user recognizes they are using Variantalker at their own risk.
