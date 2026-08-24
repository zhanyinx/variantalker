[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![Active Development](https://img.shields.io/badge/Maintenance%20Level-Actively%20Developed-brightgreen.svg)](https://gist.github.com/cheerfulstoic/d107229326a01ff0f333a1d3476e068d)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

# VarianTalker

A Nextflow pipeline for cancer variant annotation, prioritization, and biomarker extraction.

## Contents
- [Overview](#overview)
- [Features](#features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Input Format](#input-format)
- [Output](#output)
- [Modules](#modules)
- [Citations](#citations)
- [Contributing](#contributing)
- [License](#license)

## Overview

VariantAlker is a comprehensive pipeline for annotating and prioritizing genetic variants in cancer genomics. It integrates multiple annotation databases and tools to identify clinically relevant mutations and extract actionable biomarkers.

**Supported platforms:** Dragen, nf-sarek, ION-Torrent  
**Supported genomes:** hg19, hg38

## Features

### Core Annotation
- **Somatic variants** - CancerVar prioritization, ESCAT tiering and CIViC evidence levels
- **Germline variants** - InterVar classification and ReNOVo scoring
- **CNV annotation** - CNVKit (nf-sarek) and Dragen outputs
- **Multi-database integration** - Funcotator (which brings ClinVar and gnomAD), Annovar,
  CIViC, AlphaMissense
- **Pharmacogenomics** - a PharmCAT report per germline sample, on by default for hg38
  (see [Output](#output))

### Biomarker Extraction (Beta)
- Tumor Mutational Burden (TMB)
- Mutational signatures — COSMIC SBS activities collapsed into groups; the default grouping
  file, `resources/cosmic_sbs_group.csv`, defines 14 of them, among them APOBEC, tobacco and UV
- Clonal TMB (requires BAM/CRAM files)
- Gene expression (requires RNA-seq data)
- Gene CNV

Biomarkers are a separate entry point — `--analysis biomarkers`, with its own samplesheet.
See [biomarker documentation](docs/biomarkers/README.md) for details.

### Data Visualization
Filter and browse the annotated MAF interactively with the **MAFigate** app (see [`streamlit_app/`](streamlit_app/README.md)).

## Installation

### 0. What you need first

| requirement | notes |
| --- | --- |
| **Nextflow** ≥ 22.10.1 | Recent releases need `NXF_SYNTAX_PARSER=v1` — see [Quick Start](#quick-start). Checked here on **26.04.3**. |
| **Java** | Nextflow's own requirement; checked here on Java 18. |
| **Singularity/Apptainer, or Docker** | Every process runs in a container, and there is no local-software route. Pick the runtime your machine has and pass the matching `-profile`. |
| **Disk** | The annotation databases below are large: **at least 624 GiB** for the whole tree — 478 GiB of Annovar databases shared by both genome builds, plus 67 GiB (hg19) or 79 GiB (hg38) for the build you use. See [step 2](#2-download-annotation-databases). |

### 1. Clone the repository
```bash
git clone https://github.com/zhanyinx/variantalker.git
cd variantalker
```

### 2. Download annotation databases
The databases are served as a public directory tree at
[`https://repo.bioserver.ieo.it/dima/`](https://repo.bioserver.ieo.it/dima/) — anonymous HTTPS,
no credentials, reachable from outside the IEO network. It holds one directory per genome build
(`hg19/`, `hg38/`) plus the shared Annovar databases (`humandb/`), which is exactly the layout
step 3 configures, so a single recursive download puts everything where the pipeline looks:

```bash
wget -r -N --no-parent -nH --cut-dirs=1 -P path/to/public_databases https://repo.bioserver.ieo.it/dima/
```

That lands `path/to/public_databases/{hg19,hg38,humandb}`, and step 3 then needs one parameter.

Three details, each of which fails quietly if you get it wrong:

- **Keep the trailing slash on the URL.** Without it the server redirects the directory to the
  `http://` spelling, and port 80 answers **404** — which looks exactly like a wrong address.
- **`--cut-dirs` counts the directory levels in the URL, so it changes with the URL.** It is `1`
  for `/dima/` above and `2` for the per-genome URLs below. Too high a value strips a level you
  need: at `3`, `alpha_missense`, `ascat_wes_files`, `funcotator_dataSources.*` and `gwas` all
  merge into a single directory, and the download still reports success.
- **This is a very large download** — see the disk requirement in
  [step 0](#0-what-you-need-first).

If you only want one genome build, fetch it and `humandb/` separately; `humandb/` is needed
either way, and it is the bulk of the tree:

```bash
wget -r -N --no-parent -nH --cut-dirs=2 -P path/to/public_databases/hg38    https://repo.bioserver.ieo.it/dima/hg38/
wget -r -N --no-parent -nH --cut-dirs=2 -P path/to/public_databases/humandb https://repo.bioserver.ieo.it/dima/humandb/
```

For hg19 instead, swap `hg38` for `hg19` in both the destination and the URL of the first line.

The **Annovar software** is not in this tree; only the Annovar *databases* are. See
[step 3](#3-configure-the-pipeline).

> **What has and has not been checked here.** The addresses, the tree's layout and the sizes
> above were read off the live server; each `--cut-dirs` value was checked by downloading one
> tiny file from the tree and looking at where it landed, rather than by counting levels on
> paper. No database content was fetched, so nothing here says the files themselves are complete
> or current — only that these commands put them where the pipeline expects them.

### 3. Configure the pipeline
The pipeline derives every database path from a single parameter, `database_dir`, and the
layout it expects is the layout step 2 downloads into: a directory holding `hg38/`, `hg19/`
and `humandb/`. Setting that one parameter is usually the whole configuration.

Edit [nextflow.config](nextflow.config):

```groovy
params {
  database_dir = "path/to/public_databases"   // holds hg38/, hg19/ and humandb/
  fasta        = "path/to/reference.fasta"    // must match --build
  target       = "path/to/target.bed"         // exome capture regions
}
```

Each database is then looked up as follows, and each has its own parameter if your copy sits
somewhere else:

| database | default location | parameter to override it |
| --- | --- | --- |
| Funcotator, somatic | `<database_dir>/<build>/funcotator_dataSources.v1.8.2024s` | `funcotator_somatic_db` |
| Funcotator, germline | `<database_dir>/<build>/funcotator_dataSources.v1.8.2024g` | `funcotator_germline_db` |
| Annovar databases | `<database_dir>/humandb` | `annovar_db` |
| Annovar software | `<database_dir>/annovar_software` | `annovar_software_folder` |
| AlphaMissense | `<database_dir>/<build>/alpha_missense/<build>/AlphaMissense_<build>.231106.tsv.gz` | `alpha_missense` |

**The Annovar *software* is not part of step 2's download.** Annovar is not redistributable and
is not vendored here — register at
[ANNOVAR](https://annovar.openbioinformatics.org/en/latest/user-guide/download/), unpack it
wherever you like, and point `annovar_software_folder` at that directory. Step 2 downloads the
Annovar *databases* (`humandb`), not the tool that reads them.

`<build>` is the value of `--build`, which defaults to `hg38`. The default `database_dir` is
`resources/databases` inside the checkout, which does not exist in a fresh clone — so a run
that has not been configured looks for its databases there and fails on the first process that
needs one.

Two things worth knowing before the first run:

- **The run prints the paths it settled on.** Every one of the five lines above is echoed as
  `Funcotator somatic database chosen: …`, `Annovar database chosen: …` and so on, in the
  first screen of output. That listing is the fastest way to see whether `database_dir` took
  effect; it is printed whether or not the paths exist.
- **`fasta` is not checked up front.** An unset or empty `fasta` is skipped by the pipeline's
  start-up path check rather than reported, and surfaces later as a failure inside the
  variant-standardization process.

## Quick Start

> **Nextflow version / language parser**
> Recent Nextflow releases ship a new, stricter language parser as the default,
> which is **not** compatible with this pipeline's current syntax. Run with the
> legacy parser by setting `NXF_SYNTAX_PARSER=v1` (shown below), or use an older
> Nextflow release (≥ 22.10.1). All commands in this README assume the legacy
> parser. Without it, config parsing fails at `nextflow.config` before the
> pipeline starts — so the badge's `≥22.10.1` is a floor, not a promise that any
> newer release runs unaided.

### Basic annotation
```bash
NXF_SYNTAX_PARSER=v1 nextflow run main.nf -profile singularity \
  --input samplesheet.csv --build hg38 --outdir results
```

`--build` selects both the genome and the database subdirectory, and it defaults to `hg38`. It
is spelled out above because omitting it on an hg19 project is silent: the run succeeds against
the hg38 databases. Only `hg19` and `hg38` are accepted; anything else stops the run.

`-profile docker` is the equally supported alternative and takes the same command. Pick the
one your machine has — without a profile, Nextflow warns that it is running with no container
configuration and the processes will not find their tools.

> **What was verified here, and what was not.** The command above was not run to completion:
> that needs the databases from step 2. `-profile singularity` could not be exercised at all on
> the machine this page was checked from — neither singularity nor apptainer is installed there
> — while `-profile docker` was available. Both profile names resolve, and both are defined in
> [nextflow.config](nextflow.config); `conda` is not a profile of this pipeline.

### View all options
```bash
NXF_SYNTAX_PARSER=v1 nextflow run main.nf --help --show_hidden_params
```

### Biomarkers
```bash
NXF_SYNTAX_PARSER=v1 nextflow run main.nf -profile singularity --analysis biomarkers \
  --input biomarkers_samplesheet.csv --build hg38 --outdir results
```

Biomarkers take a different samplesheet from the annotation run and are still beta; the
[biomarker documentation](docs/biomarkers/README.md) is the page to follow for them.

> **Tip:** to avoid prefixing every command, export it once in your shell:
> ```bash
> export NXF_SYNTAX_PARSER=v1
> ```


## Input Format

Create a CSV samplesheet with the following columns (**header required**):

| patient  | tumor_tissue | sample_file              | sample_type |
|----------|--------------|--------------------------|-------------|
| patient1 | Lung         | /path/to/tumor.vcf.gz    | somatic     |
| patient2 | Breast       | /path/to/germline.vcf.gz | germline    |
| patient3 | Lung         | /path/to/cnv.cnr         | cnv         |

**Important:** Use absolute paths for `sample_file`, not relative paths.

All four column names must be present in the header, including `tumor_tissue` on a sheet with
no somatic rows — a missing one stops the run naming the column it wanted. An optional fifth
column, `custom_id`, is accepted: where it is set it replaces the patient name as the output
directory for that sample (see [Output](#output)), and two rows sharing a `custom_id` stop the
run.

**Only records marked `PASS` in the VCF's FILTER column are annotated.** Everything else is
dropped in the first step, before any database is consulted, so a caller that leaves FILTER
empty or `.` yields a sample with no variants — reported as header-only output files with a
warning, not as a failure.

### Sample Types

- **somatic** - Single-sample (tumor-only) or multi-sample (tumor-normal) VCF. Requires `tumor_tissue`.
- **germline** - Single-sample VCF. `tumor_tissue` is ignored; the pipeline sets it to `Germline` itself.
- **cnv** - CNVKit `.cnr` file (nf-sarek) or VCF (Dragen). `tumor_tissue` is not read.

`sample_type` is matched exactly and case-sensitively against these three words. A row whose
`sample_type` is anything else — `Somatic`, `tumour`, a stray space — is **silently dropped**:
the run reports `Pipeline completed successfully` having done nothing at all. If a run finishes
with no tasks, this column is the first thing to check.

`cnv` rows are annotated only under `--pipeline sarek` or `--pipeline dragen`. Under
`--pipeline iontorrent` they are read from the samplesheet and then go nowhere.

### Tumor Tissue Options
Adrenal_Gland, Bile_Duct, Bladder, Blood, Bone, Bone_Marrow, Brain, Breast, Cancer_all, Cervix, Colorectal, Esophagus, Eye, Head_and_Neck, Inflammatory, Intrahepatic, Kidney, Liver, Lung, Lymph_Nodes, Nervous_System, Other, Ovary, Pancreas, Pleura, Prostate, Skin, Soft_Tissue, Stomach, Testis, Thymus, Thyroid, Uterus

These 33 values are matched exactly, and an unrecognised one stops the run rather than falling
back to a default.

## Output

The pipeline generates the following directory structure:

```
results/
└── YYMMDD/                                   # the run date, e.g. 260820
    └── annotation/
        ├── somatic/
        │   └── patient_name/                 # custom_id instead, where the samplesheet sets one
        │       ├── patient.maf               # every annotated variant, with a PASS/NOPASS filter column
        │       ├── patient.vcf               # PASS variants only
        │       ├── patient.pass.tsv          # the variants that passed the filters
        │       └── patient.nopass.tsv        # the variants that did not
        ├── germline/
        │   └── patient_name/
        │       ├── patient.maf
        │       ├── patient.vcf
        │       ├── patient.pass.tsv
        │       ├── patient.nopass.tsv
        │       └── *.html                    # PharmCAT report, unless --skip_pharmgkb
        └── cnv/
            └── patient_name/
                └── patient.cnv.annotated.tsv
```

The date directory is `YYMMDD` — six digits, not an ISO date — and it is computed when the run
starts, so a run spanning midnight publishes under the day it began.

Two more things land under `results/YYMMDD/` that are not part of the per-sample annotation:

- `intermediate_files/cnvkit/cnv/` — the called CNV segments, on the nf-sarek route only.
- Under `--analysis biomarkers`, a `biomarkers/` tree instead of `annotation/`; see the
  [biomarker documentation](docs/biomarkers/README.md).

The PharmCAT report is produced for **germline** samples. On hg38 it runs by default; on hg19 it
needs `pharmgkb_hg38_fasta` and `pharmgkb_hg38_chain` for the liftover and prints a skip message
without them. `--skip_pharmgkb` turns it off. Its exact filename comes from PharmCAT rather than
from this pipeline and was not verified here.

### Default Filters

Every annotated variant carries a `filter` column of `PASS` or `NOPASS`, and the two `.tsv`
files are that split. A variant is `PASS` if **either** of the following holds.

**1. It clears every filter.** All three of these, together:

| filter | default | parameter |
| --- | --- | --- |
| Variant classification not excluded | excludes `Silent`, `IGR`, `RNA` | `filter_var_classification` |
| Depth at the site (`t_ref_count + t_alt_count`) | ≥ 50 | `filter_min_depth` |
| Variant allele frequency | **>** 0.01 somatic, **>** 0.2 germline (strictly greater) | `filter_vaf_threshold`, `filter_vaf_threshold_germline` |

…**and** at least one classification hit, which is where the OR sits — and the two arms have
different lists:

| arm | any one of | default | parameter |
| --- | --- | --- | --- |
| somatic | CancerVar | Tier_I_strong, Tier_II_potential | `filter_cancervar` |
| somatic | CIViC evidence level | A, B, C | `filter_civic_evidence_level` |
| somatic | ESCAT | IA, IB, IC, IIA, IIB | `filter_escat` |
| somatic | ClinVar significance | Pathogenic, Likely pathogenic | `filter_clinvar` |
| germline | InterVar | Pathogenic, Likely pathogenic | `filter_intervar` |
| germline | ReNOVo | LP Pathogenic, IP Pathogenic, HP Pathogenic | `filter_renovo` |
| germline | ClinVar significance | Pathogenic, Likely pathogenic | `filter_clinvar` |

A gene list may be added as a further **AND** — `filter_genes_somatic`, `filter_genes_germline`
— and no list is applied by default. CIViC drops out of the somatic list under `--skip_civic`,
and also if the MAF has no CIViC column, with a warning in either case.

**2. Or it is retained as pathogenic**, which bypasses all of the above — depth and VAF
included. The retention set is *narrower* than the filter lists: somatic keeps CancerVar
Tier_I_strong / Tier_II_potential, ClinVar pathogenic, and CIViC **A or B only**; germline keeps
InterVar Pathogenic / Likely pathogenic and ClinVar pathogenic — **not** ReNOVo.
`--skip_pathogenic_retention` turns this second route off entirely, at which point the filters
above are the only way through.

> **Note:** MAFigate applies this same filtering code — `streamlit_app/vendor/` is a
> byte-for-byte copy of it — so the thresholds you see in the app are these thresholds.

## Modules

### Database Update (`update_db/`)
Tools for updating Funcotator and Annovar databases. **Recommended before first use.**

See [update_db/README.md](update_db/README.md) for instructions, or
[update_db/README_DOCKER.md](update_db/README_DOCKER.md) to run them in a container.

### Samplesheet helpers (`utils/`)
Three shell scripts that build the annotation and biomarker samplesheets from a directory of
results. See [utils/README.md](utils/README.md).

### Data Visualization (`streamlit_app/`)
**MAFigate** — filter and browse the pipeline's annotated MAF interactively: clinical,
frequency, gene and quality filters you can change as you go, with the variants that passed
and the variants that did not side by side.

See [streamlit_app/README.md](streamlit_app/README.md) to install and run it.

## Citations

The tools and databases this pipeline builds on are listed, with their papers, in
[CITATIONS.md](CITATIONS.md).

## Contributing

Issues and pull requests are welcome. This repository is a one-way export of a private
development tree, which changes what happens to yours after you send it, and how you are
credited for it — please read [CONTRIBUTING.md](CONTRIBUTING.md) first.

## License

This software is provided for research use only. See [LICENSE](LICENSE) for details.

**Disclaimer:** VariantAlker assumes no responsibility for any injury, damage, errors, or omissions arising from its use. Users assume all risks associated with using this software.
