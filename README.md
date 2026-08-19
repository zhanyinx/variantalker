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
- [Contributing](#contributing)
- [License](#license)

## Overview

VariantAlker is a comprehensive pipeline for annotating and prioritizing genetic variants in cancer genomics. It integrates multiple annotation databases and tools to identify clinically relevant mutations and extract actionable biomarkers.

**Supported platforms:** Dragen, nf-sarek, ION-Torrent  
**Supported genomes:** hg19, hg38

## Features

### Core Annotation
- **Somatic variants** - CancerVar prioritization and CIViC evidence levels
- **Germline variants** - InterVar classification 
- **CNV annotation** - CNVKit (nf-sarek) and Dragen outputs
- **Multi-database integration** - Funcotator, Annovar, CIViC, AlphaMissense

### Biomarker Extraction (Beta)
- Tumor Mutational Burden (TMB)
- Mutational signatures (APOBEC, UV, tobacco)
- Clonal TMB (requires BAM/CRAM files)
- Gene expression (requires RNA-seq data)
- Gene CNV

See [biomarker documentation](docs/biomarkers/README.md) for details.

### Data Visualization
Filter and browse the annotated MAF interactively with the **MAFigate** app (see [`streamlit_app/`](streamlit_app/README.md)).

## Installation

### 1. Clone the repository
```bash
git clone https://github.com/zhanyinx/variantalker.git
cd variantalker
```

### 2. Download annotation databases
Separate databases are available for hg19 and hg38:

```bash
# For hg38
wget -r -N --no-parent -nH --cut-dirs=3 -P public_databases/hg38 https://bioserver.ieo.it/repo/dima/2026/hg38 

# For hg19
wget -r -N --no-parent -nH --cut-dirs=3 -P public_databases/hg19 https://bioserver.ieo.it/repo/dima/2026/hg19

# Annovar databases
wget -r -N --no-parent -nH --cut-dirs=3 -P public_databases/humandb https://bioserver.ieo.it/repo/dima/2026/humandb
```

### 3. Configure the pipeline
Edit [nextflow.config](nextflow.config) with your database paths:

```groovy
params {
  funcotator_germline_db = "path/to/public_databases/funcotator_dataSources.v1.7.20200521g"
  funcotator_somatic_db  = "path/to/public_databases/funcotator_dataSources.v1.7.20200521s"
  annovar_db             = "path/to/public_databases/humandb"
  annovar_software_folder = "path/to/annovar"
  alpha_mis_genome_basedir = "path/to/public_databases"
  fasta                  = "path/to/reference.fasta"
  target                 = "path/to/target.bed"
}
```

## Quick Start

> **Nextflow version / language parser**
> Recent Nextflow releases ship a new, stricter language parser as the default,
> which is **not** compatible with this pipeline's current syntax. Run with the
> legacy parser by setting `NXF_SYNTAX_PARSER=v1` (shown below), or use an older
> Nextflow release (≥ 22.10.1). All commands in this README assume the legacy
> parser.

### Basic annotation
```bash
NXF_SYNTAX_PARSER=v1 nextflow run main.nf -profile singularity --input samplesheet.csv --outdir results
```

### View all options
```bash
NXF_SYNTAX_PARSER=v1 nextflow run main.nf --help --show_hidden_params
```

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

### Sample Types

- **somatic** - Single-sample (tumor-only) or multi-sample (tumor-normal) VCF. Requires `tumor_tissue`.
- **germline** - Single-sample VCF. `tumor_tissue` optional.
- **cnv** - CNVKit `.cnr` file (nf-sarek) or VCF (Dragen). `tumor_tissue` optional.

### Tumor Tissue Options
Adrenal_Gland, Bile_Duct, Bladder, Blood, Bone, Bone_Marrow, Brain, Breast, Cancer_all, Cervix, Colorectal, Esophagus, Eye, Head_and_Neck, Inflammatory, Intrahepatic, Kidney, Liver, Lung, Lymph_Nodes, Nervous_System, Other, Ovary, Pancreas, Pleura, Prostate, Skin, Soft_Tissue, Stomach, Testis, Thymus, Thyroid, Uterus

## Output

The pipeline generates the following directory structure:

```
results/
└── YYYY-MM-DD/
    └── annotation/
        ├── somatic/
        │   └── patient_name/
        │       ├── patient.maf                    # Full MAF with all annotations
        │       ├── patient.vcf                    # PASS variants only
        │       ├── filtered.patient.maf.pass.tsv  # Filtered variants (passing filters)
        │       └── filtered.patient.maf.nopass.tsv # Filtered variants (not passing)
        ├── germline/
        │   └── patient_name/
        │       ├── patient.maf
        │       ├── patient.vcf
        │       ├── filtered.patient.maf.pass.tsv
        │       └── filtered.patient.maf.nopass.tsv
        └── cnv/
            └── patient_name/
                └── patient.cnv.annotated.tsv
```

### Default Filters

Variants are filtered based on the following criteria:

**Excluded variant types:**
- Silent, IGR, RNA (unless pathogenic/likely pathogenic)

**Coverage thresholds:**
- Minimum depth: 50× (unless pathogenic/likely pathogenic)

**Variant allele frequency:**
- Somatic: VAF ≥ 0.01
- Germline: VAF ≥ 0.2

**Classification filters (OR logic):**
- **InterVar** (germline): Pathogenic, Likely pathogenic
- **CancerVar** (somatic): Tier_I_strong, Tier_II_potential
- **ReNOVo**: LP Pathogenic, IP Pathogenic, HP Pathogenic
- **CIViC**: Evidence levels A, B, C

> **Note:** A variant passes if it meets at least one of the classification criteria (OR logic).

## Modules

### Database Update (`update_db/`)
Tools for updating Funcotator and Annovar databases. **Recommended before first use.**

See [update_db/README.md](update_db/README.md) for instructions.

### Data Visualization (`streamlit_app/`)
**MAFigate** — filter and browse the pipeline's annotated MAF interactively: clinical,
frequency, gene and quality filters you can change as you go, with the variants that passed
and the variants that did not side by side.

See [streamlit_app/README.md](streamlit_app/README.md) to install and run it.

## Contributing

Issues and pull requests are welcome. This repository is a one-way export of a private
development tree, which changes what happens to yours after you send it, and how you are
credited for it — please read [CONTRIBUTING.md](CONTRIBUTING.md) first.

## License

This software is provided for research use only. See [LICENSE](LICENSE) for details.

**Disclaimer:** VariantAlker assumes no responsibility for any injury, damage, errors, or omissions arising from its use. Users assume all risks associated with using this software.
