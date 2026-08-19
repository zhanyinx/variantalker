# Database Update Docker/Singularity Environment

This Docker image contains all necessary tools for updating the VariantAlker databases.

**Docker Hub:** `yinxiu/variantalker_db:latest`

## Included Software

- **tabix** (v1.19) - VCF indexing
- **bgzip** (v1.19) - Block compression
- **vt** (latest) - Variant normalization and decomposition
- **GATK** (v4.5.0.0) - Genome Analysis Toolkit
- **Python 3** with psycopg, requests, beautifulsoup4, numpy, pandas, lxml
- **Perl** - For Annovar indexing scripts
- **Java 17** - Required for GATK

## Quick Start

Go within the update_db folder in variantalker

```bash
cd path2/variantalker/update_db
```

### Option 1: Use Pre-built Docker Image

Pull the image from Docker Hub:
```
docker pull yinxiu/variantalker_db:latest
```

### Option 2: Build Locally

```bash
cd update_db
docker build -t yinxiu/variantalker_db:latest .
```

## Running with Docker

### Interactive mode
```bash
docker run -it --rm \
  -v $(pwd):/workspace \
  yinxiu/variantalker_db:latest
```

### Run Funcotator update
```bash
docker run --rm \
  -v $(pwd):/workspace \
  -v /path/to/funcotator_somatic:/somatic_db \
  -v /path/to/funcotator_germline:/germline_db \
  -v /path/to/backup:/backup \
  yinxiu/variantalker_db:latest \
  bash -c "cd /workspace && python3 scripts/update_funcotator.py \
    -sd /somatic_db \
    -gd /germline_db \
    -b /backup"
```

### Run Annovar update
```bash
docker run --rm \
  -v $(pwd):/workspace \
  -v /path/to/annovar/humandb:/humandb \
  -v /path/to/annovar:/annovar \
  yinxiu/variantalker_db:latest \
  bash -c "cd /workspace && python3 scripts/update_annovar.py \
    --annovar_db_path /humandb \
    --annovar_download_script /annovar/annotate_variation.pl \
    -v /usr/local/bin"
```

## Running with Singularity

For HPC environments that use Singularity instead of Docker:

### Build Singularity image from Docker Hub
```bash
singularity pull variantalker_db.sif docker://yinxiu/variantalker_db:latest
```

### Interactive mode
```bash
singularity shell \
  --bind $(pwd):/workspace \
  variantalker_db.sif
```

### Run Funcotator update
```bash
singularity exec \
  --bind $(pwd):/workspace \
  --bind /path/to/funcotator_somatic:/somatic_db \
  --bind /path/to/funcotator_germline:/germline_db \
  --bind /path/to/backup:/backup \
  variantalker_db.sif \
  bash -c "cd /workspace && python3 scripts/update_funcotator.py \
    -sd /somatic_db \
    -gd /germline_db \
    -b /backup"
```

### Run Annovar update
```bash
singularity exec \
  --bind $(pwd):/workspace \
  --bind /path/to/annovar/humandb:/humandb \
  --bind /path/to/annovar:/annovar \
  variantalker_db.sif \
  bash -c "cd /workspace && python3 scripts/update_annovar.py \
    --annovar_db_path /humandb \
    --annovar_download_script /annovar/annotate_variation.pl \
    -v /usr/local/bin"
```

## Environment Variables

Not required for the update scripts. The scripts work with file-based databases and do not need PostgreSQL access during updates.

## Verifying Installation

### Docker
```bash
docker run --rm yinxiu/variantalker_db:latest bash -c "
  echo 'bgzip:' && bgzip --version && \
  echo 'tabix:' && tabix --version && \
  echo 'vt:' && vt --version && \
  echo 'gatk:' && gatk --version
"
```

### Singularity
```bash
singularity exec variantalker_db.sif bash -c "
  echo 'bgzip:' && bgzip --version && \
  echo 'tabix:' && tabix --version && \
  echo 'vt:' && vt --version && \
  echo 'gatk:' && gatk --version
"
```

## Example: Complete Update Workflow

### Using Docker
```bash
# Navigate to update_db directory
cd /path/to/variantalker/update_db

# Pull latest image
docker pull yinxiu/variantalker_db:latest

# Run Funcotator update
docker run --rm \
  -v $(pwd):/workspace \
  -v /data/funcotator_somatic:/somatic_db \
  -v /data/funcotator_germline:/germline_db \
  -v /data/backup:/backup \
  yinxiu/variantalker_db:latest \
  bash -c "cd /workspace && python3 scripts/update_funcotator.py \
    -sd /somatic_db \
    -gd /germline_db \
    -b /backup"

# Run Annovar update
docker run --rm \
  -v $(pwd):/workspace \
  -v /data/annovar/humandb:/humandb \
  -v /path/to/annovar:/annovar \
  yinxiu/variantalker_db:latest \
  bash -c "cd /workspace && python3 scripts/update_annovar.py \
    --annovar_db_path /humandb \
    --annovar_download_script /annovar/annotate_variation.pl \
    -v /usr/local/bin"
```

### Using Singularity (HPC)
```bash
# Navigate to update_db directory
cd /path/to/variantalker/update_db

# Pull and convert image
singularity pull variantalker_db.sif docker://yinxiu/variantalker_db:latest

# Run Funcotator update
singularity exec \
  --bind $(pwd):/workspace \
  --bind /data/funcotator_somatic:/somatic_db \
  --bind /data/funcotator_germline:/germline_db \
  --bind /data/backup:/backup \
  variantalker_db.sif \
  bash -c "cd /workspace && python3 scripts/update_funcotator.py \
    -sd /somatic_db \
    -gd /germline_db \
    -b /backup"

# Run Annovar update
singularity exec \
  --bind $(pwd):/workspace \
  --bind /data/annovar/humandb:/humandb \
  --bind /path/to/annovar:/annovar \
  variantalker_db.sif \
  bash -c "cd /workspace && python3 scripts/update_annovar.py \
    --annovar_db_path /humandb \
    --annovar_download_script /annovar/annotate_variation.pl \
    -v /usr/local/bin"
```

## Notes

- The container uses Ubuntu 22.04 as base image
- Java 17 is installed for GATK compatibility
- htslib 1.19 provides both tabix and bgzip
- vt is compiled from source from the official repository
- All tools are available in the PATH
