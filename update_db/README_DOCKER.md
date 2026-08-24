# Database Update Docker/Singularity Environment

This Docker image contains all necessary tools for updating the VariantAlker databases.

**Docker Hub:** `yinxiu/variantalker_db:latest` (public; no login required to pull)

> **Read the warning below before you run anything.** The `/workspace` mount must be the
> **repository root**, not `update_db`, or the Annovar update cannot start (see
> [The `/workspace` mount](#the-workspace-mount)). Neither update can any longer report work it did
> not do (see [Known limitations](#known-limitations)).

The container removes the *software* prerequisites, not the *hardware* ones. These updates are
meant for an HPC system or a server: the databases need **at least 500 GB** of free disk, plus
working room for a run in progress, so the host has to have the space whatever the image contains.
See [Where this runs](README.md#where-this-runs).

## Included Software

Versions below were read out of `yinxiu/variantalker_db:latest`
(digest `sha256:1f3d0f02cc5630a573aee42beafcd437855fbf64d922b37f98b55020eade0fa9`) on 2026-08-20
with the [verification block](#verifying-installation):

- **tabix** (htslib 1.19) - VCF indexing
- **bgzip** (htslib 1.19) - Block compression
- **vt** (v0.57721) - Variant normalization and decomposition, installed at `/usr/local/bin/vt`
- **GATK** (v4.5.0.0, HTSJDK 4.1.0, Picard 3.1.1) - Genome Analysis Toolkit
- **Python 3** (3.10.12) with psycopg, requests, beautifulsoup4, numpy, pandas, lxml
- **Perl** (5.34.0) - For Annovar indexing scripts
- **Java 17** (OpenJDK 17.0.17) - Required for GATK

### Platform

The published image is **`linux/amd64` only**. It runs on Apple Silicon under Docker Desktop's
emulation, but every command prints:

```
WARNING: The requested image's platform (linux/amd64) does not match the detected host platform
(linux/arm64/v8) and no specific platform was requested
```

That warning is expected and is not an error. Pass `--platform linux/amd64` explicitly to silence it.

## Quick Start

Go to the **root** of the variantalker repository — not into `update_db`. The update scripts reach
back into `resources/` (see [The `/workspace` mount](#the-workspace-mount)), so the whole repository
has to be visible inside the container.

```bash
cd /path/to/variantalker
```

### Option 1: Use Pre-built Docker Image

Pull the image from Docker Hub:
```bash
docker pull yinxiu/variantalker_db:latest
```

### Option 2: Build Locally

> **Building under the published tag overwrites the image you just pulled.** The obvious command,
> `docker build -t yinxiu/variantalker_db:latest .`, reuses the published tag: the build succeeds,
> exits 0, says nothing about it, and every later `docker run` in this document then uses your build
> instead of the published image. On Apple Silicon it also swaps the architecture. Measured here:
>
> ```
> before:  yinxiu/variantalker_db:latest  sha256:9f9c5d650e94  amd64   (pulled)
> after:   yinxiu/variantalker_db:latest  sha256:55f33bdcf153  arm64   (built)
>          yinxiu/variantalker_db:<none>  sha256:9f9c5d650e94          (pulled image, now dangling)
> ```
>
> Build under a distinct tag unless you mean to replace it:

```bash
docker build -t variantalker_db:local update_db
```

Verified from the repository root on 2026-08-20: exit 0, and the resulting image reports the same
four tool versions as the published one. A cold build takes several minutes — it compiles htslib and
vt from source and downloads GATK.

If you do rebuild the published tag deliberately, `docker pull yinxiu/variantalker_db:latest`
afterwards to get the published image back.

## The `/workspace` mount

Mount the **repository root** at `/workspace`, and invoke the scripts as `update_db/scripts/...`.

`update_db/scripts/update_annovar.py` resolves its CancerVar and InterVar defaults relative to its
own location, two directories up:

```
scripts/../../resources/configs/config.init.CancerVar
scripts/../../resources/CancerVar/CancerVar.py
scripts/../../resources/InterVar/Intervar.py
```

`update_db/README.md` explains why those must be reachable and writable:

> since these databases are used by CancerVar and InterVar, we need to update also the corresponding
> CancerVar and InterVar. If CancerVar.py and InterVar.py and the corresponding configs are not
> provided (RECOMMENDED CASE), the tool will automatically update the CancerVar.py and InterVar.py
> and the relative configs in the resources

If you mount only `update_db` at `/workspace`, that path arithmetic escapes the mount and lands on
`/resources/...`, which does not exist in the image. The run then dies immediately, before touching
any database:

```
FileNotFoundError: [Errno 2] No such file or directory:
  '/workspace/scripts/../../resources/configs/config.init.CancerVar'
```

Mounting the repository root makes it resolve to `/workspace/resources/...`. The Funcotator script
works either way — its only reach is `scripts/../funcotator` — but use the same layout for both so
there is one thing to remember.

> **Why this differs from `update_db/README.md`, deliberately.** That document has you `cd` into
> `update_db` and run `python3 scripts/update_annovar.py`, which is correct *there*: on the host,
> `scripts/../../resources` still lands inside your checkout, because `update_db`'s parent **is** the
> repository root. In a container it is not — `/workspace`'s parent is `/`. That is the whole of the
> difference. Do not "harmonise" this document onto `scripts/...`; it would restore the
> `FileNotFoundError` above.

## What the container writes to your host

Every `-v` bind is read-write, so the container writes through to your host filesystem:

- **The database directories you mount** (`/somatic_db`, `/germline_db`, `/humandb`) are modified in
  place. The previous version is copied into whatever you mounted at `/backup`, under a
  `YYYYMMDD/` subdirectory.
- **`resources/CancerVar/CancerVar.py`, `resources/InterVar/Intervar.py` and
  `resources/configs/*`** in your checkout are rewritten in place by the Annovar update, as quoted
  above — `utils_update.py` rewrites the database names in both scripts and in both config files.
  Commit or copy your checkout first if you care about those files. (Stated from the code and from
  `update_db/README.md`; not observed directly during the audit, because a run without real Annovar
  databases stops before reaching that step.) Only databases whose new table was validated are
  rewritten, so a failed update leaves both files naming the version still on disk.
- **A log file is written into the current working directory**, i.e. into the root of your checkout:
  `YYYYMMDD_funcotator_update.log` or `YYYYMMDD_cancervar_intervar_update.log`. Read it — it is the
  only report of what happened, subject to [Known limitations](#known-limitations).
- **`update_db/scripts/__pycache__/`** appears in your checkout. Both it and `*.log` are already in
  `.gitignore`.

Docker runs as root inside the container by default. On Docker Desktop for macOS ownership is
remapped and the files come out owned by you — that is what the audit observed. Expect root-owned
files on a Linux host instead; that case was not tested here.

## Running with Docker

### Interactive mode
```bash
docker run -it --rm \
  -v $(pwd):/workspace \
  yinxiu/variantalker_db:latest
```

The image's `WORKDIR` is `/workspace`, so this drops you straight into the mounted repository.

### Run Funcotator update
```bash
docker run --rm \
  -v $(pwd):/workspace \
  -v /path/to/funcotator_somatic:/somatic_db \
  -v /path/to/funcotator_germline:/germline_db \
  -v /path/to/backup:/backup \
  yinxiu/variantalker_db:latest \
  bash -c "cd /workspace && python3 update_db/scripts/update_funcotator.py \
    -sd /somatic_db \
    -gd /germline_db \
    -b /backup"
```

`update_funcotator.py` also accepts `--build hg19`, `--build hg38` or `--build both` (default
`both`); append it to the `python3` line to restrict the update to one genome build.

### Run Annovar update
```bash
docker run --rm \
  -v $(pwd):/workspace \
  -v /path/to/annovar/humandb:/humandb \
  -v /path/to/annovar:/annovar \
  yinxiu/variantalker_db:latest \
  bash -c "cd /workspace && python3 update_db/scripts/update_annovar.py \
    --annovar_db_path /humandb \
    --annovar_download_script /annovar/annotate_variation.pl \
    -v /usr/local/bin"
```

`-v /usr/local/bin` is the directory holding `vt` **inside the container**, which is where the image
puts it — confirmed with `which vt`. Do not substitute a host path here. (The native instructions in
`update_db/README.md` say `-v /path/to/vt` because there you supply your own build of vt.)

## Known limitations

**Both updates used to be able to report work they had not done. Neither can now** — the Funcotator
route by issue 353, the CancerVar/InterVar (Annovar) route by issue 354. What follows is what they do
instead, because it changes how you read the log.

The defect was one thing with two halves. The image is Ubuntu 22.04, where `/bin/sh` is **dash**,
and every shell script under `update_db/` is `#!/bin/bash` and uses bash-only syntax — so handing
one to `sh` killed it, either on a parse error or at the first `[[ ... ]]` with every conditional
silently false. `utils_update.py` did that at seven call sites, and inspected the outcome at none of
them. All seven now run under `bash`, and both routes validate what the script produced before
anything starts using it.

Both routes run their scripts under `bash`, inspect what was produced before replacing or adopting
anything, and report one of four verdicts derived from that inspection:

```
UPDATED: clinvar hg38 20250114 -> 20260820
CURRENT: dbsnp hg38 (b156) - the published md5sum matches the installed database
SKIPPED: cosmic - update disabled, COSMIC changed its download format
FAILED:  clinvar hg19 - staged clinvar_.vcf is 0 bytes; live database untouched
RESULT:  1 updated, 1 current, 1 skipped, 1 failed
```

`SUCCESS` is gone from both routes: it had been doing four different jobs, one of which was "did not
happen at all". Two consequences worth knowing when you run either:

- **A failed step leaves the live database exactly as it was**, in both the somatic and the germline
  tree, and the run carries on to the next database instead of aborting. What used to happen is that
  a failed ClinVar update replaced *both* databases with an empty directory and logged
  `SUCCESS: clinvar updated successfully.` at exit **0**.
- **The process exits 1 if any step failed**, so a scheduled run is now checkable without parsing
  the log. Per-build skips are benign and do not affect the exit code; a step that skipped every
  build it was asked for is a failure.

`--build hg38` now updates hg38 and leaves hg19 alone, in both trees. Previously a single-build run
replaced the whole database directory, so the build you did not ask for was deleted — silently, and
no later run noticed, because the "no existing database" checks only ever looked at the build that
*was* requested. (`--build` is a Funcotator option; `update_annovar.py` has none, and does not need
one — see below.)

**The Annovar route's risk was never the same one, so read its verdicts differently.** It keeps
everything in one flat `humandb` as `<build>_<database>.txt`, with the version *in the filename*, so
every download lands on a name nothing else uses and no update can overwrite or delete the database
it replaces. What puts a database into use there is instead the two places that *name* it: the
`-protocol` list inside `CancerVar.py` / `Intervar.py`, and the `database_names` line in
`resources/configs/*`. So on that route:

- a verdict is **per database, not per build** — a name advances for both builds or for neither,
  because `-protocol` carries no build;
- a `FAILED` database means both names were left alone, so annotation carries on with the previous
  version, which is still on disk;
- a partially-downloaded new version can be left behind on disk without being adopted. That is
  harmless and the next run picks it up.

Previously this route rewrote the two annotation scripts with a GNU-only `sed -i` that nothing
checked — a silent no-op on macOS — and wrote `database_names` unconditionally at the end, so a
failure in between could leave the script naming one version and the config another. Both names now
move together, and only after the table they name has been read and found to contain data.

Separately, and by design rather than by defect, the **COSMIC** and **GENCODE** updates are disabled
in the scripts themselves and log `SKIPPED` — see the note in `update_db/README.md`.

## Running with Singularity

> **Not verified.** Neither `singularity` nor `apptainer` is installed on the machine where this
> document was audited (2026-08-20), so none of the commands in this section have been executed.
> They are given as the direct translation of the Docker commands above, which *were* executed, and
> the `/workspace` mount rule and the issue 354 limitation apply identically because they are
> properties of the image and the scripts, not of the container runtime. Verifying this section
> requires a host with Singularity or Apptainer installed.

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
  bash -c "cd /workspace && python3 update_db/scripts/update_funcotator.py \
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
  bash -c "cd /workspace && python3 update_db/scripts/update_annovar.py \
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

Run on 2026-08-20, this prints (abridged):

```
bgzip:
bgzip (htslib) 1.19
tabix:
tabix (htslib) 1.19
vt:
vt v0.57721
gatk:
The Genome Analysis Toolkit (GATK) v4.5.0.0
HTSJDK Version: 4.1.0
Picard Version: 3.1.1
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
# Navigate to the repository root
cd /path/to/variantalker

# Pull latest image
docker pull yinxiu/variantalker_db:latest

# Run Funcotator update
docker run --rm \
  -v $(pwd):/workspace \
  -v /path/to/funcotator_somatic:/somatic_db \
  -v /path/to/funcotator_germline:/germline_db \
  -v /path/to/backup:/backup \
  yinxiu/variantalker_db:latest \
  bash -c "cd /workspace && python3 update_db/scripts/update_funcotator.py \
    -sd /somatic_db \
    -gd /germline_db \
    -b /backup"

# Run Annovar update
docker run --rm \
  -v $(pwd):/workspace \
  -v /path/to/annovar/humandb:/humandb \
  -v /path/to/annovar:/annovar \
  yinxiu/variantalker_db:latest \
  bash -c "cd /workspace && python3 update_db/scripts/update_annovar.py \
    --annovar_db_path /humandb \
    --annovar_download_script /annovar/annotate_variation.pl \
    -v /usr/local/bin"
```

### Using Singularity (HPC)

Not verified — see the note under [Running with Singularity](#running-with-singularity).

```bash
# Navigate to the repository root
cd /path/to/variantalker

# Pull and convert image
singularity pull variantalker_db.sif docker://yinxiu/variantalker_db:latest

# Run Funcotator update
singularity exec \
  --bind $(pwd):/workspace \
  --bind /path/to/funcotator_somatic:/somatic_db \
  --bind /path/to/funcotator_germline:/germline_db \
  --bind /path/to/backup:/backup \
  variantalker_db.sif \
  bash -c "cd /workspace && python3 update_db/scripts/update_funcotator.py \
    -sd /somatic_db \
    -gd /germline_db \
    -b /backup"

# Run Annovar update
singularity exec \
  --bind $(pwd):/workspace \
  --bind /path/to/annovar/humandb:/humandb \
  --bind /path/to/annovar:/annovar \
  variantalker_db.sif \
  bash -c "cd /workspace && python3 update_db/scripts/update_annovar.py \
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
- `/bin/sh` in the image is **dash**, not bash — the cause of issue 339, and deliberately left that
  way: macOS `/bin/sh` is bash, so the container is the only place this class of bug is observable.
  Note that this is the *only* sense in which running natively on macOS is better: the native route
  is not supported there at all, for reasons of both disk and GNU userland — see
  [Platform support](README.md#platform-support)
