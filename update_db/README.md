[![GitHub code licence is MIT](https://img.shields.io/badge/license-MIT-brightgreen.svg)]

# Update scripts to get the most recent version of annovar and funcotator databases

## Contents
- [Contents](#contents)
- [Overview](#overview)
- [Where this runs](#where-this-runs)
- [Installation](#installation)
- [Documentation](#documentation)
- [Usage](#usage)
- [Tests](#tests)

## Overview
Regular updates are made to the resources used for annotating variants, and as a result, some variants previously classified as Variants of Uncertain Significance (VUS) may be re-annotated as Pathogenic. Therefore, it is crucial to use up-to-date databases for functional annotation. These utilities offer tools that enable automatic updates of funcotator and annovar databases.

## Where this runs

**These utilities target HPC systems and servers. A laptop or workstation is not a supported
environment, and the reason is disk rather than CPU.**

The annotation databases these scripts maintain run to several hundred gigabytes — budget **at
least 500 GB** of free space for the databases themselves, and working room on top of that for a
run in progress. The dbSNP step shows why the working room matters: it downloads the compressed
VCF, writes a **fully decompressed** rewrite of it beside the download, and only then bgzips that
back and deletes the uncompressed copy — so all three exist at once
(`funcotator/update_dbsnp/update_dbsnp.sh`, the `zcat … > $name` / `bgzip -c $name` / `rm $name`
sequence).

Both installation routes below target **Linux**. The container route is the recommended one on any
host; the manual route is **Linux-only** — see [Platform support](#platform-support) for what
breaks on macOS and why it is not being fixed.

None of this applies if you only want to read the code or run the test suite: the tests download
nothing and touch no database, so they run anywhere in about a second. See [Tests](#tests).

## Installation

### Option 1: Using Docker/Singularity (Recommended)

All required software is pre-installed in the Docker image. See [README_DOCKER.md](README_DOCKER.md) for detailed instructions.

```bash
# Using Docker
docker pull yinxiu/variantalker_db:latest

# Using Singularity (HPC)
singularity pull variantalker_db.sif docker://yinxiu/variantalker_db:latest
```

The image is built for `linux/amd64` only. On an Apple Silicon Mac it still runs, under
emulation, and every `docker run` prints a platform-mismatch warning that can be ignored.

### Option 2: Manual Installation

#### Required Software

1. **Python 3.9+** with the following libraries:
```bash
conda create -n update python=3.9 -y
conda activate update
pip install beautifulsoup4 numpy pandas requests lxml
```

`lxml` is not imported directly: it is the parser the scripts ask BeautifulSoup and
`pandas.read_html` for, so leaving it out fails at run time rather than at import.

2. **vt** - Variant normalization and decomposition tool
```bash
git clone --recursive https://github.com/atks/vt.git
cd vt
make
```

The scripts locate `vt` through the `-v/--vt` argument, not through `PATH`, and that
argument is the **directory** holding the binary — after `make` that is the clone root,
so `-v /path/to/vt` is what you pass.

*The following applies to macOS, which is **not** a supported host for this route (see
[Platform support](#platform-support)). It is kept because it is the one part of a Mac build that
is not obvious, and it is useful if you are building vt on a Mac for some other purpose.*

On macOS the link step fails with `ld: library 'crypto' not found`: vt's Makefile appends
`-lcrypto` unconditionally and macOS ships no OpenSSL for the linker to find. The
Makefile has no `LDFLAGS` hook on that rule — the library list is written into the recipe
— so the search path has to arrive through `CXXFLAGS`, which the recipe does expand:

```bash
make CXXFLAGS="-pipe -std=c++0x -O3 -I./lib -I. -I./lib/htslib -I./lib/Rmath -I./lib/pcre2 -D__STDC_LIMIT_MACROS -L$(brew --prefix openssl@3)/lib"
```

3. **GATK 4.5.0.0** - Genome Analysis Toolkit
```bash
wget https://github.com/broadinstitute/gatk/releases/download/4.5.0.0/gatk-4.5.0.0.zip
unzip gatk-4.5.0.0.zip
# Add gatk-4.5.0.0/gatk to your PATH
# Requires Java 17
```

The download is 741 MB. 4.5.0.0 is the version these scripts were tested against and the
one the Docker image ships; the only GATK tools they call are `IndexFeatureFile`,
`CreateSequenceDictionary` and `VariantsToTable`, all of which have been present
throughout GATK 4, so a newer 4.x release is very likely to work.

4. **htslib** - Provides tabix and bgzip
```bash
wget https://github.com/samtools/htslib/releases/download/1.19/htslib-1.19.tar.bz2
tar -xjf htslib-1.19.tar.bz2
cd htslib-1.19
./configure --prefix=/usr/local
make
sudo make install
```

Installing into `/usr/local` needs root. Without it, either keep the build tree and add
it to `PATH`, or configure a prefix you own — `./configure --prefix=$HOME/.local` — and
run `make install` unprivileged.

5. **Annovar** - provides `annotate_variation.pl`, which `update_annovar.py` requires

Annovar is not redistributable and is not vendored here. Register at
[ANNOVAR](https://annovar.openbioinformatics.org/en/latest/user-guide/download/); the
download link arrives by email. Unpack it anywhere and pass
`--annovar_download_script <annovar>/annotate_variation.pl`.

`--annovar_db_path` is the directory the Annovar databases are kept in — conventionally
`humandb` beside that script. A fresh unpack does not contain one, so create it before
the first run: every directory argument is validated up front and a missing one is
rejected at argument parsing (`Input value must be folder and must exist`). The same
applies to `-sd`, `-gd` and `-v`.

6. **Perl** - Required for Annovar indexing scripts (usually pre-installed)

7. **Java 17** - Required for GATK
```bash
# Ubuntu/Debian
apt-get install openjdk-17-jdk

# Conda
conda install -c conda-forge openjdk=17
```

The `apt-get` line is Debian/Ubuntu only. On any other distribution — or on an HPC system where
you cannot install packages — use the conda line, or load whatever JDK 17 module the site provides.

8. **System utilities** - `wget`, `curl` and `git`, plus **GNU** `gzip` and `sed`

The scripts call all of these by bare name and assume the GNU versions; that assumption is what
makes this route Linux-only ([Platform support](#platform-support)). A desktop Linux distribution
normally has the lot already, but a minimal container or a stripped-down image may well not — in
particular `wget` and `curl` are both used, on different scripts, so having only one of them is not
enough.

#### What must be on PATH

The Python scripts shell out to the tools below by bare name, so they have to be
findable on `PATH`:

- `gatk` (both update paths)
- `tabix` and `bgzip` (the Funcotator path)
- `perl` (the Annovar path, for Annovar's own indexing scripts)

`vt` and Annovar's `annotate_variation.pl` are the exception: they are passed in by
path, through `-v` and `--annovar_download_script`.

#### Clone Repository
```bash
git clone https://github.com/zhanyinx/variantalker.git
cd variantalker/update_db
```

### Platform support

| Host | Container route (Option 1) | Manual route (Option 2) |
|---|---|---|
| Linux | Supported | Supported |
| macOS | Runs (`linux/amd64` under emulation), subject to the disk requirement | **Not supported** |
| Windows | Not tested | Not supported |

The manual route is Linux-only because the update scripts are written against **GNU userland**, and
macOS ships BSD versions of the same commands. Three spellings break there:

- **`zcat`** — on macOS this is the `compress`-era tool. It appends `.Z` to the name it is given and
  cannot open a `.gz` at all, so a download that arrived perfectly intact still fails:
  `zcat: can't stat: clinvar.vcf.gz (clinvar.vcf.gz.Z)`. Eight call sites, across the two ClinVar
  scripts and `update_dbsnp.sh`.
- **`sed -i EXPR`** — GNU `sed` treats the expression as the script to run; BSD `sed` reads it as
  `-i`'s backup-suffix argument, exits 1 and leaves the file unchanged. Four call sites in
  `funcotator/update_clinvar/update_clinvar_funcotator.sh`.
- **`wget`** — not present on a stock macOS install at all.

This is tracked as issue 370 and is deliberately **not** being fixed. The disk requirement in
[Where this runs](#where-this-runs) rules a Mac out whatever its userland does, so making these
spellings portable would buy a route nobody could use anyway. Run the container instead.

One consequence worth knowing rather than rediscovering: **the failure is loud on both routes.**
Every live update script now runs under `set -eo pipefail` (issues 353 and 354), so the first
`zcat` aborts the run and the update reports `FAILED`. Nothing is destroyed, and nothing stale is
adopted.

That is a recent property and it is the half that matters. `date=$(zcat … | grep … | tr …)` takes
its exit status from the **last** command in the pipeline, so without `pipefail` a failed `zcat`
was invisible: `$date` came out empty and the script went on to build a database whose name carried
no version. Attempting the manual route on a Mac is a wasted run, not a corrupted database.

## Documentation

These utilities determine if updates are necessary for specific databases. For some databases, a warning will be generated if updates are required without actually updating them, while others will be automatically updated, replacing the current database version.

List of Funcotator databases:

- Achilles, cancer-gene-census, familial & simple-uniprot (from Oncotator; **check DISABLED** - the Oncotator check is commented out in `update_funcotator.py`, so these four are neither checked nor updated)
- dna_repair_genes (check only)
- oreganno (**check DISABLED** - oreganno.org has been down since November 2023, so the check is commented out)
- **cosmic (DISABLED - update script not compatible with new COSMIC download format)**
- **gencode (DISABLED - the updated database is built correctly but Funcotator does not use it properly)**
- clinvar (update if a newer version is available)
- hgnc (update if a newer version is available)
- dbsnp (update if a newer version is available)
- acmg_rec (update if a newer version is available)

List of Annovar databases:

- refGene (ignored)
- ensGene (ignored)
- knownGene (ignored)
- esp6500siv2_all (ignored)
- 1000g2015aug (ignored)
- exac03 (ignored)
- dbscsnv11 (ignored)
- dbnsfp\*_interpro (update if a newer version is available from the Annovar website)
- gnomad_genome (ignored)
- rmsk (ignored)
- avsnp (update if a newer version is available from the Annovar website)
- dbnsfp (update if a newer version is available from the Annovar website)
- clinvar (update if a newer version is available)
- **cosmic (DISABLED - update script not compatible with new COSMIC download format)**
- icgc (update if a newer version is available from the Annovar website)

The list of databases is not written into `update_annovar.py`: it is read from the
`database_names` line of the CancerVar and InterVar configs, and each entry is compared
against the ignore-list by **exact string match**. So a database is only ignored if the
config spells its name exactly as the ignore-list does; anything else — including a
version-suffixed name such as `dbnsfp47a_interpro` — falls through into the update path.

IMPORTANT for annovar databases: since these databases are used by CancerVar and InterVar, we need to update also the corresponding CancerVar and InterVar. If CancerVar.py and InterVar.py and the corresponding configs are not provided (RECOMMENDED CASE), the tool will automatically update the CancerVar.py and InterVar.py and the relative configs in the [resources](https://github.com/zhanyinx/variantalker/tree/main/resources)

## Usage

**Note:** COSMIC update is currently disabled in both Funcotator and Annovar update scripts due to changes in the COSMIC database download format and authentication method. The scripts will skip COSMIC updates and log the current version. Manual COSMIC updates are required until the update scripts are updated to support the new format. The Funcotator GENCODE update is likewise disabled (the built database is not consumed correctly by Funcotator); the current version is retained and logged.

### Using Docker/Singularity

See [README_DOCKER.md](README_DOCKER.md) for complete Docker and Singularity usage instructions.

### Manual Execution

Both commands below are run from the `update_db` directory the clone step left you in.

Funcotator update

```bash
# Update all databases for both genome builds (default)
python3 scripts/update_funcotator.py \
  -sd /path/to/funcotator_dataSources.v1.7.20200521s \
  -gd /path/to/funcotator_dataSources.v1.7.20200521g \
  -b /path/to/backup

# Update only hg38 databases
python3 scripts/update_funcotator.py \
  -sd /path/to/funcotator_dataSources.v1.7.20200521s \
  -gd /path/to/funcotator_dataSources.v1.7.20200521g \
  -b /path/to/backup \
  --build hg38

# Update only hg19 databases
python3 scripts/update_funcotator.py \
  -sd /path/to/funcotator_dataSources.v1.7.20200521s \
  -gd /path/to/funcotator_dataSources.v1.7.20200521g \
  -b /path/to/backup \
  --build hg19
```

`-b/--backup` is not optional in practice: its default is a directory on the
institutional cluster that will not exist on your machine, and the script creates the
backup directory before it does anything else. Pass a path you own.

Both `-sd` and `-gd` are required: a run missing either one stops at argument parsing,
whichever `--build` you asked for.

**A single-build run leaves the other build alone**, in both the somatic and the germline tree.
Until issue 353 that was not true: the whole database directory was replaced, so `--build hg38`
deleted the live hg19 copy, and no later run reported it missing.

### Reading the log

The Funcotator update writes `YYYYMMDD_funcotator_update.log` into the current working directory and
**exits non-zero if any step failed**, so a scheduled run can be checked without parsing it. Each
line is one of four verdicts, defined by what happened to the live database, and every one of them
follows from an inspection rather than from the code having got that far:

| verdict | meaning |
|---|---|
| `UPDATED: <db> <build> <old> -> <new>` | the live database changed, after its staged replacement passed validation |
| `CURRENT: <db> <build> (<version>)` | a version comparison ran, found nothing newer, and left the database alone |
| `SKIPPED: <db> [<build>] - <reason>` | deliberately not attempted — disabled, no credentials, not installed |
| `FAILED:  <db> [<build>] - <reason>` | attempted and did not succeed; **the live database is untouched** |
| `RESULT:  …` | once, at the end, summarising the four counts and naming any skips |

There is deliberately no `SUCCESS` — it used to mean four different things, one of which was "did
not happen at all". A per-build `SKIPPED` is benign and does not affect the exit code; a step that
skipped every build it was asked for is a `FAILED`.

Note: `--cosmic_email` and `--cosmic_password` arguments are optional but currently not used as COSMIC update is disabled.

Annovar update

```bash
python3 scripts/update_annovar.py \
  --annovar_db_path /path/to/humandb \
  --annovar_download_script /path/to/annovar/annotate_variation.pl \
  -v /path/to/vt
```

`update_annovar.py` reads the list of Annovar databases from
`https://raw.githubusercontent.com/WGLab/doc-ANNOVAR/master/docs/user-guide/download.md`
at start-up, so it needs network access to GitHub. Override the location with
`-ad/--annovar_docs` if you need to.

It writes `YYYYMMDD_cancervar_intervar_update.log` and uses **the same four verdicts and the same
non-zero exit** as the Funcotator update above (issue 354), with two differences worth knowing:

- **The verdicts are per database, not per build.** Annovar keeps everything in one flat `humandb`
  as `<build>_<database>.txt` with the version in the filename, and the `-protocol` list inside
  `CancerVar.py`/`Intervar.py` carries no build — so a database's name advances for both builds or
  for neither. There is no `--build` option here, and nothing on this route can delete the build you
  did not ask for.
- **`FAILED` means nothing started using the new version.** The two places that name a database —
  the `-protocol` lists and `database_names` in `resources/configs/*` — are only rewritten once the
  new table has been read and found to contain data, and they are rewritten together. So a failed
  database leaves annotation running on the version already installed.

An `UPDATED` line may carry an indented `detail:` note saying that a database's `.idx` is absent or
out of date. That is not a failure: Annovar checks the index against the table's size and, when it
disagrees, scans the whole table instead and says so itself — the results are the same, the lookup
is slower.

Note: `--cosmic_email` and `--cosmic_password` arguments are optional but currently not used as COSMIC update is disabled.

**Required Tools:**
- vt software can be found here: [vt](https://github.com/atks/vt)
- GATK 4.5.0.0: [GATK releases](https://github.com/broadinstitute/gatk/releases/tag/4.5.0.0)
- Annovar: [ANNOVAR download](https://annovar.openbioinformatics.org/en/latest/user-guide/download/) (registration required)

## Tests

There is a guard suite for this tooling under [tests/](tests/). Run it from the **repository
root**, not from here:

```bash
python -m pytest update_db/tests
```

It needs `pytest` and `pandas`, and nothing else — no network, no Docker, no annotation
database. It takes about a second.

Two notes on reading its output:

- **One test is an expected `xfail`.** `test_no_update_db_script_is_invoked_through_sh` records a
  known, partly-fixed defect: three call sites in `scripts/utils_update.py`, all on the Annovar
  route, still hand a bash script to `sh`, which is fatal inside the container where `/bin/sh` is
  dash. The four on the Funcotator route are fixed. `1 xfailed` is the correct result today;
  `0 xfailed` would mean either that the last three were fixed and the marker needs deleting, or
  that the guard was deleted rather than satisfied. Run it with `--runxfail` to see the offending
  call sites listed.
- **`make test` cannot reach this suite.** That target lives in `streamlit_app/` and names an
  explicit path, so it collects the application's tests only. This suite has its own CI workflow
  (`.github/workflows/update-db-contract.yml`), which is the only one that runs from the
  repository root.

## Where these scripts came from

The per-database scripts under [funcotator/](funcotator/README.md) are derived from
GATK's own Funcotator data-source scripts. That provenance is why several of them are
disabled here rather than fixed: see [funcotator/README.md](funcotator/README.md).
