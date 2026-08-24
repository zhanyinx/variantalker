# MAFigate 🧬

**Open an annotated MAF, filter it, and see which variants pass — and why.**

MAFigate is the interactive companion to the VarianTalker pipeline. It opens the annotated
MAF the pipeline writes, applies clinical, frequency, gene and quality filters you can
change as you go, and shows you which variants passed, which did not, and which criterion
decided it. You can inspect a variant, note what you think of it, and download the result.

It will open any MAF-format file, not only this pipeline's output — and when a column a
filter needs is missing, it tells you instead of filtering as though the column were there.

![Python](https://img.shields.io/badge/Python-3.11+-3776AB?style=flat-square&logo=python)
![Streamlit](https://img.shields.io/badge/Interface-Streamlit-FF6B6B?style=flat-square)

---

## What you can do

### Work on either arm

- **Somatic** — CancerVar tiers, CIViC evidence levels, ESCAT actionability
- **Germline** — InterVar and ReNOVo classifications

### Filter

- **Clinical significance** — keep or drop pathogenic, likely pathogenic and uncertain
  (VUS) calls, per classifier
- **Variant type** — exclude classes you do not want to see, such as silent, intronic,
  IGR and RNA
- **Population frequency** — an upper frequency bound, read across gnomAD (exome and
  genome), ExAC, ESP6500, 1000 Genomes and the file's own maximum-frequency column
- **Genes** — restrict to one of the panels below, or to your own gene list
- **Sequencing quality** — minimum read depth and variant allele frequency
- A variant a classifier already calls pathogenic is kept **even when it fails the other
  cuts**, so a known finding is not lost on coverage, frequency or panel membership:
  CancerVar Tier I or II, CIViC evidence A or B, or a pathogenic ClinVar classification on
  the somatic arm; InterVar pathogenic or likely pathogenic, or ClinVar, on the germline
  arm. **Auto-retain pathogenic variants** turns this off.

Four starting points — **Broad Somatic**, **Broad Germline**, **Stringent Somatic**,
**Stringent Germline** — set every threshold at once; adjust any of them afterwards. Broad
casts a wide net and keeps uncertain calls for review; Stringent keeps a short list.

### Read the result

- Passed and failed variants side by side, each with its count
- Why the report came out as it did, broken down by criterion
- Six charts: VAF distribution, variant classification, chromosome distribution, top genes,
  clinical significance passed against failed, and mutation type spectrum
- A per-variant view with the ACMG evidence behind the call, AlphaMissense, the in-silico
  predictors (REVEL, BayesDel, VEST4, MutPred2, CADD) where the MAF carries them, and links
  out to ClinVar, gnomAD, dbSNP and the UCSC browser
- A **Notes** column you can type into, and custom annotation columns you can add

### Download it

- Passed variants and failed variants as CSV, plus a summary of how the report was reached
- Your current filter settings as JSON or YAML, ready to upload again next time
  (`example_parameters.json` and `example_parameters.yaml` show the shape)

---

## How to get it

Two routes, and right now they are not equivalent: **clone + `setup.sh` works today**, and
the desktop installers are still being built. The table below is the canonical statement of
what each route asks of you and what happens to your file on it — everything else on this
page, and elsewhere in the repo, points here rather than repeating it.

|  | **Clone + `setup.sh`** | **Desktop installers** (`.dmg`, `.exe`) |
| --- | --- | --- |
| **Status** | **Works today.** The route described below, and the only one there is at the moment. | **Not built yet** — coming with the first `mafigate-v1.0.0` release. There is nothing to download until then, and no release page to send you to. |
| **Who it is for** | Collaborators on **Linux**, and anyone extending the code. Linux is not a fallback here: neither installer will target it, so this is the Linux route for good. | Everyone else — clinicians and the molecular tumour board, readers of the paper, and collaborators on macOS or Windows who would rather not open a terminal. Either way you will need an annotated MAF from the pipeline before there is anything to look at. |
| **What you get** | The full app, with no reduction — and the source, which is the other reason to choose this route. | The full app, with no reduction. The same code, packaged. |
| **What you need** | `git`, a `python3`, and a POSIX shell. Linux and macOS have all three; Windows only through a bash (Git Bash or WSL). Internet once, while `setup.sh` fills the virtual environment it builds inside the checkout — your own Python is used to build it and is not otherwise touched. | Nothing, once they ship: each brings its own Python — both the same pinned release. Internet once, on first launch, while it builds its environment. |
| **Where your file goes** | **Your MAF: nowhere — nothing leaves your machine.** It is read from your own disk, and MAFigate makes no network request of its own — the ClinVar and gnomAD links in the variant view open in your browser, they are not lookups the app performs. Streamlit's own anonymous usage reporting is off as well: the `.streamlit/config.toml` this repo ships sits beside the app, where every launch route resolves it. | Will be identical to the clone route: the same code, reading your file on your own machine, with the same config carried into the bundle. |
| **When it goes wrong** | A `python3` that cannot build a virtual environment: Debian and Ubuntu ship that separately, as `python3-venv`, and `setup.sh` stops and names it. Little else — both scripts install into and launch from that one environment, so `pip` and `streamlit` resolving to different interpreters, which used to start the app with pieces missing rather than failing outright, cannot happen. | Unsigned artifacts will be blocked on first open until you allow them, and no internet on first launch means the install completes but the app does not start. |

### Clone and run it

> **Prerequisites.** `git`, a **`python3` (3.11 or newer) with `pip`**, and
> a POSIX shell — the scripts below install MAFigate's dependencies, not Python itself. On
> Windows, run them from a bash: Git Bash or WSL.

```bash
git clone https://github.com/zhanyinx/variantalker.git
cd variantalker/streamlit_app

./setup.sh            # builds .venv/ and installs into it
./run_mafigate.sh     # opens the app at http://127.0.0.1:8501, on this machine only
```

`setup.sh` builds a **virtual environment at `.venv/` inside the checkout** with your
`python3`, installs `requirements.txt` into it, and then checks that what it asked for
actually arrived. `run_mafigate.sh` launches out of that same environment and refuses to
start if anything is missing. Both name the interpreter they are using as they start, and
neither **installs** anything outside the checkout — installing MAFigate cannot disturb the
Python your other work depends on. (`pip` still caches the downloads it makes in your home
directory, as it does for any install; nothing there belongs to an interpreter.) Re-run
`setup.sh` any time: it reuses the environment, and repairs it if the Python that built it
has moved on. Delete `.venv/` to start over.

Two details are deliberate, and each has a reason worth a line.

**One interpreter, reached as `-m` modules.** `pip` and `streamlit` are separate entries on
your `PATH` and need not belong to the same Python — on a machine with both Homebrew and
conda they usually do not. When they disagree, `pip install` succeeds into an interpreter the
app never runs in, nothing fails, and MAFigate starts up missing pieces of itself: the
symptom that prompted this was six Summary charts quietly drawing the wrong thing
(issue 162). `python -m pip` and `python -m streamlit` cannot disagree that way, and both
scripts read where the environment is from the same file, so they cannot end up in different
ones.

**In the checkout, not in your home directory.** The desktop installers keep their
environment in a shared `~/.mafigate`; two clones sharing one environment would fight over
it, so the clone route keeps its own beside the code. It is git-ignored.

To use an interpreter of your own instead — a conda environment, a specific version, a venv
you keep elsewhere — name it, and no environment is built or consulted here:

```bash
MAFIGATE_PYTHON=/path/to/python ./setup.sh
MAFIGATE_PYTHON=/path/to/python ./run_mafigate.sh
```

To check an interpreter without installing or launching anything. These two do not check
the same one, which is the point of running either: each names the interpreter it looked
at in its first line of output.

```bash
make check-deps                # the interpreter make would run the app with
python3 check_dependencies.py  # whichever python3 is first on your PATH
```

### In an environment you manage yourself

The two commands above already give you an isolated environment. If you would rather it be
one of yours — conda, pyenv, uv, a venv somewhere else — build it, and hand it to both
scripts:

```bash
conda create -n mafigate python=3.11 -y && conda activate mafigate
# ... or any other environment you keep, activated or not

MAFIGATE_PYTHON=$(command -v python) ./setup.sh
MAFIGATE_PYTHON=$(command -v python) ./run_mafigate.sh
```

Or without the scripts at all — **with that environment activated**, so that `python` is
its interpreter rather than another one:

```bash
python -m pip install -r requirements.txt
python -m streamlit run MAFigate.py --server.address 127.0.0.1
```

`python -m` rather than bare `pip` and `streamlit` for the reason given above: inside an
activated environment they are usually the same thing, and on the day they are not, `-m`
is the spelling that cannot be wrong. Outside one, `python` is whatever your `PATH` says,
and on a machine with conda installed that is usually conda's base environment — so this
pair of commands would fill the interpreter your other work depends on rather than an
isolated one, which is the whole thing `setup.sh` exists to avoid. `command -v python`
first, if you are not certain which you have.

`--server.address 127.0.0.1` because `streamlit run` without it binds **every** interface,
and this route does not get `run_mafigate.sh`'s flags for free. Left off, anyone who can
route to this machine on that port can drive MAFigate and open your MAF through it.

### What it needs

- **Operating system** — Linux, macOS, or Windows through a bash. The app itself runs
  anywhere Python and a browser do; it is the two shell scripts above that want a POSIX
  shell.
- **Browser** — any current Chrome, Firefox, Safari or Edge
- **Memory** — the whole file is held in memory while you work on it: roughly 140 MB plus
  ten times the file's size on disk, so a 100 MB MAF wants a little over 1 GB free. With
  too little, the process is killed and says nothing about why.
- **File size** — files arrive through the browser, where the upload limit is 200 MB by
  default. For a larger MAF, raise it at launch. `run_mafigate.sh` passes no extra flags
  through, so this is a direct launch and has to pin the address itself:

  ```bash
  .venv/bin/python -m streamlit run MAFigate.py \
      --server.address 127.0.0.1 --server.maxUploadSize 1024
  ```

---

## Using it

1. **Configure Parameters** — choose the arm, pick **Broad** or **Stringent** as a starting
   point, and adjust the thresholds. You can skip this and come back: the app opens with a
   working set of parameters.
2. **Open your MAF** — choose a MAF, TSV or TXT file in the **sidebar**, which offers the
   chooser on every page. The filters run on it immediately and the **Results** section
   opens with the report.
3. **Results** — sort and filter the columns, click a variant to see it in full, write
   notes, and download what you need.
4. Changed a threshold? **Re-apply Filters**, on the Load section, refreshes the report
   against the new settings.

Opening a different file replaces the one you had, and clears the notes written against it —
notes belong to the file they were written for. Download the table if you want to keep them.

The sidebar names the file you have open, how many variants it holds, how many passed the
current filters and which arm you are on — from any page, with a button back to the report.

---

## Annotations it reads

MAFigate does not query these resources. It reads the columns the pipeline has already
written into the MAF, which is why a MAF annotated without one of them will report that
column as missing.

| Annotation    | What it contributes            | Arm      |
| ------------- | ------------------------------ | -------- |
| **CancerVar** | Somatic pathogenicity tier     | Somatic  |
| **CIViC**     | Clinical evidence level        | Somatic  |
| **ESCAT**     | Actionability tier             | Somatic  |
| **InterVar**  | Germline classification        | Germline |
| **ReNOVo**    | Germline classification        | Germline |
| **ClinVar**   | Reported clinical significance | Both     |

These six are what the filters read. The per-variant view shows more where the MAF carries
it — AlphaMissense, the dbNSFP in-silico predictors and the ACMG evidence codes — and says
so rather than inventing a value when a column is absent.

The **Clinical Summary** column combines the six into one verdict per variant and says which
source it came from; see [docs/CLINICAL_SUMMARY_FEATURE.md](docs/CLINICAL_SUMMARY_FEATURE.md).

## Gene sets

Restrict the analysis to a panel or database gene set, or supply your own list:

- **MSK-IMPACT** / **MSK-IMPACT Heme** — MSK-IMPACT solid-tumour / haematologic panel genes
- **FoundationOne** / **FoundationOne Heme** — FoundationOne CDx / Heme panel genes
- **COSMIC** — COSMIC Cancer Gene Census genes
- **OncoKB** — OncoKB actionable cancer genes
- **Custom** — your own gene list: one symbol per line, or separated by commas,
  semicolons or spaces

> A set left empty behaves like "All" — no gene filtering. The lists themselves live in
> `config/gene_panels.py`.

---

## If something goes wrong

- **Help & Documentation** in the app covers the parameters, the columns and the common
  questions in more detail than this page.
- **Port already in use** — launch on another one, through the environment rather than a
  bare `streamlit`, for the reason under **One interpreter, reached as `-m` modules**
  above:
  `.venv/bin/python -m streamlit run MAFigate.py --server.address 127.0.0.1 --server.port 8502`
- **A filter seems to do nothing** — check the missing-column warnings above the report.
  A filter whose column the MAF does not carry cannot cut anything.
- **Anything else** — please
  [open an issue](https://github.com/zhanyinx/variantalker/issues). Read
  [CONTRIBUTING.md](../CONTRIBUTING.md) first if you want to know what happens to it
  afterwards; the answer is unusual and worth two minutes.

---

## For developers

`make help` lists the common targets — not quite all of them: `test-fast`, `qa`, `dev` and
`build` are in the Makefile and absent from the help text. The ones you will want:

```bash
make install          # install into .venv/ if ./setup.sh built one — it does not build it
make run              # start the app — binds every interface, unlike ./run_mafigate.sh
make test             # the test suite (pytest) — needs pytest installed; see below
make app-load-check   # boot the app and load a MAF through both load paths
make format           # black + isort (config in pyproject.toml) — see below
```

`make app-load-check` needs nothing beyond `requirements.txt`, and is the quickest evidence
that a change has not broken the app: it boots MAFigate twice through Streamlit's own test
harness — once per load path — and reports what each one loaded, filtered and drew. It runs
`tests/run_app_check.py`, which is a script rather than a pytest test, so `make test` does not
collect it and neither installer ships it. **Every command on this page works from a clone of
the public repository** — this one did not until issue 345, because the script it ran used to
live under `docs/`, which the public export strips.

**The development tools are not installed by `setup.sh`.** `requirements.txt` pins what the
app *runs* on and deliberately nothing else, so a fresh checkout has no pytest, no black, no
isort and no flake8 — and the targets that want them do not fail alike:

- `make test` fails outright, with `No module named pytest`. Install it into the environment
  `setup.sh` built: `.venv/bin/python -m pip install pytest`.
- `make format` and `make lint` print a warning, do nothing, and **exit 0** — so they are
  green whether or not they ran, and a clean `make format` is not evidence the tree is
  formatted. They also look for `black`, `isort` and `flake8` on your `PATH` rather than in
  `.venv/`, so installing them into that environment is not enough on its own: activate it,
  or call the tool directly, as `.venv/bin/black`.

  Before you let either of them write, know that this tree predates the current black:
  `black . --line-length 100 --check` reformats **101 of 133 files** on black 26, none of it
  to do with your change. Format the files you touched and no others —
  `.venv/bin/black <your files> --line-length 100` — or pin the version the tree was last
  formatted with. `--check --diff` first, always.

The filtering code under `vendor/` is a byte-for-byte copy of the pipeline's own, kept
honest by a drift guard — **read [vendor/README.md](vendor/README.md) before changing
anything near filtering**, and see `make check-vendor` and `make check-params`.

Standalone desktop installers are built from
[build/BUILD_INSTRUCTIONS.md](build/BUILD_INSTRUCTIONS.md), which is also where signing and
notarization are described. What each route offers a reader, and what it costs them, is the
table at the top of this page. This line restates none of it on purpose: a status written
down twice is a status that goes stale in one of the two places the day it changes.

Changes arrive here by a route worth knowing before you send one:
[CONTRIBUTING.md](../CONTRIBUTING.md).

---

## License

MAFigate is part of VarianTalker and shares its terms: free for non-commercial use, without
warranty. Some of the tools the annotations come from — annovar, for instance — must be
licensed separately by you. Please contact the authors for commercial use. See
[LICENSE](../LICENSE).
