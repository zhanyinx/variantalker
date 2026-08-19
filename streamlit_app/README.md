# MAFigate 🧬

**Open an annotated MAF, filter it, and see which variants pass — and why.**

MAFigate is the interactive companion to the VarianTalker pipeline. It opens the annotated
MAF the pipeline writes, applies clinical, frequency, gene and quality filters you can
change as you go, and shows you which variants passed, which did not, and which criterion
decided it. You can inspect a variant, note what you think of it, and download the result.

It will open any MAF-format file, not only this pipeline's output — and when a column a
filter needs is missing, it tells you instead of filtering as though the column were there.

![Python](https://img.shields.io/badge/Python-3.9+-3776AB?style=flat-square&logo=python)
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
| **What you need** | `git`, a `python3`, and a POSIX shell. Linux and macOS have all three; Windows only through a bash (Git Bash or WSL). | Nothing, once they ship: each will bring its own Python. Internet once, on first launch, while it builds its environment. (Today's unreleased Windows build still asks for a Python of your own; that is part of what is not finished.) |
| **Where your file goes** | **Your MAF: nowhere.** It is read from your own disk, and MAFigate makes no network request of its own — the ClinVar and gnomAD links in the variant view open in your browser, they are not lookups the app performs. One qualification, and it is Streamlit's rather than MAFigate's: Streamlit reports anonymous usage statistics by default, and this repo does not yet ship a config that turns them off. Until it does, turn them off yourself — `browser.gatherUsageStats = false` in a `.streamlit/config.toml`. | Will be identical to the clone route: the same code, reading your file on your own machine. |
| **When it goes wrong** | `pip` and `streamlit` resolving to different interpreters — the app then starts with pieces missing rather than failing outright. Both scripts below go through one interpreter so that it cannot happen. | Unsigned artifacts will be blocked on first open until you allow them, and no internet on first launch means the install completes but the app does not start. |

### Clone and run it

> **Prerequisites.** `git`, a **`python3` (3.9 or newer, 3.11 recommended) with `pip`**, and
> a POSIX shell — the scripts below install MAFigate's dependencies, not Python itself. On
> Windows, run them from a bash: Git Bash or WSL.

```bash
git clone https://github.com/zhanyinx/variantalker.git
cd variantalker/streamlit_app

./setup.sh            # installs the dependencies (equivalently: make install)
./run_mafigate.sh     # opens the app at http://127.0.0.1:8501, on this machine only
```

Both scripts work on **one interpreter**, `python3` unless you say otherwise, and both name
it as they start. `setup.sh` installs with `python3 -m pip` and then checks that what
`requirements.txt` asked for actually arrived; `run_mafigate.sh` launches with
`python3 -m streamlit` and refuses to start if anything is missing.

That is deliberate, and the reason is worth a line. `pip` and `streamlit` are separate
entries on your `PATH` and need not belong to the same Python — on a machine with both
Homebrew and conda they usually do not. When they disagree, `pip install` succeeds into an
interpreter the app never runs in, nothing fails, and MAFigate starts up missing pieces of
itself: the symptom that prompted this was six Summary charts quietly drawing the wrong
thing (issue #162).

To use a different interpreter — a virtual environment, a conda environment, a specific
version — point both scripts at it:

```bash
MAFIGATE_PYTHON=/path/to/python ./setup.sh
MAFIGATE_PYTHON=/path/to/python ./run_mafigate.sh
```

To check an interpreter without installing or launching anything:

```bash
make check-deps                                    # or: python3 check_dependencies.py
```

### In an isolated environment

```bash
git clone https://github.com/zhanyinx/variantalker.git
cd variantalker/streamlit_app

python3 -m venv .venv && source .venv/bin/activate
# ... or: conda create -n mafigate python=3.11 -y && conda activate mafigate

python -m pip install -r requirements.txt
python -m streamlit run MAFigate.py
```

`python -m` rather than bare `pip` and `streamlit` for the reason given above: inside an
activated environment they are usually the same thing, and on the day they are not, `-m`
is the spelling that cannot be wrong.

### What it needs

- **Operating system** — Linux, macOS, or Windows through a bash. The app itself runs
  anywhere Python and a browser do; it is the two shell scripts above that want a POSIX
  shell.
- **Browser** — any current Chrome, Firefox, Safari or Edge
- **Memory** — the whole file is held in memory while you work on it: roughly 140 MB plus
  ten times the file's size on disk, so a 100 MB MAF wants a little over 1 GB free. With
  too little, the process is killed and says nothing about why.
- **File size** — files arrive through the browser, where the upload limit is 200 MB by
  default. For a larger MAF, raise it at launch:
  `python3 -m streamlit run MAFigate.py --server.maxUploadSize 1024`

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
- **Port already in use** — launch on another one:
  `streamlit run MAFigate.py --server.port 8502`
- **A filter seems to do nothing** — check the missing-column warnings above the report.
  A filter whose column the MAF does not carry cannot cut anything.
- **Anything else** — please
  [open an issue](https://github.com/zhanyinx/variantalker/issues). Read
  [CONTRIBUTING.md](../CONTRIBUTING.md) first if you want to know what happens to it
  afterwards; the answer is unusual and worth two minutes.

---

## For developers

`make help` lists everything; the ones you will want:

```bash
make install          # pip install -r requirements.txt
make run              # start the app — binds every interface, unlike ./run_mafigate.sh
make test             # the test suite (pytest)
make app-load-check   # boot the app and load a MAF through both load paths
make format           # black + isort (config in pyproject.toml)
```

The filtering code under `vendor/` is a byte-for-byte copy of the pipeline's own, kept
honest by a drift guard — **read [vendor/README.md](vendor/README.md) before changing
anything near filtering**, and see `make check-vendor` and `make check-params`.

Standalone desktop installers are built from
[build/BUILD_INSTRUCTIONS.md](build/BUILD_INSTRUCTIONS.md). They are ad-hoc-signed, neither
notarized nor released yet, and what each route offers today is the table at the top of this
page rather than anything restated here.

Changes arrive here by a route worth knowing before you send one:
[CONTRIBUTING.md](../CONTRIBUTING.md).

---

## License

MAFigate is part of VarianTalker and shares its terms: free for non-commercial use, without
warranty. Some of the tools the annotations come from — annovar, for instance — must be
licensed separately by you. Please contact the authors for commercial use. See
[LICENSE](../LICENSE).
