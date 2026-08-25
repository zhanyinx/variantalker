# Changelog

Every release of **variantalker**, the Nextflow pipeline in this repository — newest first,
one section per version, saying what changed under that number.

This file exists because the pipeline prints a version number on the first screen of every
run and there was no written account of it anywhere. A reader who saw a number on that banner
had no way to learn what it succeeded, what moved, or whether it meant anything at all — and
for the two releases of the old shared line the honest answer was that it did not:
`manifest.version` read the same string at both of them, so neither bare number ever reached
the software that carried it. A number nobody can account for is worse than no number, and
this is the account.

It lives **in the tree** rather than only on a release page so that it travels: a clone reads
its own history without needing a web page, and there is one text to keep true instead of two
pages saying the same thing in slightly different words.

**It covers the pipeline only.** MAFigate — the Streamlit application under `streamlit_app/` —
is versioned separately and keeps its own record in `streamlit_app/build/RELEASES.md`. The two
numbers are independent: moving one never implies moving the other, and neither is a statement
about the other. That independence is stated here out loud because its *absence* is what went
wrong before — see the note on the shared line below. Changes to the application are not
recorded in this file.

**What holds this file honest: the release procedure, not CI.** The versioning contract is
written where the number lives — beside `manifest.version`, in `nextflow.config`'s manifest
block — and the guard it names, `.github/workflows/version-contract.yml`, refuses a release tag
that disagrees with the manifest, a tree sitting on a version already released, and any second
file stating a pipeline version. But no guard reads *this* file: keeping it current is a step
of the release act, and a release that skips that step leaves a stale record nothing will flag.
Said plainly rather than left to be discovered, because the neighbouring record in
`streamlit_app/build/RELEASES.md` *is* read by a guard, and the difference between the two
would otherwise look like an oversight instead of a chosen design.

## A note on the shared line: `1.0` and `1.1`

The `1.0` and `1.1` entries below are **not pipeline releases.** They are this repository's
former *shared* release line: one bare number covering the pipeline and the Streamlit
application together, cut before the two products were given separate identities. `1.1`'s own
notes say so on their face — one item about Annovar-versus-Funcotator annotation, one about the
application's design and interface, under a single number.

They are recorded here as history, unchanged and unrenamed. They were deliberately not
retagged into the pipeline's namespace, because that would assert a pipeline version those
releases never carried: the manifest read `1.0.0` at both. Re-labelling them as pipeline
releases is precisely the ambiguity this file was written to end, so the shape is left visible
instead of tidied away.

The pipeline's own line starts at `1.0.0` and names the pipeline in its own tag; the shared
line stops at `1.1`. **The two lines therefore both contain a `1.0`, and that is deliberate.**
The dev's call, taken when the first pipeline release was cut: this release is the version the
manuscript describes, so it is the pipeline's `1.0.0` rather than a continuation of a numbering
the pipeline never really had. The namespaces keep the two apart where it matters — the shared
line's tags were bare (`v1.0`, `v1.1`, since deleted), this line's are `variantalker-v1.0.0` —
so nothing collides, and the earlier entries stay below as history under their own numbers. An
earlier draft of this file argued the opposite, that re-using `1.0` would assert a pipeline
version those releases never carried. That reasoning is recorded here rather than deleted,
because it is the honest account of what was weighed: the ambiguity is real, and it was
accepted knowingly in exchange for a first release whose number means something to a reader of
the paper.

## Versions

### 1.0.0 — 2026-08-25

The first release cut as the pipeline's own, rather than as half of the shared line, and the
version the manuscript describes. The number restarts the count deliberately — see the note on
the shared line above.

**Read the scope of this entry before trusting it.** It is derived from the commit range since
`1.1` — 32 commits touching the pipeline's own code between 2026-05-19 and 2026-08-21 — and
every item below was checked against the tree rather than taken from a commit message. It is
**not** the result of a compatibility audit: nothing here has been tested for whether it breaks
an existing invocation, and no part of the number should be read as a finding of compatibility.
Where a change alters behaviour for an unchanged command line, it says so.

#### New and changed parameters

- The ClinVar and ESCAT terms used by the **guidelines** filter are no longer hardcoded:
  `filter_clinvar` and `filter_escat` are settable in a config file or directly on the command
  line, each a comma-separated OR filter. **This changes which variants that filter keeps by
  default**, and not only by making the list editable — the previous code enumerated every
  *composite* ClinVar significance string it accepted (`Pathogenic/Likely_pathogenic`,
  `Pathogenic|drug_response`, and dozens more), whereas the filter now splits a ClinVar cell
  into its individual calls and keeps the variant if any one of them is in the list. See the
  matching note under *Annotation*; the two are one change.

  Note what this did **not** change: the separate pathogenic-retention rule — the one that
  keeps a pathogenic variant regardless of the gene and VAF filters, unless
  `skip_pathogenic_retention` is set — still reads a built-in ClinVar term list and is not
  affected by `filter_clinvar`. Two rules read the same `ClinVar_VCF_CLNSIG` column and only
  one of them is now configurable.
- `skip_pharmgkb` skips PharmGKB annotation.
- `funcotator_ncores` sets how many parts each Funcotator chunk is split into for parallel
  annotation, independently of the overall `chunk_size`. At its default of `1` the chunk is
  annotated by a single process and no splitting happens at all.

#### Input

- The samplesheet accepts an optional **`custom_id`** column. Where it is set it replaces the
  patient name as the sample's output directory; the patient identifier is still what must
  match the VCF. Two rows sharing a `custom_id` are rejected.

#### Annotation

- Funcotator annotation moved into its own subworkflow, with a separate process for the split
  and unsplit cases, so a run that asks for no splitting no longer pays for the chunking
  machinery.
- A **balanced chunker** distributes variants evenly across chunks rather than splitting them
  unevenly, which is what made runs with few samples slow.
- VCF standardisation was reorganised to streamline the path into annotation.
- Where a ClinVar cell records more than one clinical call, the pipeline now reads it the way
  the application does: the cell is split on `|`, `/`, `;` and `,` and the variant is kept if
  any single parsed call is one being filtered for. The two no longer disagree about the same
  cell. **This changes which variants are kept** for every multi-call ClinVar cell.
- Duplicated CancerVar/InterVar annotations are reduced to the first occurrence. There is no
  information in the record on which to make an informed choice between them.
- PharmGKB annotation is restricted to germline samples.

#### Fixes

- `chrM` is translated to `MT` for Funcotator.
- Complex variants are handled correctly when the annotation database query is split.
- `PASS`-variant filtering, which was missing when input was split into chunks, is applied.
- An empty VCF, or one whose variants are all already cached, no longer fails the run.
- MAF header content and chromosome ordering are corrected.
- `*` placeholder ALT alleles are filtered during VCF standardisation.
- Germline samples are queried with the germline tissue type rather than the sample's own.
- Duplicate process tags no longer crash a run with a "multiple instances of the same file"
  error.
- `clean_final_maf.py` is faster, and pipeline resource requests were adjusted.

#### Database tooling (`update_db/`)

- The Funcotator and Annovar download routes **fail loudly instead of continuing past an
  error**, and say what failed. Previously a failed download could leave a partial database
  behind and report success.
- The manual route is documented as **Linux-only** and the claim is guarded. The databases
  require roughly 500 GB, so this tooling targets servers and HPC rather than a laptop.

#### Documentation and self-presentation

- The repository URLs the pipeline prints are built from `manifest.homePage` rather than
  assembled from the pipeline's display name.
- The pipeline presents itself as nf-core-**derived** rather than as an nf-core pipeline,
  which it is not.
- The `update_db`, Docker, vendor, biomarker and utility READMEs were walked against what
  actually runs, and corrected where they disagreed with it.
- The version number itself is now under a written contract: `manifest.version` names the
  release being prepared, its terms sit beside the value in `nextflow.config`, and a CI
  guard (`.github/workflows/version-contract.yml`) refuses a release tag, an already-used
  number, or a second in-tree statement of the pipeline's version. Before this, the manifest
  had read the same value at every release ever published.

### 1.1 — 2026-05-05 — shared line

Released as *"Stable version"*. Covered the pipeline **and** the Streamlit application under
one number. Its notes, in full:

> 1. Fix minor bugs for complex variants in annovar vs funcotator annotation
> 2. major update on streamlit design and interface

`manifest.version` read `1.0.0` at this release.

### 1.0 — 2026-01-09 — shared line

Released as *"First release"*. Its notes, in full:

> Functional release supporting databases updated in 2026

`manifest.version` read `1.0.0` at this release — the only one of the three numbers here that
the software ever actually printed for the release it named.
