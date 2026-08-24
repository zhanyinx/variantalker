# Parity fixture set

The MAFs the parity harness (issue 17) measures `streamlit_app`'s filter against `bin/filter_variants.py`
on. **Constructed end to end** by [`build_parity_fixtures.py`](build_parity_fixtures.py), which needs
nothing mounted and runs anywhere; assembled under issue 18, replaced wholesale under issue 246.

Every number below was produced by running `bin/filter_variants.py` over the committed
files, not derived on paper.

## Why the whole set was replaced

The files this set replaces were cut from the 50-sample paired GERSOM run recorded in issue 9 and had
their barcodes pseudonymised. Issue 233 measured what that left behind, and the answer was wider
than the tickets assumed: **all seven** carried real calls, not the two named as "reference".

| file | what it was |
|---|---|
| `somatic_gnomad_genome.maf` | 6 rows byte-identical to reference rows on all 427 shared columns |
| `somatic_dot_numeric.maf` | 5 reference rows, one cell overwritten with `.` |
| `somatic_synthetic.maf` | each row a real row with **4 cells** changed — 423/427 columns still the real row's |
| `germline_synthetic.maf` | ditto — 377/380 |
| `somatic_civic.maf` | ditto — 425/427 |

The `*_synthetic` files rewrote `Start_Position` to a round number, which reads as
de-identification and is not: the same row still carried

```
Genome_Change    g.chr17:37871547C>A
cDNA_Change      c.1157C>A
Protein_Change   p.A386D
```

— the real coordinate, verbatim, in an annotation column. The variant was fully recoverable.
"Synthetic" in the old set meant *one field was edited to flip a verdict*, not *no patient data*.

That is also why the public export checks recorded **provenance** rather than content: five
of those seven passed every identifier pattern the export gate can apply, because the bytes that
must not travel are, by content, indistinguishable from the bytes that must.

## The files

| file | rows | cols | what it is |
|---|---:|---:|---|
| `somatic_reference.maf` | 82 | 427 | constructed replacement for the somatic subsample |
| `germline_reference.maf` | 94 | 380 | constructed replacement for the germline subsample |
| `somatic_synthetic.maf` | 18 | 427 | shapes the reference run could not supply |
| `germline_synthetic.maf` | 17 | 380 | ditto, germline |
| `somatic_civic.maf` | 12 | 433 | the six `CIViC_*` columns grafted on |
| `somatic_gnomad_genome.maf` | 6 | 429 | `gnomAD_genome_AF*` column-presence probe |
| `somatic_dot_numeric.maf` | 5 | 427 | `.` in a numeric column — pipeline raises, app refuses |
| `genes_somatic.txt` / `genes_germline.txt` | 12 | — | pipeline-format gene lists |
| `genes_somatic_mixed_case.txt` | 12 | — | the same genes, case mangled |

The two reference files keep the per-arm shape the subsample had — 82 rows at 427 columns and 94 at
380, each behind **111 leading comment lines** — because that shape is part of what the set
exercises. `read_maf` counts `#` lines rather than assuming a fixed offset, and `writeheader` copies
the whole block into the pipeline's output. `build_parity_fixtures.py` asserts the shape at build
time, so a row added without thinking fails where it is added.

`MANIFEST.json` carries sha256, row counts, per-row expected verdicts for the constructed fixtures,
the parameter set they were measured under, and `"provenance": "generator-constructed"` per MAF.
About 1.0 MiB total.

## Baseline under the parameter contract

`nextflow.config`'s `params` block, with issue 11's one deliberate departure (underscored
`Likely_pathogenic`). Recorded in `MANIFEST.json` under `parameters`.

| fixture | invocation | PASS | NOPASS | output cols |
|---|---|---:|---:|---:|
| `somatic_reference` | somatic | 60 | 22 | 40 |
| `somatic_reference` | somatic `--skip_pathogenic` | 54 | 28 | 40 |
| `somatic_reference` | somatic + `genes_somatic.txt` | 58 | 24 | 40 |
| `somatic_reference` | somatic + `genes_somatic_mixed_case.txt` | 58 | 24 | 40 |
| `germline_reference` | germline | 59 | 35 | 39 |
| `germline_reference` | germline `--skip_pathogenic` | 51 | 43 | 39 |
| `germline_reference` | germline + `genes_germline.txt` | 57 | 37 | 39 |
| `somatic_synthetic` | somatic | 7 | 11 | 40 |
| `germline_synthetic` | germline | 7 | 10 | 39 |
| `somatic_civic` | somatic | 8 | 4 | **46** |
| `somatic_civic` | somatic `--skip_civic` | 0 | 12 | 40 |
| `somatic_gnomad_genome` | somatic | 0 | 6 | **41** |
| `somatic_dot_numeric` | somatic | — | — | pipeline raises `TypeError`; app refuses |

The twelve rows with verdict counts are transcribed into `README_BASELINE_TABLE` in
`../../parity/test_mutation_coverage.py` and asserted against a live run, which is how the
coverage instrument shows it is driving the real harness rather than a mock of it. The
thirteenth, `somatic_dot_numeric`, is asserted there too, as a pair of non-verdicts rather
than as counts. A change to these fixtures moves both tables, together and deliberately.

The gene-list rows are the point of `genes_somatic_mixed_case.txt`: 58 PASS either way,
because both sides route through the vendored `.str.upper()`. That is divergence 9
dying by routing, made observable.

**The pass rate is higher than the subsample's, and that is what construction does.** 60 of 82
somatic rows pass at contract defaults where the sampled set passed 20, because 62 of its 82 rows
carried no witness at all. Every row here is placed on purpose. One recorded consequence:
`test_param_migration::test_expanding_the_sentinel_would_saturate_the_guideline_union` used to
assert saturation as `expanded.passed > 2 * contract.passed`, which is unreachable when the
baseline is 60 of 82 — it now asserts the mechanism instead, that the expanded criteria path
reaches *exactly* the rows clearing depth, classification and VAF. That module says so in full.

## Divergence coverage

Which of the map's twelve divergences each fixture exercises. Everything below is
constructed; the column that used to say "real" is gone with the rows it described.

| # | divergence | fixture | how |
|---|---|---|---|
| 1 | `Variant_Classification` exclude vs include | `*_reference` | **6 rows per arm** whose classification the app's hardcoded vocabulary does not contain and the pipeline's exclude list does not either, plus the three excluded values as controls |
| 2 | depth: `t_alt+t_ref` vs `DP` | `*_reference`, `*_synthetic` | 4 rows per arm where the sum and `DP` disagree across `min_depth`, both directions; `missing_DP` (PASS: the pipeline never reads `DP`), `missing_t_alt_count`, `missing_both_depth_components` |
| 3 | VAF `>` vs `>=`, NaN | `*_reference`, `*_synthetic` | one row *at* each swept threshold — germline has **6 at exactly `tumor_f == 0.2`** — plus rows between them so each sweep step moves; `vaf_exactly_threshold`, `missing_tumor_f` |
| 4 | ClinVar split vs `.isin()` | `*_reference`, `*_synthetic` | every separator the vendored splitter handles, one value each, chosen to match `CLINVAR_KEEP` **only** after the split; the three `clnsig_semicolon_*` rows |
| 5 | `All` sentinel / always-OR | — | not witnessable in a MAF: the sentinel is a UI concept that never reaches the vendored filters |
| 6 | germline ESCAT filter | `germline_reference` | **6 rows via `germline_escat_only`**: in-keep `ESCAT`, Benign on InterVar, ClinVar *and* RENOVO, clearing depth on both `t_alt+t_ref` and `DP`, VAF and classification, so the row diverges on divergence 6 and on nothing else. Plus `germline_escat_outside_keep`, which exercises the column without claiming to witness the divergence |
| 7 | somatic patho-retain | `somatic_reference`, `somatic_synthetic` | `CancerVar` Tier I/II rows at sub-threshold VAF so the rescue is what decides them, the two `CLINVAR_PATHO`-only ClinVar values, and `escat_bare_I/II/III` — the values that make the app's clause dead |
| 8 | germline patho-retain | `germline_reference` | rows RENOVO-pathogenic where InterVar is not and the other way, all at sub-threshold VAF |
| 9 | gene list case + format | `genes_*.txt` | mixed case yields identical PASS sets |
| 10 | population frequency | `*_reference` | rows spanning absent, rare and high frequency, including the common-pathogenic row the exemption exists for and one common row the exemption does **not** save |
| 11 | output columns vs `KEEP` | `somatic_civic`, `somatic_gnomad_genome` | see below |
| 12 | CIViC list-repr substring vs `.isin()` | `somatic_civic` | twelve evidence-level shapes, `list_c_and_d` and `list_a_and_e` being the witnesses |

Coverage is **measured, not listed**: `../../parity/test_mutation_coverage.py` re-injects each
divergence into the app and requires some case to notice. It reports 13/13 on this set, the same
score the replaced set reached, and catches divergences 1 and 3 on 24 cases against the old set's 11 —
a constructed witness can be placed where it is visible, while a sampled one lands where the data
put it.

### Why "the column is populated" was not coverage

Divergence 6 once read **0%** on the old fixture set while running at **51%** of the reference's
attributed diverging rows, and the divergence 6 cell was not empty the whole time. It selected a
PIK3CA `ESCAT: IA` call, Benign on InterVar, ClinVar and RENOVO — the shape exactly. Its
`DP` was **46** against a `min_depth` of **50**, so the app dropped it on divergence 2
before the ESCAT clause could admit it, the pipeline rejected it on the guideline block,
and both sides answered NOPASS. Two divergences pointing opposite ways cancelled, and a
row that witnessed nothing was counted as coverage.

The lesson generalises past divergence 6, and it is the rule the generator is written to:
**a witness cell must pin the row's whole path to the verdict, not just the field under
test.** Every other divergence has to be held neutral — depth clear on *both* sides'
columns, VAF off the boundary, classification inside the app's vocabulary and outside the
exclude list, ClinVar agreeing under split and exact matching — or some other divergence
will silence the one being witnessed. A cell that names a column rather than a mechanism
cannot make that promise.

That rule has a cost, and `overlap_rows()` pays it back: a set in which every row is excluded for
exactly one reason has no row excluded for two, and the run recap's "these counts overlap" caption
then describes something the fixture does not do. Those rows are deliberately over-determined and
witness nothing.

## What the constructed fixtures add

### The CIViC matrix

The reference run set `skip_civic: true`, so all six `CIViC_*` columns were absent from
all 100 MAFs and the two paths into `filter_variants.py:376-382` were indistinguishable.
`somatic_civic.maf` separates them:

| | `skip_civic` | `CIViC_*` present | result |
|---|---|---|---|
| A | false | yes | 46 columns, six `CIViC_*` kept |
| B | true | yes | 40 columns, stripped |
| C | false | no (`somatic_reference`) | 40 columns, stripped, `:381` warning emitted |

Fixture A's output is **exactly `KEEP + ["gnomAD_exome_AF_raw"]`** — all 45 `KEEP` entries
emitted at once. Those six CIViC columns were untested before it existed.

Every row sits on the same guideline-*failing* template, so the CIViC value alone decides
the verdict. `list_c_only` PASSes (in `filter_civic`, not in `filter_patho`'s hardcoded
`["A", "B"]`), `list_a_only` PASSes through `filter_patho`, `list_d_only` and `empty_list`
do not. Under `--skip_civic` all twelve go NOPASS.

Values are list-reprs because that is the shape the annotation writes. The two shapes the old
fixture copied verbatim from a real file are not needed: what makes `list_c_and_d` and
`list_a_and_e` witnesses is the repr convention, not their provenance.

### The synthetic supplement

Every row alters exactly one field of the neutral template, and the template is chosen so
that the alteration is what flips the verdict. Two are needed:

- a **guideline-passing** one that does *not* pass via `filter_patho` — otherwise the
  unconditional pathogenic rescue absorbs every edit and nothing is observable;
- a **guideline-failing** one, so a `CLNSIG` or `ESCAT` edit is the only thing that can
  flip it.

Under the issue 11 contract the somatic guideline trigger set is a subset of `filter_patho`'s
except for the `ESCAT` clause (`filter_cancervar` is literally the same list, and
`filter_clinvar` ⊂ `CLINVAR_PATHO`), so a usable passing template must be ESCAT-driven —
germline likewise must be `RENOVO_Class`-driven. That containment is worth carrying into
the harness: **at default parameters, somatic PASS is `filter_patho` ∪ (ESCAT-in-keep ∧
common ∧ VAF ∧ genes)**, so a depth or VAF sweep has almost no power unless
`--skip_pathogenic` is also set.

Per-row expected verdicts live in `MANIFEST.json`; all 47 named rows and 59 verdicts were
measured against the pipeline, never declared.

### Missing numerics have exactly one representable shape

The reference had zero blanks in `tumor_f` / `t_alt_count` / `t_ref_count` / `DP` across
all 181,566 rows, so the dropped-NaN-tolerance ruling was unverifiable from real data.
Measured, the blank shape matters:

| written as | `read_maf` dtype | pipeline |
|---|---|---|
| `` (empty) | `float64`, value NaN | drops the row |
| `NA` | `float64`, value NaN | drops the row |
| `.` | **`object`** | **raises `TypeError`** |

`.` is the pipeline's own blank convention for *string* columns (`ESCAT` is `.` in 88,609
of 89,019 somatic rows), so this is a shape a real MAF could plausibly carry — and it
never reaches the NaN ruling at all. The synthetic fixtures use empty; `somatic_dot_numeric.maf`
pins the pipeline's `TypeError`, so the harness can assert the app fails on the same files
rather than silently coping.

Since issue 38 the two sides no longer fail *identically*, and that is the intended state:
the pipeline dies with `TypeError: '>=' not supported between instances of 'str' and 'int'`,
naming neither the column nor the value, while the app raises `UnreadableNumericColumns`
naming all three offending columns and the `.` in each. Same non-verdict, usable message —
so `harness.ParityResult.errors_in_parity` counts it as parity and `baseline.json` records
`app_refused_columns: ["t_alt_count", "t_ref_count", "tumor_f"]`. The fourth `.` in this
fixture is in `DP`, which the pipeline never reads, and it is correctly **not** named.

Note `filter_genes` defaults to `maf["tumor_f"] > -1`, which is **False for NaN** — so a
row with a missing `tumor_f` is dropped even when no gene filter is set.

### The gnomAD_genome probe

`compute_keep`'s `gnomAD_genome_AF*` branch is untaken on the reference: no such column
exists in any of the 100 MAFs. `somatic_gnomad_genome.maf` takes the branch by column
presence alone (41 output columns instead of 40). **This is not a genome-build claim** —
the rows are hg19-shaped and the values are constructed. See the scope note.

### Six RENOVO classes, not five

`IP Pathogenic` occurred on **no row** of the old `germline_reference.maf`, a gap
`test_clinical_summary` documents and works around. A constructed set simply closes it, and each
class carries a score inside its own measured `RENOVO_pls` band — the score follows the class
rather than being set beside it, because that module measures the bands over 211,744 real rows
and fails a fixture whose score sits outside its own.

## How the set is constructed

Every cell value is written by `build_parity_fixtures.py`. The only things taken from the
replaced set live in `column_profiles.json`:

* the **column names**, in order, per fixture;
* the leading **comment-line count**;
* per column, which **blank token** the annotator uses for it (`.` or `""`) and whether its
  values are numeric.

Those are Funcotator/ANNOVAR software artifacts — the same class of metadata as the column names
themselves — and no cell value is among them. The script that derived them is kept as the record
in the private notes for issue 233; it needs the replaced files and therefore cannot be re-run
here, which is exactly the status the old generator had and the reason it is not what builds the
set.

There is no RNG and no clock is read, so two builds are byte identical and `--check` holds it.

### Every column is populated, and that is deliberate

A generator that writes only the columns it reasons about leaves the rest blank, and "the rest" is
most of the file: **387 of 427** columns on the first build, against the replaced set's 7. That is
not a smaller fixture, it is a different *shape* — one in which a column being present stops being
distinguishable from its being absent.

The distinction is load-bearing. `frequency_mask` treats a column that is present-but-blank
differently from one that is missing entirely — a blank could not sink a row, while an absent
column was never in the mask's list at all — and `config/presets.py` records the measurement where
that difference decided 17 real variants. A fixture set in which nearly every column is blank
cannot witness that, and worse, cannot fail if the app stops honouring it.

So `fill_long_tail()` gives every otherwise-empty column its own constructed values on some rows
and **its own blank token** on the others, read from `column_profiles.json`. Both conventions
therefore survive the swap, which they would not if the fill picked one. No column in any of the
seven files is entirely blank.

Columns the app *parses* are excluded from that fill by being written explicitly — the ClinVar
review-status scale, the CancerVar evidence vector, `am_class` and `RENOVO_pls` among them —
because a parsed column filled with `constructed_<name>` text would read as an unrecognised value
on every row, which is a state those scales do not have.

### The two traps, both found by building

* **The pipeline drops all-`__UNKNOWN__` columns.** `filter_variants.py:452` narrows the
  NOPASS frame with `~(out == "__UNKNOWN__").all()`. Set a *whole* KEEP column to that
  placeholder and the run dies with `KeyError: ['Matched_Norm_Sample_Barcode'] not in
  index` before reaching a verdict. `__UNKNOWN__` is kept verbatim on the germline tumour
  barcode, because issue 15's finding depends on it, and the matched-normal barcode is
  deliberately *not* uniform. The old set survived this by accident.
* **`RENOVO_pls` is not free.** See above: the score follows the class, or
  `test_clinical_summary` fails.

### The gene lists carry over verbatim

Twelve HGNC symbols per arm, unchanged. A list of gene symbols is not patient data: the symbols
were *chosen* by looking at a real run, but "APC is in this panel" says nothing about any
individual, unlike a variant call, which is quasi-identifying on its own. Keeping them verbatim
also keeps `genes_somatic_mixed_case.txt`'s exact mangling, which is irregular rather than a rule
applied to the symbols, and which `test_gene_lists` asserts on by name.

## Privacy

Nothing in this directory is derived from a patient. There are no real variant calls, no sample
barcodes, no internal paths and no coordinates carried from any run — `Genome_Change` is derived
from the coordinate the generator assigns, which is the cell the replaced fixtures leaked the real
position through.

The public export does not take that on trust. Every `.maf` in an export must be recorded in
a `MANIFEST.json` **in its own directory** as `"provenance": "generator-constructed"` and match its
recorded sha256, or the export refuses.

## What was lost, so it is not re-argued

Not coverage — the **attestation that these variant shapes occur in real data at these
frequencies**. Sentences like "6 rows per arm whose classification the app's vocabulary does not
contain" describe the app's behaviour on a shape either way; what a constructed set cannot say is
that the shape *arises*.

That evidence never lived here. It lives in the full 50-sample sweep — `make parity-reference` /
`tests/parity/reference.py`, over 181,566 rows — which needs the institutional reference drive,
has never run in CI, and stays private either way. The 176 committed rows were a witness *carrier*
for the harness, not the population evidence.

The honest residual risk is narrower: a constructed set can only hold shapes someone thought of. A
future real-data shape nobody anticipated would have had some chance of being in a purposive
subsample and has none of being here. That risk sits with the private sweep, which is where it
belongs and where it already was.

## Scope: one build, one technology

The shapes here are **hg19 + iontorrent**. Parity across genome builds and sequencing
technologies is **not claimed** by this fixture set. Specifically untested:

- `compute_keep`'s `gnomAD_genome_AF` / `_raw` branch on genuinely hg38-annotated data
  (the probe above exercises the branch by column presence, nothing more);
- ~~the two gnomAD v4 names in the app's column tail (`gnomAD_exome_AF_grpmax`,
  `gnomad41_genome_AF`), which match 0 of 100 reference MAFs~~ — **no longer applicable.**
  Issue 35 dropped both from the app's tail on exactly this evidence: they match no
  reference MAF on either build, and `compute_keep` hardcodes the older spellings, so it
  could not surface them even where they existed;
- any Illumina-specific annotation shape.

This is the map's open genome-build question, not something these fixtures settle.

## Regenerating

```bash
# rebuild the MAFs and the gene lists
python streamlit_app/tests/fixtures/parity/build_parity_fixtures.py --out DIR

# verify the committed files are what the generator produces
python streamlit_app/tests/fixtures/parity/build_parity_fixtures.py --check

# re-measure the per-row verdict oracle by running the pipeline
python streamlit_app/tests/fixtures/parity/record_manifest.py

# re-record the harness baseline, then check what the fixtures can still catch
python -m tests.parity.harness --update-baseline     # from streamlit_app/
python streamlit_app/tests/parity/mutation_coverage.py
```

The generator needs **nothing mounted**, which is what makes this set reproducible in a public
checkout rather than merely present in it. `record_manifest.py` does need `bin/`, because its
verdicts are measured rather than declared.

Running `bin/filter_variants.py` against these fixtures needs three things worth knowing:

- `--config` and `--funcotator` are **required by argparse and never read**; likewise
  `--technology`, `--annovar_version`, `--genome_build`. Any string will do.
- `filter_variants.py` imports `database.database_utils`, which imports **psycopg**. It is
  not a dependency of the app. Stub it, or invoke the vendored symbols directly.
- Each invocation needs a **fresh process**: `keep = KEEP` at `filter_variants.py:352` is
  an alias, not a copy, so one germline run mutates the module-level list and the next somatic run
  in the same process raises `KeyError: "['InterVar', 'RENOVO_Class', 'RENOVO_pls'] not in index"`.
  Reproduced here; it is the bug issue 15 requires the vendored `compute_keep` to fix.

## The fixture this replaces

`streamlit_app/tests/fixtures/test_sample.maf` is deleted. It had **no consumers** — no
test in `streamlit_app/tests/` ever loaded it — and it could not have served as a parity
fixture: it has 26 columns and no `t_alt_count` / `t_ref_count`, so `common_filters`
raises `KeyError` on it before any comparison. Its values were shapes the pipeline does
not emit: `CancerVar` as `1#Tier_I_strong Evidence` (real values are bare
`Tier_I_strong`), `RENOVO_Class` as `VUS` (not one of the six real classes), bare `A`/`B`
`CIViC_Evidence_Level` (real values are list-reprs), bare `I`/`II`/`III` `ESCAT` (absent
from all 182K reference rows), and somatic *and* germline annotation columns in one file,
which no real MAF has. Those shapes survive where they are actually informative: as the
`bare_a` / `bare_c` rows of `somatic_civic.maf` and the `escat_bare_*` rows of
`somatic_synthetic.maf`, where they now witness a divergence instead of hiding one.

**Each file here carries exactly one arm's marker columns**, for the same reason that one was
deleted: arm detection keys on column presence and never on values, so a file carrying both marker
sets is ambiguous. No `InterVar`, `RENOVO_*` or `CIViC_*` on the somatic side; no `CancerVar` or
`CIViC_*` on the germline side.
