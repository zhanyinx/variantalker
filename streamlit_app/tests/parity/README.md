# The parity harness

The measurement instrument the map's destination is defined in terms of (#17). It runs
`bin/filter_variants.py` and the app's filter path over the same MAF with equivalent parameters,
diffs the PASS sets and the output column sets, and **attributes every diverging row** to one of the
map's known divergences.

```bash
cd streamlit_app

python -m pytest tests/parity/ -q                     # assert against the baseline
python tests/parity/harness.py --report               # print the baseline table
python tests/parity/harness.py --case somatic_defaults
python tests/parity/harness.py --update-baseline      # re-record deliberately

python tests/parity/mutation_coverage.py              # what the fixtures can still catch
```

Full matrix: 33 cases, ~10 s; the mutation sweep is another ~11 s, so `tests/parity/`
totals ~31 s for 252 tests. Everything that runs
the pipeline skips when `bin/` is absent — `test_baseline_integrity.py` and
`test_loader_premise.py` are what still hold in that checkout, which is why this is not
the whole story. `../README.md` is: it names every instrument in the suite, what each
owns, and which of them still assert the filter path with no pipeline in it.

| file | what it is |
|---|---|
| `contract.py` | the #11 parameter contract, its translation to each side, and the case matrix |
| `harness.py` | the runners, the diff, the attribution probes, the report |
| `baseline.json` | what the harness measures on the app **as it stands** |
| `test_parity.py` | baseline assertions + parity assertions, all green and none marked |
| `test_absent_columns.py` | the **one deliberate deviation**: a MAF missing a filter-input column, which the pipeline refuses and the app fills (#39) |
| `mutation_coverage.py` | the coverage instrument: re-injects each divergence into the app and records which cases notice (#242) |
| `test_mutation_coverage.py` | 13 of 13 divergences caught, the fixture README's table reproduced, and the instrument shown reporting a gap |
| `test_baseline_integrity.py` | what `baseline.json` must say, with **no `bin/`**: no case diverges, no diverging row is unattributed |

## How the two sides are run

**Pipeline** — a fresh subprocess per case. Fresh is not optional: `keep = KEEP` at
`filter_variants.py:352` is an alias, so one germline run mutates the module-level list
and the next somatic run in the same process emits germline columns or raises. `psycopg`
is stubbed, because `filter_variants.py:10` imports `database.database_utils` which
imports it at module scope; nothing on the filtering path touches a database.

`--config`, `--funcotator`, `--technology`, `--annovar_version` and `--genome_build` are
**required by argparse and never read** — passed as fixed placeholders and deliberately
kept out of the contract, so they cannot be mistaken for parameters under test.

**App** — `filters.variant_filters.apply_filters`, in process.

Since issue #33 that function contains no filter logic: it calls the same vendored
functions the pipeline side calls. So what this harness now measures is not "do two
implementations agree?" but "does the app's plumbing around one implementation preserve
its verdict?" — parameter reading, the gene-list-to-file adapter, and the row labelling
are all still capable of losing a verdict, which is why the pipeline subprocess stays.

**Both are handed the same DataFrame**, loaded by the vendored `pipeline_utils.read_maf`.
Feeding one frame to both sides is what makes every divergence below a divergence in
*filter logic* rather than in I/O.

That started as a workaround. The app used to load MAFs itself, building its frame from
`line.split("\t")`, so every column was a Python string and the pipeline's `common_filters`
*concatenated* `t_alt_count` with `t_ref_count` and raised `TypeError: '>=' not supported between
instances of 'str' and 'int'` (#10). The harness therefore had to do its own loading, and could not
be read as validating the app's. #32 deleted that loader — `utils.read_maf` now delegates to the
same vendored reader — so the two are the same frame by construction, asserted per fixture by
`test_harness_loads_the_frame_the_app_would_load`.

## The comparison key

`Chromosome` / `Start_Position` / `Reference_Allele` / `Tumor_Seq_Allele2`, asserted unique on all
seven MAF fixtures by `test_comparison_key_is_unique`. It cannot be extended with
`Tumor_Sample_Barcode`, which is `__UNKNOWN__` in every germline reference row; #18's generator
deduplicates on the key while pooling, which is what makes the property hold.

## Translating one contract into two parameter sets

Contract values are `nextflow.config`'s `params` block, with #11's single departure
(underscored `Likely_pathogenic`). Issue #33 made the app take the pipeline's own parameter
names, spellings and polarity, so the translation is now nearly the identity — the table
below is what is *left* of it, and none of the three is a divergence:

| contract value | app side |
|---|---|
| `vaf_threshold` | one key, holding this arm's value. The pipeline takes both thresholds on every command line and picks one by arm; the app knows its arm. |
| `filter_escat` | somatic only. `germline_filters()` takes no ESCAT argument, and a germline ESCAT clause was the largest divergence on real data. |
| `genes` (a *file*) | a list of symbols, or — in the `lines`/`spaces`/`semicolons` modes issue #34 added — one **pasted string** the app must tokenise. Its gene adapter writes the symbols back out to a file for the vendored clause, so both sides read the same symbols. `contract.app_genes` records which paste separators can witness the tokenising bug and which cannot. |

Two entries this table used to carry are gone rather than changed: `skip_civic` (the app
had no equivalent at all) and `skip_pathogenic` (inverted into `keep_pathogenic`). Both are
now the pipeline's own parameter, passed straight through.

Divergence #1 needs two readings, because the app's parameter of the same name means the
opposite thing, so each case runs in one `vc_mode`:

- **`literal`** — the same list value to both, so the app keeps exactly what the pipeline
  drops. Documents #1 at full size. Only the two `*_defaults_literal` cases use it.
- **`complement`** — the app gets its vocabulary minus the list, the closest expressible
  equivalent. Used everywhere else, so the other divergences are visible instead of
  being swamped by this one. It is not exact, and the residue is the point:
  `#1 vc_outside_app_vocabulary` counts rows the app's hardcoded vocabulary cannot
  express at all.

## What is asserted

Two layers, and the distinction is the whole design:

- **baseline** (`test_baseline_matches`, green today) — the measured numbers equal
  `baseline.json`. An unrelated change that moves parity surfaces as a diff to explain
  rather than as silence.
- **parity** (`test_row_parity`, `test_column_parity`) — the destination. Every case is
  there, and both are plain assertions: no case carries an expected-failure marker and
  there is no longer any machinery to give one.

  **The ratchet is gone, and its removal is the last step of the work it governed.** While
  divergences were open, a diverging case carried `xfail(strict=True)`, so the suite went
  red the moment the case started passing and whoever fixed it re-recorded `baseline.json`
  and dropped the marker in the same commit — progress could be neither silently lost nor
  silently claimed. #35 took `test_column_parity` off it first: column parity became
  *structural* (both sides run the same vendored `compute_keep`), so a case could not drift
  back without someone breaking that, and leaving the marker available would let a genuine
  regression be absorbed by re-recording the baseline. #41 made the same call for the rows,
  because with nothing left to be marked the marker can only ever be handed to a
  regression. A case that stops passing now fails.

There is a **third state**, and it has exactly one occupant. A MAF missing a filter-input column is
refused by the pipeline (`KeyError`) and *filled and reported* by the app, by #39's deliberate
choice — so there is no pipeline verdict for the app's to agree or disagree with. `ParityResult
.off_parity_by_construction` names it, `in_parity` stays false, and `test_absent_columns.py` asserts
it against its own expectations: that the pipeline really does refuse each of the seven columns per
arm, that the app really does report, and that a complete MAF is never in this state. It is kept out
of `errors_in_parity` on purpose — folding it in would make *any* app verdict against a pipeline
`KeyError` read as parity.

**The ratchet was proven, not assumed** — recorded here because it is the evidence the
mechanism worked, and the reason it could be retired rather than merely abandoned. Adding
the missing ClinVar clause to the app's somatic pathogenic rescue (a real one-line fix for
divergence #7) turned the suite red: four cases XPASSed and six baseline records moved,
while cases the fix does not reach stayed xfailed. The patch was reverted, and every
divergence it anticipated has since been closed for real.

Plus: neutrality of the app-only frequency extra (#16, #37) asserted **algebraically** with the `if
max_freq_population < 1.0` guard bypassed, not merely guard-deep, and end-to-end (declaring it at
1.0 equals omitting it); the harness's own exclusion from both installers; and that the mirrored
copy of the app's classification vocabulary has not drifted from `config/constants.py`.

The algebraic half calls `filters.variant_filters.frequency_mask` rather than restating the
expression. It used to keep its own copy of the loop *and* of the frequency-column list,
which made it a test of the copy — and the claim it carries got stronger with #37: the mask
now gates the pipeline's whole verdict rather than only the criteria path, so a mask that is
merely *usually* all-True at 1.0 would cost parity outright. What the extra does **below**
1.0 — the composition and the ClinVar-pathogenic exemption — has no pipeline counterpart, so
it lives in `tests/test_filter_app_extras.py` by the house rule.

## The parameter sweep

Defaults alone are a weak test, and one structural fact from #18 shapes the matrix: at
contract defaults **somatic PASS is `filter_patho` ∪ (ESCAT-in-keep ∧ common ∧ VAF ∧
genes)**, so the unconditional pathogenic rescue absorbs nearly every depth, VAF or gene
edit. Every sweep case therefore also sets `skip_pathogenic`.

That is not a theoretical concern — it was measured here. The first gene cases showed
`somatic_genes_mixed_case` reaching parity *while the app was still comparing gene
symbols case-sensitively*, because the rescue re-admitted exactly the rows the gene list
excluded. The `*_skip_patho` variants exist because of it, and under them divergence #9
costs 4 of 4 rows.

Swept: `min_depth` (0, 50, 500), somatic VAF (0.01, 0.05, 0.2), germline VAF (0.05, 0.2),
`skip_pathogenic`, `skip_civic`, all three gene lists, and empty guideline keep-lists — the
all-empty state #13 made reachable, which is pipeline-expressible (`--filter_cancervar "" …`) and
bottoms out at `filter_patho`'s unconditional rescue rather than at zero.

## The baseline

> **HISTORICAL — this table is the state *before* issue #33.** It is kept because it is
> the measurement the whole effort was argued from, and because the attribution section
> below only makes sense against a diff that existed. The live numbers are in
> `baseline.json`, which issue #33 re-recorded and issue #34 added four gene-paste cases
> to: **33 of 33 cases are in row parity** and every attribution block is empty. Issue #35
> then re-recorded it again: **33 of 33 are in column parity** too, and the suite carries
> **no `xfail`s at all**.
>
> On the full reference the same fix took `germline_defaults` from 1,461 app rows against
> 818 pipeline rows — parity on **0 of 50** samples — to **818 = 818 on 50 of 50**, and
> `somatic_defaults` to 411 = 411, with zero diverging rows across all 1,100 measurements.

Measured on the app as it stands, before any fix. `cols` is
pipeline / `config/columns.py` / `_DEFAULT_VISIBLE_COLUMNS` — the three separately
drifted lists that existed then. Issue #35 collapsed them into one resolver, so the
live report prints a single `pipeline+extras` figure instead.

| case | arm | rows | pipeline PASS | app PASS | pipeline-only | app-only | cols |
|---|---|---:|---:|---:|---:|---:|---|
| `somatic_defaults_literal` | somatic | 82 | 20 | 15 | 8 | 3 | 40/40/17 |
| `somatic_defaults` | somatic | 82 | 20 | 16 | 4 | 0 | 40/40/17 |
| `germline_defaults_literal` | germline | 94 | 31 | 35 | 1 | 5 | 39/39/18 |
| `germline_defaults` | germline | 94 | 31 | 37 | 1 | 7 | 39/39/18 |
| `somatic_synthetic` | somatic | 19 | 7 | 13 | 3 | 9 | 40/40/17 |
| `germline_synthetic` | germline | 16 | 7 | 12 | 2 | 7 | 39/39/18 |
| `civic_present` | somatic | 12 | 8 | 8 | **0** | **0** | 46/46/17 |
| `civic_skipped` | somatic | 12 | 0 | 8 | 0 | 8 | 40/40/17 |
| `gnomad_genome_present` | somatic | 6 | 0 | 0 | **0** | **0** | 41/42/17 |
| `somatic_skip_patho` | somatic | 82 | 19 | 16 | 4 | 1 | 40/40/17 |
| `germline_skip_patho` | germline | 94 | 27 | 32 | 3 | 8 | 39/39/18 |
| `somatic_genes` | somatic | 82 | 18 | 13 | 5 | 0 | 40/40/17 |
| `somatic_genes_mixed_case` | somatic | 82 | 18 | 12 | 6 | 0 | 40/40/17 |
| `germline_genes` | germline | 94 | 27 | 32 | 1 | 6 | 39/39/18 |
| `somatic_genes_skip_patho` | somatic | 82 | 4 | 3 | 1 | 0 | 40/40/17 |
| `somatic_genes_mixed_case_skip_patho` | somatic | 82 | 4 | **0** | 4 | 0 | 40/40/17 |
| `germline_genes_skip_patho` | germline | 94 | 3 | 4 | 0 | 1 | 39/39/18 |
| `somatic_depth_0` | somatic | 82 | 20 | 16 | 4 | 0 | 40/40/17 |
| `somatic_synthetic_depth_0` | somatic | 19 | 9 | 12 | 2 | 5 | 40/40/17 |
| `somatic_depth_500` | somatic | 82 | 10 | 10 | 1 | 1 | 40/40/17 |
| `germline_depth_500` | germline | 94 | 17 | 18 | 2 | 3 | 39/39/18 |
| `somatic_vaf_005` | somatic | 82 | 18 | 15 | 4 | 1 | 40/40/17 |
| `somatic_vaf_02` | somatic | 82 | 14 | 13 | 2 | 1 | 40/40/17 |
| `germline_vaf_02` | germline | 94 | 27 | 32 | 3 | 8 | 39/39/18 |
| `germline_vaf_005` | germline | 94 | 28 | 32 | 3 | 7 | 39/39/18 |
| `somatic_empty_guidelines` | somatic | 82 | 18 | 12 | 6 | 0 | 40/40/17 |
| `germline_empty_guidelines` | germline | 94 | 26 | 37 | 1 | 12 | 39/39/18 |
| `somatic_empty_guidelines_skip_patho` | somatic | 82 | 0 | 0 | **0** | **0** | 40/40/17 |
| `dot_numeric` | somatic | 5 | *raises `TypeError`* | *refuses* | — | — | — |

**3 of 29 cases in row parity. 0 of 29 in column parity.** Two of those three are
`0 == 0` and prove nothing about the filter logic; only `civic_present` (8 == 8 on the
same eight variants) is a substantive agreement.

### Attribution

Every diverging row across all 29 cases is explained by a named divergence — there is no
`unattributed` bucket. Totals count rows, and a row is counted once per case it diverges
in:

| divergence | rows |
|---|---:|
| #4 `clnsig_split_vs_exact` | 47 |
| #6 `germline_escat_unmirrored` | 39 |
| #7 `somatic_patho_no_clinvar_clause` | 29 |
| #1 `vc_exclude_vs_include` | 27 |
| #8 `germline_patho_renovo_vs_clinvar` | 26 |
| #4 `clnsig_patho_split_vs_exact` | 24 |
| #2 `depth_dp_vs_sum` | 16 |
| #1 `vc_outside_app_vocabulary` | 11 |
| #2 `depth_nan_kept_by_app` | 9 |
| `skip_civic_unmirrored` | 8 |
| #3 `vaf_nan_kept_by_app` | 6 |
| #9 `gene_case_sensitivity` | 6 |
| #3 `vaf_at_threshold` | 5 |
| #7 `escat_dead_retain_clause` | 3 |

**Probes are candidate explanations, not a partition.** Several fire on the same row —
the reference subsample was built so they overlap — so the counts do not sum to the diff
size, and a high count is not a ranking of what to fix first. Measured: the four
`somatic_defaults` rows are matched by both `#4 clnsig_split_vs_exact` and
`#7 somatic_patho_no_clinvar_clause`, and fixing **only** #7 brought the case to full
parity. Definitive attribution is "apply the fix and re-measure", which is exactly what
the ratchet supports.

**An empty `unattributed` bucket meant coverage, not completeness — until #27.** The sentence above
was true of the fixtures from the first baseline and it was misleading for the whole of that time.
Run over the 50-sample reference, #25 found **3,687** diverging rows no probe explained — the single
largest bucket — and every one of them was divergence **#6**, the ESCAT clause the app ORs into the
*germline* guideline block. There was no #6 probe, so the harness could not name it; and there was
no fixture row that diverged on it, so nothing went red. The bucket read empty because no row
reached it.

The fixture set was not simply missing the shape. It held a PIK3CA `ESCAT: IA` call, Benign on
InterVar, ClinVar *and* RENOVO — a textbook #6 witness — selected by a cell named `escat_populated`.
Its `DP` was **46** against a `min_depth` of **50**, so the app dropped it on divergence #2 before
the ESCAT clause could admit it and both sides answered NOPASS. Two divergences pointing opposite
ways cancelled, and the cell reported coverage it did not have. `germline_escat_only` replaces it as
the #6 cell and pins the row's whole path to the verdict: rejected by all three pipeline guideline
sources, rescued by neither side, clearing depth on **both** columns, VAF and classification, and
inert under every other probe. The reference holds **536** such rows; #25's ablation figure of 540
is the same population measured on the app's columns alone, the 4-row difference being rows the
pipeline drops on `t_alt+t_ref` — the same masking, in the other direction.

Two tests held the two halves of the property apart, in `test_attribution_coverage.py`:

- **completeness** — no case's attribution carries an `unattributed` key;
- **soundness** — every probe in `PROBES` is named by some case's attribution, i.e. it
  explains a row that *actually diverges*.

Soundness is the one that would have caught this, and only in that form: a check that
the predicate merely *fires* somewhere would have passed on the masked PIK3CA row.
Verified as a ratchet rather than assumed — against the pre-#27 baseline the #6 probe is
named by no case at all, so the test failed; against the post-#27 one it was named by 9.

**Issue #33 then took soundness's subject away, and issue #242 removed the module.** With
every case in row parity the recorded diff is empty, so no probe can have a witness, so
the assertion was relaxed to bind only "while divergence remains" — and in its absence
asserted instead that *no* case carries an attribution block. That is the negation of what
it was written to detect, and it passes on any fixture set that reaches row parity,
including one with no coverage at all. Completeness did not invert and survives in
[`test_baseline_integrity.py`](test_baseline_integrity.py), which still reads
`baseline.json` with no `bin/` and no reference, alongside the assertion that no case
diverges — each one paired with a test that makes its predicate fire, since a dormant
check is worth exactly what its last failure proved.

Coverage itself is now measured by **mutation** rather than by attribution:
[`mutation_coverage.py`](mutation_coverage.py) re-injects each divergence into the app
side and records which cases notice. That needs no divergence to already exist, which is
what makes it survive its own success — see
[`test_mutation_coverage.py`](test_mutation_coverage.py).

### Vocabulary: criteria path vs union PASS

Two different true numbers were both being called "the somatic baseline", which is how the
411-vs-408 gap stayed open across three tickets. They are not competing measurements of one thing,
they are two cells of #20's decomposition:

- **criteria path** — `common & specific`, the guideline route into the report:
  **408** somatic, **661** germline.
- **union PASS** — `(common & specific) | filter_patho`, what the pipeline actually
  writes: **411** somatic, **818** germline.

The difference is the unconditional pathogenic rescue: 3 somatic rows, 157 germline.
`pipeline_pass` and `app_pass` in `baseline.json` and this README are always **union
PASS** — never the criteria path — because that is the verdict the parity claim is
about. Quote the cell, not the number: "the somatic baseline is 411" and "…is 408" are
both true and neither is self-describing.

## What the harness corrected

Three of the map's premises did not survive measurement.

**Divergence #12 does not exist.** The map records the app as matching CIViC with an
exact `.isin(["A", "B"])` that "cannot match list-repr at all". It does not: both the
guideline clause (`variant_filters.py:238`) and the pathogenic clause (`:295`) call the
app's own `has_element_from_list`, a substring test behaviourally identical to the
pipeline's — including on the list-repr values the annotator actually emits. There is no
CIViC `.isin` anywhere in the app. `civic_present` reaches parity end-to-end on the
fixture built specifically to separate A/B/C.

Since #33 the claim is structural rather than measured, which is stronger: the app calls
the pipeline's `has_element_from_list`, so there is no longer a pair of implementations to
compare and `test_civic_matching_already_agrees` is gone with the app's copy. The other
real CIViC divergence — the app having no `skip_civic` at all (`civic_skipped`: pipeline 0,
app 8) — closed in the same change; the flag is a contract parameter and is passed down.

**Divergence #7 is bigger than the dead ESCAT clause.** The map's table names the app's
somatic rescue as `CancerVar | ESCAT | CIViC A/B` against the pipeline's
`CancerVar | ClinVar | CIViC A/B`, but the effect of the *missing ClinVar clause* was not
separated from the *dead ESCAT clause*. It dominates: 29 diverging rows against the dead
clause's 3, and it is the entire `somatic_defaults` gap. A ClinVar-pathogenic call whose
CancerVar is Tier III/IV — `Tier_III_Uncertain` + `Pathogenic/Likely_pathogenic` is
common in the reference — is kept by the pipeline and dropped by the app.

**`min_depth = 0` is not the unfiltered case, and the two sides differ there.** The app's
`if min_depth > 0` guard skips the depth filter entirely; the pipeline still evaluates
`t_alt + t_ref >= 0`, which is **False for NaN**, so it drops rows with a missing count
even at zero. The reference fixtures cannot witness this — #18 measured zero blanks in
the numeric columns across all 181,566 reference rows — so `somatic_synthetic_depth_0`
exists to witness it: the pipeline drops 3 of 19 rows that the app keeps.

## Caveats

- **Fixtures, not the reference.** Every number in the tables above is from the ten
  committed fixtures. The full-reference numbers live in `reference.py`, run by `make
  parity-reference` (below). #16's **411 vs 408** somatic-baseline gap **is reconciled**: it was
  definitional, not an error — see *Vocabulary* above.
- **The fixture rates are not rates.** The subsample is purposive, so a divergence's
  share of the fixture diff estimates nothing about real data, in either direction.
  Measured: #6 is **51%** of the reference's attributed diverging rows and was **0%** of
  the fixtures'; #4 is the fixtures' largest at a combined 33.6% and **3.9%** on the
  reference. The fixtures are still the only place #3's NaN cells, #7's dead ESCAT
  clause and `skip_civic` are witnessed at all — the reference contains none of them.
  Use the fixtures to ask *whether* a divergence is live, the reference to ask *how much*
  it costs.
- **`0 == 0` is not parity.** Two of the three passing cases pass vacuously. Treat the
  case list, not the headline count, as the result.
- ~~**Column parity is measured against `config/columns.py`**~~ **Closed by
  #35.** It used to be measured against the app's stated `KEEP` analogue, which never filtered any
  data, while `_DEFAULT_VISIBLE_COLUMNS` — what a user actually saw — was reported alongside (17/18
  columns) but not asserted. The counts *coincided* (40 vs 40, 39 vs 39) while the contents differed
  by −`variantalker_naive` +`gnomAD_exome_AF`: **compare the lists, never the lengths**, which is
  why `test_column_parity` asserts position by position. Both lists are gone. The app now calls one
  resolver, `config.columns.resolve_visible_columns`, which builds the pipeline's half from the
  vendored `compute_keep` and appends the app's own extras — so the assertion is the prefix
  `cols[:len(keep)] == keep` and the app is a deliberate superset. The harness *imports and calls*
  that resolver; the `ast` reader it used to need, to get the grid's list without importing
  streamlit, is deleted.
- ~~**The app does not fail where the pipeline fails.**~~ **Closed by #33, improved by
  #38.** On `somatic_dot_numeric.maf` the app used to return a verdict (0 PASS of 5) where
  the pipeline raised `TypeError`, because its own `pd.to_numeric(..., errors="coerce")`
  absorbed the `.` — a silent wrong answer, the worse direction of the two. Routing the
  decision through the vendored functions made the app raise the *same* `TypeError`.
  Issue #38 then replaced that stack trace with a refusal that names every offending column
  and value, so `baseline.json` now records `app_error: "UnreadableNumericColumns"` and
  `app_refused_columns: ["t_alt_count", "t_ref_count", "tumor_f"]`, still with
  `app_pass: null`.

  **The two exceptions differ on purpose, and that is still parity.** `in_parity` reads
  through `errors_in_parity`, which accepts a pipeline raise against an app *refusal* —
  because a refusal is the pipeline's own non-verdict with a usable message, not a second
  opinion. What it does not accept is any other app-side exception against a pipeline raise:
  "both sides failed" is not parity on its own. Which columns the refusal covers is derived
  from the vendored source by `ast` (`filters/numeric_columns.py`), so the app cannot become
  fussier or laxer than the pipeline by someone editing a list — `tests/test_numeric_columns.py`
  asserts the biconditional over every fixture and over injected values on both arms.
- **`compute_keep` is not vendored yet** (flagged by
  #19), so the harness compares against the pipeline's *emitted* header rather than against a
  vendored function.
- **One build, one technology.** hg19 + iontorrent, inherited from the fixture set. The
  `gnomad_genome_present` case takes `compute_keep`'s `gnomAD_genome_AF` branch by column
  presence only — it is not a genome-build claim.

## Running the reference sweep

`reference.py` runs this same contract and these same probes over the full 50-sample
GERSOM reference — 181,566 rows against the fixtures' 170 — and it is the only
instrument that can tell you whether the fixture set still covers what real data does.

```sh
export PARITY_REFERENCE_ROOT=/path/to/the/reference/root
make parity-reference          # ~3 min at -j8; needs the GERSOM shared drive
```

**`PARITY_REFERENCE_ROOT` has no default, and that is deliberate.** The mount path names the
institute, the unit, the project and the manuscript, so issue #247 took it out of the tree. The
alternative was to strip `reference.py` from the public export, which would have left a public
reader every rate quoted below and no way to see how it was produced. The path is a fact about
one machine; the module is the evidence. `PARITY_RESOURCE_ROOT` defaults to the reference root's
sibling `resources`, and `--reference-root` overrides the variable for a one-off run.

**It is a `make` target and deliberately not a pytest test.** The reference is a shared clinical
drive, present on a handful of machines and in no CI job, so a skip-if-absent gate would be green
everywhere by skipping — the exact failure #24 found in `test_parity.py`, whose module-level
`skipif(not PIPELINE_AVAILABLE)` removes the whole parity suite from any checkout without `bin/`
without saying so. A target you have to type cannot be mistaken for a check that ran. Nor is it a
docs-only script any more, which is how the blind spot survived: it was written, measured once, and
never run again while the fixture baseline went on reading as complete.

`--strict` (which the target passes) gates on two invariants, both of which expect the
constant **0**, so neither has a number that can go stale:

1. **No diverging row is unattributed.** Precisely the property that was false by 3,687
   rows while the fixture baseline read as complete. Proven to fail for the right
   reason: with the #6 probe removed, a two-sample sweep exits 1 on 19 unexplained rows;
   with it present, 0.
2. **The four-cell decomposition reproduces `bin/filter_variants.py` exactly** — 0
   mismatching sample-cases out of 1,100. This is a strictly stronger signal than any PASS count,
   because the cells are compared to a fresh pipeline subprocess as a *set of variant keys* per
   sample, so two runs cannot agree by cancelling errors the way equal totals can. It was already
   being computed and merely reported; asserting it costs nothing. If it fires, the vendored masks
   no longer reconstruct the pipeline's verdict and every diagnostic built on them — #20's cells,
   #21's `MAFigate_reason` — is describing a decision the pipeline did not make.

Run it whenever `PROBES` or the fixture cells change. **Leave `--strict` off while
investigating** — an unattributed bucket is the finding, not a failure, and
`CANDIDATE_PROBES` in `reference.py` is where a suspected explanation gets measured
*before* it is allowed to edit `harness.PROBES`. That separation is what made the blind
spot findable: tuning the instrument to the bucket it found would have proved nothing.
The dict is empty today and kept for the next one.

Results are written with `--out`; the committed copy in `docs/wayfinder/issue-27/` is
scrubbed on #18's precedent (samples aliased, drive path dropped, germline coordinates
reduced to counts). `--include-row-keys` writes raw variant coordinates — real germline
calls, keep them local.

## Privacy

`tests/parity/` is excluded from both installers, asserted by
`test_harness_is_excluded_from_packaged_builds`. The fixtures carry real germline variant
calls from 50 patients; see the fixture set's own privacy note before making this
repository public.
