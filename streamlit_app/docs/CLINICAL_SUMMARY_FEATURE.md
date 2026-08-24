# Clinical Summary Feature

## Overview
The Clinical Summary feature adds a new first column to the results table that provides a consolidated view of the clinical significance of each variant based on multiple clinical databases.

## How It Works

### Priority Hierarchy (Highest to Lowest)

**Two axes, not one.** The first seven are the ACMG/AMP pathogenicity ladder — how likely the
variant is to cause disease. The rest are the other things a clinical source can assert, and
are *not* points on that ladder; ClinVar itself keeps them in separate groups.

1. **🔴 Pathogenic** - Clearly disease-causing
2. **🔴 Pathogenic (low penetrance)** - Disease-associated, but many carriers are unaffected
3. **🟠 Likely Pathogenic** - Probably disease-causing
4. **🟠 Likely Pathogenic (low penetrance)** - As above, with reduced penetrance
5. **🟡 Uncertain Significance** - Unknown clinical impact
6. **🟢 Likely Benign** - Probably not disease-causing
7. **🔵 Benign** - Not disease-causing

Below the ladder, ordered by the term table's own "clinical relevance for filtering":

8. **⚠️ Disease Risk** - A susceptibility or risk allele, not a disease-causing call
9. **💊 Drug Response** - Affects response to a drug rather than disease risk
10. **🔗 Association or Trait** - A GWAS association, or a non-disease phenotype
11. **🛡️ Protective** - *Decreases* risk of a disorder
12. **⚪ No Classification** - A record exists and carries no usable classification
13. **⚪ Unrecognised Annotation** - A value present in the file that this app cannot class

And, outside the hierarchy because it is not a classification:

- **🔍 No Clinical Data** - **no source annotated this variant at all.** Said of nothing else.
  Until issue 98 it was also said of a variant whose annotation the app could not read,
  which is a different thing and false of it — that case is now item 13.

A low-penetrance class ranks directly below its full-penetrance counterpart: the term table
calls both "keep, but flag separately" and warns they must not be read as equivalent.

**The glyph is the class, in both clinical columns** (issue 104). `Pathogenicity_Overview`
draws these same glyphs, one per annotation source, so 🔵 means Benign whether it is read in
a summary label or at a position in the circle strip. Three pairs share a glyph, and the
words are what tell them apart: each low-penetrance class draws its full-penetrance
counterpart's, and items 12 and 13 both draw ⚪ — a source that recorded no usable
classification and a value this app could not read are not distinguished by the circles,
which is why the key beside the grid calls that circle *No usable classification* rather than
claiming which of the two it was.

### Data Sources Analyzed

The clinical summary reads five annotation columns. **No source outranks another**: the label
is the strongest class *any* of them asserts, which is what the Logic section below describes.
This list used to be introduced as "in order of priority", and there was never a priority —
`generate_clinical_summary` carried a 1-to-5 number per column that nothing read, deleted at
issue 108. The order below is the pipeline's, and it is presentation only.

One consequence worth stating: because the strongest class wins, reading one more source can
only move a label *up* the ladder or leave it, never down.

1. **CancerVar** - Cancer-specific variant classification
   - `Tier_I_strong` → Pathogenic
   - `Tier_II_potential` → Likely Pathogenic
   - `Tier_III_Uncertain` → Uncertain Significance
   - `Tier_IV_benign` → Benign

2. **InterVar** - ACMG guidelines-based classification
   - `Pathogenic` → Pathogenic
   - `Likely pathogenic` → Likely Pathogenic
   - `Uncertain significance` → Uncertain Significance
   - `Likely benign` → Likely Benign
   - `Benign` → Benign

3. **ClinVar** - Clinical variant database

   The values that reach this field are far more than the ACMG/AMP five — some issued by
   ClinVar, some submitter-defined — and they are **grouped by term class** rather than
   listed one by one here. The per-term list lives in
   `CLINICAL_VALUE_MAPPING` (`components/clinical_summary.py`), which is the only copy. This
   section used to transcribe it and had fallen behind by a dozen terms.

   - Pathogenicity terms → the ladder above, with the VUS sub-tiers (`VUS-high`, `VUS-mid`,
     `VUS-low`) all reading as Uncertain Significance, since the sub-tiers are
     submitter-defined and ACMG/AMP gives no criteria for them
   - `risk_factor` and the `*_risk_allele` terms → Disease Risk
   - `drug_response`, `confers_sensitivity` → Drug Response
   - `association`, `Affects` → Association or Trait
   - `protective` → Protective
   - `other`, `not_provided`, `conflicting_data_from_submitters` and the "no classification"
     aggregates → No Classification

   The classes come from the institute's ClinVar term table, cited in that module's
   docstring, which agrees with ClinVar's own grouping of its non-standard terms. Nothing
   here assigns a severity that a source did not assign: `protective` is not Benign and
   `risk_factor` is not Pathogenic, and tests refuse either.

4. **RENOVO** - a machine-learning prediction of pathogenicity, with a confidence

   RENOVO's six classes are **six contiguous bins of one score**, and the MAF carries that
   score: the pipeline copies RENOVO's `PL_score` in as `RENOVO_pls` beside the class. So
   `HP`/`IP`/`LP` is a confidence on one prediction rather than a second scale, and the two
   `LP` classes are the low-confidence pair either side of the 0.5 decision boundary —
   `LP Pathogenic` bottoms out at 0.50000, `LP Benign` tops out at 0.49985 — which is why both
   read as Uncertain Significance. The grading is the pipeline's own, from
   `bin/generate_clinical_summary.py`.

   The six are **not transcribed here**, for the reason the ClinVar section above gives:
   `CLINICAL_VALUE_MAPPING` is the only copy, and the copy that used to sit here had drifted
   from it. `tests/test_clinical_summary.py` holds the six against RENOVO's score order.

   Until issue 108 this section described a mapping that decided nothing on any real file:
   the function read `filter_renovo`, a filter *parameter* name, so RENOVO reached this column
   on no row ever rendered.

5. **CIViC** - Evidence-based clinical interpretation

   CIViC's levels grade the *evidence* behind an assertion, not the variant, so they are read
   from their own scale in `_CIVIC_LEVEL_CLASSES` rather than through `CLINICAL_VALUE_MAPPING`.
   **No level maps to a benign class**: `A` → Pathogenic, `B` → Likely Pathogenic, and `C`, `D`
   and `E` all → Uncertain Significance, because `E` means the only evidence is indirect, which
   is an absence of strong evidence rather than a benign call (issue 109).

   This section said `E` → Benign until then, and kept saying it for one commit after the code
   changed — the copy that is now a pointer rather than a transcription, as above.

### Logic
- If **any** database classifies a variant as "Pathogenic", the summary will be "🔴 Pathogenic"
- If **no** database says "Pathogenic" but **any** says "Likely Pathogenic", the summary will be "🟠 Likely Pathogenic"
- The algorithm follows the hierarchy to assign the highest significance level found across all databases
- Visual emojis make it easy to quickly identify the clinical significance at a glance

### The pipeline answers this question too

The annotation workflow runs `bin/generate_clinical_summary.py` as a process of its own and
writes its verdict into the MAF as `variantalker_naive`. That column survives filtering, so it
is in every filtered frame beside `Clinical_Summary`. The two are the same question asked twice,
and from issue 108 they agree: `tests/test_clinical_summary.py` compares them row for row on
the committed germline reference.

**The report shows one of them** (issue 117). Measured over 176 annotated real MAFs, the two
hold identical words on every somatic row and on 51 of the 52 germline files that carry both;
the one file where they part had its column written before the pipeline read RENOVO, so there
the app is the current reading. `config.columns.PIPELINE_COLUMNS_THE_APP_REPLACES` therefore
leaves the pipeline's copy out of the default view — the grid's *Show all columns* still shows
it, the column reference still describes it, and the row-for-row comparison still reads it.

Two differences are deliberate. The pipeline pools ClinVar's non-standard terms into two buckets
where the app gives them six classes from the institute's term table (issue 98); and CIViC's
level `E` reads as Uncertain Significance here where the pipeline reads it as Benign, because a
level meaning *the only evidence is indirect* is an absence of strong evidence rather than a
benign call (issue 109). The first is excluded from the row-for-row comparison explicitly; the
second no committed fixture can show, neither reference MAF carrying a CIViC column.

### Benefits
1. **Quick Assessment** - Immediately see the clinical significance without checking multiple columns
2. **Prioritized View** - Most significant classification takes precedence
3. **Visual Clarity** - Color-coded emojis for instant recognition
4. **Comprehensive** - Considers the five classification sources listed above. **Not ESCAT**,
   which grades how actionable a *target* is rather than how pathogenic a *variant* is; the
   `Pathogenicity_Overview` circles beside this column do draw it, and say so in their key
5. **Standardized** - Consistent classification scheme across different sources

### Usage in Results
- The Clinical Summary appears as the **first column** in both "Passed Filters" and "Failed Filters" tabs
- Download files include the Clinical Summary column
- Summary statistics show breakdown by clinical significance
- Sorting and filtering work with the clinical summary values

### Example
If a variant has:
- CancerVar: "Tier_II_potential" (→ Likely Pathogenic)
- ClinVar: "Pathogenic" (→ Pathogenic)
- InterVar: "Uncertain significance" (→ Uncertain Significance)

The Clinical Summary will be: **🔴 Pathogenic** (highest priority classification)
