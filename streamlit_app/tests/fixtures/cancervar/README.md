# CancerVar evidence fixtures

Two small MAFs carrying CancerVar's `CancerVar and Evidence` column, added by issue #189. Before
them no fixture in `streamlit_app/tests/` carried a CancerVar evidence string at all, so every guard
about the somatic CBP section would have had to hand-build the column name — and the column name is
one of the things that can go wrong.

## What is real and what is not

The **evidence strings are real**, byte for byte, lifted from somatic MAFs on the dev machine via
`docs/wayfinder/issue-188/cbp_shapes.json`. Everything else is synthetic: the gene symbols are
real gene names but the positions are round numbers, the read counts are invented, and no sample
barcode, patient identifier or original filename appears anywhere. An evidence string is twelve
small integers and a tier — it carries no patient data — which is why it can be checked in while
the files it came from cannot.

## The two spellings

The evidence column has two spellings in the real corpus, and both are here because
`components/variant_detail._evidence_column` matches by substring and so survives both *by luck
rather than by design* — which makes it exactly the kind of thing a guard should pin.

| file | column spelled | `CancerVar` column |
| --- | --- | --- |
| `somatic_cancervar_evidence.maf` | `' CancerVar: CancerVar and Evidence '` — padded, 121 of 124 real files | present |
| `somatic_cancervar_dotted_no_column.maf` | `'CancerVar..CancerVar.and.Evidence'` — R-mangled, 3 real files | **absent** |

The padding is real and is load-bearing: the leading and trailing spaces appear in both the
column *name* and every *value*. It sits on an interior column in both files so that no
whitespace-trimming tool can quietly repair it.

The second file also has **no `CancerVar` column**, which is the state **15 of the 109** real files
carrying an evidence vector are in (16 of 124 when #189 counted; issue #210 re-read the corpus with
the repo's own fixtures excluded from it — including this directory, whose `test_sample.maf` sibling
had been counted as a real file). On those files the tier exists only inside the evidence string,
which is what `render_cbp_evidence` badges and what issue #189 reworded `_classification_absences`
to point at rather than deny.

Issue #210 grew that file by two rows, because the state it was built for is not the only one a
file with no `CancerVar` column reaches — and on the other one the caption above the section was
**false**. See the row list below.

## The rows

`somatic_cancervar_evidence.maf`, in order:

1. **PRDM16** — the single commonest vector in the corpus (15,313 of 114,053 rows in #185): in
   COSMIC and ICGC, in a cancer gene, and ClinVar calls it benign. The row that proves a
   positive-only table would be wrong.
2. **SNRNP200** — all twelve zero, the section's empty state (7.4% of real rows).
3. **AGRN** — both criteria that *can* go negative are negative, so the printed sum carries a
   minus sign: `-2#`. The parser has to read that sign.
4. **MPL** — three criteria at 1, just over the Tier IV boundary at a sum of 3.
5. **PIK3CA p.M1043V** — Tier II, with the predictors dissenting (`CBP10 = -1`).
6. **PIK3CA p.E545G** — a genuine Tier I (48 rows of 109,416), nine of twelve criteria fired.
7. **KRAS** — a tier and a sum with **no `EVS=[...]` vector**. `parse_cancervar` returns `None`
   for it rather than badging a tier with no evidence behind it, and here the guideline row still
   shows the tier from the `CancerVar` column — because *this file has one*.

   Issue #210 found that the 73 real cells in this state **do not**: they are all in one file
   with no `CancerVar` column, so on the real rows the tier is drawn nowhere at all, which is a
   different state from the one this row models. The real one is row 3 of the dotted file below.
   The state this row models — the column present *and* the cell unparseable — has **0 rows in
   the corpus**, which is worth keeping and worth knowing.
8. **BRAF** — both columns present, both cells `.`. The variant-level absence, and rarer than
   #187 recorded: **1 row of the 107,296 carrying `CancerVar`, and 0 of the 211,634 carrying
   `InterVar`**.

`somatic_cancervar_dotted_no_column.maf`, in order — no `CancerVar` column on any of them:

1. **TP53** — a full twelve-entry vector. The tier is drawn from the string, and the caption above
   points at it.
2. **PIK3CA** — a genuine Tier I from the string alone.
3. **KRAS** — a tier and a sum with **no vector**, on a file with no `CancerVar` column: the
   73-row state exactly, and the one the caption used to lie about. Added by issue #210.
4. **EGFR** — evidence cell `.`, on a file with no `CancerVar` column. **0 real rows** are in this
   state, but the branch is reachable and its sentence has to be true, so it is pinned. Added by
   issue #210.
5. **BRAF** — a vector summing to **zero** with a real tier printed beside it
   (`0#Tier_IV_benign`), and no `CancerVar` column: **2,115 real rows in 13 files**. Added by issue
   #210 because a mutation survived without it — reading the verdict from the parsed *sum* instead
   of the *tier* passes every other guard here, and `0 or ""` is `""`, so a drawn Tier IV would be
   captioned as no tier at all.

Rows 3 and 4 are the two rows where nothing is drawn below the caption, which is why the caption
there says no tier is shown rather than pointing one out. Row 5 is the row that looks like them to
anything counting instead of reading.
