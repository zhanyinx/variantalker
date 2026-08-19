"""How a variant's clinical significance is summarised.

Reads the annotation columns a MAF carries — ClinVar, InterVar, CancerVar, RENOVO,
CIViC, ESCAT — and reduces them to the two columns the table leads with:
``Clinical_Summary``, one classification chosen by priority, and
``Pathogenicity_Overview``, one coloured circle per source.

The derivation happens once, at filter time: ``page_modules/data_loading.py`` calls
:func:`circle_sources` for the arm and the file, then :func:`add_clinical_summary_column`
on the passed and failed frames, so every later reader sees the columns already there.
:func:`format_clinical_summary_display` tidies the counts for the summary tab, and
``variant_table`` writes the key beside the grid from three things: the sources this report
drew (issue #95), plus :func:`pathogenicity_circle_legend` and :func:`circle_axis_notes` —
the glyphs, and what does not follow from a glyph (issue #100).

**Both halves of that key are derived, and they were fixed by different tickets.** #100 took
the glyphs, which named six of the eight the circles can draw; #95 took the sources, which
named all six on both arms whatever the file held. Neither half may go back to being
written out by hand — that is how ``⚪`` came to be drawn beside a key omitting it, and how
a germline reader came to be told about ``CancerVar``.

The two columns do **not** read the same sources. ``Pathogenicity_Overview`` draws a circle
per source :func:`circle_sources` returns, ESCAT among them;
:func:`generate_clinical_summary` reads five and not ESCAT, because ESCAT grades a target's
actionability and that column picks one *classification*. Circles keep it, with the key
saying which axis it is on (issue #100).

**Not every source is read off the pathogenicity ladder.** ESCAT grades how actionable a
target is and CIViC grades the *evidence* behind an assertion about the variant, so each is
mapped from its own scale — :data:`_ESCAT_GROUP_CIRCLES` and :data:`_CIVIC_LEVEL_CLASSES` —
rather than through :data:`CLINICAL_VALUE_MAPPING`, which holds the classifications a source
asserts *of the variant*. A position mapped that way draws glyphs the key labels for the
ladder, so it owes the reader a clause saying what it is really graded on;
:func:`circle_axis_notes` is which positions owe one (issues #100 and #109).

**They do speak one glyph vocabulary, since issue #104.** Reading different sources — and
grading some of them on their own scale, per the paragraph above — is what makes the two
columns worth pinning side by side; spelling the same *class* differently was only ever an
accident of their having grown up in separate functions, and it had survived in three places
— ✅/🔵 for Benign, ❓/⚪ for No Classification, and ❔/⚪ for an unreadable value, which is the
pair issue #98 opened this ticket's fourth bullet with. Both columns now read
:data:`_CLASS_GLYPHS`, so a class has one glyph wherever it is drawn and the columns can
disagree about a variant without disagreeing about a word. What still differs is what each
column has *room* for: the summary writes words after its glyph and can separate two classes
that share one, and the strip cannot.

**The pipeline answers this question too, and its answer travels in the MAF.** The annotation
workflow runs ``bin/generate_clinical_summary.py`` as a process of its own
(``modules/local/annotation/small_variants/main.nf``), downstream of RENOVO and upstream of
``filter_maf``; it writes ``variantalker_naive``, which ``compute_keep`` keeps on both arms.
This module is a re-implementation of that script — same five columns, same ClinVar split, a
hierarchy where it has seven numeric tiers — which makes it *checkable* rather than merely
similar: ``tests/test_clinical_summary.py`` compares the two row for row on the committed
germline reference, and issue #108 is what made that comparison pass.

**The report no longer shows both** (issue #117). It did until then — the pipeline's verdict
was one of the forty visible columns, underscore-spelled, a few columns along from this
module's glyphed answer to the same question, with nothing on screen saying they were twins.
``config.columns.PIPELINE_COLUMNS_THE_APP_REPLACES`` takes it out of the default view and
records the measurement that decided it: identical words on 100% of somatic rows and on 51 of
52 germline MAFs, the whole disagreement in one file whose column predates RENOVO, and 74 of
176 annotated MAFs carrying no such column at all. **The comparison is untouched by that** —
it reads the fixture's column directly, not the resolver's list, so taking the column off the
screen does not take away the oracle.

**Two divergences from it are deliberate, and neither is drift.** ClinVar's non-standard terms
get six classes of their own here where the pipeline pools all of them into ``Not_Provided`` and
``Unknown`` (issue #98, from the institute's term table); and CIViC's ``E`` reads as
``Uncertain_Significance`` where the pipeline tiers it *Benign*, because a level meaning "the
only evidence is indirect" is an absence of strong evidence rather than a benign call (issue
#109). The first is what that test's row filter excludes — failing rather than skipping when a
row falls outside what it can place, so a new case cannot pass unnoticed. **The second it cannot
see**, neither committed reference MAF carrying a CIViC column, and it is written down here for
that reason. ``docs/CLINICAL_SUMMARY_FEATURE.md`` says the same in prose for a reader who is not
reading tests, and is the copy to correct if this changes.

**Which sources the overview draws is decided once and carried, never recomputed.**
``data_loading`` stashes :func:`circle_sources`'s answer and ``variant_table`` draws the
key from that, because the two inputs do not age together: switching arms on the
parameter page resets ``filter_params`` without clearing the report
(``parameter_config.py`` warns, ``data_loading`` only *detects* the staleness), so a
key recomputed at render time would name three sources over cells holding four. Same
shape as issue #92, where ``create_data_table`` returns the column list it resolved
rather than letting the download read the widgets back out.

**Where the ClinVar classes come from (issue #98).** ``CLINICAL_VALUE_MAPPING`` groups each
ClinVar term by the *term class* it belongs to, taken from the institute's own ClinVar term
table — the authority for this repo, supplied by the dev driving the map, and the thing
issue #88 said did not exist when it declined to classify ``other``. Its classes:

===================================== ==============================================
Term class                            Terms
===================================== ==============================================
Pathogenicity / disease-causing       ``Pathogenic``, ``Likely_pathogenic``,
                                      ``Uncertain_significance``, ``VUS-high/mid/low``,
                                      ``Pathogenic,_low_penetrance``,
                                      ``Likely_pathogenic,_low_penetrance``,
                                      ``Benign``, ``Likely_benign``
Disease-risk classification           ``risk_factor``, ``Uncertain_risk_allele``,
                                      ``Likely_risk_allele``, ``Established_risk_allele``
Pharmacogenomic / drug-response       ``drug_response``
Disease-association or phenotype      ``association``, ``Affects``
Protective-effect                     ``protective``
No classification / metadata          ``other``, ``not_provided``,
                                      ``conflicting_data_from_submitters``, ``-``
===================================== ==============================================

The table agrees with ClinVar's own documentation, which lists the non-standard terms
separately from the ACMG/AMP five and states they "are not considered in the determination
of conflicting classifications" — so the grouping is the source's, not this app's. What the
table adds, and ClinVar does not, is the institute's filtering call per term, which is what
orders the classes below ``Benign`` in :data:`CLINICAL_HIERARCHY`.

The table is a table of values that reach the ``ClinVar_VCF_CLNSIG`` *field*, which is not
the same as values ClinVar itself issues: the ``VUS-high``/``VUS-mid``/``VUS-low`` sub-tiers
are submitter-defined, and the table says so — ACMG/AMP permits internal sub-tiers and
defines no criteria for them, which is why all three read as Uncertain Significance here
rather than being ranked against each other.

Two placements are **not** derived from it and are flagged where they are made:
``confers_sensitivity`` (in neither the table nor ClinVar's current page; the dev's call) and
two measured aggregate spellings the table does not list. Nothing else here is a judgement.

Column *ordering* is not here. ``PRIORITY_COLUMNS`` and ``prioritize_columns`` sat in
this file's ancestor and moved to ``config/columns.py`` with issue #76: they are pure,
streamlit-free and about which column comes first, which is that module's stated job.

Not in ``config/``, though most of it is pure pandas: :func:`add_clinical_summary_column`
reports a failed derivation through ``st.warning``, and ``config/`` is importable
without Streamlit on purpose (the parity harness depends on it). Issue #54 drew that
same line for ``page_modules/param_store.py`` — rewriting an error path is not a move.
"""

import streamlit as st
import pandas as pd

from config.columns import pipeline_output_columns
from config.missing_values import says_nothing
from config.vocabularies import CIVIC_DEFINITIONS, ESCAT_DEFINITIONS



#: The classes a variant can be put in, strongest first.
#:
#: Two axes, not one. The first seven entries are the ACMG/AMP pathogenicity ladder — how
#: likely this variant is to cause disease. The five after ``Benign`` are **not points on
#: that ladder**: they are the other things a clinical source can assert, and ClinVar keeps
#: them in separate groups for exactly that reason, stating that its non-standard values
#: "are not considered in the determination of conflicting classifications".
#:
#: That split is not new here, it was just unnamed: ``Drug_Response`` and the class formerly
#: called ``Not_Provided`` have sat below ``Benign`` since before this map, and both are
#: ClinVar non-pathogenicity terms. Issue #98 read the two as a precedent rather than an
#: oddity and gave the remaining classes the same treatment, from the term table this
#: module's docstring cites.
#:
#: **Order below ``Benign`` follows the institute's own "clinical relevance for filtering"
#: column** — review-only before usually-exclude before exclude — and *not* a claim that a
#: risk allele matters less than a benign call. Ranking the whole group under the ladder is
#: the existing convention this file already followed for ``Drug_Response``; a reader who
#: wants the other reading should move the whole group, not one member.
#:
#: **The column's own values, for the three #114 asked the dev for**, so the ordering above can
#: be checked rather than trusted: ``Disease_Risk`` is *review-only*, ``Association`` and
#: ``Protective`` are both *usually-exclude*. Consistent with the order as written — and it
#: bounds ``Drug_Response``, which sits between them here and was not asked about, to somewhere
#: from review-only to usually-exclude. **``Association`` before ``Protective`` carries no
#: authority at all**: the column ties them, so that pair's relative order is this file's
#: arbitrary tie-break and nothing may be read from it.
#:
#: What this column does **not** decide is what a preset keeps — #114 settled that separately
#: and the other way, both Broad presets keeping all three. The argument is not repeated here;
#: it lives with the decision, at ``_SOFT_CLINVAR_TERMS`` in ``config/presets.py``, together
#: with its measured row cost.
#:
#: A low-penetrance term ranks directly below its full-penetrance counterpart. The table
#: calls both "keep, but flag separately" and warns they must not be read as equivalent to a
#: fully penetrant call — so they are neither folded into it nor ranked away from it, and the
#: flag lives in the rendered label.
CLINICAL_HIERARCHY = [
    "Pathogenic",
    "Pathogenic_Low_Penetrance",
    "Likely_Pathogenic",
    "Likely_Pathogenic_Low_Penetrance",
    "Uncertain_Significance",
    "Likely_Benign",
    "Benign",
    # Below here: assertions that are not pathogenicity calls.
    "Disease_Risk",
    "Drug_Response",
    "Association",
    "Protective",
    "No_Classification",
    "Unknown",
]

#: One glyph per class, and the one place a class's glyph is chosen.
#:
#: **Both derived columns read this**, which is issue #104's answer: a glyph means the same
#: class in ``Clinical_Summary`` as in ``Pathogenicity_Overview``, and a class can no longer
#: be given one glyph in one column and a different one in the other. It could before, and
#: three of the thirteen were — measured across 183 byte-distinct real annotated MAFs and
#: 330,189 rows, by running the pre-change module beside this one:
#:
#: * ``Benign`` drew ✅ here and 🔵 there on **248,034 rows**, three quarters of the corpus;
#: * ``No_Classification`` drew ❓ here and ⚪ there on **427 rows** — the same *words* under
#:   two glyphs, since the key spells that circle out in the sentence beside the grid;
#: * ``Unknown`` drew ❔ here and ⚪ there on **18 rows**, which is the disagreement issue #98
#:   created while repairing the summary and handed to #104 by name.
#:
#: **The circles' glyphs are the ones kept**, and the reason is the ramp rather than taste.
#: 🔴🟠🟡🟢🔵 is a five-step severity scale with a distinct colour at each step; ✅ is a second
#: green beside 🟢, so folding the circles onto it would have left Likely Benign and Benign
#: reading alike in a strip whose whole job is to be scanned. The words were already the same
#: in both columns and did not move.
#:
#: **The table is deliberately not injective.** ``No_Classification`` and ``Unknown`` share
#: ⚪ — the dev's call, recorded at :data:`_UNREADABLE_CIRCLE` with what it costs — and each
#: low-penetrance class shares its full-penetrance counterpart's glyph. Nothing here may
#: assume a glyph identifies a class; what a glyph identifies is a *class's rendering*, and
#: the words beside it separate the cases the strip collapses. Which is why only one of the
#: two columns can carry the penetrance flag: it is in the words, and the strip has none.
_CLASS_GLYPHS = {
    "Pathogenic": "🔴",
    "Pathogenic_Low_Penetrance": "🔴",
    "Likely_Pathogenic": "🟠",
    "Likely_Pathogenic_Low_Penetrance": "🟠",
    "Uncertain_Significance": "🟡",
    "Likely_Benign": "🟢",
    "Benign": "🔵",
    "Disease_Risk": "⚠️",
    "Drug_Response": "💊",
    "Association": "🔗",
    "Protective": "🛡️",
    "No_Classification": "⚪",
    "Unknown": "⚪",
}

#: What the app calls each class in words, and the one place those words are written.
#:
#: Replaces an ``if``/``elif`` chain that spelled every label inline, which is why the
#: sort order beside it was a *second* hand-written copy of the same thirteen strings and
#: ``charts._strip_emoji`` was a third. Everything that needs a label, its glyph or its rank
#: now reads these two dicts, so a class cannot be added to the hierarchy and left unnamed.
#:
#: These are ``Clinical_Summary``'s words. The key beside the grid writes its own, shorter
#: ones (:data:`_KEY_WORDS`) because it labels a glyph in a legend rather than a cell in a
#: column — but it can only *shorten*, never re-name, since both are read off one glyph now.
_SUMMARY_WORDS = {
    "Pathogenic": "Pathogenic",
    "Pathogenic_Low_Penetrance": "Pathogenic (low penetrance)",
    "Likely_Pathogenic": "Likely Pathogenic",
    "Likely_Pathogenic_Low_Penetrance": "Likely Pathogenic (low penetrance)",
    "Uncertain_Significance": "Uncertain Significance",
    "Likely_Benign": "Likely Benign",
    "Benign": "Benign",
    "Disease_Risk": "Disease Risk",
    "Drug_Response": "Drug Response",
    "Association": "Association or Trait",
    "Protective": "Protective",
    "No_Classification": "No Classification",
    "Unknown": "Unrecognised Annotation",
}

#: What a ``Clinical_Summary`` cell holds: the class's glyph, then the class's words.
#:
#: Derived, so the two halves cannot drift apart and neither can drift from the circles.
_SUMMARY_LABELS = {
    name: f"{_CLASS_GLYPHS[name]} {words}" for name, words in _SUMMARY_WORDS.items()
}

#: Said of a variant no source has annotated at all — and of nothing else.
#:
#: It used to be said of a variant whose annotation this module simply could not read, which
#: is a different thing and false of it. See :func:`generate_clinical_summary`.
NO_CLINICAL_DATA = "🔍 No Clinical Data"

#: Said when the derivation itself raised. Distinct from every value above.
#:
#: ❔ identifies it, and only since issue #104: it used to be this label's glyph *and*
#: ``Unrecognised Annotation``'s, so one column said ❔ for a value it could not read and ❔
#: for a frame whose derivation had thrown — a collision inside a single column, next to the
#: one across two that #104 was opened for. The unreadable case draws ⚪ with the circles now,
#: which leaves this glyph naming one thing.
ANALYSIS_ERROR = "❔ Analysis Error"

# Mapping from various database values to standardized hierarchy
CLINICAL_VALUE_MAPPING = {
    # CancerVar mapping
    "Tier_I_strong": "Pathogenic",
    "Tier_II_potential": "Likely_Pathogenic",
    "Tier_III_Uncertain": "Uncertain_Significance",
    "Tier_IV_benign": "Benign",
    # InterVar mapping
    "Pathogenic": "Pathogenic",
    "Likely pathogenic": "Likely_Pathogenic",
    "Uncertain significance": "Uncertain_Significance",
    "Likely benign": "Likely_Benign",
    "Benign": "Benign",
    # ClinVar mapping
    "Pathogenic": "Pathogenic",
    "Likely_pathogenic": "Likely_Pathogenic",
    "Pathogenic/Likely_pathogenic": "Pathogenic",
    "Uncertain_significance": "Uncertain_Significance",
    "Likely_benign": "Likely_Benign",
    "Benign": "Benign",
    "Benign/Likely_benign": "Likely_Benign",
    # -- ClinVar, beyond the ACMG/AMP five (issue #98). ------------------------------------
    # Every entry below is placed by the *term class* column of the institute's ClinVar term
    # table, which this module's docstring cites. Nothing here invents a severity: a term
    # is grouped with the terms its own source groups it with, and the class is what the app
    # then says. Before this, all of them fell through to "Unknown", were dropped, and left
    # the variant reported as carrying no clinical data at all.
    #
    # A submitter-defined VUS sub-tier is still a VUS — the table's words — and ACMG/AMP
    # defines no criteria for the sub-tiers, so the app does not rank them against each other.
    "VUS-high": "Uncertain_Significance",
    "VUS-mid": "Uncertain_Significance",
    "VUS-low": "Uncertain_Significance",
    # Clinically relevant, but the table warns these must not be read as equivalent to a
    # fully penetrant call, so they keep their own classes and say so in the label.
    "Pathogenic,_low_penetrance": "Pathogenic_Low_Penetrance",
    "Likely_pathogenic,_low_penetrance": "Likely_Pathogenic_Low_Penetrance",
    # Disease-risk classification. `risk_factor` is the older OMIM-derived spelling of the
    # same idea as the three `*_risk_allele` terms, which is why all four share a class.
    "risk_factor": "Disease_Risk",
    # The underscored spellings of the four terms issue #103 put on offer, added for the same
    # contract the `_drug_response` note below states: every offered term is nameable here.
    # #99's rule, applied to #103's additions — the split leaves the prefix on and files differ
    # over whether the modifier keeps it. The row counts witnessing each are in
    # `config/vocabularies.py` beside the vocabulary that offers them, and deliberately not
    # repeated here: they are a fact about which terms to *offer*, and a second copy of a
    # measurement is a second thing to re-measure.
    #
    # Like `_drug_response` and `_other`, these cover the *bare-cell* case only:
    # `generate_clinical_summary` reads one piece of a composite, so `Benign|_risk_factor` is
    # decided by `Benign` whether these entries exist or not. `Affects|_association` is the
    # one measured cell where an underscored piece is load-bearing for the *filter* — 2 rows,
    # reachable by no term the control offered before #103.
    "_risk_factor": "Disease_Risk",
    "Uncertain_risk_allele": "Disease_Risk",
    "Likely_risk_allele": "Disease_Risk",
    "Established_risk_allele": "Disease_Risk",
    # Pharmacogenomic. `confers_sensitivity` is the dev's call and the one term here placed
    # without a citable definition: it is absent from the institute's table and from
    # ClinVar's current classification page, being a legacy value. Recorded as a judgement
    # rather than a derivation, so a later reader can overturn it without re-deriving it.
    "drug_response": "Drug_Response",
    # The same call as the line above, spelled as the split leaves it when the modifier
    # followed a call in the annotation (issue #99). Both spellings are offered and both Broad
    # presets keep both, so the vocabulary's contract -- every offered term is nameable by this
    # mapping -- reaches this spelling too.
    #
    # What it does **not** do is decide a real composite value, and the first version of this
    # comment claimed it did. `generate_clinical_summary` picks exactly *one* piece of a
    # composite (the first containing "pathogenic", else the first piece), so
    # `Uncertain_significance|_drug_response` is read as `Uncertain_significance` whether this
    # entry exists or not, and no cell in the measured corpus is entirely an underscored piece.
    # So this covers the bare-cell case alone, which is unobserved but offerable -- unlike
    # #88's conflicting-interpretations entry, which fixed 77 *real* rows because that spelling
    # does stand as a whole cell. Kept for the contract, not for a row it rescues.
    "_drug_response": "Drug_Response",
    "confers_sensitivity": "Drug_Response",
    "_confers_sensitivity": "Drug_Response",
    # Disease-association or phenotype annotation. `Affects` is a non-disease phenotype
    # (the table's example is lactose intolerance) rather than an association, which is why
    # the rendered label reads "Association or Trait" and not "Association".
    "association": "Association",
    "_association": "Association",
    "Affects": "Association",
    # Protective-effect annotation — the opposite end of the axis from `risk_factor`, and
    # the sharpest reason these five could not share one bucket.
    "protective": "Protective",
    "_protective": "Protective",
    # No classification / metadata. `other` lands here, which is the question issue #98 was
    # opened to answer: ClinVar's own definition is "ClinVar does not have the appropriate
    # term for your submission", and the table's reason for excluding it is that its meaning
    # depends on a free-text explanation this app never receives. That is a statement that
    # no usable classification is available — which is what this class says, and is not a
    # severity. It shares the class with `not_provided` rather than being mapped *to* it:
    # #88 refused the latter, correctly, because it would have said one ClinVar value was
    # another. A class named for the class does not make that claim.
    "not_provided": "No_Classification",
    "other": "No_Classification",
    "conflicting_data_from_submitters": "No_Classification",
    # `-` is a *term* in the table — "not submitted as a classification; the variant was
    # submitted only as part of a haplotype" — and not the same thing as an empty cell, which
    # this module spells `.`, `nan`, `None` or blank and reports as No Clinical Data. So the
    # two render differently on purpose: `-` says ClinVar holds a record with no call for
    # this variant alone, and an empty cell says ClinVar was not consulted. Worth knowing
    # that some writers use `-` as a plain missing-marker, in which case this reads it more
    # generously than they meant; no measured file does, and the table is the authority here.
    "-": "No_Classification",
    # Two spellings measured in real MAFs that the table does not list, added with the same
    # reasoning and flagged as additions beyond it. Both are ClinVar aggregate values whose
    # names state that no classification is available; the first is how the VCF spells the
    # table's `-` row, and is the second most common unreadable value in the corpus.
    "no_classification_for_the_single_variant": "No_Classification",
    "no_classifications_from_unflagged_records": "No_Classification",
    # `other` as the split leaves it when a call preceded it, the same shape #99 added
    # `_drug_response` for. #99 put both spellings on offer, so the same contract reaches
    # here — and unlike `_drug_response` this one *is* observed: `…|_other` occurs in seven
    # of the measured files. It still rescues no row on its own, because the split keeps the
    # call rather than the modifier; it is here because the term is offered.
    "_other": "No_Classification",
    "Conflicting_classifications_of_pathogenicity": "Uncertain_Significance",
    # ClinVar's name for the same call before its 2023 rename. Both spellings are offered
    # and both Broad presets keep both (issue #88), so without this entry the rows that
    # change admitted the summary reported "No Clinical Data" for a variant whose whole
    # distinguishing feature is that ClinVar holds conflicting data about it.
    "Conflicting_interpretations_of_pathogenicity": "Uncertain_Significance",
    # -- RENOVO: six classes that are six contiguous bins of one score (issue #108). --------
    # RENOVO is a machine-learning predictor, and `bin/add_renovo_to_maf.py` copies both its
    # `RENOVO_Class` and its `PL_score` (as `RENOVO_pls`) straight out of RENOVO's own output,
    # so the prefix has a numeric authority sitting in the same MAF as the class. Across the
    # **211,744 rows that carry both a class and a numeric score**, in the 64 MAFs on this
    # machine carrying `RENOVO_Class`, the six are non-overlapping intervals of that score:
    #
    #     HP Pathogenic  [0.88928 .. 1.00000]     LP Benign  [0.23501 .. 0.49985]
    #     IP Pathogenic  [0.78540 .. 0.88858]     IP Benign  [0.00920 .. 0.23494]
    #     LP Pathogenic  [0.50000 .. 0.78453]     HP Benign  [0.00000 .. 0.00919]
    #
    # So HP/IP/LP is a **confidence** on one prediction, not a second scale — and the two `LP`
    # classes are the *low-confidence pair either side of the 0.5 decision boundary*:
    # `LP Pathogenic` bottoms out at 0.50000 and `LP Benign` tops out at 0.49985. Both read as
    # Uncertain Significance here, which is the pipeline's own reading —
    # `bin/generate_clinical_summary.py` tiers both 3, with the comment "borderline /
    # conservative" — and it is what the mapping this replaced could not say, having no RENOVO
    # value landing on Uncertain Significance at all: it put the ACMG ladder's widest gap
    # (🟠 Likely Pathogenic against 🟢 Likely Benign) across RENOVO's narrowest score
    # difference, 0.00015.
    #
    # **Three of the six changed at #108, and the check is exact rather than argued.** That
    # same pipeline process writes its answer into the MAF as `variantalker_naive`, which
    # `compute_keep` keeps and the table shows, so the app can be compared against it row by
    # row. On the 52 MAFs whose column was written by a RENOVO-aware pipeline — 99,210
    # comparable rows — these six tiers agree on **every one**, where the six they replaced
    # left 974 rows (`LP Pathogenic` 662, `IP Pathogenic` 168, `LP Benign` 144) contradicting
    # a column drawn beside them on the same screen. See `tests/test_clinical_summary.py`.
    "HP Pathogenic": "Pathogenic",
    "IP Pathogenic": "Likely_Pathogenic",
    "LP Pathogenic": "Uncertain_Significance",
    "LP Benign": "Uncertain_Significance",
    "IP Benign": "Likely_Benign",
    "HP Benign": "Benign",
    # CIViC's five evidence levels used to sit here, mapping `A` to Pathogenic and `E` to
    # Benign. They are gone, and they are not rehoused in this dict under better values:
    # everything here is a *classification* a source asserted about the variant, and a CIViC
    # level is a statement about the evidence behind an assertion instead. It is read through
    # `_CIVIC_LEVEL_CLASSES` (issue #109), which is the second mapping keyed on a source's own
    # scale rather than on this one, `_ESCAT_GROUP_CIRCLES` being the first.
}

#: ClinVar terms that contain one of the separators the multi-value split uses.
#:
#: Both readers of a ClinVar cell split it on ``|``, ``/``, ``;`` and ``,`` to pull one
#: classification out of a composite. That is right for a composite — ``Benign/Likely_benign``
#: really is two calls joined — and wrong for these two, whose *own name* carries a comma:
#: the split turned ``Pathogenic,_low_penetrance`` into ``Pathogenic``, silently asserting
#: the equivalence the term table warns against, and making the mapping entries for them
#: unreachable. Checked before splitting, so an atomic term survives it.
#:
#: Deliberately this set and not "any value the mapping knows". Skipping the split whenever
#: the whole cell is mappable also resolves ``Benign/Likely_benign`` through its own entry,
#: which is a defensible reading — but it is a different one, and measured across 62 real
#: annotated MAFs it moves **2,088 rows** from Benign to Likely Benign. That is a clinical
#: reclassification seventeen times the size of the repair this ticket exists for, and not
#: issue #98's to make in passing. The inconsistency it would fix is real and left standing:
#: ``Benign/Likely_benign`` reads as ``Benign`` alone and as ``Likely_Benign`` when a
#: modifier follows it, because only then does the whole cell reach the lookup.
#:
#: **It still matches whole cells only, and since issue #285 that is enough.** The check
#: cannot be reached for a *piece*, because both terms in this set are also
#: :data:`CLINICAL_VALUE_MAPPING` keys, and :func:`main_classification` returns a piece that
#: is a known call before it could re-read it. So this set does exactly one job, at the top:
#: it stops a cell that is *nothing but* an atomic term being cut at its own comma. Without
#: it, ``Pathogenic,_low_penetrance`` splits to ``Pathogenic`` — which is a mapping key, so
#: the new check would return it, and the equivalence the term table warns against would be
#: asserted after all. Load-bearing, and easy to mistake for redundant.
#:
#: **On real data it fires on nothing, which is how a defect hid behind it for three
#: tickets.** A whole-cell atomic term occurs in **zero** of 159,573 annotated rows; the term
#: this set exists to protect always arrives embedded. The guard that covered it was written
#: against the whole-cell shape too, so both the protection and its test addressed a value
#: the corpus does not contain (issue #285).
#:
#: **The docstring here used to say the buried case "is reported as unrecognised — honest, and
#: less than the file says", and that was true of one surface and false of the other.**
#: ClinVar's own circle was honestly ⚪; the row's `Clinical_Summary` was not, because `Unknown`
#: is the bottom of :data:`CLINICAL_HIERARCHY` and any other source's `Uncertain_Significance`
#: outranked it. The row read 🟡 with nothing to say a cell had failed to parse — and the
#: *filter* kept the same row as pathogenic, since ``has_clinvar_term`` splits on ``[|/;,]`` in
#: one pass and did see the bare ``Pathogenic``. Recording the shape rather than only the fix:
#: a claim about what a value "is reported as" has to name which surface reports it.
#:
#: **The buried case does not resolve through this set, and that is the decision rather than a
#: shortfall.** ``Pathogenic/Pathogenic,_low_penetrance|other`` names *two* pathogenic calls —
#: one submitter's bare ``Pathogenic`` and another's qualified one — so the split takes the
#: first flagged piece and reports ``Pathogenic``, the stronger of two calls that are both
#: present. Widening this set to match inside a piece would have reported the weaker one.
#: Issue #98's warning is untouched by that: it is about a cell where low penetrance is the
#: **only** call, which is a different cell, and ``Benign;Pathogenic,_low_penetrance`` — where
#: ``;`` isolates the term — still resolves through this set and keeps the entry reachable.
_ATOMIC_TERMS_WITH_SEPARATOR = frozenset(
    {"Pathogenic,_low_penetrance", "Likely_pathogenic,_low_penetrance"}
)


def main_classification(value):
    """The one classification a ClinVar cell is read as.

    A cell can hold several calls joined by ``|``, ``/``, ``;`` or ``,``; this picks the
    first that mentions pathogenicity, else the first piece. An atomic term that contains a
    separator in its own name is returned whole — see
    :data:`_ATOMIC_TERMS_WITH_SEPARATOR`.

    **The piece it picks is read again when it cannot be named as it stands** (issue #285).
    The separators are tried in a fixed order, so splitting on the first one present can leave
    a piece that still holds another: ``Pathogenic/Pathogenic,_low_penetrance|other`` split on
    ``|`` yields ``Pathogenic/Pathogenic,_low_penetrance``, which was then looked up whole and
    matched nothing. Measured, that was 8 rows across 8 files reading as ``Unknown`` — and not
    visibly so, since a second source's weaker call outranked it on the same row.

    **The ``CLINICAL_VALUE_MAPPING`` check is what keeps that repair to the rows it is for,
    and it is not the same check #98 declined to make.** #98 refused to route a *whole cell*
    around the split, because ``Benign/Likely_benign`` would then resolve through its own entry
    and move thousands of rows from Benign to Likely Benign. This asks only of a piece the
    split has already produced, so a whole cell still splits exactly as it did. Recursing
    unconditionally was tried first and is wrong for the same reason in mirror image: measured
    over the corpus it also re-read the piece ``Benign/Likely_benign`` left by
    ``Benign/Likely_benign|other``, moving **68 rows across 64 files** the other way, from
    Likely Benign to Benign. Nine hand-picked cases and the pinned composites all missed it —
    the corpus did not.

    That leaves standing the inconsistency :data:`_ATOMIC_TERMS_WITH_SEPARATOR` records:
    ``Benign/Likely_benign`` reads as ``Benign`` alone and as ``Likely_Benign`` when a modifier
    follows it. Closing it is a clinical reclassification in one direction or the other, which
    is a decision for whoever asks for it and not a thing to make in passing while repairing a
    parse — the same call #98 made. ``test_a_composite_that_is_itself_a_known_call_is_not_re_read``
    fails if a later change makes it silently.

    The recursion cannot run away: a separator is removed at each step, so the value strictly
    shortens.

    One definition, one reading. Both derived columns split a cell, and they had the split
    written out separately, which is how they came to disagree about the substring that
    counts as a pathogenic mention: the summary matched ``"pathogenic"`` and the circles
    ``"athogenic"``. Issue #99 made that a *parameter* here rather than a discrepancy across
    two functions, and left choosing between them to issue #104 along with the glyphs.

    **The difference was inert, which is why removing it is a simplification and not a
    change.** ``"athogenic"`` is a substring of ``"pathogenic"`` and the piece is lowercased
    before either is looked for, so the two can only ever part over a value containing
    ``athogenic`` without ``pathogenic`` — no term in any of the six sources' vocabularies
    does, and neither does any of the 125 distinct values measured across 183 real annotated
    MAFs.

    That is a property of the *vocabularies*, not of this function, so it is the vocabularies
    a test holds: :func:`test_no_source_spells_a_call_that_would_have_split_the_two_readings`
    fails if a term arrives that would have made the old difference bite. Reinstating the
    parameter would not fail anything today, and that is the point — there was nothing to
    keep.
    """
    if value in _ATOMIC_TERMS_WITH_SEPARATOR:
        return value
    for sep in ["|", "/", ";", ","]:
        if sep in value:
            parts = value.split(sep)
            flagged = [p.strip() for p in parts if "pathogenic" in p.lower()]
            piece = flagged[0] if flagged else parts[0].strip()
            # Read again only if the piece cannot be named as it stands — see below.
            if piece in CLINICAL_VALUE_MAPPING:
                return piece
            return main_classification(piece)
    return value


#: What a CIViC cell holding no evidence items at all is spelled as.
#:
#: The annotation writes one entry per CIViC evidence item, so a variant CIViC knows and has
#: nothing on record for arrives as an empty list rather than as a blank. That is the
#: annotation saying *no evidence*, which is what ``⬜`` means since #95 — this source has no
#: call for this variant — and not what ``⚪`` means. Handled here rather than in this module's
#: general empty markers because it is one column's spelling of empty, and a bare ``[]`` in any
#: other column is a value nobody has explained.
#:
#: One spelling, the one the annotation writes. A whitespace variant was here and is gone: it
#: occurs in no measured file, and inventing spellings for a marker means the *unreadable* case
#: quietly shrinking, which is the case a reader most needs told apart from this one.
_CIVIC_CELL_WITHOUT_EVIDENCE = "[]"

#: The class each CIViC level is reported as — the second mapping keyed on a source's own scale.
#:
#: **No level maps to a benign class, and that is the decision rather than a rounding.** CIViC's
#: levels grade the study type behind an evidence item (:data:`~config.vocabularies.
#: CIVIC_DEFINITIONS`), and none of them asserts that a variant is harmless: ``E`` means the only
#: evidence is indirect, which is an absence of strong evidence and not a benign call. It used to
#: read as ``Benign`` through :data:`CLINICAL_VALUE_MAPPING`, so a variant whose CIViC evidence
#: was merely inferential drew ``🔵`` beside sources that had never said any such thing. The dev
#: chose to keep the CI position and grade it on strength, so the floor of this mapping is
#: ``Uncertain_Significance``: below ``B`` the honest reading is that CIViC does not yet know.
#:
#: ``C``, ``D`` and ``E`` therefore share a glyph, as the strip is deliberately coarse — the
#: same reason two classes may share one in :data:`_CLASS_GLYPHS`, though here it is one class
#: reached by three levels rather than one glyph lent to two classes. The distinction between a
#: case report and an in-vitro model is readable in the ``CIViC_Evidence_Level`` column itself,
#: and the help page defines each level.
#:
#: What this mapping cannot fix is that the glyphs it lands on are named for the pathogenicity
#: ladder — ``🔴`` is labelled *Pathogenic* in the key. :data:`CIVIC_CIRCLE_NOTE` is what the key
#: says about that, and it is the same compromise #100 struck for ESCAT.
_CIVIC_LEVEL_CLASSES = {
    "A": "Pathogenic",
    "B": "Likely_Pathogenic",
    "C": "Uncertain_Significance",
    "D": "Uncertain_Significance",
    "E": "Uncertain_Significance",
}

#: CIViC's levels strongest first, read off the cited table rather than written out again.
_CIVIC_LEVELS_STRONGEST_FIRST = tuple(CIVIC_DEFINITIONS)


def civic_levels(value):
    """Every CIViC evidence level a cell mentions, strongest first.

    **The cell holds one entry per CIViC evidence item, not one value.** The annotation writes
    the whole list — ``"['B', 'C', 'D']"`` — and every CIViC column on the row is a list of the
    same length, one element per evidence item. So a variant with thirty-five CIViC evidence
    items arrives with thirty-five levels in this column, and reading it as a single
    classification is reading the first character of a list literal: before issue #109 this went
    through :func:`main_classification`, which split on ``,`` and returned ``"['B'"`` — bracket,
    quote and all — a value no mapping holds, so **133 of 140 CIViC-annotated rows across 34
    real MAFs drew the circle for a value the app could not read**, beside the ⬜ the key defines
    as "No data".

    **The rule is the vendored filter's, deliberately.** ``vendor.pipeline_filters.
    has_element_from_list`` asks ``element in s`` — a plain substring test against the raw cell —
    so the pipeline's answer to *which levels does this cell mention* is exactly this, and CIViC
    filtering has always worked on the files where this column's reading failed. The two cannot
    share code: the matcher is in frozen ``vendor/``, and it answers a different question
    (*is any level I selected mentioned*) than the report does (*what is the strongest level
    mentioned*). What they can share is the rule, so it is applied here per level rather than
    reimplemented as a list parse.

    That is not the only reading that works — ``ast.literal_eval`` recovers the list — and the
    substring rule is chosen for two reasons. It agrees with the filter **by construction**
    rather than by coincidence, which is what came apart here; and it reads any shape the
    annotation writes, the bare ``A`` of one measured file included, where a list parse needs a
    branch per shape. Measured over all 140 CIViC-annotated rows of the corpus, the two agree on
    every cell — 0 disagreements — so nothing is lost by taking the one the filter defines.

    Returns ``[]`` for a cell mentioning no level at all, which is both an empty evidence list
    and a value this app cannot read; :func:`civic_reports_no_evidence` is what tells those
    apart.
    """
    return [level for level in _CIVIC_LEVELS_STRONGEST_FIRST if level in value]


def civic_class(value):
    """The class a CIViC cell is reported as, in three outcomes both readers need.

    * a class from :data:`_CIVIC_LEVEL_CLASSES`, for the strongest level the cell mentions;
    * ``None`` when the cell is the annotation saying *no evidence items* — see
      :data:`_CIVIC_CELL_WITHOUT_EVIDENCE`;
    * ``"Unknown"``, this module's existing name for a value it cannot read, when the cell
      mentions no level and is not that.

    The three are separated here rather than at each call site because **both** derived columns
    read this column, and the whole defect this repairs is the two of them reading it
    differently: the circle strip and ``Clinical_Summary`` each had their own copy of the split,
    which is how one came to say ⚪ where the other said a class. Each reader still decides what
    to *draw* for the three outcomes — ⬜ against ⚪ is a circle's distinction, and the summary
    simply skips a source that said nothing — but what the cell *says* is decided once.
    """
    if value.strip() == _CIVIC_CELL_WITHOUT_EVIDENCE:
        return None
    levels = civic_levels(value)
    return _CIVIC_LEVEL_CLASSES[levels[0]] if levels else "Unknown"


#: The annotation columns the summary reads, in the pipeline's own order.
#:
#: **The names are the pipeline's, and one of them was not.** Until issue #108 the fourth entry
#: read ``filter_renovo``, which is the name of a filter *parameter* — the RENOVO classes a run
#: keeps — and not a column any MAF carries: measured at 0 of 297 byte-distinct MAFs on this
#: machine. ``compute_keep`` emits ``RENOVO_Class`` and has never emitted ``filter_renovo``, so
#: the entry matched nothing and RENOVO reached this column on no row ever rendered. It read as
#: an unfinished rename rather than a decision, and the authority says so plainly: the pipeline
#: runs ``bin/generate_clinical_summary.py`` as a live annotation process, this function is a
#: re-implementation of it, and that file names ``RENOVO_Class`` in this exact slot.
#:
#: **A tuple, because the numbers it used to hold were dead.** This was a dict mapping each
#: column to a priority 1–5, and nothing read them: the priorities went into a local list that
#: was appended to and never inspected, while the label is chosen by walking
#: :data:`CLINICAL_HIERARCHY` over the classes collected. So no source outranks another —
#: the *strongest class asserted by any of them* wins — and issue #108's question of "what
#: priority should RENOVO have" had no answer because there was no priority. Deleted rather
#: than wired up: a number that has never decided anything is not a decision to preserve, and
#: the pipeline's own copy of this function carries the same dead dict.
#:
#: One consequence worth stating, since it bounds every change to this tuple: because the
#: hierarchy takes the strongest, **adding a source can only strengthen a label or leave it
#: alone.** It can never weaken one, and it can only shrink :data:`NO_CLINICAL_DATA`.
_SUMMARY_SOURCE_COLUMNS = (
    "CancerVar",
    "InterVar",
    "ClinVar_VCF_CLNSIG",
    "RENOVO_Class",
    "CIViC_Evidence_Level",
)


def generate_clinical_summary(row):
    """Name this variant with the strongest class any of its sources asserts.

    Reads the five columns of :data:`_SUMMARY_SOURCE_COLUMNS`, maps each value to a class in
    :data:`CLINICAL_HIERARCHY`, and renders the strongest. ESCAT is deliberately not among
    them (see the module docstring).

    **RENOVO reaches this column from issue #108, and did not before it.** What that costs is
    measured rather than estimated, because the pipeline's own answer to this question travels
    in the MAF as ``variantalker_naive``: on the 52 real germline MAFs whose column is
    RENOVO-aware — 99,210 comparable rows — this function agreed with it on 81,448 before the
    fix and on **all 99,210** after. Across the 218,349 rows of the 64 MAFs carrying
    ``RENOVO_Class`` the label moves on 44,828 — almost all of them 🔵 Benign → 🟢 Likely Benign
    or 🟡 Uncertain Significance — and ``🔍 No Clinical Data`` falls from **2,033 of those rows
    to 1**, the one survivor being a row where RENOVO is empty too. It also moves the grid's
    opening sort, which reads :func:`get_clinical_summary_sort_order` over exactly this label.

    Those three totals are three different populations of one corpus and are written out as such
    deliberately: 218,349 rows in the files, 218,348 of them carrying a RENOVO value, and 211,744
    carrying both a value and a numeric score. The suite pins the *mechanism* — that a RENOVO
    call alone names a class — and not the corpus counts, which no committed fixture can witness:
    ``germline_reference.maf`` has no ``🔍 No Clinical Data`` row either side of this change.

    **A value this module cannot read is no longer discarded** (issue #98). It used to be:
    an unmapped value became ``"Unknown"`` and was dropped, so a variant whose *only*
    annotation was one the mapping did not hold came out as :data:`NO_CLINICAL_DATA` —
    telling a clinician there was nothing on record about a variant that ClinVar had, for
    instance, called ``protective``. It is now kept and ranks last, so the column says the
    annotation is unrecognised rather than absent, and :data:`NO_CLINICAL_DATA` means what
    it says: no source annotated this variant at all.

    That repair is not ClinVar-only, because the drop was not. Measured across 62 real
    annotated MAFs, it moves **89 rows**: 82 carrying a ClinVar term this app could not
    read, and 7 carrying a ``CancerVar`` value of ``1``, which is an annotation artefact —
    and "unrecognised annotation" is true of that too, where "no clinical data" was not.
    """

    clinical_values = []

    for col in _SUMMARY_SOURCE_COLUMNS:
        if col in row and not says_nothing(row[col]):
            value = str(row[col]).strip()

            # Handle complex ClinVar values (extract main classification)
            if col == "ClinVar_VCF_CLNSIG":
                value = main_classification(value)

            # Map to a class. An unreadable value is kept as "Unknown" rather than dropped:
            # it is evidence that this source said *something*, which is what separates an
            # unrecognised annotation from an absent one.
            if col == "CIViC_Evidence_Level":
                # CIViC is read on its own scale, and nothing in this column reaches
                # `CLINICAL_VALUE_MAPPING` any more (issue #109). A cell holding no evidence
                # items is not an annotation at all, so this source said nothing and is skipped
                # -- which keeps `NO_CLINICAL_DATA` true of a variant whose only CIViC cell is
                # an empty list.
                mapped_value = civic_class(value)
                if mapped_value is None:
                    continue
            else:
                mapped_value = CLINICAL_VALUE_MAPPING.get(value, "Unknown")
            clinical_values.append(mapped_value)

    # No source carried a value at all — the one state this sentence is true of.
    if not clinical_values:
        return NO_CLINICAL_DATA

    for priority_level in CLINICAL_HIERARCHY:
        if priority_level in clinical_values:
            return _SUMMARY_LABELS[priority_level]

    return _SUMMARY_LABELS["Unknown"]


# Every annotation source the overview can draw a circle for, in circle order.
# (abbreviation, display name, MAF column).
#
# Public because the grid draws the key for this column from it (see
# ``variant_table._render_aggrid_with_detail``): the circles and the key explaining
# them have to name the same sources in the same order, and the only way to guarantee
# that is for both to read this list. It was private while the two sat in one file,
# which made a shared fact look like a local one.
#
# This is the *catalogue*, not what any one report draws — :func:`circle_sources`
# narrows it to the arm and the file. Nothing should read this list directly to decide
# what a report shows.
#
# The **display name** is here rather than in the key's markup because the key spelled
# "(ClinVar, InterVar, CancerVar, RENOVO, CIViC, ESCAT)" out by hand beside
# abbreviations it read from this list: half derived, half transcribed, and only the
# transcribed half free to drift.
#
# ``RN``'s column is ``RENOVO_Class``, with no fallback. It used to be ``filter_renovo``
# with ``RENOVO_Class`` behind it, and the fallback was doing all the work:
# ``filter_renovo`` is the name of a filter *parameter*, and a column of that name
# appears in 0 of 60 distinct MAFs on this machine. ``compute_keep`` emits
# ``RENOVO_Class`` and has never emitted ``filter_renovo``, so under
# :func:`circle_sources` the old primary could not be selected at all. Issue #95.
#
# ``generate_clinical_summary`` read ``filter_renovo`` too, with no fallback behind it, so
# RENOVO reached ``Clinical_Summary`` on no row ever rendered. Issue #108 fixed that at
# :data:`_SUMMARY_SOURCE_COLUMNS`, so **both derived columns now read RENOVO by the same name**
# — which is what makes the six RENOVO entries of :data:`CLINICAL_VALUE_MAPPING` shared rather
# than one column's private business, and why re-grading three of them there moved this
# column's ``RN`` glyph on 4,485 of the 218,348 rows carrying a RENOVO value as well.
CLINICAL_CIRCLE_SOURCES = [
    ("CV", "ClinVar", "ClinVar_VCF_CLNSIG"),
    ("IV", "InterVar", "InterVar"),
    ("CA", "CancerVar", "CancerVar"),
    ("RN", "RENOVO", "RENOVO_Class"),
    ("CI", "CIViC", "CIViC_Evidence_Level"),
    ("ES", "ESCAT", "ESCAT"),
]

#: The words the key writes beside each glyph — and, by its keys, which rows the key has.
#:
#: The glyphs are not here. Both this and the drawing read them off :data:`_CLASS_GLYPHS`,
#: which is what stops the key from explaining one glyph while the column draws another; the
#: key used to be a hand-written HTML string in ``variant_table`` naming six of the eight
#: glyphs the circles emit, and that is exactly how ``💊`` and ``⚪`` came to be drawn beside a
#: key that omitted them (issue #100). What is *not* structural is the other direction —
#: that every glyph drawn has a row here — because a class can be added to the glyph table
#: and left out of this one. :func:`test_the_key_names_every_glyph_the_circles_can_draw`
#: holds it, by drawing real values rather than by reading either dict.
#:
#: The order is the order the key lists them in: strongest first, then the classes that are
#: not points on that scale.
#:
#: Issue #98 added four of these. They are the ClinVar term classes that used to have no
#: glyph of their own and so drew ⚪ by fallthrough — ``protective`` and ``risk_factor``
#: among them, opposite claims rendered as the same blank. Because the key derives from this
#: dict, naming them here is what puts them in the key; nothing had to be transcribed.
#:
#: **A class is missing from here when it shares another's glyph**, and that is the whole
#: rule for which rows the key has — a key lists glyphs, so a glyph gets one row however many
#: classes draw it. The three absentees are the two low-penetrance classes, whose flag the
#: circles cannot carry and the ``Clinical_Summary`` words can, and ``Unknown``, which shares
#: ⚪ with ``No_Classification`` (see :data:`_UNREADABLE_CIRCLE`). Since issue #104 the words
#: on the shared row have to be true of every class that draws it, which is why ⚪'s row no
#: longer says "No classification": that is a claim about what the *source* recorded, and
#: **18 of the 1,005 positions** drawing ⚪ across the measured corpus are values this app
#: could not read, of which the source had recorded nothing of the kind. A small minority,
#: and enough — the words are either true of both classes or false of one of them.
_KEY_WORDS = {
    "Pathogenic": "Pathogenic",
    "Likely_Pathogenic": "Likely Path.",
    "Uncertain_Significance": "VUS",
    "Likely_Benign": "Likely Benign",
    "Benign": "Benign",
    "Disease_Risk": "Disease risk",
    "Drug_Response": "Drug response",
    "Association": "Association / trait",
    "Protective": "Protective",
    "No_Classification": "No usable classification",
}

#: Every glyph the key lists, with the words it lists it by. Derived from the shared table.
_CIRCLE_LEGEND = {name: (_CLASS_GLYPHS[name], words) for name, words in _KEY_WORDS.items()}

#: What a source with nothing to say about this variant draws: the column is absent, or its
#: value is one of this module's empty markers. Named because the key has to explain it and
#: it is the one glyph that is not a classification, so it is not in the map above.
_ABSENT_CIRCLE = ("⬜", "No data")

#: The name for what a value this module cannot classify draws: ``Unknown``'s glyph, ⚪.
#:
#: **Nothing reaches its circle through this constant any more.** Until issue #104 it was the
#: ``.get`` fallback in :func:`pathogenicity_circle_glyphs`, because ``Unknown`` was the one
#: class with no glyph of its own; it has one now, so an unreadable value is looked up in
#: :data:`_CLASS_GLYPHS` like every other class and the drawing has no special case at all.
#: What is left here is a *name* — for the tests that ask what an unreadable value draws
#: without wanting to know which class lends it the glyph.
#:
#: **The collision is kept, and it is now a decision rather than an inheritance** (issue
#: #104, the dev's call). ⚪ means both "the source declined to classify this" and "the
#: source said something this app could not read". What changed is that it is no longer
#: *borrowed*: ``Unknown`` has its own entry in :data:`_CLASS_GLYPHS` which happens to hold
#: the same glyph, so the sharing is written down where a reader can see both sides of it,
#: rather than expressed as this constant reaching into ``No_Classification``'s.
#:
#: Measured across 183 byte-distinct real annotated MAFs, 330,189 rows. On ``dev_marco`` ⚪
#: was drawn at **1,108 positions — 987 by a term that really is a no-classification record,
#: 121 by a value that could not be read**. What decided it against splitting them is the
#: second number's *shape*, not its size: after issue #98 there is no clinical term left that
#: this app declines to classify, so every one of the 121 was a value it could not **parse**,
#: and a glyph of its own would have been a permanent key row for other tickets' unfinished
#: parsing.
#:
#: **Issue #109 then proved that point by taking 103 of them away.** This branch is based on
#: it, and re-measuring over the same corpus with its CIViC reading underneath left **18**:
#: 10 rows of ``Pathogenic/Pathogenic,_low_penetrance|other`` (the composite limit
#: :data:`_ATOMIC_TERMS_WITH_SEPARATOR` recorded), 7 ``CancerVar`` values of ``1``, and one
#: ``ClinVar_VCF_CLNSIG`` of ``3``. Not one is a classification anybody issued.
#:
#: **Issue #285 then took the first of those three away**, which was the largest contributor:
#: that cell is read now, so what is left unreadable is the ``CancerVar`` ``1``\ s and the
#: single ``3`` — annotation artefacts in both cases. The conclusion is unchanged and better
#: supported: every remaining ⚪-as-unreadable is a value no source *issued*, so the key still
#: may not call it "No classification", because the words have to be true of all of them
#: rather than of most rows.
#:
#: **The totals above are left as the figures each ticket measured, and are not restated
#: here, because this corpus is not byte-stable between runs.** #285 counted that cell at 8
#: rows in 8 files where this paragraph records 10 in 10; two runs of one enumeration twenty
#: minutes apart, over the canonical rule, saw 187 files and then 161, because the OneDrive
#: mirror hydrates and dehydrates underneath it (the ``3`` is in a file that was present for
#: one run and absent for the other). What is durable is the *composition* — which values are
#: unreadable and why — so that is what the words beside the glyph are answerable to.
#:
#: **The two meanings have never shared a cell**: 0 rows of the corpus draw ⚪ from both
#: causes at once, so the collision is between rows and files rather than inside one strip.
#:
#: What the sharing costs is paid in :data:`_KEY_WORDS` instead. The row now reads "No usable
#: classification", which is true of a source that declined and of a value this app could not
#: read; "No classification" was false of the second, as "Not provided" before it was false
#: of both. The claim the key must not make is the one about *whose* limit it was.
_UNREADABLE_CIRCLE = _CLASS_GLYPHS["Unknown"]

#: An ESCAT level's circle, keyed by the group ESCAT itself puts the level in.
#:
#: Derived rather than pattern-matched. This branch used to read the tier as a *string* —
#: ``val.startswith("I") and not val.startswith("II")`` — which had two consequences neither
#: of them intended. ``V`` matched nothing and fell through to ⚪, a glyph the key did not
#: list, on the app's most reliably annotated gene: the only ``TP53`` row in
#: ``resources/escat_tiering.csv`` is a tier-``V`` wildcard for head and neck, so every
#: ``TP53`` SNP in an HNSCC sample lands there. And ``IV`` would have drawn ``🔴``, the
#: strongest circle, for ESCAT's *weakest* actionable level, unreachable only because that
#: same file happens to assign none — a fact about a data file, not a check this function
#: made.
#:
#: The groups are :data:`ESCAT_DEFINITIONS`' own ``group`` field, which is Table 2's
#: grouping, so a level's circle now follows from what ESCAT says the level *is*. Keying on
#: the group rather than the level is also what makes the mapping short enough to read: eight
#: levels, four groups, and the group is the part a circle can express.
#:
#: ``V`` draws ``💊`` because that is what its definition says — a matched drug produces
#: objective responses, without improving outcome. That is a statement about a drug, not a
#: point on a benign-to-pathogenic scale, and ⚪/⬜ would both have said the app had nothing
#: to report about a variant ESCAT calls a real target under combination development. It is
#: the same glyph ClinVar's ``drug_response`` draws, and that is not the ⚪ collision in a
#: new place: ⚪ shares one *position* between two meanings, while 💊 says the same thing —
#: *there is a drug story here* — at two positions the key names separately. Which drug
#: story is the ``ESCAT`` column's own to tell, and the help page defines all eight levels.
#:
#: **The ES circle is graded on a different axis from the other five, and the key says so.**
#: ESCAT ranks how actionable a *target* is; the rest classify how pathogenic a *variant* is.
#: The two disagree on real files, and not rarely: across the annotated MAFs measured for
#: issue #100, 34 of 48, 32 of 39 and 34 of 41 ESCAT-annotated rows draw 🔴 or 🟠 here while
#: every other populated source on the row says benign. ``SF3B1 p.V1219V`` — synonymous,
#: ClinVar *Benign*, InterVar *Benign*, RENOVO *HP Benign* — draws ``🔵🔵⬜🔵⬜🔴``, because
#: the tiering table matches gene plus tissue with a wildcard on the variant, so any SNP in
#: an ESCAT target gene inherits that gene's tier. Correct as an ESCAT statement, and not a
#: pathogenicity call. Keeping the circle was the dev's decision against removing it; saying
#: what it means is the part this mapping owes the reader.
_ESCAT_GROUP_CIRCLES = {
    "Ready for routine use": "Pathogenic",
    "Investigational": "Likely_Pathogenic",
    "Hypothetical target": "Uncertain_Significance",
    "Combination development": "Drug_Response",
}

#: The clause the key adds for the ES position, since it is not read off the same scale.
ESCAT_CIRCLE_NOTE = "ES grades how actionable the target is, not how pathogenic the variant is."

#: And the CI position's, which is not on that scale either — nor on ESCAT's.
#:
#: CIViC's levels rank the *evidence* behind an assertion about the variant, from a proven
#: association in human medicine to an indirect one. So the strip has three axes on it: four
#: positions classifying the variant, one grading a target, and one grading evidence. Two of the
#: six therefore need a sentence, and both are drawn only when their position is (issue #109).
CIVIC_CIRCLE_NOTE = "CI grades how strong the CIViC evidence is, not how pathogenic the variant is."

#: The positions whose glyphs are not read on the pathogenicity ladder, by abbreviation.
#:
#: A mapping rather than a special case in the grid, because the grid used to test for ESCAT by
#: hand — ``any(source[0] == "ES" …)`` — which is one transcription per exceptional source, in a
#: module that draws no circles. Which positions need a clause is decided here, beside the
#: mappings that make them exceptional.
_CIRCLE_AXIS_NOTES = {
    "ES": ESCAT_CIRCLE_NOTE,
    "CI": CIVIC_CIRCLE_NOTE,
}


def pathogenicity_circle_legend():
    """Every glyph the circles can draw, in key order, as ``(glyph, label)`` pairs.

    Read by the grid to write its key. Derived from the same constants the drawing reads,
    which is the point: the key it replaced was transcribed, and named six of the eight
    glyphs :func:`generate_pathogenicity_circles` emits.
    """
    return list(_CIRCLE_LEGEND.values()) + [_ABSENT_CIRCLE]


def circle_axis_notes(sources):
    """The clauses the key owes these positions, in circle order.

    A position graded on something other than pathogenicity needs a sentence saying so, and a
    position the column does not draw needs none — a MAF whose ESCAT column its arm does not
    emit has no ES circle for the clause to be about. See :data:`_CIRCLE_AXIS_NOTES`.
    """
    return [
        _CIRCLE_AXIS_NOTES[abbreviation]
        for abbreviation, _display, _column in sources
        if abbreviation in _CIRCLE_AXIS_NOTES
    ]


def circle_sources(sample_type, available_columns, skip_civic=False):
    """The sources the overview draws for this arm and this file, in circle order.

    A source is drawn when the arm's ``compute_keep`` emits its column **and** this MAF
    carries it — so the overview can only ever summarise a column the user can open in the
    table beside it and check.

    That used to read "the same test :func:`~config.columns.resolve_visible_columns`
    applies", and since issue #117 it is not the same test: the resolver applies this one
    **and then** subtracts
    :data:`~config.columns.PIPELINE_COLUMNS_THE_APP_REPLACES`. The guarantee above
    therefore holds because no circle source is in that list, which is a fact about the
    list rather than about this function — so it is pinned by
    ``test_no_circle_source_is_a_column_the_app_takes_off_the_screen`` rather than left to
    be true by luck. If a source is ever added to it, this function must subtract it too
    or stop claiming the reader can check the column.

    Both halves are load-bearing, and each rules out a different false claim. Measured
    across 167 distinct annotated MAFs on this machine (issue #95):

    * **the arm**, because a germline MAF may still *carry* ``CancerVar`` — 1 of 64 did.
      ``compute_keep``'s germline branch removes it from the output, so reading the file
      alone would draw a circle summarising a column that appears in neither the table
      nor the export, with nothing on screen to check it against.
    * **the file**, because ``compute_keep`` names most of these columns unconditionally.
      It emits ``InterVar`` for every germline run, and 2 of 64 germline MAFs do not carry
      it; the same holds for ``RENOVO_Class``, ``ESCAT`` and ``ClinVar_VCF_CLNSIG``.
      Reading the arm alone would name a source over a column of ``⬜``.

    **CIViC is not what the presence test is for**, though it looks like the obvious case:
    ``compute_keep`` guards those columns itself — it is the one filter input the pipeline
    checks — so an absent CIViC block is already gone from ``emitted``. It matters anyway,
    because it is why a somatic MAF with the annotation (34 of 103) draws a position that
    one without does not, on the same arm.

    ``skip_civic`` is taken for the same reason the grid takes it. The resolver is a pure
    function of the arm, that flag and the columns present, so dropping any one of the
    three is how this stops agreeing with the table it sits above — and the CIViC position
    is exactly the one the flag decides. It is ``False`` on every live path today,
    ``parameter_config`` stripping it from loaded parameter files as a deprecated name,
    but a caller that passes ``True`` must get the answer the grid would give.

    Returns:
        A sublist of :data:`CLINICAL_CIRCLE_SOURCES`, order preserved. May be empty: a
        MAF carrying none of its arm's sources gets no overview column at all, which is
        :func:`add_clinical_summary_column`'s branch rather than this one's.
    """
    emitted = set(
        pipeline_output_columns(
            sample_type, skip_civic=skip_civic, available_columns=available_columns
        )
    )
    present = set(available_columns)
    return [
        source
        for source in CLINICAL_CIRCLE_SOURCES
        if source[2] in emitted and source[2] in present
    ]


def pathogenicity_circle_glyphs(row, sources):
    """One glyph per entry of ``sources``, in that order, as a list.

    Read this rather than indexing the joined column value. **A glyph is not always one code
    point**: ``⚠️`` and ``🛡️`` are a symbol plus U+FE0F, so from issue #98 the strip's length
    stopped equalling the number of sources, and ``circles[2]`` stopped meaning "the third
    source" the moment a two-code-point glyph appeared before it. Iterating the string has
    the same fault — it yields the variation selector as if it were a glyph of its own.

    The joined string is still what the column stores; this is what anything asking about a
    *particular* source uses.

    The glyphs and what they mean are :func:`pathogenicity_circle_legend`'s, so this
    docstring does not list them: it used to, and its list had drifted from the key beside
    the grid, which had drifted from the glyphs below (issue #100).

    ``sources`` comes from :func:`circle_sources`, so every column named in it is one this
    frame carries. That is what makes ``⬜`` mean one thing here: *this source has no call
    for this variant*. It used to mean that **or** "this arm never consults this source",
    the second being a whole column of ``⬜`` the reader had to scroll to identify — two
    facts a clinician weighs differently, drawn identically (issue #95).
    """
    circles = []
    for entry in sources:
        label, _display, col = entry
        val = None
        if col in row.index:
            v = row[col]
            if not says_nothing(v):
                val = str(v).strip()

        if val is None:
            circles.append(_ABSENT_CIRCLE[0])
            continue

        # A CIViC cell lists one level per evidence item, so it is read as a *set* of levels
        # and drawn for the strongest — the reading the vendored filter already uses, which
        # is why filtering worked on the files where this circle did not (issue #109). An
        # empty list is CIViC holding no evidence, which is ⬜'s meaning and not ⚪'s.
        # A cell this app cannot read draws ⚪, which is #109's question and #104's answer:
        # `Unknown` is a class with a glyph of its own now, and that glyph is the one
        # `No_Classification` draws. So this needs no fallback and states no exception — it
        # is the same lookup the other five positions make.
        if label == "CI":
            mapped = civic_class(val)
            circles.append(_ABSENT_CIRCLE[0] if mapped is None else _CLASS_GLYPHS[mapped])
            continue

        # Handle ClinVar multi-value (e.g. "Benign/Likely_benign") — the same reading the
        # summary makes of the same cell, since issue #104.
        val = main_classification(val)

        # An ESCAT level's circle follows its group in `ESCAT_DEFINITIONS`, not the shape of
        # its name. A level the scale does not define — `IV`, `X`, or an annotation this app
        # cannot read — is unreadable rather than silently the strongest tier.
        if label == "ES":
            meaning = ESCAT_DEFINITIONS.get(val)
            mapped = _ESCAT_GROUP_CIRCLES[meaning.group] if meaning else "Unknown"
        else:
            mapped = CLINICAL_VALUE_MAPPING.get(val, "Unknown")

        # Looked up rather than defaulted. Until #104 an unreadable value reached its circle
        # through a `.get` fallback, because `Unknown` was the one class with no glyph of its
        # own; it has one now, so every class this reaches is a key here and a default could
        # only hide a class someone forgot to name.
        # `test_every_class_a_value_can_map_to_has_a_glyph` is what makes that safe.
        circles.append(_CLASS_GLYPHS[mapped])
    return circles


def generate_pathogenicity_circles(row, sources):
    """The ``Pathogenicity_Overview`` cell: one circle per source in ``sources``, joined.

    What the column stores and the grid shows. To ask about one source, read
    :func:`pathogenicity_circle_glyphs` instead — see there for why the joined form cannot
    be indexed.
    """
    return "".join(pathogenicity_circle_glyphs(row, sources))


#: Rendered label -> rank, derived from the hierarchy rather than written out beside it.
#:
#: This was a hand-maintained second copy of the label list, and it had already fallen behind
#: the first: every class the ``if``/``elif`` chain rendered through its ``else`` branch was
#: missing here and sorted to the default rank. Deriving it means a class cannot be added to
#: :data:`CLINICAL_HIERARCHY` and left unsortable.
#:
#: The two rows that are not classes come last, in the order they were already in: a variant
#: with an unreadable annotation sorts above one with no annotation, which sorts above a row
#: whose derivation raised.
_SORT_ORDER = {
    label: rank
    for rank, label in enumerate(
        [_SUMMARY_LABELS[name] for name in CLINICAL_HIERARCHY] + [NO_CLINICAL_DATA, ANALYSIS_ERROR]
    )
}


def get_clinical_summary_sort_order(clinical_summary_value):
    """Get sort order for clinical summary values (lower numbers = higher priority)."""

    # Anything unrecognised sorts after everything named, as it always has.
    return _SORT_ORDER.get(clinical_summary_value, len(_SORT_ORDER))


def add_clinical_summary_column(data: pd.DataFrame, sources) -> pd.DataFrame:
    """Add the clinical columns as the frame's first columns, sorted by significance.

    Args:
        data: the labelled frame, every column the MAF carried still on it.
        sources: what the overview draws, from :func:`circle_sources`. **Empty means no
            overview column is added at all** — an overview of nothing is a pinned
            180px column of empty strings, and ``Clinical_Summary`` already says
            "🔍 No Clinical Data" for that MAF, which is the honest sentence and needs
            no blank column beside it. The resolver drops absent app extras silently,
            so nothing downstream has to be told. Issue #95.
    """

    if len(data) == 0:
        return data

    # Create a copy to avoid modifying original data
    data_with_summary = data.copy()

    try:
        # Generate clinical summary for each row
        data_with_summary["Clinical_Summary"] = data_with_summary.apply(
            generate_clinical_summary, axis=1
        )

        # Generate pathogenicity overview circles
        if sources:
            data_with_summary["Pathogenicity_Overview"] = data_with_summary.apply(
                generate_pathogenicity_circles, axis=1, args=(sources,)
            )

        # Sort by clinical significance priority (pathogenic first, then likely pathogenic, etc.)
        data_with_summary["_clinical_sort_order"] = data_with_summary[
            "Clinical_Summary"
        ].apply(get_clinical_summary_sort_order)
        data_with_summary = data_with_summary.sort_values(
            ["_clinical_sort_order", "Hugo_Symbol"], ascending=[True, True]
        )
        # Remove the temporary sort column
        data_with_summary = data_with_summary.drop(columns=["_clinical_sort_order"])

        return _clinical_columns_first(data_with_summary)

    except Exception as e:
        # If clinical summary generation fails, return original data with a default column
        st.warning(f"⚠️ Could not generate clinical summary: {str(e)}")
        data_with_summary["Clinical_Summary"] = ANALYSIS_ERROR
        if sources:
            data_with_summary["Pathogenicity_Overview"] = "—"

        return _clinical_columns_first(data_with_summary)


def _clinical_columns_first(data: pd.DataFrame) -> pd.DataFrame:
    """``data`` with whichever clinical columns it has moved to the front, in order.

    Both paths through :func:`add_clinical_summary_column` end here, and both used to
    spell the reorder out. Neither may name ``Pathogenicity_Overview`` unconditionally
    now that a MAF carrying none of its arm's sources does not get one.
    """
    leading = [
        col
        for col in ("Clinical_Summary", "Pathogenicity_Overview")
        if col in data.columns
    ]
    rest = [col for col in data.columns if col not in leading]
    return data[leading + rest]


def strip_summary_glyph(text):
    """Return a ``Clinical_Summary`` label without its leading glyph.

    Shared by the summary tab's counts and the dashboard's clinical-significance axis, which
    each had their own copy — one slicing two characters off the front, the other matching
    against a hand-written tuple of nine emoji.

    Both broke on the glyphs issue #98 added, and in different ways. ``⚠️`` and ``🛡️`` are two
    code points each (the symbol plus U+FE0F), so ``label[2:]`` left a leading space; and a
    glyph absent from the tuple was not stripped at all, which on the chart meant the label
    no longer matched ``CLINICAL_COLORS`` and lost its place in the axis order. Splitting on
    the first space is blind to how many code points the glyph is, so a new one cannot
    reintroduce either failure.
    """
    if not isinstance(text, str):
        return text
    glyph, sep, rest = text.partition(" ")
    if sep and glyph and not glyph[0].isalnum():
        return rest
    return text


def format_clinical_summary_display(clinical_summary_counts):
    """Format clinical summary counts for display."""

    formatted_counts = {}

    for clinical_sig, count in clinical_summary_counts.items():
        clean_sig = strip_summary_glyph(clinical_sig)

        # Format for better readability
        clean_sig = clean_sig.replace("_", " ")
        formatted_counts[clean_sig] = count

    return formatted_counts
