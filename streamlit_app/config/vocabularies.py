"""The values a filter control may offer, and how a stored value is cleaned against them.

Two things that only make sense together: the vocabularies themselves — the arms, the
classification exclude list and the six guideline sources' keep-lists — and the
normaliser that turns whatever a cache, an upload or a hand-edited JSON holds into a
clean list of terms from one of them.

Two readers, not one: the parameter page renders a control from each list, and the help page
documents each list. ``page_modules/data_loading.py`` looked like a third — it imported all
seven — but it used none of them and those imports are gone. So the case for a module of its
own is not "many consumers"; it is that both readers want the *same* list and neither is its
owner, and that anything living in a page is unreadable to a bare ``pytest``.

Streamlit-free, and that is load-bearing. ``config/`` and ``filters/`` are importable
without Streamlit so the parity harness and a bare ``pytest`` can read the app's
vocabulary without booting a UI; :func:`validate_multiselect_params` and
:func:`filter_terms` moved here out of ``page_modules/parameter_config.py`` on that
basis, and nothing in this module may acquire an ``st`` call.
"""

from typing import NamedTuple

from config.param_migration import LEGACY_SENTINEL

#: The two arms, as the arm selector offers them.
#:
#: A vocabulary rather than a preset, which is why it is here and not in ``config/presets.py``
#: where it first landed: it is the set of values ``sample_type`` may take, exactly as
#: :data:`ESCAT_OPTIONS` is the set ``filter_escat`` may take. ``pipeline_params.ARMS`` is the
#: contract's own copy of the same two names, and ``tests/test_param_contract.py`` asserts the
#: two agree — a copy with a guard, per the rule that copies need one.
SAMPLE_TYPES = ["somatic", "germline"]

# Variant classifications the exclude control can offer.
#
# No catch-all sentinel, here or in any of the six guideline lists below (issue #36).
# `All` read as "no restriction", but was measured to be identical to the bug it looked
# like the way out of: three guideline sources match a listed value in 100% of reference
# rows, so one source reading as unrestricted made the OR true for every row — 68,364
# somatic rows against a parity baseline of 411. An *empty* multiselect is the pipeline's
# own empty keep-list, and `isin([])` drops that source out of the union instead.
#
# This list is the whole of the vocabulary the control can offer, and it is deliberately
# **not** the whole of what a MAF can carry: 211 reference rows hold five classifications
# missing from it, two of them minted by the pipeline itself. That is survivable only
# because the parameter is an *exclude* list — the unknown ones are kept, as the pipeline
# keeps them. It could never be fixed by listing more values, because the parameter page
# runs before any MAF is loaded, so the vocabulary can never be data-driven.
VARIANT_CLASSIFICATIONS = [
    # Protein-coding variants
    "Missense_Mutation",
    "Nonsense_Mutation",
    "Silent",
    "Frame_Shift_Del",
    "Frame_Shift_Ins",
    "In_Frame_Del",
    "In_Frame_Ins",
    "Nonstop_Mutation",
    "Translation_Start_Site",
    # Splicing and regulatory region variants
    "Splice_Site",
    "Splice_Region",
    "5'UTR",
    "3'UTR",
    "5'Flank",
    "3'Flank",
    "Intron",
    # Non-coding / other
    "IGR",
    "RNA",
]

# The six guideline sources' options. Each feeds a *keep* list: a value listed here is
# kept, an empty list places no constraint. There is no catch-all — see the note on
# VARIANT_CLASSIFICATIONS above for what it cost.
#
# `CIVIC_RATING_OPTIONS` used to sit between CIViC and ClinVar, offering evidence ratings
# 1-5. It was imported and never used, and there is no parameter behind it: the
# pipeline's only CIViC setting is `filter_civic_evidence_level`. Deleted with the
# sentinel, as a control that cannot affect the report.
CANCERVAR_OPTIONS = [
    "Tier_I_strong",
    "Tier_II_potential",
    "Tier_III_Uncertain",
    "Tier_IV_benign",
]
CIVIC_OPTIONS = ["A", "B", "C", "D", "E"]

#: The ClinVar pathogenicity calls — the first of two controls, and the head of the
#: vocabulary :data:`CLINVAR_OPTIONS` assembles.
#:
#: This block documents all three constants below it, ClinVar's vocabulary being one subject.
#: Every term in any of them is **atomic** — a term the filter can actually match.
#:
#: ``ClinVar_VCF_CLNSIG`` is a composite field: a variant carrying two calls is written
#: ``Pathogenic/Likely_pathogenic``, and one carrying a call plus a modifier is written
#: ``Benign|_risk_factor``. ``vendor.pipeline_utils.has_clinvar_term`` splits *that value*
#: on ``[|/;,]`` and asks whether any resulting piece is in this list — it does not split
#: the list. So a keep-term containing a separator can never equal a piece, and offering
#: one is offering a control that cannot fire.
#:
#: This list carried two of them until issue #88: ``Pathogenic/Likely_pathogenic`` and
#: ``Benign/Likely_benign``. They are gone, and the removal changed no row anywhere,
#: because a value spelled ``Benign/Likely_benign`` is kept by ``Benign`` or by
#: ``Likely_benign`` — the atoms it splits into, both still on offer. Measured before
#: removal on 8,859 rows of real annotated MAF: ``Benign/Likely_benign`` occurs **344
#: times as a value** and matched **zero** rows as a keep-term. That is the shape of the
#: bug — the term looks like the most obviously useful thing in the list.
#:
#: **Both spellings of the conflicting-classification term are offered, deliberately.**
#: ClinVar renamed ``Conflicting_interpretations_of_pathogenicity`` to
#: ``Conflicting_classifications_of_pathogenicity`` in 2023. Which one reaches a MAF is a
#: property of the ClinVar release the annotation step used, and this repo pins no such
#: release — so both are live spellings of one classification depending on when a file was
#: annotated. Offering only the current name is what the app did before #88, and on every
#: real MAF measured there it kept nothing while 79 rows carried the old name. Neither
#: spelling is safe to drop; a user cannot re-annotate their file to suit a dropdown.
#:
#: **The underscore-prefixed spellings are offered too, for the same reason** (issue #99).
#: A modifier that follows a call is written with a leading underscore in some files and
#: without one in others — ``Conflicting_interpretations_of_pathogenicity|_other`` against
#: ``Conflicting_interpretations_of_pathogenicity|other``, the same ClinVar content twice.
#: The split leaves the prefix on, so ``_other`` and ``other`` are two pieces and a keep-list
#: naming one of them misses the other.
#:
#: Measured over 48 byte-distinct real annotated MAFs, 16,776 ClinVar-annotated rows: **no
#: file carries both spellings of** ``other``. Two files carry ``_other`` 7 times each (14
#: rows) and ``other`` never; seven others carry ``other`` 1–3 times (14 rows) and ``_other``
#: never. So on those two files the offered ``other`` kept **nothing at all** — not a partial
#: shortfall but the whole kind #88 removed, arriving by a different route, and dead exactly as
#: ``Conflicting_classifications_of_pathogenicity`` is dead on a pre-2023 file.
#: ``drug_response`` is the partial case: 19 bare against 2 underscored in the same file.
#:
#: **The limit of that evidence, stated rather than left to be found:** those two files are two
#: annotation runs of *one sample*, so every underscored row measured comes from one patient's
#: MAF. What keeps this from resting on that sample alone is the pipeline's own
#: ``--filter_clinvar`` help, which advertises ``_low_penetrance`` — an underscored piece from a
#: vocabulary nobody here chose — so the artefact is one the pipeline's authors met
#: independently. A wider corpus would settle how common it is; it cannot make it go away.
#:
#: This is a repair rather than a leak, and specifically **not** a divergence from the
#: pipeline. The term the app filters with is the term it writes into the parameter file, and
#: the pipeline given that file runs this very matcher over the same value — so they agree
#: exactly, by construction rather than by coincidence. What would break decision 3 of map
#: #50 is splitting the *keep-term*, which #88 ruled out and which stays ruled out: the
#: matcher is in ``vendor/``, a byte-for-byte copy of the pipeline's code, and is not edited
#: here. Nor is the list guessing at a vocabulary this repo does not pin: it offers a second
#: spelling only for terms already on offer, each one measured in real data.
#:
#: **Nothing is owed upstream beyond this record.** A pipeline user typing
#: ``--filter_clinvar other`` has the identical shortfall, and that is the pipeline's to keep or
#: fix: changes to it are out of map #50's scope, and its help text shows its authors already
#: met these pieces and chose to expose them raw. So no issue is filed against it and no
#: behaviour here depends on one being; what this map owes is that the finding is written where
#: the list lives, which is here.
#:
#: **The list is now the whole of what the report can name** (issue #103), and that is the
#: rule rather than a coincidence. Until #103 the vocabulary held ten classifications while
#: ``ClinVar_VCF_CLNSIG`` carried twenty-odd, so a term could be *rendered* by
#: ``components.clinical_summary`` and be unselectable here — the mirror image of #88, where a
#: term was selectable and unrenderable. Both are the same defect: two ClinVar vocabularies in
#: one app, disagreeing. There is one now, and ``tests/test_app_defaults.py`` holds it in
#: **both** directions.
#:
#: **The authority is the institute's own ClinVar term table**, by way of issue #98, which
#: adopted it and classed every term with it — not the pipeline's ``--filter_clinvar`` help,
#: which #103 read and found is not a vocabulary anyone chose: it advertises an **empty term**
#: and a bare ``_low_penetrance``, neither of which is a classification at all, while omitting
#: ``Conflicting_interpretations_of_pathogenicity``, which #88 measured live in 498 rows.
#: Mirroring it would have imported those defects and still missed real data. Nothing here is a
#: clinical claim of this app's: a term is offered because the table names it and the matcher
#: can match it.
#:
#: What is **not** an argument against that help is that some of its terms are unobserved here,
#: and the first draft of this note made it before stating the opposite rule two paragraphs
#: below. ``Established_risk_allele`` occurs in no measured row and is offered, deliberately —
#: see the next paragraph. The help's problem is that it lists things that are not terms, not
#: that it lists terms this corpus lacks.
#:
#: **Reach does not decide membership; structure does.** A term is offered when the report can
#: name it *and* it can survive the split — so the four terms whose own names carry a
#: separator (``Pathogenic/Likely_pathogenic``, ``Benign/Likely_benign``,
#: ``Pathogenic,_low_penetrance``, ``Likely_pathogenic,_low_penetrance``) stay out, being
#: #88's kind of dead, while ``Established_risk_allele`` at **zero** observed rows and ``-``
#: at zero are in. That asymmetry is the line #88 and #99 already drew and is worth restating,
#: because it is easy to collapse: #88 deleted the composites as impossible *on every file for
#: ever*, and kept ``Conflicting_classifications_of_pathogenicity`` at zero rows because this
#: corpus is 2022–2024 vintage and the term is live upstream. Absence from files of this age is
#: not death.
#:
#: **What #103 measured, on a wider corpus than #99's.** An unbounded sweep of ``$HOME``,
#: deduplicated by content hash and excluding agent fixtures, reaches **176 byte-distinct MAFs
#: with a populated** ``ClinVar_VCF_CLNSIG`` **and 159,580 annotated rows** — against the 48
#: files and 16,776 rows recorded above, which came from a ``-maxdepth 6`` sweep. Verified
#: through ``vendor.pipeline_utils.has_clinvar_term`` itself rather than a re-implementation of
#: the split: **1,215 rows could be kept by no selection the old list offered**, 0.76% of the
#: corpus, and the two arms were hit alike (germline 0.73%, somatic 0.84%).
#:
#: The single largest term in that shortfall was one the ticket did not know about:
#: ``no_classification_for_the_single_variant`` at **854 rows across 102 files**, seven times
#: the reach of ``risk_factor`` and four times that of the ``other`` which #88 and #99 each
#: spent a ticket on. It is also the sharpest argument that the old list was a transcription
#: and not a shortlist: the app offered ``not_provided`` while refusing this, and both are
#: ClinVar saying it holds no classification for the variant — the same class, and now the
#: same treatment. The rest, for the record: ``risk_factor`` 247 rows in 102 files,
#: ``association`` 152 in 93, ``protective`` 107 in 92, ``Affects`` 35 in 24,
#: ``confers_sensitivity`` 15 in 5.
#:
#: **``_low_penetrance`` is deliberately not offered**, and it is the one place #99's
#: every-spelling rule is not applied. It is not a classification — it is the tail of
#: ``Pathogenic,_low_penetrance``, left behind because that term's own name carries a comma —
#: so there is no honest class for it and ``CLINICAL_VALUE_MAPPING`` does not name it.
#:
#: **It costs nothing to refuse, and the fact is sharper than "it sits beside a real call"**
#: (issue #278 measured it; #284 re-measured and tightened this sentence, which used to claim
#: only the weaker thing). Every occurrence in the corpus is the **same single cell value** —
#: ``Pathogenic/Pathogenic,_low_penetrance|other`` — carried **once each by 8 real MAFs**, and
#: on all 8 it is **one variant**: ``SERPINA1`` chr14:94847262 T>A, the same row
#: :func:`filters.variant_filters.frequency_mask`'s pathogenic exemption is documented against.
#: (The two committed parity fixtures carry the value once each as well, at constructed loci of
#: their own — ``NOTCH3`` and ``PDS5B`` — which is why a sweep including them reports 10.) So
#: there is no cell anywhere in which this piece is the whole value, and refusing it cannot lose
#: a row: the ``Pathogenic`` that shares every one of those cells is offered and is what keeps
#: them. ``Likely_pathogenic,_low_penetrance`` occurs in **zero** rows, so its tail is not even
#: reachable. Measured over 198 byte-distinct real MAFs, 183 with a populated
#: ``ClinVar_VCF_CLNSIG``, 159,819 annotated rows.
#:
#: The pipeline's help advertises it; that is the vocabulary nobody chose, not a reason. Issue
#: #278 revisited the whole question — the dev first chose to give it a class, then reversed once
#: the measurement above was in — and the outcome is that it stays out of this list, out of
#: ``CLINICAL_VALUE_MAPPING``, and out of every preset. Its decision is not a default nobody
#: examined.
#:
#: ``-`` is offered because the table names it as a term (*not submitted as a classification;
#: the variant was submitted only as part of a haplotype*) and the matcher can match it. Worth
#: knowing that some writers use ``-`` as a plain missing-marker, in which case selecting it
#: would keep their unannotated rows; no measured file does, and #98 drew this line first when
#: it gave ``-`` a class distinct from an empty cell.
#:
#: Split in two because the whole vocabulary in one multiselect buries ``Pathogenic`` (#103).
#: The boundary is **#98's class boundary, which is the source's and not this app's**: these
#: are the classifications, and :data:`CLINVAR_OTHER_ASSERTION_TERMS` is everything ClinVar
#: asserts that is not a pathogenicity call. Both controls write one ``filter_clinvar`` list,
#: so the parameter file the pipeline reads is exactly what it was — decision 3 of map #50
#: holds untouched, and the split is a fact about the screen only.
CLINVAR_PATHOGENICITY_TERMS = [
    "Pathogenic",
    "Likely_pathogenic",
    "Uncertain_significance",
    # Submitter-defined sub-tiers of Uncertain Significance. ACMG/AMP permits them and
    # defines no criteria for them, which is why #98 reads all three as one class rather
    # than ranking them, and why they are offered beside the term they subdivide.
    "VUS-high",
    "VUS-mid",
    "VUS-low",
    "Likely_benign",
    "Benign",
    "Conflicting_classifications_of_pathogenicity",
    "Conflicting_interpretations_of_pathogenicity",
]

#: Everything ClinVar asserts that is not a pathogenicity call — the second control.
#:
#: Five of #98's classes, in the order that module ranks them: disease risk, drug response,
#: association or trait, protective, no classification. A user who wants only pathogenic
#: variants never opens this control; a user chasing a pharmacogenomic or risk annotation
#: finds it in one place instead of eight rows below ``Benign``.
#:
#: Every term here is in ``CLINICAL_VALUE_MAPPING``, so the report has a name for each. The
#: four underscored spellings are #99's rule applied to the terms #103 adds: the split leaves
#: the prefix on, files differ over whether the modifier keeps it, and all four are observed
#: — ``_risk_factor`` 4 rows, ``_association`` 4, ``_protective`` 2, ``_confers_sensitivity``
#: 2. ``Affects|_association`` is a real cell that **neither** spelling of ``association``
#: reached before this.
CLINVAR_OTHER_ASSERTION_TERMS = [
    # Disease risk. `risk_factor` is the older OMIM-derived spelling of the same idea as the
    # three `*_risk_allele` terms; #98 gives all four one class on the table's authority.
    "risk_factor",
    "_risk_factor",
    "Uncertain_risk_allele",
    "Likely_risk_allele",
    "Established_risk_allele",
    # Pharmacogenomic.
    "drug_response",
    "_drug_response",
    "confers_sensitivity",
    "_confers_sensitivity",
    # Disease association or non-disease phenotype.
    "association",
    "_association",
    "Affects",
    # Protective effect — the opposite end of the axis from `risk_factor`, which is why #98
    # could not put the two in one bucket.
    "protective",
    "_protective",
    # No classification. `not_provided`, `other` and `_other` were already offered; the rest
    # of the class arrives with #103, and the class is the one both Broad presets already
    # keep — see `_SOFT_CLINVAR_TERMS` in config/presets.py.
    "not_provided",
    "other",
    "_other",
    "conflicting_data_from_submitters",
    "-",
    "no_classification_for_the_single_variant",
    "no_classifications_from_unflagged_records",
]

#: The whole ClinVar vocabulary, and the name every other module reads.
#:
#: One list from the two the controls draw, so anything asking *what may a ClinVar keep-list
#: hold* asks one name and needs not know it is rendered as two widgets. Concatenated rather
#: than sorted: the order is the order each control offers, and a reader comparing a saved
#: cache against this list should see the screen's order.
#:
#: **Its readers are few, and worth naming rather than gesturing at** — the first draft of this
#: note claimed ``filter_terms``, ``param_migration`` and the parity harness among them and
#: none of the three reads it. What does: the Help page's underscore caption, which walks it,
#: and the guards in ``tests/test_app_defaults.py``, which hold the vocabulary's contract with
#: ``CLINICAL_VALUE_MAPPING`` in both directions. ``filter_terms`` takes its options from the
#: caller, and since #103 the callers pass the two halves; ``param_migration`` names this
#: constant only in a comment.
#:
#: That makes it thinner than it looks, and it is kept for two reasons that are not "someone
#: reads it": splitting a vocabulary into two lists with no name for the whole is how the two
#: come to disagree about what a ClinVar term is, and the contract #103 exists to enforce —
#: *the control offers exactly what the report can name* — is a statement about the whole and
#: cannot be written without it.
CLINVAR_OPTIONS = CLINVAR_PATHOGENICITY_TERMS + CLINVAR_OTHER_ASSERTION_TERMS
INTERVAR_OPTIONS = [
    "Pathogenic",
    "Likely pathogenic",
    "Uncertain significance",
    "Likely benign",
    "Benign",
]
RENOVO_OPTIONS = [
    "LP Pathogenic",
    "IP Pathogenic",
    "HP Pathogenic",
    "LP Benign",
    "IP Benign",
    "HP Benign",
]
ESCAT_OPTIONS = ["IA", "IB", "IC", "IIA", "IIB", "IIIA", "IIIB", "V"]


#: The paper that defines the scale, so the app can cite its glosses rather than assert them.
#:
#: The levels are defined here; the *assignments* — which gene, in which cancer, gets which
#: level — come from disease-specific recommendations, which is what the pipeline's
#: `resources/escat_tiering.csv` transcribes: all 172 of its rows carry `official: yes` and
#: cite one of eight papers, seven of them in Annals of Oncology — ESMO's own journal — and
#: one in JCO Precision Oncology, with 105 rows from Mosele F et al., Annals of Oncology
#: 2020; 31(11):1491-1505 alone. So the institute reads the published scale, and this is its
#: defining table — not a second opinion about it.
ESCAT_SOURCE = (
    "Mateo J, Chakravarty D, Dienstmann R, et al. A framework to rank genomic alterations "
    "as targets for cancer precision medicine: the ESMO Scale for Clinical Actionability "
    "of molecular Targets (ESCAT). Annals of Oncology 2018;29(9):1895-1902."
)
ESCAT_SOURCE_DOI = "https://doi.org/10.1093/annonc/mdy263"


class ESCATLevel(NamedTuple):
    """One ESCAT level: the group it belongs to, what earns it, and what it implies.

    Three fields because Table 2 of :data:`ESCAT_SOURCE` gives three, and the third is the
    one a clinician reads for — a level is worth knowing because of what it says to do.
    """

    group: str
    evidence: str
    implication: str


#: What each level the control offers actually means, condensed from ESCAT's own table.
#:
#: This exists because the app used to gloss these levels from nowhere, and was wrong.
#: `IIIA`/`IIIB` were described as "Resistance mechanism with clinical/preclinical evidence"
#: and **ESCAT has no resistance tier at all** — Table 1 of the same paper attributes
#: resistance-biomarker grading to OncoKB, whose `R1`/`R2` this app does not read. `V` read
#: "Not actionable" on the help page, "case reports" in the parameter tooltip, and mapped to
#: `Unknown` in the Pathogenicity Overview: three surfaces, three stories, and ESCAT's `V` is
#: none of them. "Lack of evidence for actionability" is the scale's **`X`**, a different
#: level. #79 deleted the invented prose rather than rewrite eight clinical definitions from
#: memory; this is the sourced replacement (issue #89).
#:
#: Two of those three surfaces read this constant now. The third does not:
#: ``components/clinical_summary.py`` still maps ESCAT tiers into the pathogenicity hierarchy
#: by string prefix, so `V` still draws as `Unknown` — which this text now contradicts in
#: writing, since it says the drug is active. That mapping is issue #100's, deliberately: what
#: a level *means* and what a circle should *say* are different questions, and the second is
#: decided on a screen issue #95 is also editing.
#:
#: The keys are :data:`ESCAT_OPTIONS`, in that order, which is the paper's own — strongest
#: evidence first. Both facts are asserted rather than trusted, so a level cannot go
#: undescribed and a description cannot outlive its level.
#:
#: The scale also defines `IV` (preclinical evidence only) and `X` (no evidence of
#: actionability). Neither is described here because neither is offered, and neither is
#: offered because the annotation never carries one: `resources/escat_tiering.csv`, the only
#: source of the `ESCAT` column, assigns exactly these eight across its 172 rows. A control
#: offering `IV` or `X` could keep nothing — the shape of issue #88, not a gap to fill.
ESCAT_DEFINITIONS = {
    "IA": ESCATLevel(
        "Ready for routine use",
        "Prospective randomised trials show that matching a drug to this alteration in "
        "this tumour type gives a clinically meaningful improvement in a survival end "
        "point.",
        "Access to the treatment should be considered standard of care.",
    ),
    "IB": ESCATLevel(
        "Ready for routine use",
        "Prospective non-randomised trials show clinically meaningful benefit in this "
        "tumour type, as measured by the ESMO Magnitude of Clinical Benefit Scale 1.1.",
        "Access to the treatment should be considered standard of care.",
    ),
    "IC": ESCATLevel(
        "Ready for routine use",
        "Trials across several tumour types, or basket trials, show benefit for the "
        "alteration-drug match, of similar size in each tumour type.",
        "Access to the treatment should be considered standard of care.",
    ),
    "IIA": ESCATLevel(
        "Investigational",
        "Retrospective studies show clinically meaningful benefit from the matched drug "
        "in this tumour type, against patients who do not carry the alteration.",
        "Treatment to be considered preferable in the context of evidence collection — "
        "either a prospective registry or a prospective trial.",
    ),
    "IIB": ESCATLevel(
        "Investigational",
        "Prospective trials show the matched drug is more likely to get a response in "
        "this tumour type, but no survival data are available yet.",
        "Treatment to be considered preferable in the context of evidence collection — "
        "either a prospective registry or a prospective trial.",
    ),
    "IIIA": ESCATLevel(
        "Hypothetical target",
        "Benefit is established for this alteration, as tier I or II, but in a different "
        "tumour type — with limited or no clinical evidence in this one.",
        "Clinical trials to be discussed with the patient.",
    ),
    "IIIB": ESCATLevel(
        "Hypothetical target",
        "The alteration is predicted to act like a tier I abnormality in the same gene or "
        "pathway, but has no supporting clinical data of its own.",
        "Clinical trials to be discussed with the patient.",
    ),
    "V": ESCATLevel(
        "Combination development",
        "Prospective studies show the matched drug produces objective responses, but "
        "this does not lead to improved outcome.",
        "Trials of drug combinations could be considered.",
    ),
}

#: The strongest level, read off the ordering above rather than written out again.
#:
#: Two readers, and both used to spell `IA` out by hand: the parameter page's tooltip and the
#: help page's opening sentence, each saying which way the scale runs. That is one claim about
#: the ordering made on two surfaces, which is the shape that drifted before — so it is asked
#: rather than repeated.
ESCAT_STRONGEST = ESCAT_OPTIONS[0]


#: CIViC's own documentation of its evidence levels, so the app can cite them rather than
#: gloss them from memory — #89's line, paid a second time.
#:
#: The distinction this citation exists to hold is what a level *grades*. CIViC's levels rank
#: the **study type behind an evidence item**, from a proven association in human medicine down
#: to indirect inference. They are not a pathogenicity scale, and the app used to treat them as
#: one: ``A`` read as *Pathogenic* and ``E`` as *Benign* in the Pathogenicity Overview, which
#: says a variant is benign when what CIViC said is that the only evidence about it is indirect
#: (issue #109).
#:
#: What an evidence item *asserts* is a different field — ``CIViC_Entity_Significance``, which
#: holds ``RESISTANCE``, ``SENSITIVITYRESPONSE``, ``POOR_OUTCOME``. Measured over the rows of
#: real annotated MAF that still carry it, **784 of 960 evidence items assert resistance to a
#: therapy** and 4 assert likely pathogenicity; 26 carry a ``DOES_NOT_SUPPORT`` direction, so
#: they are evidence *against* their own assertion. And the pipeline's output keeps six CIViC
#: columns of forty, ``CIViC_Entity_Significance`` not among them — so on a filtered MAF the app
#: cannot know what CIViC asserted, only how strongly. That is the whole reason the level is
#: read as a strength and reported as one.
CIVIC_SOURCE = (
    "CIViC (Clinical Interpretation of Variants in Cancer) documentation, "
    "Evidence Level. Griffith M, Spies NC, Krysiak K, et al. CIViC is a community "
    "knowledgebase for expert crowdsourcing the clinical interpretation of variants "
    "in cancer. Nature Genetics 2017;49(2):170-174."
)
CIVIC_SOURCE_URL = "https://civic.readthedocs.io/en/latest/model/evidence/level.html"


class CivicLevel(NamedTuple):
    """One CIViC evidence level: its own name, and the study type that earns it.

    Two fields where :class:`ESCATLevel` has three, because CIViC's table gives two and the
    third has no counterpart: an ESCAT level implies a clinical action, while a CIViC level
    implies nothing on its own — what to do about the variant depends on what the evidence
    item asserted, which is a different column and one the pipeline's output drops.
    """

    name: str
    evidence: str


#: What each level the control offers means, in CIViC's own words.
#:
#: The keys are :data:`CIVIC_OPTIONS`, in that order, which is CIViC's own — strongest evidence
#: first — and both facts are asserted rather than trusted, so a level cannot go undescribed and
#: a description cannot outlive its level.
#:
#: This replaces a hand-written dict in ``page_modules/help.py``. That one was the last
#: transcribed vocabulary on its tab and it survived #79's check on the strength of being
#: *roughly* right, which was enough while nothing else read it and is not enough now that the
#: circles are graded from the same scale. Its ``A`` gloss — "Validated - Strong clinical
#: significance" — is the drift this constant exists to stop: CIViC's ``A`` is a proven
#: association, and "clinical significance" is the name of the other axis.
CIVIC_DEFINITIONS = {
    "A": CivicLevel(
        "Validated association",
        "A proven or consensus association in human medicine.",
    ),
    "B": CivicLevel(
        "Clinical evidence",
        "A clinical trial or other primary patient data supports the association.",
    ),
    "C": CivicLevel(
        "Case study",
        "Individual case reports from clinical journals.",
    ),
    "D": CivicLevel(
        "Preclinical evidence",
        "In vivo or in vitro models support the association.",
    ),
    "E": CivicLevel(
        "Inferential association",
        "An indirect association.",
    ),
}


#: Every filter parameter drawn from a fixed vocabulary: the classification exclude list
#: and the six guideline keep-lists.
#:
#: Private, which it could not be while it sat in a page beside its only reader. Nothing
#: outside this module names it — :func:`validate_multiselect_params` is what callers want,
#: and *which* parameters it sweeps is this module's business.
_MULTISELECT_PARAMS = (
    "filter_variant_classification",
    "filter_cancervar",
    "filter_civic",
    "filter_clinvar",
    "filter_intervar",
    "filter_renovo",
    "filter_escat",
)

#: The guideline sources each arm's report is actually built from — the clauses the
#: vendored ``somatic_filters`` and ``germline_filters`` OR together. ESCAT is absent from
#: germline because ``germline_filters()`` takes no ESCAT argument.
#:
#: Used for one thing: knowing when *every* source has been emptied, which is the state the
#: parameter page warns about. Written down rather than derived because deriving it means
#: importing the vendored signatures, which would pull pandas into every reader of the app's
#: vocabulary — including the bare-pytest guards that exist to run without it.
#: ``tests/test_param_contract.py`` derives the same relation from those signatures and is
#: what keeps this copy honest.
GUIDELINE_SOURCES = {
    "somatic": ("filter_cancervar", "filter_civic", "filter_escat", "filter_clinvar"),
    "germline": ("filter_intervar", "filter_renovo", "filter_clinvar"),
}


def validate_multiselect_params(params):
    """Normalise the shape of every multiselect parameter, and nothing else.

    It used to rewrite an empty selection to ``["All"]``, which is how the catch-all
    sentinel became reachable without a user ever choosing it. That is deleted (issue
    #36): an empty list is a real, expressible value — it is what the pipeline's own
    ``--filter_cancervar ""`` denotes — and the vendored code reads it as "this source
    places no constraint" rather than as "keep everything". Backfilling it was the
    opposite of neutral: one source reading as unrestricted made the guideline union true
    for every row.

    What is left is shape repair, which a hand-edited JSON still needs: a bare string
    becomes a one-element list, and anything that is neither becomes empty.

    It also strips the legacy sentinel, and this is the *only* place that catches every
    key. The widgets clean a value on its way in, so they never see a key the current arm
    does not render — a germline session leaves a stale ``filter_cancervar: ["All"]``
    untouched, and the auto-save at the end of the parameter page would write it straight
    back out to the cache and to any export the user downloads. This runs over the whole
    dict, after every tab.
    """
    for param in _MULTISELECT_PARAMS:
        if param in params:
            params[param] = filter_terms(params[param])

    # The panel dropdown deliberately gets no default here. It is UI state that resolves
    # to `filter_genes`, never a filter parameter of its own (issue #34) — see
    # gene_panel_state_key. Injecting it made the panel *name* part of every saved
    # parameter file and every cache entry.
    return params


def filter_terms(value, options=None):
    """One multiselect parameter's value, as a clean list of terms.

    The app's single normaliser for such a value. Three callers needed the same three steps —
    coerce to a list, drop the pre-#36 sentinel, and (for a widget) keep only what is on
    offer — and writing them out three times is how the three could come to disagree
    about what a parameter value is.

    ``options`` is what separates the two jobs. Without it this is *storage* hygiene: fix
    the shape and drop the sentinel, keeping every term the user meant, including ones no
    control offers. With it this is a *widget default*: also drop anything not on offer,
    because Streamlit raises on a default outside its options.

    Dropping the sentinel is behaviour-preserving rather than merely convenient, and only
    because issue #33 got there first — nothing branches on it any more, so ``["All"]``
    and ``[]`` reach the vendored code as the same restriction.
    """
    if isinstance(value, str):
        value = [value]
    if not isinstance(value, list):
        return []
    terms = [item for item in value if item != LEGACY_SENTINEL]
    if options is None:
        return terms
    return [item for item in terms if item in options]
