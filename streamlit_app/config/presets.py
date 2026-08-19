"""The four presets the parameter page offers.

Two on one axis — **Broad** and **Stringent** — each per arm. Nothing more: this module
holds the settings a user can *choose*, and none of them is what the app opens with.

**The settings the app opens with are not here.** They are the contract, in
``config/pipeline_params.py``, and both places that seed a session — ``MAFigate.py`` and
``page_modules/parameter_config.py`` — call ``pipeline_params("somatic")`` directly. This
module used to also export a ``DEFAULT_PARAMS`` naming itself the app's default, which no
app module ever read; issue #77 deleted it, because a second expression of the default
consulted only by its own tests goes on saying "somatic" whatever the app does.

The constants are still spelled ``SOFT_*`` and ``CLINICAL_*`` while the labels the user
reads say Broad and Stringent. Renaming them is held back deliberately: whether four
presets on that axis is the right *set* is still open, and the names are asserted across
two test modules, so the rename is worth doing once — after the set is settled, not
before.
"""

#: The classifications SOFT reported, as the exclude list the parameter now is.
#:
#: SOFT and CLINICAL used to state this key as *include* lists of 16 and 9 values, from
#: back when the app re-implemented the pipeline's filter and inverted it. Issue #33
#: routed the filter through the pipeline's own code, which reads the parameter the
#: pipeline's way, so an unconverted preset would say "drop every missense call". Each is
#: the complement of the include list it replaces, within the vocabulary the control can
#: offer — with one deliberate difference the inversion buys: a classification neither
#: list has heard of is now *kept* rather than silently dropped, which is what the
#: pipeline does with the 211 reference rows that carry five such values.
_SOFT_EXCLUDED_CLASSIFICATIONS = ["IGR", "RNA"]

#: CLINICAL's, likewise. It reported protein-coding and canonical splice-site calls only,
#: so everything regulatory, intronic and synonymous is excluded.
_CLINICAL_EXCLUDED_CLASSIFICATIONS = [
    "Silent",
    "Splice_Region",
    "5'UTR",
    "3'UTR",
    "5'Flank",
    "3'Flank",
    "Intron",
    "IGR",
    "RNA",
]

#: The ClinVar terms both Broad presets keep — one definition, read by each arm.
#:
#: Atomic terms only, and every term in **both** the spellings a MAF can carry it in: the
#: conflicting-classification pair, renamed by ClinVar in 2023 (issue #88), and the
#: underscore-prefixed modifiers (issue #99). See ``CLINVAR_OPTIONS`` in
#: ``config/vocabularies.py`` for why either matters. Without both spellings, whether a preset
#: reports a variant depends on which ClinVar vintage annotated the file rather than on what
#: the file says.
#:
#: **What the underscore pair costs this preset, measured:** nothing yet. Through the app's own
#: filter on four real MAFs the Broad report is unchanged, +0 rows — the depth, VAF and
#: population-frequency clauses drop those variants anyway and the pathogenic rescue admits
#: others regardless. With the ClinVar clause isolated the pair is worth 0 → 7 rows on each of
#: the two files that spell it with the underscore, and +0 on a file spelled the other way. So
#: this changes what the preset *says* it keeps, and does not currently change what it emits —
#: the same trade #88 recorded when it removed the composites at zero row cost.
#:
#: Shared rather than written out twice because it has now been edited identically twice, and
#: the second copy is where the drift would live: #88 added the rename pair to both arms and
#: #99 added the underscore pair to both, neither for a reason that distinguishes somatic from
#: germline — a MAF's spelling is a fact about the annotation, not about the arm. Should an arm
#: ever need a different ClinVar list, inline this back into both rather than parameterising
#: it; two lists that differ on purpose are two lists.
#: **What issue #103 added, and why it is not a new judgement.** This preset already kept
#: ``not_provided``, ``other`` and ``_other`` — three members of the class issue #98 named
#: **No Classification**. It did not keep the rest of that class only because the control did
#: not offer them. #103 offers them, so they arrive here on the reasoning already in this
#: list rather than on a fresh one: a preset that keeps a class keeps the class, exactly as
#: #88 and #99 each added a second *spelling* of a term already kept.
#:
#: Measured at the ClinVar clause over 176 real MAFs and 159,580 annotated rows, this is the
#: preset's one row-level change: **+855 rows**, taking the clause from 10,490 to 11,345, and
#: 854 of the 855 are ``no_classification_for_the_single_variant`` — the term #103 found the
#: old vocabulary could not select at all, in 102 of those files. Both Broad arms, since a
#: MAF's ClinVar spelling is a fact about the annotation and not about the arm.
#:
#: **And unlike #88 and #99, this one moves the report.** Both of those measured the clause and
#: then measured the whole filter and found ``+0`` — the depth, VAF and frequency clauses drop
#: those rows anyway — so each changed what a preset *said* rather than what it emitted. #103
#: does not get that for free, which is why it was measured the same way rather than assumed to
#: land the same place: through ``filters.variant_filters.apply_filters`` on twelve real MAFs
#: that carry an added term, six per arm, the Broad report goes **1,455 → 1,463 rows, +8**,
#: with eight of the twelve files moving by exactly one. 103 of the 176 measured files carry a
#: term this adds, so the effect is ordinary rather than a property of one sample.
#:
#: Small, then, but real: a variant ClinVar holds a record for and no classification of now
#: reaches a discovery report where it did not. That is the intended consequence — it is the
#: state ``not_provided`` was already admitting — and it is stated here rather than left for a
#: reader to discover from a row count that moved.
#:
#: **One of the four carries a caveat, and it belongs here rather than only where the term is
#: offered.** ``-`` is a ClinVar *term* on the institute table's authority — the variant was
#: recorded only as part of a haplotype, so there is no call for it alone — and #98 gave it a
#: class distinct from an empty cell on that basis. But some writers use a bare dash as a plain
#: missing-marker, and this preset now *selects* it, so on such a file Broad would keep rows
#: ClinVar never said anything about. No measured file does that: ``-`` occurs as a whole cell
#: in **zero** of the 176 files, which is also why it costs nothing here. The exposure is a
#: property of a MAF this corpus does not contain, so no test can hold it and it is written
#: down instead — a reader meeting an implausibly large Broad report on an unfamiliar
#: annotation pipeline should look here first.
#:
#: **What issue #114 added, and the rule it made explicit.** The four classes #103 left
#: selectable and unselected are all here now — ⚠️ Disease Risk, 🔗 Association, 🛡️ Protective,
#: and ``confers_sensitivity`` completing 💊 Drug Response. Which means this list is no longer
#: a set of classes chosen one at a time: **Broad keeps every ClinVar class except 🟢 Likely
#: Benign and 🔵 Benign**, and that is the whole of it, pinned by
#: ``test_broad_keeps_every_clinvar_class_except_the_benign_ones``. The axis is benign against
#: not-benign, not pathogenicity against other-assertion, which is what the preset's own
#: description has said since #51 — *you would rather review a variant than miss one*.
#:
#: **The institute term table's *clinical relevance for filtering* column, recorded here
#: because #114 is where it entered this repository** (the dev holds the table; #98 read it for
#: ``CLINICAL_VALUE_MAPPING`` and used this column to order
#: ``components.clinical_summary.CLINICAL_HIERARCHY``, but the calls themselves were never
#: written down): ⚠️ Disease Risk is **review-only**; 🔗 Association is **usually-exclude**;
#: 🛡️ Protective is **usually-exclude**. ``confers_sensitivity`` is in no row of it — the one
#: placement #98 made with no citable definition, a legacy value ClinVar's current
#: classification page does not list either.
#:
#: **Broad keeps all three anyway, and that is not the column being overruled.** The column
#: says what to do with a term in a *report*; Broad is the cut a reviewer makes *before* there
#: is a report, and Stringent is the report. The proof that this mapping was never a lookup is
#: already in the list above: Broad has kept ❓ No Classification whole since #103, and that
#: class is the column's **exclude** end. So a rule of the form *keep/review-only means Broad
#: selects it* would have contradicted the preset as it already stood, which is exactly why
#: #103 declined to derive one and left the judgement to be made here.
#:
#: **Measured, and unlike #88 (+77 → +0), #99 (+0) and #103 (+855 → +8) the report number is
#: not small.** The corpus is **289 byte-distinct real MAFs** swept from ``$HOME`` — name ending
#: in ``.maf``, hash-deduplicated, committed fixtures excluded, which is the definition #103's
#: 176 and #108's 297 use and is worth stating because a looser one is not comparable to them:
#: matching any name *containing* ``.maf`` pulls in ``*.maf.nopass.tsv``, the pipeline's own
#: already-filtered rejects. Of those, **158 carry a ClinVar column, an arm-identifying column
#: and survive the filter** — 97 somatic, 61 germline, 322,032 rows, 157,590 ClinVar-annotated.
#:
#: Through ``filters.variant_filters.apply_filters`` on those 158: the ClinVar clause goes
#: 11,007 → 11,418 (**+411**) and the Broad report goes 21,490 → **21,662, +172** — somatic +82,
#: germline +90, with **102 of the 158 files** moving. So this is not a preset made consistent,
#: it is a preset widened, and per class through the whole filter:
#:
#: * ⚠️ Disease Risk — clause +133, report **+2** (``HLA-A`` chr6:29913298, twice). The cheap one.
#: * 🔗 Association — clause +165, report **+84**, of which **82 are one variant**.
#: * 🛡️ Protective — clause +102, report **+86**, **all one variant**.
#: * 💊 ``confers_sensitivity`` — clause +14, report **+0**. 17 rows reach the term at all, every
#:   one of them ``FUT2``, in 5 files, and the 14 the class newly reaches are all germline; 3 of
#:   the 17 were already kept through ``other``.
#:
#: **The two variants are the finding, and they are why the clause number was not an upper bound
#: the other clauses spend here.** Association's 84 rows are ``HMGCR`` chr5:74648603 on 82 files
#: plus ``HSPA1A`` twice; Protective's 86 are all ``NOS3`` chr7:150690079 — two common
#: polymorphisms, ``Freq_1000g2015aug_all`` of **0.436** and **0.766**, surviving a 1%
#: population-frequency cut.
#:
#: Not a defect of *this* preset, and the mechanism was worth getting exactly right because
#: ``FUT2`` showed the same filter doing the opposite. At the time
#: :func:`filters.variant_filters.frequency_mask` ORed ``(value <= threshold) | value.isna()``
#: over the frequency columns **the frame has**, so *one blank kept a row whatever the populated
#: columns said*. Every file carrying ``HMGCR`` or ``NOS3`` leaves ``gnomAD_exome_AF`` empty and
#: spells ``Freq_esp6500siv2_all`` as ``.`` — two blanks, either sufficient. ``FUT2`` is common in
#: *every* panel its files fill in (gnomAD 0.38, ESP 0.49, ExAC 0.39, 1000G 0.32) and so had no
#: blank to be rescued by, which is why all 17 were dropped. It was **not** that the ``FUT2``
#: files carry gnomAD values and the others do not: ``gnomAD_exome_AF`` is absent from the file
#: altogether on 9 of those 17 rows, and an absent column was not in the mask's list at all — so
#: a column present but empty kept a variant where a column missing entirely could not.
#:
#: **Issue #122 changed that rule, and this whole widening went with it.** A blank is no longer
#: counted as zero: it still cannot sink a row, but it no longer rescues one the populated
#: columns call common, and a row no column measured is still kept. Re-measured on the same
#: corpus, these twelve terms are now worth **+0 report rows on both arms** — against the +172
#: below, reproduced exactly as a control before the change. ``HMGCR``, ``NOS3`` and Disease
#: Risk's two ``HLA-A`` rows were the entire delta, and every one of them was a common
#: polymorphism reaching the report through a blank rather than through this preset's reach.
#:
#: So the numbers below are **#114's measurement, kept as the record of what the blank rule was
#: doing** rather than as a description of the report today. What survives unchanged is #114's
#: *rule* — Broad keeps every ClinVar class except 🟢 Likely Benign and 🔵 Benign — which was
#: never the thing costing the rows. And the sentence a reviewer needed then is now the opposite
#: one: a Broad report **no longer** opens with the same recurrent polymorphism in every patient.
#:
#: **What the report then calls these rows reads like a contradiction and is not.** Measured on
#: all 172: ``Clinical_Summary`` labels them **🟢 Likely Benign (86) and 🔵 Benign (86)** — the
#: two together account for every one, and never 🛡️, 🔗 or ⚠️.
#: :func:`components.clinical_summary.generate_clinical_summary` names the
#: *strongest* class across five sources and ``CLINICAL_HIERARCHY`` ranks both benign classes
#: above every non-pathogenicity one, so a common polymorphism InterVar or RENOVO calls benign
#: is labelled benign whatever modifier ClinVar attaches. The Pathogenicity Overview's ``CV``
#: position, pinned immediately beside it, draws ClinVar's *own* class: 🛡️ for ``protective``,
#: 🔗 for ``association``, ⚠️ for ``risk_factor``, 💊 for ``confers_sensitivity``. So a row
#: reading 🛡️ at one pinned column and 🔵 Benign at the next is possible and is **not** a
#: collision: that is #100's and #104's design — the strip answers *what did each source say*,
#: the summary *what is the strongest call*. #114 recorded it here because that preset was about
#: to make it ordinary on nearly every report; #122 has made it rare again, these rows being
#: exactly the ones the blank rule was carrying. The explanation is kept rather than deleted
#: because the pairing is still reachable and still reads like a contradiction when it is met.
#:
#: Which gave the honest description of what this widening admitted, and it was **not** "risk
#: and protective findings": it was common benign polymorphisms that happen to carry a
#: non-pathogenicity ClinVar assertion. The app did not disagree with the *pipeline* about any
#: of them either — ``variantalker_naive`` is present for 162 of the 172 and agrees on all 162 —
#: so #117's measurement of those two columns is untouched by this. That description is also
#: what made #122 answerable: once the admitted rows are named as common polymorphisms rather
#: than as findings, the question is no longer whether this preset is too wide but whether the
#: frequency layer should have been keeping them at all. It should not have been.
#:
#: **Some cells carry two of these classes at once, so adding the four separately over-counts
#: adding them together — at the clause only.** The four clause deltas sum to +414 against
#: **+411** together, the difference being **3 rows**:
#: ``Benign|_association|_confers_sensitivity`` ×2 and ``Uncertain_risk_allele|protective`` ×1.
#: Through the whole filter there is **no** such gap — separately and together both give +172 —
#: because none of those three survives the depth, VAF and frequency clauses. Which is the same
#: caution as everywhere else on this list: a number measured at the clause does not carry to
#: the report, not even a correction to one.
#:
#: Completing 💊 Drug Response also **empties the one documented exception** to *a preset keeps
#: a ClinVar class whole or not at all* — ``_PARTIALLY_KEPT_BY_DESIGN`` in
#: ``tests/test_app_defaults.py`` is deleted rather than left standing empty, so there is no
#: set for the next inconsistency to be added to quietly.
#:
#: Stringent keeps ``Pathogenic``/``Likely_pathogenic`` and is untouched by all of it: nothing
#: here argues a report you would act on should carry a non-pathogenicity assertion, and #103
#: left it alone on the same ground.
_SOFT_CLINVAR_TERMS = [
    "Pathogenic",
    "Likely_pathogenic",
    "Conflicting_classifications_of_pathogenicity",
    "Conflicting_interpretations_of_pathogenicity",
    # The Drug Response class, whole (issue #114). `confers_sensitivity` is the term the
    # institute's table does not list at all, so the reason it is here is the class rather than
    # the term: keeping `drug_response` without it made whether a pharmacogenomic variant is
    # reported depend on which of two terms a submitter chose. **+0 rows** through the filter.
    "confers_sensitivity",
    "_confers_sensitivity",
    "drug_response",
    "_drug_response",
    # The Disease Risk class (issue #114), the table's *review-only* call. **Four calls in five
    # entries**: the three `*_risk_allele` tiers, and `risk_factor` — the older OMIM-derived
    # spelling of the same idea, which #98 put in one class with them on the table's authority —
    # in both the spellings #99's rule requires. **+2 rows** through the filter against +133 at
    # the clause: the frequency, depth and VAF clauses spend this one almost entirely, unlike
    # the two below — and **+0 since issue #122**, which spends the last two as well.
    "risk_factor",
    "_risk_factor",
    "Uncertain_risk_allele",
    "Likely_risk_allele",
    "Established_risk_allele",
    # The Association class (issue #114), the table's *usually-exclude* call. Measured at +84
    # rows then, of which 82 were one `HMGCR` variant at chr5:74648603 with a 1000 Genomes
    # frequency of 0.436 — **+0 since issue #122**, that variant having reached the report only
    # through a blank frequency column.
    "association",
    "_association",
    "Affects",
    # The Protective class (issue #114), likewise *usually-exclude*. Measured at +86 rows then,
    # every one of them one `NOS3` variant at chr7:150690079 with a 1000 Genomes frequency of
    # 0.766 — **+0 since issue #122**, for the same reason as Association above. Kept
    # because Broad's line is drawn at benign, and `protective` is not a benign call — ClinVar
    # groups it apart from the ACMG/AMP five and the app names it 🛡️, not 🔵.
    "protective",
    "_protective",
    "Uncertain_significance",
    # The Uncertain Significance class, whole (issue #103) — the same completion as the block
    # below, on a class this preset also already kept part of. A submitting laboratory may
    # record a VUS in its own sub-tier, and #98 reads all three as that one call because
    # ACMG/AMP permits the sub-tiers and defines no criteria for them. Keeping the plain term
    # and not these would make whether a variant is reported depend on how its submitter chose
    # to write the same call — which is #99's defect exactly, in a new place. **+0 rows**: no
    # sub-tier occurs anywhere in the 159,580-row corpus, so this is a preset made consistent
    # rather than widened, and the reason it is here at all is that reach does not decide.
    "VUS-high",
    "VUS-mid",
    "VUS-low",
    # The No Classification class, whole (issue #103). The first three were here already.
    "not_provided",
    "other",
    "_other",
    "conflicting_data_from_submitters",
    "-",
    "no_classification_for_the_single_variant",
    "no_classifications_from_unflagged_records",
]

# Soft filter parameters
SOFT_SOMATIC_PARAMS = {
    "sample_type": "somatic",
    "min_depth": 50,
    "vaf_threshold": 0.05,
    # The one gene-list parameter, unprefixed and empty: no gene filtering.
    # The panel dropdown is UI state that resolves to this list, never a parameter of
    # its own -- see page_modules/parameter_config.gene_panel_state_key (issue #34).
    "filter_genes": [],
    "filter_variant_classification": list(_SOFT_EXCLUDED_CLASSIFICATIONS),
    "filter_cancervar": ["Tier_I_strong", "Tier_II_potential", "Tier_III_Uncertain"],
    "filter_civic": ["A", "B", "C", "D"],
    "filter_clinvar": list(_SOFT_CLINVAR_TERMS),
    "filter_escat": ["IA", "IB", "IC", "IIA", "IIB", "IIIA", "IIIB"],
    "max_freq_population": 0.01,  # 5% for soft somatic (more lenient for discovery)
    "keep_pathogenic": True,  # Auto-retain pathogenic variants by default
}

SOFT_GERMLINE_PARAMS = {
    "sample_type": "germline",
    "min_depth": 50,
    "vaf_threshold_germline": 0.2,
    # The one gene-list parameter, unprefixed and empty: no gene filtering.
    # The panel dropdown is UI state that resolves to this list, never a parameter of
    # its own -- see page_modules/parameter_config.gene_panel_state_key (issue #34).
    "filter_genes": [],
    "filter_variant_classification": list(_SOFT_EXCLUDED_CLASSIFICATIONS),
    "filter_intervar": ["Pathogenic", "Likely pathogenic", "Uncertain significance"],
    "filter_renovo": ["LP Pathogenic", "IP Pathogenic", "HP Pathogenic"],
    "filter_clinvar": list(_SOFT_CLINVAR_TERMS),
    "filter_escat": ["IA", "IB", "IC", "IIA", "IIB", "IIIA", "IIIB"],
    "max_freq_population": 0.01,  # 5% for soft germline (more lenient for discovery)
    "keep_pathogenic": True,  # Auto-retain pathogenic variants by default
}

# Clinical somatic parameters
CLINICAL_SOMATIC_PARAMS = {
    "sample_type": "somatic",
    "min_depth": 50,
    "vaf_threshold": 0.05,
    # The one gene-list parameter, unprefixed and empty: no gene filtering.
    # The panel dropdown is UI state that resolves to this list, never a parameter of
    # its own -- see page_modules/parameter_config.gene_panel_state_key (issue #34).
    "filter_genes": [],
    "filter_variant_classification": list(_CLINICAL_EXCLUDED_CLASSIFICATIONS),
    "filter_cancervar": ["Tier_I_strong", "Tier_II_potential"],
    "filter_civic": ["A", "B"],
    "filter_clinvar": ["Pathogenic", "Likely_pathogenic"],
    "filter_escat": ["IA", "IB", "IC"],
    "max_freq_population": 0.005,
    "keep_pathogenic": True,  # Auto-retain pathogenic variants by default
}

# Clinical germline parameters
CLINICAL_GERMLINE_PARAMS = {
    "sample_type": "germline",
    "min_depth": 50,
    "vaf_threshold_germline": 0.2,
    # The one gene-list parameter, unprefixed and empty: no gene filtering.
    # The panel dropdown is UI state that resolves to this list, never a parameter of
    # its own -- see page_modules/parameter_config.gene_panel_state_key (issue #34).
    "filter_genes": [],
    "filter_variant_classification": list(_CLINICAL_EXCLUDED_CLASSIFICATIONS),
    "filter_intervar": ["Pathogenic", "Likely pathogenic"],
    "filter_renovo": ["HP Pathogenic"],
    "filter_clinvar": ["Pathogenic", "Likely_pathogenic"],
    "filter_escat": ["IA", "IB", "IC"],
    "max_freq_population": 0.005,
    "keep_pathogenic": True,  # Auto-retain pathogenic variants by default
}
