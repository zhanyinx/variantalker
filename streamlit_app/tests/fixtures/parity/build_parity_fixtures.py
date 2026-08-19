"""Build the parity fixture set. Constructed end to end; nothing mounted, nothing real.

WHAT REPLACED WHAT, AND WHY IT WAS ALL SEVEN FILES. The generator this replaces cut rows
out of the 50-sample paired GERSOM run and pseudonymised the barcodes. Issue #233 measured
what that left behind, and the answer reframed the problem: **all seven** MAFs carried real
calls, not the two the tickets named.

* ``somatic_gnomad_genome.maf`` — 6 rows byte-identical to reference rows on all 427
  shared columns.
* ``somatic_dot_numeric.maf`` — 5 reference rows with one cell overwritten by ``.``.
* ``somatic_synthetic.maf``, ``germline_synthetic.maf``, ``somatic_civic.maf`` — each row a
  real reference row with two to four cells overwritten, so 423-425 of 427 columns were
  still the real row's.

The ``*_synthetic`` files rewrote ``Start_Position`` to a round number, which reads as
de-identification and is not: the same row still carried ``g.chr17:37871547C>A`` in
``Genome_Change`` and ``p.A386D`` in ``Protein_Change``, so the real coordinate survived
verbatim in an annotation column and the variant was fully recoverable. "Synthetic" in the
old set meant *one field was edited to flip a verdict*, not *no patient data*.

So this builds all of it. The only things taken from the old set are in
``column_profiles.json``: column **names**, the leading comment-line **count**, and per
column which blank token the annotator uses for it and whether its values are numeric.
Those are Funcotator/ANNOVAR software artifacts — the same class of metadata as the column
names themselves — and no cell value is among them. Every value below is written here.

THE CONSTRUCTION RULE, which the fixture README paid for the hard way:

    a witness cell must pin the row's whole path to the verdict, not just the field
    under test.

So every witness row holds every *other* divergence neutral — depth clear on both
``t_alt+t_ref`` and ``DP``, VAF off every swept boundary, classification inside the app's
vocabulary and outside the exclude list, ClinVar agreeing under split and exact matching —
and the rows that witness a *rescue* divergence are given a sub-threshold VAF, because the
rescue is OR'd in after the VAF conjunction and is otherwise invisible: at contract
defaults a ``Tier_I_strong`` row passes the guideline block anyway, so removing it from the
rescue changes nothing. That is the trap divergence #6 fell into, in a different clause.

That rule has a cost this file has to pay back deliberately, and ``overlap_rows()`` is
where: a set in which every row is excluded for exactly one reason has no row excluded for
two, and the run recap's "these counts overlap" caption then describes something the
fixture does not do. Witness rows stay single-reason; the overlap rows are separate and
witness nothing.

Usage::

    python streamlit_app/tests/fixtures/parity/build_parity_fixtures.py --out DIR
    python streamlit_app/tests/fixtures/parity/build_parity_fixtures.py --check

``--check`` rebuilds into a temporary directory and compares against the committed files,
which is what makes "the generator reproduces the fixture set" a checked claim rather than
a remembered one. There is no RNG anywhere and no clock is read, so two builds are byte
identical.

The per-row expected verdicts in ``MANIFEST.json`` are **not** written here. They are
measured by ``record_manifest.py`` in a second step, by running the pipeline — a verdict
this file declared would be an assertion about the pipeline rather than an observation of
it.
"""

from __future__ import annotations

import argparse
import filecmp
import json
import tempfile
import zlib
from pathlib import Path

HERE = Path(__file__).resolve().parent
PROFILES = json.loads((HERE / "column_profiles.json").read_text())
COLUMN_PROFILES = PROFILES["columns"]
FIXTURES = PROFILES["fixtures"]

#: The reference set's leading comment-line count, reproduced. ``read_maf`` counts ``#``
#: lines rather than assuming a fixed offset, and ``writeheader`` copies the whole block
#: into the pipeline's output, so the block's *size* is part of what the fixtures exercise.
COMMENT_LINES = 111

# The parts of the contract this file actually builds rows against, mirrored so the
# generator is readable on its own. Only the ones something below *uses*: an unreferenced
# mirror is a copy of the contract with nothing holding it in step, and the parity suite
# reads tests/parity/contract.py rather than trusting any of this. Where a docstring needs
# to name a value it does not build with — `min_depth`, `filter_clinvar` — it names the
# contract's key rather than a constant kept alive to be quoted.
VAF_SOMATIC = 0.01
VAF_GERMLINE = 0.2
EXCLUDED_VC = ["Silent", "IGR", "RNA"]
CANCERVAR_KEEP = ["Tier_II_potential", "Tier_I_strong"]
ESCAT_KEEP = ["IA", "IB", "IC", "IIA", "IIB"]
RENOVO_KEEP = ["LP Pathogenic", "IP Pathogenic", "HP Pathogenic"]

#: Classifications the app's hardcoded vocabulary does **not** contain, and which are not
#: in the pipeline's exclude list either — the only shape that witnesses divergence #1.
#: Real Funcotator classification names; none is patient data.
VC_OUTSIDE_APP_VOCABULARY = [
    "Start_Codon_SNP",
    "Start_Codon_Del",
    "Start_Codon_Ins",
    "De_novo_Start_InFrame",
    "De_novo_Start_OutOfFrame",
    "Stop_Codon_Del",
]

#: ClinVar strings that match the contract's ``filter_clinvar`` **only after the separator
#: split**, and are absent from ``CLINVAR_PATHO`` as literals too, so one value witnesses
#: both legs of divergence #4. Every separator the vendored splitter handles is represented.
CLNSIG_SPLIT_ONLY = [
    "Pathogenic|not_provided",
    "Likely_pathogenic;other",
    "Pathogenic/Likely_pathogenic|not_provided",
    "not_provided,Pathogenic",
]

#: In ``CLINVAR_PATHO`` but matching no term of the contract's ``filter_clinvar`` even
#: under the split, so the row can only ever pass through the pathogenic rescue. Witnesses
#: the rescue's ClinVar
#: leg — and it is the only ClinVar value that does, which is why these rows are kept even
#: though they plant a two-code-point ``⚠️`` on the overview strip. Issue #241 fixed the
#: legend guard to count glyphs rather than dropping the rows to satisfy it.
CLNSIG_RESCUE_ONLY = ["Established_risk_allele", "Likely_risk_allele"]

#: The gene lists, **carried over verbatim** rather than reinvented.
#:
#: A list of twelve HGNC symbols is not patient data. The symbols were *chosen* by looking
#: at a real run, but "APC is in this panel" says nothing about any individual — unlike a
#: variant call, which is quasi-identifying on its own. Keeping them verbatim also keeps
#: ``genes_somatic_mixed_case.txt``'s exact mangling, which is irregular rather than a rule
#: applied to the symbols, and which ``test_gene_lists`` asserts on by name.
GENES_SOMATIC = [
    "APC", "CCDC107", "CUL3", "DROSHA", "ERBB3", "HLA-A",
    "KMT2C", "MAP3K4", "MUC6", "PDS5B", "PMS2", "SCAF4",
]
GENES_GERMLINE = [
    "APOB", "AXIN2", "CDA", "DPYD", "EZH2", "FH",
    "KLK3", "MSH2", "NOTCH3", "PIK3CA", "ROS1", "SRSF2",
]
GENES_SOMATIC_MIXED_CASE = [
    "Apc", "ccdc107", "Cul3", "drosha", "Erbb3", "hla-a",
    "Kmt2C", "map3k4", "Muc6", "pds5b", "Pms2", "scaf4",
]

#: Symbols the gene lists carry, and symbols they do not. Sliced from the lists above
#: rather than written out again: a second copy would go on looking in-list after the panel
#: changed under it, and "this row survives the gene filter" would quietly stop being true.
#: The default template gene is in-list on both arms, matching the old set's balance —
#: there, a gene restriction moved 20 PASS to 18 rather than to nearly zero, and a fixture
#: where the list excludes almost everything would make the gene cases test something much
#: coarser.
IN_LIST_SOMATIC = GENES_SOMATIC[:6]
IN_LIST_GERMLINE = GENES_GERMLINE[:6]
OFF_LIST = ["MYC", "JAK2", "NRAS", "IDH1"]

#: The three symbols ``test_param_migration``'s legacy somatic parameter file names, and
#: the two its germline one does. They are in no gene list this fixture set ships, which is
#: the point: that module migrates a *legacy* file and asks whether the migrated
#: restriction reaches the filter at all.
LEGACY_GENES_SOMATIC = ["TP53", "KRAS", "EGFR"]
LEGACY_GENES_GERMLINE = ["BRCA1", "BRCA2"]

#: RENOVO's six classes and a score inside each one's measured ``RENOVO_pls`` interval.
#: The intervals are ``tests/test_clinical_summary.py``'s ``RENOVO_TIERS``, measured over
#: 211,744 real rows; a fixture whose score sits outside its class's band fails that
#: module. All six appear in the germline fixture on purpose — ``IP Pathogenic`` occurred on
#: **no row** of the old ``germline_reference.maf``, a gap that module documents and works
#: around, and a constructed set can simply close.
RENOVO_SCORES = {
    "HP Pathogenic": 0.95,
    "IP Pathogenic": 0.83,
    "LP Pathogenic": 0.65,
    "LP Benign": 0.35,
    "IP Benign": 0.12,
    "HP Benign": 0.005,
}

#: ``am_class``'s two boundaries on ``am_pathogenicity``, from ``test_alphamissense``:
#: below 0.34 benign, above 0.564 pathogenic, ambiguous between. The class is *derived*
#: here rather than set beside the score, because that module asserts the function.
AM_BENIGN_BELOW = 0.34
AM_PATHOGENIC_ABOVE = 0.564

#: Every population-frequency column the app reads
#: (``filters.variant_filters.FREQUENCY_COLUMNS``), less the two ``gnomAD_genome_*`` names,
#: which are a column-*presence* probe and must stay absent from the other fixtures.
FREQUENCY_COLUMNS = (
    "gnomAD_exome_AF",
    "gnomAD_exome_AF_raw",
    "Freq_ExAC_ALL",
    "Freq_esp6500siv2_all",
    "Freq_1000g2015aug_all",
)


def _frequency(value) -> dict:
    """Set every frequency column to one value, so a row is unambiguously common or rare."""
    return {column: value for column in FREQUENCY_COLUMNS}


def _no_frequency() -> dict:
    """A row with no population frequency at all — the all-absent branch."""
    return {column: "" for column in FREQUENCY_COLUMNS}


# ---------------------------------------------------------------------------
# The comment block
# ---------------------------------------------------------------------------


def comment_block(arm: str, fixture: str) -> list[str]:
    """Exactly ``COMMENT_LINES`` well-formed VCF header lines.

    The pipeline's ``read_maf`` counts every line starting with ``#`` anywhere in the file
    and passes that count as ``header=``, so the block must sit at the top and nothing
    below may start with ``#``. ``writeheader`` copies it into the pipeline's output, so it
    has to be well-formed as well as counted.

    The count is reproduced rather than shortened because the block's size is part of what
    these fixtures exercise: a fixture with a three-line header would leave "counts ``#``
    lines rather than assuming a fixed offset" an untested claim.

    The filler lines are labelled as filler rather than dressed up as reference metadata.
    Inventing 80 plausible-looking contig and command lines would make the block *look*
    like a real run's, which is the one property this set exists not to have.
    """
    lines = [
        "##fileformat=VCFv4.2",
        "##source=build_parity_fixtures.py",
        f"##fixture={fixture}",
        f"##arm={arm}",
        "##provenance=constructed by streamlit_app/tests/fixtures/parity/"
        "build_parity_fixtures.py -- no real variant calls, no real sample identifiers",
        '##INFO=<ID=CONSTRUCTED,Number=0,Type=Flag,Description="Constructed fixture row">',
        '##FILTER=<ID=PASS,Description="All filters passed">',
    ]
    for n in list(range(1, 23)) + ["X", "Y", "M"]:
        lines.append(f"##contig=<ID=chr{n},length=100000000>")
    index = 0
    while len(lines) < COMMENT_LINES:
        index += 1
        lines.append(
            f"##fixture.headerFiller={index:02d},Description="
            '"padding to the reference set\'s comment-line count; carries no information"'
        )
    return lines[:COMMENT_LINES]


# ---------------------------------------------------------------------------
# Row construction
# ---------------------------------------------------------------------------


def neutral(arm: str) -> dict:
    """A row that fails every guideline block and every rescue, and diverges on nothing.

    Read this as the definition of "held neutral": depth clears the contract's ``min_depth`` on the sum
    *and* on ``DP``; VAF is above every swept threshold (0.01, 0.05, 0.2) and on none of
    them; the classification is inside the app's vocabulary and outside the exclude list;
    the ClinVar string matches under neither exact nor split comparison.
    """
    row = {
        "Hugo_Symbol": "APC" if arm == "somatic" else "APOB",
        "Variant_Classification": "Missense_Mutation",
        "Variant_Type": "SNP",
        "Reference_Allele": "C",
        "Tumor_Seq_Allele1": "C",
        "Tumor_Seq_Allele2": "A",
        "t_alt_count": 40,
        "t_ref_count": 60,
        "DP": 100,
        "n_alt_count": 0,
        "n_ref_count": 80,
        "tumor_f": 0.5,
        "ClinVar_VCF_CLNSIG": "Benign",
        "ESCAT": ".",
        "ESCAT_TISSUE": ".",
        "ESCAT_CANCER": ".",
        "am_pathogenicity": 0.1,
        "Annotation_Transcript": "ENST00000000001.1",
        "cDNA_Change": "c.100C>A",
        "Codon_Change": "c.(100-102)CCA>ACA",
        "Protein_Change": "p.P34T",
        "Transcript_Exon": 3,
        "AAChange.refGene": "SYNTH:ENST00000000001.1:exon3:c.100C>A:p.P34T",
        "variantalker_naive": "constructed",
        "Otherinfo": ".",
        "tumor_tissue": "constructed_tissue",
        "cosmic": ".",
        "Freq_ExAC_ALL": 0.001,
        "Freq_esp6500siv2_all": 0.001,
        "Freq_1000g2015aug_all": 0.001,
        "gnomAD_exome_AF": 0.001,
        "gnomAD_exome_AF_raw": 0.001,
        "Tumor_Sample_Barcode": "SYNTH_01",
        "Matched_Norm_Sample_Barcode": "SYNTH_01_N",
    }
    if arm == "somatic":
        row["CancerVar"] = "Tier_III_Uncertain"
    else:
        # ``__UNKNOWN__`` verbatim on the tumour barcode: issue #15's finding is about this
        # value, and it is a pipeline placeholder rather than an identifier.
        #
        # The matched-normal barcode must NOT be ``__UNKNOWN__`` as well, and finding that
        # out cost a build: ``filter_variants.py:452`` narrows the NOPASS frame with
        # ``~(out == "__UNKNOWN__").all()``, dropping every column whose *entire* value set
        # is that placeholder — and ``Matched_Norm_Sample_Barcode`` is in the germline
        # ``keep``, so the run dies on ``KeyError: ['Matched_Norm_Sample_Barcode'] not in
        # index`` before producing a verdict. A fixture that pseudonymises a KEEP column to
        # one constant crashes the pipeline; the old set avoided it by accident, because
        # only the column germline ``keep`` removes was uniform.
        row["Tumor_Sample_Barcode"] = "__UNKNOWN__"
        row["Matched_Norm_Sample_Barcode"] = "SYNTH_NORMAL_01"
        row["InterVar"] = "Benign"
        row["RENOVO_Class"] = "HB Benign"
        row["RENOVO_pls"] = 0.05
        row["CLINSIG"] = "Benign"
    return row


def somatic_rows() -> list[tuple[str, dict]]:
    """The somatic reference-replacement rows, one ``(cell, overrides)`` pair each."""
    rows: list[tuple[str, dict]] = []

    def add(name, **over):
        rows.append((name, over))

    # -- Guideline-pass carriers. ESCAT-driven on purpose: under the #11 contract the
    # other somatic guideline triggers are a subset of the rescue, so an ESCAT row is the
    # only kind whose verdict a depth/VAF/gene sweep can move without --skip_pathogenic.
    add("escat_pass_depth60", ESCAT="IA", t_alt_count=25, t_ref_count=35, DP=60)
    add("escat_pass_depth600", ESCAT="IB", t_alt_count=200, t_ref_count=400, DP=600,
        Hugo_Symbol="CCDC107")
    add("escat_pass_depth700", ESCAT="IIA", t_alt_count=300, t_ref_count=400, DP=700,
        Hugo_Symbol="CUL3")
    add("escat_pass_depth520", ESCAT="IIB", t_alt_count=220, t_ref_count=300, DP=520,
        Hugo_Symbol="DROSHA")
    # VAF boundaries: one row *at* each swept threshold, which is what witnesses the
    # ``>`` / ``>=`` divergence, plus rows between them so each sweep step moves.
    add("vaf_exactly_001", ESCAT="IA", tumor_f=0.01, Hugo_Symbol="ERBB3")
    add("vaf_exactly_005", ESCAT="IA", tumor_f=0.05, Hugo_Symbol="HLA-A")
    add("vaf_exactly_02", ESCAT="IB", tumor_f=0.2, Hugo_Symbol="KMT2C")
    add("vaf_between_005_02", ESCAT="IIA", tumor_f=0.12, Hugo_Symbol="CCDC107")
    add("vaf_above_02", ESCAT="IIB", tumor_f=0.35, Hugo_Symbol="CUL3")
    # Gene-filter carriers: off-list rows the gene cases must drop, in-list rows they keep.
    add("escat_pass_off_list", ESCAT="IA", Hugo_Symbol="MYC")
    add("escat_pass_off_list_2", ESCAT="IB", Hugo_Symbol="JAK2")
    add("escat_pass_in_list", ESCAT="IIA", Hugo_Symbol="MAP3K4")
    add("escat_pass_in_list_2", ESCAT="IA", Hugo_Symbol="MUC6")

    # -- #1: classification outside the app's vocabulary but not in the exclude list.
    for i, vc in enumerate(VC_OUTSIDE_APP_VOCABULARY):
        add(
            f"vc_outside_app_vocabulary_{i + 1}",
            Variant_Classification=vc,
            ESCAT=ESCAT_KEEP[i % len(ESCAT_KEEP)],
            Hugo_Symbol=(IN_LIST_SOMATIC + OFF_LIST)[i % 10],
        )
    # -- #1 control: the exclude list itself, guideline-passing so the drop is the filter's.
    for vc in EXCLUDED_VC:
        add(f"vc_excluded_{vc.lower()}", Variant_Classification=vc, ESCAT="IA")

    # -- #2: DP and the sum disagree across the contract's min_depth, both directions.
    add("dp_high_sum_low", ESCAT="IA", t_alt_count=10, t_ref_count=20, DP=100)
    add("dp_low_sum_high", ESCAT="IB", t_alt_count=40, t_ref_count=40, DP=30)
    add("dp_high_sum_low_2", ESCAT="IIA", t_alt_count=15, t_ref_count=25, DP=500)
    add("dp_low_sum_high_2", ESCAT="IIB", t_alt_count=60, t_ref_count=60, DP=10)

    # -- #4: ClinVar matching only after the split, on both the guideline and rescue legs.
    for i, value in enumerate(CLNSIG_SPLIT_ONLY):
        add(f"clnsig_split_only_{i + 1}", ClinVar_VCF_CLNSIG=value)
    add("clnsig_exact_pathogenic", ClinVar_VCF_CLNSIG="Pathogenic")
    add("clnsig_exact_likely", ClinVar_VCF_CLNSIG="Likely_pathogenic")

    # -- #7, rescue leg: in CLINVAR_PATHO, matching no keep term even after the split, and
    # sub-threshold VAF so the guideline conjunction cannot mask the rescue.
    for i, value in enumerate(CLNSIG_RESCUE_ONLY):
        add(f"patho_rescue_only_{i + 1}", ClinVar_VCF_CLNSIG=value, tumor_f=0.005)
    # -- #7, CancerVar rescue: same trick, so --skip_pathogenic actually moves them.
    add("cancervar_tier1_rescue_only", CancerVar="Tier_I_strong", tumor_f=0.005)
    add("cancervar_tier2_rescue_only", CancerVar="Tier_II_potential", tumor_f=0.005)
    add("cancervar_tier1_guideline", CancerVar="Tier_I_strong", Hugo_Symbol="DROSHA")
    # -- #7, the dead clause: bare ESCAT tiers the app's retain clause used to admit.
    for value in ["I", "II", "III"]:
        add(f"escat_bare_{value}", ESCAT=value)
    add("escat_outside_keep", ESCAT="IIIA")

    # -- #10 and the app-only frequency filter. The app's population-frequency filter is
    # its own extra, so no parity case exercises it — but ``test_filter_app_extras``,
    # ``test_filter_notes`` and ``test_stale_banners`` all refuse to run without rows it
    # can actually move, and say so rather than passing vacuously. What they need:
    #   * rows above the 0.01 preset that pass on criteria alone, so the filter removes
    #     something and the report fires;
    #   * a *common pathogenic* row — high frequency and ClinVar-pathogenic — passing on
    #     criteria, so the exemption has something to exempt;
    #   * rows with every frequency column missing, so the all-absent branch is taken.
    add("gnomad_absent", ESCAT="IA", **_no_frequency())
    add("gnomad_absent_2", ESCAT="IB", Hugo_Symbol="CCDC107", **_no_frequency())
    add("gnomad_rare", ESCAT="IB", **_frequency(0.00004))
    add("gnomad_high_patho_exempt", ClinVar_VCF_CLNSIG="Pathogenic", ESCAT="IA",
        **_frequency(0.35))
    add("gnomad_high_patho_exempt_2", ClinVar_VCF_CLNSIG="Likely_pathogenic",
        ESCAT="IIA", Hugo_Symbol="CUL3", **_frequency(0.28))
    # See the germline twin: the exemption needs a split-only witness, not only exact
    # CLINVAR_PATHO matches.
    add("gnomad_high_patho_exempt_split_only",
        ClinVar_VCF_CLNSIG="Pathogenic/Pathogenic,_low_penetrance|other",
        Hugo_Symbol="PDS5B", **_frequency(0.41))
    # A common row that reaches the report on the *rescue* path and is **not** exempt, so
    # the frequency mask is shown reaching the pathogenic union rather than only the
    # criteria path. CancerVar rescues it; the exemption reads ClinVar, which is Benign
    # here, so nothing saves it. Sub-threshold VAF keeps it off the criteria path.
    add("gnomad_high_rescue_dropped", CancerVar="Tier_I_strong", tumor_f=0.005,
        Hugo_Symbol="SCAF4", **_frequency(0.55))
    for i, value in enumerate([0.02, 0.11, 0.4, 0.62, 0.85]):
        add(f"gnomad_common_not_exempt_{i + 1}", ESCAT=ESCAT_KEEP[i % 5],
            Hugo_Symbol=(IN_LIST_SOMATIC + OFF_LIST)[i], **_frequency(value))
    # ``.`` in a frequency column: the pipeline never reads these columns, but the app
    # coerces them, and every fixture is required to exercise that coercion.
    add("frequency_dot", ESCAT="IA", Freq_esp6500siv2_all=".",
        Freq_1000g2015aug_all=".")

    # One row per gene in ``genes_somatic.txt``, so a gene-list absence warning that names
    # a listed gene is a real finding rather than a fixture gap.
    for gene in GENES_SOMATIC:
        add(f"gene_row_{gene}", Hugo_Symbol=gene, ESCAT="IB")

    rows.extend(review_status_rows("somatic"))
    rows.extend(cancervar_evidence_rows())
    rows.extend(legacy_gene_rows("somatic"))
    rows.extend(overlap_rows("somatic"))

    add("control_guideline_fail")
    return rows


def germline_rows() -> list[tuple[str, dict]]:
    rows: list[tuple[str, dict]] = []

    def add(name, **over):
        rows.append((name, over))

    # -- Guideline-pass carriers, RENOVO-driven: ``RENOVO_Class`` is the one germline
    # guideline trigger absent from the rescue, so these are the rows a sweep can move.
    add("renovo_lp_depth60", RENOVO_Class="LP Pathogenic", t_alt_count=25,
        t_ref_count=35, DP=60)
    add("renovo_ip_depth600", RENOVO_Class="IP Pathogenic", t_alt_count=200,
        t_ref_count=400, DP=600, Hugo_Symbol="AXIN2")
    add("renovo_hp_depth700", RENOVO_Class="HP Pathogenic", t_alt_count=300,
        t_ref_count=400, DP=700, Hugo_Symbol="CDA")
    add("renovo_lp_depth520", RENOVO_Class="LP Pathogenic", t_alt_count=220,
        t_ref_count=300, DP=520, Hugo_Symbol="DPYD")

    # -- #3: six rows at exactly the 0.2 germline threshold, matching the distributional
    # witness the old set supplied. Constructed rather than sampled, which is the point.
    for i in range(6):
        add(
            f"vaf_exactly_02_{i + 1}",
            RENOVO_Class=RENOVO_KEEP[i % 3],
            tumor_f=0.2,
            Hugo_Symbol=(IN_LIST_GERMLINE + OFF_LIST)[i % 10],
        )
    add("vaf_between_005_02", RENOVO_Class="LP Pathogenic", tumor_f=0.12)
    add("vaf_exactly_005", RENOVO_Class="IP Pathogenic", tumor_f=0.05)
    add("vaf_above_02", RENOVO_Class="HP Pathogenic", tumor_f=0.45)

    add("renovo_pass_off_list", RENOVO_Class="LP Pathogenic", Hugo_Symbol="MYC")
    add("renovo_pass_off_list_2", RENOVO_Class="IP Pathogenic", Hugo_Symbol="NRAS")
    add("renovo_pass_in_list", RENOVO_Class="HP Pathogenic", Hugo_Symbol="MSH2")
    add("renovo_pass_in_list_2", RENOVO_Class="LP Pathogenic", Hugo_Symbol="EZH2")

    # -- #1, germline.
    for i, vc in enumerate(VC_OUTSIDE_APP_VOCABULARY):
        add(
            f"vc_outside_app_vocabulary_{i + 1}",
            Variant_Classification=vc,
            RENOVO_Class=RENOVO_KEEP[i % 3],
            Hugo_Symbol=(IN_LIST_GERMLINE + OFF_LIST)[i % 10],
        )
    for vc in EXCLUDED_VC:
        add(f"vc_excluded_{vc.lower()}", Variant_Classification=vc,
            RENOVO_Class="LP Pathogenic")

    # -- #2, germline.
    add("dp_high_sum_low", RENOVO_Class="LP Pathogenic", t_alt_count=10,
        t_ref_count=20, DP=100)
    add("dp_low_sum_high", RENOVO_Class="IP Pathogenic", t_alt_count=40,
        t_ref_count=40, DP=30)
    add("dp_high_sum_low_2", RENOVO_Class="HP Pathogenic", t_alt_count=15,
        t_ref_count=25, DP=500)
    add("dp_low_sum_high_2", RENOVO_Class="LP Pathogenic", t_alt_count=60,
        t_ref_count=60, DP=10)

    # -- #4, germline: guideline leg (VAF clear) and rescue leg (VAF sub-threshold).
    for i, value in enumerate(CLNSIG_SPLIT_ONLY):
        add(f"clnsig_split_only_{i + 1}", ClinVar_VCF_CLNSIG=value)
    add("clnsig_split_only_rescue", ClinVar_VCF_CLNSIG="Pathogenic|not_provided",
        tumor_f=0.05)
    add("clnsig_exact_pathogenic", ClinVar_VCF_CLNSIG="Pathogenic")

    # -- #6: the largest divergence on real data. Everything neutral, Benign on InterVar,
    # ClinVar *and* RENOVO, depth clear on both columns, VAF off every boundary, so the row
    # can diverge on the unmirrored ESCAT clause and on nothing else. Six of them, as the
    # old set had.
    for i in range(6):
        add(
            f"germline_escat_only_{i + 1}",
            ESCAT=ESCAT_KEEP[i % len(ESCAT_KEEP)],
            ESCAT_TISSUE="constructed_tissue",
            ESCAT_CANCER="constructed_cancer",
            Hugo_Symbol=(IN_LIST_GERMLINE + OFF_LIST)[i % 10],
        )
    # ESCAT populated but outside the keep list: exercises the column without claiming to
    # witness the divergence. The README's ``#6 shape`` distinction, kept explicit.
    add("germline_escat_outside_keep", ESCAT="IIIA")

    # -- #8: the germline rescue's key. Sub-threshold VAF on all of them, so the rescue is
    # the only thing that can decide the row.
    add("rescue_only_intervar_patho", InterVar="Pathogenic", tumor_f=0.05)
    add("rescue_only_intervar_likely", InterVar="Likely pathogenic", tumor_f=0.05)
    add("rescue_only_intervar_patho_2", InterVar="Pathogenic", tumor_f=0.03,
        Hugo_Symbol="FH")
    # No ``RENOVO_pls`` here: the score follows the class (see ``materialise``). Setting it
    # beside the class is exactly how the first build broke ``test_clinical_summary`` — an
    # ``LP Pathogenic`` row carrying 0.88, which is ``IP Pathogenic``'s band.
    add("renovo_patho_intervar_benign", RENOVO_Class="HP Pathogenic", tumor_f=0.05)
    add("renovo_patho_intervar_benign_2", RENOVO_Class="LP Pathogenic", tumor_f=0.05,
        Hugo_Symbol="KLK3")
    add("rescue_only_clnsig_patho", ClinVar_VCF_CLNSIG="Pathogenic", tumor_f=0.05)
    # In CLINVAR_PATHO, no keep term under the split: the rescue's ClinVar leg alone.
    for i, value in enumerate(CLNSIG_RESCUE_ONLY):
        add(f"patho_rescue_only_{i + 1}", ClinVar_VCF_CLNSIG=value, tumor_f=0.005)

    add("gnomad_absent", RENOVO_Class="LP Pathogenic", **_no_frequency())
    add("gnomad_absent_2", RENOVO_Class="IP Pathogenic", Hugo_Symbol="FH",
        **_no_frequency())
    add("gnomad_rare", RENOVO_Class="IP Pathogenic", **_frequency(0.00004))
    add("gnomad_high_patho_exempt", ClinVar_VCF_CLNSIG="Pathogenic",
        RENOVO_Class="LP Pathogenic", **_frequency(0.35))
    add("gnomad_high_patho_exempt_2", ClinVar_VCF_CLNSIG="Likely_pathogenic",
        RENOVO_Class="HP Pathogenic", Hugo_Symbol="CDA", **_frequency(0.28))
    # The exemption must be witnessed by a *split-only* ClinVar string as well as by exact
    # matches, or ``test_filter_app_extras`` reports that the fixture has stopped covering
    # the low-penetrance case.
    add("gnomad_high_patho_exempt_split_only",
        ClinVar_VCF_CLNSIG="Pathogenic/Pathogenic,_low_penetrance|other",
        Hugo_Symbol="NOTCH3", **_frequency(0.41))
    # The germline twin of the rescue-path drop: InterVar rescues it, ClinVar is Benign so
    # the exemption does not, and the sub-threshold VAF keeps it off the criteria path.
    add("gnomad_high_rescue_dropped", InterVar="Pathogenic", tumor_f=0.05,
        Hugo_Symbol="SRSF2", **_frequency(0.55))
    for i, value in enumerate([0.02, 0.11, 0.4, 0.62, 0.85]):
        add(f"gnomad_common_not_exempt_{i + 1}", RENOVO_Class=RENOVO_KEEP[i % 3],
            Hugo_Symbol=(IN_LIST_GERMLINE + OFF_LIST)[i], **_frequency(value))
    add("frequency_dot", RENOVO_Class="LP Pathogenic", Freq_esp6500siv2_all=".",
        Freq_1000g2015aug_all=".")

    # All six RENOVO classes, each with a score inside its own measured band. The old
    # fixture carried five: ``IP Pathogenic`` appeared on no row of it.
    for renovo_class in RENOVO_SCORES:
        add(f"renovo_class_{renovo_class.replace(' ', '_')}", RENOVO_Class=renovo_class)

    for gene in GENES_GERMLINE:
        add(f"gene_row_{gene}", Hugo_Symbol=gene, RENOVO_Class="IP Pathogenic")

    # The two AlphaMissense columns are different annotations, and the app reads exactly
    # one of them. See ``alphamissense_disagreement`` for why they are set apart.
    rows.extend(alphamissense_rows())
    rows.extend(review_status_rows("germline"))
    rows.extend(legacy_gene_rows("germline"))
    rows.extend(overlap_rows("germline"))

    add("control_guideline_fail")
    return rows


#: ClinVar's review-status vocabulary, weakest assertion first. This is a **separate,
#: ordered** scale from ``CLNSIG`` — a call and how well reviewed that call is — and two of
#: its spellings share a level, which is why the terms are written out rather than
#: generated. All are ClinVar's own published terms; none is patient data.
CLINVAR_REVIEW_STATUS = [
    "no_assertion_criteria_provided",
    "criteria_provided,_single_submitter",
    "criteria_provided,_conflicting_classifications",
    "criteria_provided,_multiple_submitters,_no_conflicts",
    "reviewed_by_expert_panel",
    "practice_guideline",
]


def review_status_rows(arm: str) -> list[tuple[str, dict]]:
    """Rows spanning ClinVar's review-status scale.

    ``ClinVar_VCF_CLNREVSTAT`` is parsed, not displayed verbatim, so it is one of the
    columns that must carry real vocabulary rather than be left to the long-tail fill —
    a column populated with ``constructed_<name>`` strings would read as an unrecognised
    review status on every row, which is a state the scale does not have.

    Zero stars is not the same as absent, so both are here: the first row asserts a call
    with no criteria behind it, while every other row in the file leaves the column blank.
    The ``CLNSIG`` value is held at the neutral Benign throughout, so these rows say
    nothing about pathogenicity and cannot move a verdict.
    """
    genes = IN_LIST_SOMATIC if arm == "somatic" else IN_LIST_GERMLINE
    # The somatic arm carries the whole scale; the germline arm carries the three levels a
    # germline report actually turns on — a conflicted call, an expert panel's, and a
    # practice guideline's — because those are the ones whose ordering has to hold.
    wanted = (
        CLINVAR_REVIEW_STATUS
        if arm == "somatic"
        else [CLINVAR_REVIEW_STATUS[2], CLINVAR_REVIEW_STATUS[4], CLINVAR_REVIEW_STATUS[5]]
    )
    return [
        (
            f"clnrevstat_{status.split(',')[0]}_{i + 1}",
            {
                "ClinVar_VCF_CLNREVSTAT": status,
                "Hugo_Symbol": genes[i % len(genes)],
            },
        )
        for i, status in enumerate(wanted)
    ]


#: The column CancerVar writes its whole evidence vector into. The leading and trailing
#: spaces are in the **name and every value**, which is real and load-bearing: the app
#: matches this column by substring and so survives the padding by luck rather than by
#: design.
CANCERVAR_EVIDENCE_COLUMN = " CancerVar: CancerVar and Evidence "


def cancervar_evidence_rows() -> list[tuple[str, dict]]:
    """Somatic rows carrying CancerVar's evidence vector beside its tier.

    The tier is printed *inside* the evidence string as well as standing in its own
    column, and the app reads the string. Left to the long-tail fill this column would
    carry ``constructed_`` text on a quarter of the rows, which parses as no tier at all —
    so the shapes are written out: a vector whose criteria sum is negative, one summing to
    zero beside a real tier, a genuine Tier I, and a tier with no vector behind it.

    The ``CancerVar`` column agrees with each string, because these rows are not about the
    two disagreeing — that state is measured where it belongs, on the dedicated fixtures
    in ``tests/fixtures/cancervar/``. The tier-with-no-vector shape lives there too, and is
    deliberately not repeated here.
    """
    shapes = [
        ("cancervar_evidence_negative_sum", "Tier_IV_benign",
         " CancerVar: -2#Tier_IV_benign EVS=[0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -1, 0] "),
        ("cancervar_evidence_zero_sum", "Tier_IV_benign",
         " CancerVar: 0#Tier_IV_benign EVS=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] "),
        ("cancervar_evidence_tier_i", "Tier_I_strong",
         " CancerVar: 9#Tier_I_strong EVS=[1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0] "),
    ]
    return [
        (
            name,
            {
                "CancerVar": tier,
                CANCERVAR_EVIDENCE_COLUMN: evidence,
                # Sub-threshold VAF on the two rows whose tier is in the keep list, so the
                # tier decides the row through the rescue rather than through a guideline
                # block these rows are not meant to be about.
                **({"tumor_f": 0.005} if tier in CANCERVAR_KEEP else {}),
                "Hugo_Symbol": IN_LIST_SOMATIC[i % len(IN_LIST_SOMATIC)],
            },
        )
        for i, (name, tier, evidence) in enumerate(shapes)
    ]


def alphamissense_rows() -> list[tuple[str, dict]]:
    """Rows where ``am_pathogenicity`` and ``AlphaMissense_score`` disagree.

    They are different annotations — Funcotator's AlphaMissense call and dbNSFP's — and
    ``components.variant_detail`` reads ``am_pathogenicity`` and must not read the other.
    A builder that derives one from the other destroys that witness silently: every guard
    about which column is read still passes, because the two agree everywhere.

    So the scores are set apart deliberately, and in both directions, and one row is left
    with dbNSFP's column blank while Funcotator's is populated — which is the state a file
    annotated by only one of the two is in.
    """
    return [
        ("am_disagrees_dbnsfp_lower",
         {"am_pathogenicity": 0.91, "AlphaMissense_score": 0.12,
          "AlphaMissense_pred": "B", "AlphaMissense_rankscore": 0.2,
          "Hugo_Symbol": "MSH2"}),
        ("am_disagrees_dbnsfp_higher",
         {"am_pathogenicity": 0.08, "AlphaMissense_score": 0.97,
          "AlphaMissense_pred": "P", "AlphaMissense_rankscore": 0.99,
          "Hugo_Symbol": "ROS1"}),
        ("am_only_funcotator",
         {"am_pathogenicity": 0.44, "AlphaMissense_score": "",
          "Hugo_Symbol": "PIK3CA"}),
    ]


def legacy_gene_rows(arm: str) -> list[tuple[str, dict]]:
    """Rows carrying the symbols a *legacy* parameter file names.

    ``test_param_migration`` migrates a pre-parity parameter file and asks whether the
    migrated gene restriction reaches ``filter_genes`` at all. It measures that on the
    criteria path rather than on the union, because a gene restriction here does not
    shrink the union — it moves rows across it, and a test written on the total would pass
    with the migration deleted.

    That needs two things from the fixture, and neither is a count:

    * the legacy symbols must carry criteria-path rows, or the restricted cell is zero and
      the measurement is "the restriction reached nothing" rather than "the restriction
      reached these";
    * every criteria-path row **not** on the legacy list must also be pathogenic, so the
      rescue re-admits exactly what the restriction drops and the union coincides.

    The second is why these rows are ``Pathogenic`` on ClinVar as well as guideline-
    passing: a row that were only guideline-passing would leave the union moving with the
    restriction, and the docstring's argument about totals would stop describing this set.
    """
    genes = LEGACY_GENES_SOMATIC if arm == "somatic" else LEGACY_GENES_GERMLINE
    trigger = {"ESCAT": "IA"} if arm == "somatic" else {"RENOVO_Class": "HP Pathogenic"}
    return [
        (
            f"legacy_gene_{gene}",
            {
                "Hugo_Symbol": gene,
                "ClinVar_VCF_CLNSIG": "Pathogenic",
                # Above the legacy file's 0.05 VAF threshold and its 30 min_depth, and
                # under its 0.01 population-frequency cut, so the row is on the criteria
                # path under the *migrated* parameters rather than only the contract ones.
                "tumor_f": 0.4,
                **_frequency(0.0005),
                **trigger,
            },
        )
        for gene in genes
    ]


def overlap_rows(arm: str) -> list[tuple[str, dict]]:
    """Rows excluded by more than one reason at once.

    The run recap tells the reader its exclusion counts sum to more than the number left
    out, because a row can be excluded twice over. That caption is load-bearing — a reader
    who adds the counts and overshoots has otherwise been told something that looks wrong —
    and ``test_run_recap`` refuses to assert it against a fixture where the counts do not
    actually overlap.

    Every witness row in this set is single-reason by construction, which is the rule that
    makes the mutation instrument work and is exactly what leaves no overlap behind. These
    rows are therefore separate and deliberately over-determined: each fails depth *and*
    carries an excluded classification, so it lands in two buckets at once. They witness
    no divergence and are not meant to.
    """
    trigger = {"ESCAT": "IA"} if arm == "somatic" else {"RENOVO_Class": "LP Pathogenic"}
    return [
        (
            f"excluded_twice_{vc.lower()}",
            {
                "Variant_Classification": vc,
                "t_alt_count": 4,
                "t_ref_count": 6,
                "DP": 10,
                **trigger,
            },
        )
        for vc in EXCLUDED_VC
    ]


def synthetic_supplement(arm: str) -> list[tuple[str, dict]]:
    """The shapes a reference run could not supply — missing numerics above all.

    The reference had zero blanks in ``tumor_f`` / ``t_alt_count`` / ``t_ref_count`` /
    ``DP`` across all 181,566 rows, so these rows are the *only* witness for the two
    NaN-tolerance divergences. A constructed set has no trouble here: this is the half of
    the fixture set that was always built rather than cut.
    """
    rows: list[tuple[str, dict]] = []

    def add(name, **over):
        rows.append((name, over))

    trigger = {"ESCAT": "IA"} if arm == "somatic" else {"RENOVO_Class": "LP Pathogenic"}
    vaf_threshold = VAF_SOMATIC if arm == "somatic" else VAF_GERMLINE

    add("control_guideline_pass", **trigger)
    add("control_guideline_fail")
    # A missing depth component: the sum is NaN, `>=` is False, the pipeline drops the row
    # even at min_depth 0 — which is what makes the *_depth_0 case more than "unfiltered".
    add("missing_t_alt_count", t_alt_count="", **trigger)
    add("missing_t_ref_count", t_ref_count="", **trigger)
    add("missing_both_depth_components", t_alt_count="", t_ref_count="", **trigger)
    # ``DP`` absent from the row: PASS, because the pipeline never reads it. The row that
    # tells a DP-based depth filter apart from the sum-based one.
    add("missing_DP", DP="", **trigger)
    add("missing_tumor_f", tumor_f="", **trigger)
    add("vaf_exactly_threshold", tumor_f=vaf_threshold, **trigger)
    add("dp_high_sum_low", t_alt_count=10, t_ref_count=20, DP=100, **trigger)
    add("dp_low_sum_high", t_alt_count=40, t_ref_count=40, DP=30, **trigger)
    # The ``;`` branch of the ClinVar splitter, which never occurs in the reference.
    add("clnsig_semicolon_pathogenic", ClinVar_VCF_CLNSIG="Pathogenic;other")
    add("clnsig_semicolon_likely", ClinVar_VCF_CLNSIG="other;Likely_pathogenic")
    add("clnsig_semicolon_benign", ClinVar_VCF_CLNSIG="Benign;other")
    # ``.`` in a frequency column, in every fixture: ``test_numeric_columns`` requires each
    # one to exercise the coercion, and refuses to pass if it cannot.
    add("frequency_dot", Freq_esp6500siv2_all=".", Freq_1000g2015aug_all=".", **trigger)

    if arm == "somatic":
        for value in ["I", "II", "III"]:
            add(f"escat_bare_{value}", ESCAT=value)
        add("vc_outside_app_vocabulary", Variant_Classification="Start_Codon_SNP",
            **trigger)
    else:
        add("escat_only", ESCAT="IA")
        # ``RENOVO_Class`` is a *guideline* trigger, never a rescue one, so this row stays
        # inside the supplement's rescue-free invariant while still being a one-field edit.
        add("renovo_rescue_candidate", RENOVO_Class="HP Pathogenic", tumor_f=0.05)
        add("vc_outside_app_vocabulary", Variant_Classification="Start_Codon_SNP",
            **trigger)
    # No rescued row in either supplement, deliberately: ``test_filter_entry_point``
    # records these two as unable to separate the criteria path from the union, and a
    # fixture that quietly gained a rescued row would silently invalidate that record.
    # The rescue-only shapes live in the reference-replacement files instead.
    return rows


def civic_rows() -> list[tuple[str, dict]]:
    """The CIViC matrix: every row on the same guideline-*failing* template, so the CIViC
    value alone decides the verdict.

    Values are list-reprs because that is the shape the annotation writes. ``list_c_and_d``
    and ``list_a_and_e`` are the witnesses for divergence #12: a substring match finds the
    kept level inside the repr, an ``.isin()`` on the whole string does not.
    """
    values = [
        ("list_a_only", "['A']"),
        ("list_b_only", "['B']"),
        ("list_c_only", "['C']"),
        ("list_d_only", "['D']"),
        ("list_a_and_b", "['A', 'B']"),
        ("list_c_and_d", "['C', 'D']"),
        ("list_a_and_e", "['A', 'E']"),
        ("list_d_and_e", "['D', 'E']"),
        ("empty_list", "[]"),
        ("bare_a", "A"),
        ("bare_c", "C"),
        ("blank", ""),
    ]
    return [
        (
            name,
            {
                # One row carries the frequency ``.`` this fixture is also required to
                # exercise; it sits on ``blank``, whose CIViC verdict is unaffected.
                **({"Freq_esp6500siv2_all": ".", "Freq_1000g2015aug_all": "."}
                   if name == "blank" else {}),
                "CIViC_Evidence_Level": value,
                "CIViC_Evidence_Rating": 3,
                "CIViC_Entity_Disease": "constructed_disease",
                "CIViC_Variant_URL": "https://example.invalid/variant",
                "CIViC_Entity_URL": "https://example.invalid/entity",
                "CIViC_Entity_Status": "accepted",
            },
        )
        for name, value in values
    ]


def gnomad_genome_rows() -> list[tuple[str, dict]]:
    """Six rows that take ``compute_keep``'s ``gnomAD_genome_AF*`` branch by column
    presence alone. Not a genome-build claim: the values are constructed, the rows are
    hg19-shaped, and the branch is taken because the columns exist.
    """
    return [
        (
            f"gnomad_genome_{i + 1}",
            {
                "gnomAD_genome_AF": value,
                "gnomAD_genome_AF_raw": value,
                "Hugo_Symbol": (IN_LIST_SOMATIC + OFF_LIST)[i],
                **({"Freq_esp6500siv2_all": ".", "Freq_1000g2015aug_all": "."}
                   if i == 0 else {}),
            },
        )
        for i, value in enumerate([0.0, 0.00002, 0.001, 0.02, 0.3, 0.75])
    ]


def dot_numeric_rows() -> list[tuple[str, dict]]:
    """``.`` in a numeric column: the pipeline raises ``TypeError``, the app refuses.

    One ``.`` per column, in the three the pipeline reads plus ``DP``, which it never
    reads — so the app's refusal must name exactly ``t_alt_count``, ``t_ref_count`` and
    ``tumor_f``, and not ``DP``.
    """
    return [
        ("dot_t_alt_count", {"t_alt_count": "."}),
        ("dot_t_ref_count", {"t_ref_count": "."}),
        ("dot_tumor_f", {"tumor_f": "."}),
        ("dot_DP", {"DP": "."}),
        # A ``.`` in a *frequency* column alongside the numeric ones, so this fixture also
        # shows that the app's refusal is about filter inputs and nothing else: the
        # frequency columns are coerced, never refused.
        ("clean", {"ESCAT": "IA", "Freq_esp6500siv2_all": ".",
                   "Freq_1000g2015aug_all": "."}),
    ]


# ---------------------------------------------------------------------------
# Writing
# ---------------------------------------------------------------------------

CHROMS = [f"chr{n}" for n in list(range(1, 23)) + ["X"]]

#: Columns the long-tail fill must leave alone, and why each one.
#:
#: ``gnomAD_genome_AF*`` are a column-*presence* probe: they exist only in
#: ``somatic_gnomad_genome.maf`` and their values there are set per row. The frequency
#: columns are excluded because ``_no_frequency()`` rows depend on being blank *on those
#: rows*, and a fill that reads "this column is blank everywhere" would not touch them
#: anyway — they are listed so that stays true if a future row set changes.
NEVER_FILLED = set(FREQUENCY_COLUMNS) | {
    "gnomAD_genome_AF",
    "gnomAD_genome_AF_raw",
}


def _filler(column: str, index: int) -> str:
    """A deterministic constructed value for a column nothing else populates.

    Deterministic from the column name alone — ``zlib.crc32`` rather than ``hash()``,
    which is salted per process and would make two builds differ. Numeric columns get a
    number so the frame's dtype stays what the annotator's own values would make it;
    everything else gets a string that says what it is.
    """
    seed = zlib.crc32(column.encode()) + index
    if COLUMN_PROFILES.get(column, {}).get("kind") == "numeric":
        return str(round((seed % 10_000) / 10_000, 4))
    return f"constructed_{column}_{seed % 100:02d}"


def fill_long_tail(columns: list[str], rows: list[dict]) -> None:
    """Populate every column no row would otherwise put a value in.

    THE CASE THIS EXISTS FOR. A constructed builder that writes only the columns it
    reasons about leaves the rest blank, and "the rest" is most of the file: 387 of 427
    columns, against the reference set's 7. That is not a smaller fixture, it is a
    different *shape* — one in which a column being present stops being distinguishable
    from its being absent.

    The distinction is load-bearing rather than cosmetic. ``frequency_mask`` treats a
    column that is present-but-blank differently from one that is missing entirely — a
    blank could not sink a row, while an absent column was never in the mask's list at all
    — and ``config/presets.py`` records the measurement where that difference decided 17
    real variants. A fixture set in which nearly every column is blank cannot witness it,
    and worse, cannot fail if the app stops honouring it.

    So each such column gets its own constructed values on some rows and its own blank
    token on the others. Which token is the annotator's choice, not this file's: it is read
    from ``column_profiles.json``, where ``.`` and ``""`` are recorded per column exactly
    as the reference software emitted them. Both conventions therefore survive the swap,
    which they would not if the fill picked one.

    **Density is the annotator's, not this file's.** An earlier draft spread each column
    over a stride derived from its own name, which reproduced the blanking convention while
    inventing the texture — a column the annotation fires on 3% of rows and one it fires on
    97% of them came out looking alike, and "the coverage per arm is reproduced" would have
    been a claim about nothing. ``column_profiles.json`` records the measured ``fill`` per
    column instead, and the rows chosen are evenly spaced across the file so a fixture of 5
    rows and one of 94 both land on it.

    Which rows is deterministic, so the files stay reproducible, and every filled column is
    guaranteed at least one populated row even where the measured fill rounds to nothing —
    a column populated on no row would be the very thing this is preventing.
    """
    if not rows:
        return
    for column in columns:
        if column in NEVER_FILLED:
            continue
        if any(row.get(column) not in ("", ".", None) for row in rows):
            continue
        profile = COLUMN_PROFILES.get(column, {})
        blank = profile.get("blank", "")
        wanted = max(1, round(profile.get("fill", 0.25) * len(rows)))
        # Evenly spaced rather than the first N, so a column's values are spread through
        # the file the way an annotation firing on some variants and not others is.
        step = len(rows) / wanted
        populated = {min(len(rows) - 1, int(n * step)) for n in range(wanted)}
        for index, row in enumerate(rows):
            row[column] = _filler(column, index) if index in populated else blank


def materialise(columns: list[str], arm: str, specs: list[tuple[str, dict]]) -> list[dict]:
    """Turn ``(cell, overrides)`` pairs into full rows over ``columns``.

    Coordinates are assigned by index rather than carried from anywhere: a stable spread
    over chr1..chrX at 10 kb spacing, which keeps the four-column variant key unique
    without any real position appearing. ``Genome_Change`` is derived from the assigned
    coordinate, not from a template — the cell the old fixtures leaked the real position
    through.
    """
    out = []
    for index, (cell, overrides) in enumerate(specs):
        row = {column: "" for column in columns}
        row.update(neutral(arm))
        chrom = CHROMS[index % len(CHROMS)]
        position = 1_000_000 + index * 10_000
        row["Chromosome"] = chrom
        row["Start_Position"] = position
        row["End_Position"] = position
        row.update(overrides)
        # ``RENOVO_pls`` is not free: ``test_clinical_summary`` measures each class's score
        # band over 211,744 real rows and fails a fixture whose score sits outside its own
        # band, so the score follows the class rather than being set beside it.
        if row.get("RENOVO_Class") in RENOVO_SCORES and "RENOVO_pls" not in overrides:
            row["RENOVO_pls"] = RENOVO_SCORES[row["RENOVO_Class"]]
        # ``am_class`` is an exact function of ``am_pathogenicity`` at 0.34 / 0.564, which
        # ``test_alphamissense`` asserts. Derived here for the same reason ``RENOVO_pls``
        # is: a class set beside its score is a fixture that disagrees with itself.
        if "am_class" not in overrides:
            row["am_class"] = _am_class(row.get("am_pathogenicity"))
        ref, alt = row["Reference_Allele"], row["Tumor_Seq_Allele2"]
        row["Genome_Change"] = f"g.{chrom}:{position}{ref}>{alt}"
        # The cell name travels in a column the filters never read, so a diff can say
        # which constructed shape moved. ``Otherinfo`` is in KEEP, so it survives to the
        # pipeline's output and stays readable in a failure.
        row["Otherinfo"] = f"constructed_cell={cell}"
        out.append({column: row.get(column, "") for column in columns})
    return out


def _am_class(score) -> str:
    """AlphaMissense's own three-way call from the score, at the published boundaries."""
    try:
        value = float(score)
    except (TypeError, ValueError):
        return ""
    if value < AM_BENIGN_BELOW:
        return "likely_benign"
    if value > AM_PATHOGENIC_ABOVE:
        return "likely_pathogenic"
    return "ambiguous"


def fill_pipeline_verdict(rows: list[dict]) -> None:
    """Write ``variantalker_naive`` as the pipeline's clinical summary would.

    Not a constant, and not a guess. ``test_clinical_summary`` asserts that this column
    equals the app's own ``generate_clinical_summary`` on every row — that is the whole
    point of the column: it is the pipeline's answer, held beside the app's so the two can
    be compared. So the generator computes it by *calling the app*, then inverting the
    label back to the pipeline's ladder name.

    A row the app can say nothing about would have no ladder name to write. None occurs
    here — the neutral template is Benign on all three sources — and if one ever does, this
    raises rather than writing a placeholder the module would then have to excuse.

    Called after ``fill_long_tail`` on purpose: the summary reads the row, so a column the
    fill populates has to be visible to it or the committed value would describe a row that
    is not the one shipped.
    """
    import sys

    sys.path.insert(0, str(HERE.parent.parent.parent))
    import pandas as pd
    from components.clinical_summary import _SUMMARY_LABELS, generate_clinical_summary

    inverse = {label: name for name, label in _SUMMARY_LABELS.items()}
    for row in rows:
        series = pd.Series({k: (pd.NA if v == "" else v) for k, v in row.items()})
        label = generate_clinical_summary(series)
        name = inverse.get(label)
        if name is None:
            raise SystemExit(
                f"row {row.get('Otherinfo')} produces the label {label!r}, which is not one "
                "of the pipeline's ladder names — give it a clinical source"
            )
        row["variantalker_naive"] = name


def write_maf(path: Path, arm: str, columns: list[str], rows: list[dict]) -> None:
    lines = comment_block(arm, path.name)
    lines.append("\t".join(columns))
    for row in rows:
        lines.append("\t".join(str(row[column]) for column in columns))
    path.write_text("\n".join(lines) + "\n")


PLAN = [
    ("somatic_reference.maf", "somatic", somatic_rows),
    ("germline_reference.maf", "germline", germline_rows),
    ("somatic_synthetic.maf", "somatic", lambda: synthetic_supplement("somatic")),
    ("germline_synthetic.maf", "germline", lambda: synthetic_supplement("germline")),
    ("somatic_civic.maf", "somatic", civic_rows),
    ("somatic_gnomad_genome.maf", "somatic", gnomad_genome_rows),
    ("somatic_dot_numeric.maf", "somatic", dot_numeric_rows),
]

#: The per-arm shape the reference set had, reproduced. Asserted at build time rather than
#: left to a test, so a row added without thinking is caught where it is added.
EXPECTED_SHAPE = {
    "somatic_reference.maf": (82, 427),
    "germline_reference.maf": (94, 380),
}


def build(out_dir: Path) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)

    manifest = {}
    for name, arm, make_rows in PLAN:
        columns = FIXTURES[name]["columns"]
        specs = make_rows()
        rows = materialise(columns, arm, specs)
        fill_long_tail(columns, rows)
        fill_pipeline_verdict(rows)
        expected = EXPECTED_SHAPE.get(name)
        if expected and (len(rows), len(columns)) != expected:
            raise SystemExit(
                f"{name}: built {len(rows)} rows at {len(columns)} columns, but the "
                f"reference set's shape is {expected[0]} rows at {expected[1]} columns"
            )
        write_maf(out_dir / name, arm, columns, rows)
        manifest[name] = {
            "rows": len(rows),
            "columns": len(columns),
            "cells": [cell for cell, _ in specs],
        }

    (out_dir / "genes_somatic.txt").write_text("\n".join(GENES_SOMATIC) + "\n")
    (out_dir / "genes_germline.txt").write_text("\n".join(GENES_GERMLINE) + "\n")
    (out_dir / "genes_somatic_mixed_case.txt").write_text(
        "\n".join(GENES_SOMATIC_MIXED_CASE) + "\n"
    )
    # ``BUILD.json``, not ``MANIFEST.json``: the manifest the suite reads is written by
    # ``record_manifest.py`` in a second step, because its per-row verdicts have to be
    # *measured* by running the pipeline rather than declared here. Two writers of one
    # filename made the build look non-deterministic when only the second step had run.
    (out_dir / "BUILD.json").write_text(json.dumps(manifest, indent=1) + "\n")
    return manifest


BUILT_FILES = [name for name, _, _ in PLAN] + [
    "genes_somatic.txt", "genes_germline.txt", "genes_somatic_mixed_case.txt",
    "BUILD.json",
]


def check(committed: Path) -> int:
    """Rebuild into a temporary directory and compare, file by file."""
    with tempfile.TemporaryDirectory() as tmp:
        rebuilt = Path(tmp)
        build(rebuilt)
        differing = [
            name for name in BUILT_FILES
            if not (committed / name).exists()
            or not filecmp.cmp(committed / name, rebuilt / name, shallow=False)
        ]
    if differing:
        print("the committed set does not match what the generator produces:")
        for name in differing:
            print(f"  {name}")
        return 1
    print(f"{len(BUILT_FILES)} files match what the generator produces")
    return 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", type=Path, help="write the fixture set here")
    parser.add_argument(
        "--check",
        action="store_true",
        help="rebuild into a temporary directory and compare against the committed set",
    )
    args = parser.parse_args(argv)
    if args.check:
        return check(args.out or HERE)
    if args.out is None:
        parser.error("one of --out or --check is required")
    manifest = build(args.out)
    for name, record in manifest.items():
        print(f"{name}: {record['rows']} rows, {record['columns']} cols")
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
