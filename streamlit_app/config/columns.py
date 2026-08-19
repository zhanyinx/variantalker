"""Everything the app knows about columns: which ones a filtered table has and in what
order, which ones a MAF must carry, and what each one means.

Two orderings live here and they answer different questions. :func:`resolve_visible_columns`
decides what the grid *opens* with — selection and order together. :data:`PRIORITY_COLUMNS`
and :func:`prioritize_columns` order a frame that already has its columns, and reach the
user only through the grid's "Show all columns" checkbox and the "Add columns" dropdown.
The second arrived from ``components/`` with issue #76, as the one pure column-ordering
decision that was sitting in a UI module; every other ordering in ``components/`` is a
*viewport* decision — pixel widths, which columns freeze at the left edge — and stays with
the widget it configures, for the reason ``variant_table._PINNED_LEFT`` gives.

The first is the module's original job and the rest arrived with it. :data:`REQUIRED_COLUMNS`,
:data:`OPTIONAL_COLUMNS` and :data:`COLUMN_DESCRIPTIONS` used to live in
``config/constants.py`` and be imported straight back here — the accessors that read them
were already here, and so was the address consumers used. Holding the vocabulary next to the
resolver is one responsibility stated twice, not two responsibilities: both answer "what are
this app's columns", one by listing them and one by deciding which of them to show.

:func:`resolve_visible_columns` is the app's single answer. It is a pure function of
the arm, the CIViC skip flag and the columns actually present, and it drives every
consumer: the grid's default selection and the help pages' schema listing. Since issue
#92 the "shown columns" export is no longer a third consumer beside the grid — it reads
the list the grid resolved, so what it inherits from here it inherits through the grid.
It is deliberately streamlit-free so the parity harness can import and call it rather
than reading a list back out of the source with ``ast``.

The list is built in three parts:

* the **pipeline's** columns, from the vendored :func:`vendor.pipeline_filters.compute_keep`
  — a verbatim copy of the statements in ``bin/filter_variants.py:main()``, not a
  transcription of them;
* minus :data:`PIPELINE_COLUMNS_THE_APP_REPLACES`, the one place the app shows *less* than
  the pipeline emits, and the only subtraction permitted here;
* the **app's** extras, :data:`APP_EXTRA_COLUMNS`, appended after.

So the app's list opens with the pipeline's as an exact prefix **once the replaced columns
are taken out of both**, and is otherwise a deliberate superset. Assert that prefix
element-by-element, never by length: the three lists this function replaced were once
measured at 40 columns against the pipeline's 40 while differing by a substitution, and a
length check is precisely the check that misses it.

The subtraction is new with issue #117 and is deliberately a *named list rather than a
mechanism*. Nothing may drop a pipeline column for being empty, or narrow, or ugly; a name
goes in that list only when the app draws its own answer to the same question, and the
argument for each entry is written beside it.

The lists replaced, all now deleted:

* ``get_filtered_columns`` here — a hand-maintained transcription of ``bin/utils.py:KEEP``
  which had drifted, and which only ``help.py`` ever called;
* ``_DEFAULT_VISIBLE_COLUMNS``, then in ``components/ui_components.py`` — the 25 columns
  that actually reached the user, two of them gnomAD v4 names matching no reference MAF
  on either genome build. That file is gone too: issue #76 split it, and the grid it
  drove is ``components/variant_table.py``;
* a ``KEEP`` in ``utils/main_utils.py``, deleted earlier by issue #32.
"""

import re
import warnings

import pandas as pd


class MissingColumnsWarning(UserWarning):
    """A MAF is missing columns the pipeline's output would carry.

    The pipeline reselects with ``out[keep]`` and raises ``KeyError`` on such a MAF.
    The app is interactive and the user has already waited for the load, so it drops
    the absent columns and says which ones instead of failing the render.
    """


#: The app's own columns, appended after the pipeline's. Every entry is either derived
#: at display time or written by the app's filter, so none of them can come from
#: ``compute_keep``, and none of them may be a name ``compute_keep`` already produces.
#:
#: What is *not* here: ``gnomAD_exome_AF_grpmax`` and ``gnomad41_genome_AF``. Both were
#: in ``_DEFAULT_VISIBLE_COLUMNS``, and both are fiction — they appear in no reference
#: MAF on either genome build, and ``compute_keep`` hardcodes the older
#: ``gnomAD_{exome,genome}_AF[_raw]`` spellings, so it could not surface them even if a
#: MAF carried them. Also not here: ``MAX_AF``, which the deleted ``get_filtered_columns``
#: appended and the pipeline has never emitted.
#: Which of these the app *makes* and which the *file* carries is not narrated here any
#: more: it is declared in :data:`COLUMNS_MAFIGATE_ADDS` and
#: :data:`COLUMNS_THE_MAF_CARRIES_OUTSIDE_THE_KEEP_LIST` below, which a guard holds against
#: this list. Two comments used to say it, nothing checked them, and being the only thing
#: that knew the split is precisely how six columns came to be called "Optional" on a
#: clinician's screen (issue #124). An unguarded restatement of it here would rebuild that.
APP_EXTRA_COLUMNS = [
    "Clinical_Summary",
    "Pathogenicity_Overview",
    "Notes",
    "Sample_Name",
    "dbSNP_ID",
    # The app's verdict and its reason, appended last. Spelled out rather than imported
    # from filters.variant_filters, which would pull pandas into a module that is
    # otherwise pandas-free (see create_column_info_table). The two spellings are held
    # together by test_verdict_and_reason_are_the_last_two_columns, which asserts this
    # tail against MAFIGATE_FILTER and MAFIGATE_REASON themselves.
    "MAFigate_filter",
    "MAFigate_reason",
]


#: The subset of :data:`APP_EXTRA_COLUMNS` that **no MAF ever carries** — the app makes
#: every one of them, on every report, from columns the file does supply.
#:
#: The split this records was a comment until issue #124, and the comment was the only
#: thing that knew it: :func:`column_source_status` classified by membership of
#: :data:`REQUIRED_COLUMNS`/:data:`OPTIONAL_COLUMNS`, these names are in neither, so all
#: six fell off the end of its chain and the column reference called them **"Optional"** —
#: a claim about the *user's file*, false in both directions about a column their file
#: never has and the app always makes.
#:
#: ``dbSNP_ID`` is deliberately **not** here, and it is the whole reason this is a subset
#: rather than a rename of :data:`APP_EXTRA_COLUMNS`. The app appends it because
#: ``compute_keep`` does not, not because it makes it — the *file* supplies it — so
#: "Optional" is exactly right for that one: #120 measured 167 of 290 real MAFs carrying
#: it, and an absent app extra is dropped from the view in silence.
COLUMNS_MAFIGATE_ADDS = [
    "Clinical_Summary",
    "Pathogenicity_Overview",
    "Notes",
    "Sample_Name",
    "MAFigate_filter",
    "MAFigate_reason",
]

#: The rest of :data:`APP_EXTRA_COLUMNS`: names the *file* carries that the pipeline's keep
#: list leaves out, so the app puts them back.
#:
#: Declared rather than inferred as "everything not in :data:`COLUMNS_MAFIGATE_ADDS`",
#: because inferring it is what lets a new app extra join the table classified by accident.
#: ``test_every_app_extra_column_is_declared_on_one_side_or_the_other`` asserts the two
#: halves partition :data:`APP_EXTRA_COLUMNS`, so a seventh entry fails there until whoever
#: added it says which kind it is — the same derivation #120 used for the descriptions, and
#: for the same reason: a hand-listed set at the point of *use* is the thing that did not
#: get updated six times.
COLUMNS_THE_MAF_CARRIES_OUTSIDE_THE_KEEP_LIST = [
    "dbSNP_ID",
]


#: Pipeline output columns the app leaves out of its default view, because it draws its
#: own answer to the same question on the same screen.
#:
#: **The only subtraction from the pipeline's list, and the first this map has made.**
#: Every other decision here adds; this one takes away, so it is a named list with the
#: argument written down rather than a rule anything can invoke. The bar is not "empty",
#: "rarely useful" or "cluttered" — it is that a column of the app's own answers the same
#: question, pinned, in the same table.
#:
#: ``variantalker_naive`` is the pipeline's clinical verdict, written by
#: ``bin/generate_clinical_summary.py`` as a live annotation process and kept by
#: ``compute_keep`` on both arms. ``components/clinical_summary.py`` re-implements that
#: script and draws its answer as ``Clinical_Summary``, pinned first and glyphed — so the
#: report carried the same verdict twice, once with glyphs at the front and once
#: underscore-spelled in the middle, with nothing on screen saying they were twins.
#:
#: **Measured before deciding** (issue #117), over 297 byte-distinct real MAFs on the
#: developer's machine, 176 of them carrying clinical annotation, read with
#: ``vendor.pipeline_utils.read_maf``:
#:
#: * the two columns hold **identical words on 100% of somatic rows** — 89,019 of 89,019,
#:   across the 50 somatic MAFs carrying both — and on 51 of the 52 germline MAFs;
#: * the whole of the disagreement is **one file**, ``1669-01_N.maf``, at 7,508 of its
#:   32,162 annotated rows (22.7%), whose column was written by a pipeline that did not
#:   yet read RENOVO. There the app is the *current* answer and the file's column the
#:   stale one, which is the wrong way round for a tie-breaker;
#: * what the app says and the pipeline cannot — issue #98's six classes below the ladder
#:   — reaches **8 rows in 329,955** (3 Drug Response, 2 No Classification, 1 Protective,
#:   1 Disease Risk, 1 Unrecognised Annotation, every one of them somatic, none germline);
#: * and **74 of the 176 annotated MAFs carry no ``variantalker_naive`` at all**, so on
#:   108,769 rows it is the app's column or nothing.
#:
#: So the pipeline's copy adds nothing on the files that have it, is missing from 42% of
#: the files that do not, and is stale on the one file where it differs. It is a duplicate
#: of a pinned column, which is what puts it here.
#:
#: **Left out of the view, not out of the app.** The column is still in the frame, so the
#: grid's *Show all columns* checkbox shows it, *Add columns* offers it, and #92's
#: download carries it when either does — the file's own content is never dropped, only
#: its place in the columns the report opens with. It keeps its entry in
#: :data:`COLUMN_DESCRIPTIONS` for the same reason: a reader who meets it needs to be able
#: to look it up, and that entry is now the only thing the app says about it.
#:
#: What this is **not**: it is not a claim the pipeline is wrong, and it does not touch
#: :func:`pipeline_output_columns`, which must go on answering what the pipeline emits —
#: the drift guard, the missing-column account and ``circle_sources`` all depend on that
#: staying the pipeline's own answer.
PIPELINE_COLUMNS_THE_APP_REPLACES = ("variantalker_naive",)


class _ColumnsShim(frozenset):
    """The ``out.columns`` that ``compute_keep``'s verbatim statements read.

    They test membership two ways — ``col in out.columns`` and
    ``name in out.columns.values`` — because that is what ``main()`` does. A set covers
    the first in constant time; ``values`` returning self covers the second. Order is
    not needed: ``compute_keep`` only ever asks whether a name is present, and the order
    of the result comes from ``KEEP``, never from the frame. Building a real DataFrame
    here would work too, but would make this module import pandas, which it otherwise
    does not need.
    """

    @property
    def values(self):
        return self


class _AllColumnsPresent:
    """Stands in for ``out.columns`` when the caller has no MAF in hand.

    Reports every name as present, so the help pages get the full schema rather than
    a list shaped by whichever frame happened to be loaded.

    Deliberately not iterable: an infinite set of names has no meaningful iteration, and
    ``compute_keep`` only ever tests membership. If it grows a loop over the columns,
    that must fail loudly here rather than quietly seeing none.
    """

    def __contains__(self, name):
        return True

    @property
    def values(self):
        return self




# The column vocabulary. It lived in `config/constants.py` and was imported straight back
# here, which put the data one module away from the accessors below that read it — and
# `tests/test_config.py` was already importing REQUIRED_COLUMNS and OPTIONAL_COLUMNS *from
# this module* rather than from where they were defined. This module was the address
# consumers already used; now it is also where the data is.
#
# Not private, and not funnelled through the accessors: `page_modules/data_loading.py`
# reads the two category maps directly to decide what it can refuse and what it must warn
# about, and `page_modules/help.py` reads the descriptions directly to build its schema
# table. Routing those through a function would be a behaviour change dressed as tidying,
# so the three dicts are part of this module's surface alongside the accessors.

# Two jobs, and only the first of them is a requirement (issue #127).
#
# ``core`` is what `validate_required_columns` refuses a file over, and the only category the
# word "Required" is true of: everything else here is tolerated, filled or blank, and the app
# carries on. The other three are the *display* columns each arm expects — the ones the same
# function names in its "Not in this MAF" note, after subtracting the filter's own inputs —
# and they are still spelled per arm because that note is drawn for one arm at a time.
#
# What they are no longer is the source of the Help page's "In your MAF?" answers. That job
# went to `_classify_column_source`, which reads what an absence *does* rather than which list
# a column sits in; three of the old answers claimed a requirement nothing enforces, and two
# of them contradicted the Sample Types column on their own row.
REQUIRED_COLUMNS = {
    # Absolutely required for basic analysis
    "core": [
        "Hugo_Symbol",
        "Chromosome",
        "Start_Position",
        "End_Position",
        "Variant_Classification",
        "Variant_Type",
        "Reference_Allele",
        "Tumor_Sample_Barcode",
    ],
    # Both arms read these; `validate_required_columns` subtracts the arm's filter inputs, so
    # what survives from here is `DP` alone.
    #
    # `DP` sat under `# For depth filtering` and does no filtering: the pipeline's depth gate
    # is `(t_alt_count + t_ref_count) >= coverage` and the vendored body reads `DP` on no path
    # of either arm. Nor is it that sum by another name — across 322,913 rows of 157 real MAFs
    # the two numbers differ on 72.7% of them, `DP` being the larger on 57.3%, which is what
    # its own description says it is (every sample, not the tumour). The stale comment had a
    # live source: the widget over the gate was labelled `Minimum Depth (DP)`, naming a column
    # the gate does not read. Issue #127 relabelled it.
    "filtering": [
        "DP",  # shown in the variant dialog as "Depth"; read by no filter
        "tumor_f",  # a filter input, and therefore filled rather than blank when absent
    ],
    # Not "required for somatic analysis": the somatic filter reads `t_alt_count` and
    # `t_ref_count` and fills them when absent, and the two allele columns are read by no
    # filter on either arm. What makes this the *somatic* list is that these are the columns
    # the somatic arm's note asks about — the pipeline emits all four on both arms, which is
    # why the reference table reports them as "Both".
    "somatic": ["Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "t_alt_count", "t_ref_count"],
    # The same, on the other arm, and the whole list is display-only: nothing reads the normal
    # sample's read counts — no filter, no chart, no panel — so they reach the user as table
    # columns and nothing else.
    "germline": ["n_alt_count", "n_ref_count"],
}

# A `clinical` category listed the five annotations the pipeline's filter reads, as a mirror of
# `filters/absent_columns.REQUIRED_INPUTS`, kept in step by a drift guard. Issue #127 deleted
# it: `_classify_column_source` was its only reader, and it now reads the derivation itself, so
# the mirror had nothing left to do. The reason the guard gave for the mirror existing — that
# this module "cannot import the derivation" — was no longer true either; nothing in
# `filters/absent_columns.py` imports `config/`, and the classifier takes it function-locally
# exactly as `pipeline_output_columns` takes `compute_keep`.
#
# What the category was *for* survives in the derivation and is stronger there: issue #39
# measured that filling `CancerVar` neutrally drops 70% of the somatic report and
# `RENOVO_Class` 39% of the germline one, so "optional" was the opposite of true for them —
# and now no filter input can be called optional whatever anyone writes in a list here, which
# `test_no_filter_input_is_called_optional` holds.

# Optional columns that enhance analysis but aren't required
OPTIONAL_COLUMNS = {
    # The one genuinely optional filter input, and the only one the *pipeline* guards: its
    # own check_civic_column_exists skips the CIViC clause when the column is absent, so a
    # MAF without it is filtered by both tools alike and stays on parity.
    "civic": [
        "CIViC_Evidence_Level",
    ],
    # No `population_frequency` entry, deliberately, since issue #126. It held a
    # hand-written second copy of `filters.variant_filters.FREQUENCY_COLUMNS` -- the same
    # eight names, derived from nothing -- and its only live reader was
    # `page_modules/data_loading.py`'s "does this MAF carry any frequency column?" check,
    # which now reads the filter's own list. So the load-time account and the filter cannot
    # disagree about what a frequency column is, which is exactly what dropping `MAX_AF`
    # from one of two copies would have made them do.
    #
    # What the *reader* is told is unmoved by the deletion, and that had to be re-checked
    # against #124 rather than assumed: `_classify_column_source` iterates this dict, so the
    # five of these names that reach the reference table (`Freq_ExAC_ALL`,
    # `Freq_esp6500siv2_all`, `Freq_1000g2015aug_all` and the two `_raw` columns) move from
    # its `optional` token to its `uncategorised` one. #124 made those two deliberately
    # separate facts, but `_COLUMN_SOURCE_PROSE` words both as **"Optional"**, so the column
    # reference renders exactly what it rendered before -- measured, not argued.
    #
    # And #124's distinction is not hollowed out by the move: `optional` still has five
    # genuine members on screen, `CIViC_Evidence_Level` and the four `annotation` names below.
    # Worth stating, because a guard whose subjects quietly leave is this repo's recurring way
    # of going vacuous.
    #
    # One list of these names does survive, and it is not a copy of the filter's: the
    # category chain in `create_column_info_table` spells three `Freq_*` names to file them
    # under "Population Frequency". That answers *what kind of column is this*, for a column
    # already on screen -- a display question, not a claim about what the filter reads -- so
    # it is not the drift this deletion closes. It would be, if it grew a name the filter
    # does not read, which is why `MAX_AF` left it too.
    "annotation": [
        "AAChange.refGene",
        "cDNA_Change",
        "Protein_Change",
        "Transcript_Exon",
    ],
}


# Column descriptions for MAF file and annotation fields
COLUMN_DESCRIPTIONS = {
    # Sample and project information
    "Tumor_Sample_Barcode": "Unique identifier for the tumor sample aliquot from which the mutation was called",
    "Matched_Norm_Sample_Barcode": "Unique identifier for the matched normal sample aliquot used as reference",
    "project_id": "Project or study identifier for data organization and tracking",
    # Gene and transcript information
    "Hugo_Symbol": "HUGO Gene Nomenclature Committee approved gene symbol (always in caps, 'Unknown' for non-gene regions)",
    "Annotation_Transcript": "Transcript identifier used for annotation, typically Ensembl or RefSeq ID",
    # Genomic coordinates
    "Chromosome": "Chromosome on which the variant is located (e.g., chr1, chr2, etc.)",
    "Start_Position": "Genomic start coordinate of the variant (1-based coordinate system)",
    "End_Position": "Genomic end coordinate of the variant (1-based coordinate system, same as start for SNVs)",
    # Variant classification and type
    "Variant_Classification": "Functional effect of the variant on the transcript (e.g., Missense_Mutation, Nonsense_Mutation, Silent)",
    "Variant_Type": "Type of sequence alteration (SNP=single nucleotide, DNP=dinucleotide, TNP=trinucleotide, INS=insertion, DEL=deletion)",
    # Allele information
    "Reference_Allele": "Reference genome allele at this position ('+' strand)",
    "Tumor_Seq_Allele1": "First allele observed in tumor sequencing data",
    "Tumor_Seq_Allele2": "Second allele observed in tumor sequencing data",
    # Protein and transcript changes
    "AAChange.refGene": "Amino acid change annotation from RefGene database including transcript and protein information",
    "cDNA_Change": "Description of the variant in coding DNA sequence using HGVS nomenclature",
    "Codon_Change": "Change in DNA codon sequence showing reference and variant codons",
    "Protein_Change": "Amino acid change in protein sequence using standard three-letter codes",
    "Transcript_Exon": "Exon number where the variant is located within the transcript",
    # Sequencing metrics
    "tumor_f": "Variant allele frequency (VAF) in tumor sample (ratio of variant reads to total reads)",
    "DP": "Total read depth at the variant position across all samples",
    "t_alt_count": "Number of reads supporting the variant allele in tumor sample",
    "t_ref_count": "Number of reads supporting the reference allele in tumor sample",
    "n_alt_count": "Number of reads supporting the variant allele in normal sample",
    "n_ref_count": "Number of reads supporting the reference allele in normal sample",
    # The two clinical verdicts, and the whole of what the app says to tell them apart
    # (issue #117). One is derived here and shown; the other is written into the MAF by
    # the pipeline and no longer shown by default, so this pair of lines is the only place
    # a reader can find out that they answer the same question — and which one they are
    # looking at when "Show all columns" puts both on screen.
    #
    # Both are written for a clinician: no issue numbers, no module names, and the
    # divergence is described by *when* each was worked out rather than by which source
    # each reads, because "the same five columns" is not a difference a reader can act on
    # and "annotated earlier" is.
    "Clinical_Summary": "MAFigate's clinical verdict for this variant: the strongest classification any of its annotation sources gives it, worked out when you loaded the file",
    "variantalker_naive": "The same clinical verdict as Clinical_Summary, as written into the MAF when it was annotated. Normally identical; it can differ on a file annotated by an earlier version of the pipeline, and then Clinical_Summary is the current reading",
    # Clinical significance databases
    "ClinVar_VCF_CLNSIG": "Clinical significance annotation from ClinVar database (Pathogenic, Likely_pathogenic, Benign, etc.)",
    # Undescribed until issue #219, and by then the variant panel was drawing it: issue #204 put
    # this column on ClinVar's badge as star glyphs beside the significance. A reader who saw
    # `Pathogenic ★★☆☆` had nowhere to look it up — the column is not in the resolver's set
    # either, so it is not a row of the reference table, which is why the help page's
    # `_KEY_CLINICAL_COLUMNS` shortlist names it.
    #
    # "How many stars" is the whole point of the description, and it is the panel's authority:
    # the star hierarchy is ClinVar's own published scale, dev-confirmed as ordered (issue #200),
    # and `components/clinical_badges.CLINVAR_REVIEW_STARS` is the mapping. Zero stars is named
    # explicitly because that is the reading most likely to look like a bug — a variant showing
    # `Pathogenic ☆☆☆☆` is ClinVar reporting a classification submitted with no assertion
    # criteria, not MAFigate failing to read one.
    "ClinVar_VCF_CLNREVSTAT": "How ClinVar's classification was reviewed, on ClinVar's own star scale - the variant panel shows it as stars beside the significance, where zero stars means ClinVar holds a classification submitted without assertion criteria",
    "CancerVar": "Cancer-specific variant interpretation using AI-based CancerVar classification (Tier_I_strong to Tier_IV_benign)",
    "ESCAT": "ESMO Scale for Clinical Actionability of molecular Targets - clinical evidence level (IA to V)",
    "ESCAT_TISSUE": "Tissue or cancer type context for ESCAT classification",
    "ESCAT_CANCER": "Specific cancer diagnosis for ESCAT clinical actionability assessment",
    # CIViC database annotations
    # Not a per-level gloss: the five levels are defined once, from CIViC's own documentation,
    # in `vocabularies.CIVIC_DEFINITIONS`, and the help page renders them from there. This line
    # said "A=Validated, B=Clinical, C=Case study, D=Preclinical" — four of five, `E` missing —
    # and called the scale a clinical evidence level, which is what invited reading it as a
    # significance call (issue #109). One cell holds one entry per CIViC evidence item, so the
    # description says so: a reader matching this against their own header row meets a list.
    "CIViC_Evidence_Level": "How strong the CIViC evidence is, from A (a proven association) to E (an indirect one) - one entry per CIViC evidence item",
    "CIViC_Evidence_Rating": "Quality rating of CIViC evidence (1-5 stars, 5 being highest quality)",
    "CIViC_Entity_Disease": "Disease or cancer type associated with the CIViC evidence",
    "CIViC_Variant_URL": "Direct URL link to the variant page in CIViC database",
    "CIViC_Entity_URL": "Direct URL link to the evidence item in CIViC database",
    "CIViC_Entity_Status": "Status of the CIViC evidence entry (accepted, submitted, rejected)",
    # AlphaMissense, named (issue #219). "Annotation class from external annotation tools" and
    # "Pathogenicity prediction from annotation tools" named no tool at all, on the two columns
    # the variant panel's AlphaMissense section is built from — so a reader who saw
    # "AlphaMissense (in silico missense predictor)" on the panel could not find either column
    # here by that name, and the FAQ that sends them to this table promised a detailed
    # description.
    #
    # Three facts issue #203 measured, each one a way a shorter line would mislead:
    #
    # * `am_pathogenicity` is **not** `AlphaMissense_score`, which is also in these MAFs. Of
    #   55,353 rows carrying both, only 6.3% agree to 1e-9 and 346 rows across 84 files fall in
    #   different AlphaMissense classes — so a reader matching one against the other must be told
    #   they are different annotations rather than left to assume a rounding difference.
    # * `am_class` is an **exact function** of `am_pathogenicity` at AlphaMissense's own 0.34 and
    #   0.564 cutoffs, so it carries no information the score lacks.
    # * 100 of 139 real files spell `am_class` `benign`/`pathogenic` for a tool that publishes
    #   only the `likely_` form. The panel names the band from the score rather than echoing the
    #   class for exactly that reason, and a description that presented the class as
    #   AlphaMissense's own vocabulary would contradict the panel beside it.
    "am_class": "AlphaMissense's own class for this missense variant, derived from am_pathogenicity at its 0.34 and 0.564 cutoffs - some files spell it benign/pathogenic where AlphaMissense itself publishes only likely_benign/likely_pathogenic",
    "am_pathogenicity": "AlphaMissense pathogenicity score (0-1), the value the variant panel's AlphaMissense section reads - not the same annotation as AlphaMissense_score, which some files also carry with a different value",
    "Otherinfo": "Additional annotation information and comments",
    "tumor_tissue": "Tissue type or anatomical site of the tumor sample",
    "cosmic": "COSMIC database identifier for known somatic mutations in cancer",
    # Population frequency data
    "Freq_ExAC_ALL": "Allele frequency in ExAC database across all populations (now part of gnomAD)",
    "Freq_esp6500siv2_all": "Allele frequency in NHLBI Exome Sequencing Project (ESP) across all populations",
    "Freq_1000g2015aug_all": "Allele frequency in 1000 Genomes Project across all populations",
    # None of these four says "recommended" any more, and the raw pair no longer calls
    # itself high-quality (issue #126). MAFigate reads all four the same way, taking the
    # lowest frequency any of them reports, so a description that ranked them was describing
    # advice the filter does not act on -- and "high-quality population reference" was the
    # exact opposite of what a pre-QC frequency is.
    "gnomAD_exome_AF_raw": "Allele frequency in gnomAD exomes, before gnomAD's own genotype quality filters",
    "gnomAD_exome_AF": "Allele frequency in gnomAD exomes, over the genotypes that passed gnomAD's quality filters",
    "gnomAD_genome_AF_raw": "Allele frequency in gnomAD genomes, before gnomAD's own genotype quality filters",
    "gnomAD_genome_AF": "Allele frequency in gnomAD genomes, over the genotypes that passed gnomAD's quality filters",
    # No `MAX_AF` entry, since issue #126 stopped the filter reading it. The description that
    # stood here called it "useful for conservative filtering", a claim about this app that
    # the same ticket made false -- and there is no surface left to correct it on. `MAX_AF` is
    # in neither the resolver's set nor `APP_EXTRA_COLUMNS`, so it has never appeared in the
    # column reference table (verified: 57 rows, none of them this one); the note at
    # `APP_EXTRA_COLUMNS` already records that only the deleted `get_filtered_columns` ever
    # appended it and the pipeline has never emitted it. Help's frequency list is derived from
    # `FREQUENCY_COLUMNS` now, so that reader is gone too. A description nothing can render is
    # dead weight, and correcting one is worse than deleting it: it reads as a live surface.
    # Germline-specific annotations
    "InterVar": "Automated ACMG/AMP guidelines-based variant interpretation for germline variants",
    "RENOVO_Class": "RENOVO algorithm classification for germline variant pathogenicity assessment",
    "RENOVO_pls": "RENOVO pathogenicity likelihood score for germline variants",
    # Quality and filtering
    "filter": "Quality filter flags applied during variant calling and processing (PASS indicates high quality)",
    # Variant identity carried by the MAF but outside the pipeline's keep list.
    #
    # "an rs number" would have been false. Measured over the 290 byte-distinct real MAFs on
    # the developer's machine, the 167 of them that carry this column, 313,257 rows: **5.3%
    # of cells are empty** and **0.26% hold several rs IDs separated by `|`**, which is what
    # a position with more than one dbSNP entry looks like. A reader matching this against
    # their own header row meets both, so the description names both.
    #
    # The first sweep for this said 95 files and 131,925 rows and was wrong in a way worth
    # recording: it split `find` output on whitespace, and these MAFs live under a OneDrive
    # path containing spaces, so two thirds of the corpus was shredded into unreadable
    # fragments and silently skipped. Split on NUL. Both percentages moved; neither
    # qualitative claim did.
    "dbSNP_ID": "dbSNP reference SNP identifier (rs number) for this variant. Empty where the variant is not in dbSNP; a position with more than one dbSNP entry carries them separated by |",
    # The five undescribed columns the app invents, that no MAF supplies and no annotation
    # source explains (issue #120). `Clinical_Summary` is a sixth of that kind and was
    # described by #117; what these five had in common is that nobody had described them.
    #
    # Being the app's own is what makes this table the only place they can be looked up —
    # every other column here is at least a name a user can search for outside the app — and
    # the reference table renders `get_column_description`, so it said "Description not
    # available for Pathogenicity_Overview" about the column it pins **second**.
    #
    # Every one of them borrows wording the app already uses for the same fact rather than
    # inventing a second one, which is the discipline `Clinical_Summary` above and the ESCAT
    # definitions both follow: the reason cells are spelled the way the run report spells
    # them, and the note the way the variant dialog's caption does. A glossary that says it
    # differently is a second story to keep true.
    #
    # What `Notes` deliberately does **not** carry is the dialog's *"do not include patient
    # identifiers"*. That is an instruction, and it belongs where the writing happens; this
    # is the fact the instruction follows from, and it belongs where the column is looked up.
    # Two copies of one safety sentence is how a safety sentence comes to have two wordings.
    #
    # `Pathogenicity_Overview` **points at the key rather than repeating it**. The key beside
    # the grid already names every glyph it can draw and every source it drew them for, and
    # both halves are derived per arm and per file — so a copy here could not be right for
    # every report even on the day it was written, and a transcription of that same key is
    # how `docs/CLINICAL_SUMMARY_FEATURE.md` came to hold three lines that were wrong. It
    # says "marker" and not "circle": since #98 and #104 the strip also draws ⚠️ 💊 🔗 🛡️ ⬜,
    # so "one coloured circle per source" was false of five of the eleven glyphs in the key.
    "Pathogenicity_Overview": "One marker per annotation source for this variant, read left to right. The key drawn above your results table names each position and says what its colours and symbols mean - which sources appear depends on whether you are running somatic or germline, and on what your file carries",
    "Notes": "Your own note about this variant, typed in the variant details panel. A note is attached to the variant rather than to the file, so it appears on any file carrying this variant. Notes live for this session only: download the table to keep them",
    # "or the matched normal's" is true only where the file has that column: the fallback is
    # `Matched_Norm_Sample_Barcode` where it exists and an empty string where it does not,
    # and where `Tumor_Sample_Barcode` itself is absent the column is never derived at all.
    # So the entry says what a blank cell means rather than implying it cannot happen.
    "Sample_Name": "The sample this variant was called in: the tumor sample barcode, or the matched normal's where the tumor one is missing or a placeholder - blank if your file carries neither",
    "MAFigate_filter": "MAFigate's verdict on this variant under the parameters you set: PASS if it is in your report, NOPASS if it was filtered out",
    "MAFigate_reason": "Why this variant got that verdict: criteria - it passed on the filter criteria; both - it met the criteria and is also pathogenic; pathogenic_rescue - it is kept only by pathogenic retention, having not met the criteria; rejected - it did not pass",
}


#: What :func:`get_column_description` returns for a name nobody has described.
#:
#: **Kept rather than raised**, because the accessor is called with names from lists this
#: module does not own — ``page_modules/help.py``'s curated shortlists among them — and a
#: ``KeyError`` there would take a whole page out over one missing sentence. What it must
#: not do is stay *quiet*, and it did: the column reference rendered a plausible sentence
#: for six columns nobody had described, and nothing said so until the table was measured
#: (issue #120).
#:
#: So the loudness lives in the suite rather than on a clinician's screen.
#: ``test_every_column_the_reference_table_lists_has_a_real_description`` builds the table
#: and fails on any row holding this string — and the table is built from the *resolver*,
#: so a new entry in :data:`APP_EXTRA_COLUMNS` joins it automatically and fails there until
#: it is described. That derivation is the guard; a hand-listed set of names would need
#: updating by whoever added the column, which is the thing that did not happen.
#:
#: A named constant, and not an f-string inline below, so that guard can construct the exact
#: fallback rather than sniff for "not available" in prose a description could legitimately
#: contain.
_MISSING_DESCRIPTION = "Description not available for {column}"


def get_column_description(column_name):
    """Get description for a specific column"""
    return COLUMN_DESCRIPTIONS.get(
        column_name, _MISSING_DESCRIPTION.format(column=column_name)
    )


#: The header the column reference draws over :func:`column_source_status`, and the whole
#: of what makes its answers one answer.
#:
#: It said **"Required"** until issue #124, and four of its five values answered that: what
#: your MAF must carry, in degrees. A column the app *derives* is not an input at all, so
#: under that header there was nothing true to say about it and it said "Optional" — which
#: reads as *your MAF may or may not carry this and MAFigate copes either way*, the exact
#: opposite of the truth for a column no MAF has and the app always makes.
#:
#: A fifth value under the old header would have answered a question the header did not
#: ask. Widening the *question* instead is what lets one column carry both kinds of row:
#: every value below is an answer to "is this column in my file, and must it be?", which
#: is what a reader holding their own header row against this table is asking.
#:
#: **The frame's column key, not just its label.** ``st.column_config`` matches on the
#: frame's column name and ignores a key that matches nothing *in silence*, so a rename
#: spelled as a literal at one end would quietly drop the other end's width and label
#: rather than fail. :func:`create_column_info_table` and ``page_modules/help.py`` both
#: read this constant, and ``test_the_reference_table_configures_only_columns_it_actually_has``
#: proves the two agree rather than merely having been written together.
COLUMN_SOURCE_HEADER = "In your MAF?"

#: Prose for each token :func:`_classify_column_source` returns. Split from the classifier
#: on purpose (issue #124): the old function decided and worded in one pass, so "Optional"
#: was returned from **two** places — a genuine :data:`OPTIONAL_COLUMNS` member, and the
#: fallthrough for a name in no category at all — and two different facts spelled the same
#: way is what let the derived columns hide among them for as long as they did.
#:
#: ``optional`` and ``uncategorised`` still render **the same sentence**, and that is the
#: decision rather than an oversight: for a column the *file* supplies, "you do not have to
#: carry this" is true whether or not anyone has filed it under a category. What changed is
#: that the two are now distinguishable in code, so a guard can say that no column the app
#: derives reaches ``uncategorised`` — which no assertion over the rendered string could.
#:
#: **Three values went at issue #127, because they claimed a requirement nothing enforces.**
#: Only ``core`` refuses a file: ``validate_required_columns`` errors and returns ``False``
#: for a missing core column, and for every other category it says "Filtering is unaffected"
#: and carries on. Two things made a per-category answer unfixable as wording:
#:
#: * ``filtering`` and ``somatic`` each held columns of **two** kinds — ``tumor_f``,
#:   ``t_alt_count`` and ``t_ref_count`` are filled by the filter, while ``DP`` and the two
#:   allele columns beside them are read on no path at all — so one string per category could
#:   not be true of both.
#: * ``Required — somatic reports`` and ``Required — germline reports`` **contradicted the
#:   Sample Types column on their own row.** This table is the union of both arms, and every
#:   one of those six columns is emitted on both, so the row read ``Both`` beside a status
#:   naming one arm. The arm question was already answered one column to the left, and
#:   answered differently.
#:
#: So the status is derived from what an absence *does*, per column: refused (``core``),
#: filled by the filter (:data:`~filters.absent_columns.REQUIRED_INPUTS`, both arms), or
#: neither — in which case the column is one the file may or may not carry and the app copes,
#: which is what "Optional" has always said. The five that landed there — ``DP``,
#: ``Tumor_Seq_Allele1``, ``Tumor_Seq_Allele2``, ``n_alt_count``, ``n_ref_count`` — join it
#: on #124's own line, and ``cosmic`` is the precedent measured there: its absence *warns*
#: and its answer is still "Optional", because whether an absence is noticed is a different
#: question from whether the file has to carry the column.
_COLUMN_SOURCE_PROSE = {
    "core": "Required — the file is refused without it",
    # An absent filter input is filled neutrally so the app can still report, and the fill
    # costs the report rows — up to 70% of it (issue #39). See filters/absent_columns.py.
    # Was ``clinical``, keyed off a hand-written mirror of the derivation; #127 reads the
    # derivation itself, which is how ``tumor_f`` and the two tumour read counts came to give
    # the same answer as the five annotations they behave exactly like.
    "filled": "Filled if absent",
    "optional": "Optional",
    "uncategorised": "Optional",
    "derived": "Not yours — MAFigate adds it",
}


def _classify_column_source(column_name):
    """Which fact about the user's file this column stands for. A token, not prose.

    Kept separate from :data:`_COLUMN_SOURCE_PROSE` so the two facts that share the word
    "Optional" stay two facts. Returns one of that mapping's keys.

    The order is the order of the questions, and each step is a different kind of fact:

    1. **Does any MAF supply it?** The derived columns first, because no membership of a list
       describing MAFs can be right about a column no MAF has.
    2. **Is the file refused without it?** ``REQUIRED_COLUMNS["core"]`` — the one category
       ``validate_required_columns`` actually refuses on. ``Variant_Classification`` is both
       core and a filter input and is answered here, which is the true order: the file is
       refused before there is anything to fill.
    3. **Does the filter fill it?** Read out of the vendored source rather than from a
       category (issue #127), so this answer cannot drift from what the filter does. The union
       of both arms, because this table is the union of both arms.
    4. Otherwise the file may or may not carry it, and ``OPTIONAL_COLUMNS`` membership only
       decides *which token* says so — see :data:`_COLUMN_SOURCE_PROSE` on why both spell the
       same sentence, and ``test_no_filter_input_is_called_optional`` on what the distinction
       buys.

    The import is function-local, following :func:`pipeline_output_columns` right below:
    ``filters/absent_columns.py`` parses the vendored source at import time and this module is
    the leaf every page imports, so the cost is paid by the one caller that needs it. There is
    no cycle to dodge — ``filters/absent_columns.py`` imports nothing from ``config/`` — but
    the layering is worth not inverting for a function the Help page calls once per render.
    """
    from filters.absent_columns import REQUIRED_INPUTS

    if column_name in COLUMNS_MAFIGATE_ADDS:
        return "derived"

    if column_name in REQUIRED_COLUMNS["core"]:
        return "core"

    if any(column_name in columns for columns in REQUIRED_INPUTS.values()):
        return "filled"

    for columns in OPTIONAL_COLUMNS.values():
        if column_name in columns:
            return "optional"

    return "uncategorised"


def column_source_status(column_name):
    """What the column reference says under :data:`COLUMN_SOURCE_HEADER` for this column.

    Named ``get_column_requirement_status`` until issue #124, when the column above it
    stopped being called "Required": a function of that name returning "Not yours —
    MAFigate adds it" carries the same mismatch one layer down. Its only caller is
    :func:`create_column_info_table`.
    """
    return _COLUMN_SOURCE_PROSE[_classify_column_source(column_name)]


def pipeline_output_columns(
    sample_type="somatic", skip_civic=False, available_columns=None
):
    """The columns the *pipeline* would emit, in the pipeline's order.

    Delegates to the vendored :func:`vendor.pipeline_filters.compute_keep`, which is
    ``bin/filter_variants.py:main()``'s own statements. Nothing here transcribes the
    column list; drift between this and the pipeline is caught by
    ``tests/test_vendor_drift.py``.

    Args:
        sample_type: ``"somatic"`` or ``"germline"`` — the arm.
        skip_civic: whether the CIViC annotation step was skipped.
        available_columns: the columns actually present, which decide the CIViC and
            gnomAD branches. ``None`` means "assume everything is present".

    Returns:
        The ordered column list. No entry is filtered for presence here — that is
        :func:`resolve_visible_columns`'s job, so it can report what went missing.
    """
    from types import SimpleNamespace

    from vendor.pipeline_filters import compute_keep

    columns = (
        _AllColumnsPresent()
        if available_columns is None
        else _ColumnsShim(available_columns)
    )
    args = SimpleNamespace(sample_type=sample_type, skip_civic=skip_civic)
    return compute_keep(args, SimpleNamespace(columns=columns))


def resolve_visible_columns(
    sample_type="somatic", skip_civic=False, available_columns=None
):
    """The ordered columns the app shows and exports. The app's single answer.

    The pipeline's list first, less :data:`PIPELINE_COLUMNS_THE_APP_REPLACES`, then
    :data:`APP_EXTRA_COLUMNS` — so the return value opens with
    :func:`pipeline_output_columns` as an exact prefix once the replaced columns are
    taken out of both.

    That subtraction is the one place the app shows less than the pipeline emits, it is
    a list of names rather than a rule, and the argument for each name is written beside
    it. Nothing here decides that a column is not worth showing; the constant records
    that the app already answers the same question in a column of its own.

    Args:
        sample_type: ``"somatic"`` or ``"germline"`` — the arm.
        skip_civic: whether the CIViC annotation step was skipped.
        available_columns: the columns actually present in the frame this list will
            index. ``None`` means "assume everything is present", which is what the
            help pages want: they describe the schema with no MAF in hand.

    Returns:
        Columns present in ``available_columns``, in order. Callers index a frame with
        this, so every name in it is safe to select.

    Warns:
        MissingColumnsWarning: a column the pipeline's output would carry is absent
            from ``available_columns``. Absent *app extras* are dropped silently —
            they are derived at display time, so most callers legitimately run before
            they exist, and warning about them would drown the real signal.

            The wording puts the absence on the **file** — *this MAF does not carry*,
            not *left out of the table* — because since issue #93 the app says this
            with the grid's *Show all columns* ticked as well, and in that state the
            table is showing every column the file has. Nothing is being left out
            then; the columns were never there. Both halves have to stay true in both
            states, which the old phrasing did not: it read as MAFigate doing the
            omitting.
    """
    # The subtraction happens here and not in `pipeline_output_columns`, which must go on
    # answering what the *pipeline* emits. It is applied before the missing-column check
    # below on purpose: a column the app has decided not to show is not a column the file
    # is short of, and warning "this MAF does not carry variantalker_naive … it is in
    # neither the table nor the export" would name an absence the app was going to arrange
    # anyway. Measured: that sentence names it today on 74 of the 176 annotated real MAFs.
    pipeline = [
        col
        for col in pipeline_output_columns(sample_type, skip_civic, available_columns)
        if col not in PIPELINE_COLUMNS_THE_APP_REPLACES
    ]
    ordered = pipeline + APP_EXTRA_COLUMNS

    if available_columns is None:
        return ordered

    present = set(available_columns)
    missing = [col for col in pipeline if col not in present]
    if missing:
        warnings.warn(
            f"This MAF does not carry {len(missing)} column(s) the {sample_type} "
            f"pipeline emits: {', '.join(missing)}. They are in neither the table nor "
            "the export; every other column is unaffected.",
            MissingColumnsWarning,
            stacklevel=2,
        )

    return [col for col in ordered if col in present]


def create_column_info_table():
    """Create a comprehensive column information table with descriptions"""
    try:
        import pandas as pd
    except ImportError:
        raise ImportError(
            "pandas is required for create_column_info_table(). Install with: pip install pandas"
        )

    # Get all columns from the resolver — the same list the grid and the export use —
    # plus the columns the app deliberately keeps out of that view.
    #
    # The glossary is wider than the default table on purpose, and issue #117 is why:
    # a column the app stops *showing* is still in the user's file and still one tick of
    # "Show all columns" away, so this table is where they can look it up. It is the only
    # thing the app now says about `variantalker_naive`, which is what makes leaving it
    # out of the view honest rather than silent. The arm membership is read from the
    # pipeline's own list for exactly these columns — reading it from the resolver would
    # report them as absent from both arms and the branch below would then label them
    # "Germline", which is a fact about the if/elif chain rather than about the column.
    all_somatic_cols = resolve_visible_columns("somatic")
    all_germline_cols = resolve_visible_columns("germline")
    for arm, resolved in (("somatic", all_somatic_cols), ("germline", all_germline_cols)):
        emitted = pipeline_output_columns(arm)
        resolved.extend(
            col for col in PIPELINE_COLUMNS_THE_APP_REPLACES if col in emitted
        )

    # Combine and get unique columns
    all_columns = list(set(all_somatic_cols + all_germline_cols))

    # Create table data
    table_data = []
    for col in sorted(all_columns):
        description = get_column_description(col)

        # Determine which sample types use this column
        in_somatic = col in all_somatic_cols
        in_germline = col in all_germline_cols

        if in_somatic and in_germline:
            sample_types = "Both"
        elif in_somatic:
            sample_types = "Somatic"
        else:
            sample_types = "Germline"

        # Categorize column type
        if col in [
            "Tumor_Sample_Barcode",
            "Matched_Norm_Sample_Barcode",
            # Derived from the two above it, and falling through to "Other" until issue
            # #120 described it: a user filtering this table by Sample Information got the
            # barcodes and not the column the grid actually shows them.
            "Sample_Name",
            "project_id",
        ]:
            category = "Sample Information"
        elif col in ["Hugo_Symbol", "Annotation_Transcript"]:
            category = "Gene Information"
        elif col in ["Chromosome", "Start_Position", "End_Position"]:
            category = "Genomic Coordinates"
        elif col in [
            "Variant_Classification",
            "Variant_Type",
            "Reference_Allele",
            "Tumor_Seq_Allele1",
            "Tumor_Seq_Allele2",
        ]:
            category = "Variant Information"
        elif col in [
            "AAChange.refGene",
            "cDNA_Change",
            "Codon_Change",
            "Protein_Change",
            "Transcript_Exon",
        ]:
            category = "Functional Annotation"
        elif col in [
            "tumor_f",
            "DP",
            "t_alt_count",
            "t_ref_count",
            "n_alt_count",
            "n_ref_count",
        ]:
            category = "Sequencing Metrics"
        elif col in [
            "ClinVar_VCF_CLNSIG",
            "CancerVar",
            "ESCAT",
            "ESCAT_TISSUE",
            "ESCAT_CANCER",
            # The two verdicts of issue #117. They belong beside the sources they read,
            # not in "Other", which is where the chain's fallthrough had been putting
            # `Clinical_Summary` — the column the table pins first.
            "Clinical_Summary",
            "variantalker_naive",
            # And its twin, pinned immediately to its right, which #117 left in "Other"
            # because it was only moving the columns it was describing. Same argument, same
            # pair: these two summarise the same sources for the same variant and are read
            # side by side, so a category filter that returns one and not the other splits
            # them where the screen puts them together (issue #120).
            "Pathogenicity_Overview",
        ]:
            category = "Clinical Significance"
        elif col.startswith("CIViC"):
            category = "CIViC Database"
        # `MAX_AF` left this branch with its description (issue #126). It has never reached
        # this table -- the table is the resolver's set plus the app's extras, and `MAX_AF` is
        # in neither -- so the entry was classifying a column that cannot arrive.
        elif col in [
            "Freq_ExAC_ALL",
            "Freq_esp6500siv2_all",
            "Freq_1000g2015aug_all",
        ] or col.startswith("gnomAD"):
            category = "Population Frequency"
        elif col in ["InterVar", "RENOVO_Class", "RENOVO_pls"]:
            category = "Germline Assessment"
        else:
            category = "Other"

        table_data.append(
            {
                "Column Name": col,
                "Category": category,
                "Sample Types": sample_types,
                # The frame's key is the header the reader sees. See COLUMN_SOURCE_HEADER
                # for why both ends read the constant.
                COLUMN_SOURCE_HEADER: column_source_status(col),
                "Description": description,
            }
        )

    return pd.DataFrame(table_data)


#: The order columns are shown in when the user asks to see *all* of them.
#:
#: Not the default view — that is :func:`resolve_visible_columns`, which selects and
#: orders what the grid opens with. This list reaches the user in exactly two places:
#: the grid's "Show all columns" checkbox, and the order of the "Add columns" dropdown.
#: So it is a presentation order over a frame that already has its columns, which is why
#: it can sit here beside the resolver without competing with it.
#:
#: That two-places claim was false when it was written, and issue #90 is what made it
#: true: the Summary tab's ``Column Analysis`` block was a third reader, and it used this
#: list to answer *which columns is my MAF missing?* — a question only the resolver can
#: answer, because only the resolver knows the arm and asks ``compute_keep`` what it emits
#: for the columns actually present. Read as a membership set this list says a complete
#: germline MAF is missing nine columns. It is an order; do not measure absence with it.
#:
#: Named ``KEY_COLUMNS`` until issue #76, in ``components/ui_components.py``. Renamed on
#: the way here because the old name was wrong twice over: beside
#: :data:`REQUIRED_COLUMNS`, :data:`OPTIONAL_COLUMNS` and :data:`APP_EXTRA_COLUMNS` it
#: reads as a fourth *membership* set when it is an *ordering*, and ``tests/parity``
#: already has a ``KEY_COLUMNS`` meaning something else entirely — the four columns that
#: identify a variant when joining the app's output against the pipeline's.
PRIORITY_COLUMNS = [
    "Clinical_Summary",  # New summary column (first)
    "Tumor_Sample_Barcode",
    "Matched_Norm_Sample_Barcode",
    "project_id",
    "Hugo_Symbol",
    "Annotation_Transcript",
    "Chromosome",
    "Start_Position",
    "End_Position",
    "Variant_Classification",
    "Variant_Type",
    "Reference_Allele",
    "Tumor_Seq_Allele1",
    "Tumor_Seq_Allele2",
    "AAChange.refGene",
    "cDNA_Change",
    "Codon_Change",
    "Protein_Change",
    "Transcript_Exon",
    "tumor_f",
    "DP",
    "t_alt_count",
    "t_ref_count",
    "n_alt_count",
    "n_ref_count",
    "ClinVar_VCF_CLNSIG",
    "CancerVar",
    "ESCAT",
    "ESCAT_TISSUE",
    "ESCAT_CANCER",
    "CIViC_Evidence_Level",
    "CIViC_Evidence_Rating",
    "CIViC_Entity_Disease",
    "CIViC_Variant_URL",
    "CIViC_Entity_URL",
    "CIViC_Entity_Status",
    "am_class",
    "am_pathogenicity",
    "Otherinfo",
    "tumor_tissue",
    "cosmic",
    "Freq_ExAC_ALL",
    "Freq_esp6500siv2_all",
    "Freq_1000g2015aug_all",
    "filter",
]


def prioritize_columns(data: pd.DataFrame) -> pd.DataFrame:
    """Reorder columns to show key columns first, followed by others."""

    available_key_cols = [col for col in PRIORITY_COLUMNS if col in data.columns]
    other_cols = [col for col in data.columns if col not in PRIORITY_COLUMNS]

    # Create new column order: key columns first, then others
    new_column_order = available_key_cols + sorted(other_cols)

    return data[new_column_order]


#: Which characters R's ``make.names`` and pandas' own mangling replace with a dot. Every one of
#: them is legal in an ANNOVAR/dbNSFP column name and illegal in an R identifier, which is how a
#: column arrives spelled two ways.
#:
#: A dot is **not** in this set, because a dot is legal in an R identifier: ``make.names`` leaves
#: ``AAChange.refGene``, ``Func.refGene`` and ``ExonicFunc.refGene`` exactly as they are. That is
#: what makes the character class a usable rule rather than a rough filter — it is the set of
#: characters whose presence in a name means the name *can* arrive spelled two ways, and whose
#: absence means it cannot. Issue #212 measured it: of the **58** column names the variant panel
#: reads, exactly **3** contain one of these characters.
_MANGLED_TO_DOT = re.compile(r"[^A-Za-z0-9_.]")


def is_mangling_exposed(name: str) -> bool:
    """Whether ``make.names`` would rewrite this column name (issue #212).

    Args:
        name: a column name as the annotator documents it.

    Returns:
        bool: True when the name carries a character illegal in an R identifier, so a file that
        passed through ``make.names`` spells it differently and an exact-name lookup misses it.

    This is the rule ``tests/test_column_spelling.py`` holds the variant panel to: a name this
    returns True for must reach the row through :func:`spelled_in` rather than by equality. It is
    a statement about the character class, so a column added tomorrow is covered without anyone
    re-measuring the corpus.
    """
    stripped = str(name).strip()
    return _normalised_name(stripped) != stripped


def spelled_in(columns, name: str) -> "str | None":
    """The name of a column as *this file* spells it, or ``None`` if it has no such column.

    Args:
        columns: any iterable of column names — a ``DataFrame.columns``, a ``Series.index``, or a
            plain list.
        name: the canonical spelling, as the annotator documents it.

    Returns:
        str: the file's own spelling. The canonical name where the file uses it, and otherwise a
        dot-mangled equivalent — ``'GERP.._RS'`` for ``'GERP++_RS'`` — where exactly one column
        mangles to the same shape. ``None`` when the file has neither.

    **Why this exists.** Somewhere between the pipeline and the MAF a step passes the header
    through R's ``make.names`` (or pandas' equivalent), which replaces every character illegal in
    an identifier with a dot. Issue #189 met it on the evidence column —
    ``'CancerVar..CancerVar.and.Evidence'`` on 3 of 124 files — and worked around it by matching on
    a substring. Issue #190 met it again on ``GERP++_RS``, spelled ``'GERP.._RS'`` on **2 of 167**
    real MAFs, and a substring match is no use there: ``GERP`` alone would also match
    ``GERP++_NR``, a different score.

    So the comparison is *made* on the mangled shape rather than searched for in it: both sides are
    reduced to the same normal form, and a name matches when the normal forms are equal. The
    canonical spelling is preferred whenever it is present, so a file carrying both spellings — none
    does — gets the documented one rather than an arbitrary one.

    **Surrounding whitespace is stripped before normalising, and that is what lets this serve every
    read on the panel.** Issue #212 measured the resolver as issue #190 left it and found it could
    not: the evidence columns arrive space-padded — ``' CancerVar: CancerVar and Evidence '`` — and
    a space is a character ``make.names`` rewrites, so the padding normalised to leading and
    trailing dots and the comparison failed on **93 of the 96** files carrying the CancerVar
    evidence column and **56 of 56** carrying InterVar's. Wiring the resolver in at that read
    without the strip would have resolved the 3 dot-mangled files and lost the 93 padded ones.
    With it, the resolver and the substring match it replaced agree on **all 152** files that carry
    either column, with **0** disagreements.

    The strip is applied to the *comparison* only. What comes back is the file's own spelling,
    padding included, because that string is the key the caller indexes the row with.

    **A tie is refused.** If two columns mangle to the same shape and neither is the canonical
    spelling, there is no way to tell which the caller meant, and returning either would be a guess
    about a number a clinician reads. ``None`` is the honest answer, and the caller's
    already-required *absent* path handles it. No real MAF holds such a tie: across 167 files and
    727 distinct column names, **0** normalised shapes carry more than one spelling within a single
    file.
    """
    names = list(columns)
    if name in names:
        return name
    target = _normalised_name(name)
    matches = [column for column in names if _normalised_name(column) == target]
    if len(matches) == 1:
        return matches[0]
    return None


def _normalised_name(name) -> str:
    """A column name reduced to the shape ``make.names`` would leave it in.

    Stripped first, so a padded spelling and its bare equivalent share a normal form — see
    :func:`spelled_in` for the measurement that made the strip necessary.
    """
    return _MANGLED_TO_DOT.sub(".", str(name).strip())
