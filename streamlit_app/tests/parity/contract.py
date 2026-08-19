"""The parameter contract, its translation to both sides, and the case matrix.

One dict of *contract* values (issue #11) is translated twice: into
``bin/filter_variants.py`` command-line arguments, and into the ``params`` dict the
app's ``filters.variant_filters.apply_filters`` consumes. Keeping the translation in one
place is what makes "equivalent parameters" a checkable claim rather than a
judgement call — every divergence the harness reports is then a divergence in
*filter logic*, not in how the harness chose to configure the two sides.

Since issue #33 the app-side translation is nearly the identity: the app's parameter names
*are* the pipeline's long options, with the pipeline's ``skip_*`` polarity and the
pipeline's meanings. What remains in :func:`app_params` is three shape differences, each
documented there, and none of them a divergence.

Contract values are ``nextflow.config``'s ``params`` block, not argparse defaults:
``modules/local/annotation/small_variants/main.nf:316`` always passes every
argument, so the argparse defaults are unreachable in production. See issue #11.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from pathlib import Path

FIXTURE_DIR = Path(__file__).resolve().parent.parent / "fixtures" / "parity"

# ---------------------------------------------------------------------------
# The contract
# ---------------------------------------------------------------------------

# nextflow.config lines 79-96, with issue #11's single deliberate departure:
# filter_clinvar carries underscored "Likely_pathogenic" rather than the config's
# spaced form, which matches 0 of 181,566 reference rows.
CONTRACT = {
    "min_depth": 50,                                    # filter_min_depth
    "vaf_threshold": 0.01,                              # filter_vaf_threshold
    "vaf_threshold_germline": 0.2,                      # filter_vaf_threshold_germline
    "filter_variant_classification": ["Silent", "IGR", "RNA"],  # filter_var_classification
    "filter_cancervar": ["Tier_II_potential", "Tier_I_strong"],
    "filter_civic": ["A", "B", "C"],                    # filter_civic_evidence_level
    "filter_clinvar": ["Pathogenic", "Likely_pathogenic"],
    "filter_escat": ["IA", "IB", "IC", "IIA", "IIB"],
    "filter_intervar": ["Pathogenic", "Likely pathogenic"],
    "filter_renovo": ["LP Pathogenic", "IP Pathogenic", "HP Pathogenic"],
    "skip_civic": False,                                # skip_civic
    "skip_pathogenic": False,                           # skip_pathogenic_retention
    "genes": None,                                      # filter_genes_{somatic,germline} = null
    # App-only extra. 1.0 is the neutral value issue #16 requires PIPELINE_PARAMS
    # to carry; the pipeline has no counterpart at all.
    "max_freq_population": 1.0,
}

# The app's hardcoded classification vocabulary, mirrored here rather than imported
# so the harness records the *baseline* vocabulary even after issue #12's fix moves
# or deletes it. Kept in step by test_parity.py::test_app_vocabulary_unchanged.
APP_VARIANT_CLASSIFICATIONS = [
    "Missense_Mutation",
    "Nonsense_Mutation",
    "Silent",
    "Frame_Shift_Del",
    "Frame_Shift_Ins",
    "In_Frame_Del",
    "In_Frame_Ins",
    "Nonstop_Mutation",
    "Translation_Start_Site",
    "Splice_Site",
    "Splice_Region",
    "5'UTR",
    "3'UTR",
    "5'Flank",
    "3'Flank",
    "Intron",
    "IGR",
    "RNA",
]

# The four-column variant identity the two outputs are joined on.
KEY_COLUMNS = ["Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2"]


# ---------------------------------------------------------------------------
# Cases
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Case:
    """One (fixture, arm, parameter set) the two sides are compared under."""

    name: str
    fixture: str
    arm: str
    overrides: dict = field(default_factory=dict)
    # How the contract's *exclude* list is handed to the app, whose parameter of
    # the same name is an *include* list (divergence #1):
    #   "literal"    -- the same list value, so the app keeps exactly what the
    #                   pipeline drops. Documents the divergence at full size.
    #   "complement" -- the app gets APP_VARIANT_CLASSIFICATIONS minus the list,
    #                   the closest expressible equivalent. Lets the other
    #                   divergences be seen instead of being swamped by this one.
    vc_mode: str = "complement"
    # How the gene list reaches the app, where the pipeline always gets a file path. The
    # separator is the variable because the separator is where the bug lived:
    #   "tokens"     -- pre-split, one symbol per list element. What a saved parameter
    #                   file holds, and what the harness has always used.
    #   "lines"      -- the gene file's raw text, one symbol per line.
    #   "spaces"     -- the symbols separated by single spaces.
    #   "semicolons" -- the symbols separated by "; ".
    # The last three are all one string the app must tokenise. See `app_genes` for which
    # of them can actually witness the bug, and why "lines" cannot.
    genes_as: str = "tokens"
    # Set when neither side is expected to produce a verdict. Records the *pipeline's*
    # exception only; what the app does with the same file is measured, not declared here.
    # Documentation: nothing reads this field. See `harness.ParityResult.errors_in_parity`.
    expect_error: str | None = None
    notes: str = ""

    @property
    def contract(self) -> dict:
        merged = dict(CONTRACT)
        merged.update(self.overrides)
        return merged

    @property
    def maf_path(self) -> Path:
        return FIXTURE_DIR / self.fixture

    @property
    def genes_path(self) -> Path | None:
        genes = self.contract["genes"]
        return None if genes is None else FIXTURE_DIR / genes


_EMPTY_GUIDELINES_SOMATIC = {
    "filter_cancervar": [],
    "filter_civic": [],
    "filter_clinvar": [],
    "filter_escat": [],
}

_EMPTY_GUIDELINES_GERMLINE = {
    "filter_intervar": [],
    "filter_renovo": [],
    "filter_clinvar": [],
}


CASES: list[Case] = [
    # -- Contract defaults, both readings of divergence #1 -------------------
    Case(
        "somatic_defaults_literal",
        "somatic_reference.maf",
        "somatic",
        vc_mode="literal",
        notes="divergence #1 at full size: the same list, opposite semantics",
    ),
    Case(
        "somatic_defaults",
        "somatic_reference.maf",
        "somatic",
        notes="the headline case",
    ),
    Case(
        "germline_defaults_literal",
        "germline_reference.maf",
        "germline",
        vc_mode="literal",
        notes="divergence #1 at full size, germline",
    ),
    Case(
        "germline_defaults",
        "germline_reference.maf",
        "germline",
        notes="the headline case, germline",
    ),
    # -- The constructed shapes --------------------------------------------
    Case(
        "somatic_synthetic",
        "somatic_synthetic.maf",
        "somatic",
        notes="one-field edits: depth, VAF, NaN, ';' CLNSIG, bare ESCAT",
    ),
    Case(
        "germline_synthetic",
        "germline_synthetic.maf",
        "germline",
        notes="ditto, germline",
    ),
    # -- CIViC, the only fixture that separates skip_civic from column absence
    Case(
        "civic_present",
        "somatic_civic.maf",
        "somatic",
        notes="divergence #12: substring vs isin against list-repr values",
    ),
    Case(
        "civic_skipped",
        "somatic_civic.maf",
        "somatic",
        overrides={"skip_civic": True},
        notes="the app has no skip_civic at all",
    ),
    Case(
        "gnomad_genome_present",
        "somatic_gnomad_genome.maf",
        "somatic",
        notes="column-set probe: compute_keep's gnomAD_genome branch",
    ),
    # -- skip_pathogenic, which is what gives a depth/VAF sweep any power ----
    # Issue #18: at contract defaults somatic PASS is filter_patho union
    # (ESCAT-in-keep & common & VAF & genes), so the unconditional pathogenic
    # rescue absorbs nearly every depth or VAF edit.
    Case(
        "somatic_skip_patho",
        "somatic_reference.maf",
        "somatic",
        overrides={"skip_pathogenic": True},
    ),
    Case(
        "germline_skip_patho",
        "germline_reference.maf",
        "germline",
        overrides={"skip_pathogenic": True},
    ),
    # -- Gene lists (divergence #9) -----------------------------------------
    Case(
        "somatic_genes",
        "somatic_reference.maf",
        "somatic",
        overrides={"genes": "genes_somatic.txt"},
    ),
    Case(
        "somatic_genes_mixed_case",
        "somatic_reference.maf",
        "somatic",
        overrides={"genes": "genes_somatic_mixed_case.txt"},
        notes="divergence #9: same genes, case mangled",
    ),
    Case(
        "germline_genes",
        "germline_reference.maf",
        "germline",
        overrides={"genes": "genes_germline.txt"},
    ),
    # The gene cases above have almost no power on their own: the unconditional
    # pathogenic rescue re-admits the rows a gene list excludes, so a case-sensitivity
    # divergence can be present and still leave the PASS sets equal. Measured, not
    # assumed — a trial fix to the pathogenic clause brought `somatic_genes_mixed_case`
    # to parity while the app was still comparing gene symbols case-sensitively. These
    # variants add `skip_pathogenic` so divergence #9 actually reaches the verdict.
    Case(
        "somatic_genes_skip_patho",
        "somatic_reference.maf",
        "somatic",
        overrides={"genes": "genes_somatic.txt", "skip_pathogenic": True},
    ),
    Case(
        "somatic_genes_mixed_case_skip_patho",
        "somatic_reference.maf",
        "somatic",
        overrides={"genes": "genes_somatic_mixed_case.txt", "skip_pathogenic": True},
        notes="divergence #9 with the pathogenic rescue out of the way",
    ),
    Case(
        "germline_genes_skip_patho",
        "germline_reference.maf",
        "germline",
        overrides={"genes": "genes_germline.txt", "skip_pathogenic": True},
    ),
    # -- Gene lists as a user pastes them (issue #34) ------------------------
    # The cases above hand the app a pre-split list, so they cannot see the tokeniser at
    # all. These hand it one string and make it do the splitting, while the pipeline reads
    # the same symbols out of the file.
    #
    # `skip_pathogenic` on every one of them, and that is not incidental: without it the
    # unconditional pathogenic rescue re-admits exactly the rows a gene filter excludes,
    # so these cases sit at apparent parity while the bug is fully present. Measured, on
    # the variants issue #34's "Watch out" names.
    #
    # `app_genes` records which separators can witness the bug and which cannot — the
    # newline shape cannot, for a reason worth reading before adding more cases here.
    Case(
        "somatic_genes_pasted_skip_patho",
        "somatic_reference.maf",
        "somatic",
        overrides={"genes": "genes_somatic.txt", "skip_pathogenic": True},
        genes_as="lines",
        notes="the one-per-line paste; already at parity via a file round-trip",
    ),
    Case(
        "somatic_genes_spaced_skip_patho",
        "somatic_reference.maf",
        "somatic",
        overrides={"genes": "genes_somatic.txt", "skip_pathogenic": True},
        genes_as="spaces",
        notes="a spreadsheet row: the separator the old comma-only split lost",
    ),
    Case(
        "somatic_genes_mixed_case_semicolons_skip_patho",
        "somatic_reference.maf",
        "somatic",
        overrides={
            "genes": "genes_somatic_mixed_case.txt",
            "skip_pathogenic": True,
        },
        genes_as="semicolons",
        notes="both bugs at once: an unhandled separator *and* case-mangled symbols",
    ),
    Case(
        "germline_genes_spaced_skip_patho",
        "germline_reference.maf",
        "germline",
        overrides={"genes": "genes_germline.txt", "skip_pathogenic": True},
        genes_as="spaces",
    ),
    # -- Depth sweep --------------------------------------------------------
    Case(
        "somatic_depth_0",
        "somatic_reference.maf",
        "somatic",
        overrides={"min_depth": 0, "skip_pathogenic": True},
        notes="min_depth 0: the app's `if min_depth > 0` guard skips the filter "
        "entirely, the pipeline still evaluates `sum >= 0` (False for NaN)",
    ),
    Case(
        "somatic_synthetic_depth_0",
        "somatic_synthetic.maf",
        "somatic",
        overrides={"min_depth": 0, "skip_pathogenic": True},
        notes="min_depth 0 is not the unfiltered case: the pipeline still evaluates "
        "`t_alt + t_ref >= 0`, which is False for NaN, and drops 3 of these 19 rows. "
        "The reference fixtures cannot witness it -- they have no missing numerics.",
    ),
    Case(
        "somatic_depth_500",
        "somatic_reference.maf",
        "somatic",
        overrides={"min_depth": 500, "skip_pathogenic": True},
    ),
    Case(
        "germline_depth_500",
        "germline_reference.maf",
        "germline",
        overrides={"min_depth": 500, "skip_pathogenic": True},
    ),
    # -- VAF sweep ----------------------------------------------------------
    Case(
        "somatic_vaf_005",
        "somatic_reference.maf",
        "somatic",
        overrides={"vaf_threshold": 0.05, "skip_pathogenic": True},
        notes="3 reference rows sit at exactly 0.05 (issue #18)",
    ),
    Case(
        "somatic_vaf_02",
        "somatic_reference.maf",
        "somatic",
        overrides={"vaf_threshold": 0.2, "skip_pathogenic": True},
        notes="4 reference rows sit at exactly 0.2",
    ),
    Case(
        "germline_vaf_02",
        "germline_reference.maf",
        "germline",
        overrides={"skip_pathogenic": True},
        notes="6 germline rows sit at exactly the 0.2 contract threshold",
    ),
    Case(
        "germline_vaf_005",
        "germline_reference.maf",
        "germline",
        overrides={"vaf_threshold_germline": 0.05, "skip_pathogenic": True},
    ),
    # -- Empty keep-lists (issue #13's now-reachable all-empty state) --------
    Case(
        "somatic_empty_guidelines",
        "somatic_reference.maf",
        "somatic",
        overrides=_EMPTY_GUIDELINES_SOMATIC,
        notes="pipeline-expressible; bottoms out at filter_patho's rescue, not zero",
    ),
    Case(
        "germline_empty_guidelines",
        "germline_reference.maf",
        "germline",
        overrides=_EMPTY_GUIDELINES_GERMLINE,
    ),
    Case(
        "somatic_empty_guidelines_skip_patho",
        "somatic_reference.maf",
        "somatic",
        overrides={**_EMPTY_GUIDELINES_SOMATIC, "skip_pathogenic": True},
        notes="nothing rescues anything: the pipeline's floor is 0",
    ),
    # -- Error parity -------------------------------------------------------
    Case(
        "dot_numeric",
        "somatic_dot_numeric.maf",
        "somatic",
        expect_error="TypeError",
        notes=(
            "'.' makes the column object dtype; the pipeline raises (issue #18) and "
            "the app refuses, naming t_alt_count, t_ref_count and tumor_f (issue #38). "
            "The fourth '.' is in DP, which the pipeline never reads, and is not named."
        ),
    ),
]

CASES_BY_NAME = {case.name: case for case in CASES}


# ---------------------------------------------------------------------------
# Translation
# ---------------------------------------------------------------------------


def pipeline_args(case: Case) -> list[str]:
    """Contract values as ``bin/filter_variants.py`` command-line arguments.

    ``--config``, ``--funcotator``, ``--technology``, ``--annovar_version`` and
    ``--genome_build`` are required by argparse and never read (issue #18), so they
    are passed as fixed placeholders and are deliberately not part of the contract.
    """
    c = case.contract
    args = [
        "--maf", str(case.maf_path),
        "--sample_type", case.arm,
        "--min_depth", str(c["min_depth"]),
        "--vaf_threshold", str(c["vaf_threshold"]),
        "--vaf_threshold_germline", str(c["vaf_threshold_germline"]),
        "--filter_variant_classification", ",".join(c["filter_variant_classification"]),
        "--filter_cancervar", ",".join(c["filter_cancervar"]),
        "--filter_civic", ",".join(c["filter_civic"]),
        "--filter_clinvar", ",".join(c["filter_clinvar"]),
        "--filter_escat", ",".join(c["filter_escat"]),
        "--filter_intervar", ",".join(c["filter_intervar"]),
        "--filter_renovo", ",".join(c["filter_renovo"]),
        # argparse-required, never read:
        "--config", "unused",
        "--funcotator", "unused",
        "--technology", "iontorrent",
        "--annovar_version", "unused",
        "--genome_build", "hg19",
    ]
    if c["skip_civic"]:
        args.append("--skip_civic")
    if c["skip_pathogenic"]:
        args.append("--skip_pathogenic")
    if c["genes"] is not None:
        flag = "--filter_genes_somatic" if case.arm == "somatic" else "--filter_genes_germline"
        args += [flag, str(case.genes_path)]
    return args


def read_gene_tokens(path: Path) -> list[str]:
    """Gene symbols from a pipeline-format gene file (one per line)."""
    return [line.strip() for line in path.read_text().splitlines() if line.strip()]


#: How each ``Case.genes_as`` paste mode joins the file's symbols back into one string.
_PASTE_SEPARATORS = {"lines": "\n", "spaces": " ", "semicolons": "; "}


def app_genes(case: Case):
    """The gene list as the app is handed it — a list of symbols, or one pasted string.

    The pipeline is always handed a *path*, so this is the one place the two sides are
    configured from genuinely different shapes, and ``Case.genes_as`` records which shape
    each case documents. Both must arrive at the same symbols; that is the claim.

    Which modes can witness the bug, measured
    -----------------------------------------
    ``lines`` **cannot**, and the reason is worth writing down because it is not obvious
    and it would otherwise be mistaken for coverage. The pre-#34 reader split on commas
    only, so a newline-delimited paste became a single token — but the gene adapter writes
    tokens to a file verbatim, and the newlines *inside* that token became line breaks in
    the file, which ``pd.read_csv`` then split back into the right symbols. Measured
    against the recorded pipeline PASS set with the old reader restored: 4 of 4 somatic and
    3 of 3 germline, identical. The headline paste shape was therefore already at parity
    after issue #33 — accidentally, via a file round-trip — and ``lines`` is kept as the
    case that pins that, not as the case that proves the fix.

    ``spaces`` and ``semicolons`` **do**: neither separator survives being written into a
    one-line file, so the old reader restricts to a symbol that matches no row and the
    somatic criteria path falls to 0. Those are the cases with power, and both are shapes a
    user really produces — a row copied out of a spreadsheet, or a semicolon-delimited
    export.

    Deliberately no mode for an *invalid* gene list. The app drops letterless tokens and
    applies no gene filter where the pipeline crashes on the same file (``AttributeError``
    out of ``.str.upper()`` on a numeric column), so such a case could only ever be a
    permanent xfail recording a deviation the app wants. It is asserted where it belongs
    instead — ``tests/test_gene_lists.py``, against the vendored function directly.
    """
    if case.contract["genes"] is None:
        return []
    tokens = read_gene_tokens(case.genes_path)
    if case.genes_as in _PASTE_SEPARATORS:
        return _PASTE_SEPARATORS[case.genes_as].join(tokens)
    return tokens


def app_variant_classifications(case: Case) -> list[str]:
    """The classification list as the app is handed it — now an *exclude* list.

    Before issue #33 the app's parameter of this name was an *include* list, and
    ``vc_mode`` chose between two ways of translating the contract's exclude list into
    it: ``literal`` (the same values, so the app kept exactly what the pipeline dropped)
    and ``complement`` (the app's vocabulary minus those values, the closest expressible
    equivalent). Routing the decision through the vendored ``common_filters`` made the
    app's parameter mean what the pipeline's means, so there is nothing left to translate
    and both modes now yield the same value.

    ``vc_mode`` is kept on :class:`Case` as the record of which reading each case was
    built to document — the two ``*_literal`` cases are now duplicates of their
    siblings, deliberately retained so the case set and its baseline records stay
    comparable across the fix.
    """
    return list(case.contract["filter_variant_classification"])


def app_params(case: Case) -> dict:
    """Contract values as the ``params`` dict :func:`filters.apply_filters` consumes.

    The translation is now almost the identity, which is the point of issue #33: the app
    takes the pipeline's own parameter names, with the pipeline's ``skip_*`` polarity and
    the pipeline's meanings. Three shape differences remain, none of them a divergence:

    * ``genes`` — the pipeline takes a *file path*, the app a list of symbols or the raw
      text of a paste (see :func:`app_genes` and ``Case.genes_as``). The app tokenises
      whichever it gets and its gene adapter writes the symbols back out to a file for the
      vendored clause to read, so both sides end up reading the same symbols.
    * ``vaf_threshold`` — one key in the app, whose value is the arm's. The pipeline
      takes both thresholds on every command line and picks one by arm.
    * ``filter_escat`` — somatic only. ``germline_filters`` has no ESCAT argument, so the
      germline arm is handed no such key; a germline ESCAT clause was the largest
      divergence on real data.

    Kept in step with ``config/pipeline_params.py`` — the app's own statement of the same
    contract — by ``tests/test_param_contract.py``.
    """
    c = case.contract
    params = {
        "sample_type": case.arm,
        "min_depth": c["min_depth"],
        "filter_variant_classification": app_variant_classifications(case),
        "filter_clinvar": list(c["filter_clinvar"]),
        "max_freq_population": c["max_freq_population"],
        "skip_civic": c["skip_civic"],
        "skip_pathogenic": c["skip_pathogenic"],
        "filter_genes": app_genes(case),
    }
    if case.arm == "somatic":
        params.update(
            {
                "vaf_threshold": c["vaf_threshold"],
                "filter_cancervar": list(c["filter_cancervar"]),
                "filter_civic": list(c["filter_civic"]),
                "filter_escat": list(c["filter_escat"]),
            }
        )
    else:
        params.update(
            {
                "vaf_threshold": c["vaf_threshold_germline"],
                "filter_intervar": list(c["filter_intervar"]),
                "filter_renovo": list(c["filter_renovo"]),
            }
        )
    return params


def with_overrides(case: Case, **overrides) -> Case:
    """A copy of ``case`` with extra contract overrides — used by neutrality tests."""
    merged = dict(case.overrides)
    merged.update(overrides)
    return replace(case, overrides=merged)
