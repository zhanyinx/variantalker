"""Coverage by mutation: which parity divergences the fixture set can still catch.

``test_parity.py`` asserts that the app and the pipeline agree. This measures something
the agreement cannot tell you: whether the fixtures would **notice if they stopped
agreeing**. Agreement is the cheapest property to hold — a fixture whose rows neither
side does anything interesting with agrees trivially — so an empty fixture set and a rich
one both read green there.

The instrument this replaces measured coverage by *attribution* (issue #33's
``test_attribution_coverage.py``): every probe had to explain a row the two sides
disagreed about. That worked while divergences were open and inverted the moment the last
one closed — with nothing diverging, no probe has a witness, so its soundness assertion
relaxed into asserting the **absence** of attribution, and any fixture set reaching row
parity passed it. Issue #242 replaced it with this.

How it measures
---------------
Each of the map's divergences is re-injected into the **app** side, one at a time, and
every case is re-run against a pipeline result computed once and cached. The pipeline is
untouched *by construction*, not by promise: it runs ``bin/filter_variants.py`` in a
subprocess, while a mutation replaces one of the three vendored symbols
``filters.variant_filters`` imports — which is the exact surface issue #33 closed the
divergences on. A fixture set **covers** a divergence when some case notices it.

Mutations are written against the vendored *signatures*, not as textual patches of the old
app: the point is to change the app's decision in the direction the divergence had, not to
reconstruct a historical diff.

What a "catch" is, and what it is not
------------------------------------
A catch is a **verdict-level** difference: both sides reached a PASS set and the sets no
longer match. An app that merely *crashed* under mutation is recorded separately
(:attr:`MutationOutcome.caught_by_crash`) and does not count towards coverage — a fixture
set that only notices a divergence because the app blew up has not shown that the rows
discriminate.

Three states are recorded rather than counted as silence, because each one means "this
case could not have witnessed it" rather than "this case failed to":

* the mutation's target is the other arm's filter (a germline mutation over a somatic case
  never runs);
* the unmutated app already produced no verdict — ``dot_numeric``, where the app refuses
  the file before any filter runs;
* the mutation itself raised :class:`NotApplicable` (no ``DP`` column, no gene list, the
  rescue disabled by ``skip_pathogenic``).

Two of the twelve divergences cannot be scored this way at all, and are named as
:data:`UNSCOREABLE` rather than quietly dropped from the denominator: **#5** (the ``All``
sentinel is a UI concept with no MAF shape) and **#10** (population frequency is the app's
own extra with no pipeline counterpart, and every parity case runs at the neutral
``max_freq_population = 1.0``, so no case exercises it). **#11** is measured directly, by
the harness's own column-prefix rule — see :func:`column_parity`.

Runtime: **~11 s** measured — 33 pipeline subprocesses, run once, plus one in-process app
run per (mutation, applicable case). Deliberately **not** marked ``slow``: ``make
test-fast`` deselects that marker, and this is the check someone wants precisely when they
are about to change the fixtures.

Standalone, to compare two fixture sets on the same mutations::

    python streamlit_app/tests/parity/mutation_coverage.py
    python streamlit_app/tests/parity/mutation_coverage.py --fixture-dir DIR --json out.json
    python streamlit_app/tests/parity/mutation_coverage.py --without-case civic_present

``--fixture-dir`` is what makes the coverage of a *candidate* fixture set measurable
against the coverage of the committed one, which is the question issue #233 asked and
issue #246 acts on. ``--without-case`` withholds a case, which is how the instrument is
shown able to report a gap. The CLI **exits non-zero when any scoreable divergence is
uncovered**, so it can gate a fixture swap from a shell as well as from the suite.
"""

from __future__ import annotations

import argparse
import json
import sys
from contextlib import contextmanager
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

import pandas as pd

HERE = Path(__file__).resolve().parent
STREAMLIT_APP = HERE.parent.parent
if str(STREAMLIT_APP) not in sys.path:
    sys.path.insert(0, str(STREAMLIT_APP))

from tests.parity import contract as C  # noqa: E402
from tests.parity import harness as H  # noqa: E402
from vendor.pipeline_utils import CLINVAR_PATHO, has_clinvar_term  # noqa: E402

#: The three vendored symbols ``filters.variant_filters`` imports and calls by name. A
#: mutation replaces exactly one of them, or ``TARGET_BOTH_ARMS`` for a divergence that
#: lived in a clause both arms share. Prefixed, because a mutation also declares *arms* and
#: the two vocabularies are one line apart in the table below: a target is a symbol to
#: patch, an arm is a case the patch can reach.
TARGET_COMMON = "pipeline_common_filters"
TARGET_SOMATIC = "pipeline_somatic_filters"
TARGET_GERMLINE = "pipeline_germline_filters"
TARGET_BOTH_ARMS = "both_arms"

SOMATIC_ARM = ("somatic",)
GERMLINE_ARM = ("germline",)
BOTH_ARMS = ("somatic", "germline")


class NotApplicable(Exception):
    """The mutation cannot act on this case — recorded, never scored as silence.

    Raised from inside a mutated filter, so it surfaces through
    ``harness.run_app``'s catch-all as a recorded error string rather than propagating.
    :func:`_classify` reads it back off that string; there is no other route out.
    """


# ---------------------------------------------------------------------------
# The mutations
# ---------------------------------------------------------------------------


def _sum_depth(maf: pd.DataFrame, coverage: float) -> pd.Series:
    """The pipeline's depth test, so a mutation can keep it while changing something else."""
    return (maf["t_alt_count"] + maf["t_ref_count"]) >= coverage


def _blanking(real, column: str, is_loose_only, guard=None):
    """Injection by data: blank the cells only the *loose* reading of a match would keep.

    Three of the divergences below are differences in how one string column is matched, and
    the matching lives inside the vendored body. Rather than reimplement that body — which
    would make the mutation a second implementation to get wrong — the frame is edited: the
    cells the loose test matches and the strict test does not are blanked, so the vendored
    function decides exactly as the strict version would have.

    Args:
        real: the vendored function being wrapped.
        column: the column whose cells are blanked.
        is_loose_only: ``(value, kwargs) -> bool``, true for a cell the divergence turns on.
        guard: optional ``(maf, kwargs)`` precondition, raising :class:`NotApplicable` when
            the case cannot host the injection at all.
    """

    def mutated(maf, **kwargs):
        if guard is not None:
            guard(maf, kwargs)
        patched = maf.copy()
        values = patched[column]
        loose_only = values.apply(lambda value: is_loose_only(value, kwargs))
        patched[column] = values.where(~loose_only, "")
        return real(patched, **kwargs)

    return mutated


def mut_vc_include(real):
    """#1 — the app's *include* vocabulary instead of the pipeline's *exclude* list.

    Witnessed by any row whose ``Variant_Classification`` is neither in the exclude list
    nor in the app's eighteen-value vocabulary: the pipeline keeps it, the app drops it.
    ``variant_classification_filter`` is deliberately unread — that is the divergence.
    """

    def mutated(maf, coverage, variant_classification_filter):
        return maf["Variant_Classification"].isin(
            C.APP_VARIANT_CLASSIFICATIONS
        ) & _sum_depth(maf, coverage)

    return mutated


def mut_depth_dp(real):
    """#2 — depth read off ``DP`` instead of ``t_alt_count + t_ref_count``."""

    def mutated(maf, coverage, variant_classification_filter):
        if "DP" not in maf.columns:
            raise NotApplicable("DP absent")
        return (~maf["Variant_Classification"].isin(variant_classification_filter)) & (
            pd.to_numeric(maf["DP"], errors="coerce") >= coverage
        )

    return mutated


def mut_depth_nan_kept(real):
    """#2 — a missing depth kept rather than dropped (NaN comparisons are False)."""

    def mutated(maf, coverage, variant_classification_filter):
        depth = maf["t_alt_count"] + maf["t_ref_count"]
        return (~maf["Variant_Classification"].isin(variant_classification_filter)) & (
            (depth >= coverage) | depth.isna()
        )

    return mutated


def _vaf_variant(real, op):
    """Shift the values, not the comparison — the vendored body owns the ``>``.

    The cheapest faithful injection for a boundary divergence is to make the affected
    values compare the other way at exactly the threshold, which is indistinguishable from
    having written the other operator.
    """

    def mutated(maf, **kwargs):
        vaf = kwargs["vaf"]
        original = maf["tumor_f"]
        patched = maf.copy()
        patched["tumor_f"] = original.where(~op(original, vaf), vaf + 1e-12)
        return real(patched, **kwargs)

    return mutated


def mut_vaf_ge(real):
    """#3 — ``>=`` instead of ``>``: a row at exactly the threshold passes."""
    return _vaf_variant(real, lambda s, v: s == v)


def mut_vaf_nan_kept(real):
    """#3 — a missing VAF kept rather than dropped.

    Filled with 1.0 rather than with the threshold, so the row clears the gene clause's
    ``tumor_f > -1`` too: that clause is also False for NaN, so a threshold-only fill
    would leave the row dropped for the other reason and the mutation would read silent.
    """

    def mutated(maf, **kwargs):
        patched = maf.copy()
        patched["tumor_f"] = patched["tumor_f"].fillna(1.0)
        return real(patched, **kwargs)

    return mutated


def _clnsig_exact(real, keep_key):
    """#4 — ``.isin()`` on the whole ClinVar string instead of splitting it.

    The cells that only match *after* the split are blanked, so the vendored
    ``has_clinvar_term`` behaves as an exact-match test would on this frame.
    """

    def split_only(value, kwargs):
        keep = kwargs[keep_key] if keep_key else CLINVAR_PATHO
        return has_clinvar_term(value, keep) and value not in list(keep)

    return _blanking(real, "ClinVar_VCF_CLNSIG", split_only)


def mut_clnsig_split(real):
    """#4 — in the guideline block, against the case's ``clinvar_keep``."""
    return _clnsig_exact(real, "clinvar_keep")


def mut_clnsig_patho_split(real):
    """#4 — in the pathogenic rescue, against the vendored ``CLINVAR_PATHO``."""
    return _clnsig_exact(real, None)


def mut_germline_escat(real):
    """#6 — the ESCAT clause the app used to OR into the germline guideline block.

    The largest divergence on real data: 51% of attributed diverging rows on the
    50-sample reference sweep, while reading 0% on the fixtures until issue #27 built a
    row whose whole path was pinned to the verdict.
    """

    def mutated(maf, **kwargs):
        keep, rescue = real(maf, **kwargs)
        if "ESCAT" not in maf.columns:
            raise NotApplicable("ESCAT absent")
        escat = maf["ESCAT"].isin(["IA", "IB", "IC", "IIA", "IIB"])
        vaf_ok = maf["tumor_f"] > kwargs["vaf"]
        return (keep | (escat & vaf_ok)), rescue

    return mutated


def mut_escat_dead_clause(real):
    """#7 — bare ``I``/``II``/``III`` ESCAT values retained by the somatic rescue.

    The clause is dead on 182K reference rows (real values are ``IA``, ``IB``, ...), which
    is why the synthetic fixtures construct the bare values: a dead clause can only be
    witnessed by a row the reference does not contain.
    """

    def mutated(maf, **kwargs):
        keep, rescue = real(maf, **kwargs)
        if "ESCAT" not in maf.columns:
            raise NotApplicable("ESCAT absent")
        return keep, (rescue | maf["ESCAT"].isin(["I", "II", "III"]))

    return mutated


def mut_somatic_patho_no_clinvar(real):
    """#7's other half — the somatic rescue without its ClinVar leg.

    Injected by blanking the ``CLNSIG`` values that only the rescue reads: the rescue
    matches ``CLINVAR_PATHO`` while the guideline block matches the case's
    ``clinvar_keep``, so removing exactly the values in the former and not the latter
    leaves the guideline block's decision untouched. The containment that makes that true
    is asserted rather than assumed — a case whose ``clinvar_keep`` reached outside
    ``CLINVAR_PATHO`` would have the mutation quietly editing the guideline block as well.
    """

    def rescue_only(maf, kwargs):
        keep = kwargs["clinvar_keep"]
        if not set(keep) <= set(CLINVAR_PATHO):
            raise NotApplicable(
                f"clinvar_keep {sorted(set(keep) - set(CLINVAR_PATHO))} reaches outside "
                "CLINVAR_PATHO, so this injection is not rescue-only"
            )

    def only_the_rescue_reads_it(value, kwargs):
        return has_clinvar_term(value, CLINVAR_PATHO) and not has_clinvar_term(
            value, kwargs["clinvar_keep"]
        )

    return _blanking(
        real, "ClinVar_VCF_CLNSIG", only_the_rescue_reads_it, guard=rescue_only
    )


def mut_germline_patho_renovo(real):
    """#8 — the germline rescue keyed on ``RENOVO_Class`` instead of ``ClinVar``.

    Both implementations rescue on ``InterVar``, so the substitution only reaches a row
    that InterVar does not already rescue.
    """

    def mutated(maf, **kwargs):
        keep, rescue = real(maf, **kwargs)
        if kwargs.get("skip_pathogenic"):
            raise NotApplicable("rescue disabled by skip_pathogenic")
        renovo = maf["RENOVO_Class"].isin(["LP Pathogenic", "HP Pathogenic"])
        intervar = maf["InterVar"].isin(["Pathogenic", "Likely pathogenic"])
        return keep, ((rescue & ~intervar) | renovo)

    return mutated


def mut_gene_case(real):
    """#9 — the gene list compared case-sensitively.

    Reads the gene file the app's own adapter just wrote, so the symbols compared are the
    ones the app tokenised rather than a second copy of the fixture's list.
    """

    def mutated(maf, **kwargs):
        key = "somatic_genes" if "somatic_genes" in kwargs else "germline_genes"
        path = kwargs[key]
        if path == "null" or not Path(path).exists():
            raise NotApplicable("no gene list")
        genes = pd.read_csv(path, header=None)[0]
        keep, rescue = real(maf, **kwargs)
        case_sensitive = maf["Hugo_Symbol"].isin(genes.values)
        case_insensitive = maf["Hugo_Symbol"].str.upper().isin(genes.str.upper().values)
        lost = case_insensitive & ~case_sensitive
        return (keep & ~lost), rescue

    return mutated


def mut_civic_isin(real):
    """#12 — CIViC matched with ``.isin()`` instead of a substring of the list-repr.

    The map recorded this divergence as open; it was measured not to exist (both sides ran
    a substring test), and issue #33 removed the app's copy entirely. It stays here as the
    mutation that says so: the fixture that separates ``skip_civic`` from column absence is
    the only thing standing between the app and an exact match silently returning.
    """

    def civic_is_read_here(maf, kwargs):
        if "CIViC_Evidence_Level" not in maf.columns:
            raise NotApplicable("CIViC_Evidence_Level absent")
        if kwargs.get("skip_civic"):
            raise NotApplicable("skip_civic")

    def substring_only(value, kwargs):
        keep_levels = list(kwargs["civic_keep"]) + ["A", "B"]
        return (
            pd.notna(value)
            and any(level in str(value) for level in keep_levels)
            and str(value) not in keep_levels
        )

    return _blanking(
        real, "CIViC_Evidence_Level", substring_only, guard=civic_is_read_here
    )


@dataclass(frozen=True)
class Mutation:
    """One divergence, injectable into the app side.

    Attributes:
        name: the divergence number and a short mechanism, matching the vocabulary
            ``harness.PROBES`` used — so a reader comparing the two instruments sees the
            same names.
        divergence: the map's number, for the accounting in :data:`THE_TWELVE`.
        target: which vendored symbol the injection replaces, or :data:`TARGET_BOTH_ARMS`.
        arms: the arms whose cases can witness it. A germline mutation over a somatic case
            is not silence — the patched function is never called.
        make: ``(real) -> mutated``, closing over the symbol it replaces. A mutation that
            replaces a whole function rather than wrapping it ignores ``real``, which is
            what an *include*-list or a different depth column amounts to.
    """

    name: str
    divergence: str
    target: str
    arms: tuple[str, ...]
    make: Callable[[Callable], Callable]

    @property
    def why(self) -> str:
        """The mutation's first docstring line, so the report explains itself."""
        doc = (self.make.__doc__ or "").strip().splitlines()
        return doc[0] if doc else ""


MUTATIONS: list[Mutation] = [
    Mutation("#1 vc_exclude_vs_include", "#1", TARGET_COMMON, BOTH_ARMS, mut_vc_include),
    Mutation("#2 depth_dp_vs_sum", "#2", TARGET_COMMON, BOTH_ARMS, mut_depth_dp),
    Mutation("#2 depth_nan_kept_by_app", "#2", TARGET_COMMON, BOTH_ARMS, mut_depth_nan_kept),
    Mutation("#3 vaf_at_threshold", "#3", TARGET_BOTH_ARMS, BOTH_ARMS, mut_vaf_ge),
    Mutation("#3 vaf_nan_kept_by_app", "#3", TARGET_BOTH_ARMS, BOTH_ARMS, mut_vaf_nan_kept),
    Mutation("#4 clnsig_split_vs_exact", "#4", TARGET_BOTH_ARMS, BOTH_ARMS, mut_clnsig_split),
    Mutation(
        "#4 clnsig_patho_split_vs_exact",
        "#4",
        TARGET_BOTH_ARMS,
        BOTH_ARMS,
        mut_clnsig_patho_split,
    ),
    Mutation(
        "#6 germline_escat_unmirrored",
        "#6",
        TARGET_GERMLINE,
        GERMLINE_ARM,
        mut_germline_escat,
    ),
    Mutation(
        "#7 escat_dead_retain_clause",
        "#7",
        TARGET_SOMATIC,
        SOMATIC_ARM,
        mut_escat_dead_clause,
    ),
    Mutation(
        "#7 somatic_patho_no_clinvar_clause",
        "#7",
        TARGET_SOMATIC,
        SOMATIC_ARM,
        mut_somatic_patho_no_clinvar,
    ),
    Mutation(
        "#8 germline_patho_renovo_vs_clinvar",
        "#8",
        TARGET_GERMLINE,
        GERMLINE_ARM,
        mut_germline_patho_renovo,
    ),
    Mutation("#9 gene_case_sensitivity", "#9", TARGET_BOTH_ARMS, BOTH_ARMS, mut_gene_case),
    Mutation(
        "#12 civic_isin_vs_substring", "#12", TARGET_SOMATIC, SOMATIC_ARM, mut_civic_isin
    ),
]

#: Divergences with no MAF-level injection at all. Named, and **excluded from the
#: denominator explicitly** rather than by omission: a divergence that silently left the
#: table would read as full coverage of a smaller map.
UNSCOREABLE = {
    "#5": (
        "the All sentinel is a UI concept with no MAF shape — the fixture README says so, "
        "and the sentinel never reaches the vendored filters, so no MAF can witness it"
    ),
    "#10": (
        "population frequency is the app's own extra, which the pipeline has no "
        "counterpart for. Every parity case runs at the neutral max_freq_population=1.0, "
        "so no case exercises it and no fixture set can be scored on it here. Its "
        "neutrality at 1.0 is asserted directly, in test_parity.py"
    ),
}

#: Divergences measured by the harness's own rule rather than by mutation.
MEASURED_DIRECTLY = {
    "#11": "output columns against the pipeline's keep — see column_parity()",
}

#: The map's twelve divergences, titled as the fixture README's coverage table titles
#: them. Every one must be scored, unscoreable-by-name, or measured directly — asserted by
#: ``test_mutation_coverage.py``, which is what stops the table above from quietly
#: shrinking.
THE_TWELVE = {
    "#1": "Variant_Classification exclude vs include",
    "#2": "depth: t_alt+t_ref vs DP",
    "#3": "VAF > vs >=, NaN",
    "#4": "ClinVar split vs .isin()",
    "#5": "All sentinel / always-OR",
    "#6": "germline ESCAT filter",
    "#7": "somatic patho-retain",
    "#8": "germline patho-retain",
    "#9": "gene list case + format",
    "#10": "population frequency",
    "#11": "output columns vs KEEP",
    "#12": "CIViC list-repr substring vs .isin()",
}


# ---------------------------------------------------------------------------
# Installing a mutation
# ---------------------------------------------------------------------------


@contextmanager
def mutated(mutation: Mutation):
    """Patch ``filters.variant_filters`` for the body, and restore it afterwards.

    Restoration is ``finally``-guaranteed and additionally asserted by
    ``test_mutation_coverage.py``: a leaked mutation would not fail here, it would fail
    somewhere else in the session, as a parity regression in a module that never asked for
    one.
    """
    import filters.variant_filters as V

    targets = (
        (TARGET_SOMATIC, TARGET_GERMLINE)
        if mutation.target == TARGET_BOTH_ARMS
        else (mutation.target,)
    )
    saved = [(attr, getattr(V, attr)) for attr in targets]
    try:
        for attr, real in saved:
            setattr(V, attr, mutation.make(real))
        yield
    finally:
        for attr, real in saved:
            setattr(V, attr, real)


# ---------------------------------------------------------------------------
# Measuring
# ---------------------------------------------------------------------------

#: One pipeline subprocess per (fixture directory, case), for the life of the process.
#: Keyed on the directory as well as the case because ``--fixture-dir`` repoints the whole
#: case set at files of the same names.
_pipeline_cache: dict[tuple[str, str], H.PipelineResult] = {}


def pipeline_result(case: C.Case) -> H.PipelineResult:
    """The pipeline's answer for one case, computed at most once per fixture directory."""
    key = (str(C.FIXTURE_DIR), case.name)
    if key not in _pipeline_cache:
        _pipeline_cache[key] = H.run_pipeline(case)
    return _pipeline_cache[key]


@contextmanager
def fixtures_from(directory: Path | None):
    """Point the whole case set at another directory of fixtures, and put it back.

    ``Case.maf_path`` resolves against a module global, so repointing the case set means
    rebinding it in two modules and clearing the loader's cache. Restoring it is not
    tidiness: this module is importable inside pytest, and a repointing left in place would
    silently hand every later parity assertion in the session a different set of MAFs.
    """
    if directory is None:
        yield
        return

    previous = C.FIXTURE_DIR
    C.FIXTURE_DIR = H.FIXTURE_DIR = directory
    H._read_maf_cached.cache_clear()
    try:
        yield
    finally:
        C.FIXTURE_DIR = H.FIXTURE_DIR = previous
        H._read_maf_cached.cache_clear()


def agrees(pipeline: H.PipelineResult, app: H.AppResult) -> bool:
    """Whether the two sides said the same thing, by the harness's own rule.

    ``ParityResult.errors_in_parity`` rather than a local "both failed" test, and that
    matters in both directions: it accepts the app's typed refusal against the pipeline's
    ``TypeError`` on ``dot_numeric`` (the same non-verdict), and it refuses to accept any
    *other* app-side exception against a pipeline raise — which is what stops a mutation
    that merely crashes the app from reading as agreement.
    """
    verdict = H.ParityResult(
        case="",
        arm="",
        fixture="",
        rows=0,
        pipeline_error=pipeline.error,
        app_error=app.error,
        app_refused_columns=app.refused_columns,
        app_filled_columns=app.filled_columns,
    )
    if not verdict.errors_in_parity:
        return False
    if pipeline.error is not None or app.error is not None:
        return True
    return pipeline.passed == app.passed


def _not_applicable_reason(app: H.AppResult) -> str | None:
    """The mutation's own ``NotApplicable``, read back off the recorded error string.

    ``harness.run_app`` converts every exception to ``"TypeName: message"``, so this is
    the only route out of a mutated filter — there is no exception left to catch.
    """
    if app.error and app.error.startswith(f"{NotApplicable.__name__}:"):
        return app.error.split(":", 1)[1].strip()
    return None


@dataclass
class MutationOutcome:
    """What one mutation did to every case."""

    mutation: str
    divergence: str
    why: str
    caught_by: list[str] = field(default_factory=list)
    #: Cases where the mutated app raised something other than :class:`NotApplicable`.
    #: Recorded, and deliberately **not** counted as coverage — see the module docstring.
    caught_by_crash: list[str] = field(default_factory=list)
    silent_on: list[str] = field(default_factory=list)
    not_applicable: dict[str, str] = field(default_factory=dict)

    @property
    def covered(self) -> bool:
        return bool(self.caught_by)

    def to_json(self) -> dict:
        return {
            "divergence": self.divergence,
            "why": self.why,
            "covered": self.covered,
            "caught_by": sorted(self.caught_by),
            "caught_by_crash": sorted(self.caught_by_crash),
            "silent_on": sorted(self.silent_on),
            "not_applicable": dict(sorted(self.not_applicable.items())),
        }


@dataclass
class CoverageReport:
    """The measurement. Every count in it is derived, so nothing can drift apart."""

    fixture_dir: str
    cases: list[str]
    #: Cases the two sides already disagree about *before* any mutation. Any entry here
    #: invalidates the whole run: a case that is off parity to begin with cannot report
    #: whether it noticed a mutation, so it is excluded and named.
    disagreeing_unmutated: dict[str, str] = field(default_factory=dict)
    outcomes: dict[str, MutationOutcome] = field(default_factory=dict)
    #: Per-case PASS / NOPASS / output columns, straight off the pipeline — the fixture
    #: README's table, measured.
    table: dict[str, dict] = field(default_factory=dict)
    column_parity: dict[str, dict] = field(default_factory=dict)

    @property
    def gaps(self) -> list[str]:
        """The scoreable divergences no case noticed. The suite's verdict."""
        return sorted(n for n, o in self.outcomes.items() if not o.covered)

    @property
    def score(self) -> str:
        return f"{len(self.outcomes) - len(self.gaps)}/{len(self.outcomes)}"

    @property
    def widths(self) -> list[int]:
        """The distinct pipeline output widths the case set exercises."""
        return sorted({v["pipeline_columns"] for v in self.column_parity.values()})

    def to_json(self) -> dict:
        return {
            "fixture_dir": self.fixture_dir,
            "cases": sorted(self.cases),
            "disagreeing_unmutated": dict(sorted(self.disagreeing_unmutated.items())),
            "score": self.score,
            "gaps": self.gaps,
            "mutations": {n: o.to_json() for n, o in sorted(self.outcomes.items())},
            "unscoreable": UNSCOREABLE,
            "measured_directly": MEASURED_DIRECTLY,
            "table": dict(sorted(self.table.items())),
            "column_parity": dict(sorted(self.column_parity.items())),
            "output_widths": self.widths,
        }


def _case_table_row(case: C.Case, pipeline: H.PipelineResult, app: H.AppResult) -> dict:
    if pipeline.error is not None:
        return {
            "pipeline_error": pipeline.error,
            "app_error": app.error,
            "app_refused_columns": app.refused_columns,
            "agree": agrees(pipeline, app),
        }
    passed = len(pipeline.passed)
    return {
        "pass": passed,
        "nopass": len(pipeline.verdicts) - passed,
        "output_columns": len(pipeline.output_columns),
        "agree": agrees(pipeline, app),
    }


def column_parity(cases, pipeline_results: dict[str, H.PipelineResult]) -> dict:
    """#11 — the app's column list against the pipeline's output, by the harness's rule.

    Not a mutation, and not a length comparison: ``ParityResult.columns_in_parity`` checks
    that the app's list *opens with* the pipeline's expected prefix position by position,
    which is the check that caught a 40-against-40 substitution before issue #35. The
    **widths** exercised are the coverage — a set that only ever reaches 40 columns has
    lost the ``KEEP``-completeness the CIViC fixture bought.

    Runs off the cached pipeline results rather than calling ``harness.compare``, which
    would spend a second subprocess per case to re-measure what is already in hand.
    """
    out = {}
    for case in cases:
        pipeline = pipeline_results[case.name]
        if pipeline.error is not None or not pipeline.output_columns:
            continue
        frame = H.load_frame(case)
        result = H.ParityResult(
            case=case.name,
            arm=case.arm,
            fixture=case.fixture,
            rows=len(frame),
            pipeline_columns=pipeline.output_columns,
            app_columns=H.app_columns(case, frame),
        )
        out[case.name] = {
            "pipeline_columns": len(pipeline.output_columns),
            "expected_prefix": len(result.expected_prefix),
            "app_columns": len(result.app_columns),
            "in_parity": bool(result.columns_in_parity),
        }
    return out


def measure(cases=None, mutations=None, fixture_dir: Path | None = None) -> CoverageReport:
    """Run every mutation over every applicable case.

    Args:
        cases: the cases to measure over. Defaults to the whole contract case set;
            withholding one is how the instrument is shown able to report a gap.
        mutations: the mutations to inject. Defaults to all of :data:`MUTATIONS`.
        fixture_dir: repoint the case set at a different directory of fixtures with the
            same file names, for the duration of this measurement only — see
            :func:`fixtures_from`, which puts it back.
    """
    with fixtures_from(fixture_dir):
        return _measure(
            list(C.CASES if cases is None else cases),
            list(MUTATIONS if mutations is None else mutations),
        )


def _measure(cases: list, mutations: list) -> CoverageReport:
    report = CoverageReport(
        fixture_dir=str(C.FIXTURE_DIR), cases=[c.name for c in cases]
    )

    unmutated: dict[str, H.AppResult] = {}
    for case in cases:
        pipeline = pipeline_result(case)
        app = H.run_app(case)
        unmutated[case.name] = app
        report.table[case.name] = _case_table_row(case, pipeline, app)
        if not agrees(pipeline, app):
            report.disagreeing_unmutated[case.name] = (
                f"pipeline_error={pipeline.error!r} app_error={app.error!r} "
                f"pipeline_pass={len(pipeline.passed)} app_pass={len(app.passed)}"
            )

    for mutation in mutations:
        outcome = MutationOutcome(
            mutation=mutation.name, divergence=mutation.divergence, why=mutation.why
        )
        with mutated(mutation):
            for case in cases:
                if case.name in report.disagreeing_unmutated:
                    outcome.not_applicable[case.name] = "off parity before mutation"
                    continue
                if case.arm not in mutation.arms:
                    outcome.not_applicable[case.name] = (
                        f"{case.arm} case, mutation targets {'/'.join(mutation.arms)}"
                    )
                    continue
                if unmutated[case.name].error is not None:
                    outcome.not_applicable[case.name] = (
                        f"the app reaches no verdict on this file unmutated "
                        f"({unmutated[case.name].error})"
                    )
                    continue

                app = H.run_app(case)
                reason = _not_applicable_reason(app)
                if reason is not None:
                    outcome.not_applicable[case.name] = reason
                elif app.error is not None:
                    outcome.caught_by_crash.append(case.name)
                elif agrees(pipeline_result(case), app):
                    outcome.silent_on.append(case.name)
                else:
                    outcome.caught_by.append(case.name)
        report.outcomes[mutation.name] = outcome

    report.column_parity = column_parity(cases, {c.name: pipeline_result(c) for c in cases})
    return report


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------


def render(report: CoverageReport) -> str:
    """The report as a human reads it: one line per divergence, the two excuses named."""
    lines = [
        f"fixtures: {report.fixture_dir}",
        f"cases:    {len(report.cases)}",
    ]
    bad = report.disagreeing_unmutated
    lines.append(
        "unmutated baseline: all cases agree"
        if not bad
        else f"unmutated baseline: {len(bad)} case(s) DISAGREE -> {sorted(bad)}"
    )
    lines.append("")
    lines.append(f"{'divergence':38s} {'covered':8s} caught_by")
    lines.append("-" * 100)
    for name, outcome in report.outcomes.items():
        caught = sorted(outcome.caught_by)
        shown = ", ".join(caught[:4]) + (f" (+{len(caught) - 4})" if len(caught) > 4 else "")
        crashed = (
            f"   [crash-only: {sorted(outcome.caught_by_crash)}]"
            if outcome.caught_by_crash and not caught
            else ""
        )
        lines.append(
            f"{name:38s} {'yes' if outcome.covered else 'NO':8s} {shown}{crashed}"
        )
    lines.append("")
    for number, why in {**UNSCOREABLE, **MEASURED_DIRECTLY}.items():
        lines.append(f"{number} {THE_TWELVE[number]:34s} {'n/a':8s} {why}")
    lines.append("")
    mismatched = sorted(k for k, v in report.column_parity.items() if not v["in_parity"])
    lines.append(
        f"column parity (#11): {len(report.column_parity)} case(s) compared, "
        f"{len(mismatched)} mismatched" + (f" -> {mismatched}" if mismatched else "")
    )
    lines.append(f"  distinct pipeline output widths exercised: {report.widths}")
    lines.append("")
    lines.append(f"score: {report.score} scoreable divergences caught")
    if report.gaps:
        lines.append(f"GAPS:  {report.gaps}")
    return "\n".join(lines)


def main(argv=None) -> int:
    """Measure one fixture set from a shell, and exit non-zero if anything is uncovered.

    Two flags, each with a caller: ``--fixture-dir`` is how a candidate fixture set is
    measured against the same mutations (issues #233 and #246), and ``--without-case``
    is how a reader reproduces the falsifiability demonstration by hand.
    """
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--fixture-dir", type=Path, default=None)
    parser.add_argument("--json", type=Path, default=None)
    parser.add_argument(
        "--without-case",
        action="append",
        default=[],
        help="withhold a case, to show the instrument reporting a gap",
    )
    args = parser.parse_args(argv)

    if not H.PIPELINE_AVAILABLE:
        print("bin/ is absent: the pipeline side cannot run", file=sys.stderr)
        return 2

    withheld = set(args.without_case)
    unknown = withheld - set(C.CASES_BY_NAME)
    if unknown:
        print(f"unknown case(s): {sorted(unknown)}", file=sys.stderr)
        return 2
    cases = [c for c in C.CASES if c.name not in withheld]

    report = measure(cases=cases, fixture_dir=args.fixture_dir)
    print(render(report))
    if args.json:
        args.json.write_text(json.dumps(report.to_json(), indent=2, sort_keys=True) + "\n")
    return 1 if report.gaps else 0


if __name__ == "__main__":
    raise SystemExit(main())
