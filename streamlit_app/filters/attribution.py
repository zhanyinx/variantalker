"""Why a particular variant is not in the report, and which setting keeps the most out.

Issue #147. The report already says *what the cut was* — every setting live on this arm, in
the parameter page's own control labels, snapshotted at filter time (issue #137's recap) —
and it already says a row was **rejected**, because ``MAFigate_reason`` carries the cell of
the union that decided it. What neither says is *which of those settings* a particular
variant fell outside, which is the question a clinician asks when a variant they expected is
absent. The recap makes that gap sharper rather than smaller: a user can read eleven
settings that all look right and still not know which one dropped the row they came for.

The mechanism, and why it re-implements nothing
-----------------------------------------------
The obvious implementation is to keep each clause's mask instead of their conjunction. It is
not available. ``vendor/pipeline_filters.py``'s ``somatic_filters`` and ``germline_filters``
build their per-clause masks as **local variables** and return only
``filter_guidelines & filter_vaf & filter_genes`` plus the rescue mask; ``common_filters``
likewise returns one conjunction of its two clauses. ``vendor/`` may not be edited — the
drift guard is an ``ast.dump()`` comparison against ``bin/`` — so keeping per-clause masks
would mean recomputing the clauses here, which is the re-implementation
:mod:`filters.variant_filters` exists to have deleted and which issue #16 measured eleven
divergences in, the largest worth 540 rows.

**The cost of that shape is parity, not memory**, which is worth recording because the
ticket priced it the other way round. Measured over 161 real MAFs
(``docs/wayfinder/issue-147/measure147.py``), the per-clause masks are **0.26%** of the
frame typically and **3.6%** worst case — against a frame :func:`~filters.variant_filters.
_label`'s ``assign`` already duplicates in full on every run.

So this module asks the shipped filter instead. Every claim it makes is
:func:`~filters.variant_filters.apply_filters`'s own verdict, re-run with one clause
**neutralised by a parameter value the app can already hold**. Nothing is recomputed, no
vendored signature is touched, and there is no second opinion about the filter here — which
is what makes the result safe to render under issue #136's line that a message states what
MAFigate did.

Two directions, two different questions
---------------------------------------
* **leave-one-in** — everything neutral *except* one clause, and pathogenic retention off.
  A row that fails now is a row that fails **that clause**. This is the set the ticket
  asks for, and :func:`explain_variant` reports it.
* **leave-one-out** — one clause neutral, the rest as the user set them. A rejected row that
  passes now is a row *relaxing that clause alone would bring back*.

They are not the same set, and the difference decides the copy. For AND-composed clauses the
leave-one-out is empty whenever more than one clause fails, and **a rejected row usually
fails more than one**: exactly one clause in 10% of germline and 30% of somatic rejected rows
at Broad, 5% and 8% at Stringent, over 322,238 rows of 161 real MAFs. A surface built on the
leave-one-out alone would therefore go silent on 70–95% of rejected rows — exactly the ones a
user is most likely to be asking about. So the answer is a **list of every closed path**, and
the ticket's "which clause rejected it" is reported as *which clauses*.

Neutral values, and the one clause that has none
------------------------------------------------
The two quality gates use the vendored code's own neutral values: ``>= 0`` over a sum of
non-negative counts, and ``> -1`` over a fraction, which is literally what ``somatic_filters``
writes for its own all-True gene default (``maf["tumor_f"] > -1``).

The guideline OR is the exception and it matters: an **empty** keep-list makes ``isin([])``
all-*False* and narrows the report, so no constant neutralises it. It is neutralised by
**widening** each source to the values the file's own column holds, and a single source is
isolated by **emptying** the others — both vendored string predicates return False on an
empty list, ``has_clinvar_term`` by an explicit ``if not clinvar_keep`` guard and
``has_element_from_list`` by never entering its loop.

Widening cannot rescue a row the file is **blank** on, and that is a measured property rather
than a caveat: 415 somatic rows are admitted by widening ESCAT and not by widening CancerVar,
because CancerVar says nothing at those positions. So a value is always *shown* beside its
setting and never compared to it here — the failure claim comes from the re-run, the value is
only what the file says.

Cost
----
:func:`explain_variant` runs seven single-row filters: **median 9.6ms, p95 10.1ms, max
16.0ms** over 1,920 sampled rejected rows, and *flat in file size*, because the cost is
per-call overhead times clause count rather than row count. Every clause in the filter is
row-wise, so a one-row frame reaches the same verdict — asserted on real data at
**1,920 of 1,920** for both the verdict and the failing set.

:func:`attribute_report` runs six whole-frame filters, 27–68ms per preset on real MAFs. It is
therefore called **once by the filter run** and snapshotted, like the recap beside it, rather
than on every rerender.

Streamlit-free, like the rest of ``filters/``: the parity harness imports this package
without Streamlit installed, and both entry points take a frame and a parameter dict.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import pandas as pd

from config.missing_values import says_nothing
from config.param_labels import label_of
from filters.absent_columns import plan_fills
from filters.variant_filters import (
    FREQUENCY_COLUMNS,
    MAFIGATE_FILTER,
    MAFIGATE_REASON,
    PASS,
    REASON_BOTH,
    REASON_CRITERIA,
    REASON_RESCUE,
    _Settings,
    apply_filters,
)

#: One entry per **thing the user set**, in the order a sentence should name them. The
#: guideline OR is one clause and not four: the block is False only when *every* source is
#: False, so "which source rejected it" has no answer — they all did, together.
CLAUSES: tuple[str, ...] = (
    "exclude_classifications",
    "depth",
    "vaf",
    "guidelines",
    "gene_list",
    "frequency",
)

#: What to call the guideline block on screen. **The one rendered string in this module that
#: is not a control's own label**, and it needs a reason: the block is one clause of the
#: filter drawn by four or five separate controls, so naming it by its controls produces a
#: five-label chain in the middle of a sentence, and naming it by any one of them would be
#: false. The per-source lines underneath carry the control labels, so nothing is hidden —
#: this is a heading over them, not a substitute for them.
CLASSIFICATION_GROUP = "Clinical classification filters"

#: Which parameter carries each clause's value, and what the interface calls it. ``vaf`` is
#: per-arm because the two arms' controls are two widgets with two labels.
CLAUSE_LABEL_IDS: dict[str, object] = {
    "exclude_classifications": "exclude_classifications",
    "depth": "min_depth",
    "vaf": {"somatic": "vaf_threshold", "germline": "vaf_threshold_germline"},
    "gene_list": "gene_selection",
    "frequency": "max_freq_population",
}

#: The columns whose value is worth showing beside a clause. Shown, never compared.
CLAUSE_COLUMNS: dict[str, tuple[str, ...]] = {
    "exclude_classifications": ("Variant_Classification",),
    "depth": ("t_alt_count", "t_ref_count"),
    "vaf": ("tumor_f",),
    "gene_list": ("Hugo_Symbol",),
    "frequency": FREQUENCY_COLUMNS,
}

#: Each clause's neutral parameter value. There is deliberately no ``guidelines`` entry; see
#: the module docstring for why a constant cannot neutralise an OR of keep-lists.
NEUTRAL: dict[str, dict] = {
    "exclude_classifications": {"filter_variant_classification": []},
    "depth": {"min_depth": 0},
    "vaf": {"vaf_threshold": -1.0, "vaf_threshold_germline": -1.0},
    "gene_list": {"filter_genes": []},
    "frequency": {"max_freq_population": 1.0},
}


@dataclass(frozen=True)
class GuidelineSource:
    """One term of the guideline OR: its column, its control, and its name on screen.

    ``name`` is not free vocabulary — it is the word the control's own label already uses
    ("Keep **InterVar** Classifications"), which is why
    ``tests/test_attribution.py::test_source_names_come_from_their_control`` requires each
    one to appear inside :func:`~config.param_labels.label_of` of its ``label_id``. A name
    that stopped matching its control would be a second vocabulary, which is what issue #79
    found this repository is punished for.
    """

    #: The MAF column the vendored clause reads.
    column: str
    #: The :mod:`config.param_labels` row whose control sets this source's keep-list. Not the
    #: same thing as :attr:`param_key`: the label table is keyed by *control*, and ClinVar is
    #: drawn by two controls that partition one parameter (issue #103).
    label_id: str
    #: The key the parameter dict carries this source's keep-list under.
    param_key: str
    #: What to call this source on screen. Must appear inside its control's label.
    name: str


#: The guideline sources per arm, in the order the vendored function ORs them. Germline has
#: no ESCAT term and no CIViC term — that arm's function does not take them — so naming them
#: there would describe a clause that does not exist.
CLASSIFICATION_SOURCES: dict[str, tuple[GuidelineSource, ...]] = {
    "somatic": (
        GuidelineSource("CancerVar", "cancervar_keep", "filter_cancervar", "CancerVar"),
        GuidelineSource(
            "ClinVar_VCF_CLNSIG", "clinvar_pathogenicity", "filter_clinvar", "ClinVar"
        ),
        GuidelineSource("ESCAT", "escat_keep", "filter_escat", "ESCAT"),
        GuidelineSource(
            "CIViC_Evidence_Level", "civic_keep", "filter_civic", "CIViC"
        ),
    ),
    "germline": (
        GuidelineSource("InterVar", "intervar_keep", "filter_intervar", "InterVar"),
        GuidelineSource(
            "ClinVar_VCF_CLNSIG", "clinvar_pathogenicity", "filter_clinvar", "ClinVar"
        ),
        GuidelineSource("RENOVO_Class", "renovo_keep", "filter_renovo", "RENOVO"),
    ),
}

#: What to say where the file has no such column at all — distinct from a blank value,
#: because widening a setting cannot rescue a row the file is silent on either way, but only
#: one of the two is something the user could fix by re-annotating.
NOT_IN_FILE = "not in your file"
SAYS_NOTHING = "nothing"


def shown_value(row: pd.Series, column: str) -> str:
    """What the file says at ``column``, or why it says nothing there.

    The blank test is :func:`config.missing_values.says_nothing` — the app's one answer for
    display — rather than a set written here. This module had its own fifth copy of that set
    briefly, including ``-``, which that module deliberately does *not* treat as a sentinel.
    """
    if column not in row.index:
        return NOT_IN_FILE
    return SAYS_NOTHING if says_nothing(row[column]) else str(row[column])


def label_for(clause: str, arm: str) -> str:
    """The words the interface uses for ``clause`` on ``arm``."""
    if clause == "guidelines":
        return CLASSIFICATION_GROUP
    entry = CLAUSE_LABEL_IDS[clause]
    row_id = entry[arm] if isinstance(entry, dict) else entry
    return label_of(row_id)


@dataclass(frozen=True)
class ClauseOutcome:
    """One setting, and how this variant stands against it."""

    clause: str
    label: str
    #: ``(name, value)`` pairs to show beneath the label. For the guideline block these are
    #: the sources and what the file says at each; for the others, the columns the clause
    #: reads. Only columns the file carries appear, except for the guideline sources, where
    #: an absent column is itself the answer.
    values: tuple[tuple[str, str], ...] = ()
    #: Whether neutralising this clause alone would bring the variant back.
    recoverable: bool = False


@dataclass(frozen=True)
class VariantExplanation:
    """Why one variant is, or is not, in the report."""

    in_report: bool
    #: The cell of the union that put it there — ``criteria``, ``both`` or
    #: ``pathogenic_rescue``. ``None`` when the variant is not in the report.
    reason: str | None = None
    #: Every setting the variant falls outside, in control order. Empty when it is in the
    #: report.
    failing: tuple[ClauseOutcome, ...] = ()
    #: True where the variant is rejected even with every setting neutral and retention off,
    #: which means it holds no usable depth or VAF value and **no setting explains it**.
    #: 13 rows in 322,238 across the measured corpus.
    unreachable: bool = False
    #: Whether pathogenic retention was on and did not apply. False where the user turned it
    #: off, because then it is not a path that closed — it is one they closed.
    retention_declined: bool = False
    #: Filter-input columns this MAF does not carry that a failing clause reads. Non-empty
    #: means part of the answer is about a column the file lacks rather than a setting the
    #: user chose.
    filled_inputs: tuple[str, ...] = field(default_factory=tuple)

    @property
    def sole_reason(self) -> ClauseOutcome | None:
        """The one setting responsible, where there is exactly one."""
        return self.failing[0] if len(self.failing) == 1 else None

    @property
    def recoverable(self) -> tuple[ClauseOutcome, ...]:
        """The settings that would each, alone, bring this variant back."""
        return tuple(outcome for outcome in self.failing if outcome.recoverable)


@dataclass(frozen=True)
class ReportAttribution:
    """How many variants each setting keeps out of one report.

    :attr:`excluded_by` **overlaps** — a row outside two settings is counted by both — so it
    sums to more than :attr:`left_out`. That is the framing issue #147 chose, because it is
    the one that answers *why is my report so small*: on real MAFs the classification block
    and the population-frequency cut each account for 80–98% of the rejected rows, and they
    overwhelmingly do it together. The partitioning alternative — "the one setting
    responsible" — is arithmetically tidier and puts 70–95% of the rows in a single "outside
    more than one setting" bucket, which is honest and nearly useless.
    """

    rows: int
    in_report: int
    #: ``(label, count)``, descending, for the settings that exclude at least one row.
    excluded_by: tuple[tuple[str, int], ...] = ()

    @property
    def left_out(self) -> int:
        return self.rows - self.in_report


def _widened(frame: pd.DataFrame, arm: str) -> dict:
    """Every guideline keep-list widened to the values the file's own columns hold.

    A column the file does not carry is skipped: it is one :mod:`filters.absent_columns`
    fills, and its filled value matches no keep-list, so widening it would claim a term that
    cannot fire.
    """
    override: dict = {}
    for source in CLASSIFICATION_SOURCES[arm]:
        if source.column not in frame.columns:
            continue
        override[source.param_key] = sorted(
            {str(v) for v in frame[source.column].dropna().unique() if str(v).strip()}
        )
    return override


def _all_neutral(frame: pd.DataFrame, arm: str, params) -> dict:
    """Parameters under which no clause excludes anything and the rescue is off.

    Retention off is load-bearing: with it on, ``PASS`` is ``criteria | rescue``, so a row
    the rescue carries would pass whatever the isolated clause says and every clause would
    read as satisfied.
    """
    relaxed = dict(params)
    for override in NEUTRAL.values():
        relaxed.update(override)
    relaxed.update(_widened(frame, arm))
    relaxed["skip_pathogenic"] = True
    relaxed.pop("keep_pathogenic", None)
    return relaxed


def _isolate(frame: pd.DataFrame, arm: str, params, clause: str) -> dict:
    """Everything neutral except ``clause``, restored to what the user set."""
    isolated = _all_neutral(frame, arm, params)
    if clause == "guidelines":
        for source in CLASSIFICATION_SOURCES[arm]:
            key = source.param_key
            if key in params:
                isolated[key] = params[key]
            else:
                isolated[key] = []
        return isolated
    for key in NEUTRAL[clause]:
        if key in params:
            isolated[key] = params[key]
        else:
            isolated.pop(key, None)
    return isolated


def _relax(frame: pd.DataFrame, arm: str, params, clause: str) -> dict:
    """``params`` with ``clause`` alone neutralised."""
    relaxed = dict(params)
    if clause == "guidelines":
        relaxed.update(_widened(frame, arm))
        return relaxed
    relaxed.update(NEUTRAL[clause])
    return relaxed


def _passes(frame: pd.DataFrame, params) -> pd.Series:
    labelled, _ = apply_filters(frame, params)
    return labelled[MAFIGATE_FILTER] == PASS


def _values_for(clause: str, row: pd.Series, arm: str) -> tuple[tuple[str, str], ...]:
    """What to show beneath a clause's label.

    For the guideline block, every source — including ones the file does not carry, because
    an absent source *is* the answer for those rows. For the others, only columns the file
    carries: the frequency layer reads seven names and a real MAF holds two or three, so
    listing the absent ones would bury the answer in "not in your file", and a layer skipping
    a panel it has no data for is not a finding.
    """
    if clause == "guidelines":
        return tuple(
            (source.name, shown_value(row, source.column))
            for source in CLASSIFICATION_SOURCES[arm]
        )
    return tuple(
        (column, shown_value(row, column))
        for column in CLAUSE_COLUMNS.get(clause, ())
        if column in row.index
    )


def explain_variant(one_row: pd.DataFrame, params) -> VariantExplanation:
    """Why the variant in ``one_row`` is, or is not, in the report these ``params`` made.

    Args:
        one_row: a **one-row frame**, sliced from the frame the filter ran on — not a
            ``Series`` transposed back into one. The vendored comparisons raise on object
            dtypes, and a ``Series`` built from a display frame has been through
            ``fillna("")`` and ``astype(str)``, so it would either raise or answer about
            different data. Slicing with ``frame.loc[[index]]`` preserves both dtypes and
            index.
        params: the parameters the run used — ``RunRecap.params``, not the live dict, so the
            explanation describes the report on screen rather than the controls as they now
            stand.

    Returns:
        A :class:`VariantExplanation`. When the variant is in the report the failing set is
        empty and :attr:`~VariantExplanation.reason` names the cell that admitted it — the
        detail dialog is shared by the passed and failed tabs, so the same block has to say
        something true either way.

    Raises:
        ValueError: if ``one_row`` does not hold exactly one row. Every caller slices a
            single row, and a silently-averaged answer over several would be worse than a
            refusal.
    """
    if len(one_row) != 1:
        raise ValueError(
            f"explain_variant takes a one-row frame; given {len(one_row)} rows"
        )

    settings = _Settings.from_params(params)
    arm = settings.arm
    row = one_row.iloc[0]

    labelled, _ = apply_filters(one_row, dict(params))
    if bool((labelled[MAFIGATE_FILTER] == PASS).iloc[0]):
        return VariantExplanation(
            in_report=True, reason=str(labelled[MAFIGATE_REASON].iloc[0])
        )

    failing = []
    for clause in CLAUSES:
        if bool(_passes(one_row, _isolate(one_row, arm, params, clause)).iloc[0]):
            continue
        failing.append(
            ClauseOutcome(
                clause=clause,
                label=label_for(clause, arm),
                values=_values_for(clause, row, arm),
                recoverable=bool(
                    _passes(one_row, _relax(one_row, arm, params, clause)).iloc[0]
                ),
            )
        )

    unreachable = not bool(_passes(one_row, _all_neutral(one_row, arm, params)).iloc[0])

    # Derived here rather than threaded in from the run's ``Diagnostics``, and it cannot
    # disagree with that run: the plan is a function of the frame's columns, the arm and
    # ``skip_pathogenic`` alone, and ``one_row`` is a slice of the very frame the filter ran
    # on, so it carries the same columns. A caller passing its own copy would be a second
    # answer to a question ``absent_columns`` already answers once.
    plan = plan_fills(one_row, arm, skip_pathogenic=settings.skip_pathogenic)
    reads = {
        column
        for outcome in failing
        for column in CLAUSE_COLUMNS.get(outcome.clause, ())
    }
    reads |= {source.column for source in CLASSIFICATION_SOURCES[arm]}
    return VariantExplanation(
        in_report=False,
        failing=tuple(failing),
        unreachable=unreachable and not failing,
        retention_declined=not settings.skip_pathogenic,
        filled_inputs=tuple(sorted(set(plan.filled).intersection(reads))),
    )


def attribute_report(frame: pd.DataFrame, params) -> ReportAttribution:
    """How many of this report's excluded variants each setting is responsible for.

    Six whole-frame filters, 27–68ms per preset on real MAFs, so this is called **once by the
    filter run** and the result snapshotted — the same clock as issue #137's recap, and for
    the same reason: a count re-derived at render time would describe the controls as they
    now stand rather than the cut on screen.

    Counted over the rows this filter **left out**, and overlapping by design — see
    :class:`ReportAttribution`.
    """
    settings = _Settings.from_params(params)
    arm = settings.arm

    labelled, _ = apply_filters(frame, dict(params))
    rejected = labelled[MAFIGATE_FILTER] != PASS
    left_out = int(rejected.sum())

    counts: list[tuple[str, int]] = []
    if left_out:
        for clause in CLAUSES:
            outside = int(
                (rejected & ~_passes(frame, _isolate(frame, arm, params, clause))).sum()
            )
            if outside:
                counts.append((label_for(clause, arm), outside))
    counts.sort(key=lambda pair: (-pair[1], pair[0]))

    return ReportAttribution(
        rows=len(frame),
        in_report=len(frame) - left_out,
        excluded_by=tuple(counts),
    )


#: ``MAFigate_reason`` values that mean the variant is in the report, and what to say about
#: each. Read from :mod:`filters.variant_filters`'s own constants so a renamed cell cannot
#: leave this table naming one that no longer exists.
IN_REPORT_BECAUSE = {
    REASON_CRITERIA: "it met your filter criteria",
    REASON_BOTH: "it met your filter criteria, and it is also called pathogenic",
    REASON_RESCUE: (
        "it is called pathogenic — it did not meet your filter criteria"
    ),
}
