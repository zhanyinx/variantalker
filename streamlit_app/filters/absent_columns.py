"""Which filter-input columns a MAF must carry, and what the app does when it doesn't.

The sibling of ``filters/numeric_columns.py``, and the two split the same question along a
mechanical line issue #28 draws:

* a filter-input column that is **present but unreadable** — a refusal, next door;
* a filter-input column that is **absent** — a neutral fill and a warning, here.

The line is mechanical rather than a matter of taste. An unreadable column means *this file
says something the pipeline cannot interpret*; an absent one means *this file never said it*.
The first has an offending value to name and no verdict to reach; the second has no value at
all, and a report is still reachable if the app decides what the missing column would have
said.

Off parity, and recorded as such
--------------------------------
The pipeline raises ``KeyError`` on these files. It does not fill anything, so **there is no
pipeline verdict for a filled run to match** — the app is inventing a behaviour rather than
reproducing one. That is a deliberate deviation, taken because the app is interactive and a
clinician who has already waited for a load is better served by a report they are told to
distrust than by a stack trace. It is recorded as intentional in the parity suite
(``tests/parity/test_absent_columns.py``) rather than quietly excluded from it.

The risk direction inverts intuition
------------------------------------
This is the part that decides how the warnings are built, and it is the opposite of what the
shapes suggest. Measured against the 411 / 818 reference baselines (issue #20, `measure20c.py`):

===========================  =========  ==========  ==============================
column filled                neutral    vs PASS     ``filter_patho``
===========================  =========  ==========  ==============================
somatic ``CancerVar``        ``.``      −287 (−70%) 388 → **100**
somatic ``ClinVar_VCF_CLNSIG`` NaN      −54         388 → 334
somatic ``ESCAT``            ``.``      −23         388
somatic ``Variant_Classification`` unlisted  **+86**  388
germline ``RENOVO_Class``    ``.``      −322 (−39%) 496
germline ``InterVar``        ``.``      −147        496 → **90**
germline ``ClinVar_VCF_CLNSIG`` NaN     −5          496 → 475
germline ``t_alt_count``     numeric    **+76**     496
germline ``tumor_f``         numeric    **+27**     496
===========================  =========  ==========  ==============================

Filling a **guideline source** — which looks harmless, since dropping a term from an OR is
algebraically clean — *shrinks* a clinical report, and the rows it removes are invisible to
whoever reads it. Filling a **quality gate** — which looks alarming, since it disables a
threshold outright — only ever *adds* rows, and an extra row is on screen to be judged. So
the family that looks safe is the dangerous one. Issue #14 settled this effort's principle as
"extra rows are visible, missing rows are not"; :data:`NEUTRAL_FILLS` records which direction
each fill pushes so the warning can say which one the user is looking at.

The escalation
--------------
Four columns are read by the pipeline's **pathogenic retention** as well as by the guideline
OR. Filling one of those neutral drops a term from an OR that runs *beside* the criteria
rather than inside them, which is why the same algebraic move costs so much more there —
somatic 388 → 100 for ``CancerVar`` and germline 496 → 90 for ``InterVar``. Those get a
distinct, escalated warning rather than the same quiet dropout line, and filtering still
proceeds.

The escalation is by **which rule reads the column**, not by how much the fill costs, and
issue #136 measured the gap that opens between the two: ``ClinVar_VCF_CLNSIG`` escalates on
both arms and costs 388 → 334 and 496 → 475, a −14% and a −4% against the other two's −74%
and −82%. The warning therefore states the mechanism and not the magnitude, since one
sentence has to be true at both ends — and no per-column severity is recorded, which would be
a fourth measured field to keep in step with the reference for a distinction the escalation
does not itself draw.

Issue #136 also counted how this warning is actually reached. On the file's **own** arm it is
nearly unreachable: 2 of 173 placeable real MAFs, both germline files missing ``InterVar``.
Every other time it fires, it fires because the *arm* is wrong — 63 of 63 germline files draw
it, and nothing else, when filtered as somatic. That is why the sentence blaming the file had
to go, and why nothing here branches on the arm to do it: the mismatch notice drawn directly
above (issue #135) supplies the cause, and the false clause was false on the arm-correct case
too. :data:`PATHOGENIC_INPUTS` is derived, not written down, and
``tests/test_absent_columns.py`` asserts the trigger set exactly — including the two
counter-intuitive members of the boundary:

* ``RENOVO_Class`` is **not** escalated even though it is the single largest source loss at
  −322 rows (−39%), because germline pathogenic retention has no RENOVO clause. Its fill costs
  guideline rows, which is bad; it does not disarm the safety net, which is worse.
* ``CIViC_Evidence_Level`` **is** read by somatic pathogenic retention, and is still never
  escalated — because it is never filled. It is the one filter input the pipeline itself
  guards (``check_civic_column_exists``), so an absent CIViC column is real pipeline behaviour
  and the app is on parity by doing nothing about it. It is also the only filter-input column
  actually absent in the reference: 100 of 100 files.

Both are asserted by name, because both read as bugs to anyone who has not seen the
measurements and will be "fixed" otherwise.

The lists are derived, not written down
---------------------------------------
:data:`REQUIRED_INPUTS` and :data:`PATHOGENIC_INPUTS` are computed at import by parsing the
same vendored source the app actually calls, exactly as ``numeric_columns`` derives its three.
A hand-written list would be a fifth copy of a pipeline decision in a codebase whose entire
parity problem began with lists someone copied correctly and then stopped updating.

What "required" means here is precise and is what makes the CIViC carve-out fall out rather
than being special-cased: a column is required when the vendored body reads it **on every path
through the function**. ``CIViC_Evidence_Level`` appears only inside one branch of the
pipeline's own ``if``; ``Hugo_Symbol`` only inside the gene-file branch; ``CancerVar`` appears
inside branches too, but in *both* arms of them, so it survives the intersection. That is the
same distinction as "would an absent column raise", reached from the source rather than by
dropping columns one at a time and watching.

The **values** are not derived and could not be. Which value is neutral for a column is a
clinical judgement about what the missing annotation would have said, measured per column
against the reference; only the *set of columns needing one* comes from the source. The two
are held together by an import-time check that every derived column has a fill.
"""

from __future__ import annotations

import ast
from dataclasses import dataclass
from typing import Iterable, Sequence

import pandas as pd

from filters.notes import Note
from filters.numeric_columns import vendored_source
from filters.vendored_ast import frame_parameters, functions as functions_in, is_frame_column

#: Which vendored functions decide each arm — the same three names
#: ``filters/variant_filters.py`` imports and calls. Written here because a source parser has
#: no way to know which arm a function serves; a rename in the vendored code breaks that
#: module's import before it could silently empty a set here.
ARM_FILTERS: dict[str, tuple[str, ...]] = {
    "somatic": ("common_filters", "somatic_filters"),
    "germline": ("common_filters", "germline_filters"),
}

#: The pipeline's own name for its pathogenic-retention mask, in every arm's body. The one
#: identifier this module has to know, and the escalation is the whole reason: the columns that
#: reach *this* variable are the ones whose absence disarms the safety net, and nothing in the
#: shape of the code distinguishes them from the guideline block otherwise.
PATHOGENIC_MASK = "filter_patho"

#: Rows removed, rows added — which way a fill moves the report. The distinction the module
#: docstring is about, carried per column so the warning can state it rather than the user
#: having to infer it from the column's name.
REMOVES = "removes"
ADDS = "adds"

#: The phrase that marks an escalated warning, and the whole of :func:`is_escalated`'s
#: implementation. It is a phrase in the text rather than a marker smuggled alongside it
#: because it is what the sentence *says* — a caller re-deriving which warning is which by
#: matching against ``degraded_columns``, say, can get it wrong on a MAF where a plain fill
#: happens to name a pathogenic column in passing, and this cannot.
#:
#: Since issue #151 it is no longer what the app branches on: the level travels with the
#: message as :attr:`filters.notes.Note.level`. See :func:`is_escalated` for why the text test
#: survives anyway, and :mod:`filters.notes` for why it was not extended to a second level.
ESCALATION_MARKER = "PATHOGENIC RETENTION DEGRADED"

def _not_a_complete_result() -> str:
    """The statement every fill warning has to carry, spelled in one place.

    Both notes below make it, so two copies is one copy too many. It replaced "off parity by
    construction" when that word left the interface (issue #51): the user-facing point was
    never the comparison, it was that a report built on a stand-in value must not be read as
    a whole one.

    Issue #136 found the replacement had kept the comparison after all. It read *"a file
    without this column would normally be refused rather than filtered"*, and **MAFigate never
    refuses on absence** — absence fills and warns, here; only an unreadable *value* refuses,
    next door. So "normally" could only mean *in the pipeline*, which is the sales pitch
    decision 2 of the map retired, and the sentence made a claim about the user's file where
    the app's own behaviour was the subject.

    That reading is what makes this a wrong-arm fix without a wrong-arm branch. A germline MAF
    on the somatic arm is a complete file being told it would be refused, which is false and
    sends the user to their annotation pipeline instead of the Sample Type control; a genuinely
    somatic MAF missing ``CancerVar`` is told the same thing, which is *also* not what MAFigate
    would do with it. One sentence true of every run fixes both, and ``filters/`` keeps its
    freedom from Streamlit and from :mod:`filters.arm_detection` — the mismatch notice drawn
    directly above this (issue #135) is what supplies the *why*.

    One column makes the old sentence true and is unreachable through the app:
    ``Variant_Classification`` is also in ``REQUIRED_COLUMNS["core"]``, so a MAF without it is
    stopped at load and never filtered. The other eight, including both that escalate, are the
    ones a user can actually meet.
    """
    return (
        "MAFigate filtered your file anyway, so **this report is not a complete result**."
    )


def is_escalated(text: str) -> bool:
    """Whether this warning is one of the escalated ones, read off the text alone.

    **No longer what the app renders from.** Since issue #151 the level travels with the
    message as :attr:`filters.notes.Note.level`, set by whichever builder writes the sentence,
    because a second level would have meant testing the leading emoji and an emoji is copy an
    editorial pass may change. This survives for the callers that only have the string: the
    tests that assert the escalated sentence is the one that says so, and any reader checking a
    recorded message after the fact.

    It stays honest for the same reason it always worked — the phrase it looks for *is* the
    copy, not a marker smuggled alongside it — so the two answers cannot disagree while
    :func:`_degraded_note` is the only builder that writes :data:`ESCALATION_MARKER`, which
    ``tests/test_filter_notes.py`` asserts directly.
    """
    return ESCALATION_MARKER in text


@dataclass(frozen=True)
class Fill:
    """One column's neutral value, the direction it moves the report, and why it is neutral.

    Attributes:
        value: what the mask-computation copy carries for this column. Never written into a
            frame anyone can see — see :meth:`FillPlan.frame_for_masks`.
        direction: :data:`REMOVES` or :data:`ADDS`, measured, not reasoned about.

    What makes each particular value neutral is recorded as a comment beside it in
    :data:`NEUTRAL_FILLS` rather than as a third field. That text is read by people and never
    by code, and a field nothing reads is one more thing to keep in step for no gain.
    """

    value: object
    direction: str


#: The neutral value per filter-input column. Every entry is a measured decision (issue #20);
#: the module docstring carries the table of what each one costs.
#:
#: Two of them are worth reading twice.
#:
#: ``tumor_f`` is filled with ``1.0`` rather than something larger because a VAF *has* a domain
#: maximum: the pipeline's gate is ``tumor_f > vaf`` with ``vaf`` a fraction, so 1.0 is the
#: largest value a real measurement could have taken and therefore neutral wherever neutrality
#: is reachable at all. At a threshold of exactly 1.0 nothing passes — but nothing passes there
#: with a *real* column either, so the fill is not what emptied the report.
#:
#: The depth columns get ``inf``, which no measurement could be, because read depth has no
#: domain maximum to reach for. Issue #20 measured this fill at 1e9 and any finite choice is a
#: bet that no user sets a larger threshold; ``inf`` is that bet's limit and settles it. It
#: reaches only ``(t_alt_count + t_ref_count) >= coverage``, where it is neutral for every
#: threshold rather than every plausible one.
NEUTRAL_FILLS: dict[str, Fill] = {
    # -- Guideline sources: the fill makes the source contribute nothing to the OR --------
    #
    # isin(keep) is False for a value no keep-list carries, so the clause contributes nothing
    # to the guideline union.
    "CancerVar": Fill(".", REMOVES),
    "InterVar": Fill(".", REMOVES),
    # The same rule, and the largest single-source loss on the reference — see the module
    # docstring on why it is still not escalated.
    "RENOVO_Class": Fill(".", REMOVES),
    # has_clinvar_term returns False on a missing value, which is *exactly* neutral rather than
    # merely unmatched: the vendored function tests pd.notna before looking at any term.
    "ClinVar_VCF_CLNSIG": Fill(float("nan"), REMOVES),
    # The pipeline's own blank for this column, in 88,609 of 89,019 reference rows — the only
    # one of these fills that reproduces an encoding the pipeline itself emits.
    "ESCAT": Fill(".", REMOVES),
    # -- Quality gates: the fill makes the gate a no-op -----------------------------------
    #
    # The gate is ~isin(excluded), so any value the user cannot have excluded passes every row.
    # A sentinel rather than a real classification because the exclude list is the *user's*,
    # and a real name could be in it — which would empty the report instead of widening it.
    "Variant_Classification": Fill("__MAFIGATE_ABSENT__", ADDS),
    # The depth gate is (t_alt_count + t_ref_count) >= coverage, which inf satisfies at every
    # threshold; the second is the other half of the same sum.
    "t_alt_count": Fill(float("inf"), ADDS),
    "t_ref_count": Fill(float("inf"), ADDS),
    # The gate is tumor_f > vaf with vaf a fraction, so 1.0 is the largest value a real VAF
    # could have taken.
    "tumor_f": Fill(1.0, ADDS),
}


# ---------------------------------------------------------------------------
# Reading the contract out of the vendored source
# ---------------------------------------------------------------------------


def derive_required_inputs(source: str) -> dict[str, tuple[str, ...]]:
    """Per arm, the MAF columns ``source``'s filters read on **every** path, sorted.

    Args:
        source: Python source text — normally :func:`~filters.numeric_columns.vendored_source`,
            and mutated copies of it in the tests, which is what makes this checkable as a
            derivation rather than as a remembered answer.

    Returns:
        ``{"somatic": (...), "germline": (...)}``. A column here is one whose absence makes the
        vendored call raise ``KeyError`` whatever the parameters are, which is precisely the
        set this module has to fill.

    Every-path rather than every-mention, because the difference *is* the CIViC carve-out and
    the gene clause. A branch-local read is a read the pipeline has already decided it can do
    without.
    """
    return _per_arm(
        source, lambda function, frames: _always_read(function.body, frames)
    )


def derive_pathogenic_inputs(source: str) -> dict[str, tuple[str, ...]]:
    """Per arm, the MAF columns that reach the pipeline's pathogenic-retention mask, sorted.

    Every path, not every-path: a column that decides the rescue on *any* branch is a column
    whose absence can degrade the rescue, and the branches here are ``skip_civic`` and
    ``skip_pathogenic`` — settings, not properties of the file. Being conservative in the other
    direction than :func:`derive_required_inputs` is deliberate; this set only ever escalates a
    warning, so over-inclusion costs a louder message and under-inclusion costs a silent one.

    Returns:
        ``{"somatic": (...), "germline": (...)}``. The union over both arms is the four columns
        issue #28 names.
    """
    return _per_arm(source, _pathogenic_read)


def _per_arm(source: str, collect) -> dict[str, tuple[str, ...]]:
    """Run ``collect`` over each arm's vendored functions and union what it finds, sorted.

    The two derivations above differ only in *which* columns of a function they want; the
    walk from ``source`` to "one sorted tuple per arm" is the same both times, and was written
    out twice before this. ``collect`` takes ``(function, frame_parameter_names)``.
    """
    functions = _functions_by_name(source)
    return {
        arm: tuple(
            sorted(
                set().union(
                    *(
                        collect(functions[name], frame_parameters(functions[name]))
                        for name in names
                    )
                )
            )
        )
        for arm, names in ARM_FILTERS.items()
    }


def _functions_by_name(source: str) -> dict[str, ast.FunctionDef | ast.AsyncFunctionDef]:
    """Every function ``source`` defines, by name, checked to cover what the app calls."""
    functions = {node.name: node for node in functions_in(ast.parse(source))}
    missing = sorted(
        name for names in ARM_FILTERS.values() for name in names if name not in functions
    )
    if missing:
        raise RuntimeError(
            f"the vendored filter module no longer defines {', '.join(missing)}, so which "
            "columns each arm requires cannot be derived. The app calls these functions by "
            "name; fix filters/absent_columns.ARM_FILTERS to match."
        )
    return functions


def _always_read(body: Iterable[ast.stmt], frames: set[str]) -> set[str]:
    """Columns read on every path through a statement list.

    The only interesting case is ``if``: a column read in the test is read on every path, and a
    column read in a branch is read on every path only if *every* branch reads it. An ``if``
    with no ``else`` therefore contributes nothing from its body, which is how ``Hugo_Symbol``
    and ``CIViC_Evidence_Level`` fall out — and ``CancerVar``, read in both arms of the
    ``skip_civic`` branch, survives.

    Loops and ``try`` bodies are treated as *not* guaranteed, which is conservative in the
    direction that matters: a column this misses is a column the app does not fill and the
    vendored code raises on, which is today's behaviour rather than a new silent verdict. The
    vendored filters contain neither today.
    """
    always: set[str] = set()
    for statement in body:
        if isinstance(statement, ast.If):
            always |= _columns_in(statement.test, frames)
            always |= _always_read(statement.body, frames) & _always_read(
                statement.orelse, frames
            )
        elif isinstance(statement, (ast.For, ast.AsyncFor, ast.While, ast.Try, ast.With)):
            continue
        else:
            always |= _columns_in(statement, frames)
    return always


def _pathogenic_read(function: ast.FunctionDef, frames: set[str]) -> set[str]:
    """Columns reaching any assignment to :data:`PATHOGENIC_MASK` in ``function``.

    The mask is built by straight-line assignment in every branch of both arms, so reading the
    assignments is enough; nothing here follows a name through a later rebinding. A vendored
    body that started composing the rescue some other way would derive an empty set, which
    :data:`PATHOGENIC_INPUTS`'s import-time check turns into a startup failure rather than an
    escalation that quietly stops firing.
    """
    found: set[str] = set()
    for node in ast.walk(function):
        if isinstance(node, ast.Assign) and any(
            isinstance(target, ast.Name) and target.id == PATHOGENIC_MASK
            for target in node.targets
        ):
            found |= _columns_in(node.value, frames)
    return found


def _columns_in(node: ast.AST, frames: set[str]) -> set[str]:
    """Every ``maf["literal"]`` anywhere under ``node``.

    A full walk, unlike ``numeric_columns._numeric_operands``, which descends only through
    arithmetic. The question there is what a column's values are *compared as*, so a call in
    between changes the answer; the question here is only whether the column is looked up at
    all, and ``maf["X"].apply(f)`` raises ``KeyError`` on an absent ``X`` exactly as
    ``maf["X"] > 0`` does.
    """
    return {
        child.slice.value
        for child in ast.walk(node)
        if isinstance(child, ast.Subscript) and is_frame_column(child, frames)
    }


#: The columns whose absence makes the vendored call raise, per arm. Derived at import, once —
#: a source the parser cannot make sense of is then a startup failure with a clear cause rather
#: than an app that silently stops filling anything.
REQUIRED_INPUTS: dict[str, tuple[str, ...]] = derive_required_inputs(vendored_source())

#: The columns that feed pathogenic retention, per arm. The union across arms is four:
#: ``CancerVar``, ``ClinVar_VCF_CLNSIG``, ``CIViC_Evidence_Level``, ``InterVar``.
PATHOGENIC_INPUTS: dict[str, tuple[str, ...]] = derive_pathogenic_inputs(vendored_source())

_undecided = sorted(
    {column for columns in REQUIRED_INPUTS.values() for column in columns}
    - set(NEUTRAL_FILLS)
)
if _undecided:  # pragma: no cover - a pipeline change, caught at import
    raise RuntimeError(
        "the vendored filters now read "
        f"{', '.join(_undecided)} on every path, and no neutral fill is recorded for "
        "them. What value stands in for a missing annotation is a clinical judgement with a "
        "measured cost, not something this module can guess: add an entry to "
        "filters/absent_columns.NEUTRAL_FILLS, with the direction it moves the report."
    )

if not any(PATHOGENIC_INPUTS.values()):  # pragma: no cover - a pipeline change
    raise RuntimeError(
        "no pathogenic-retention columns could be derived from the vendored filters. Either "
        f"the pipeline no longer builds a {PATHOGENIC_MASK!r} mask — in which case the "
        "escalated warning has nothing left to warn about — or the derivation has stopped "
        "understanding the source, in which case the loudest warning the app has would "
        "silently never fire again."
    )


# ---------------------------------------------------------------------------
# The plan
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class FillPlan:
    """What this MAF is missing, what stands in for it, and what to tell the user.

    Built once per filter call and consulted three times — for the frame the masks are
    computed on, for the warnings, and for the structured fields on
    ``filters.variant_filters.Diagnostics``. One object so those three cannot come to disagree
    about which columns were filled.
    """

    arm: str
    #: The absent required columns, in derived order. Empty on any MAF the pipeline would
    #: accept, which is the case worth keeping cheap.
    filled: tuple[str, ...] = ()
    #: The subset of :attr:`filled` that feeds pathogenic retention. Escalated, not merely
    #: reported.
    degraded: tuple[str, ...] = ()

    def frame_for_masks(self, maf: pd.DataFrame) -> pd.DataFrame:
        """``maf`` with the neutral fills, for computing masks and nothing else.

        Returns ``maf`` **itself** when nothing is missing. That identity is the whole
        performance story of this module: the happy path is every MAF the pipeline would
        accept, and a 427-column frame duplicated on every re-filter click to add no columns
        would be a real cost paid for a rare case. ``tests/test_absent_columns.py`` asserts the
        identity and measures the allocation on a large frame.

        On the unhappy path it *is* a copy, deliberately and unavoidably: ``assign`` returns a
        new frame, which is the point. The caller labels the original, so the fill reaches the
        masks and nothing else — a ``tumor_f`` of 1.0 in the grid would be the app lying to the
        user about their own data, which is the failure that made "fill and carry on" nearly
        unacceptable in the first place. The missing column stays missing everywhere the user
        can see it.
        """
        if not self.filled:
            return maf
        return maf.assign(**{column: NEUTRAL_FILLS[column].value for column in self.filled})

    def warnings(self) -> tuple[Note, ...]:
        """What to tell the user: the escalations first, then the ordinary dropout line.

        Escalations first because they are the ones that mean *rows are missing from your
        report*, and they are per column rather than pooled — there are at most two, and each
        names a specific clinical consequence that a joined sentence would blur.

        The ordinary line is a single message covering every other filled column, grouped by
        which way it moved the report. One message and not one per column: issue #20 measured
        the app emitting eight undifferentiated yellow boxes per arm, and a warning nobody
        finishes reading is not a warning.

        Every note here is at least a ``WARNING``: this class only speaks when a filter input
        was absent and a stand-in was used, which is by definition a report that is not the one
        the user asked for. Nothing it produces is ever ``INFO``, and
        ``tests/test_filter_notes.py`` asserts that rather than leaving it to be read off the
        two builders.
        """
        notes = [_degraded_note(column) for column in self.degraded]
        ordinary = [column for column in self.filled if column not in self.degraded]
        if ordinary:
            notes.append(_dropout_note(ordinary))
        return tuple(notes)


def plan_fills(maf: pd.DataFrame, arm: str, skip_pathogenic: bool = False) -> FillPlan:
    """Which of ``arm``'s required inputs ``maf`` does not carry.

    Args:
        maf: any frame. Only its column names are read.
        arm: ``"somatic"`` or ``"germline"``.
        skip_pathogenic: the pipeline's flag, in the pipeline's polarity. When it is set the
            rescue mask is already all-``False``, so filling a pathogenic-retention column
            **cannot** degrade a safety net that the user has switched off — the fill costs
            guideline rows and nothing more, and those get the ordinary dropout line.

            Escalating anyway would be the more cautious-looking choice and the wrong one. The
            escalated warning says *variants that would have been retained as pathogenic are
            absent from this report*, and under ``skip_pathogenic`` no variant would have been
            retained: the sentence names a loss that did not happen, in the app's loudest
            voice. A warning that cries wolf on a setting the user chose deliberately is how
            the next real one gets ignored.

    Raises:
        KeyError: for an arm with no derived contract. The caller has already refused an
            unknown arm by the time this runs; a silent empty plan would mean an unknown arm
            filled nothing and filtered anyway.
    """
    required = REQUIRED_INPUTS[arm]
    present = set(maf.columns)
    filled = tuple(column for column in required if column not in present)
    escalated = set() if skip_pathogenic else set(PATHOGENIC_INPUTS[arm])
    return FillPlan(
        arm=arm,
        filled=filled,
        degraded=tuple(column for column in filled if column in escalated),
    )


def _degraded_note(column: str) -> Note:
    """The escalated warning, for a column that feeds pathogenic retention.

    Says which way the fill moves the report, because a user reading a shorter one has no way
    to tell a strict filter from a disarmed safety net, and the two call for opposite
    responses. It holds for every column that reaches here: all three are :data:`REMOVES`.

    It says it as a property of the **fill** — *a stand-in value can only take rows out* —
    rather than as the bare comparison *"rows are missing, not added"* it used to make. That
    phrasing named no referent, and ``components/arm_notice.py`` had already recorded the
    collision it caused and handed it here. Directly above this box the mismatch notice
    prices the *other arm*, and on the germline reference filtered as somatic it reports 21
    kept against 31, with 14 in both — so seven variants are on screen that the correct arm
    would not have returned. *Rows added*, in the only comparison the user can see, under a
    sentence saying rows are not added. Both were true of their own counterfactual — this
    one's is the same run with a real column — and nothing on screen said which. Naming the
    fill as the subject removes the collision without the module learning anything about arms.

    It also carries the not-a-complete-result statement, rather than leaving that to the dropout
    line. An escalated column is *excluded* from the pooled line — that is what makes the
    escalation distinct — so a MAF whose only missing column is this one would otherwise be told
    its safety net is gone and never told the report cannot be read as a whole one.

    What issue #136 rewrote, and why the old mechanism clause was false
    ------------------------------------------------------------------
    It said the fill *"does not just drop one term from an OR — it collapses the safety net"*.
    Read against ``vendor/pipeline_filters.py`` that is exactly backwards. ``filter_patho`` **is
    an OR**, over two sources or three: somatic ``CancerVar | ClinVar_VCF_CLNSIG`` — plus
    ``CIViC_Evidence_Level`` only where that column exists, which on the reference is 0 of 100
    files — and germline ``InterVar | ClinVar_VCF_CLNSIG``. Filling one member neutral drops
    one term and the rest survive; nothing collapses. The distinction the sentence was reaching
    for is real but is about **which** OR, not about whether a term was dropped: this one runs
    beside the criteria rather than inside them, so the same algebraic move costs far more.

    The overstatement was not free, because ``ClinVar_VCF_CLNSIG`` escalates on *both* arms and
    drew that same paragraph. Against the module docstring's #20 baselines its fill costs
    388 → 334 somatic and 496 → 475 germline — −14% and −4% — while ``CancerVar``'s costs
    388 → 100 and ``InterVar``'s 496 → 90. One sentence saying "collapses" covered a −4% loss.
    The replacement states the mechanism and not the magnitude, so it is true at both ends, and
    no per-column severity is invented to be kept in step with the reference.

    *"...keeps variants that did not meet your filter criteria"* is lifted from what the run
    report already renders — ``_decomposition_summary``'s *"kept only by pathogenic retention
    (these did not meet the criteria)"* — so the warning and the report tell one story rather
    than a fourth. It replaces the word *unconditional*, which over-claimed in a second way:
    the rescue clears the criteria, but the app's own population-frequency filter is applied to
    it as well, exempted only by the narrower ClinVar-only :func:`pathogenic_exemption`.
    """
    return Note.error(
        f"❌ **{ESCALATION_MARKER} — `{column}` column not found.** This is one of the "
        "sources the pathogenic-retention rule reads, and that rule keeps variants that did "
        "not meet your filter criteria. With a stand-in value it matches nothing, so variants "
        f"retained on `{column}` evidence alone are **absent from this report**. A stand-in "
        "value can only take rows out of your report, never put them in. "
        + _not_a_complete_result()
    )


def _dropout_note(columns: Sequence[str]) -> Note:
    """The ordinary warning: what was filled, and which way each fill moved the report.

    Both directions are named even when only one occurred, in the sense that the sentence
    always says which one this is. "Filled neutrally" alone reads as *nothing changed*, which
    is true of the parameter and false of the report.

    A ``Sequence`` rather than an ``Iterable``, because the columns are read twice — once to
    group them by direction and once to name them all — and a generator would silently produce
    a message naming nothing.

    Two things issue #136 changed here, both consequences of the shared sentence above rather
    than decisions of their own. The opening no longer ends *"so filtering could run"*: that
    said the same thing :func:`_not_a_complete_result` now says in the very next sentence, and
    it only read as new information while the sentence beside it was claiming the file would
    have been refused. And the clause verbs agree with what they are about — this note is
    pooled, so it is reached with one column as readily as with five, and *"``RENOVO_Class``
    ... now match nothing"* was the commonest single wrong-arm rendering the app had. It is
    what a somatic MAF filtered on the germline arm draws beside the escalated ``InterVar``
    note, in 110 of 110 somatic-detected files in the measured corpus.
    """
    by_direction: dict[str, list[str]] = {}
    for column in columns:
        by_direction.setdefault(NEUTRAL_FILLS[column].direction, []).append(column)

    clauses = []
    if REMOVES in by_direction:
        removed = by_direction[REMOVES]
        clauses.append(
            f"{_names(removed)} fed the guideline filter, and now "
            f"{_matches(removed)} nothing — so rows are **missing** from this report, and a "
            "missing row is not visible in the table"
        )
    if ADDS in by_direction:
        added = by_direction[ADDS]
        clauses.append(
            f"{_names(added)} gated on values this MAF does not carry, and now "
            f"{_passes(added)} every row — so this report may contain rows that would "
            "otherwise have been dropped, which are at least visible in the table"
        )

    return Note.warning(
        "⚠️ **Missing filter-input column(s) filled with a stand-in value**: "
        f"{_names(columns)}. "
        + _not_a_complete_result()
        + " "
        + "; ".join(clauses)
        + ". The stand-in was used only to compute the filter: the columns stay absent in the "
        "table and in the export."
    )


def _derive_civic_guard_column(source: str) -> str:
    """The column the pipeline's own CIViC guard tests for.

    Derived rather than written down, for the reason this module derives everything else:
    a hand-written ``"CIViC_Evidence_Level"`` here would be one more copy of a pipeline
    decision, in a codebase whose parity problem began with copies someone made correctly
    and then stopped updating.

    Found by shape rather than by name — the ``<str> in <frame>.columns`` test inside
    ``check_civic_column_exists`` — so it follows the guard if the column is renamed, and
    raises here rather than silently reporting "not applied" if the guard is restructured
    into something this cannot read.
    """
    for function in functions_in(ast.parse(source)):
        if function.name != "check_civic_column_exists":
            continue
        for node in ast.walk(function):
            if (
                isinstance(node, ast.Compare)
                and len(node.ops) == 1
                and isinstance(node.ops[0], ast.In)
                and isinstance(node.left, ast.Constant)
                and isinstance(node.left.value, str)
            ):
                return node.left.value
    raise RuntimeError(
        "cannot find the column check_civic_column_exists guards on, so the app cannot "
        "say whether the CIViC clause ran"
    )


#: The column whose absence makes the pipeline skip CIViC entirely.
CIVIC_GUARD_COLUMN = _derive_civic_guard_column(vendored_source())


def civic_clause_applies(skip_civic: bool, available_columns: Iterable[str]) -> bool:
    """Whether the CIViC clause ran at all, on the pipeline's own two conditions.

    ``somatic_filters`` opens ``if skip_civic or not civic_column_exists`` and drops the
    CIViC term from **both** the criteria OR and the pathogenic-retention OR. So a report
    filtered from a MAF without that column was not cut on CIViC at any point, whatever the
    user's CIViC selection says.

    This is not an edge case and that is why it is worth a function: ``CIViC_Evidence_Level``
    is absent in **100 of 100** reference files (see this module's header), so a surface
    that reports the user's CIViC selection as a restriction that applied would be wrong on
    essentially every MAF this app has been pointed at. It is the one filter input the
    pipeline guards for itself, which is why nothing fills it and why
    :func:`filled_input_note` never mentions it.
    """
    return not skip_civic and CIVIC_GUARD_COLUMN in set(available_columns or ())


def filled_input_note(
    filled: Sequence[str], degraded: Sequence[str] = ()
) -> str | None:
    """The same fact as the warnings above, short enough to sit inside a filter recap.

    ``None`` when the file carried every filter input, which is the ordinary case: on the
    measured corpus the escalated warning fires on the file's own arm for 2 of 173 placeable
    MAFs (issue #136).

    **Why a second spelling of one fact is not a second copy of it.** The warnings above are
    drawn once, by the render that follows a filter run, and they are gone from the screen
    after it — a user who scrolls back to a report an hour later has no record that a column
    was filled at all. This one is carried with the report by
    :class:`~filters.run_recap.RunRecap` and is on screen for as long as the report is. What
    keeps the two from drifting is that the sentence they share is
    :func:`_not_a_complete_result`, spelled here as it is there, and that both take their
    column lists from the same :class:`~filters.variant_filters.Diagnostics`.

    It is deliberately shorter than the warnings rather than a summary of them. The warning
    has to explain the mechanism to a user meeting it for the first time; this one appears
    beside the settings that produced the report, where the question it answers is the
    narrow *is this report whole?*

    A plain string and **not** a :class:`~filters.notes.Note`, which is not an oversight: the
    recap draws this itself, in its own place, and never routes it through the data page's slot.
    A level would say which of three boxes to draw it in, and there are no boxes here.
    """
    if not filled:
        return None

    note = (
        f"⚠️ This MAF did not carry {_names(filled)}, filled with a stand-in value to "
        f"filter it. {_not_a_complete_result()}"
    )
    if degraded:
        # The one clause worth a second sentence: a fill that reaches pathogenic retention
        # means rows are *absent*, which is the only thing here the user cannot see for
        # themselves by reading the table.
        note += (
            f" {_names(degraded)} also {_feeds(degraded)} the pathogenic-retention rule, "
            "so variants it would have kept are missing."
        )
    return note


def _feeds(columns: Sequence[str]) -> str:
    """``feeds`` for one column, ``feed`` for several. The sibling of :func:`_matches`."""
    return "feeds" if len(columns) == 1 else "feed"


def _matches(columns: Sequence[str]) -> str:
    """``matches`` for one column, ``match`` for several.

    The clause's subject is :func:`_names`' list, so the verb has to follow how many columns
    landed in *this* direction rather than how many were filled — a MAF missing one guideline
    source and two quality gates puts one name in front of this verb and two in front of
    :func:`_passes`.
    """
    return "matches" if len(columns) == 1 else "match"


def _passes(columns: Sequence[str]) -> str:
    """``passes`` for one column, ``pass`` for several. The sibling of :func:`_matches`."""
    return "passes" if len(columns) == 1 else "pass"


def _names(columns: Sequence[str]) -> str:
    """Column names for a message, always spelled out.

    Never truncated with an "and N more": the derived contract is seven columns per arm, so the
    longest possible list is seven names, and a user whose MAF is that incomplete needs all of
    them at once rather than one per reload.
    """
    return ", ".join(f"`{column}`" for column in columns)


