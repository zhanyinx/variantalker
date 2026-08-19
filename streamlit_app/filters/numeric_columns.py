"""Which MAF columns the pipeline reads as numbers, and what to do when they aren't.

Two things live here, and they are one thing really: a **contract read out of the
pipeline's own source**, and the refusal that enforces it.

The contract
------------
``bin/filter_variants.py`` compares three columns arithmetically —
``(t_alt_count + t_ref_count) >= coverage`` and ``tumor_f > vaf`` — and reads every other
column it touches as text, through ``.isin``, ``.apply`` or ``.str.upper``. When one of the
three cannot be read as a number, the comparison raises ``TypeError`` and the pipeline dies
with a traceback. That is the pipeline's non-verdict, and the app's job is to reach the same
non-verdict with a message a clinician can act on.

That list of three is **derived, not written down**. :data:`NUMERIC_COLUMNS` is computed at
import by parsing ``vendor/pipeline_filters.py`` — the very module the app calls — and
collecting the columns whose own values reach an arithmetic or ordering operator. Writing
the three names into this file would have been shorter by twenty lines and would have been
a fourth copy of a pipeline decision, which is the failure mode the whole vendoring effort
exists to prevent: the app's re-implementation of this filter drifted from the pipeline in
twelve distinct ways, and every one of them started as a list someone copied correctly.

Why *numeric context* rather than *every mention*
-------------------------------------------------
The distinction is load-bearing in both directions.

``compute_keep`` names ``gnomAD_exome_AF`` and its three siblings, but only to ask whether
they are in ``out.columns``. Every committed fixture carries ``.`` in a frequency column —
MAF's own blank marker — so a derivation that collected every mention would refuse all seven
of them, including the six the pipeline filters happily. And ``maf["A"].apply(f) > 0``
compares the *result of a call*, not the column, so following the compare into the call would
refuse a text column for being text.

Conversely, a plain dtype check would be fussier than the pipeline in the other direction:
``pd.Series([1.0, 2.0], dtype=object) > 0.5`` is a perfectly good mask, and object dtype
alone is not a reason to refuse anything. The trigger is a **value that is not a number**.

The biconditional
-----------------
:func:`require_readable_numerics` refuses **if and only if** the vendored code would raise.
``tests/test_numeric_columns.py`` asserts that over all seven fixtures and over injected
values covering both arms and both answers, with the vendored functions themselves as the
oracle. The rule that makes it hold is exactly:

    a present derived column is unreadable when it holds an entry that is neither
    missing nor a Python number

which tracks pandas' own behaviour cell for cell — ``"." > 0.5`` and ``"5" + 3`` both raise,
while ``None``, ``NaN``, ``pd.NA``, ``bool`` and ``Decimal`` all compare and add fine. Note
what that means for a numeric *string*: it coerces cleanly and is still a refusal, because
the pipeline concatenates it rather than adding it. An implementation built on
``pd.to_numeric(errors="coerce")`` gets that case wrong in the dangerous direction — it
returns a confident verdict on a file the pipeline crashes on.

What this deliberately does not own
-----------------------------------
* **Absent columns.** A derived column that is not in the frame is skipped. The pipeline
  raises ``KeyError`` there, and what the app should do about it — a neutral fill, an
  escalated warning — is issue #39's ruling, with its own measurements. Annexing it here
  would decide it by accident.
* **Frequency columns.** They are coerced to missing by the app's own frequency mask and
  are never refused, because a blank in a population panel means *not seen in this panel*
  rather than *this call cannot be assessed*. The carve-out is not a special case in this
  module — it falls out of the derivation, since the pipeline never compares them.
* **Row-level rescue.** There is none available. ``read_maf`` fixes a column's dtype once,
  when it reads the file, so one bad cell leaves the whole column object-typed and dropping
  the offending rows afterwards leaves it that way. The refusal is whole-file because the
  damage is whole-column.
"""

from __future__ import annotations

import ast
import inspect
import numbers
from typing import Mapping

import pandas as pd

from filters.vendored_ast import frame_parameters, functions, is_frame_column
from vendor import pipeline_filters as _pipeline_filters

#: Operators that force a numeric reading of both sides. Ordering only: ``==`` and ``in``
#: compare object-wise and never raise, which is why ``maf["CancerVar"] == x`` is not a
#: numeric context even though it is a comparison.
_ORDERING_OPS = (ast.Lt, ast.LtE, ast.Gt, ast.GtE)

#: Arithmetic operators, for a column that reaches one without being compared — and for
#: ``(t_alt_count + t_ref_count)``, which is where the depth columns are actually found.
_ARITHMETIC_OPS = (
    ast.Add,
    ast.Sub,
    ast.Mult,
    ast.Div,
    ast.FloorDiv,
    ast.Mod,
    ast.Pow,
    ast.MatMult,
)

class UnreadableNumericColumns(ValueError):
    """A MAF whose depth or VAF columns cannot be read as numbers.

    Raised rather than reported through :class:`~filters.variant_filters.Diagnostics`
    because there is **no verdict to accompany**. A warning describes a report that exists;
    here the app cannot produce one, and a frame handed back alongside a warning would have
    to contain a decision nothing was able to make. The pipeline's answer to this file is a
    traceback; the app's is this exception, which is the same non-verdict with a message.

    A ``ValueError`` subclass so that a caller which has not heard of this type still
    handles it as bad input rather than as a bug — the Streamlit call site catches broadly,
    and the message is what the user needs either way.

    Attributes:
        offenders: column name -> ``{rendered value: row count}``, for a caller that wants
            to present the refusal its own way instead of parsing :func:`str`.
    """

    def __init__(self, offenders: Mapping[str, Mapping[str, int]]):
        self.offenders = {column: dict(values) for column, values in offenders.items()}
        super().__init__(_message(self.offenders))


# ---------------------------------------------------------------------------
# Reading the contract out of the vendored source
# ---------------------------------------------------------------------------


def vendored_source() -> str:
    """The source text of the vendored filter module the app calls.

    Read through the **import system**, not off a path this module constructs. Two reasons,
    and the second is the one that bites:

    * it is the module actually imported, so the contract cannot be derived from some other
      copy of the filter code — a stale file, an installed egg, or ``bin/``, which no
      installer ships;
    * a packaged build has no ``bin/`` and may have no ``__pycache__``, and asking the
      loader for source works whether the module came from a plain ``.py``, a zip import, or
      a tree that has never been byte-compiled.

    The loader is asked *first* because it answers the question actually being asked — "what
    source did this module come from" — while ``inspect.getsource`` re-derives a *filename*
    and reads that, which is one more step that can point somewhere else. It is a fallback
    rather than an error path: measured under zipimport, both routes return the same text
    (``tests/test_numeric_columns.py`` says so and how it was checked), so the ordering is a
    preference for directness and not a claim that ``inspect`` fails anywhere tested.
    """
    loader = getattr(_pipeline_filters, "__loader__", None)
    get_source = getattr(loader, "get_source", None)
    if get_source is not None:
        source = get_source(_pipeline_filters.__name__)
        if source:
            return source

    try:
        return inspect.getsource(_pipeline_filters)
    except OSError as exc:  # pragma: no cover - a build that ships no source at all
        raise RuntimeError(
            "cannot read the source of "
            f"{_pipeline_filters.__name__}, so the numeric-column contract cannot be "
            "derived. The app must ship the vendored filter source, not bytecode alone."
        ) from exc


def derive_numeric_columns(source: str) -> tuple[str, ...]:
    """The MAF columns ``source`` compares or combines arithmetically, sorted.

    Args:
        source: Python source text — normally :func:`vendored_source`, and a mutated copy
            of it in the tests, which is what makes the derivation checkable as a
            derivation rather than as a remembered answer.

    Returns:
        Column names, sorted for a stable identity across runs.

    Which subscripts count is decided by :func:`_numeric_operands`; which names denote the
    MAF is decided by :func:`~filters.vendored_ast.frame_parameters`. Both are deliberately
    conservative: a shape this parser does not understand yields *nothing* rather than a
    guess, so the failure mode is the app being no fussier than it is today rather than the
    app refusing a file for a reason nobody can explain. The self-test that
    :data:`NUMERIC_COLUMNS` is non-empty guards the one way that could go silently wrong.
    """
    found: set[str] = set()
    for function in functions(ast.parse(source)):
        frames = frame_parameters(function)
        if not frames:
            continue
        for node in ast.walk(function):
            if isinstance(node, ast.Compare) and any(
                isinstance(op, _ORDERING_OPS) for op in node.ops
            ):
                for operand in [node.left, *node.comparators]:
                    found |= _numeric_operands(operand, frames)
            elif isinstance(node, ast.BinOp) and isinstance(node.op, _ARITHMETIC_OPS):
                found |= _numeric_operands(node, frames)
    return tuple(sorted(found))


def _numeric_operands(node: ast.AST, frames: set[str]) -> set[str]:
    """Columns whose **own values** reach ``node`` as a number.

    Descends through arithmetic and parentheses only. It deliberately does *not* descend
    into calls or attribute access, and that restraint is the whole difference between this
    and a naive walk: ``maf["ClinVar_VCF_CLNSIG"].apply(f) > 0`` compares whatever ``f``
    returned, and treating it as a numeric read of ``ClinVar_VCF_CLNSIG`` would refuse
    every real MAF for holding text in a text column.
    """
    if isinstance(node, ast.Subscript):
        return {node.slice.value} if is_frame_column(node, frames) else set()
    if isinstance(node, ast.BinOp):
        return _numeric_operands(node.left, frames) | _numeric_operands(
            node.right, frames
        )
    if isinstance(node, ast.UnaryOp):
        return _numeric_operands(node.operand, frames)
    return set()


#: The columns the vendored pipeline code reads as numbers. Derived at import, once.
#:
#: Import time rather than call time so that a source the parser cannot make sense of is a
#: startup failure with a clear cause, not a filter that silently validates nothing on the
#: first MAF a user loads.
NUMERIC_COLUMNS: tuple[str, ...] = derive_numeric_columns(vendored_source())

if not NUMERIC_COLUMNS:  # pragma: no cover - a parser that has stopped parsing
    raise RuntimeError(
        "no numeric columns could be derived from "
        f"{_pipeline_filters.__name__}. Either the pipeline's filters no longer compare "
        "any column arithmetically — in which case delete this module — or the derivation "
        "has stopped understanding the source it is given, in which case the app would "
        "validate nothing at all."
    )


# ---------------------------------------------------------------------------
# The refusal
# ---------------------------------------------------------------------------


def unreadable_values(maf: pd.DataFrame) -> dict[str, dict[str, int]]:
    """Every derived column of ``maf`` the pipeline could not compare, and why.

    Args:
        maf: any frame. Columns of :data:`NUMERIC_COLUMNS` that are absent are **skipped**
            — see the module docstring on why the missing-column decision is not this
            module's.

    Returns:
        ``{}`` when nothing is wrong. Otherwise one entry per refused column, in
        :data:`NUMERIC_COLUMNS` order, mapping the **distinct values that cannot be read as
        a number at all** to how many rows hold them. Distinct rather than per-row so a
        column of ten thousand ``.`` is one line rather than ten thousand.

        The mapping is **empty** when the column refuses without any single value being at
        fault — a column of numeric *strings*, which the pipeline concatenates rather than
        adds. See :func:`_message` for what is said in that case. That shape is not
        reachable from a file: a token ``pd.to_numeric`` accepts is one ``read_csv`` accepted
        first, so the column would have arrived numeric and never been refused. It is a
        hand-built frame, and ``test_no_token_that_poisons_a_column_goes_unnamed`` pins the
        gap shut over both readings so this stays a statement of fact rather than a belief.

    Which values get named is a separate question from whether the column refuses, and
    conflating them makes the message useless on real data. Refusal is decided by "is every
    present value a number?", because that is what pandas raises on. But once one cell holds
    ``.``, ``read_maf`` leaves the *whole* column as strings — so on a 180,000-row MAF
    "every value that is not a number" is 150,000 distinct numeric strings, and a message
    listing them names everything and points at nothing. What the user can act on is the
    handful of tokens that are not numbers in any reading, so those are what is named.

    A numeric dtype short-circuits: pandas has already proved every value is a number, and
    the per-cell scan below is the expensive part on a large MAF.
    """
    offenders: dict[str, dict[str, int]] = {}
    for column in NUMERIC_COLUMNS:
        if column not in maf.columns:
            continue

        series = maf[column]
        if pd.api.types.is_numeric_dtype(series):
            continue

        present = series[series.notna()]
        not_numbers = present[~present.map(_is_number)]
        if not_numbers.empty:
            # Object dtype holding nothing but real numbers and missing values. The
            # pipeline compares this frame without complaint, so refusing it would make
            # the app fussier than the pipeline on a file the pipeline reports on.
            continue

        rendered = not_numbers.astype(str)
        illegible = rendered[pd.to_numeric(rendered, errors="coerce").isna()]
        counts = illegible.value_counts()
        offenders[column] = {
            str(value): int(count)
            for value, count in sorted(
                counts.items(), key=lambda item: (-item[1], str(item[0]))
            )
        }
    return offenders


def _is_number(value: object) -> bool:
    """Whether the pipeline could compare and add ``value``.

    ``numbers.Number`` covers ``int``, ``float``, ``bool``, ``Decimal``, ``Fraction`` and
    every NumPy scalar, all of which pandas compares and adds without complaint. ``complex``
    is excluded because it is the one ``Number`` that has no ordering — ``1j > 0`` raises,
    so a complex cell is as unreadable as a string.
    """
    return isinstance(value, numbers.Number) and not isinstance(value, complex)


def require_readable_numerics(maf: pd.DataFrame) -> None:
    """Refuse ``maf`` if any derived column holds a value that is not a number.

    Raises:
        UnreadableNumericColumns: naming every offending column and every distinct
            offending value in one message.
    """
    offenders = unreadable_values(maf)
    if offenders:
        raise UnreadableNumericColumns(offenders)


def _message(offenders: Mapping[str, Mapping[str, int]]) -> str:
    """One message, naming every column and every value.

    One and not one-per-column: the columns are unreadable *together*, because they are one
    file. Told about ``t_alt_count`` alone, a user fixes it, reloads, and is told about
    ``t_ref_count`` — three round trips through a clinical drive to learn what one sentence
    could have said.

    It also says what *would* work, because the near miss is the common case: a blank cell
    and ``NA`` read as missing, the row then fails the depth or VAF gate, and that is a
    verdict rather than a crash. A user whose export writes ``.`` for missing has a
    one-column fix and no way to guess it from "unreadable".
    """
    columns = ", ".join(offenders)
    detail = "; ".join(
        f"{column}: " + _render(values) for column, values in offenders.items()
    )
    plural = "column" if len(offenders) == 1 else "columns"
    return (
        f"This MAF cannot be filtered: the {plural} {columns} hold values that are not "
        "numbers, and the depth and VAF gates compare them arithmetically — so this file "
        "is refused outright rather than filtered wrongly. One such value makes the whole "
        "column unreadable, so removing the affected rows does not help. "
        f"Offending values — {detail}. "
        "An empty cell or 'NA' would read as missing and the row would simply fail the "
        "depth or VAF threshold; '.' and other text will not."
    )


def _render(values: Mapping[str, int]) -> str:
    """One column's offending values, or a plain statement when none is at fault.

    The empty case is not a fallback for a shape that cannot happen: a column of numeric
    strings refuses correctly — the pipeline concatenates ``"5"`` and ``3`` rather than
    adding them — while no individual cell is wrong. Saying "no single value is at fault"
    is the honest report, where naming every string in the column would name the whole file.
    """
    if not values:
        return "no single value is at fault; the column holds text rather than numbers"
    return ", ".join(f"{value!r} ({count} row(s))" for value, count in values.items())
