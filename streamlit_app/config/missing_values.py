"""What *this row does not say* means, for everything the interface draws.

The display half of a question `filters/absent_columns.py` owns the filter half of. That
module decides what a missing **annotation** would have said to the filter, per column and
measured against a reference report. This one decides whether a value is worth putting on a
clinician's screen at all. They must not be folded together: a fill exists to keep the filter
running, and is never shown to anyone (``FillPlan.frame_for_masks``), while the answers here
are only ever about rendering and reach no verdict.

Two questions, not one
----------------------
Issue #131 charted this expecting one set of strings and found the app holding four. Measured
against 303 byte-distinct real MAFs (350,003 rows), the honest answer is **two** predicates,
and the second contains the first:

* :func:`is_blank` — *the cell is empty*. Nothing was written there.
* :func:`says_nothing` — *the row does not say*. Blank, **or** an annotator wrote a sentinel
  meaning it could not assign a value.

The split is load-bearing rather than tidy-minded, and the note key is what makes it so. A
variant whose ``Hugo_Symbol`` is ANNOVAR's ``Unknown`` — 2,605 rows across **136 of the 303
files** — is still a variant identified by its chromosome, position and alleles, so it can
carry a note. A variant with an *empty* key column cannot, because two rows differing only in
that component would collapse onto one key (issue #127). So ``components/variant_table``
gates its note key on :func:`is_blank` and never on :func:`says_nothing`, and applying the
wrong one there would silently withdraw notes from nearly half of real files.

Pandas has already answered the first one
-----------------------------------------
The loader is the pipeline's — ``pd.read_csv(file, header=c, sep="\\t", low_memory=False)``,
with default missing-value handling and no nullable dtypes anywhere in the app. So pandas
converts ``""``, ``None``, ``NaN``, ``NA``, ``N/A``, ``null`` and ``<NA>`` to ``NaN`` *before
the app sees the value*, and ``str(NaN)`` is ``"nan"``. Measured over the corpus, every one of
those spellings occurs **zero** times as literal text in a file. That is why the sets this
module replaces were mostly ceremony: of the six members of the old ``_NO_VALUE``, four could
not arrive through the loader at all, and ``""`` reaches a check only via
``row.get(name, "")`` when the *column* is absent.

The spellings that do occur
---------------------------
=================  ============  ==========  =============================================
spelling           rows          columns     what wrote it
=================  ============  ==========  =============================================
``.``              39,423,699    422         ANNOVAR / dbNSFP: no annotation for this row
``__UNKNOWN__``    1,789,697     87          Funcotator: value not determined
``UNKNOWN``        4,043         7           ANNOVAR: no transcript consequence
``unknown``        3,854         8           ANNOVAR: exonic function not assigned
``Unknown``        2,702         3           ANNOVAR: no gene assigned
``NONE``           76            1           ANNOVAR: ``Ref.Gene`` unassigned
``-``              204,251       18          **a real value** — see below
=================  ============  ==========  =============================================

``-`` is deliberately absent from both predicates. Every column carrying it that the app reads
means something by it: the six allele columns and ``OREF``/``OALT`` spell the absent side of an
indel that way, and ``Transcript_Strand`` means the minus strand. The one column where it may
stand for *no value* is ``MutPred_score`` (46 rows across 15 files), which this app does not
read — ``CLINGEN_SVI_TABLE`` wires no dbNSFP column for MutPred2 — and which would in any
case fall out of ``reference_scales.parse_score``'s ``float()``. The rest are ``Unnamed: N`` columns,
the artefact of a ragged row. Issue #127 declined to read ``-`` as missing in
``Reference_Allele`` and ``Tumor_Seq_Allele2``; this extends that to every column, on
measurement rather than on assumption (issue #131).

``.`` is read as missing, and that **agrees with the filter** rather than overruling it: it is
``NEUTRAL_FILLS``'s own value for ``CancerVar``, ``InterVar``, ``RENOVO_Class`` and ``ESCAT``,
where it is neutral precisely because it is the value no keep-list carries. The charting
session suspected the opposite — that hiding a ``.`` would conceal something the filter had
valued — and the fill table says otherwise.

The sentinels are matched case-insensitively, because ANNOVAR is not consistent about case:
``Unknown``, ``UNKNOWN`` and ``unknown`` all occur, in different columns, in the same corpus.
The cost is that a literal value spelled ``unknown`` would be hidden; across the corpus it
occurs only in columns that mean it as a sentinel — the three ``ExonicFunc.*``, the three
``AAChange.*``, ``HGNC_Locus_Type`` and ``Hugo_Symbol`` — and the alternative is a list of case
variants that grows every time an annotator is added. The *blank* spellings are matched
exactly, and the asymmetry is deliberate — see :data:`_BLANK`.

The failure direction
---------------------
Issue #14's principle — extra rows are visible, missing rows are not — applies to fields as
well as rows: a value read as missing leaves the screen with no warning attached. That is why
this module is small, exact, and measured, and why :func:`shown` renders a dash rather than
dropping the field's label: a reader can tell an em dash from a value that was never drawn.
"""

from __future__ import annotations

import pandas as pd

#: Spellings that mean *nothing was written in this cell*, matched **exactly**.
#:
#: Exactly, because every member is the fixed rendering of one of Python's own missing
#: objects — ``str(float("nan"))``, ``str(None)``, ``str(pd.NA)`` — and their case is not a
#: matter of opinion the way an annotator's is. ``tests/test_missing_values.py`` derives them
#: from those objects rather than trusting this list.
#:
#: ``"nan"`` does nearly all the work, being what the loader leaves behind for every missing
#: cell; ``""`` arrives from ``row.get(name, "")`` when the column is absent rather than from
#: a file. ``"None"`` and ``"<NA>"`` are kept from issue #127: neither can survive the loader,
#: but a value reaching the interface from somewhere else — session state, a default — can
#: still be one. ``"NaT"`` is here because the guard that derives this set from the objects
#: asked for it: no column the loader produces is a datetime, since the vendored reader passes
#: no ``parse_dates``, so nothing observed renders it — but it is a missing value's rendering
#: and the four sets this module replaces all lacked it.
_BLANK = frozenset({"", "nan", "NaN", "NaT", "None", "none", "<NA>"})

#: Spellings an annotator writes to say *I could not assign a value here*, matched **case
#: insensitively**, because ANNOVAR is not consistent about case: ``Unknown``, ``UNKNOWN`` and
#: ``unknown`` all occur, in different columns, in the same corpus. Every member is measured in
#: the module docstring's table; ``-`` is deliberately not one.
#:
#: ``"none"`` is here as well as in :data:`_BLANK`, and the overlap is the point: ANNOVAR
#: writes ``NONE`` in ``Ref.Gene`` for a gene it could not assign, which is a sentinel and not
#: an empty cell, while ``None`` is a Python object that reached a render. Folding the two
#: matches together would have made those indistinguishable — caught by the guard that asserts
#: a sentinel is *not* blank.
_WROTE_NOTHING = frozenset({".", "unknown", "__unknown__", "none"})


def is_blank(value) -> bool:
    """Whether this cell is empty — nothing was written there.

    Args:
        value: any cell value, of any dtype, or a missing one.

    Returns:
        bool: ``True`` when the cell holds no value at all.

    Use this where an *identity* is being assembled and a hole would collapse two rows onto
    one, which in this app means the note key alone. Everything a user reads wants
    :func:`says_nothing` instead — a sentinel is still not worth rendering, but it does still
    identify the variant it sits on.
    """
    return str(value).strip() in _BLANK


def says_nothing(value) -> bool:
    """Whether this row says anything here worth putting on screen.

    Args:
        value: any cell value, of any dtype, or a missing one.

    Returns:
        bool: ``True`` when the cell is blank *or* holds an annotator's sentinel.

    This is the app's one answer for display. It is deliberately not the filter's — see
    ``filters/absent_columns.NEUTRAL_FILLS``, which answers a different question per column
    and whose values are never rendered.
    """
    rendered = str(value).strip()
    return rendered in _BLANK or rendered.lower() in _WROTE_NOTHING


def says_nothing_over(values: pd.Series) -> pd.Series:
    """:func:`says_nothing` over a whole column, for the frame-shaped call sites.

    Args:
        values: any column, of any dtype.

    Returns:
        pandas.Series: a boolean mask, index-aligned with ``values``.

    Separate from :func:`says_nothing` because ``Series.apply`` over a MAF-sized frame is not
    free and the two call sites that need this run on every render. Reads the same two sets,
    with the same two matching rules, so the row-wise and column-wise answers cannot drift
    apart — asserted value for value in ``tests/test_missing_values.py``.
    """
    rendered = values.astype(str).str.strip()
    return rendered.isin(_BLANK) | rendered.str.lower().isin(_WROTE_NOTHING)


def shown(value, dash: str = "—") -> str:
    """The value as the interface should render it, or ``dash`` when it says nothing.

    Args:
        value: any cell value, of any dtype, or a missing one.
        dash: what to render instead. The default em dash is what the detail panel already
            passed to ``row.get`` as a default, which only ever defended against an absent
            *column* — never against a column that was present and empty.

    Returns:
        str: the value, stripped, or ``dash``.

    For the fields that are drawn unconditionally — a gene, a position, an allele pair — where
    the choice is not *whether* to draw but *what*. A field that should disappear entirely when
    empty asks :func:`says_nothing` and skips itself.
    """
    return dash if says_nothing(value) else str(value).strip()
