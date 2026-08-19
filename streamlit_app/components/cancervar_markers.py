"""The clinical markers behind a CBP criterion, resolved from CancerVar's own table.

Issue #198. CancerVar writes three columns of *supporting evidence* into the MAF beside its twelve
CBP scores — ``Therap_list`` (CBP1), ``Diag_list`` (CBP2) and ``Prog_list`` (CBP3) — and each holds
nothing but comma-separated integers::

    ,8989,8941,8931,8999,8980,8987,8930,8960

They are **0-based line offsets** into CancerVar's marker table, vendored at
``vendor/cancervar_markers.txt``. ``CancerVar.py:193`` reads that file with
``list(csv.reader(fh, delimiter="\\t"))`` and never skips a header, so index 0 *is* the
header row and index 8989 is the 8,990th line. ``:1090`` then indexes straight into the
result. Nothing about the mapping is heuristic, and it was measured rather than assumed:
across 109 real somatic MAFs, **112,328 indices resolved with 0 out of range and 0 whose
``Evidence_type`` disagreed with the criterion that cited them.**

This module is pure — no Streamlit, no pandas — so that
:mod:`components.cbp_evidence` keeps its one Streamlit seam and the resolution is
testable without a session. It reads a file, which no other ``components`` module does;
that is why it is its own module rather than more of ``cbp_evidence``.

Four things the renderer must not get wrong
-------------------------------------------

1. **A marker list is not tissue-filtered, and most of it is for other cancers.**
   ``CancerVar.py:1092`` appends every matching row to ``out_list`` *before* the
   tumour-type test at ``:1104``, which only demotes the criterion's ``level``. So the
   list spans tumour types by construction. Measured over the 94 files that carry
   ``tumor_tissue``: only **10.6%** of therapy markers match the sample's own tissue
   (10,408 of 97,861), and **80.1%** of populated cells contain no in-tissue marker at
   all. Level A markers run **781 in-tissue against 6,542 off-tissue**. A renderer that
   presents these as "the therapy for this variant" is wrong on four rows in five, which
   is why :attr:`Marker.in_tissue` exists and why the renderer groups on it.

2. **A populated list does not imply a non-zero score.** ``level`` is reset *inside* the
   per-transcript loop while ``out_list`` accumulates across all of them (issue #185
   §4), so the last transcript parsed decides the score. Rare but real: 45 of 9,758
   populated ``Therap_list`` cells and 103 of 3,905 ``Prog_list`` cells sit behind a
   criterion scored 0.

3. **The response is not decoration.** ``Resistance`` (2,475), ``Resistant`` (2,447),
   ``no benefit`` (2,841) and ``Poor Outcome`` (5,289) are ~11k of the ~101k resolved
   rows. A drug named without its response inverts the claim on those.

4. **The evidence level is CancerVar's ``Evidence_level`` (column 10), not ``Level``
   (column 6).** Column 6 holds the source's own free text — ``"Phase II"``,
   ``"Preclinical - Cell line xenograft"`` — while column 10 is the A/B/C/D that
   ``CancerVar.py:1093`` actually grades on. Overall distribution on resolved rows:
   D 40,243 · C 40,206 · B 22,163 · A 9,716, so most markers are the weaker two.
"""

import csv
import re
from functools import lru_cache
from pathlib import Path
from typing import NamedTuple, Optional

from config.missing_values import says_nothing

#: The vendored marker table. Under ``vendor/`` and not ``resources/`` because it is a
#: copy of the pipeline's file under a drift guard — see ``vendor/README.md``. Resolved
#: from this module's own location so the packaged app finds it too.
MARKERS_PATH = Path(__file__).resolve().parent.parent / "vendor" / "cancervar_markers.txt"

#: Column offsets in CancerVar's marker table, from its header row. Named because
#: ``row[10]`` and ``row[6]`` are both "level" to a reader and only one of them is the
#: one CancerVar grades on — see the module docstring.
_GENE = 0
_CANCER = 3
_DRUG = 4
_RESPONSE = 5
_PMIDS = 7
_CANCERVAR_TYPE = 8
_EVIDENCE_TYPE = 9
_EVIDENCE_LEVEL = 10

#: How many columns a usable marker row has. A shorter row is a truncated file, not a
#: marker; :func:`_load` drops it rather than reading ``None`` into a drug name.
_MIN_COLUMNS = 11

#: Which criterion each MAF column stands behind, and the ``Evidence_type`` its rows must
#: carry. The ``Evidence_type`` is the load-bearing half: it is what makes a resolved row
#: *checkable* rather than merely present, and a table of the wrong vintage announces
#: itself here — a shifted index lands on a row whose type no longer matches.
#:
#: ``CancerVar.py`` filters on exactly this in ``check_Thera`` (``:1091``),
#: ``check_Diagno`` (``:1238``) and ``check_Progno`` (``:1386``).
CRITERION_COLUMNS = (
    ("CBP1", "Therap_list", "Therapeutic", "therapy marker"),
    ("CBP2", "Diag_list", "Diagnostic", "diagnostic marker"),
    ("CBP3", "Prog_list", "Prognostic", "prognostic marker"),
)

#: CancerVar's A/B/C/D evidence levels, strongest first, and what each is worth to the
#: score. ``CancerVar.py:1095-1102`` grades A and B alike (both reach 2) and C and D
#: alike (both reach 1), so the pairs are not four distinct strengths.
_LEVEL_ORDER = ("A", "B", "C", "D")

#: A response that *reverses* the claim a drug name makes on its own. Matched by
#: substring, lowercased, because the table spells the same idea several ways
#: (``Resistance`` / ``Resistant`` / ``no benefit`` / ``Poor Outcome``).
_ADVERSE_RESPONSES = ("resist", "no benefit", "poor outcome", "no response", "decreased")

#: How *this table* spells "nothing here", which is not how a MAF cell spells it.
#: :func:`config.missing_values.says_nothing` is the app's measured answer for MAF cells
#: and is deliberately not widened: ``N/A`` is not a MAF sentinel, and ``-`` is
#: explicitly excluded there because it is a real value in some columns.
#:
#: CancerVar's table has its own two conventions, and they pair exactly — of 10,442
#: marker rows, **1,119 have an empty drug and 119 have the literal ``N/A``**, and those
#: same 119 carry ``not applicable`` as their response. That pair is a marker with no
#: drug at all (a diagnostic or prognostic one), not a marker whose drug is unknown, so
#: rendering either string would be noise dressed as data.
_TABLE_BLANKS = frozenset({"", ".", "-", "n/a", "na", "not applicable", "none", "null"})


def _table_value(value) -> str:
    """A marker-table cell as text, or ``""`` where the table says nothing there."""
    text = str(value).strip() if value is not None else ""
    if text.lower() in _TABLE_BLANKS or says_nothing(text):
        return ""
    return text

_DIGITS_RE = re.compile(r"\d+")


class Marker(NamedTuple):
    """One row of CancerVar's marker table, as cited by one MAF cell.

    ``in_tissue`` is the app's own answer, not the file's: it replays CancerVar's test
    (``re.findall(cancer_type, row[8], flags=re.IGNORECASE)`` — the tissue string used as
    a **regular expression**) so the renderer can group the sample's own tumour type
    first. ``None`` where the sample's tissue is unknown, which is a third state and not
    a ``False``: 15 of 109 CancerVar-bearing files carry no ``tumor_tissue`` column, and
    on those "not in this tumour type" would be a claim the app cannot make.
    """

    index: int
    gene: str
    cancer: str
    drug: str
    response: str
    level: str
    pmids: tuple
    in_tissue: Optional[bool]

    @property
    def is_adverse(self) -> bool:
        """Whether the response reverses what the drug name implies on its own."""
        lowered = self.response.lower()
        return any(word in lowered for word in _ADVERSE_RESPONSES)


class MarkerSet(NamedTuple):
    """Every marker one MAF cell cites, plus what could not be resolved.

    ``cited`` is the count taken from the MAF cell itself, so it is knowable even when
    the table is unreadable — which is what lets the renderer say *"7 therapy markers,
    which cannot be named"* rather than falling silent. Issue #187 settled that absence
    is named per state; an unresolvable marker is one of those states.
    """

    criterion: str
    noun: str
    cited: int
    markers: tuple
    unresolved: int
    table_missing: bool

    @property
    def best_level(self) -> Optional[str]:
        """The strongest ``Evidence_level`` among the resolved markers, or ``None``."""
        levels = [m.level for m in self.markers if m.level in _LEVEL_ORDER]
        if not levels:
            return None
        return min(levels, key=_LEVEL_ORDER.index)

    @property
    def in_tissue(self) -> tuple:
        """Markers matching the sample's own tumour type, strongest evidence first."""
        return tuple(m for m in self.markers if m.in_tissue)

    @property
    def other_tissue(self) -> tuple:
        """Markers for other tumour types — or all of them, where the tissue is unknown."""
        return tuple(m for m in self.markers if not m.in_tissue)

    @property
    def tissue_known(self) -> bool:
        """Whether the sample's tumour type was available to compare against."""
        return any(m.in_tissue is not None for m in self.markers)


@lru_cache(maxsize=1)
def _load(path: str) -> tuple:
    """The marker table as a tuple of rows, read once per process.

    Cached on the path rather than on nothing so a test can point at a fixture and get
    its own entry. Returns ``()`` when the file is absent — the caller turns that into
    :attr:`MarkerSet.table_missing`, because "no table" and "no markers" are different
    facts and the renderer says different things about them.

    Rows shorter than :data:`_MIN_COLUMNS` are kept as ``None`` rather than dropped: the
    index of every *later* row is the whole point, so removing a malformed line would
    silently re-base the rest.
    """
    handle = Path(path)
    if not handle.is_file():
        return ()
    with handle.open(newline="", encoding="utf-8", errors="replace") as fh:
        return tuple(
            tuple(row) if len(row) >= _MIN_COLUMNS else None
            for row in csv.reader(fh, delimiter="\t")
        )


def _cited_indices(cell) -> tuple:
    """The distinct indices in a MAF marker cell, in ascending order.

    Real cells begin with a comma (``,8989,8941``) because ``CancerVar.py:1092`` prepends
    each hit to a string that already ends in one, and then de-duplicates through a
    ``set`` — which keeps the empty leading field. Ascending rather than the file's own
    order because that order is a ``set``'s, i.e. arbitrary.
    """
    if says_nothing(cell):
        return ()
    return tuple(sorted({int(match) for match in _DIGITS_RE.findall(str(cell))}))


def _matches_tissue(tissue: str, cancervar_type: str) -> Optional[bool]:
    """CancerVar's own tumour-type test, replayed.

    ``check_Thera`` uses the tissue string as a **regular expression** against the
    table's ``Cancervar_type`` column, case-insensitively. Replayed here rather than
    re-invented so the grouping agrees with the demotion CancerVar applied to the score.

    Returns ``None`` where the sample's tissue is unknown. A tissue that is not a valid
    regular expression also yields ``None`` — the app declines to guess rather than
    reporting every marker as off-tissue, which would read as a finding.
    """
    if not tissue:
        return None
    try:
        return bool(re.findall(tissue, cancervar_type, flags=re.IGNORECASE))
    except re.error:
        return None


def resolve_markers(get, tissue: str = "", path=None) -> dict:
    """Resolve every marker cell on one variant row.

    Args:
        get: a callable taking a column name and returning that cell, or ``None``.
            Deliberately not a ``pandas.Series``: this module stays pandas-free, and the
            only thing it needs of a row is ``row.get``.
        tissue: the sample's own tumour type, from the ``tumor_tissue`` column. Empty
            where the file has no such column, which makes every ``in_tissue`` ``None``.
        path: the marker table to read. Defaults to :data:`MARKERS_PATH`.

    Returns:
        dict: ``{criterion_code: MarkerSet}``, holding only the criteria whose column is
        present *and* cites at least one index. A criterion absent from this dict has
        nothing to disclose; one present with ``markers=()`` has something to disclose
        that could not be named.

    Markers are ordered strongest evidence first, then by index so the order is stable.
    """
    table = _load(str(path or MARKERS_PATH))
    resolved = {}

    for code, column, evidence_type, noun in CRITERION_COLUMNS:
        indices = _cited_indices(get(column))
        if not indices:
            continue

        markers = []
        unresolved = 0
        for index in indices:
            row = table[index] if 0 <= index < len(table) else None
            # A row that is absent, malformed, or of the wrong Evidence_type is all one
            # thing to the reader: this index did not resolve to a marker for this
            # criterion. The wrong-type case is the one that catches a bad vintage.
            if row is None or row[_EVIDENCE_TYPE] != evidence_type:
                unresolved += 1
                continue
            markers.append(
                Marker(
                    index=index,
                    gene=_table_value(row[_GENE]),
                    cancer=_table_value(row[_CANCER]),
                    drug=_table_value(row[_DRUG]),
                    response=_table_value(row[_RESPONSE]),
                    level=_table_value(row[_EVIDENCE_LEVEL]),
                    pmids=tuple(
                        piece
                        for piece in (
                            _table_value(p) for p in re.split(r"[;,]", row[_PMIDS])
                        )
                        if piece.isdigit()
                    ),
                    in_tissue=_matches_tissue(tissue, row[_CANCERVAR_TYPE]),
                )
            )

        markers.sort(
            key=lambda m: (
                _LEVEL_ORDER.index(m.level) if m.level in _LEVEL_ORDER else len(_LEVEL_ORDER),
                m.index,
            )
        )
        resolved[code] = MarkerSet(
            criterion=code,
            noun=noun,
            cited=len(indices),
            markers=tuple(markers),
            unresolved=unresolved,
            table_missing=not table,
        )

    return resolved
