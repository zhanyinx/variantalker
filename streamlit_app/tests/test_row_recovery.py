"""Which row the variant dialog is handed, and what happens when it cannot be sure.

Issue #310. ``_full_row`` recovers the pristine row behind the partial one the grid hands
back, and it had two ways of answering confidently and wrongly:

* **It read ``__pandas_index`` as a label.** It is not one. ``st_aggrid`` builds the field as
  ``[str(i) for i in range(len(frame))]`` (``st_aggrid/AgGrid.py:41``, streamlit-aggrid
  1.1.6), so it is the row's **position** in the frame the grid was handed. Resolving it with
  ``.loc`` is right only while the frame's labels happen to be ``0..n-1``, and the frames
  this app draws are boolean-mask slices of the report — the passed and failed splits — whose
  labels are neither contiguous nor zero-based. On such a frame the lookup lands on a
  different variant and the dialog says so with a straight face.

* **The three-column fallback answered when it had nothing to match on.** The mask started
  all-True and was narrowed only by columns present on the partial row, so a partial carrying
  none of ``Hugo_Symbol``/``Chromosome``/``Start_Position`` left it untouched and
  ``matches.iloc[0]`` returned row 0. One of the three is worse than none in a way: it
  narrows to a *set* and then silently takes the first of it. Measured on this machine's
  corpus, 37 of 104 real MAFs repeat at least one ``(gene, chromosome, position)`` triple,
  and one pooled MAF repeats it on 1711 of 2361 rows — so the ambiguous branch is the
  ordinary case, not a corner.

The two are one bug from the reader's side: the dialog opens a variant they did not pick.
They are separated here because only the first is reachable through the app today, and only
the second was what #310 was filed about.

Every frame below is indexed the way a report's passed/failed split leaves it, and every
``__pandas_index`` is the value the component would really send for that row — which is what
the eight tests in ``test_double_click_dialog.py`` could not see, one of them asserting the
label reading outright.
"""

import sys
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

import pandas as pd
import pytest

STREAMLIT_APP = Path(__file__).resolve().parent.parent
if str(STREAMLIT_APP) not in sys.path:
    sys.path.insert(0, str(STREAMLIT_APP))

from components import variant_table  # noqa: E402

STAMP = variant_table._DOUBLE_CLICK_COLUMN
INDEX_KEY = variant_table._AGGRID_INDEX_KEY


# --- what the component actually sends ----------------------------------------------------


def test_the_component_sends_a_position_and_not_a_label():
    """Read off the installed ``st_aggrid``, so the fiction cannot come back quietly.

    This is the fact the whole module turns on, and it is a fact about a third-party package
    that a comment cannot hold. If a future version sends labels instead, this fails and the
    lookup below has to be revisited rather than silently starting to miss.

    Mutation: none needed — it asserts someone else's code. It is here to fail on *upgrade*.
    """

    from st_aggrid.AgGrid import __parse_row_data as parse_row_data

    frame = _masked_frame()
    rows, _ = parse_row_data(frame)

    assert [row[INDEX_KEY] for row in rows] == ["0", "1", "2", "3"], (
        "st_aggrid no longer numbers rows by position; `_grid_position` reads it as one"
    )
    assert list(frame.index) == [2, 3, 5, 7], "the fixture stopped being a masked frame"


# --- fixtures -----------------------------------------------------------------------------


def _masked_frame() -> pd.DataFrame:
    """Four variants with the labels a boolean mask leaves behind, and no repeated triple.

    Distinguishable by gene alone, so a failure names the variant that was opened.
    """

    return pd.DataFrame(
        {
            "Hugo_Symbol": ["APC", "TP53", "KRAS", "BRCA1"],
            "Chromosome": ["chr5", "chr17", "chr12", "chr17"],
            "Start_Position": [112000, 7670000, 25245000, 43050000],
        },
        index=[2, 3, 5, 7],
    )


def _repeated_triple_frame() -> pd.DataFrame:
    """One variant called in three samples — all the fallback's three columns agree.

    The shape 37 of this machine's 104 real MAFs carry. ``Tumor_Sample_Barcode`` is what
    tells the three rows apart, and it is not one of the columns the fallback matches on.
    """

    return pd.DataFrame(
        {
            "Hugo_Symbol": ["TP53", "TP53", "TP53"],
            "Chromosome": ["chr17", "chr17", "chr17"],
            "Start_Position": [7670000, 7670000, 7670000],
            "Tumor_Sample_Barcode": ["first", "second", "third"],
        },
        index=[4, 9, 12],
    )


def _double_click_response(row: dict, stamp: str = "1:1"):
    """A grid response shaped like the one a stamped double-click produces."""

    return SimpleNamespace(
        event_data={
            "streamlitRerunEventTriggerName": "cellValueChanged",
            "type": "cellValueChanged",
            "data": {**row, STAMP: stamp},
            "newValue": stamp,
        },
        selected_rows=None,
        data=[],
    )


def _sent_row(frame: pd.DataFrame, position: int) -> dict:
    """The JSON the component sends when the row at ``position`` is double-clicked."""

    row = frame.iloc[position].to_dict()
    row[INDEX_KEY] = str(position)
    return row


def _opened(monkeypatch, response, full_data, *, button: bool = False):
    """Render the grid with that response and return the rows the dialog was handed."""

    opened = []
    fake_st = MagicMock()
    fake_st.session_state = {}
    fake_st.button.return_value = button
    monkeypatch.setattr(variant_table, "st", fake_st)
    monkeypatch.setattr(variant_table, "GridOptionsBuilder", MagicMock(), raising=False)
    monkeypatch.setattr(variant_table, "AgGrid", lambda *a, **k: response, raising=False)
    monkeypatch.setattr(
        variant_table, "_show_variant_dialog", lambda row: opened.append(row)
    )

    variant_table._render_aggrid_with_detail(
        full_data, full_data, 500, "passed_variants", list(full_data.columns)
    )
    return opened


# --- the double-click opens the row that was double-clicked --------------------------------


@pytest.mark.parametrize("position", [0, 1, 2, 3])
def test_a_double_click_opens_the_row_at_that_position(monkeypatch, position):
    """On a frame whose labels are not its positions — which is every report frame.

    Mutation: resolve the value with ``full_data.loc`` instead of ``iloc`` and positions 2
    and 3 open ``APC`` and ``TP53``, the rows at labels 2 and 3. Positions 0 and 1 pass
    either way, which is why all four are here: half the parametrisation would have called
    the old reading correct.
    """

    frame = _masked_frame()
    expected = frame.iloc[position]["Hugo_Symbol"]

    opened = _opened(
        monkeypatch, _double_click_response(_sent_row(frame, position)), frame
    )

    assert len(opened) == 1, f"the dialog opened {len(opened)} times"
    assert opened[0]["Hugo_Symbol"] == expected, (
        f"double-clicked the row at position {position} ({expected}) and the dialog opened "
        f"{opened[0]['Hugo_Symbol']}"
    )


def test_position_zero_is_a_position_and_not_an_absence(monkeypatch):
    """The first row of the frame, whose ordinal is the one falsy value in the range.

    Guarded on its own because the reading has to be ``is not None``: a truthiness test on
    the position sends exactly this row down the fallback, where it would be answered
    correctly by accident on any frame and wrongly on a frame with a repeated triple.

    Mutation: ``if grid_position:`` in place of ``if grid_position is not None:``. The frame
    here has a repeated triple, so the fallback cannot cover for it.
    """

    frame = _repeated_triple_frame()

    opened = _opened(monkeypatch, _double_click_response(_sent_row(frame, 0)), frame)

    assert len(opened) == 1, f"the dialog opened {len(opened)} times"
    assert opened[0]["Tumor_Sample_Barcode"] == "first"


@pytest.mark.parametrize("position", [0, 1, 2])
def test_a_double_click_picks_the_sample_and_not_just_the_variant(monkeypatch, position):
    """Three rows the fallback's three columns cannot tell apart.

    Mutation: any of them. Before the fix every one of these opened ``first``, because the
    three-column match returned all three rows and ``iloc[0]`` took the first of them.
    """

    frame = _repeated_triple_frame()
    expected = frame.iloc[position]["Tumor_Sample_Barcode"]

    opened = _opened(
        monkeypatch, _double_click_response(_sent_row(frame, position)), frame
    )

    assert len(opened) == 1, f"the dialog opened {len(opened)} times"
    assert opened[0]["Tumor_Sample_Barcode"] == expected, (
        f"double-clicked sample {expected} and the dialog opened "
        f"{opened[0]['Tumor_Sample_Barcode']}"
    )


# --- the button route resolves the same way -----------------------------------------------
#
# It reaches `_full_row` from the other side: `st_aggrid` has already turned
# `__pandas_index` into the index of `selected_rows` (`AgGridReturn.py:266`), so the position
# arrives as the row's *name* rather than as a field. Both routes must come to the same
# answer about which variant the reader picked — that is what `_full_row` exists for.


def _selection_response(frame: pd.DataFrame, position: int):
    selected = pd.DataFrame([_sent_row(frame, position)]).set_index(INDEX_KEY)
    return SimpleNamespace(event_data=None, selected_rows=selected, data=[])


@pytest.mark.parametrize("position", [0, 1, 2, 3])
def test_the_button_opens_the_row_at_that_position(monkeypatch, position):
    """Mutation: as above. The button's ordinal arrives as a string, so before the fix it
    missed the integer labels outright and every press went to the three-column fallback —
    which is why this route was the one that looked right.
    """

    frame = _masked_frame()
    expected = frame.iloc[position]["Hugo_Symbol"]

    opened = _opened(
        monkeypatch, _selection_response(frame, position), frame, button=True
    )

    assert len(opened) == 1, f"the button opened the dialog {len(opened)} times"
    assert opened[0]["Hugo_Symbol"] == expected


@pytest.mark.parametrize("position", [0, 1, 2])
def test_the_button_picks_the_sample_and_not_just_the_variant(monkeypatch, position):
    """The button on a repeated triple, which is where its fallback stopped being enough.

    Mutation: any of the fix. Before it, all three presses opened ``first``.
    """

    frame = _repeated_triple_frame()
    expected = frame.iloc[position]["Tumor_Sample_Barcode"]

    opened = _opened(
        monkeypatch, _selection_response(frame, position), frame, button=True
    )

    assert len(opened) == 1, f"the button opened the dialog {len(opened)} times"
    assert opened[0]["Tumor_Sample_Barcode"] == expected


# --- and when it cannot be sure, it does not guess -----------------------------------------
#
# `_full_row` is called directly from here down. The position is what the app always sends,
# so these states are not reachable through either route as the module stands — which is the
# point: the fallback is now the answer to a question nothing asks, and the cost of it
# refusing is therefore nothing. It is kept, and kept honest, because #147 recorded that the
# component's identifier is not promised and the fallback is what is left when it is absent.


def test_nothing_to_match_on_is_not_answered_with_row_zero():
    """The bug #310 was filed about.

    Mutation: drop the `applied` count and the mask stays all-True, so `matches` is the whole
    frame and `iloc[0]` hands back ``APC`` — the first row of the report, presented as the
    variant the reader picked.
    """

    frame = _masked_frame()
    partial = pd.Series({"Variant_Classification": "Missense_Mutation"})

    recovered = variant_table._full_row(partial, frame)

    assert "Hugo_Symbol" not in recovered.index, (
        f"a row with nothing to match on was answered with {recovered.to_dict()}"
    )
    assert recovered["Variant_Classification"] == "Missense_Mutation", (
        "the partial row the browser sent is what is left to draw"
    )


def test_an_ambiguous_match_is_not_resolved_to_the_first_of_it():
    """One clause applied, three rows satisfying it.

    Mutation: relax `len(matches) == 1` to `len(matches) > 0` and this returns ``first``,
    which is the same silent guess by a shorter route.
    """

    frame = _repeated_triple_frame()
    partial = pd.Series({"Hugo_Symbol": "TP53"})

    recovered = variant_table._full_row(partial, frame)

    assert "Tumor_Sample_Barcode" not in recovered.index, (
        f"an ambiguous match was resolved to one row: {recovered.to_dict()}"
    )


def test_an_unambiguous_match_is_still_answered():
    """The fallback is made strict, not vestigial: one row satisfying every clause is an
    answer, and returning the partial row instead would lose the columns the reader hid.

    Mutation: return `partial` unconditionally and this loses ``Tumor_Sample_Barcode``.
    """

    frame = _masked_frame().assign(Tumor_Sample_Barcode=["a", "b", "c", "d"])
    partial = pd.Series(
        {"Hugo_Symbol": "KRAS", "Chromosome": "chr12", "Start_Position": 25245000}
    )

    recovered = variant_table._full_row(partial, frame)

    assert recovered["Tumor_Sample_Barcode"] == "c"


@pytest.mark.parametrize("raw", ["4", "-1", "", "not-a-number", None, "1.5"])
def test_a_position_the_frame_cannot_hold_is_not_a_position(raw):
    """Out of range, negative, unparseable, absent. Each must fall through rather than raise
    or wrap around — ``iloc[-1]`` is a valid pandas call and the last row of the report.

    Mutation: drop the bounds check and ``"-1"`` opens ``BRCA1``, the last variant in the
    frame, for a reader who clicked something else.
    """

    frame = _masked_frame()

    assert variant_table._grid_position(raw, frame) is None


def test_the_dialog_is_never_handed_the_machinery_on_the_honest_path():
    """The refusal draws the partial row, which arrives carrying both internal names.

    ``test_double_click_dialog.py`` guards this through the render path. It is guarded here
    too because the fallback's return is now reached by a second route — the refusal — and a
    row that is *returned as itself* is exactly the one that could carry them.

    Mutation: drop the `partial.drop(...)` and the dialog lists `_double_click_stamp` and
    `__pandas_index` among the variant's own fields.
    """

    frame = _masked_frame()
    partial = pd.Series(
        {"Variant_Classification": "Silent", STAMP: "1:1", INDEX_KEY: "99"}
    )

    recovered = variant_table._full_row(partial, frame)

    assert STAMP not in recovered.index
    assert INDEX_KEY not in recovered.index
