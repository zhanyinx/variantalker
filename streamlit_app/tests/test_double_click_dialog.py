"""Opening the variant dialog by double-clicking its row, and keeping it shut (issue #159).

The dev's request was one sentence — *double-clicking a row should open the variant dialog*
— and the whole difficulty is in the second half of it: **staying shut afterwards**. Three
measured facts from ``docs/wayfinder/issue-159/`` decide the shape of what is guarded here,
and each of them was reproduced in a real browser because ``AppTest`` never runs
``st_aggrid``'s JavaScript at all:

* a component's value **persists across reruns**, so opening the dialog on a standing "the
  last event was a double-click" condition re-opens it on every later interaction — measured
  at 8 openings for 5 double-clicks;
* the payload of a *genuine* second double-click on the same row is **byte-identical** to
  that persisted value, all 249KB of it, so no comparison of the event's content can tell
  them apart — a guard built that way suppresses the re-open a user asks for (measured:
  ``[True, False, False]`` over three consecutive double-clicks on one row);
* dismissing an ``st.dialog`` reruns nothing, so there is no close event to clear a token on.

Hence the stamp: the double-click writes a value carrying ``Date.now()`` into a hidden cell,
AG Grid dispatches a fresh ``cellValueChanged`` carrying it, and the stamp is consumed the
moment the row is handed to the dialog. Every guard below is mutation-checked — a test that
passes on both the standing condition and the one-shot token is not a guard.
"""

from __future__ import annotations

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


def _frame() -> pd.DataFrame:
    """Three variants, indexed the way a report's passed/failed split leaves them."""

    return pd.DataFrame(
        {
            "Hugo_Symbol": ["TP53", "BRCA1", "EGFR"],
            "Chromosome": ["chr17", "chr17", "chr7"],
            "Start_Position": [7577120, 41244936, 55249071],
        },
        index=[4, 7, 11],
    )


def _double_click(row: dict, stamp: str):
    """A grid response shaped like the one a stamped double-click actually produces.

    The field names are the ones measured off the running component: ``eventData`` carries
    ``data`` (the row, plus ``__pandas_index`` and the stamped cell) and the trigger's name.
    """

    return SimpleNamespace(
        event_data={
            "streamlitRerunEventTriggerName": "cellValueChanged",
            "type": "cellValueChanged",
            "data": {**row, STAMP: stamp},
            "newValue": stamp,
            "rowIndex": 1,
        },
        selected_rows=None,
        data=[],
    )


def _row(index: int, gene: str = "BRCA1", position: int = 41244936) -> dict:
    return {
        "Hugo_Symbol": gene,
        "Chromosome": "chr17",
        "Start_Position": position,
        INDEX_KEY: str(index),
    }


# --- the one-shot token -------------------------------------------------------------------
#
# `_fresh_double_click` is where "opens once and stays shut" lives, so it is guarded directly
# as well as through the render path: it is a pure function of the response and a state dict,
# which is the only reason this behaviour is testable without a browser.


def test_a_stamped_double_click_is_a_double_click():
    state = {}
    row = variant_table._fresh_double_click(
        _double_click(_row(7), "1786974405748:1"), state, "k"
    )

    assert row is not None, "a stamped double-click was not recognised at all"
    assert row["Hugo_Symbol"] == "BRCA1"


def test_the_same_payload_arriving_again_is_not_a_second_double_click():
    """The measured hazard: the value persists, so it arrives again on every later rerun.

    Mutation: return the row whenever the event carries one, and this fails on the second
    call — which is what a user sees as the dialog re-opening by itself.
    """

    state = {}
    response = _double_click(_row(7), "1786974405748:1")

    assert variant_table._fresh_double_click(response, state, "k") is not None
    for _ in range(3):
        assert variant_table._fresh_double_click(response, state, "k") is None, (
            "the same payload was read as a fresh double-click, so the dialog re-opens "
            "on every interaction after it"
        )


def test_the_same_row_double_clicked_again_does_open_again():
    """The other half, and the one a content comparison fails.

    A user who closes the dialog and double-clicks the same variant again must get it back.
    The row, its index and its position are all identical here — only the stamp differs, so
    a guard keyed on which row was clicked cannot pass this.

    Mutation: key the token on `__pandas_index` or on the row's identity and this fails.
    """

    state = {}
    row = _row(7)

    assert variant_table._fresh_double_click(_double_click(row, "1:1"), state, "k")
    assert variant_table._fresh_double_click(_double_click(row, "2:1"), state, "k"), (
        "a second double-click on the same row was suppressed"
    )
    assert variant_table._fresh_double_click(_double_click(row, "3:1"), state, "k")


def test_the_two_tabs_do_not_consume_each_other_s_stamps():
    """One state key per grid. Both results tabs render through this same function."""

    state = {}
    response = _double_click(_row(7), "1786974405748:1")

    assert variant_table._fresh_double_click(response, state, "passed") is not None
    assert variant_table._fresh_double_click(response, state, "failed") is not None
    assert variant_table._fresh_double_click(response, state, "passed") is None


def test_each_grid_remembers_its_own_stamp(monkeypatch):
    """The guard above tests the helper; this one tests the key the render path passes it.

    Both results tabs render through this same function, so one shared key would let a
    double-click in one tab consume the other tab's stamp.

    Mutation: pass a constant key instead of one per ``key_suffix`` and this fails.
    """

    opened = []
    fake_st = MagicMock()
    fake_st.session_state = {}
    fake_st.button.return_value = False
    monkeypatch.setattr(variant_table, "st", fake_st)
    monkeypatch.setattr(variant_table, "GridOptionsBuilder", MagicMock(), raising=False)
    monkeypatch.setattr(
        variant_table,
        "AgGrid",
        lambda *a, **k: _double_click(_row(7), "1:1"),
        raising=False,
    )
    monkeypatch.setattr(
        variant_table, "_show_variant_dialog", lambda row: opened.append(row)
    )

    frame = _frame()
    variant_table._render_aggrid_with_detail(
        frame, frame, 500, "failed_variants", list(frame.columns)
    )

    assert opened, "the double-click did not open the dialog through the render path"
    remembered = [key for key in fake_st.session_state if "failed_variants" in str(key)]
    assert remembered, (
        f"no per-tab state remembers the stamp: {list(fake_st.session_state)!r}"
    )


@pytest.mark.parametrize(
    "response",
    [
        MagicMock(),
        SimpleNamespace(event_data=None),
        SimpleNamespace(event_data={"streamlitRerunEventTriggerName": "selectionChanged"}),
        SimpleNamespace(event_data={"data": {"Hugo_Symbol": "TP53"}}),
        SimpleNamespace(event_data={"data": {"Hugo_Symbol": "TP53", STAMP: ""}}),
    ],
    ids=["mock", "none", "selection-only", "unstamped-row", "blank-stamp"],
)
def test_nothing_else_counts_as_a_double_click(response):
    """Including a ``MagicMock``.

    Not a hypothetical: several tests in this suite stand ``AgGrid`` up as a ``MagicMock``,
    whose ``event_data`` is a truthy mock whose ``.get`` answers with another truthy mock.
    Read loosely, every one of those tests would open the variant dialog.

    Mutation: swap either `isinstance` for a truthiness test and the ``mock`` case fails.
    """

    assert variant_table._fresh_double_click(response, {}, "k") is None


# --- what the grid is configured to do ----------------------------------------------------


def _captured_grid_call(monkeypatch, frame):
    """Render the grid with a real ``GridOptionsBuilder`` and capture the ``AgGrid`` call."""

    if not variant_table.AGGRID_AVAILABLE:  # pragma: no cover - st_aggrid is a requirement
        pytest.skip("streamlit-aggrid is not installed")

    calls = {}

    def _capture(*args, **kwargs):
        calls["args"] = args
        calls["kwargs"] = kwargs
        return SimpleNamespace(event_data=None, selected_rows=None, data=[])

    fake_st = MagicMock()
    fake_st.session_state = {}
    fake_st.button.return_value = False
    monkeypatch.setattr(variant_table, "st", fake_st)
    monkeypatch.setattr(variant_table, "AgGrid", _capture, raising=False)

    variant_table._render_aggrid_with_detail(
        frame, frame, 500, "passed_variants", list(frame.columns)
    )
    assert calls, "AgGrid was never called"
    return calls


def test_the_grid_asks_to_hear_about_the_cell_the_double_click_writes(monkeypatch):
    """``update_on``, and ``update_mode`` still beside it.

    ``st_aggrid`` reports any AG Grid event named in ``update_on`` and nothing else, so
    without this the double-click never reaches Python. ``SELECTION_CHANGED`` stays because
    the ``🔍 View details`` button reads the selection.

    Mutation: drop either argument and one half of the table stops working.
    """

    calls = _captured_grid_call(monkeypatch, _frame())

    assert "cellValueChanged" in calls["kwargs"]["update_on"], (
        f"the grid does not subscribe to the stamped cell: {calls['kwargs'].get('update_on')!r}"
    )
    assert calls["kwargs"]["update_mode"] == "SELECTION_CHANGED", (
        "the grid stopped reporting selection, which the View details button reads"
    )
    assert calls["kwargs"]["allow_unsafe_jscode"] is True, (
        "the double-click handler is JsCode and does not run without this"
    )


def test_the_double_click_handler_writes_the_stamp(monkeypatch):
    """The handler itself, read off the built grid options.

    It has to write a *cell*: a JsCode handler cannot stamp the outgoing event, because
    ``st_aggrid`` has already cloned and sent it by the time the handler runs — measured, with
    the handler shouting on the browser console to prove it fired.

    Mutation: stamp the event (`event.toIndex = Date.now()`) instead and no double-click ever
    reaches Python at all.
    """

    calls = _captured_grid_call(monkeypatch, _frame())
    options = calls["kwargs"]["gridOptions"]

    handler = options.get("onRowDoubleClicked")
    assert handler is not None, "no double-click handler is installed on the grid"
    source = getattr(handler, "js_code", handler)
    assert "setDataValue" in str(source), (
        f"the handler does not write a cell, so nothing new reaches Python: {source!r}"
    )
    assert STAMP in str(source), f"the handler does not write the stamp column: {source!r}"
    assert "Date.now()" in str(source), (
        f"the handler writes no per-click value, so a second double-click on the same row "
        f"cannot be told from the first: {source!r}"
    )


def test_the_stamp_column_is_hidden_and_cannot_be_unhidden(monkeypatch):
    """It is machinery, not data. The side bar must not offer to show it.

    Mutation: drop ``suppressColumnsToolPanel`` and the columns panel lists
    ``_double_click_stamp`` beside the user's own columns.
    """

    calls = _captured_grid_call(monkeypatch, _frame())
    options = calls["kwargs"]["gridOptions"]

    stamp_def = [c for c in options["columnDefs"] if c.get("field") == STAMP]
    assert stamp_def, f"the stamp column has no column definition: {options['columnDefs']!r}"
    assert stamp_def[0].get("hide") is True, f"the stamp column is visible: {stamp_def[0]!r}"
    assert stamp_def[0].get("suppressColumnsToolPanel") is True, (
        f"the stamp column can be unhidden from the side bar: {stamp_def[0]!r}"
    )
    assert stamp_def[0].get("editable") is False, (
        f"the stamp column is editable: {stamp_def[0]!r}"
    )


def test_the_stamp_reaches_the_grid_and_not_the_frame(monkeypatch):
    """The download reads the column list this module returns, so the stamp must stay out.

    Mutation: add the column to ``display_data`` itself rather than to a copy, and the file
    offered beneath the grid grows a column named ``_double_click_stamp``.
    """

    frame = _frame()
    calls = _captured_grid_call(monkeypatch, frame)

    handed_to_grid = calls["args"][0]
    assert STAMP in handed_to_grid.columns, "the grid was never given the stamp column"
    assert STAMP not in frame.columns, (
        "the stamp was added to the caller's frame, which is what the download reads"
    )


# --- which row the dialog is handed -------------------------------------------------------


def _opened_row(monkeypatch, response, full_data):
    """Render with the given grid response and return the row the dialog was handed."""

    opened = []
    fake_st = MagicMock()
    fake_st.session_state = {}
    fake_st.button.return_value = False
    monkeypatch.setattr(variant_table, "st", fake_st)
    monkeypatch.setattr(variant_table, "GridOptionsBuilder", MagicMock(), raising=False)
    monkeypatch.setattr(variant_table, "AgGrid", lambda *a, **k: response, raising=False)
    monkeypatch.setattr(
        variant_table, "_show_variant_dialog", lambda row: opened.append(row)
    )

    frame = full_data
    variant_table._render_aggrid_with_detail(
        frame, frame, 500, "passed_variants", list(frame.columns)
    )
    return opened


def test_a_double_click_opens_the_dialog_on_the_row_that_was_clicked(monkeypatch):
    """And on the pristine row, found by the index the component sends.

    ``__pandas_index`` arrives as a string while a MAF's index is integer, so the raw value
    alone misses and falls through to matching on three columns.

    Mutation: stop converting it and this still finds BRCA1 by column match — so the guard
    below checks the *index* was used, by giving the frame a second row the column match
    would also accept.
    """

    full_data = pd.DataFrame(
        {
            "Hugo_Symbol": ["BRCA1", "BRCA1"],
            "Chromosome": ["chr17", "chr17"],
            "Start_Position": [41244936, 41244936],
            "Tumor_Sample_Barcode": ["first", "second"],
        },
        index=[4, 7],
    )

    opened = _opened_row(monkeypatch, _double_click(_row(7), "1:1"), full_data)

    assert len(opened) == 1, f"the dialog opened {len(opened)} times"
    assert opened[0]["Tumor_Sample_Barcode"] == "second", (
        "the dialog was handed a different variant than the one double-clicked — the index "
        "the component sent was not used"
    )


def test_the_dialog_is_never_handed_the_machinery(monkeypatch):
    """No internal name reaches the interface, including where the row cannot be matched.

    An unmatchable row falls back to the partial row the browser sent, which carries both the
    stamp and ``__pandas_index``. The variant dialog renders the fields it is given.

    Mutation: drop the `partial.drop(...)` in `_full_row` and the dialog lists
    ``_double_click_stamp`` and ``__pandas_index`` among the variant's own fields.
    """

    empty = pd.DataFrame(columns=["Hugo_Symbol", "Chromosome", "Start_Position"])
    opened = _opened_row(monkeypatch, _double_click(_row(7), "1:1"), empty)

    assert len(opened) == 1, f"the dialog opened {len(opened)} times"
    assert STAMP not in opened[0].index, f"the stamp reached the dialog: {list(opened[0].index)}"
    assert INDEX_KEY not in opened[0].index, (
        f"the component's index key reached the dialog: {list(opened[0].index)}"
    )


def test_the_selection_route_still_opens_the_same_row(monkeypatch):
    """The ``🔍 View details`` button was not disturbed.

    What happens to that button is issue #160's to settle; until it does, both routes must
    reach the dialog, and both must go through the same row recovery so they cannot come to
    different answers about which variant the user picked.
    """

    full_data = _frame()
    selected = pd.DataFrame([{**_row(7)}]).set_index(INDEX_KEY)

    opened = []
    fake_st = MagicMock()
    fake_st.session_state = {}
    fake_st.button.return_value = True
    monkeypatch.setattr(variant_table, "st", fake_st)
    monkeypatch.setattr(variant_table, "GridOptionsBuilder", MagicMock(), raising=False)
    monkeypatch.setattr(
        variant_table,
        "AgGrid",
        lambda *a, **k: SimpleNamespace(
            event_data=None, selected_rows=selected, data=[]
        ),
        raising=False,
    )
    monkeypatch.setattr(
        variant_table, "_show_variant_dialog", lambda row: opened.append(row)
    )

    variant_table._render_aggrid_with_detail(
        full_data, full_data, 500, "passed_variants", list(full_data.columns)
    )

    assert len(opened) == 1, f"the button opened the dialog {len(opened)} times"
    assert opened[0]["Hugo_Symbol"] == "BRCA1"
