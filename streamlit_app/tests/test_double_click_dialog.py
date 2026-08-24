"""Opening the variant dialog by double-clicking its row, and keeping it shut (issue #159).

The dev's request was one sentence — *double-clicking a row should open the variant dialog*
— and the whole difficulty is in the second half of it: **staying shut afterwards**. Three
measured facts from issue #159 decide the shape of what is guarded here,
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


def _row(
    row_position: int,
    gene: str = "BRCA1",
    position: int = 41244936,
    *,
    of_rows: int | None = None,
    impossible: bool = False,
) -> dict:
    """The JSON for the row at ``row_position`` in the frame the grid was drawn from.

    The first argument was called ``index`` and read as a label. It is not one: ``st_aggrid``
    numbers the rows it sends by position (issue #310, ``tests/test_row_recovery.py``), so a
    call passing one of ``_frame()``'s labels was describing a payload the component cannot
    produce — and one such call *asserted the bug* while all 22 tests here passed.

    Hence the bounds check, which is issue #311's ruling on what stops that happening again.
    A comment was tried and is not enough: ``variant_detail.py`` has stated this exact truth
    since issue #147 and the bug shipped anyway, because the comment lived at a downstream
    consumer while the code reading the field had none of it. A helper that **refuses** the
    impossible fails at the moment the false premise is typed, rather than passing loudly.

    ``of_rows`` is how many rows the frame in the test has, defaulting to :func:`_frame`'s
    three; pass it when the test builds its own. ``impossible=True`` is the deliberate escape,
    for the one test that needs a position no frame could have sent — and being a keyword it
    cannot be reached by accident.
    """

    limit = len(_frame()) if of_rows is None else of_rows
    if not impossible and not 0 <= row_position < limit:
        raise AssertionError(
            f"position {row_position} is not a row of a {limit}-row frame, so `st_aggrid` "
            f"could not have sent it: it numbers rows `[str(i) for i in range(len(frame))]`. "
            f"If you meant one of the frame's index *labels*, that is the premise issue #310 "
            f"was filed about. Pass `impossible=True` only to exercise the refusal itself."
        )

    return {
        "Hugo_Symbol": gene,
        "Chromosome": "chr17",
        "Start_Position": position,
        INDEX_KEY: str(row_position),
    }


def test_the_helper_refuses_a_payload_the_component_cannot_send():
    """The instrument issue #311 chose over a comment, guarded so it cannot quietly go away.

    ``7`` is the value that was written here five times, because it is a *label* of
    :func:`_frame` — and a three-row frame has no position 7. The escape stays available and
    is asserted too, so a stricter helper cannot break the one test that needs it.

    Mutation: drop the bounds check and this test fails, which is the whole of its job.
    """

    with pytest.raises(AssertionError, match="could not have sent it"):
        _row(7)
    with pytest.raises(AssertionError, match="could not have sent it"):
        _row(-1)
    with pytest.raises(AssertionError, match="could not have sent it"):
        _row(0, of_rows=0)

    assert _row(7, impossible=True)[INDEX_KEY] == "7"
    assert _row(2)[INDEX_KEY] == "2"


# --- the one-shot token -------------------------------------------------------------------
#
# `_fresh_double_click` is where "opens once and stays shut" lives, so it is guarded directly
# as well as through the render path: it is a pure function of the response and a state dict,
# which is the only reason this behaviour is testable without a browser.


def test_a_stamped_double_click_is_a_double_click():
    state = {}
    row = variant_table._fresh_double_click(
        _double_click(_row(1), "1786974405748:1"), state, "k"
    )

    assert row is not None, "a stamped double-click was not recognised at all"
    assert row["Hugo_Symbol"] == "BRCA1"


def test_the_same_payload_arriving_again_is_not_a_second_double_click():
    """The measured hazard: the value persists, so it arrives again on every later rerun.

    Mutation: return the row whenever the event carries one, and this fails on the second
    call — which is what a user sees as the dialog re-opening by itself.
    """

    state = {}
    response = _double_click(_row(1), "1786974405748:1")

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
    row = _row(1)

    assert variant_table._fresh_double_click(_double_click(row, "1:1"), state, "k")
    assert variant_table._fresh_double_click(_double_click(row, "2:1"), state, "k"), (
        "a second double-click on the same row was suppressed"
    )
    assert variant_table._fresh_double_click(_double_click(row, "3:1"), state, "k")


def test_the_two_tabs_do_not_consume_each_other_s_stamps():
    """One state key per grid. Both results tabs render through this same function."""

    state = {}
    response = _double_click(_row(1), "1786974405748:1")

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
        lambda *a, **k: _double_click(_row(1), "1:1"),
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
    """And on the pristine row, found by the identifier the component sends.

    Rewritten by issue #310, which found this test asserting the bug. It sent
    ``__pandas_index: "7"`` for a two-row frame labelled ``[4, 7]`` and required the row at
    *label* 7 — a payload ``st_aggrid`` cannot produce, since it numbers rows by position
    (``[str(i) for i in range(len(frame))]``). Passing it meant resolving that value with
    ``.loc``, which is what opened a different variant on every report frame whose labels are
    not its positions — that is, on every passed/failed split.

    What it guards is unchanged and still worth guarding: the identifier is used, not the
    column match. The frame's two rows agree on all three matched columns, so the column
    match cannot tell them apart and only the identifier can.

    Mutation: read the value with ``full_data.loc`` instead of ``iloc`` and this opens
    ``first`` — label 1 is absent, so it falls to the column match and takes the first of two.

    It is no longer the only cover for that, and is kept as the render-path instance in the
    double-click's own file: ``test_row_recovery.py`` parametrises the same guard over four
    positions of a masked frame and over both routes, deriving each payload from the frame so
    the premise cannot be misstated at all. Read that file first when this one fails.
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

    opened = _opened_row(
        monkeypatch, _double_click(_row(1, of_rows=len(full_data)), "1:1"), full_data
    )

    assert len(opened) == 1, f"the dialog opened {len(opened)} times"
    assert opened[0]["Tumor_Sample_Barcode"] == "second", (
        "the dialog was handed a different variant than the one double-clicked — the "
        "identifier the component sent was not used"
    )


def test_the_dialog_is_never_handed_the_machinery(monkeypatch):
    """No internal name reaches the interface, including where the row cannot be matched.

    An unmatchable row falls back to the partial row the browser sent, which carries both the
    stamp and ``__pandas_index``. The variant dialog renders the fields it is given.

    Mutation: drop the `partial.drop(...)` in `_full_row` and the dialog lists
    ``_double_click_stamp`` and ``__pandas_index`` among the variant's own fields.
    """

    empty = pd.DataFrame(columns=["Hugo_Symbol", "Chromosome", "Start_Position"])
    # The one place `impossible=True` is warranted: an empty frame has no positions at all, so
    # every payload is one the component could not have sent. That is the point — it is how the
    # unmatchable row is reached, and the unmatchable row is what carries the machinery.
    opened = _opened_row(monkeypatch, _double_click(_row(7, impossible=True), "1:1"), empty)

    assert len(opened) == 1, f"the dialog opened {len(opened)} times"
    assert STAMP not in opened[0].index, f"the stamp reached the dialog: {list(opened[0].index)}"
    assert INDEX_KEY not in opened[0].index, (
        f"the component's index key reached the dialog: {list(opened[0].index)}"
    )


def test_the_selection_route_still_opens_the_same_row(monkeypatch):
    """The ``🔍 View details`` button was not disturbed.

    Issue #160 has since settled that the button stays — it is the discoverable route for a
    reader who does not guess that rows are double-clickable — so both routes reach the dialog,
    and both go through the same row recovery so they cannot come to different answers about
    which variant the user picked. They *did* come to different answers, on screen, in one
    Chromium run (issue #309): this route missed the identifier outright because
    ``selected_rows`` leaves the ordinal a string, and resolved by the column match instead.
    """

    full_data = _frame()
    # Position 1 of `_frame()`, which is BRCA1. It was `_row(7)` — the label — and passed on
    # the column match while the position path was never reached; #310 made the two routes
    # share one reader for this value, so the payload here has to be one the component sends.
    selected = pd.DataFrame([{**_row(1)}]).set_index(INDEX_KEY)

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
