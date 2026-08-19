"""The grid-less path's route into the variant dialog (issue #160).

``create_data_table`` has two exits. The one everybody sees draws ``st_aggrid``; the other is
reached whenever ``st_aggrid`` is absent, draws ``st.dataframe``, and offers a ``st.selectbox``
of variants. It cannot offer the double-click issue #159 built — ``st.dataframe`` dispatches no
such event — and that difference is deliberate. What was **not** deliberate is what this file
guards:

    the selectbox was given no ``index``, so it answered ``0`` rather than ``None``, and the
    guard beneath it — ``if row_idx is not None`` — was therefore true on *every render*.

Measured on the code before the fix, three untouched renders of a three-variant frame opened
the dialog three times, on a variant nobody had chosen. This path is not inside a fragment
either, so any interaction anywhere on the page reran it and the dialog came back. It is the
standing-condition failure ``tests/test_double_click_dialog.py`` exists to prevent on the grid,
which had been sitting in the fallback the whole time.

The fix is not ``index=None`` on its own — that only postpones it, since a chosen variant keeps
being answered on every later rerun — but moving the opening behind a **press**, which happens
once. So the two paths differ on the gesture and agree on the rule: no dialog appears unbidden.

The first guard below runs against **real** Streamlit in bare mode rather than a ``MagicMock``,
and that is the point of it. The defect lived entirely in what ``st.selectbox`` answers when it
is handed no ``index``; a fake whose return value the test chooses cannot see it, and would pass
just as happily on the broken code.
"""

from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import MagicMock

import pandas as pd

STREAMLIT_APP = Path(__file__).resolve().parent.parent
if str(STREAMLIT_APP) not in sys.path:
    sys.path.insert(0, str(STREAMLIT_APP))

from components import variant_table  # noqa: E402


def _frame() -> pd.DataFrame:
    """Three variants, carrying enough for the dialog's key to be whole."""

    return pd.DataFrame(
        {
            "Hugo_Symbol": ["TP53", "BRCA1", "EGFR"],
            "Chromosome": ["chr17", "chr17", "chr7"],
            "Start_Position": [7577120, 41244936, 55249071],
            "Reference_Allele": ["C", "A", "G"],
            "Tumor_Seq_Allele2": ["T", "G", "A"],
            "Variant_Classification": ["Missense_Mutation"] * 3,
        }
    )


def _opened(monkeypatch) -> list:
    """Record every variant handed to the dialog, and draw none of it."""

    opened = []
    monkeypatch.setattr(
        variant_table,
        "_show_variant_dialog",
        lambda row: opened.append(row.get("Hugo_Symbol", "?")),
    )
    monkeypatch.setattr(variant_table, "AGGRID_AVAILABLE", False)
    return opened


def _fallback_st() -> MagicMock:
    """A stand-in ``st`` for the questions that are about wiring rather than widget defaults.

    ``button`` answers ``False`` unless a test says otherwise: a bare ``MagicMock`` returns a
    truthy mock from every call, so left alone it would report the button as pressed on every
    render and this file's central claim would pass for the wrong reason.
    """

    fake_st = MagicMock()
    fake_st.checkbox.return_value = False
    fake_st.multiselect.return_value = []
    fake_st.button.return_value = False
    return fake_st


class _RecordingStreamlit:
    """Real Streamlit, with the buttons it is asked to draw written down.

    Not a ``MagicMock``: the whole defect lived in what ``st.selectbox`` answers when handed
    no ``index``, so the selectbox has to be the real one and compute its own answer. Only
    ``button`` is intercepted, and it answers ``False`` — the honest reading of a render in
    which nobody pressed anything.
    """

    def __init__(self, real):
        self._real = real
        self.buttons = []

    def button(self, label, *args, **kwargs):
        self.buttons.append(label)
        return False

    def __getattr__(self, name):
        return getattr(self._real, name)


def test_the_grid_less_path_offers_nothing_on_its_own(monkeypatch):
    """Render it three times, touch nothing: no dialog, and no button either.

    Both halves are needed, and the second is the one with teeth. Once the opening sits behind
    a press, a selectbox that quietly answers ``0`` no longer opens anything unbidden — it
    draws a *View details* button for the first variant in the file, which nobody chose, and
    a guard that only counted dialogs would sail past it. Mutations: drop ``index=None`` and
    this records a button naming ``TP53``; drop the ``st.button`` gate as well and it records
    the three openings the code produced before issue #160.
    """

    assert getattr(variant_table.st, "__name__", None) == "streamlit", (
        "another test left a stand-in behind on `variant_table.st`; this guard measures what "
        "the real widget answers and is worthless against a mock"
    )

    opened = _opened(monkeypatch)
    recording = _RecordingStreamlit(variant_table.st)
    monkeypatch.setattr(variant_table, "st", recording)
    data = _frame()

    for _ in range(3):
        variant_table.create_data_table(data, key_suffix="passed_variants")

    assert opened == [], (
        f"the grid-less path opened the variant dialog {len(opened)} times with nobody "
        f"choosing anything: {opened!r}"
    )
    assert recording.buttons == [], (
        "the grid-less path offered a variant nobody had chosen: "
        f"{recording.buttons!r}"
    )


def test_choosing_a_variant_draws_the_button_but_does_not_open_it(monkeypatch):
    """A choice is not a press.

    Mutation: call `_show_variant_dialog` on the selection rather than inside the `st.button`
    branch, and this fails — which is the shape of the original defect.
    """

    opened = _opened(monkeypatch)
    fake_st = _fallback_st()
    fake_st.selectbox.return_value = 1
    monkeypatch.setattr(variant_table, "st", fake_st)

    variant_table.create_data_table(_frame(), key_suffix="passed_variants")

    assert fake_st.button.called, "no button was drawn for the chosen variant"
    assert opened == [], f"the dialog opened without the button being pressed: {opened!r}"


def test_pressing_the_button_opens_the_variant_that_was_chosen(monkeypatch):
    """And opens it once, with the row the selector was pointing at — not the first row."""

    opened = _opened(monkeypatch)
    fake_st = _fallback_st()
    fake_st.selectbox.return_value = 1
    fake_st.button.return_value = True
    monkeypatch.setattr(variant_table, "st", fake_st)

    variant_table.create_data_table(_frame(), key_suffix="passed_variants")

    assert opened == ["BRCA1"], (
        f"pressing the button did not open the chosen variant exactly once: {opened!r}"
    )


def test_no_variant_chosen_means_no_button_at_all(monkeypatch):
    """The empty state matches the grid's, where no selected row means no button.

    Left on ``index=0`` the button would stand there permanently offering the first variant in
    the file, which nobody chose.
    """

    opened = _opened(monkeypatch)
    fake_st = _fallback_st()
    fake_st.selectbox.return_value = None
    monkeypatch.setattr(variant_table, "st", fake_st)

    variant_table.create_data_table(_frame(), key_suffix="passed_variants")

    assert not fake_st.button.called, (
        "a View details button was drawn with no variant chosen: "
        f"{[call.args for call in fake_st.button.call_args_list]!r}"
    )
    assert opened == [], f"the dialog opened with nothing chosen: {opened!r}"


def test_both_paths_label_the_button_the_same_way(monkeypatch):
    """Two buttons, one spelling — read off the paths that draw them, not off the helper.

    Issue #149's ``chrchr17`` was a `chr` literal written into two places that had drifted
    apart, and issue #160 added a second button. Mutation: re-inline the label at either call
    site and this fails as soon as the two spellings differ by so much as a space.
    """

    variant = {
        "Hugo_Symbol": "TP53",
        "Chromosome": "chr17",
        "Start_Position": 7577120,
    }

    # The grid path: a selected row, and a button drawn under it.
    grid_st = _fallback_st()
    monkeypatch.setattr(variant_table, "st", grid_st)
    monkeypatch.setattr(variant_table, "GridOptionsBuilder", MagicMock(), raising=False)
    grid_response = MagicMock()
    grid_response.selected_rows = [variant]
    monkeypatch.setattr(
        variant_table, "AgGrid", lambda *a, **k: grid_response, raising=False
    )
    frame = _frame()
    variant_table._render_aggrid_with_detail(
        frame, frame, 500, "passed_variants", list(frame.columns)
    )
    grid_labels = [call.args[0] for call in grid_st.button.call_args_list if call.args]

    # The grid-less path: the same variant, chosen in the selectbox.
    _opened(monkeypatch)
    fallback_st = _fallback_st()
    fallback_st.selectbox.return_value = 0
    monkeypatch.setattr(variant_table, "st", fallback_st)
    variant_table.create_data_table(frame, key_suffix="passed_variants")
    fallback_labels = [
        call.args[0] for call in fallback_st.button.call_args_list if call.args
    ]

    assert grid_labels, "the grid path drew no button to compare"
    assert fallback_labels, "the grid-less path drew no button to compare"
    assert grid_labels == fallback_labels, (
        f"the two paths label the same variant differently: {grid_labels!r} vs "
        f"{fallback_labels!r}"
    )
    assert "chrchr" not in grid_labels[0], (
        f"the button label doubled the chromosome prefix: {grid_labels[0]!r}"
    )
