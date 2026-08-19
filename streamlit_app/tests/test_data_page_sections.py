"""The section switch: fixed order, and a jump that fires once (issue #59).

Called *Load Data* / *Results* when this module was written. The first section is now the
**filter run** (issue #65): once the uploader moved to the sidebar (issue #64) and the load's
banners moved up to the page (issue #59), what was left on it was how the report was reached
and the button that re-reaches it — no loading at all. The switch itself is drawn only when
a file is open, since with nothing open both sections say the same thing.

The two sections of the data page used to be ``st.tabs``, rendered in one order before a file was
filtered and in the *other* order after — because Streamlit always opens the first tab, and Results
was easy to miss as the second one. Reordering surfaced Results at the cost of moving the page's
primary control out from under the user, which is what the walkthrough in #53 reported.

The order is now fixed and Results is reached by *navigation*: the page selects it once,
when a filter run produces a report. That is only expressible with a widget whose value
can be written from code, and ``st.tabs`` is not one — it takes neither a key nor a
default — so the switch is an ``st.segmented_control``.

Three claims, and the third is the one with teeth:

1. the page opens on the filter run, with both sections offered in a fixed order;
2. a filter run that produces a report selects Results;
3. **a rerender that merely re-meets the same file does not.** The uploader hands the
   same file back on every rerender, so a jump wired to "a file is present" rather than to
   "a file has arrived" re-fires forever, and the user who clicks back to the first section
   is dragged to Results again — a worse trap than the swapping tabs it replaces. The page
   guards this with a token, and this module exists mostly to keep that guard honest.

Why every test renders exactly once
-----------------------------------
``AppTest`` cannot drive this widget across runs, and the failure is silent rather than
loud, so it is written down here.

``AppTest`` models a segmented control with its ``ButtonGroup`` class, which was written
for ``st.feedback`` and assumes the widget's value is a *list*. In single-selection mode
``st.segmented_control`` returns a scalar, so ``ButtonGroup.indices`` — ``[options.index(
format_func(v)) for v in self.value]`` — iterates the string ``"filter_run"`` one character
at a time and hands the next run a selection that has nothing to do with what was selected.
A second ``app.run()`` therefore arrives on the wrong section, having reported no error.

Nothing about that reaches the running app, where the browser returns the selected
*index*; it is a limit on this harness only. The response is one render per test, seeded
into the state under test rather than clicked into it — which is why the trap test seeds a
stamped token instead of selecting the first section and rerunning.

A unit module: nothing here has a pipeline counterpart — the pipeline has no tabs — and
none of it needs ``bin/``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

STREAMLIT_APP = Path(__file__).resolve().parent.parent
if str(STREAMLIT_APP) not in sys.path:
    sys.path.insert(0, str(STREAMLIT_APP))

FIXTURES = STREAMLIT_APP / "tests" / "fixtures" / "parity"
MAF = FIXTURES / "somatic_reference.maf"

#: A MAF whose depth and VAF columns hold ``.``. The filter refuses it rather than
#: guessing, so it is the load that produces no report — see ``test_numeric_columns.py``.
UNREADABLE_MAF = FIXTURES / "somatic_dot_numeric.maf"

#: The page under test, handed a chosen file the way the sidebar's chooser hands one over
#: (issue #64): a name plus bytes from ``getvalue()``, parked in ``PENDING_UPLOAD_KEY``, which is
#: all the load path reads. There is no uploader on this page to stub any more. Whether a
#: file has been chosen is read from session state, so a test can render the page with or
#: without one without touching the page itself.
_SCRIPT = """
import os
import sys
sys.path.insert(0, {app!r})

import streamlit as st
from components.sidebar import PENDING_UPLOAD_KEY
from config.pipeline_params import pipeline_params
from page_modules import data_loading


class UploadStub:
    def __init__(self, path):
        self.name = os.path.basename(path)
        self._path = path

    def getvalue(self):
        with open(self._path, "rb") as handle:
            return handle.read()


def _chosen_file():
    try:
        chosen = st.session_state["upload"]
    except KeyError:
        chosen = False
    if not chosen:
        return None
    return UploadStub(chosen if isinstance(chosen, str) else {maf!r})


_pending = _chosen_file()
if _pending is not None:
    st.session_state[PENDING_UPLOAD_KEY] = _pending

if "filter_params" not in st.session_state or st.session_state.filter_params is None:
    st.session_state.filter_params = pipeline_params("somatic")

data_loading.show_data_loading_page()
"""


def _render(**seed):
    """The data page, rendered once, over a seeded session."""
    pytest.importorskip("streamlit", reason="streamlit not installed")
    from streamlit.testing.v1 import AppTest

    app = AppTest.from_string(_SCRIPT.format(app=str(STREAMLIT_APP), maf=str(MAF)))
    for key, value in seed.items():
        app.session_state[key] = value
    app.run(timeout=180)
    assert not app.exception, [str(e.value) for e in app.exception]
    return app


def _state(app, key, default=None):
    """``AppTest.session_state`` has no ``get``."""
    try:
        return app.session_state[key]
    except (KeyError, AttributeError):
        return default


def _section(app):
    return _state(app, "data_page_section")


def _labels(app):
    """What the switch reads on screen. A leading emoji becomes the segment's icon."""
    groups = app.get("button_group")
    assert groups, "the data page rendered no section control"
    return [f"{option.content_icon} {option.content}".strip() for option in groups[0].options]


def test_no_switch_is_offered_before_a_file_is_open():
    """The page draws the way in, not a choice between two empty sections (issue #65).

    Both sections say *load a file first* in this state, so the switch was the page arguing
    with itself — and once the first section stopped being the load (the uploader went to
    the sidebar in issue #64), a segment labelled *Filter run* offered before any file
    exists names something the user cannot do.
    """
    app = _render()

    assert not app.get("button_group"), "a section switch was drawn with no file open"
    assert _state(app, "maf_data") is None
    assert any("Open a MAF file" in header.value for header in app.subheader)


def test_both_sections_are_offered_filter_run_first():
    """The order is fixed, across the two states the old tabs code reordered between.

    Those states are *a file with no report* and *a file with a report* — which is what the
    ``st.tabs`` version keyed its ordering off, and why Results used to move out from under
    the user. Before issue #65 this was checked as before-upload against after-upload; the
    no-file render now draws no switch at all, so the comparison moves to the two states
    that actually differ while a switch is on screen.
    """
    without_report = _labels(_render(upload=str(UNREADABLE_MAF)))
    with_report = _labels(_render(upload=True))

    assert without_report[0] == with_report[0] == "🔍 Filter run"
    assert without_report[1].startswith("📊 Results")
    assert with_report[1].startswith("📊 Results")


def test_the_filter_run_section_holds_only_the_run():
    """The six blocks are one (issue #65).

    Data Overview's three ``st.metric`` tiles, and the two ``st.dataframe`` expanders —
    Data Preview and Column Information — were the bulk of this section and are gone. So is
    the Configure Parameters button, which was a third door to a page the sidebar's nav
    radio names on this same screen. Asserted by *kind* rather than by label: what this
    pins is that the section stopped being a place blocks accumulate.
    """
    app = _render(upload=str(UNREADABLE_MAF))

    assert _section(app) == "filter_run"
    assert not app.get("metric"), "the filter run section drew metric tiles"
    assert not app.get("arrow_data_frame"), "the filter run section drew a table"

    labels = [button.label for button in app.button]
    assert labels == ["🔍 Re-apply Filters"], labels


def test_a_filter_run_selects_results():
    app = _render(upload=True)

    assert _section(app) == "results"
    assert _state(app, "filtered_data") is not None


def test_the_results_label_carries_the_report_size():
    app = _render(upload=True)

    passed = len(app.session_state["filtered_data"])
    assert _labels(app)[1] == f"📊 Results ({passed})"


def test_the_jump_is_spent_once_it_has_been_honoured():
    """A jump left set would be honoured again on the next rerender."""
    app = _render(upload=True)

    assert _state(app, "jump_to_results", "<absent>") in ("<absent>", False)


def test_meeting_the_same_file_again_neither_reloads_it_nor_moves_the_user():
    """The trap this ticket exists to avoid, and the re-read that came with it.

    A stamped token and a user on Load Data is the state after they have jumped to Results
    once and navigated back. The chooser still holds the file, so the load path is re-entered
    on this render — and must stay quiet.

    Since issue #64 the token guards the *load* as well as the jump: the page body runs on
    every render of either section, so an unguarded read would re-read and re-filter the
    whole MAF on every click anywhere in the app. ``maf_data`` staying unset with a file
    still on offer is what says the read did not happen.

    The section is no longer asserted here. With no file open the page draws no switch at
    all (issue #65), so a section assertion in this seeded state would be reading back the
    seed rather than the page. *The user was not moved* is carried by the jump instead, and
    carried better: nothing can move them if no jump was ever requested.
    """
    app = _render(
        upload=True,
        data_page_section="filter_run",
        last_upload_token=(MAF.name, None),
    )

    assert _state(app, "maf_data") is None
    assert _state(app, "jump_to_results", "<absent>") in ("<absent>", False)


def test_a_different_file_gets_its_own_jump():
    """The guard is per file, not once per session."""
    app = _render(
        upload=True,
        data_page_section="filter_run",
        last_upload_token=("something_else.maf", None),
    )

    assert _section(app) == "results"
    assert _state(app, "last_upload_token") == (MAF.name, None)


def test_a_refused_maf_does_not_open_an_empty_report():
    """No report, no jump.

    Jumping on a refusal would land the user on a Results section with nothing in it, and
    take them away from the banner saying why. The banner is asserted against the running
    app by ``make app-load-check``; what is asserted here is that the page stayed put.
    """
    app = _render(upload=str(UNREADABLE_MAF))

    assert _state(app, "filtered_data") is None
    assert _section(app) == "filter_run"
    assert _state(app, "jump_to_results", "<absent>") in ("<absent>", False)


def test_a_refused_maf_still_says_what_is_in_the_file():
    """The counts survive the state that has no report (issue #65).

    Samples and Genes moved off the load page onto the results summary — but that summary
    lives inside the report's own tab, and the results section returns before it whenever
    there is no report. A refused MAF is exactly that state, and it is the one where the old
    Data Overview still had something to show, so rehousing them into the summary alone
    would have dropped them precisely where they were still wanted.

    The MAF here is refused by the *filter*, not by validation, so the file is loaded and
    ``maf_data`` holds it — which is what lets the counts be asserted against the frame they
    describe rather than against numbers copied into this file.
    """
    app = _render(upload=str(UNREADABLE_MAF), data_page_section="results")

    assert _state(app, "filtered_data") is None

    maf = _state(app, "maf_data")
    assert maf is not None, "the refusal was a validation failure, not a filter refusal"

    tiles = {tile.label: tile.value for tile in app.get("metric")}
    assert tiles == {
        "Samples": str(maf["Tumor_Sample_Barcode"].nunique()),
        "Genes": str(maf["Hugo_Symbol"].nunique()),
    }, tiles


def test_the_refusal_is_on_screen_and_not_behind_the_other_section():
    """Banners belong to the page.

    Every load path filters silently and stashes what it has to say, so a banner rendered
    by a section is a banner the user can navigate away from — or, once the page started
    jumping to Results on its own, one they are moved off before it is drawn.
    """
    app = _render(upload=str(UNREADABLE_MAF))

    banners = [banner.value for banner in app.error]
    assert any("cannot be filtered" in banner for banner in banners), banners
