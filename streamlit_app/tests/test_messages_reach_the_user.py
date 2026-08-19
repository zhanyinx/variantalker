"""The two sentences issue #140 moved, asserted where a user would read them.

``tests/test_discarded_frames.py`` next door holds the rule — that nothing is drawn into a
frame ``st.rerun`` throws away — by reading the source. A rule about what the app must *not*
do cannot say whether it still says anything at all: delete the drain call from the data
page and that guard stays green while *Saved for this session.* vanishes exactly as it did
before, which is the failure this ticket exists to repair. Found by review, and the reason
this module exists.

So these two boot the app and read the elements back:

* the note confirmation, parked by the variant dialog and drawn as a toast on the run after
  the rerun that closes it;
* a page that raises, reported **in place** — named as the sidebar names it, with the app
  still on the page the user asked for rather than bounced to Home.

The instrument has a known blind spot in this exact area, which is why the claims are
phrased the way they are. ``AppTest`` reports an element drawn before a rerun when it sits
inside a container (issue #133, measured with a one-element probe), so *the element exists*
is not by itself evidence a user saw it. Neither claim below rests on that: the toast is
asserted on a run with no rerun in it, together with the stash being **drained**, and the
page failure is asserted together with the app *not* having moved.
"""

from __future__ import annotations

import sys
from pathlib import Path
from unittest import mock

import pytest

STREAMLIT_APP = Path(__file__).resolve().parent.parent
if str(STREAMLIT_APP) not in sys.path:
    sys.path.insert(0, str(STREAMLIT_APP))

#: The sentence issue #67 chose, unchanged by issue #140 — only where it is said changed.
SAVED = "Saved for this session."

#: The data page alone, over a seeded session. The variant dialog is not driven here: it
#: cannot be, since ``AppTest`` cannot open a modal and click inside it, and it is not what
#: is under test either. What is under test is the seam between the two — that a sentence
#: parked by anything reaches the screen on the following run — so the test parks it the way
#: the dialog does, through the same public function.
_DATA_PAGE = """
import sys
sys.path.insert(0, {app!r})

import streamlit as st
from config.pipeline_params import pipeline_params
from page_modules import data_loading

if "filter_params" not in st.session_state or st.session_state.filter_params is None:
    st.session_state.filter_params = pipeline_params("somatic")

data_loading.show_data_loading_page()
"""

#: The whole app, so the router runs with a real sidebar and a real page beneath it.
_WHOLE_APP = """
import sys
sys.path.insert(0, {app!r})

import MAFigate

MAFigate.main()
"""


def _render_data_page(**seed):
    """The data page, rendered once, over a seeded session."""
    pytest.importorskip("streamlit", reason="streamlit not installed")
    from streamlit.testing.v1 import AppTest

    app = AppTest.from_string(_DATA_PAGE.format(app=str(STREAMLIT_APP)))
    for key, value in seed.items():
        app.session_state[key] = value
    app.run(timeout=180)
    assert not app.exception, [str(e.value) for e in app.exception]
    return app


def _render_whole_app(page, failing=None, **seed):
    """The app on ``page``, optionally with one of its calls replaced by one that raises.

    ``failing`` names any attribute of the ``MAFigate`` module, not only a page function.
    Issue #140 only ever broke a page, because ``route_to_page`` was what it was repairing;
    issue #144 needs the calls **outside** the router — ``render_header`` and its neighbours —
    since those are the only things that can still reach ``main()``'s own handler now that the
    router reports its own failures and lets nothing escape.

    The cache functions are stubbed for the reason ``test_app_defaults.py`` gives at length:
    the app consults ``~/.mafigate`` before falling back to the contract, so an unstubbed
    run is decided by whatever arm this developer last used — and ``discard_stale_cache``
    really does move that file, which would make this suite have a side effect on the
    machine it runs on.

    The call is replaced with ``mock.patch.object`` rather than by assignment inside the
    script, and that is not style: the script runs in *this* process, so an assignment there
    would outlive the test and leave the app permanently broken for everything that ran
    afterwards.
    """
    pytest.importorskip("streamlit", reason="streamlit not installed")
    from streamlit.testing.v1 import AppTest

    import MAFigate
    from page_modules import param_store

    def explode():
        raise RuntimeError("the page exploded")

    patches = [
        mock.patch.object(MAFigate, "load_parameters_from_cache", lambda: None),
        mock.patch.object(MAFigate, "get_cache_info", lambda: None),
        mock.patch.object(param_store, "discard_stale_cache", lambda: None),
    ]
    if failing is not None:
        patches.append(mock.patch.object(MAFigate, failing, explode))

    for patch in patches:
        patch.start()
    try:
        app = AppTest.from_string(_WHOLE_APP.format(app=str(STREAMLIT_APP)))
        app.session_state["current_page"] = page
        for key, value in seed.items():
            app.session_state[key] = value
        app.run(timeout=180)
        assert not app.exception, [str(e.value) for e in app.exception]
        return app
    finally:
        for patch in reversed(patches):
            patch.stop()


def _state(app, key, default=None):
    """``AppTest.session_state`` has no ``get``."""
    try:
        return app.session_state[key]
    except (KeyError, AttributeError):
        return default


# ---------------------------------------------------------------------------
# The note confirmation
# ---------------------------------------------------------------------------


def test_a_parked_note_confirmation_is_said_on_the_next_run():
    """The claim the whole change is for: the sentence reaches a user.

    Fails if the drain is dropped from the data page — which is the specific way this could
    silently regress, the guard next door being blind to a message that is never drawn at
    all.

    Seeded through :data:`components.variant_table._NOTE_CONFIRMATION` rather than a literal,
    so the test cannot go on passing against a key the app has renamed away from. The dialog
    itself is not driven: ``AppTest`` cannot open a modal and press a button inside it, and
    the seam under test is the one between the park and the draw.
    """
    from components.variant_table import _NOTE_CONFIRMATION

    app = _render_data_page(**{_NOTE_CONFIRMATION: SAVED})
    assert [toast.value for toast in app.toast] == [SAVED]


def test_the_confirmation_is_drained_so_it_is_not_said_twice():
    """Said once, on the run after the save, and not on every run thereafter.

    A toast redrawn on each later render would follow the user around the app, and a stash
    that is read rather than drained does exactly that.
    """
    from components.variant_table import _NOTE_CONFIRMATION

    app = _render_data_page(**{_NOTE_CONFIRMATION: SAVED})
    assert [toast.value for toast in app.toast] == [SAVED]
    assert _state(app, _NOTE_CONFIRMATION) in (None, "")

    app.run(timeout=180)
    assert not app.exception, [str(e.value) for e in app.exception]
    assert [toast.value for toast in app.toast] == []


def test_nothing_is_said_when_no_note_was_saved():
    """The ordinary render draws no toast, so the one above is evidence of something."""
    app = _render_data_page()
    assert [toast.value for toast in app.toast] == []


# ---------------------------------------------------------------------------
# A page that fails
# ---------------------------------------------------------------------------


def test_a_page_that_raises_says_so_where_the_user_is():
    """Named as the sidebar names it, with the detail a bug report needs one click away."""
    app = _render_whole_app("help", failing="show_help_page")

    reported = [element.value for element in app.error]
    assert any("❓ Help & Documentation" in text for text in reported), reported
    assert any("could not be opened" in text for text in reported), reported

    # The internal identifier is what the old message printed, and what the map's Style rule
    # keeps off the screen. It may appear in the traceback, which is folded; it may not be
    # the name the sentence uses.
    assert not any("'help'" in text for text in reported), reported

    assert any("RuntimeError" in block.value for block in app.code), [
        block.value for block in app.code
    ]


def test_a_page_that_raises_does_not_move_the_user_to_home():
    """The bounce is gone, so the sidebar and the page body still name the same page.

    The nav radio is drawn *before* the router, so a reset here could only take effect on a
    later run — which is why the old code needed a rerun, and why the rerun is what threw
    the explanation away.
    """
    app = _render_whole_app("help", failing="show_help_page")
    assert _state(app, "current_page") == "help"


def test_a_page_that_works_reports_nothing():
    """So the two above are evidence of the failure rather than of the harness."""
    app = _render_whole_app("home")
    assert not [
        text for text in (element.value for element in app.error) if "could not be opened" in text
    ]


# ---------------------------------------------------------------------------
# The app itself failing (issue #144)
# ---------------------------------------------------------------------------
#
# `render_header` is broken rather than a page, because a page can no longer reach
# `main()`'s handler at all — issue #140 gave `route_to_page` its own, and nothing escapes
# it. What is left is the four calls before the router and `render_into_status_slot` in the
# `finally` after it, and `render_header` is one of them.
#
# The unit tests in `test_components.py` own what the rescue *contains*; `AppTest` serves the
# button's payload over a URL and cannot read it back. These own the half those cannot: that
# a user in front of a broken app is offered it at all.


def _notes_for(**values):
    """A seeded ``variant_notes`` store, keyed the way the app keys one.

    Built through ``_variant_key`` rather than by writing the string out, so a test cannot go
    on passing against a spelling the app has moved away from.
    """
    import pandas as pd

    from components.variant_table import _variant_key

    row = pd.Series(
        {
            "Hugo_Symbol": "TP53",
            "Chromosome": "17",
            "Start_Position": "7577120",
            "Reference_Allele": "C",
            "Tumor_Seq_Allele2": "T",
        }
    )
    return {_variant_key(row): values.get("note", "hotspot, discussed at board")}


def test_a_broken_app_offers_the_users_writing_back():
    """The claim the ticket exists for.

    The advice of last resort used to be *refresh the page*, which starts the session over —
    and a note lives in session state and nowhere else (issue #67), so the app was telling the
    user to throw their writing away at the moment they were most likely to do it.
    """
    app = _render_whole_app("home", failing="render_header", variant_notes=_notes_for())

    offered = [button.label for button in app.get("download_button")]
    assert any("what you have written" in label.lower() for label in offered), offered

    said = [element.value for element in app.markdown]
    assert any("lives in this session only" in text for text in said), said
    assert any("download it first" in text for text in said), said


def test_a_broken_app_with_nothing_written_offers_no_download_and_no_alarm():
    """The gate, from the other side.

    A user who has written nothing loses nothing a refresh can destroy — the report re-derives
    from a file still on their disk — so they are told what a refresh does and not warned about
    it. Without this the pair above would pass over a handler that draws the button
    unconditionally.
    """
    app = _render_whole_app("home", failing="render_header")

    assert [button.label for button in app.get("download_button")] == []

    said = [element.value for element in app.markdown]
    assert any("starts MAFigate over" in text for text in said), said
    assert not any("lives in this session only" in text for text in said), said


def test_a_broken_app_no_longer_says_the_two_things_it_used_to():
    """`contact support` names nothing this repo has, and `str(e)` is Python on a clinician's
    screen — the traceback below already ends with that exact sentence."""
    app = _render_whole_app("home", failing="render_header")

    reported = [element.value for element in app.error] + [
        element.value for element in app.markdown
    ]
    assert not any("contact support" in text for text in reported), reported
    assert not any("An unexpected error occurred" in text for text in reported), reported
    assert not any("Application Error" in text for text in reported), reported

    assert any("MAFigate ran into a problem" in text for text in reported), reported


def test_a_broken_app_folds_the_traceback_away_for_a_bug_report():
    """The only part of this written in Python, and the only part that can get it fixed."""
    app = _render_whole_app("home", failing="render_header", variant_notes=_notes_for())

    assert any("RuntimeError" in block.value for block in app.code), [
        block.value for block in app.code
    ]


def test_an_app_that_works_says_none_of_this():
    """So the four above are evidence of the failure rather than of the harness."""
    app = _render_whole_app("home")

    said = [element.value for element in app.markdown] + [
        element.value for element in app.error
    ]
    assert not any("ran into a problem" in text for text in said), said
