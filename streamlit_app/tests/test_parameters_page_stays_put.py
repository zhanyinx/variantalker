"""Changing a filter on the parameters page must not move you to another page.

Issue #277 named it, issue #283 fixed it. The report: *selecting any item on the multiple
choice menu from any filters on the parameters page jumps to the data loading and analysis
page without recording the item*.

**What the fix was.** ``components/sidebar.py`` no longer calls ``st.rerun`` anywhere the file
chooser has not been drawn yet. Both reruns inside ``create_sidebar_navigation`` are gone —
neither had work to do, as the comments there now record — and the one in ``_nav_button``, which
is drawn *above* the chooser, is parked in ``NAV_RERUN_PENDING`` and taken by
``render_into_status_slot`` once the column is complete. Nothing is pruned, so the chooser's
state survives, so the browser's next message reports a file the session already holds and no
``on_change`` fires. That is the whole fix, and it is the third of the routes issue #283 listed:
it removes the class rather than the symptom, which is why ``_keep_a_section_selected`` on the
data page — the app's only other ``on_change``, and exposed to the same shape — is out of reach
of it too.

**The mechanism that was, named.** It was two of the ticket's suspects acting as one chain, and
neither half was wrong on its own:

1. ``create_sidebar_navigation`` called ``st.rerun()`` (for a nav-radio click, and for the
   reconcile pop) — and it did so *before* the sidebar's file chooser had been drawn, because
   the chooser lives in the status slot that ``render_into_status_slot`` fills at the end of
   the render. That call sits in ``main()``'s ``finally``, and the ``finally`` wraps
   ``route_to_page()`` alone, not the sidebar above it. So the aborted run never registered
   ``sidebar_maf_upload``. (Only the sidebar's reruns were exposed: ``RerunException``
   subclasses ``BaseException``, so a page's own ``st.rerun`` passes through
   ``route_to_page``'s ``except Exception`` and the ``finally`` draws the chooser anyway.)
2. Streamlit does not treat ``st.rerun`` as a premature stop
   (``scriptrunner/exec_code.py``: ``premature_stop = False`` on ``RerunException``), so it
   still runs ``on_script_finished(ctx.widget_ids_this_run)``, which drops the state of every
   widget the run did not register. The chooser's state went with it. The programmatic rerun
   that followed carried no widget states from the browser at all — ``st.rerun``'s
   ``RerunData`` has none — so the chooser re-registered at its initial value, ``None``.
3. The browser was still holding the file. The next message it sent — *any* widget
   interaction — reported that file against a stored value of ``None``, so
   ``SessionState._widget_changed`` was True and Streamlit fired the chooser's ``on_change``:
   ``_hand_the_chosen_file_to_the_data_page``, which sets ``current_page = "data_loading"``.
   Callbacks run before the script body, so the parameters page never rendered on that run and
   the user's selection was never written to ``filter_params``.

So the selection was lost **because of** the navigation, not independently of it — and it was
not specific to multiselects. ``test_a_number_change_does_not_move_you`` drives a
``number_input`` and it jumped the same way, losing the number. That is also why this is a
*different* defect from issue #276: this one lost numbers as readily as lists, and #276's
evidence is a cache in which every number survived and every list was emptied.

**The precondition.** The browser must be holding a file in the sidebar chooser — i.e. the
MAF was opened through that chooser rather than through the macOS ``MAFIGATE_OPEN_FILE``
route, which never populates it. ``test_no_file_held_no_jump`` is the control: with nothing
held, the same interactions stay put.

**Why the file is faked at the proto layer.** ``AppTest`` has no ``file_uploader`` widget, so
it never carries one in the widget states it sends. That is exactly the state under test —
what the *browser* keeps reporting between runs — so it is supplied here the way a browser
supplies it, as one more entry in the ``WidgetStates`` message, and the run is driven through
``AppTest._run`` with that message. Everything else in the message is ``AppTest``'s own.
``MemoryUploadedFileManager.get_files`` is patched to hand back the file those ids name,
which is the one piece of the upload path a test process has no way to have populated.

Two cautions this suite has already paid for elsewhere and honours here: ``filter_params`` is
seeded rather than inherited (the app boots on whatever ``~/.mafigate`` last held), and the
cache functions are patched so a test run cannot write to this developer's machine.
"""

from __future__ import annotations

import sys
from pathlib import Path
from unittest import mock

import pytest

STREAMLIT_APP = Path(__file__).resolve().parent.parent
if str(STREAMLIT_APP) not in sys.path:
    sys.path.insert(0, str(STREAMLIT_APP))

#: The label the sidebar's nav radio gives the parameters page, and the page behind it.
CONFIGURE = "⚙️ Configure Parameters"
PARAMETERS_PAGE = "parameter_config"

#: One row is enough. The file is never filtered here — what matters is that the browser is
#: holding *a* file, not what is in it.
MAF_TEXT = (
    "Hugo_Symbol\tChromosome\tStart_Position\tVariant_Classification\t"
    "t_alt_count\tt_ref_count\n"
    "TP53\t17\t7577120\tMissense_Mutation\t40\t60\n"
)

_WHOLE_APP = """
import sys
sys.path.insert(0, {app!r})

import MAFigate

MAFigate.main()
"""


class Browser:
    """The app, plus the file the browser's uploader is still holding.

    ``AppTest`` models the session; this models the one thing on the far side of the wire
    that ``AppTest`` cannot express. ``run()`` sends what a browser sends: every widget
    ``AppTest`` is showing, and — once ``holds_file`` is set — the chooser's file alongside
    them, under the element id Streamlit registered for it on the previous run.
    """

    FILE_ID = "issue-277-probe-file"

    def __init__(self, app):
        self.app = app
        self.holds_file = False

    def _uploader_id(self):
        from components.sidebar import UPLOAD_CHOOSER_KEY

        mapper = self.app.session_state._state._key_id_mapper
        return mapper.get_id_from_key(UPLOAD_CHOOSER_KEY, None)

    def _widget_states(self):
        from streamlit.proto.WidgetStates_pb2 import WidgetStates

        states = WidgetStates()
        states.widgets.extend(self.app._tree.get_widget_states().widgets)
        widget_id = self._uploader_id()
        if self.holds_file and widget_id is not None:
            held = states.widgets.add()
            held.id = widget_id
            info = held.file_uploader_state_value.uploaded_file_info.add()
            info.file_id = self.FILE_ID
            info.name = "issue_277.maf"
            info.size = len(MAF_TEXT.encode())
        return states

    def run(self):
        self.app._run(self._widget_states(), timeout=300)
        assert not self.app.exception, [str(e.value) for e in self.app.exception]
        return self.app

    @property
    def page(self):
        return self.app.session_state["current_page"]


@pytest.fixture
def browser():
    """The app booted on Home with nothing open, and the patches that keep it hermetic."""
    pytest.importorskip("streamlit", reason="streamlit not installed")

    from streamlit.runtime.memory_uploaded_file_manager import MemoryUploadedFileManager
    from streamlit.runtime.uploaded_file_manager import UploadedFileRec
    from streamlit.testing.v1 import AppTest

    import MAFigate
    from config.pipeline_params import pipeline_params
    from page_modules import param_store
    from page_modules import parameter_config

    record = UploadedFileRec(
        file_id=Browser.FILE_ID,
        name="issue_277.maf",
        type="text/plain",
        data=MAF_TEXT.encode(),
    )
    state = {"holds_file": lambda: False}

    def get_files(self, session_id, file_ids):
        return [record] if state["holds_file"]() else []

    patches = [
        mock.patch.object(MAFigate, "load_parameters_from_cache", lambda: None),
        mock.patch.object(MAFigate, "get_cache_info", lambda: None),
        mock.patch.object(param_store, "discard_stale_cache", lambda: None),
        mock.patch.object(parameter_config, "load_parameters_from_cache", lambda: None),
        mock.patch.object(parameter_config, "get_cache_info", lambda: None),
        mock.patch.object(parameter_config, "save_parameters_to_cache", lambda p: True),
        mock.patch.object(MemoryUploadedFileManager, "get_files", get_files),
    ]
    for patch in patches:
        patch.start()
    try:
        app = AppTest.from_string(_WHOLE_APP.format(app=str(STREAMLIT_APP)))
        # Seeded, never inherited: the real app opens on whatever the on-disk cache last
        # held, so an unseeded run measures this developer's last session instead.
        app.session_state["current_page"] = "home"
        app.session_state["filter_params"] = pipeline_params("somatic")
        app.run(timeout=300)
        assert not app.exception, [str(e.value) for e in app.exception]

        driver = Browser(app)
        state["holds_file"] = lambda: driver.holds_file
        yield driver
    finally:
        for patch in reversed(patches):
            patch.stop()


def _open_a_file_through_the_sidebar(browser):
    """Choose a MAF in the sidebar chooser and let the app settle on the data page.

    One run is the upload, which legitimately navigates (issue #64). The two after it are
    idle — the browser reporting the widgets it is already showing, which is what it sends on
    any interaction — and they are not padding: before issue #283 the chooser's ``on_change``
    fired a *second* time on the run after the upload, because that upload's own navigation
    rerun inside ``create_sidebar_navigation`` took the chooser's state with it. That second
    firing was harmless here only because it named the page already on screen. It is the same
    firing that moves you off the parameters page, so an idle run that stays put is worth
    asserting where it is cheapest to read.
    """
    browser.holds_file = True
    browser.run()
    assert browser.page == "data_loading"
    for _ in range(2):
        browser.run()
        assert browser.page == "data_loading"


def _go_to_the_parameters_page(browser):
    """Click the sidebar's nav radio, and confirm the parameters page is on screen."""
    browser.app.radio[0].set_value(CONFIGURE)
    browser.run()
    assert browser.page == PARAMETERS_PAGE
    assert "⚙️ Parameter Configuration" in [t.value for t in browser.app.title]


# Both symptom tests below carried `pytest.mark.xfail(strict=True)` while this module was the
# diagnosis of issue #277 and nothing had been fixed. Deleting the marker is part of issue
# #283's fix rather than tidying after it: strict xfail *fails* the moment the behaviour is
# repaired, so a fix that left it in place would report XPASS(strict) and turn a working app
# into a red suite. These are ordinary tests now, and they are what stops the defect coming
# back.


def test_a_filter_selection_does_not_move_you(browser):
    """The reported symptom: pick one item in a filter, and the app leaves the page."""
    _open_a_file_through_the_sidebar(browser)
    _go_to_the_parameters_page(browser)

    control = browser.app.multiselect[0]
    label, before = control.label, list(control.value)
    chosen = next(option for option in control.options if option not in before)
    control.select(chosen)
    browser.run()

    assert browser.page == PARAMETERS_PAGE, (
        f"selecting {chosen!r} in {label!r} moved the app to {browser.page!r}"
    )
    assert chosen in browser.app.session_state["filter_params"][
        "filter_variant_classification"
    ], f"{chosen!r} was not recorded"


def test_a_number_change_does_not_move_you(browser):
    """Not a multiselect defect: a number input goes the same way, and loses its number."""
    _open_a_file_through_the_sidebar(browser)
    _go_to_the_parameters_page(browser)

    control = browser.app.number_input[0]
    label = control.label
    control.set_value(51)
    browser.run()

    assert browser.page == PARAMETERS_PAGE, (
        f"changing {label!r} moved the app to {browser.page!r}"
    )
    assert browser.app.session_state["filter_params"]["min_depth"] == 51


def test_the_nav_radio_names_the_page_the_body_renders(browser):
    """What the deleted reconcile rerun was for, kept without it.

    Choosing a file navigates from a callback, so ``current_page`` moves without the radio
    being touched — and the radio ignores ``index`` once its widget key holds a value, which
    is why ``create_sidebar_navigation`` deletes that key when it finds the two disagreeing.
    It used to delete and *rerun*, and that rerun is half of issue #277. Deleting alone works
    because ``SessionState.__delitem__`` drops the browser's value for the key as well, so the
    registration further down the same run has nothing left to prefer over ``index``.

    That last sentence is a fact about Streamlit rather than about this app, which is exactly
    why it is asserted here: if a future version keeps the sent value, this app's sidebar and
    page body would quietly start naming different pages, and nothing else in the suite would
    notice.
    """
    from components.sidebar import PAGE_LABELS

    _open_a_file_through_the_sidebar(browser)
    assert browser.app.radio[0].value == PAGE_LABELS["data_loading"], (
        "the file chooser navigated to the data page, but the nav radio still reads "
        f"{browser.app.radio[0].value!r}"
    )

    _go_to_the_parameters_page(browser)
    assert browser.app.radio[0].value == PAGE_LABELS[PARAMETERS_PAGE]


def test_the_chooser_hands_its_file_over_once(browser):
    """The mechanism under the two tests above, counted rather than inferred.

    Choosing a file is one event, so its callback runs once — not once per interaction the
    user makes afterwards. This is what the two symptom tests observe from the outside, and
    counting it here is what tells the next reader *why* they failed: a second firing is a
    pruned widget re-arriving, never a second choice.

    The wrapper is installed before an idle run rather than straight before the upload,
    because Streamlit calls the callback recorded when the widget was last *drawn*. Patching
    and uploading in one go would let the previous run's recording answer the first call and
    the count would be short by one — the reading that would make this test lie.
    """
    import components.sidebar as sidebar

    real = sidebar._hand_the_chosen_file_to_the_data_page
    fired = []

    def counting_the_hand_over():
        fired.append(browser.page)
        real()

    with mock.patch.object(
        sidebar, "_hand_the_chosen_file_to_the_data_page", counting_the_hand_over
    ):
        browser.run()
        assert fired == [], "nothing was chosen, so nothing should have been handed over"

        _open_a_file_through_the_sidebar(browser)
        assert len(fired) == 1, f"one file chosen, {len(fired)} hand-overs"

        _go_to_the_parameters_page(browser)
        control = browser.app.multiselect[0]
        control.select(next(o for o in control.options if o not in control.value))
        browser.run()

        assert len(fired) == 1, (
            f"the chooser handed its file over again on run {len(fired)}, from "
            f"{fired[-1]!r} — the file was chosen once and nothing has chosen another"
        )


def test_the_route_back_still_navigates_from_below_the_chooser(browser):
    """The third rerun site, which could not be deleted — only moved.

    ``_nav_button`` draws the sidebar's route back to your file, and it is drawn *above* the
    chooser inside the status block. Its rerun is the one this fix could not do without: the
    button is drawn after the page has already rendered, so setting ``current_page`` alone
    would leave the new page unshown until the user's next interaction. So the request is
    parked in ``NAV_RERUN_PENDING`` and taken at the end of ``render_into_status_slot``, and
    the two things worth asserting are that it is still taken — one click, the data page on
    screen — and that the parking key does not survive the run that honours it.

    A file is seeded as open rather than opened, because the button is drawn only in that
    state and this module's MAF is one row wide on purpose: it is refused by the data page's
    column validation, which is all the navigation under test needs and less than the status
    block needs.
    """
    import pandas as pd

    from components.sidebar import NAV_RERUN_PENDING

    browser.app.session_state["maf_data"] = pd.DataFrame({"Hugo_Symbol": ["TP53"]})
    browser.app.session_state["maf_source_name"] = "seeded.maf"
    browser.app.session_state["current_page"] = PARAMETERS_PAGE
    browser.run()

    back = [button for button in browser.app.button if "Back to your" in button.label]
    assert len(back) == 1, (
        f"expected the sidebar's one route back with a file open, found {len(back)}: "
        f"{[b.label for b in back]}"
    )

    back[0].click()
    browser.run()

    assert browser.page == "data_loading"
    assert "📊 Data Loading & Analysis" in [t.value for t in browser.app.title], (
        "the button set the page but nothing put it on screen — the parked rerun was not "
        "taken, so the user is looking at the page they asked to leave"
    )
    assert NAV_RERUN_PENDING not in browser.app.session_state, (
        "the parked rerun request outlived the run that honoured it, which would rerun every "
        "render from here on"
    )


def test_no_file_held_no_jump(browser):
    """The control case that fixes the precondition: with no file held, nothing moves.

    Not a weaker copy of the two above — it is what makes their failure mean *the chooser's
    callback*, rather than anything else the parameters page does on a widget change.
    """
    browser.app.radio[0].set_value(CONFIGURE)
    browser.run()
    assert browser.page == PARAMETERS_PAGE

    control = browser.app.multiselect[0]
    chosen = next(o for o in control.options if o not in control.value)
    control.select(chosen)
    browser.run()

    assert browser.page == PARAMETERS_PAGE
    assert chosen in browser.app.session_state["filter_params"][
        "filter_variant_classification"
    ]
