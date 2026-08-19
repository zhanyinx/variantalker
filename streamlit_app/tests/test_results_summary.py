"""What the Summary tab says about the user's columns (issue #90).

The tab used to end with a **Column Analysis** block: two columns of ``✅``/``❌`` MAF
column names, ten at a time, under an *Available Key Columns (39)* heading. It was deleted,
and this module holds the reason so it cannot quietly come undone.

**It measured against the wrong list.** The block read
:data:`config.columns.PRIORITY_COLUMNS`, which is a presentation *order* — what the grid
shows when the user ticks "Show all columns" — and not
:func:`config.columns.resolve_visible_columns`, the list that selects the table and the
export. The difference is not academic, because the resolver asks the *pipeline* what it
would emit for this arm and this file, via the vendored ``compute_keep``, and the ordering
list asks nothing: it is one flat list of somatic-flavoured names for both arms.

So the block reported absence where the pipeline had made a *choice*. A MAF with no CIViC
annotation is not an incomplete MAF — ``compute_keep`` sees the CIViC columns are absent
and drops them from its own output list, which is why the pipeline prints
``CIViC columns will be automatically excluded from output`` rather than complaining. The
block knew nothing of that branch and drew six red ``❌`` anyway. The germline arm is worse
again: ``CancerVar``, ``cosmic``, ``Freq_ExAC_ALL`` and the six ``CIViC_*`` are columns the
germline pipeline does not emit at all, so a complete germline MAF was told it was missing
nine "key" columns.

That is what :func:`test_the_reference_maf_is_reported_complete` pins, and it is the whole
case in one line: **on both of this repo's reference MAFs the app's own verdict is that
nothing is missing, and the block drew six and nine ``❌`` at them respectively.** Every
marker it drew on the files this suite is built from was a false alarm.

Its count was rigged in the other direction too, though only the code shows it now:
``Clinical_Summary`` heads the ordering list and is derived by the app at filter time, so
it counted as present whatever the file said — one guaranteed ``✅`` in every denominator.

Nothing replaced the block. The app answers the question correctly on the same screen when there is
something to answer — :func:`test_a_genuinely_absent_column_is_named` renders a MAF with one — and
``filters/absent_columns.py`` escalates separately when a neutral fill would move the verdict. A
fresh line here would also have prejudged #93, which is deciding how many times that one sentence
should be said.

Why the warning is asserted by presence, not by count
-----------------------------------------------------
``_visible_columns`` surfaces the resolver's warning on every call and the results section
calls it several times — which is exactly what #93 owns. This module asserts only that the
sentence is said **at least once**, naming the column and the arm. Whichever way #93 goes,
the claim tested here still holds and this file does not have to move with it.

A unit module: nothing here has a pipeline counterpart — the pipeline draws no tabs — and
none of it needs ``bin/``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

STREAMLIT_APP = Path(__file__).resolve().parent.parent
if str(STREAMLIT_APP) not in sys.path:
    sys.path.insert(0, str(STREAMLIT_APP))

MAF = STREAMLIT_APP / "tests" / "fixtures" / "parity" / "germline_reference.maf"

#: The germline reference is used rather than the somatic one every other page test
#: reaches for, because it is the file the deleted block was worst on: nine ``❌`` against
#: a file the app considers complete.
ARM = "germline"

#: A column the germline pipeline emits unconditionally, dropped from the MAF in the second
#: render so there is a real absence to report. AlphaMissense's class rather than
#: ``InterVar`` or ``RENOVO_Class``, which are guideline sources feeding pathogenic
#: retention: those take ``filters/absent_columns.py``'s escalated path as well, and this
#: module is about the *table and export* report, not the escalation.
DROPPED = "am_class"

#: Wording the deleted block owned. ``Column Analysis`` was its heading; ``Key Columns``
#: appeared in both subheadings and is a phrase the interface defines nowhere — which was
#: the copy half of the complaint, separate from the block being wrong.
BLOCK_PHRASES = ("Column Analysis", "Key Columns", "key columns")

#: Names the block drew ``❌`` against on a complete germline MAF. Used to recognise a
#: roll-call independently of whatever heading it might come back under.
FALSELY_FLAGGED = ("CancerVar", "cosmic", "Freq_ExAC_ALL", "CIViC_Evidence_Level")

_SCRIPT = """
import os
import sys
sys.path.insert(0, {app!r})

import streamlit as st

from config.pipeline_params import pipeline_params
from page_modules import data_loading
from utils import read_maf

maf = read_maf({maf!r})
drop = {drop!r}
if drop:
    assert drop in maf.columns, drop
    maf = maf.drop(columns=[drop])

st.session_state.filter_params = pipeline_params({arm!r})
st.session_state.maf_source_name = os.path.basename({maf_path!r})
st.session_state.maf_data = maf
data_loading.apply_filters_to_data(show_messages=False)

if {whole_page!r}:
    st.session_state["data_page_section"] = "results"
    data_loading.show_data_loading_page()
else:
    data_loading.show_results_section()
"""


def _render(drop=None, whole_page=False):
    """The results section over the germline reference, rendered once.

    Streamlit renders every tab's body whether or not it is the open one, so the Summary
    tab's content is in this tree without the harness having to select it — just as well,
    since ``AppTest`` cannot drive the page's section control across runs
    (``test_data_page_sections.py`` records why).

    ``whole_page`` draws the page around the section instead, which one test below needs:
    since issue #93 the missing-column sentence is collected by the grids and said once by
    the *page*, so a section rendered on its own no longer draws it anywhere.
    """
    pytest.importorskip("streamlit", reason="streamlit not installed")
    from streamlit.testing.v1 import AppTest

    app = AppTest.from_string(
        _SCRIPT.format(
            app=str(STREAMLIT_APP),
            maf=str(MAF),
            maf_path=str(MAF),
            arm=ARM,
            drop=drop,
            whole_page=whole_page,
        )
    )
    app.run(timeout=180)
    assert not app.exception, [str(e.value) for e in app.exception]

    # The page turns its own render failures into `st.error`, so an exception-free run is
    # not by itself a rendered results section.
    assert not app.error, [element.value for element in app.error]
    return app


@pytest.fixture(scope="module")
def complete():
    """The reference MAF as it ships — a file the app considers complete."""
    return _render()


@pytest.fixture(scope="module")
def missing_a_column():
    """The same file with one column the germline pipeline emits taken out of it.

    Drawn with the page around it, unlike its sibling above, because the one thing this
    fixture exists to show is now said by the page (issue #93).
    """
    return _render(drop=DROPPED, whole_page=True)


def _texts(app):
    """Every string the section drew, across the element kinds it draws them with."""
    values = []
    for kind in ("markdown", "text", "warning", "info", "caption", "subheader"):
        values.extend(str(element.value) for element in getattr(app, kind))
    return values


def test_the_summary_tab_actually_rendered(complete):
    """The anchor under every absence test below, which is worthless without it.

    "The block is not on screen" passes just as well when *nothing* is on screen — a
    section that returned early, a fixture that filtered to nothing, a harness that drew a
    blank page. So the tab is pinned by two things it does still draw: the metrics that
    open it, and the run report beneath them, the only artifact anywhere recording the
    parameters a cut was made with.
    """
    labels = [str(element.label) for element in complete.metric]
    assert "Total Variants" in labels, labels
    assert "Pass Rate" in labels, labels

    reports = [
        element
        for element in complete.get("download_button")
        if "report" in str(element.label).lower()
    ]
    assert reports, "the Summary tab drew no run report download"


def test_the_summary_tab_inventories_no_columns(complete):
    """The deletion itself: no ``✅``/``❌`` roll-call of column names, and no heading.

    Checked against every string the section drew rather than against the source of
    ``results_view``, because a guard that reads a module cannot see what the module
    renders — issue #71's ``not hasattr(help, "APP_VERSION")`` passed over a version
    written out by hand two lines away, and issue #79 made rendering the page the rule.
    """
    drawn = _texts(complete)
    for phrase in BLOCK_PHRASES:
        offenders = [text for text in drawn if phrase in text]
        assert not offenders, f"{phrase!r} is back on the Summary tab: {offenders}"

    roll_call = [
        text
        for text in drawn
        if ("✅" in text or "❌" in text)
        and any(name in text for name in FALSELY_FLAGGED)
    ]
    assert not roll_call, f"a column roll-call is back on the Summary tab: {roll_call}"


def test_the_reference_maf_is_reported_complete(complete):
    """The finding that condemned the block, stated as the app's own verdict.

    This file is complete as far as the germline pipeline is concerned: ``compute_keep``
    asks what the arm emits for the columns actually present, so the absent CIViC block is
    a branch it takes rather than a loss it suffers. The app says nothing, correctly. The
    deleted block drew **nine** red ``❌`` at this same file — every one of them a name the
    germline pipeline never emits.

    So this is not a test that the app is quiet; it is a test that the app is quiet *here*,
    which is the state the block could not represent.

    Scoped to the *missing-column* claim rather than to the somatic names themselves. When
    this was written the grid's Pathogenicity Overview key named all six sources —
    ``CancerVar`` and ``CIViC`` included — on a germline run, where two of its glyph
    positions could only ever be ``⬜``, so a broader assertion would have failed on a key
    that was not making this block's mistake. Issue #95 has since narrowed that column to
    the sources the arm emits and the file carries, and ``test_overview_legend.py`` holds
    it there; the scoping stays because it is the right scope for *this* claim, not
    because the key is still arm-blind.
    """
    complaints = [
        text
        for text in _texts(complete)
        if "missing" in text.lower() and "column" in text.lower()
    ]
    assert not complaints, (
        "the reference MAF is reported as missing columns, which it is not: "
        f"{complaints}"
    )


def test_a_genuinely_absent_column_is_named(missing_a_column):
    """Nothing true was lost — the claim the deletion rests on, tested rather than argued.

    With a column the germline pipeline really does emit taken out of the file, the
    resolver says so by name and names the arm it measured against. That is the report the
    block was a worse copy of: arm-correct and untruncated.

    Issue #90 predicted this test would survive issue #93 "either way", on the grounds
    that it asserts presence rather than count. Half right — the count is untouched, but
    #93 moved *where* the sentence is drawn, out of the tabs and up to the page-level slot,
    and presence is read off a screen. So the fixture draws the page now. The prediction
    was about the assertion and the thing that moved was the harness under it.
    """
    # Recognised by the column it names, not by a word from the sentence. It used to match
    # on "missing", which #93's rewording dropped — and a predicate transcribing the copy
    # fails when the copy is improved rather than when the claim stops being made. Nothing
    # else on this screen can name `am_class`: it is emitted by the pipeline and read by no
    # filter, which is why this module picked it.
    said = [text for text in _texts(missing_a_column) if DROPPED in text]
    assert said, (
        f"{DROPPED} was removed from the MAF and dropped from the table, and nothing on "
        "screen says so"
    )
    assert ARM in said[0], (
        "the missing-column warning does not name the arm it measured against, which is "
        f"the whole of why the deleted block was wrong here: {said[0]!r}"
    )


def test_the_fixture_is_the_file_this_module_claims_it_is(complete):
    """The docstring's germline facts, checked instead of trusted.

    A fixture quietly regenerated with different columns would leave every test above
    passing while testing nothing of what they describe — so the property the reason for
    deletion rests on is asserted here: the somatic-only names the block drew ``❌``
    against really are absent from this file, and really are not the file's fault.
    """
    columns = set(complete.session_state["maf_data"].columns)
    present = sorted(name for name in FALSELY_FLAGGED if name in columns)
    assert not present, (
        "the germline fixture now carries columns the block used to flag, so this "
        f"module no longer demonstrates what it describes: {present}"
    )
