"""What a download actually carries (issue #67).

A note the user types is never written to disk. That is a decision, not an oversight: the
real store is institute-wide and hosted, and a per-laptop one built first would have the
wrong reach and leave notes on N machines for the hosted store to gather. So the app says
plainly that notes live for the session and points at the download — which makes the
download the *only* way a note survives the tab being closed, and makes a download that
silently drops one the single failure this app cannot afford.

Three of the seven download buttons on the results page were dropping user-authored
columns when this module was written:

* the two under **Download Results** serialised ``st.session_state.filtered_data``
  directly, and ``Notes`` is built onto a *copy* at display time, so that frame has never
  carried it — those two exported every pipeline column and none of the user's;
* the two *shown columns* buttons called the column resolver alone, while the grid called
  it and then appended the user's invented annotation columns, so a column the user named
  was on screen and missing from the CSV offered directly beneath it.

``test_components.py`` owns the two seams that fixed this (``with_user_columns`` and
``_shown_columns``). What this module owns is the *wiring*: the page really does route its
buttons through them. That was the whole of the bug — both seams existed and the page did
not use them — so a test of the seams alone would have passed against the broken app.

Why the payload is read from a stub
-----------------------------------
``AppTest`` exposes a download button but not its bytes: Streamlit hands the payload to
its media file manager and the element carries a URL. Rather than reach into that, the
script below replaces ``st.download_button`` with a recorder, exactly as it already
replaces ``st.file_uploader`` — the same technique this suite's page tests use, for the
same reason. The recorded ``data`` is the string the user's browser would save.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

STREAMLIT_APP = Path(__file__).resolve().parent.parent
if str(STREAMLIT_APP) not in sys.path:
    sys.path.insert(0, str(STREAMLIT_APP))

MAF = STREAMLIT_APP / "tests" / "fixtures" / "parity" / "somatic_reference.maf"

#: The note text, and the name of a column the user invented. Both are searched for
#: verbatim in the CSV, so they are deliberately unlike anything a MAF carries.
NOTE = "discussed at molecular board, keep"
INVENTED_COLUMN = "Board_Outcome"
INVENTED_VALUE = "escalate"

#: The results section, rendered **once**, over a session seeded with a filtered MAF and a
#: note on every variant.
#:
#: Rendered once because Streamlit keys widgets per run: calling the section twice — load
#: and filter, then annotate and redraw — re-registers ``download_passed_shown`` and the
#: duplicate-key error is *caught* by the page's own ``except`` around the results table, so
#: the run ends early having reported nothing. The load is therefore done directly rather
#: than through the page, which needs no uploader stub either.
#:
#: A note on *every* variant, because the passed and failed tables are different frames and
#: only one of them holds any given variant — noting one would leave the other's downloads
#: asserting nothing. That the notes land on the right variants rather than smeared across
#: all of them is ``test_components.py``'s job, where one row is noted and the next is
#: checked to be empty.
_SCRIPT = """
import os
import sys
sys.path.insert(0, {app!r})

import streamlit as st

from components.variant_table import _variant_key
from config.pipeline_params import pipeline_params
from page_modules import data_loading
from utils import read_maf

# Unwrapped, not read straight off the module: `import streamlit` is cached in
# sys.modules, so a second AppTest in the same pytest process starts with the *previous*
# run's recorder already installed. Taken plainly, each run would wrap the last and every
# button would be recorded once per run so far.
_real_download_button = getattr(st.download_button, "_wraps", st.download_button)


def _recording_download_button(label, data=None, **kwargs):
    # Keyed by the button's own `key` where it has one: the two tabs reuse their labels,
    # and only the key tells "shown columns" on the passed table from the failed one.
    st.session_state.setdefault("recorded_downloads", []).append(
        {{
            "identity": kwargs.get("key") or label,
            "label": label,
            "file_name": kwargs.get("file_name"),
            "data": data,
        }}
    )
    return _real_download_button(label, data=data, **kwargs)


_recording_download_button._wraps = _real_download_button
st.download_button = _recording_download_button

# What the grid decided to show, per tab, recorded off `create_data_table`'s return value
# (issue #92). Patched on `results_view`, which imported the name at import time, so this
# is the object that module actually calls. The test compares this against the header of
# the CSV beside it: two code paths that must agree, rather than the test recomputing what
# the grid ought to have shown and checking the app agrees with the test.
from components import results_view as _results_view

_real_create_data_table = getattr(
    _results_view.create_data_table, "_wraps", _results_view.create_data_table
)


def _recording_create_data_table(data, height=None, key_suffix=None, *args, **kwargs):
    columns = _real_create_data_table(data, height, key_suffix, *args, **kwargs)
    st.session_state.setdefault("grid_columns", {{}})[key_suffix] = list(columns)
    return columns


_recording_create_data_table._wraps = _real_create_data_table
_results_view.create_data_table = _recording_create_data_table

st.session_state.filter_params = pipeline_params("somatic")
st.session_state.maf_source_name = os.path.basename({maf!r})
st.session_state.maf_data = read_maf({maf!r})
data_loading.apply_filters_to_data(show_messages=False)

maf = st.session_state.maf_data
keys = [_variant_key(row) for _, row in maf.iterrows()]
st.session_state["variant_notes"] = {{key: {note!r} for key in keys}}
st.session_state["custom_annotations"] = {{{column!r}: {{key: {value!r} for key in keys}}}}

{seed}

data_loading.show_results_section()
"""

#: Ticking *Show all columns* on the passed table, before it is drawn.
#:
#: The grid's two controls are keyed widgets, so seeding the key is the same lever the
#: user pulls, and it happens in the first render — no ``AppTest`` interaction and re-run,
#: which would draw the results section twice in one run and trip the duplicate-widget-key
#: error the page's own ``except`` swallows. One render per state, one ``AppTest`` each;
#: what the runs share is the ``streamlit`` module, which is why the stub above unwraps.
#:
#: Only the passed table is ticked. The failed one is left alone deliberately, so the same
#: render shows a collapsed pair beside an uncollapsed one and the assertion is about the
#: state rather than about the page.
_SEED_SHOW_ALL = """
st.session_state["show_all_passed_variants"] = True
"""

#: Adding one column through *Add columns* on the passed table.
#:
#: The column is chosen off the frame rather than written down: it must be one the MAF
#: carries and the resolver does not return, and naming a real MAF column here would make
#: this test a second copy of the resolver's list. ``added_column`` is read back out of
#: session state by the test.
_SEED_ADD_COLUMN = """
from components.variant_table import shown_columns, with_user_columns

_exported = with_user_columns(st.session_state.filtered_data)
_default = shown_columns(_exported)
_pool = [c for c in _exported.columns if c not in _default]
assert _pool, "every column is already shown; this fixture has nothing to add"
st.session_state["added_column"] = _pool[0]
st.session_state["extra_cols_passed_variants"] = [_pool[0]]
"""


def _run(seed: str = ""):
    """Render the results section once, with ``seed`` executed just before it."""
    pytest.importorskip("streamlit", reason="streamlit not installed")
    from streamlit.testing.v1 import AppTest

    app = AppTest.from_string(
        _SCRIPT.format(
            app=str(STREAMLIT_APP),
            maf=str(MAF),
            note=NOTE,
            column=INVENTED_COLUMN,
            value=INVENTED_VALUE,
            seed=seed,
        )
    )
    app.run(timeout=180)
    return app


def _collect(app):
    """Every download the section offered, as ``{identity: {label, file_name, data}}``."""
    assert not app.exception, [str(e.value) for e in app.exception]

    # The page catches its own render errors and turns them into `st.error`, so an
    # exception-free run is not by itself a rendered results table.
    assert not app.error, [element.value for element in app.error]

    recorded = app.session_state["recorded_downloads"]
    assert recorded, "the results section offered no downloads at all"
    names = [button["identity"] for button in recorded]
    assert len(names) == len(set(names)), f"two downloads share an identity: {names}"
    return {button["identity"]: button for button in recorded}


def _buttons(seed: str = ""):
    return _collect(_run(seed))


@pytest.fixture(scope="module")
def buttons():
    return _buttons()


@pytest.fixture(scope="module")
def downloads(buttons):
    """Just the payloads, as ``{identity: bytes-or-str the browser would save}``."""
    return {name: button["data"] for name, button in buttons.items()}


def _variant_tables(downloads):
    """The downloads that are tables of variants, found by their header.

    Selected by content rather than by label: the page also offers a plain-text summary,
    and picking the variant tables out by name would need the same list of button labels
    this module exists to stop trusting.
    """
    return {
        name: payload
        for name, payload in downloads.items()
        if isinstance(payload, str) and "Hugo_Symbol" in payload.splitlines()[0]
    }


def test_every_variant_download_carries_the_users_note(downloads):
    """The claim the app now makes in the variant dialog, tested where it is kept.

    "Notes live for this session only: download the table to keep them" is a promise about
    *the* download, and the page offers several. A user who reads that line and presses the
    nearest button has kept their notes or has not; which button it was cannot be part of
    the answer.
    """
    tables = _variant_tables(downloads)
    # Four: shown/all × passed/failed. An equality since issue #83, which made how many
    # doors lead to a table a decision this module owns — see the test below.
    assert len(tables) == 4, f"expected the page's variant tables, got {list(tables)}"

    missing = [name for name, payload in tables.items() if NOTE not in payload]

    assert not missing, (
        f"these downloads dropped the user's note: {missing}. Every one of them must "
        "route through components.variant_table.with_user_columns — see issue #67"
    )


def test_every_variant_download_names_the_notes_column(downloads):
    """The header, not just the text — an empty ``Notes`` column is still a lost note.

    The passed and failed tables are different frames, and only one of them holds the
    annotated variant. Searching for the note text alone would let the other export the
    column away entirely and still pass by never being asked.
    """
    headerless = [
        name
        for name, payload in _variant_tables(downloads).items()
        if "Notes" not in payload.splitlines()[0].split(",")
    ]

    assert not headerless, f"these downloads have no Notes column: {headerless}"


def test_a_column_the_user_invented_is_downloaded_too(downloads):
    """The narrower half of the bug, and the one the earlier fix created.

    Pinning the user's own annotation columns into the grid (issue #60) made them visible
    without making them exportable: the resolver cannot name them — it is streamlit-free
    and they live in session state — so the "shown columns" download, which called the
    resolver alone, showed one thing and wrote another.
    """
    missing = [
        name
        for name, payload in _variant_tables(downloads).items()
        if INVENTED_COLUMN not in payload.splitlines()[0].split(",")
    ]

    assert not missing, (
        f"these downloads dropped the user's own {INVENTED_COLUMN} column: {missing}"
    )


# --------------------------------------------------------------------------------------
# How many doors lead to the same file (issue #83)
# --------------------------------------------------------------------------------------


def test_the_results_section_offers_exactly_these_downloads(buttons):
    """Five buttons, and the list is the decision — not an inventory of what accreted.

    There were seven. Two of them, under a `💾 Download Results` block below both tabs,
    served the *same bytes* as the tabs' "all columns" buttons and differed only in the
    name the file landed under — so a user wanting their passed variants met three doors
    to two files. Naming the survivors here means a sixth arrives by decision.
    """
    assert set(buttons) == {
        "download_passed_shown",
        "download_passed_all",
        "download_failed_shown",
        "download_failed_all",
        "download_summary",
    }, f"the results section's downloads have changed: {sorted(buttons)}"


def test_no_two_downloads_hand_over_the_same_bytes(downloads):
    """The defect itself, rather than the count that happened to expose it.

    Counting buttons would pass a page that grew an eighth button duplicating a seventh.
    What made the old block indefensible was that two of its buttons were byte-for-byte
    another two — measured, not argued: identical SHAs over 351,331 and 463,492 bytes.
    """
    by_payload = {}
    for name, payload in downloads.items():
        by_payload.setdefault(payload, []).append(name)

    duplicates = [names for names in by_payload.values() if len(names) > 1]

    assert not duplicates, (
        f"these downloads carry identical bytes and differ only in name: {duplicates}. "
        "Two buttons offering one file is what issue #83 removed."
    )


def test_each_label_states_the_column_count_it_delivers(buttons):
    """The number in the label is measured off the frame, so it cannot be wrong — unless
    it is written down somewhere, which is what this stops.

    Nothing on screen told a user that "all columns" was nine times wider than "shown
    columns"; the counts are in the labels for that reason (issue #83), which only helps
    if they are the counts of the file that actually arrives.
    """
    for name, button in buttons.items():
        if not name.startswith("download_passed") and not name.startswith(
            "download_failed"
        ):
            continue
        delivered = len(button["data"].splitlines()[0].split(","))
        assert str(delivered) in button["label"], (
            f"{name} is labelled {button['label']!r} but hands over {delivered} columns"
        )


# Issue #83 also asserted here that no variant label contains the word "shown". That ban
# is gone, and deliberately: it was right for the app it was written against, where the
# grid's set was the user's to change and the download was always the resolver's list, so
# "Download shown columns" was a lie for exactly as long as nobody checked. Issue #92 made
# the download follow the grid, which makes the word true — and replaces a ban on saying
# it with a test that it *holds*, below. Policing the wording was the weaker guard: it
# would have passed a label reading "Download the 48 columns you can see".


def test_every_download_is_named_after_the_loaded_file(buttons):
    """A clinician works several samples in a sitting.

    Four files called `passed_variants_all.csv` are four files they cannot tell apart in
    a Downloads folder, and the browser's `(1)` names nothing. The app has recorded the
    source file's name since issue #58; issue #83 is what made the downloads read it.
    """
    unnamed = [
        name
        for name, button in buttons.items()
        if not button["file_name"].startswith("somatic_reference_")
    ]

    assert not unnamed, (
        f"these downloads are not named after the loaded MAF: "
        f"{[(name, buttons[name]['file_name']) for name in unnamed]}"
    )


@pytest.mark.parametrize("dataset", ["passed", "failed"])
def test_the_shown_columns_download_is_narrower_than_the_all_columns_one(
    downloads, dataset
):
    """The two buttons must still mean different things.

    A tempting fix for the above is to hand every button the full frame, which would pass
    all three tests by making "shown columns" a lie. The distinction is the point: one
    download is the columns in front of the user, the other is everything the MAF carries.

    Issue #83 kept both buttons on exactly this ground — 48 columns against 434 is a real
    choice, where the block below the tabs was offering the same 434 twice — so this now
    runs over both tables rather than the passed one alone, and checks the narrow file is
    a *subset* rather than merely shorter.
    """
    shown_header = downloads[f"download_{dataset}_shown"].splitlines()[0].split(",")
    all_header = downloads[f"download_{dataset}_all"].splitlines()[0].split(",")

    assert len(shown_header) < len(all_header)
    assert set(shown_header) <= set(all_header)
    # ...and the narrower one still keeps what the user wrote.
    assert "Notes" in shown_header
    assert INVENTED_COLUMN in shown_header


# --------------------------------------------------------------------------------------
# Whether the download follows the grid or the report (issue #92)
# --------------------------------------------------------------------------------------


def test_a_column_added_through_the_multiselect_reaches_the_csv():
    """The complaint this ticket opened on, driven rather than argued.

    A user searches *Add columns* for a column they need, adds it, sees it in the grid,
    presses the button directly beneath — and before issue #92 did not find it in the
    file. That is issue #60's shape one layer out: fixed there for the columns a user
    *invents* in the variant dialog, which ``shown_columns`` appends, and not for the
    pipeline columns they *add* through the multiselect, which nothing did.

    The added column is chosen off the frame by the seed, not named here, so this cannot
    quietly become a test that one particular MAF column is exported.
    """
    app = _run(_SEED_ADD_COLUMN)
    buttons = _collect(app)
    added = app.session_state["added_column"]

    header = buttons["download_passed_shown"]["data"].splitlines()[0].split(",")

    assert added in header, (
        f"the user added {added!r} through 'Add columns' and the download beneath the "
        f"grid does not carry it: {header[:5]}…"
    )
    # The failed tab's grid was not touched, so its download must not have moved either —
    # the two tables keep their own column choices, as their widget keys already do.
    failed_header = buttons["download_failed_shown"]["data"].splitlines()[0].split(",")
    assert added not in failed_header, (
        "adding a column to the passed grid changed the failed table's download"
    )


@pytest.mark.parametrize("seed", ["", _SEED_ADD_COLUMN], ids=["as-opened", "column-added"])
@pytest.mark.parametrize("dataset", ["passed", "failed"])
def test_the_shown_download_carries_exactly_the_grids_columns(seed, dataset):
    """The guard that replaces issue #83's ban on the word "shown".

    The label may now say what is on screen, so something has to hold that true — and this
    is the stronger form. It compares the CSV's header against the list
    ``create_data_table`` returned for that tab on that render: two paths through the app
    that must agree, where the ban only checked that a word did not appear and would have
    passed a label reading "the 48 columns you can see".

    Neither side is recomputed here. Asking ``resolve_visible_columns`` what the grid shows
    would be this module answering the question it exists to check — and it cannot anyway,
    since the user's own annotation columns live in a session state a bare test process
    does not have.

    An exact list equality, so it holds the *order* too: the grid resolves the pipeline's
    order with additions appended, and issue #60 keeps the export in that order by putting
    ``_PINNED_LEFT`` downstream of the resolver. Both states are driven because the file
    the app opens with and the file after a user adds a column fail differently — drop the
    ``grid_columns`` argument and only the second notices.
    """
    app = _run(seed)
    buttons = _collect(app)

    grid = app.session_state["grid_columns"][f"{dataset}_variants"]
    header = buttons[f"download_{dataset}_shown"]["data"].splitlines()[0].split(",")

    assert header == grid, (
        f"the {dataset} download and the grid above it disagree: only in the file "
        f"{[c for c in header if c not in grid]}, only on screen "
        f"{[c for c in grid if c not in header]}"
    )


def test_the_pipelines_own_verdict_is_off_the_default_screen_and_one_tick_away():
    """Issue #117's decision, read off what renders rather than off the resolver.

    ``variantalker_naive`` is the pipeline's answer to the question ``Clinical_Summary``
    answers, written into the MAF by ``bin/generate_clinical_summary.py`` and kept by
    ``compute_keep`` on both arms — so the report used to show the same verdict twice,
    once glyphed at the front and once underscore-spelled in the middle, with nothing
    saying they were twins. It is out of the default view now and still in the frame.

    **Both halves are the decision, and neither is safe alone.** Taking it off the screen
    without leaving a way back would drop a column of the user's own file; leaving it in
    the default view is the duplication. So this asserts the pair: absent from the grid the
    report opens with, present once the user asks for everything.

    Driven rather than read, because :func:`~config.columns.resolve_visible_columns`
    returning the right list is not the claim — the claim is that the grid the user sees
    does. A test against the resolver would pass over a call site that had gone on using
    its own list, which is the fault issue #92 found between the grid and the download.

    ``somatic_reference.maf`` carries the column on all 82 rows, so both halves are live
    on this fixture.

    **Both halves read the grid's own list**, the one ``create_data_table`` returned, and
    not a download's header. That distinction is not pedantry: a first draft asserted the
    second half against the *all columns* CSV, which serialises the whole frame and so
    carries the column whatever the grid does — mutating ``create_data_table`` to drop it
    from the *Show all columns* branch left that assertion green, and only the collapse
    test next door noticed. An assertion about what the user can get on screen has to read
    what went on screen.
    """
    as_opened = _run().session_state["grid_columns"]["passed_variants"]
    with_show_all = _run(_SEED_SHOW_ALL).session_state["grid_columns"]["passed_variants"]

    assert "variantalker_naive" not in as_opened, (
        "the report opens with the pipeline's verdict beside the app's own — the "
        "duplication issue #117 removed"
    )
    assert "variantalker_naive" in with_show_all, (
        "'Show all columns' no longer reaches the pipeline's verdict, so leaving it out "
        "of the default view has dropped a column of the user's file instead of moving it"
    )
    # And the app's own verdict is on screen in both states — the half that stays.
    assert "Clinical_Summary" in as_opened
    assert "Clinical_Summary" in with_show_all


def test_show_all_columns_collapses_the_pair_to_one_button():
    """A checkbox must not turn two downloads into one file offered twice.

    Tick *Show all columns* and the grid's set becomes the whole frame, so a download that
    follows the grid delivers what the *all columns* button beside it already delivers.
    Issue #83 spent a ticket removing exactly that — seven buttons carrying five files —
    and its guard could not have caught this one: the grid frame is
    ``prioritize_columns``-reordered and the download frame is not, so the two payloads
    would have differed by column *order* alone and
    ``test_no_two_downloads_hand_over_the_same_bytes`` would have passed over them.

    Only the passed table is ticked, so the failed table's uncollapsed pair is the control
    in the same render: this is a fact about the state, not about the page.
    """
    buttons = _buttons(_SEED_SHOW_ALL)

    assert set(buttons) == {
        "download_passed_all",
        "download_failed_shown",
        "download_failed_all",
        "download_summary",
    }, (
        f"with 'Show all columns' ticked on the passed table the downloads are "
        f"{sorted(buttons)}; the passed pair should have collapsed to one"
    )
    assert "shown" not in buttons["download_passed_all"]["label"].lower(), (
        "the surviving button delivers every column, so it must not describe itself as "
        f"the ones shown: {buttons['download_passed_all']['label']!r}"
    )
