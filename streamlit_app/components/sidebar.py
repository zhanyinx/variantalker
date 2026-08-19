"""The sidebar: what file is open, where you are, and how to open another.

Three things that share one column and one question — *what am I looking at?* The
status block answers it, the navigation moves between the pages that act on it, and
the file chooser is the app's only uploader (issue #64).

Called from ``MAFigate.py`` and from nowhere else. It holds no reference to the table,
the charts or the derived columns: opening a file and reading its results are separate
jobs, and this module owns only the first.

A *Recently opened* list was drawn here too, inside the chooser (issue #158), with the
sentences explaining what it kept and a button to clear it. All of it is deleted, along
with the store behind it: the column's job is to say what is open, and a fold-out list of
names that nothing could click cost more of it than it returned.
"""

import streamlit as st

from utils.main_utils import settings_have_moved


LOAD_STATUS_SLOT = "_load_status_slot"

#: Set by a control drawn *inside* the status block, drained by `render_into_status_slot`.
#:
#: The block is the last thing the app renders and the file chooser is the last thing in the
#: block, so a control above the chooser that called `st.rerun` itself would end the run with
#: the chooser undrawn. That is not merely a discarded frame. `st.rerun` is not a premature
#: stop (`scriptrunner/exec_code.py` sets `premature_stop = False` on `RerunException`), so
#: Streamlit still runs `on_script_finished` and **prunes the state of every widget the run
#: did not register** — the chooser included. The rerun that follows carries no widget states,
#: so the chooser comes back holding `None` while the browser is still holding the file, and
#: the next interaction *anywhere in the app* reports that file against `None`, reads as a
#: fresh choice, and fires `_hand_the_chosen_file_to_the_data_page` — which navigates, before
#: the page whose control the user actually touched has run (issue #277).
#:
#: So a rerun wanted from inside the block is parked here and taken one line later, with the
#: whole column drawn. Same shape as the parameter page's parked confirmations (issue #140)
#: and for a related reason: what a control needs to happen cannot always happen where the
#: control is.
NAV_RERUN_PENDING = "_sidebar_nav_rerun"

#: Every page the app has, and the name the user knows it by.
#:
#: The nav radio built this dict inline as ``label -> page`` and was its only reader, which
#: was enough until the router had to *name* a page it could not open (issue #140). It said
#: ``Error loading page 'data_loading'`` — the internal identifier, on a clinician's screen,
#: for a page the sidebar calls *📊 Load & Analyze Data*. Written this way round because a
#: page is what the app holds and the label is what it shows, and because a router looking
#: up a label from a page is the direction that now has two callers.
PAGE_LABELS = {
    "home": "🏠 Home",
    "parameter_config": "⚙️ Configure Parameters",
    "data_loading": "📊 Load & Analyze Data",
    "help": "❓ Help & Documentation",
}

#: The sidebar's file chooser, and the file it is holding for the data page to load.
#:
#: Two keys rather than one, and the second is not a convenience. The chooser is drawn into
#: the status slot, which is filled *after* the page has run (see `render_into_status_slot`),
#: while the page reads the file near the top of its own body — so the file has to reach the
#: page by a route that does not run through the widget. Widget state belongs to the widget
#: and is only guaranteed while that widget keeps rendering, too: Streamlit prunes the state
#: of every widget a completed run did not register, so a run ending before the slot is
#: filled takes the chosen file with it. The sidebar no longer ends a run that way (issue
#: #283, `NAV_RERUN_PENDING`), which is the defect that made the point; the plain session key
#: holds the file either way, because nothing prunes it.
UPLOAD_CHOOSER_KEY = "sidebar_maf_upload"
PENDING_UPLOAD_KEY = "pending_maf_upload"


def _hand_the_chosen_file_to_the_data_page() -> None:
    """Carry the sidebar's chosen file to the one page that can load it, and go there.

    The chooser is on every page (issue #64), but a MAF is read, validated and filtered by
    the data page, which is where the account of that — the missing columns, the refusals,
    the filter's own warnings — has room to be read. So choosing a file is a navigation as
    well as a load: you land on the report of the file you just opened, which is what the
    upload path inside the page already did (issues #58, #59).

    A callback, because it has to run before the page routes. Streamlit reruns the script
    after it, so this does not ask for a rerun of its own.
    """
    chosen = st.session_state.get(UPLOAD_CHOOSER_KEY)
    if chosen is None:
        # The ✕ inside the uploader **withdraws** the file: it is no longer offered to the
        # page, so a file the page could not read stops being retried and stops re-drawing its
        # error. What it does not do is close the file you are reading — the open MAF, its
        # report and its notes are untouched, because clearing a chooser is not that request.
        st.session_state.pop(PENDING_UPLOAD_KEY, None)
        return

    st.session_state[PENDING_UPLOAD_KEY] = chosen
    st.session_state.current_page = "data_loading"


def _render_file_chooser(ui, a_file_is_open: bool) -> None:
    """Draw the app's one file chooser, sized to the state it is in.

    Prominent when there is nothing open, because it is then the only thing to do. Folded
    into an expander once a file is open, because the sidebar's first job is to say what
    you are looking at — a dropzone above that would push the answer down the column on
    every page.
    """
    box = ui.expander("📂 Open a different file") if a_file_is_open else ui
    # The uploader is the whole of it. A *Recently opened* list was drawn under it (issue
    # #158) and is deleted: nothing in it was clickable, so it named files without offering
    # a way back to any of them, and it made the one state where this chooser is unfolded —
    # no file open — the longest the column ever got.
    box.file_uploader(
        "Choose a MAF file",
        type=["txt", "tsv", "maf"],
        key=UPLOAD_CHOOSER_KEY,
        on_change=_hand_the_chosen_file_to_the_data_page,
        help="MAF, TSV or TXT — the annotated output of the pipeline",
        label_visibility="collapsed" if a_file_is_open else "visible",
    )


def _nav_button(ui, label: str, page: str, to_results: bool = False) -> None:
    """A sidebar button that goes to ``page``.

    A rerun matters, and *where it is taken* matters as much. This button is drawn after the
    page has run, into the status slot, so setting `current_page` alone would leave the new
    page unrendered until the next interaction — and it would leave the nav radio pointing at
    where you used to be, since `create_sidebar_navigation` rebuilds the radio from
    `current_page` only at the top of a run. But the chooser is drawn *below* this button, so
    rerunning from here would end the run with the chooser unregistered and its state pruned
    (issue #277 — see `NAV_RERUN_PENDING`). The request is parked instead, and
    `render_into_status_slot` takes it once the column, chooser included, is complete.

    ``to_results`` asks the data page for its Results section as well as for the page.
    Naming the page is not enough since issue #59: the data page keeps a fixed section
    order and remembers which section you left it on, so a button labelled *back to your
    results* pressed after a visit to Load Data would land on Load Data. Two features that
    are each right alone — a route back, and an order that stops moving — and the label
    is the thing that ends up lying.

    One caller since issue #161, and that argument is why: the nav radio lists every page in
    `PAGE_LABELS`, so a sidebar button naming one of them and nothing else is a second door
    onto a place the column already offers. What earns a button here is a destination the
    radio cannot express — a *section*, as ``to_results`` asks for — which is the rule
    ``tests/test_sidebar_doors.py`` now holds.
    """
    if ui.button(label, use_container_width=True):
        st.session_state.current_page = page
        if to_results:
            # Imported at point of use: `page_modules.data_loading` imports this module,
            # so a module-level import would close the cycle. One writer, one key.
            from page_modules.data_loading import request_results_view

            request_results_view()
        st.session_state[NAV_RERUN_PENDING] = True


def render_load_status(target):
    """Show which file is open, offer the way back to it, and open another.

    Lives above the navigation, so it is on screen on every page. Leaving the results
    for the parameters does not discard the loaded file — but before this block, nothing
    on screen said so and nothing pointed back, so the transition read as losing your
    data (issue #58). The fix is additive: state plus a route, never a guard on the
    transition itself.

    The file chooser sits here too, and this is the app's only one (issue #64): opening a
    sample used to mean navigating to the data page and finding the right section first.
    State and the control that changes it belong in one place — the block already answers
    *which file*, so *choose another* reads as part of the same answer rather than as a
    second feature bolted beside it.

    Three states, and the middle one is the reason this is not a one-liner: a file can be
    open with no results behind it (the filters have not run, or they raised). That state
    must not read as "no file loaded".

    Drawn *after* the page, into the slot `create_sidebar_navigation` reserved above the
    nav — see `render_into_status_slot`. The sidebar is built before the page it sits
    beside, and both load paths read their file from inside that page, so a status block
    rendered in call order would announce "No file open" on the very render that opens it,
    and stay wrong until the next interaction.
    """
    ui = target
    file_name = st.session_state.get("maf_source_name")
    maf_data = st.session_state.get("maf_data")
    filtered_data = st.session_state.get("filtered_data")
    on_data_page = st.session_state.get("current_page") == "data_loading"

    if maf_data is None:
        ui.markdown("📄 **No file open**")
        # No route button here any more. It existed to carry you to the page that held the
        # uploader; the uploader is now the next thing in this column, so a button whose
        # whole job was to take you somewhere else to do this would be the longer way round.
        _render_file_chooser(ui, a_file_is_open=False)
        ui.markdown("---")
        return

    # Both load paths record the name; one that somehow arrived without it still gets a
    # truthful line rather than a blank.
    ui.markdown(f"📄 **{file_name or 'Your MAF file'}**")
    ui.caption(f"{len(maf_data):,} variants in this file")

    if filtered_data is None:
        # This state is usually reached by the filters *failing*, not by their not having
        # run — both load paths filter immediately, and a refusal leaves no results behind
        # while stashing its banner for the data page. So the line says what is true of
        # either route, and does not claim the filters never ran.
        ui.caption("No results for this file yet")
        back_label = "📊 Back to your file"
    else:
        # "current" only holds while the settings behind these results are the settings now
        # in force. The data page warns about exactly this gap; saying "current" across it
        # would have the two halves of one interface disagree.
        stale = settings_have_moved(
            st.session_state.get("filter_params"),
            st.session_state.get("data_params_hash"),
        )
        if stale:
            ui.caption(f"{len(filtered_data):,} passed your last filter run")
            ui.caption("⚠️ Settings have changed since")
        else:
            ui.caption(f"{len(filtered_data):,} passed your current filters")
        back_label = "📊 Back to your results"

    # The one fact the old Session Info panel held that this line did not: which arm the
    # filters are set to.
    sample_type = (st.session_state.get("filter_params") or {}).get("sample_type")
    if sample_type:
        ui.caption(f"{sample_type.title()} analysis")

    if not on_data_page:
        # Only the results route asks for the Results section; the middle state's button
        # goes back to a file that has no report behind it yet.
        _nav_button(ui, back_label, "data_loading", to_results=filtered_data is not None)

    _render_file_chooser(ui, a_file_is_open=True)

    ui.markdown("---")


def render_into_status_slot():
    """Draw the load status into the slot the sidebar reserved, once the page has run.

    Called at the end of the render, so the block reflects a file the page has just
    opened. A no-op until a sidebar has reserved a slot — `create_sidebar_navigation`
    reserves a fresh one on every render, ahead of this call.
    """
    slot = st.session_state.get(LOAD_STATUS_SLOT)
    if slot is None:
        return
    # `st.empty()` holds one element, so `.container()` fills it with the whole block —
    # and the status is drawn exactly once per run, which keeps its button's widget id
    # unique.
    #
    # Drained in a `finally`, and honoured outside it. A request parked by the button above
    # is answered on the run that parked it or not at all: a status block that raises partway
    # has already lost the frame the rerun was for, and a flag left set would fire on the
    # *next* render instead — a page the user did not ask for, arriving one interaction late.
    # So the key is cleared either way and the rerun is taken only on the way out.
    try:
        render_load_status(slot.container())
    finally:
        rerun_wanted = st.session_state.pop(NAV_RERUN_PENDING, False)

    # The last line of the render, and the only place in this module a rerun is taken. Every
    # widget the run draws — this column's chooser last of all — is registered by now, so
    # nothing is pruned when the run ends (issue #277).
    if rerun_wanted:
        st.rerun()


def create_sidebar_navigation():
    """Create sidebar navigation and update session state."""

    st.sidebar.title("🧬 MAFigate Navigation")

    # What is open, and the way back to it, above the nav rather than below it: the state
    # you are in is read before the place you are going. Only the space is claimed here;
    # `render_into_status_slot` fills it after the page has run, for the reason
    # `render_load_status` gives.
    st.session_state[LOAD_STATUS_SLOT] = st.sidebar.empty()

    # Page selection with radio buttons
    page_options = {label: page for page, label in PAGE_LABELS.items()}

    # Find the current page index to set as default
    current_page = st.session_state.get("current_page", "home")
    page_keys = list(page_options.keys())
    page_values = list(page_options.values())

    # Find the index of the current page
    try:
        current_index = page_values.index(current_page)
    except ValueError:
        current_index = 0  # Default to home if page not found

    # The radio uses a widget key, so Streamlit ignores `index` once that key
    # exists and keeps whatever the radio last held. If navigation was triggered
    # programmatically (e.g. a page button changed `current_page` since the last
    # render), drop the stale widget state so the radio is rebuilt cleanly from
    # `index` below. We only delete the key (never assign to it, which corrupts
    # the radio's index-based state); genuine user radio clicks are handled by
    # the reconcile step after the widget.
    #
    # Deleting is enough, and it takes effect in *this* run: `SessionState.__delitem__`
    # drops the key from `_new_widget_state` and `_old_state` alike, so the value the
    # browser sent for the radio goes with it and the registration below has nothing left
    # to prefer over `index`. An `st.rerun()` stood here to force that rebuild, and it was
    # half of issue #277 — it ended the run before the sidebar's file chooser was drawn,
    # which prunes the chooser's state and makes the *next* interaction anywhere in the app
    # fire its `on_change` (see `NAV_RERUN_PENDING`). Nothing between here and the radio
    # reads the key, so there was never a second run's worth of work to do.
    prev_page = st.session_state.get("_nav_prev_page")
    if prev_page is not None and prev_page != current_page:
        st.session_state.pop("navigation_radio", None)
    st.session_state["_nav_prev_page"] = current_page

    selected_page_key = st.sidebar.radio(
        "Navigate to:", page_keys, index=current_index, key="navigation_radio"
    )

    # User selected a different page via the radio -> make it the current page.
    #
    # No rerun, and none is needed: `MAFigate.main` calls this function before
    # `route_to_page`, which reads `current_page` from session state, so the page the user
    # just picked is the page this same run renders. The `st.rerun()` that stood here was
    # the other half of issue #277 — it ended the run before the file chooser below was
    # drawn, pruning its state (see `NAV_RERUN_PENDING`) — and all it bought was rendering
    # the same page a second time.
    selected_page = page_options[selected_page_key]
    if selected_page != current_page:
        st.session_state.current_page = selected_page
        st.session_state["_nav_prev_page"] = selected_page

    # The `📊 Session Info` expander stood here. It is absorbed by `render_load_status`
    # above (issue #58): of the four things it reported, the variant counts and the
    # analysis type are now on screen unfolded, and "Current Page" only ever restated the
    # nav radio sitting directly above it. A collapsed panel was the wrong home for the
    # one fact the user needed at a glance.
    #
    # A `❓ Need Help? Click here` button stood below it, and a rule above that fenced the
    # two apart. Both are deleted (issue #161). It was a survivor rather than an accident —
    # issue #58 absorbed the expander and explicitly kept the button — and the argument for
    # keeping it was that a radio entry reads as a *place* while a button reads as an
    # *offer*, which is the argument #58 made for its own back-to-results button.
    #
    # What separates the two: the back-to-results button asks for something the radio cannot
    # express — a section of the data page, `to_results` below — and this one asked for
    # `help` and nothing else, the same page the radio's fourth entry opens, two elements
    # above it and behind the same ❓ glyph. Measured on the running app: with a file open
    # the column drew eighteen elements, and the only state in which this button behaved
    # differently from that radio entry was the Help page itself, where it offered a trip to
    # the page already on screen — the one thing `render_load_status` guards against for its
    # own button (`if not on_data_page`) and this one did not. A door that is either the same
    # as its neighbour or wrong is not a second invitation.
    #
    # Nothing else routed to `help` expecting it: the nav radio and this button were the two
    # writers of `current_page = "help"` in the app, and the ⋮ menu carries *Get Help* and
    # *Report a bug* for the routes that leave the app (`setup_page_config`). Re-aiming it at
    # a named Help tab was considered and is not available: `st.tabs` cannot be selected from
    # Python. The closest thing to that mechanism was a `help_tab_focus` session key on the
    # Help page, which only printed a hint naming the tab to click and had no writer anywhere
    # in the app; issue #167 deleted it, so nothing on this route survives to be re-aimed.

    # A three-line footer stood here — the version, "🧬 Advanced MAF File Analysis", and
    # "Made with ❤️ using Streamlit" — greyed and centred at the bottom of every page. It is
    # deleted (issue #71): the version now renders in the About dialog alone, the tagline
    # was a truncation of the one `render_header` draws above it on the same screen, and the
    # credit is a fact about how the app was built rather than about the user's file. The
    # sidebar is not spare room either — `render_load_status` put what is open and the way
    # back to it at the top of this same column, and a footer restating the app's name
    # competes with it.
    return selected_page
