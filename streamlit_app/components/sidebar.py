"""The sidebar: what file is open, where you are, and how to open another.

Three things that share one column and one question — *what am I looking at?* The
status block answers it, the navigation moves between the pages that act on it, and
the file chooser is the app's only uploader (issue #64).

Called from ``MAFigate.py`` and from nowhere else. It holds no reference to the table,
the charts or the derived columns: opening a file and reading its results are separate
jobs, and this module owns only the first.

The file *history* is drawn here too, inside the chooser (issue #158). Every sentence the
history shows is written in this module; :mod:`config.file_history` owns only the file it
is stored in.
"""

import streamlit as st

from config.file_history import (
    HISTORY_LIMIT,
    NAME_KEY,
    OPENED_KEY,
    clear_history,
    format_opened,
    read_history,
)
from utils.main_utils import settings_have_moved


LOAD_STATUS_SLOT = "_load_status_slot"

#: Set by the clear button's callback, drained by the render that follows it.
#:
#: A callback rather than a branch on the button's return value, and the ordering is the
#: reason: the list is drawn *above* its own clear control, so a clear handled after the
#: button had returned would leave the just-cleared names on screen until the next
#: interaction. Streamlit runs an ``on_click`` before the rerun, so the render that follows
#: reads an already-empty history. The confirmation cannot be drawn from inside the callback
#: either — nothing is rendering yet — so it is parked here and drawn below, which is the
#: same shape as the parameter page's parked confirmations and for the same reason
#: (issue #140).
HISTORY_CLEARED = "_file_history_cleared"

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
#: so the run that loads a file is a later run than the one that chose it — and on the way
#: there `create_sidebar_navigation` may call `st.rerun` to rebuild the nav radio. Widget
#: state belongs to the widget and is only guaranteed while that widget keeps rendering; a
#: plain session key is not, so the chosen file is copied into one that nothing prunes.
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


def _clear_the_file_history() -> None:
    """Forget the names of the files opened, and park the sentence that says so."""
    clear_history()
    st.session_state[HISTORY_CLEARED] = True


def _render_file_history(box) -> None:
    """Draw the files opened before, under the chooser that opened them.

    Placed here rather than in the sidebar column itself, and that is the decision issue
    #158 turned on. The column then carried four things — the status block, this chooser, the
    nav radio and a Help button — each put there for a good reason and none of them reviewed
    as one column; a fifth *unfolded* item is the question that ticket declined to answer
    here. The Help button has since gone (issue #161), so the count is three, and the reason
    holds unchanged: it was the *unfolded* item this list would have added, not the button,
    that #158 was avoiding. Inside the chooser it costs the column nothing when a file is
    open, because
    the chooser is already folded away in that state, and it sits beside the one control it
    is about: *what have I opened* and *open another* are one question asked twice.

    It follows that the list is on screen unfolded in exactly one state — no file open —
    which is the state where it is most worth reading and where the column is at its
    shortest anyway.

    **Nothing here is clickable.** No entry re-opens a file, on either route. The chooser
    hands the app an in-memory buffer and the app keeps only its name and size, so
    re-opening an uploaded MAF would mean copying patient variant data to disk; the OS
    "Open With" route does have a real path, and still gets no click, because a list that
    behaved differently depending on how the file arrived would be two lists.
    """
    history = read_history()

    # Before the empty check: clearing the list is exactly the act that empties it, so a
    # confirmation drawn only when there is something left to show would never be drawn.
    if st.session_state.pop(HISTORY_CLEARED, False):
        box.caption("The list of files you have opened has been cleared.")

    if not history:
        # Nothing is said in the empty state. A first-time user has no history to explain,
        # and a paragraph about what the app would keep, in a column whose job is to say
        # what is open, would be the app introducing a feature to someone who has not used
        # it yet.
        return

    box.markdown("---")
    box.markdown("**Recently opened**")

    lines = []
    for entry in history:
        when = format_opened(entry.get(OPENED_KEY, ""))
        # A time the stored stamp cannot supply leaves the name standing alone rather than
        # putting "Unknown" beside it: the name is the entry, and the time is what makes the
        # list easy to read.
        lines.append(f"📄 {entry[NAME_KEY]}" + (f" · {when}" if when else ""))
    # One element rather than one per file: ten captions in a sidebar column read as ten
    # separate statements, and this is a single list. Markdown hard breaks keep it tight.
    box.caption("  \n".join(lines))

    # What the app keeps, said once, here — beside the list rather than beside the button,
    # because it is a fact about the list and the reader meets the list first. Two copies of
    # it would be two places to keep true.
    #
    # The cap is named because a name silently dropping off the bottom is the app losing
    # something the user could see; the outliving is named because it is the one thing about
    # this list that is not obvious from looking at it — the file can be gone and the name
    # is still here.
    box.caption(
        f"MAFigate keeps the last {HISTORY_LIMIT} file names on this computer so you can "
        "see what you have looked at. A name stays in this list even after the file itself "
        "is gone. Nothing else about those files is kept, and no name here re-opens one."
    )
    box.button(
        "🗑️ Clear this list",
        on_click=_clear_the_file_history,
        # Named so it stays this button across reruns. The chooser is drawn once per render
        # into a fresh slot, so an unkeyed button here would take its identity from its
        # position in a column whose contents change with the state above it.
        key="clear_file_history",
        help="Forget these file names. Your saved filter settings are not affected.",
        use_container_width=True,
    )


def _render_file_chooser(ui, a_file_is_open: bool) -> None:
    """Draw the app's one file chooser, sized to the state it is in.

    Prominent when there is nothing open, because it is then the only thing to do. Folded
    into an expander once a file is open, because the sidebar's first job is to say what
    you are looking at — a dropzone above that would push the answer down the column on
    every page.
    """
    box = ui.expander("📂 Open a different file") if a_file_is_open else ui
    box.file_uploader(
        "Choose a MAF file",
        type=["txt", "tsv", "maf"],
        key=UPLOAD_CHOOSER_KEY,
        on_change=_hand_the_chosen_file_to_the_data_page,
        help="MAF, TSV or TXT — the annotated output of the pipeline",
        label_visibility="collapsed" if a_file_is_open else "visible",
    )
    _render_file_history(box)


def _nav_button(ui, label: str, page: str, to_results: bool = False) -> None:
    """A sidebar button that goes to ``page``.

    The rerun matters: `create_sidebar_navigation` rebuilds the nav radio from
    `current_page` when it has been changed from outside the radio, so a button that set
    the page without rerunning would leave the radio pointing at where you used to be.

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
        st.rerun()


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
    render_load_status(slot.container())


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
    # render), drop the stale widget state and re-render so the radio is rebuilt
    # cleanly from `index` below. We only delete the key (never assign to it,
    # which corrupts the radio's index-based state); genuine user radio clicks
    # are handled by the reconcile step after the widget.
    prev_page = st.session_state.get("_nav_prev_page")
    if prev_page is not None and prev_page != current_page:
        st.session_state.pop("navigation_radio", None)
        st.session_state["_nav_prev_page"] = current_page
        st.rerun()
    st.session_state["_nav_prev_page"] = current_page

    selected_page_key = st.sidebar.radio(
        "Navigate to:", page_keys, index=current_index, key="navigation_radio"
    )

    # User selected a different page via the radio -> make it the current page.
    selected_page = page_options[selected_page_key]
    if selected_page != current_page:
        st.session_state.current_page = selected_page
        st.session_state["_nav_prev_page"] = selected_page
        st.rerun()

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
