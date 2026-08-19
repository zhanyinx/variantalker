"""
Data loading and analysis page for MAFigate application.
"""

import streamlit as st
import pandas as pd
import os
from datetime import datetime

# Both loaders below are the *pipeline's* reader with different input plumbing (see
# utils/__init__.py). Every load path on this page goes through one of them, so no path
# can hand the filters a differently-typed frame.
from utils import read_maf, read_maf_without_comment_lines
from utils.main_utils import params_hash, settings_have_moved

# One import per component module rather than one per package: `components/` re-exports
# nothing since issue #76, so an import here says which module owns the name.
from components.arm_notice import render_ambiguous_notice, render_mismatch_notice
from components.clinical_summary import add_clinical_summary_column, circle_sources
from components.results_view import create_enhanced_data_table, render_file_facts
from components.sidebar import PENDING_UPLOAD_KEY
from components.variant_table import (
    clear_missing_column_reports,
    drain_missing_column_reports,
    draw_parked_note_confirmation,
)
from filters.numeric_columns import UnreadableNumericColumns
from filters.absent_columns import REQUIRED_INPUTS
from filters.notes import ERROR, INFO, WARNING
from filters.arm_detection import price_other_arm, read_arm_evidence
from filters.attribution import attribute_report
from filters.run_recap import describe_run
from filters.variant_filters import (
    FREQUENCY_COLUMNS,
    MAFIGATE_FILTER,
    MAFIGATE_REASON,
    PASS,
    apply_filters,
    decomposition,
)

# The arm-change route, joined rather than re-invented. Issue #133 made "no path changes
# the arm without saying so" one rule with four call sites; the switch offered by the
# mismatch notice is a fifth, and saying it in this page's own words would be the fifth
# vocabulary that rule exists to prevent. The only import this page takes from another
# page module — `components/sidebar.py` already reaches the other way, so the edge is not
# new, but it is worth naming rather than leaving to be discovered.
from page_modules.parameter_config import adopt_parameters, show_parameter_notice
# This page reads the column categories and nothing else out of `config/`. It used to
# import nineteen more names — every vocabulary list, all four presets, the gene panels,
# the sample types — and used none of them. They were left behind by the filtering work
# that moved the decision into `filters/`, and they made this page look like a consumer of
# constants it has not read in some time.
from config.columns import REQUIRED_COLUMNS
from config.chromosome_spelling import normalise_chromosome_spelling
from config.contaminated_columns import (
    SESSION_KEY as _CONTAMINATED_COLUMNS,
    contaminated_columns,
)
from config.file_history import record_open_file
from config.pipeline_params import pipeline_params

#: Where :func:`_open_the_file_just_read` parks the fact that it re-spelled a chromosome
#: column, for :func:`_show_stashed_banners` to draw. Stashed rather than drawn on the
#: spot for the reason every other load-time message here is: the load paths all filter
#: with ``show_messages=False`` and the page jumps to Results the moment that run
#: finishes, so a message drawn during the load is drawn on the screen the user is being
#: moved off.
_CHROMOSOMES_RESPELLED = "chromosomes_respelled"

#: Where a *silent* filter run parks its account of itself, for the same slot to draw: the
#: notes it produced, or the refusal that stopped it producing any.
#:
#: Named rather than spelled inline at the five sites that touch them, which is the shape
#: issue #155's fix needs: the writer, the drain and the shared load tail all have to agree
#: about the key, and a sixth site spelling it by hand is how one of the three quietly stops
#: clearing it. The two are one account and are always written as a pair — see
#: :func:`_stash_the_runs_account`.
_FILTER_NOTES = "filter_notes"
_FILTER_ERROR = "filter_error"


def validate_required_columns(data):
    """Refuse a MAF the app cannot open at all, and say what the display will be missing.

    What this does **not** do any more, since issue #39: tell the user that a filter has been
    disabled. No filter is ever disabled by a missing column. The vendored filters index every
    input they read, so the app fills an absent one with a neutral value and reports the fill —
    with an escalation where the column also fed pathogenic retention. That report is made by
    the filter, at filter time, from the columns it actually read, and it is the accurate one.

    A second account here would be at best a duplicate and at worst a contradiction: the
    warnings this replaces named ``DP`` as a filter input, which the pipeline has never read —
    it filters on ``t_alt_count + t_ref_count`` — so a user with no ``DP`` column was told a
    depth filter had been switched off when nothing of the kind had happened.

    What survives is the part the filter cannot speak for: the core columns, whose absence
    stops the app before there is anything to filter, and the display columns, which are about
    what the user will be able to *see* rather than what was decided.
    """
    missing_warnings = []

    missing_core = [col for col in REQUIRED_COLUMNS["core"] if col not in data.columns]

    # Show errors for missing core columns
    if missing_core:
        st.error(f"❌ **Missing Core Required Columns**: {', '.join(missing_core)}")
        st.error(
            "⚠️ **Analysis cannot proceed without these columns.** Please check your MAF file format."
        )
        return False

    # The columns the *display* reads: the allele columns, the normal-sample read counts, and
    # DP. The arm's filter inputs are subtracted rather than listed as exceptions, so the claim
    # "filtering is unaffected" is true by construction on any future edit to either list.
    #
    # It said "display columns and **charts** will be blank", and no chart reads any column
    # this note can name (issue #127): the only chart over a sequencing metric is the VAF
    # histogram, which reads `tumor_f` — a filter input, and therefore one of the columns
    # subtracted above. So the sentence promised a consequence for the wrong set and stayed
    # silent about the one case it was true of. It no longer names charts, and it no longer
    # calls these the arm's columns either: `DP` is not arm-specific, and the pipeline emits
    # every one of them on both arms. What is per-arm is which of them this note asks about.
    sample_type = st.session_state.filter_params.get("sample_type", "somatic")
    filter_inputs = set(REQUIRED_INPUTS[sample_type])
    missing_display = [
        col
        for col in REQUIRED_COLUMNS["filtering"] + REQUIRED_COLUMNS[sample_type]
        if col not in data.columns and col not in filter_inputs
    ]
    if missing_display:
        # "will be empty in the table" was the first replacement for the charts claim and was
        # false in the same direction: `resolve_visible_columns` **drops** a column the file
        # does not carry — `[col for col in ordered if col in present]` — so it is not an
        # empty column, it is not a column. The phrasing is the resolver's own, which #93
        # settled for the same reason: the absence belongs to the file, not to MAFigate.
        note = (
            f"ℹ️ **Not in this MAF**: {', '.join(missing_display)} — these are in neither "
            "the table nor anything you download; every other column is unaffected. "
            "Filtering is unaffected too; the filter reports on its own inputs when you "
            "run it."
        )
        # One of them costs more than a column, so the note does not stop there. Notes *and*
        # the user's own annotation values are stored against the variant's gene, position and
        # both alleles, so without the alternate allele two variants at one position cannot be
        # told apart — the app offers neither field on such a row rather than show one user's
        # writing against both, and a user who finds them missing should be told why here
        # rather than have to guess. Both are named because both go: the dialog returns before
        # either.
        if "Tumor_Seq_Allele2" in missing_display:
            note += (
                " Without `Tumor_Seq_Allele2` a variant cannot be identified precisely "
                "enough to write against, so neither the note field nor your own annotation "
                "columns are offered for it."
            )
        missing_warnings.append(note)

    # The filter's own list, not a copy of it (issue #126). This used to read
    # `OPTIONAL_COLUMNS["population_frequency"]`, eight names written out a second time in
    # `config/columns.py`; the two could drift, and this warning is a claim about the very
    # thing `frequency_mask` decides -- so a name in one list and not the other would have
    # the load say "no frequency columns" over a filter that found one, or the reverse.
    available_freq = [col for col in FREQUENCY_COLUMNS if col in data.columns]
    if not available_freq:
        missing_warnings.append(
            "ℹ️ **No Population Frequency Columns Found** - Population frequency filtering will be skipped"
        )

    # Show warnings
    if missing_warnings:
        for warning in missing_warnings:
            st.warning(warning)

    return True


def current_arm() -> str:
    """The arm MAFigate is set to, read the way every other reader reads it.

    ``.get("sample_type", "somatic")``, as in ``MAFigate.initialize_session_state``,
    ``parameter_config.adopt_parameters`` and :func:`validate_required_columns` above: a
    parameter set that states no arm *is* somatic to the app, and a second spelling here
    would let the mismatch notice disagree with the filter about which arm the run was made
    on.
    """
    return (st.session_state.get("filter_params") or {}).get("sample_type", "somatic")


def arm_evidence():
    """What the open file's header says about the arm it was annotated for.

    ``None`` where no file is open — which is not the same as the header being unable to
    place one, and the callers here need to tell those apart.
    """
    maf = st.session_state.get("maf_data")
    if maf is None:
        return None
    return read_arm_evidence(maf.columns)


def file_and_arm_disagree() -> bool:
    """Whether the open file was annotated for an arm other than the selected one.

    Read by the notice *and* by the load paths deciding whether to jump to Results, so the
    screen the user lands on and the sentence they read there cannot come apart. False when
    the header cannot place the file: *cannot tell* is not *disagrees*.
    """
    evidence = arm_evidence()
    return bool(evidence and evidence.disagrees_with(current_arm()))


def _kept_index():
    """The index of the report on screen, or ``None`` when there is no report.

    The distinction matters to the notice. An empty frame is a report that kept nothing and
    can be counted; ``None`` is a run that was refused, where "your current settings keep 0"
    would describe a cut that never happened.
    """
    filtered = st.session_state.get("filtered_data")
    return filtered.index if isinstance(filtered, pd.DataFrame) else None


def _switch_arm_and_refilter(arm: str) -> None:
    """Move to ``arm``, re-cut the open file on it, and open the report.

    Re-seeded from the contract, exactly as the Sample Type control does — and exactly what
    :func:`~filters.arm_detection.price_other_arm` was handed, so the counts the notice
    quoted are the counts this delivers.

    Through ``adopt_parameters`` rather than beside it. This path really can move the arm,
    so it is the kind of writer issue #133's guard exists to catch, and hand-rolling the
    park-and-rerun here would leave the announcement to whoever next edits this function.
    The filtering and the navigation ride in on that helper's ``then`` hook because both
    read the parameters being installed.

    The jump the load withheld is granted here: the mismatch that justified withholding it
    is the thing this button has just resolved, so there is no longer a reason to keep the
    user off their results. A run that produces no report leaves it withheld — jumping to an
    empty Results section would hide the banner explaining why.
    """

    def refilter_and_open_results():
        if apply_filters_to_data(show_messages=False):
            request_results_view()

    adopt_parameters(
        pipeline_params(arm),
        f"Parameters reset to the settings MAFigate opens with for {arm} analysis.",
        then=refilter_and_open_results,
    )


def _show_arm_notice() -> None:
    """Say it when the file and the arm disagree, or when the file cannot be placed at all.

    Re-derived every render rather than stashed, which is a third clock in this slot beside
    the filter's stashed banners and the grids' per-render collection. A mismatch is a
    *standing condition*, not an event: it stays true until the arm changes or another file
    opens, so a one-shot message would vanish on the first navigation while the condition it
    described was still there. Deriving it also means there is nothing to clear and nothing
    that can go stale — the notice cannot outlive the file or the setting it describes.

    Nothing is said for a file carrying **neither** arm's markers, which is every one of the
    12 real unplaceable files — 10 of them, once the 2 the app cannot parse at all are set
    aside, since those never get this far. The filter's own fill warning is already on that
    screen saying the report is not a complete result: escalated where pathogenic retention
    is on, and the ordinary dropout line where the user has switched it off, but never
    silent, because every one of these files is missing a filter input by construction. A
    second box beside it is the wall issue #93 avoided. That
    silence is not a claim the arm is right; this surface only ever speaks on a positive
    finding, and nothing in the app says the arm has been checked.
    """
    evidence = arm_evidence()
    if evidence is None:
        return

    arm = current_arm()

    if evidence.ambiguous:
        render_ambiguous_notice(evidence, arm)
        return

    if not evidence.disagrees_with(arm):
        return

    kept = _kept_index()
    comparison = (
        price_other_arm(
            st.session_state.maf_data, pipeline_params(evidence.detected), kept
        )
        if kept is not None
        else None
    )

    render_mismatch_notice(
        evidence=evidence,
        current_arm=arm,
        comparison=comparison,
        # The switch re-seeds, so it discards — but only where there is something to
        # discard. On the reported session the arm came from the cache and the settings
        # were untouched, so nothing is lost and the notice says nothing about losing it.
        customised=st.session_state.get("filter_params") != pipeline_params(arm),
        on_switch=_switch_arm_and_refilter,
    )


def _open_the_file_just_read(data, name: str) -> bool:
    """Make a frame that has just been read the open file, and filter it.

    The tail every load path shares, in one copy: the incoming name is recorded, the frame
    becomes ``maf_data``, everything the *previous* file's filter run left behind is dropped,
    the columns are checked, and the filters run silently. Returns whether a report was
    produced — the callers use it to decide whether to open Results, and a refused MAF has no
    report to open.

    What is **not** dropped is what the user wrote. Notes and annotation columns survive a
    change of file, because a note is about the variant rather than about the sample in front
    of you (issue #67); issue #64 cleared them here and that step is gone.

    There are three doors into this — the sidebar's chooser, and the OS "Open With" handler's
    two readers — and before issue #64 each carried its own copy of the sequence. What must
    not depend on which door the file came through is the answer to *what is still true of the
    last file*: ``failed_data`` and ``data_params_hash`` both outlived the file they describe,
    so a new MAF whose filter run then failed left the session holding the old file's rejects
    and the old run's stamp — the stamp that decides whether the sidebar calls your results
    *current*.

    The name is recorded here rather than on arrival (which is what issue #58 did) because
    here the read has already succeeded: a file the app could not open never becomes the file
    the sidebar names.
    """
    st.session_state.maf_source_name = name

    # Before the frame becomes `maf_data`, because everything downstream reads the
    # chromosome as a *settled* value: the four sites that render it, and the note key,
    # which builds a variant's identity out of it and would otherwise file the same
    # variant under two names depending on how its file spelled it (issue #149). The
    # pipeline settles this one step upstream, while annotating, and this is the same
    # act at the same point in the sequence — see `config/chromosome_spelling.py`.
    #
    # Whether to *say so* is a separate question with a later answer, below: the sentence
    # promises something about a download, and a MAF that is about to be refused has no
    # report to download.
    respelled = normalise_chromosome_spelling(data)

    # And immediately after it, because this verdict *is* a comparison against `Chromosome`
    # (issue #194): a predictor column ANNOVAR did not supply holds the chromosome rather than
    # a score, and asking whether it does before the spelling is settled would ask it of one of
    # two spellings. Taken once per file rather than per rendered cell — the measurement behind
    # that choice is in `config/contaminated_columns.py`.
    #
    # Assigned rather than stashed-and-drained: the panel reads it on every variant a user
    # opens, so it has to survive the render that first uses it. Replaced here, in the tail
    # that replaces `maf_data`, so it describes the frame beside it and never the last one.
    st.session_state[_CONTAMINATED_COLUMNS] = contaminated_columns(data)

    st.session_state.maf_data = data
    st.session_state.filtered_data = None
    st.session_state.failed_data = pd.DataFrame()
    st.session_state.data_params_hash = None
    # In the same tail for the same reason: which sources the overview draws is a fact
    # about the file that has just been replaced, so it must not outlive it either. A new
    # MAF whose filter run then fails would otherwise leave the key describing the old one.
    st.session_state.overview_sources = []
    # And the recap of the run that produced the old report (issue #137), which is a
    # statement about a cut that no longer has a table beside it. Its `note` names columns
    # *that* file did not carry, so leaving it standing would attach one MAF's missing
    # columns to the next one's report.
    st.session_state.filter_run = None
    # Its other half (issue #147), for the same reason plus one it does not share: the recap
    # is wrong about the replaced file's *columns*, and this is wrong about its **row
    # counts** — "33,535 variants read" over a MAF holding 200 is the kind of claim a reader
    # checks the report against rather than doubts.
    st.session_state.filter_attribution = None
    # And this file's chromosome notice, for the plainest form of the same reason: it is a
    # statement about the file being replaced. Cleared *before* it can be set below, so a
    # bare-spelled MAF whose notice was never drawn — no render of this page between the
    # two loads — cannot have it drawn against the prefixed file that followed it (issue
    # #149).
    st.session_state.pop(_CHROMOSOMES_RESPELLED, None)

    # The last file's filter run, on the same rule and for the same reason. Issue #149 put
    # the line above here and called it "the only stash in this tail that a later render
    # rather than a later run consumes" — true of the keys the tail *writes*, and misleading
    # about the tail as a whole, which is the mistake issue #155 exists to undo. These two
    # are stashed one step later, by the filter run this tail ends with, and drawn by the
    # same render that draws the notice above; so a MAF whose account was never drawn used
    # to leave it standing for the next file to be blamed for.
    #
    # Not made redundant by the writer clearing them (`_stash_the_runs_account`), and this
    # is the case that says so: a MAF **refused** by `validate_required_columns` below
    # returns before the filter runs at all, so no run happens to replace anything. Without
    # these two pops the previous file's notes are drawn beside the refusal of the next.
    st.session_state.pop(_FILTER_NOTES, None)
    st.session_state.pop(_FILTER_ERROR, None)

    if not validate_required_columns(data):
        st.session_state.maf_data = None
        return False

    # Here, and this is the same judgement the name above is recorded by, one step further on
    # (issue #158). The sidebar's history means *files you have looked at*, so it is written
    # once the file has been read **and** accepted: a MAF the app could not read never gets
    # this far, and one refused for a missing column has just been unloaded on the line
    # above. The account of either is the error on the page, which has room for the reason;
    # a name in a recents list has none.
    #
    # Not moved below the filter run either. A file whose filters refuse is still open — the
    # sidebar names it, and its middle state exists precisely for that — so it is a file the
    # user has looked at.
    record_open_file(name)

    # After the refusal, not before it. The sentence ends "*and your download will spell
    # them that way too*", and a refused MAF produces no report and offers no download —
    # so stashed any earlier it is drawn beside an error, promising something about a file
    # the app has just declined to open. It is deliberately *not* also behind the filter
    # run below: a file that loads and then fails to filter has still been re-spelled, and
    # the user can still reach the columns it changed.
    if respelled:
        st.session_state[_CHROMOSOMES_RESPELLED] = True

    with st.spinner("Applying default filters..."):
        return apply_filters_to_data(show_messages=False)


def _auto_load_file_from_path(file_path):
    """Load a MAF file directly from a file path (e.g., macOS 'Open With')."""
    name = os.path.basename(file_path)

    try:
        with st.spinner(f"Loading {name}..."):
            data = read_maf(file_path)

    except Exception as error:
        # Fallback: same reader, reached with the comment lines removed so it locates the
        # header itself instead of inheriting the pipeline's count of them.
        try:
            with open(file_path, "r", encoding="utf-8") as handle:
                content = handle.read()
            data = read_maf_without_comment_lines(content)

        except Exception as fallback_error:
            st.error(
                f"❌ Error loading file: {str(error)}. Fallback: {str(fallback_error)}"
            )
            return False

    return _open_the_file_just_read(data, name)


# The two sections of this page, keyed by identity rather than by their labels: the
# Results label carries a live count, and a keyed widget whose stored selection stops
# matching its own options loses that selection.
#
# The first section is the *filter run*, not the load, since issue #65. Issue #64 moved the
# uploader to the sidebar and issue #59 moved the load's banners up to the page, which
# between them left a section called "Load Data" with no loading on it — what remained was
# how the report was reached and the button that re-reaches it.
_FILTER_RUN = "filter_run"
_RESULTS = "results"
_SECTION_KEY = "data_page_section"
_SECTION_PREV = "data_page_section_prev"
_JUMP_TO_RESULTS = "jump_to_results"
_LAST_UPLOAD = "last_upload_token"


def show_data_loading_page():
    """Display the data loading and analysis page."""

    # `Load your MAF file and apply filtering parameters.` stood here. Both halves had
    # stopped being true of this page: the file is chosen in the sidebar (issue #64) and the
    # parameters are set on their own page. What a corrected line would say — which sections
    # this page offers — the section control directly beneath it already says.
    st.title("📊 Data Loading & Analysis")

    # Initialize session state if needed
    if "maf_data" not in st.session_state:
        st.session_state.maf_data = None
    if "filtered_data" not in st.session_state:
        st.session_state.filtered_data = None
    if (
        "filter_params" not in st.session_state
        or st.session_state.filter_params is None
    ):
        # Deep copy, for the reason MAFigate.py gives at the same assignment: a shallow
        # one shares the nested keep-lists with the module constant.
        st.session_state.filter_params = pipeline_params("somatic")

    # Auto-load file if opened via "Open With" (macOS / Windows). This runs before the
    # section control exists, so the jump it asks for is honoured in this very run — and
    # `auto_load_file` is popped, so it happens once per opened file.
    if "auto_load_file" in st.session_state:
        file_path = st.session_state.pop("auto_load_file")
        if os.path.isfile(file_path):
            st.info(f"Opening: **{os.path.basename(file_path)}**")
            # Withheld on a mismatch, like the sidebar's load below: a report cut on the
            # wrong arm is not one to drop the user into.
            if _auto_load_file_from_path(file_path) and not file_and_arm_disagree():
                request_results_view()

    # A file chosen in the sidebar is loaded here, in the page body, for the reason
    # `_load_pending_upload` gives — and *before* the jump is honoured, which is what keeps
    # everything the load has to say on the screen the user ends up looking at. A jump asked
    # for here is applied by the pop below, in this same run, so no `st.rerun` discards the
    # frame the missing-column warnings were drawn into. The upload block this replaces sat
    # after the section control, could only ask for the jump *next* run, and so lost those
    # warnings on every successful first load.
    _load_pending_upload()

    # A slot for the filter run's banners, claimed above the control and filled at the end
    # of this function. They belong to the page rather than to a section — a filter run
    # stashes them from whichever section is on screen, and the page may be about to move
    # the user off it — but they have to be *written* after that section has run, because
    # that is when the run that stashes them happens. The slot buys the second without
    # giving up the first.
    banners = st.empty()

    # The same slot now also carries what the *grids* found missing (issue #93), and the
    # collection is emptied here — before the section draws into it — rather than at the
    # drain. Emptying at the drain would leave it holding this render's sentences for the
    # next render to redraw, on a page whose whole point is that the sentence is said once
    # per render; emptying here means what the slot shows was resolved by the section
    # immediately above it, on this frame, or is not shown at all.
    clear_missing_column_reports()

    if st.session_state.maf_data is None:
        # No switch when there is nothing to switch between (issue #65). Both sections say
        # "load a file first" in this state, so offering the choice is the page arguing with
        # itself; and since issue #64 the way in is a sentence pointing at the sidebar, which
        # is not a section of anything.
        #
        # A pending jump is deliberately *not* consumed here. Nothing can request one while
        # no file is open — every caller of `request_results_view` has either just loaded a
        # file or is the sidebar's route back, which is gated on there being a report — but
        # popping it on a run that draws no control would discard it in the one place it
        # could never be honoured, and a jump that goes missing is harder to see than one
        # that waits.
        _show_the_way_in()
    else:
        # Honoured *before* the control is instantiated: a widget's key cannot be written
        # once the widget exists in the same run. Popping it is what makes the jump fire
        # once — a jump that survived into the next rerender would pull the user back to
        # Results every time they left it.
        if st.session_state.pop(_JUMP_TO_RESULTS, False):
            st.session_state[_SECTION_KEY] = _RESULTS

        section = _section_control()

        if section == _RESULTS:
            show_results_section()
        else:
            show_filter_run_section()

    # A jump requested by the *section* that just rendered — the Re-apply button — cannot
    # take effect until the next run, because the control above it already exists. This is
    # that run. The flag is only ever set by a genuinely new filter run, so this cannot loop.
    #
    # Before the banners, and this is load-bearing: `st.rerun` throws away everything this
    # run drew, so consuming the stash first would consume it into a discarded frame. Left
    # stashed, it is rendered by the run the user actually sees.
    if st.session_state.get(_JUMP_TO_RESULTS):
        st.rerun()

    with banners.container():
        _show_stashed_banners()

    # After the rerun check above, for the same reason the banners are: a note saved in the
    # variant dialog parks its confirmation and reruns, and draining it before a rerun would
    # drain it into a frame nobody sees (issue #140). Not inside the banner slot — a toast
    # is not drawn into the page, and this sentence belongs over the grid it just changed.
    draw_parked_note_confirmation()


def _section_label(section: str) -> str:
    """The user-facing label for a section, with the report's size when there is one.

    The count survives the fixed order. It was introduced to draw the eye, which the jump
    now does, but it says something the eye cannot get anywhere else on this page: how big
    the report is, readable from the filter-run side without leaving it.
    """
    if section == _FILTER_RUN:
        return "🔍 Filter run"

    filtered = st.session_state.get("filtered_data")
    if isinstance(filtered, pd.DataFrame):
        return f"📊 Results ({len(filtered)})"
    return "📊 Results"


def _keep_a_section_selected() -> None:
    """Undo a deselection.

    ``st.segmented_control`` in single-selection mode lets a click on the *active* segment
    clear it, which would leave this page with no section chosen. Restoring the previous
    value makes that click a no-op. It has to happen in a callback: that runs before the
    script, so the widget does not yet exist and its key is still writable.
    """
    if st.session_state.get(_SECTION_KEY) is None:
        st.session_state[_SECTION_KEY] = st.session_state.get(_SECTION_PREV, _FILTER_RUN)


def _section_control() -> str:
    """Render the Filter-run/Results switch and return the chosen section.

    Drawn only when a file is open (issue #65) — see `show_data_loading_page`.

    A segmented control rather than ``st.tabs``, for two reasons that are really one.
    Streamlit's tabs cannot be selected from code — ``st.tabs`` takes neither a key nor a
    default — so the only way to surface Results used to be to *reorder* the tabs, which
    moved the page's primary control out from under the user. And tab labels are plain
    text, so enlarging them would have meant a CSS block pinned to Streamlit's internal
    class names; a segmented control is a button group at full width already.
    """
    # Seeded rather than passed as `default=`, which Streamlit warns about once the key
    # also lives in session state — and it always does here, because the jump writes it.
    if _SECTION_KEY not in st.session_state:
        st.session_state[_SECTION_KEY] = _FILTER_RUN

    section = st.segmented_control(
        "Section",
        options=[_FILTER_RUN, _RESULTS],
        format_func=_section_label,
        key=_SECTION_KEY,
        on_change=_keep_a_section_selected,
        label_visibility="collapsed",
        width="stretch",
    )
    section = section or _FILTER_RUN
    st.session_state[_SECTION_PREV] = section
    return section


def request_results_view() -> None:
    """Ask the data page to show Results, once, at the next opportunity to honour it.

    Public because the sidebar's route back to your results calls it too (issue #58): the
    section order is fixed and remembered, so naming the page is no longer enough to land
    on the report.
    """
    st.session_state[_JUMP_TO_RESULTS] = True


def _load_pending_upload() -> None:
    """Load the file the sidebar's chooser is holding, once, and open its report.

    Called from the page body rather than from the load section, and that is the whole
    reason this function exists rather than the upload block it replaces. The chooser is in
    the sidebar now (issue #64), so a file can arrive while the user is reading Results, or
    the help page, or the parameters — and a load wired into the load section would sit
    unread until they happened to visit it.

    Guarded by the same ``(name, size)`` token the jump used to carry alone. The chooser
    hands the same file back on every rerender, so "a file is here" and "a file has arrived"
    are different events, and this now runs on every render of every section: without the
    guard a 500 MB MAF would be re-read and re-filtered on every click anywhere in the app.
    Stamping it only after a successful filter run keeps the old retry behaviour — a MAF
    refused for missing a column the *other* arm requires loads on the next attempt, once
    the arm is switched, with no need to choose the file again.

    A file that cannot be **read** is withdrawn rather than retried, and that asymmetry is the
    point: the arm-switch retry above changes the verdict, while re-reading the same bytes with
    the same two readers cannot. Left on offer it would re-raise and re-draw its banner on
    every render of the page, with the chooser's own ✕ the only way out — and this is a session
    where nothing else went wrong, so there must be a way out that is not a page refresh.

    A read that fails also leaves the *open* file alone: nothing has been written to session
    state at that point, so the MAF the user was reading, its report and its notes all survive
    a mis-picked file. Only a file that reads and then fails validation replaces them.
    """
    upload = st.session_state.get(PENDING_UPLOAD_KEY)
    if upload is None:
        return

    token = (upload.name, getattr(upload, "size", None))
    if st.session_state.get(_LAST_UPLOAD) == token:
        return

    try:
        # read_maf takes the upload buffer directly; spilling it to a temporary file for the
        # pipeline's path-taking reader is its job, not this page's.
        with st.spinner("Loading MAF file..."):
            data = read_maf(upload)
    except Exception as error:
        # Fallback: same reader, reached with the comment lines removed so it locates the
        # header itself instead of inheriting the pipeline's count of them.
        try:
            data = read_maf_without_comment_lines(upload.getvalue())
        except Exception as fallback_error:
            st.session_state.pop(PENDING_UPLOAD_KEY, None)
            st.error(
                f"❌ Error loading file: {str(error)}. Fallback error: {str(fallback_error)}"
            )
            return

    if _open_the_file_just_read(data, upload.name):
        st.session_state[_LAST_UPLOAD] = token
        # The jump is withheld when the file and the arm disagree (issue #135), so the user
        # meets the mismatch notice on the Filter run section rather than arriving at a
        # screen of variants cut on the wrong arm. The *token* is stamped either way: it
        # guards the re-read, not the navigation, and leaving it unstamped would re-read and
        # re-filter the file on every render for exactly the files this notice is about.
        if not file_and_arm_disagree():
            request_results_view()


def _show_the_way_in():
    """What the page says when no file is open.

    Not a *section* since issue #65 — the page draws this instead of the switch, because
    both sections say "load a file first" in this state. No uploader of its own either,
    since issue #64: the app has exactly one file chooser and it lives in the sidebar,
    where it is reachable from every page. Two uploaders for one job would be two widget
    states, two accounts of which file is current, and a user who can hold a different file
    in each.
    """
    st.subheader("📁 Open a MAF file")
    st.markdown(
        "Choose your file in the sidebar, on the left. If the sidebar is hidden, the "
        "**»** arrow at the top left brings it back."
    )
    st.caption(
        "MAFigate reads the annotated `.maf`, `.tsv` or `.txt` output of the pipeline."
    )


def _results_are_stale() -> bool:
    """Whether the report on screen was produced by settings that have since changed.

    The load section's warning and the recap above the results are the same question asked
    from two places, so they ask it once — and through
    :func:`utils.main_utils.settings_have_moved`, which the sidebar asks it through too.
    """
    return settings_have_moved(
        st.session_state.get("filter_params"),
        st.session_state.get("data_params_hash"),
    )


def _reapply_from_results() -> None:
    """Re-run the filters for the recap's stale branch, and stay where the user is.

    Silent, like every other call site: a warning rendered from inside the results section
    would be drawn below the report rather than in the page-level slot that collects them,
    and the stash is what puts it above the section switch on the next render.

    No ``request_results_view`` here, unlike the load section's button — the press came
    from Results, so there is nowhere to travel to. The rerun is what redraws the report
    and the recap that now describes it.
    """
    apply_filters_to_data(show_messages=False)
    st.rerun()


def show_filter_run_section():
    """How the open file's report was reached, and the button that reaches it again.

    Everything else that stood here is gone (issue #65), and the reason each went is the
    same one: it was drawn beside its own twin.

    - *"To open a different file, use the sidebar"* pointed at a chooser **on the same
      screen**, labelled ``📂 Open a different file``. It was written when the uploader had
      just moved (issue #64) and read as a note to users who remembered the old place; it
      rendered on every visit, forever, next to the thing it described.
    - *Configure Parameters* was a third door to a page the sidebar's nav radio already
      offers by name, on this same screen — and the parameter page's own *Go to Data
      Analysis* button is the return leg, so the round trip survives its deletion.
    - *Data Overview*'s **Total Variants** repeats the sidebar's ``N variants in this
      file``, again on the same screen. Its other two tiles, Samples and Genes, said
      something no other surface did, so they moved to the Results summary rather than
      being dropped — and they now read ``N/A`` rather than ``0`` where the column they
      count is absent, which is what the file actually justifies.
    - *Data Preview* showed ten raw rows of a file whose every row is on the Results grid.
    - *Column Information* tabulated all 427 columns of a reference MAF with their pandas
      dtypes — ``float64``, ``object`` — which is implementation vocabulary reaching the
      interface, and the Results summary's Column Analysis already answers the question a
      user actually has, which is whether the columns that matter are present.

    What is left is one job: this is the report's provenance and the control that refreshes
    it, which is why the section is no longer called *Load Data*.
    """

    # Warn if filter parameters changed since last filter run
    if _results_are_stale():
        st.warning(
            "⚠️ Filter parameters have changed since the last filter run. "
            "Click **Re-apply Filters** to update results."
        )

    # What the last filter run decided, decomposed into the four cells.
    #
    # This replaces an echo of the parameter dict. The echo was broken by the
    # parameter contract in four separate ways — it read a catch-all sentinel that no
    # longer means anything, the redundant germline VAF key, and the classification
    # list as though it were an include list — and an echo of the inputs never told
    # the user the thing they need, which is how the report was reached. The cells
    # come off the masks that made the verdict, so they cannot disagree with the rows
    # on screen.
    #
    # The heading is drawn only when there is a report to have reached. On a refused MAF
    # this section is most of what is on screen, which is how the pair came to read *How
    # this report was reached: No filters have been applied yet.* — a heading promising a
    # report directly above a line saying there is none.
    if st.session_state.get("filtered_data") is not None:
        st.markdown("**How this report was reached:**")
    st.info(_decomposition_summary())

    # Sized to its text rather than stretched across the page. It was half of a two-button
    # row; alone at `use_container_width=True` it would be a full-width bar, which is more
    # emphasis than it had before this section was thinned, not less.
    if st.button("🔍 Re-apply Filters", type="primary"):
        # Silent, like every other call site, because this one now ends on the
        # Results section too: a warning rendered inline here would be rendered
        # onto the screen the user is being moved off. Stashing it puts it above
        # the section control on arrival. The success line is not stashed with it —
        # the refreshed report, and the count beside "Results", say the same thing
        # in the place the user is looking.
        if apply_filters_to_data(show_messages=False):
            request_results_view()


def _show_stashed_banners() -> None:
    """Render what the last silent filter run had to say, and what this render resolved.

    Page-level rather than inside the first section, which is where it used to sit. Every
    load path filters with ``show_messages=False`` and stashes instead, and the page now
    jumps to Results as soon as that run finishes — so a banner rendered by the Load
    section is a banner shown on the one screen the user has just been moved off.

    Four things arrive here, on three clocks. The filter's refusal and its run notes
    are stashed by a *filter run* and read by whichever later render draws this page. The
    chromosome notice is stashed by the *load*, one step earlier (issue #149), and so
    survives a file that then failed to filter — which is right: the re-spelling happened
    either way and the export carries it either way. The missing-column sentence is
    collected by *this* render's grids, a few frames down, and is gone again at the end of
    it — same seam, same slot, because the user cannot tell one clock from the other and
    all of them are statements about the file they just opened.

    That last clause is the one with a rule under it. Everything drawn here is a statement
    about **the open file's latest run**, so none of it may outlive either: the three stashed
    keys are cleared by the shared load tail on a change of file, and the two written by the
    filter run are replaced — cleared when there is nothing to say — on every silent run
    (issue #155). A drain alone cannot keep that promise, because it only runs when this page
    renders, and the whole reason these are stashed is that the load does not wait for one.
    """
    # What a wholesale parameter replacement parked for us (issue #133). Drawn here as well
    # as on the parameter page because the mismatch notice's switch button is a fifth site
    # for that event, and it fires from *this* page — a confirmation drawn only where the
    # parameter page draws it would reach the user on some later visit, announcing a switch
    # they made minutes ago somewhere else. Only one page renders per run, so the two
    # drains cannot both fire.
    show_parameter_notice()

    # Before the filter's account, because it explains it. A germline MAF on the somatic arm
    # draws `❌ PATHOGENIC RETENTION DEGRADED — CancerVar column not found`, which is true of
    # the somatic arm and blames the file; read second, it is a consequence of the sentence
    # above rather than a defect in the user's data. What that warning should *say* once a
    # mismatch is known is issue #136's, deliberately not settled here.
    _show_arm_notice()

    # Display any filter errors that occurred during silent application. The banner is
    # shown verbatim: `apply_filters_to_data` stores the framing along with the message,
    # because it is the only place that knows which of the two it hit. This used to
    # prefix whatever it found with "Error applying filters", which meant a refused MAF
    # (#38) was framed as a filter error on the silent path and as a refusal on the
    # loud one — the same file, two different accounts of what went wrong.
    if _FILTER_ERROR in st.session_state:
        st.error(st.session_state[_FILTER_ERROR])
        del st.session_state[_FILTER_ERROR]

    # Same seam, the other outcome: the file *was* filtered, but a filter-input column
    # was missing and had to be filled (#39). The load paths all filter silently, so
    # without this the user meets a report that is not a complete result — and possibly
    # much smaller than it should be — with nothing on screen saying so.
    if _FILTER_NOTES in st.session_state:
        _show_filter_notes(st.session_state[_FILTER_NOTES])
        del st.session_state[_FILTER_NOTES]

    # Under both of the above, because it is the least consequential thing this slot says:
    # nothing about the verdict moves on it — the pipeline's filters never read the
    # chromosome — and it is the one message here that reports something the app *did*
    # rather than something the file lacks. It is said at all because it reaches the
    # download: `filtered_data` descends from `maf_data`, so a report taken from a
    # bare-spelled file carries the new spelling out of the app with it (issue #149).
    if st.session_state.pop(_CHROMOSOMES_RESPELLED, False):
        st.info(
            "ℹ️ This file numbers its chromosomes 1, 2, 3. MAFigate reads them as "
            "chr1, chr2, chr3, and your download will spell them that way too."
        )

    # Last: what the grids found this render (issue #93). Under the filter's
    # account deliberately, and **not merged with it**, which was the dev's call. The two
    # answer different questions about the same file — the filter's is *your verdict may
    # be wrong*, this is *your table and your export are narrower* — and their column sets
    # only partly overlap: a MAF missing `am_class` and `cosmic` alone gets nothing from
    # the filter, those being output columns it never reads, and this sentence is then the
    # only thing the app says. Adjacent and slightly overlapping beats one message that
    # has to be true of both a filter run and a render.
    for report in drain_missing_column_reports():
        st.warning(report)


def show_results_section():
    """Display the results section."""

    st.subheader("📊 Analysis Results")

    if st.session_state.maf_data is None:
        st.info("📁 Please load a MAF file first.")
        return

    if st.session_state.filtered_data is None:
        # Not "please apply filters first", which is the claim issue #58 corrected in the
        # sidebar and left standing here: both load paths filter the moment a file opens, so
        # this state is normally reached by the filters *running and failing*, not by their
        # never having run. What went wrong is already on screen, in the page-level banner
        # above the section switch.
        st.info("ℹ️ No results for this file yet — see the message above the sections.")

        # The file still loaded, so what is in it is still worth saying. Without this the
        # counts rehoused by issue #65 would live only on the summary tab, which this branch
        # returns before ever reaching — so a refused MAF, the one case where the old Data
        # Overview still had something to show, would show them nowhere.
        render_file_facts(st.session_state.maf_data)
        return

    # Get failed data if available, otherwise create empty DataFrame
    failed_data = getattr(st.session_state, "failed_data", pd.DataFrame())

    # Defensive check - ensure data integrity
    try:
        # Verify data is still valid
        if not isinstance(st.session_state.filtered_data, pd.DataFrame):
            st.error("❌ Filtered data is corrupted. Please reapply filters.")
            return

        if not isinstance(failed_data, pd.DataFrame):
            failed_data = pd.DataFrame()

        # Guarded on its own, inside the integrity check that precedes it and separately
        # from the tables that follow. A report that will not assemble must say so without
        # taking the variant tables down with it — which is what the old block below the
        # tabs achieved by keeping a separate `except` around each of its buttons — and
        # must not be attempted at all over a frame already known to be corrupt.
        try:
            summary_report = create_summary_report()
        except Exception as e:
            summary_report = None
            st.error(f"❌ Error preparing the run report: {str(e)}")

        # Display enhanced results with passed/failed variants. Every download the
        # results section offers is drawn inside this call, under the tab it belongs to:
        # the two CSVs beneath each variant table, the run report beneath the summary of
        # the run. The `💾 Download Results` block that used to follow is gone (issue
        # #83) — its two variant buttons served the same bytes as the tabs' `all columns`
        # buttons, differing only in the name the file landed under.
        create_enhanced_data_table(
            st.session_state.filtered_data,
            failed_data,
            summary_report=summary_report,
            # What the filters were, taken when they ran (issue #137). Staleness is
            # computed here, from the same stamp and the same recipe the sidebar and the
            # load section use — three surfaces saying one thing, not three answers.
            recap=st.session_state.get("filter_run"),
            stale=_results_are_stale(),
            on_reapply=_reapply_from_results,
            # And what that cut did (issue #147), taken on the same clock as the recap it
            # is drawn beneath.
            attribution=st.session_state.get("filter_attribution"),
        )

    except Exception as e:
        st.error(f"❌ Error displaying results: {str(e)}")
        st.info("Please try reapplying the filters to refresh the data.")
        return


def apply_filters_to_data(show_messages=True) -> bool:
    """Apply filters to the MAF data. Returns whether a report was produced.

    The return value is what the callers use to decide whether to open Results. A refusal
    or a crash leaves no report to open, and jumping to an empty Results section would
    hide the very banner explaining why.
    """

    if st.session_state.maf_data is None:
        if show_messages:
            st.error("❌ No MAF data loaded.")
        return False

    try:
        with st.spinner("Applying filters..."):
            labelled, diagnostics = apply_filters(
                st.session_state.maf_data, st.session_state.filter_params
            )

        # The entry point returns every row labelled; the split is the caller's to make.
        # Deriving it here rather than having the filter return two frames is what keeps
        # one place where the decision is made — and the verdict column travels with both
        # frames, so a row on screen carries the reason it is there.
        verdict = labelled[MAFIGATE_FILTER] == PASS
        passed_data = labelled[verdict]
        failed_data = labelled[~verdict]
        notes = list(diagnostics.notes)

        # Which sources the Pathogenicity Overview can speak for is a fact about this arm
        # and this file, so it is settled here, once, and carried — `variant_table` draws
        # the key from `overview_sources` rather than working it out again beside the
        # grid. Recomputing there would read the *current* arm, which a visit to the
        # parameter page can change without clearing this report, and the key would then
        # name a different set of sources from the one the circles were drawn for.
        sources = circle_sources(
            st.session_state.filter_params.get("sample_type", "somatic"),
            list(st.session_state.maf_data.columns),
            skip_civic=bool(st.session_state.filter_params.get("skip_civic", False)),
        )
        st.session_state.overview_sources = sources

        # Enrich with clinical summary once, then store
        st.session_state.filtered_data = add_clinical_summary_column(
            passed_data, sources
        )
        st.session_state.failed_data = add_clinical_summary_column(failed_data, sources)

        # Stamp the params hash so we can detect stale data later. The sidebar's load
        # status compares against this same stamp, via the same recipe (issue #58).
        st.session_state.data_params_hash = params_hash(st.session_state.filter_params)

        # And, on the same clock, what those parameters *say* — the recap the results view
        # draws above its tabs (issue #137). Taken here rather than read live at render
        # time for the reason the stamp beside it exists: the parameters can move away from
        # the report on screen, and a recap that followed them would describe a cut the
        # table beside it was not made with. `circle_sources` above is settled here for the
        # same reason (issue #95). It also carries the diagnostics' account of any filled
        # column, which exists only at this moment — the warning below is drawn once and is
        # gone from the screen after that render.
        # The file's columns as well as the parameters, because two clauses can be dropped
        # whole by a column the MAF does not carry — CIViC, which the pipeline guards for
        # itself, and the frequency layer, which skips when no frequency column is present.
        # Neither fills anything, so neither reaches `filled_columns`, and a recap built
        # from parameters alone would report both as restrictions that applied.
        st.session_state.filter_run = describe_run(
            st.session_state.filter_params,
            list(st.session_state.maf_data.columns),
            diagnostics,
        )

        # And what the cut *did*, on the same clock and for the same reason (issue #147):
        # how many of the excluded variants each setting is responsible for. Taken here so it
        # cannot be re-derived at render time against parameters that have since moved and
        # end up contradicting the recap it is drawn beneath.
        #
        # Six extra whole-frame filter runs, measured at 27-68ms per preset on real MAFs
        # against a single run's 6-16ms — so it roughly quadruples the work behind a spinner
        # that already says "Applying filters...", and reaches ~600ms on the largest MAF in
        # the measured corpus (33,535 rows). Paid eagerly rather than on first open of the
        # expander, because a lazy read would pair whatever frame is loaded *then* with this
        # run's parameters, and `filter_run` outlives a file the way `filtered_data` does not.
        #
        # It re-derives the rejected mask from the same parameters it isolates against rather
        # than taking `labelled` from above, which is the one thing that makes it unable to
        # pair one cut's counts with another cut's settings — the defect #137 found in the
        # downloaded report.
        st.session_state.filter_attribution = attribute_report(
            st.session_state.maf_data, st.session_state.filter_params
        )

        # Three levels, not one. A missing column that feeds pathogenic retention means rows
        # are *absent* from this report — the one thing on this page the user cannot see for
        # themselves and cannot recover by re-reading their parameters — so it is an error
        # box, above the success line, not the eighth yellow box in a column of them. And a
        # note about a run that did exactly what was asked is blue, so that the yellow ones
        # are worth reading; issue #151 measured that the slot's only output on a complete
        # MAF was a yellow box about a filter that had worked. Which level each note is, is
        # the answer of whichever module wrote its sentence, not a substring test here.
        #
        # And they are shown on *both* paths. Every load call site filters with
        # `show_messages=False`, because it filters before the results tab exists — so a
        # warning shown only when `show_messages` is true is a warning the user never sees on
        # the load that produced it, and only meets later if they happen to press Apply. That
        # is exactly the seam #38 found for the refusal banner; the answer is the same one:
        # stash it and let the page render it on the next rerender.
        if show_messages:
            _show_filter_notes(notes)
        else:
            _stash_the_runs_account(notes)

        if show_messages:
            st.success(
                f"✅ Filters applied successfully! {len(passed_data)} variants passed, {len(failed_data)} variants failed."
            )

        return True

    except UnreadableNumericColumns as refusal:
        # A refusal, not a crash — so it is not framed as "Error applying filters", which
        # would send the user off to change a filter. Nothing about the parameters is wrong:
        # this MAF cannot be filtered by the pipeline either, and the message already names
        # every offending column and value.
        _report_filter_failure(f"🛑 {refusal}", show_messages)
        return False

    except Exception as e:
        _report_filter_failure(f"❌ Error applying filters: {str(e)}", show_messages)
        return False


#: How each of :data:`filters.notes.LEVELS` is drawn. A mapping and not a chain of ``if``s so
#: that a level added to ``filters.notes`` and forgotten here is a *missing key*, which
#: ``tests/test_filter_notes.py`` asserts against ``LEVELS`` directly rather than trusting this
#: to be read alongside it.
#:
#: **Streamlit function *names*, not the functions.** Holding ``st.error`` itself would bind the
#: real callable at import time, and the suite's stubs work by replacing this module's ``st``
#: global — eight sites do it — so a frozen reference silently bypasses the stub and calls the
#: live Streamlit. The old single-predicate renderer resolved ``st.warning`` through the global
#: on every call and so was interceptable by accident; this keeps that property on purpose.
_NOTE_RENDERERS = {
    ERROR: "error",
    WARNING: "warning",
    INFO: "info",
}


def _show_filter_notes(notes) -> None:
    """Draw each of the filter's notes at the level it carries.

    Three levels since issue #151, and no test on the message text at all — the module that
    writes the sentence states what kind of note it is, and this reads
    :attr:`filters.notes.Note.level`. What that replaced was a single predicate,
    ``filters.absent_columns.is_escalated``, which put everything non-escalated in a yellow
    warning box: on both reference MAFs, on their own arm, under all four shipped presets, the
    only box that slot ever drew was a yellow one reporting a population-frequency filter that
    had worked perfectly, so yellow meant nothing by the time a real warning arrived.

    Deliberately **not** grouped by severity. The order is the filter's own — fills, then the
    gene list, then frequency — because a user reading the boxes is reading an account of one
    run in the order it happened, and two of these notes explain each other when adjacent
    (an absent ``Hugo_Symbol`` beside the requested symbols it could not look for).
    """
    for note in notes:
        getattr(st, _NOTE_RENDERERS[note.level])(note.text)


def _stash_the_runs_account(notes) -> None:
    """Leave this run's account in the slot, and nothing of any earlier run's.

    A **replacement**, not an addition, and that is issue #155's fix. This was
    ``elif notes:``, so a run with nothing to say neither wrote the stash nor cleared it: the
    previous run's notes stayed standing and were drawn by whichever render came next,
    against a file — or a cut — they were not about. Silence is a statement too. It says the
    filter had nothing to report, and the slot has to be able to say it.

    **Zero notes is ordinary, not a corner**, which is why the hole mattered: it is better than
    a third of real filter runs, and it is what the app's own opening state does, because
    ``pipeline_params`` ships ``max_freq_population = 1.0`` and ``apply_filters`` skips the whole
    frequency layer above 1.0. That set is what the app opens with, what *Reset to Defaults*
    installs and what the mismatch notice's arm switch installs. Issue #151's "4 of 4" measured
    the four shipped presets on the two references — the one configuration that is reliably
    loud — which is what made this read as the rare branch.

    The corpus, the run counts and the rates are recorded once, in
    ``tests/test_stale_banners.py``'s module docstring, beside the tests that hold them. Not
    repeated here: two copies of one measurement are two things to keep true, and this one is
    already load-bearing in a test.

    The refusal is cleared alongside, because the two keys are one account: a run either
    produced a report to describe or was stopped before it could, never both.
    :func:`_report_filter_failure` clears this one from the other side for the same reason.

    One outcome of a silent run reaches neither of them, and that is deliberate rather than an
    omission: :func:`apply_filters_to_data` returns early with no file open, touching no key.
    Every one of its four callers is gated on an open file in the same run — this one runs after
    ``validate_required_columns`` succeeded, both *Re-apply* buttons are drawn only where
    ``maf_data`` is not ``None``, and the arm switch is instantiated by :func:`_show_arm_notice`,
    whose :func:`arm_evidence` returns ``None`` without a file. Clearing there would be a line
    no test could make fail.
    """
    if notes:
        st.session_state[_FILTER_NOTES] = list(notes)
    else:
        st.session_state.pop(_FILTER_NOTES, None)
    st.session_state.pop(_FILTER_ERROR, None)


def _report_filter_failure(banner: str, show_messages: bool) -> None:
    """Show ``banner`` now, or leave it for the page to show on the next rerender.

    ``show_messages=False`` is the silent load path, which filters before the results tab
    exists and so has nowhere to put a banner yet. The **framed** text is what gets stored,
    not the bare message: the framing is a judgement about what went wrong, this is the only
    place that has made that judgement, and storing the message alone is what let the two
    paths describe one file two different ways.

    The notes go with it (issue #155). A refusal is raised by ``apply_filters`` before this
    run reaches :func:`_stash_the_runs_account`, so without clearing them here the slot would
    draw an *earlier* run's notes underneath this run's refusal — an account of a report
    beside the news that no report was produced.
    """
    if show_messages:
        st.error(banner)
    else:
        st.session_state[_FILTER_ERROR] = banner
        st.session_state.pop(_FILTER_NOTES, None)


def _decomposition_summary() -> str:
    """The four cells of the last filter run, as one line of prose.

    Re-derived from the reason column carried by the two stored frames rather than from a
    new session-state entry, so there is nothing that can go stale against what is on
    screen: if a frame is displayed, its own rows are what is being counted.

    Worth naming for a clinician reading it: the rescue-only cell is the number of rows in
    the report that did **not** meet the quality thresholds and are present solely because
    they are pathogenic. Nothing in the app reported that before — on the reference it is
    19.2% of the germline report against 0.7% of the somatic.
    """
    columns = [
        frame[MAFIGATE_REASON]
        for frame in (
            st.session_state.get("filtered_data"),
            st.session_state.get("failed_data"),
        )
        if frame is not None and MAFIGATE_REASON in frame.columns
    ]
    if not columns:
        return "No filters have been applied yet."

    cells = decomposition(pd.concat(columns))
    total = sum(cells.values())
    if not total:
        return "No variants to report."

    return (
        f"{total} variants: "
        f"{cells['criteria_only']} passed on the filter criteria, "
        f"{cells['both']} met the criteria and are also pathogenic, "
        f"{cells['rescue_only']} kept only by pathogenic retention "
        f"(these did not meet the criteria), "
        f"{cells['rejected']} rejected."
    )


def _parameter_provenance() -> str:
    """The parameters that produced this cut, listed verbatim.

    A downloaded report that does not say what it was filtered with cannot be checked
    against anything later, which matters more here than on screen — the file outlives the
    session. So the export keeps parameter provenance even though issue #28 deletes all
    eight of the app's on-screen parameter echoes.

    Deliberately a **verbatim dump of the dict**, not a rendered interpretation. The
    function this replaces tried to describe the parameters in prose and was broken by the
    parameter contract in four separate ways — it read a catch-all sentinel that no longer
    means anything, the superseded germline VAF key, and the classification list as though
    it were an include list. A dump cannot misdescribe what it prints, and it survived
    #36's move of these keys to the pipeline's contract without being touched.

    **The run's parameters, not the session's** (issue #137). This read
    ``filter_params`` live, so a user who changed a setting and downloaded the report
    without re-filtering got a file whose counts came from one cut and whose parameter list
    came from another — the defect that ticket exists to fix on screen, in the one artifact
    that outlives the session and so cannot be checked against anything later. It reads the
    snapshot the filter run took instead, and falls back to the live dict only where there
    is no run to describe.
    """
    run = st.session_state.get("filter_run")
    params = (run.params if run else None) or st.session_state.get("filter_params") or {}
    if not params:
        return "  (none recorded)"
    return "\n".join(f"  {key}: {params[key]!r}" for key in sorted(params))


def create_summary_report():
    """Create a summary report of the filtering results."""

    if st.session_state.filtered_data is None:
        return "No filtered data available."

    original_count = len(st.session_state.maf_data)
    filtered_count = len(st.session_state.filtered_data)
    failed_count = len(getattr(st.session_state, "failed_data", pd.DataFrame()))

    report = f"""
MAF Analysis Summary Report
Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

Input Data:
- Original variants: {original_count}
- Passed filters: {filtered_count}
- Failed filters: {failed_count}
- Pass rate: {(filtered_count/original_count*100):.1f}%

How this report was reached:
{_decomposition_summary()}

Filter Parameters:
{_parameter_provenance()}

Top Genes (in passed variants):
"""

    if (
        "Hugo_Symbol" in st.session_state.filtered_data.columns
        and len(st.session_state.filtered_data) > 0
    ):
        top_genes = (
            st.session_state.filtered_data["Hugo_Symbol"].value_counts().head(10)
        )
        for gene, count in top_genes.items():
            report += f"- {gene}: {count} variants\n"

    # Add failed variants analysis if available
    failed_data = getattr(st.session_state, "failed_data", pd.DataFrame())
    if len(failed_data) > 0 and "Hugo_Symbol" in failed_data.columns:
        report += "\nTop Genes (in failed variants):\n"
        top_failed_genes = failed_data["Hugo_Symbol"].value_counts().head(10)
        for gene, count in top_failed_genes.items():
            report += f"- {gene}: {count} variants\n"

    return report
