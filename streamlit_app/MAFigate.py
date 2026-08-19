"""
MAFigate: Advanced MAF File Analysis and Variant Filtering

MAFigate is a comprehensive Streamlit application designed for clinical researchers
and bioinformatics specialists to analyze Mutation Annotation Format (MAF) files.
The application provides sophisticated filtering capabilities for both somatic and
germline variants, over the clinical annotations already written into the file.

"Clinical database integration" is what this said until issue #219, and the same claim had already
been deleted from the About dialog further down this file — the comment at ``menu_items`` records
why — on a ground that applies here unchanged: MAFigate **integrates with nothing**. It queries no
service and reads the annotations your file already carries. So the wording below is that
correction's, not a new one. *Database* is false of most of the sources the bullet named, too: map
#199 established that RENOVO is an in-silico predictor this pipeline ran (issue #201, where passing
*computed here* is exactly the leg a third-party database fails), that ``InterVar`` and
``CancerVar`` are guideline verdicts computed here (issue #187), and that AlphaMissense is a
predictor on its own scale (issue #203). ClinVar and CIViC were the only two of the five the word
fitted. The line below names the sources and claims no category, which is the one framing no member
can falsify.

Key Features:
- Multi-format file support (MAF, TSV, TXT)
- Clinical and research parameter presets
- Gene-based filtering with predefined panel gene sets (MSK-IMPACT, FoundationOne, COSMIC, OncoKB)
- Reads the clinical annotations already in your file (ClinVar, CIViC, ESCAT, CancerVar,
  InterVar, RENOVO, AlphaMissense) — MAFigate queries nothing
- Population frequency filtering (gnomAD, ExAC, ESP, 1000 Genomes)
- Real-time filtering with detailed statistics
- Export capabilities in multiple formats

Architecture:
- Modular design with separation of concerns
- Page-based navigation system
- Reusable UI components
- Comprehensive error handling
- Performance optimization for large datasets

Author: Development Team
License: See repository LICENSE file

The version is deliberately not written here. It was, twice — in the opening line and in a
``Version:`` field at the end — and because a docstring renders nowhere, issue #71's sweep
over the app's *surfaces* could not see either copy while both went stale by a major number.
Neither is quoted back here on purpose: a note explaining a stale literal by reprinting it
leaves the literal in the file. ``config/constants.py``'s ``APP_VERSION`` is the one place
the number lives, the About dialog is the one place it renders, and ``build/version.py`` is
how it reaches an installer (issue #260).
"""

import streamlit as st
import sys
import os
import traceback

# Add the current directory to Python path for imports
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

# Import application components
from config.build_identity import build_identity
from config.constants import APP_NAME, APP_TAGLINE, APP_VERSION
from config.pipeline_params import pipeline_params
from components.sidebar import (
    PAGE_LABELS,
    create_sidebar_navigation,
    render_into_status_slot,
)
from components.variant_table import written_notes
from page_modules.home import show_home_page
from page_modules.parameter_config import show_parameter_config_page

# The cache is restored and the discard banner rendered before any page exists, which is
# why these three come from a module of their own rather than out of the parameter page.
from page_modules.param_store import (
    get_cache_info,
    load_parameters_from_cache,
    show_discarded_cache_banner,
)
from page_modules.data_loading import show_data_loading_page
from page_modules.help import show_help_page

#: Where the two clickable menu items land, and the only two URLs in this app a user can reach
#: by clicking (issue #243). ``GET_HELP_URL`` is the app's **own** README, one directory into
#: the public repository, rather than the repository root: the root README is the *Nextflow
#: pipeline's*, and a clinician who clicks for help would land on DSL2 badges and a channel
#: diagram with MAFigate one line down the contents. Named here rather than written inline
#: because the first of them does not fit on the line it is used on.
GET_HELP_URL = "https://github.com/zhanyinx/variantalker/blob/main/streamlit_app/README.md"
REPORT_A_BUG_URL = "https://github.com/zhanyinx/variantalker/issues"


def initialize_session_state() -> None:
    """
    Initialize Streamlit session state variables.

    Sets up default values for all session state variables used throughout
    the application, including navigation state, filter parameters, data storage,
    and parameter tracking for change detection.

    Session State Variables:
        current_page (str): Current active page identifier
        filter_params (Dict[str, Any]): Current filtering parameters
        maf_data (Optional[pd.DataFrame]): Loaded MAF data
        maf_source_name (Optional[str]): File name behind maf_data, shown in the sidebar
        filtered_data (Optional[pd.DataFrame]): Results after filtering (with Clinical_Summary)
        data_params_hash (Optional[str]): Hash of filter_params at last filter run (stale data detection)
        overview_sources (list): Sources the Pathogenicity Overview drew, settled at filter time
    """

    # Page navigation
    if "current_page" not in st.session_state:
        st.session_state.current_page = "home"

    # Filter parameters - try loading from cache first, fallback to defaults
    if "filter_params" not in st.session_state:
        # Only a cache carrying this version's parameter-format stamp is restored. An
        # older one is set aside by `show_discarded_cache_banner` in the header, which
        # runs before any page: the app cannot tell a pre-parity cache from a current one
        # by looking at it, and restoring one silently is what made "the app opens at
        # parity" untrue for every user who had opened it before (issue #40).
        cached_params = load_parameters_from_cache()
        if cached_params:
            st.session_state.filter_params = cached_params
            # Store cache loading info for display
            if "cache_loaded" not in st.session_state:
                st.session_state.cache_loaded = True
                cache_info = get_cache_info()
                if cache_info:
                    st.session_state.cache_timestamp = cache_info["timestamp"]
        else:
            # `pipeline_params` returns a fresh deep copy, and both halves matter: the
            # pages edit this dict in place and the keep-lists are nested, so anything
            # shallower would share those lists with the contract itself and let one
            # session redefine the app's default for the rest of the process.
            st.session_state.filter_params = pipeline_params("somatic")
            st.session_state.cache_loaded = False

    # Data storage
    if "maf_data" not in st.session_state:
        st.session_state.maf_data = None

    if "filtered_data" not in st.session_state:
        st.session_state.filtered_data = None

    # Name of the file behind `maf_data`, for the sidebar's status block (issue #58).
    # Only meaningful while `maf_data` is set — see `_auto_load_file_from_path`.
    if "maf_source_name" not in st.session_state:
        st.session_state.maf_source_name = None

    # Hash of filter_params when filters were last applied (stale data detection)
    if "data_params_hash" not in st.session_state:
        st.session_state.data_params_hash = None

    # Which annotation sources the Pathogenicity Overview drew circles for, settled at
    # filter time from the arm and the file (issue #95). Empty until a report exists, and
    # empty is a real answer: a MAF carrying none of its arm's sources gets no overview.
    if "overview_sources" not in st.session_state:
        st.session_state.overview_sources = []

    # Auto-load file from environment (e.g., macOS "Open With")
    if "auto_load_checked" not in st.session_state:
        st.session_state.auto_load_checked = True
        open_file = os.environ.get("MAFIGATE_OPEN_FILE", "")
        if open_file and os.path.isfile(open_file):
            st.session_state.auto_load_file = open_file
            st.session_state.current_page = "data_loading"


def setup_page_config() -> None:
    """
    Configure Streamlit page settings and metadata.

    Sets up the page title, icon, layout, sidebar state, and menu items
    for the MAFigate application. This function should be called once
    at the beginning of the application execution.

    Configuration:
        - Wide layout for better data visualization
        - Expanded sidebar for easy navigation
        - Custom menu items with help and bug report links
        - DNA emoji (🧬) as page icon

    This dialog is the **only** place the version number renders (issue #71). It was in
    six: here, the browser tab, the header line under the title, the sidebar footer, the
    Help page's subtitle, and — found by review, after the first four had been dealt with —
    a worked citation example in the Help page's FAQ, where it was written out by hand
    rather than read from the constant. A release number is a fact a bug report needs and a
    clinician does not, and the About menu is where a bug reporter already looks for it.

    Since issue #263 it names the **channel** and the **build** as well, on the same line and
    for the same reader. A version alone cannot identify a build: the same ``APP_VERSION``
    reaches a user as a .dmg, as a Windows .exe or as a clone, and with no update check ever
    (#229) this dialog is the only thing that tells the maintainer which. The two extra facts
    are added *here*, in the one surface already justified, rather than to a second one —
    ``tests/test_app_identity.py`` runs the tab, the header, the sidebar and the Help page to
    prove all three of them render nowhere else, and that sweep now covers these two strings
    as well as the version. See ``config/build_identity.py`` for where they come from.
    """
    st.set_page_config(
        page_title=APP_NAME,
        page_icon="🧬",
        layout="wide",
        initial_sidebar_state="expanded",
        menu_items={
            "Get Help": GET_HELP_URL,
            "Report a bug": REPORT_A_BUG_URL,
            # The five-bullet Features list that stood here is gone rather than corrected.
            # Two bullets were untrue of the app a user can actually drive — it advertised
            # "Export capabilities in multiple formats" when every reachable variant
            # download is CSV, and "Clinical database integration" when the app integrates
            # with nothing and reads the annotations already in your file — and a third
            # named four of the six gene panels. Home's expanders say all of this
            # accurately already, so a corrected list here would be a fourth copy to keep
            # true. About points at them instead.
            #
            # One line, three facts, in the order a bug report needs them: what the app is,
            # which release, and which build of that release. `build_identity()` is called
            # here rather than read at import so that a stamp written between two runs is
            # picked up by the next launch rather than by the next reinstall.
            "About": f"""
            {APP_NAME} v{APP_VERSION} · {build_identity()}

            {APP_TAGLINE}

            What {APP_NAME} reads, and what each setting does, is on the Home page and
            under Help & Documentation.
            """,
        },
    )


def render_header() -> None:
    """
    Render the application header with title and navigation.

    Displays the MAFigate title and its one-sentence description at the top of each page
    for consistent branding and user orientation.

    The sentence is :data:`config.constants.APP_TAGLINE` rather than a literal, because
    this is not the only surface that says what the app is — the About dialog reads the
    same constant. On Home it renders directly above the two buttons with nothing between,
    so it is the only description a user landing on the app is given, which is why it says
    what the app does to a file instead of leading with a version number (issue #71).
    """
    st.title(f"🧬 {APP_NAME}")
    st.markdown(APP_TAGLINE)

    # A parameter cache written before the format was stamped is set aside here, before
    # any page renders, and the user is told what it held (issue #40).
    show_discarded_cache_banner()

    # Show cache loading notification if parameters were loaded from cache
    if (
        st.session_state.get("cache_loaded", False)
        and "cache_timestamp" in st.session_state
    ):
        try:
            from datetime import datetime

            cache_time = datetime.fromisoformat(
                st.session_state.cache_timestamp
            ).strftime("%Y-%m-%d %H:%M")
            st.success(
                f"💾 Restored your previous parameter settings from {cache_time}"
            )
            # Only show this once per session
            st.session_state.cache_loaded = False
        except Exception:
            # Silently ignore datetime parsing errors
            pass


def route_to_page() -> None:
    """
    Route to the appropriate page based on current selection.

    Uses the current_page session state variable to determine which
    page component to render. Each page is implemented as a separate
    module in the page_modules directory.

    A page that fails is **said so where the user is standing**, and the app stays there
    (issue #140). It used to draw ``st.error``, set the page to Home and rerun — and
    ``st.rerun`` discards the frame it had just drawn into, so a user whose page failed was
    moved to Home with no explanation at all, which is an error they cannot report. There
    is nothing to rerun *for*: the nav radio is drawn before this function, so leaving it
    alone is what keeps the sidebar and the page body naming the same page, and Home is not
    a destination the user asked for.

    The unknown-page branch raises rather than reporting itself. It is unreachable —
    every writer of ``current_page`` writes one of these four literals — and a second,
    unreachable account of a broken page is a second sentence to keep true. Nothing escapes
    this function, which is why the docstring no longer claims a ``ValueError`` does.
    """
    current_page = st.session_state.get("current_page", "home")

    try:
        if current_page == "home":
            show_home_page()
        elif current_page == "parameter_config":
            show_parameter_config_page()
        elif current_page == "data_loading":
            show_data_loading_page()
        elif current_page == "help":
            show_help_page()
        else:
            raise ValueError(f"Unknown page: {current_page}")
    except Exception:
        # Named as the sidebar names it. `Error loading page 'data_loading'` put an
        # internal identifier on screen for a page the user knows as `📊 Load & Analyze
        # Data`; the label comes from the nav's own table so the two cannot drift.
        label = PAGE_LABELS.get(current_page, current_page)
        st.error(f"❌ **{label} could not be opened.**")
        st.write("You can pick another page in the sidebar.")
        # Folded, and labelled for the job it does. The traceback is the only part of this
        # that can get a defect fixed, and it is the only part written in Python rather
        # than in the user's terms — so it is a click away rather than absent. Refreshing
        # is deliberately not suggested: every note in the session lives in memory alone
        # (issue #67), so that advice would quietly throw the user's writing away.
        with st.expander("Details, for a bug report"):
            st.code(traceback.format_exc())


def main() -> None:
    """
    Main application entry point.

    Orchestrates the entire MAFigate application by:
    1. Setting up page configuration
    2. Initializing session state
    3. Rendering the header and navigation
    4. Routing to the appropriate page
    5. Filling the sidebar's load-status slot, now that the page has run

    This function serves as the central coordinator for all application
    components and handles the overall application flow.
    """
    try:
        # Setup page configuration
        setup_page_config()

        # Initialize session state
        initialize_session_state()

        # Create sidebar navigation (this updates session state)
        create_sidebar_navigation()

        # Render header
        render_header()

        # Add separator
        st.markdown("---")

        # Route to appropriate page
        try:
            route_to_page()
        finally:
            # Fill the sidebar's reserved status slot last: both load paths open their file
            # from inside the page above, so this is the first point in the render where
            # what is loaded is actually known (issue #58). In a `finally`, because the
            # render that most needs to say your file is still open — and offer the way
            # back to it — is the one where the page itself has just failed.
            render_into_status_slot()

    except Exception:
        # One box rather than two, and no `str(e)` line. The traceback folded away below ends
        # with exactly that sentence, so drawing it here as well was the same fact twice —
        # and the copy on the clinician's screen was the one written in Python.
        #
        # `ran into a problem` rather than `could not start`, because this handler covers two
        # unalike moments. Since issue #140 `route_to_page` reports its own failures and lets
        # nothing escape, so what reaches here is `setup_page_config`, `initialize_session_state`,
        # `create_sidebar_navigation` and `render_header` — none of which has drawn a page yet —
        # or `render_into_status_slot`, which runs in the `finally` **after** the page, and can
        # fail with the whole app already on screen. `could not start` is false in that second
        # case, and it is the case where the user can simply carry on.
        st.error("❌ **MAFigate ran into a problem.**")
        _say_what_a_refresh_costs()
        # `contact support` is gone with nothing in its place. There is no support address in
        # this repo to name — issue #55's line, no promise the app does not keep — and the ⋮
        # menu already carries *Report a bug*, set by `setup_page_config`, which runs first and
        # has therefore succeeded in every failure that reaches here bar its own.
        with st.expander("Details, for a bug report"):
            st.code(traceback.format_exc())


def _say_what_a_refresh_costs() -> None:
    """The advice of last resort, and the one act that makes writing survive it.

    A refresh is sometimes the only way out of here — if `create_sidebar_navigation` is what
    failed there is no navigation left to use — and it works *precisely because* it throws the
    session away. That is also what destroys the user's writing: a note and every custom
    annotation value live in `st.session_state` and nowhere else (issue #67), which the variant
    dialog says in as many words. So the old line, *Please refresh the page or contact support*,
    told the user to discard their writing at the moment they were most likely to do it, from a
    surface that never asked whether there was any.

    It asks now. `written_notes` is empty exactly when there is nothing a refresh can destroy
    that the user cannot get back — the report re-derives from a file still on their disk, the
    writing does not — so it is both the gate and the payload, and a user with nothing written
    is not alarmed about it.

    Deliberately not guarded by a `try`. The rescue reads two dicts out of session state and
    consults no frame, which is why issue #144 chose it over serialising the report; wrapping it
    would add a second account of failure to a handler whose whole job is to give the first one.
    """
    written = written_notes()
    if written.empty:
        st.write("Refreshing the page starts MAFigate over.")
        return

    # Offered before the sentence that explains it: this is the act, and the sentence is why.
    st.download_button(
        label="📥 Download what you have written (CSV)",
        data=written.to_csv(index=False),
        file_name="mafigate_notes.csv",
        mime="text/csv",
        key="download_writing_after_error",
    )
    st.write(
        "Refreshing the page starts MAFigate over. What you have written lives in this "
        "session only — download it first."
    )


if __name__ == "__main__":
    main()
