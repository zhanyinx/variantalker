"""The parameter configuration page: six tabs of controls over ``filter_params``.

What this module is *not* responsible for, and used to be:

* where those parameters live between sessions — that is
  :mod:`page_modules.param_store`, which the app's boot sequence needs before any page
  exists and so had no business being reachable only through a page;
* what values a control may offer and how a stored value is cleaned against them — that is
  :mod:`config.vocabularies`, which is Streamlit-free and readable by the parity harness.

What is left is the widgets: the tabs, their layout, and the session-state keys that carry
each arm's gene selection. It is long because Streamlit pages are long, not because it
holds more than one kind of thing.
"""

import copy
import json
from datetime import datetime

import streamlit as st
import yaml

from config.gene_panels import GENE_LISTS, PREDEFINED_GENE_SETS
from config.param_labels import label_of
from config.param_migration import (
    PARAM_SCHEMA_VERSION,
    SCHEMA_VERSION_KEY,
    migrate_params,
)
from config.pipeline_params import pipeline_params
from config.presets import (
    CLINICAL_GERMLINE_PARAMS,
    CLINICAL_SOMATIC_PARAMS,
    SOFT_GERMLINE_PARAMS,
    SOFT_SOMATIC_PARAMS,
)
from config.vocabularies import (
    CANCERVAR_OPTIONS,
    CIVIC_OPTIONS,
    CLINVAR_OTHER_ASSERTION_TERMS,
    CLINVAR_PATHOGENICITY_TERMS,
    ESCAT_OPTIONS,
    ESCAT_STRONGEST,
    GUIDELINE_SOURCES,
    INTERVAR_OPTIONS,
    RENOVO_OPTIONS,
    SAMPLE_TYPES,
    VARIANT_CLASSIFICATIONS,
    filter_terms,
    validate_multiselect_params,
)
from filters.gene_lists import (
    GENE_LIST_EXTENSIONS,
    GeneSelection,
    panel_symbols,
    parse_gene_list,
)
from page_modules.param_store import (
    PARAMS_CACHE_FILE,
    SUPERSEDED_CACHE_SUFFIX,
    clear_parameters_cache,
    export_text,
    get_cache_info,
    load_parameters_from_cache,
    save_parameters_to_cache,
    show_discarded_cache_banner,
)


def gene_panel_state_key(sample_type):
    """Where the gene panel dropdown's choice lives — session state, per arm.

    Deliberately **not** ``filter_params``. The panel is UI state that *resolves* to a gene
    list; the list is the parameter. Storing the panel name as a parameter — which this page
    used to do, on both arms — wrote it into every saved parameter file and every cache
    entry, where it read as a filter setting that nothing applied.
    ``config/pipeline_params.py`` states the rule, and
    ``tests/test_gene_lists.py::test_the_panel_choice_is_never_a_filter_parameter`` checks
    this file against it.

    The widget's own ``key`` *is* the storage, and that is the point rather than an
    incidental detail. Streamlit persists a keyed widget's value across reruns by itself, so
    mirroring it into a second dict and feeding that back as ``index=`` gives the widget two
    sources of truth — and when they disagree ``index`` wins. Measured on the first draft of
    this page: selecting a panel and then selecting ``All`` again left the dropdown reading
    ``All`` while the mirror still said ``MSK_IMPACT``, so 505 genes went on filtering a
    report the user had just widened. One source of truth, no ``index``.
    """
    return f"param_{sample_type}_gene_set"


def gene_text_state_key(sample_type):
    """Where the typed gene list lives — the text box's own key, per arm.

    Named for the same reason :func:`gene_panel_state_key` is: it is read in one place and
    written in another, and an inlined f-string that has to match across those two places is
    a typo away from a box whose contents are silently never read.
    """
    return f"gene_text_{sample_type}"


def gene_upload_state_key(sample_type):
    """Where the uploaded gene file lives — the uploader's own key, per arm."""
    return f"gene_file_{sample_type}"


def _pathogenic_retention_on(params):
    """Whether pathogenic retention is on, read the way the filter reads it.

    The page must not answer this question differently from
    ``filters.variant_filters._skip_pathogenic``, and the two spellings are exactly where
    it could: the contract says ``skip_pathogenic`` and the SOFT and CLINICAL presets
    still say ``keep_pathogenic``. The filter prefers ``skip_pathogenic`` whenever a dict
    carries it — so a checkbox that read and wrote the *other* key had no effect at all
    once the contract became the default, which is precisely what happened. Measured on
    somatic_reference.maf: ticked and unticked both gave the same 18 rescued rows, where
    the real setting moves the union 20 → 19 and empties both rescue cells.
    """
    if "skip_pathogenic" in params:
        return not params["skip_pathogenic"]
    return bool(params.get("keep_pathogenic", True))


#: Where a wholesale parameter replacement parks what to say about itself, for the render
#: that follows the rerun. See :func:`adopt_parameters`.
PARAM_NOTICE_KEY = "last_parameter_notice"


def switched_to(arm):
    """How the app says the arm moved — one definition, both routes.

    The Sample Type control and :func:`show_parameter_notice` are the same event
    reached differently, so they open with the same words and then diverge: that control
    really does reset to the settings MAFigate opens with, and a preset brings its own.
    Sharing the opening is the point; sharing the sentence would make one of them false.
    """
    return f"🔄 Switched to {arm}:"


def park_confirmation(confirmation, arm=None):
    """Park one sentence for the render after the rerun, because this one is discarded.

    ``st.rerun`` raises immediately and the run's output is thrown away, so a message
    written on the way to it reaches nobody. That is not a hazard this page merely avoids;
    it is a defect found across it — ``✅ … parameters loaded!``,
    ``✅ Parameters reloaded from cache!``, ``✅ Parameters saved to cache!``,
    ``✅ Cache cleared!`` and ``✅ Parameters reset to the settings MAFigate opens with``
    were each drawn immediately before a rerun. ``AppTest`` reports such an element anyway
    when it sits inside a container, which is how they all read as working.

    Args:
        confirmation: one sentence saying what happened, ending in a full stop.
        arm: the arm the user has just been moved to, or ``None`` where nothing moved.
    """
    st.session_state[PARAM_NOTICE_KEY] = {
        "confirmation": confirmation,
        "switched_to": arm,
    }


def adopt_parameters(params, confirmation, then=None):
    """Replace the session's parameters, and say so — naming the arm if it moved.

    One rule with five call sites, not five fixes. The presets tab, an uploaded custom
    preset, ``📥 Reload from Cache``, *Reset to Defaults* and the mismatch notice's arm
    switch all replace ``filter_params`` *whole*. All but the fourth can land a set that
    carries a different ``sample_type``, and so move the user between the germline and
    somatic arms; that one is handed the arm the page is on and cannot. **No path changes
    the arm without saying so**, and none replaces the parameters without saying that
    either (issue #133).

    That mattered because the arm is not one setting among many: it decides which
    guideline sources are read, which thresholds apply, which gene panel is offered and
    which columns the filter requires. A germline MAF filtered on the somatic arm draws
    ``❌ PATHOGENIC RETENTION DEGRADED`` for a ``CancerVar`` column it was never supposed
    to carry — a message that is correct about the somatic arm and blames the file.

    The notice is *parked* rather than drawn here, for the reason
    :func:`_adopt_uploaded_parameters` already had to learn: ``st.rerun`` discards the
    frame, so anything written on the way to it is thrown away. That was not a hazard this
    change avoided but a defect it found — the ``✅ … parameters loaded!`` and
    ``✅ Parameters reloaded from cache!`` lines were both drawn immediately before a
    rerun, so **neither ever reached the user**. ``AppTest`` reports them anyway when they
    sit inside a container, which is how they read as working.

    Args:
        params: the parameter set to adopt, used as-is. A preset is a complete set for its
            own arm, so it is not re-seeded from that arm's contract the way the Sample
            Type control re-seeds: re-seeding would discard the very preset the user asked
            for by name.
        confirmation: one sentence from the call site saying what was loaded, ending in a
            full stop. The arm clause in front of it is this function's, so the paths
            cannot drift into five vocabularies for one event.
        then: optional work to do once the new parameters are in place and before the
            rerun. Issue #135's arm switch needs it: that path re-cuts the open file and
            opens its report, and both read ``filter_params``, so they have to run after
            the replacement — while the rerun below still has to be the *last* thing, since
            it discards the frame. Without this hook that call site would have to replace
            the dict itself, which is the one thing this function exists to stop, and
            ``test_only_the_sanctioned_places_replace_the_parameters_wholesale`` would have
            been asked to sanction a fifth writer that really can move the arm.
    """
    # Both arms are read the way every other reader reads them — `.get("sample_type",
    # "somatic")`, as in `MAFigate.initialize_session_state`,
    # `data_loading.current_arm` and `show_parameter_config_page` below. Named by function
    # rather than by line: this comment carried three line numbers, and issue #135 shifted
    # two of them by adding code above. A set
    # that states no arm *is* somatic to the app, so comparing the raw values would let a
    # cache with no `sample_type` move a germline session to somatic under a plain "loaded"
    # line: the arm changing without saying so, by the one route that could still do it.
    previous = st.session_state.get("filter_params")
    previous_arm = (
        previous.get("sample_type", "somatic")
        if isinstance(previous, dict) and previous
        else None
    )
    st.session_state.filter_params = params

    # No previous arm is not a switch. There is nothing to have moved away from, and
    # announcing one tells a user they were somewhere they have never been.
    arm = params.get("sample_type", "somatic")
    moved = previous_arm is not None and arm != previous_arm

    # After the replacement, because the work reads the parameters it just installed; before
    # the parking and the rerun, because the rerun discards this frame and anything the work
    # wants shown has to be stashed rather than drawn.
    if then is not None:
        then()

    park_confirmation(confirmation, arm=arm if moved else None)
    st.rerun()


def show_parameter_notice():
    """Draw what a replacement parked for us, once, on the run after the rerun.

    The arm clause borrows the Sample Type route's own ``🔄 Switched to {arm}`` opening
    rather than inventing a second vocabulary: a preset switch is the same event reached
    differently. It cannot borrow that route's *sentence*, which promises "parameters
    reset to the settings MAFigate opens with" — true of the selectbox, false of a preset,
    which brings its own.

    The consequence sentence names the open file only when there is one. The parameter
    page is reachable with nothing loaded, and "MAFigate now filters this file" would be a
    claim about a file that is not there.
    """
    notice = st.session_state.pop(PARAM_NOTICE_KEY, None)
    if not notice:
        return

    arm = notice["switched_to"]
    if not arm:
        st.success(f"✅ {notice['confirmation']}")
        return

    if st.session_state.get("maf_data") is not None:
        consequence = f"MAFigate now filters this file as {arm}."
    else:
        consequence = f"MAFigate now filters as {arm}."
    st.info(f"{switched_to(arm)} {notice['confirmation']} {consequence}")


def show_parameter_config_page():
    """Display the parameter configuration page."""

    st.title("⚙️ Parameter Configuration")
    st.markdown("Configure filtering parameters for your MAF analysis.")

    show_discarded_cache_banner()
    show_parameter_notice()
    show_migration_report()

    # Initialize parameters - try loading from cache first
    if (
        "filter_params" not in st.session_state
        or st.session_state.filter_params is None
    ):
        # Only a cache this version wrote is restored; an older one has already been set
        # aside by the banner above, which is why this cannot silently reopen off parity.
        cached_params = load_parameters_from_cache()
        if cached_params:
            st.session_state.filter_params = cached_params
            # Show info message about loaded cache
            cache_info = get_cache_info()
            if cache_info:
                cache_time = datetime.fromisoformat(cache_info["timestamp"]).strftime(
                    "%Y-%m-%d %H:%M"
                )
                st.info(f"💾 Loaded cached parameters from {cache_time}")
        else:
            # `pipeline_params` returns a fresh deep copy, and both halves matter: this
            # page edits the dict in place and the keep-lists are nested, so anything
            # shallower would let a widget append to the contract's own list and redefine
            # the app's default for the rest of the process.
            st.session_state.filter_params = pipeline_params("somatic")

    # A `params_cached` flag was initialised to `False` here under the comment "Initialize
    # cache tracking", and those two lines were its only mentions in the whole repository —
    # nothing set it when parameters were cached and nothing read it to find out. Deleted in
    # passing at issue #167, which was about the same shape one page over: a session key read
    # by a branch that nothing else writes. What actually tracks the cache is
    # `page_modules/param_store.py`, which stamps `cached_parameters.json` on disk.

    # Sample type selection
    st.subheader("🧬 Sample Type")
    previous_arm = st.session_state.filter_params.get("sample_type", "somatic")
    sample_type = st.selectbox(
        "Select analysis type:",
        SAMPLE_TYPES,
        index=safe_index(SAMPLE_TYPES, previous_arm, "somatic"),
        help="Choose between somatic (cancer) or germline (hereditary) variant analysis",
    )

    if sample_type != previous_arm:
        # Switching arms re-seeds from *that arm's* contract, so the app opens at parity
        # on whichever arm the user is on — which is what "a freshly loaded MAF is at
        # parity without the user doing anything" has to mean for a germline MAF.
        #
        # Carrying the old arm's dict across was silently off parity, because the two
        # arms share almost no guideline sources: the somatic default carries no
        # `filter_intervar` or `filter_renovo`, so both arrived at the germline tab
        # absent, became empty, and were written back. Measured on
        # germline_reference.maf, that is a criteria path of 13 against the contract's
        # 27 — and the all-empty warning cannot catch it, because `filter_clinvar` is
        # shared and stays populated.
        #
        # Re-seeding discards edits made on the other arm. That is the right trade and
        # not merely the easy one: the arms' parameters are mostly disjoint, so there is
        # nothing coherent to carry, and the alternative is the silent narrowing above.
        # A deviation is one click away in the Presets tab.
        st.session_state.filter_params = pipeline_params(sample_type)
        # Drawn inline rather than parked: this route does not rerun, so its frame is the
        # one the user sees. The opening is shared with `show_parameter_notice`
        # and the sentence deliberately is not — this control really does reset to the
        # settings MAFigate opens with, and a preset brings its own.
        st.info(
            f"{switched_to(sample_type)} parameters reset to the settings MAFigate "
            f"opens with for {sample_type} analysis."
        )
    else:
        st.session_state.filter_params["sample_type"] = sample_type

    # Create tabs for different parameter categories
    tab1, tab2, tab3, tab4, tab5 = st.tabs(
        [
            "📋 Presets",
            "📊 Basic Filters",
            "🧬 Clinical Filters",
            "🎯 Population Filters",
            "💾 Cache Settings",
        ]
    )

    with tab1:
        show_parameter_presets_tab(sample_type)

    with tab2:
        show_basic_filters_tab(sample_type)

    with tab3:
        show_clinical_filters_tab(sample_type)

    with tab4:
        show_population_frequency_tab(sample_type)

    with tab5:
        show_cache_management_tab()

    # Validate multiselect parameters
    st.session_state.filter_params = validate_multiselect_params(
        st.session_state.filter_params
    )

    # Auto-save parameters to cache when they change
    if "last_params_hash" not in st.session_state:
        st.session_state.last_params_hash = None

    # Calculate hash of current parameters to detect changes
    current_params_hash = str(hash(str(sorted(st.session_state.filter_params.items()))))

    # Save to cache if parameters have changed
    if st.session_state.last_params_hash != current_params_hash:
        # Save current state as new hash before showing message
        prev_hash = st.session_state.last_params_hash
        st.session_state.last_params_hash = current_params_hash

        # Save to cache
        if save_parameters_to_cache(st.session_state.filter_params):
            # Only show message if this isn't the first load
            if prev_hash is not None:
                st.sidebar.success("💾 Parameters auto-saved", icon="💾")

    # Save parameters
    st.markdown("---")
    st.subheader("💾 Save & Export Parameters")

    # Current parameters are automatically saved to session state
    st.info(
        "ℹ️ Your current parameters are automatically saved. Use the export buttons below to download them as files."
    )

    # Export buttons
    col_json, col_yaml, col_reset = st.columns(3)

    with col_json:
        st.download_button(
            label="📄 Download JSON",
            data=export_text(st.session_state.filter_params, "json"),
            file_name="mafigate_parameters.json",
            mime="application/json",
            use_container_width=True,
            help="Download current parameters as JSON file",
        )

    with col_yaml:
        st.download_button(
            label="📄 Download YAML",
            data=export_text(st.session_state.filter_params, "yaml"),
            file_name="mafigate_parameters.yaml",
            mime="text/yaml",
            use_container_width=True,
            help="Download current parameters as YAML file",
        )

    with col_reset:
        if st.button(
            "🔄 Reset to Defaults",
            use_container_width=True,
            help="Put every parameter back to the settings MAFigate opens with.",
        ):
            reset_parameters(sample_type)

    # Navigation
    st.markdown("---")
    if st.button("📊 Go to Data Analysis", use_container_width=True, type="primary"):
        st.session_state.current_page = "data_loading"
        st.rerun()

    # Display current parameters
    with st.expander("🔍 Current Parameters Preview"):
        st.json(st.session_state.filter_params)


def show_parameter_presets_tab(sample_type):
    """Show parameter preset loading options."""

    st.subheader("📋 Parameter Presets")
    st.markdown("Load predefined parameter sets for common analysis scenarios.")

    # Four presets, two per arm, described by the cut each one makes rather than by how far
    # it sits from a named baseline. The two that reloaded the app's own opening settings are
    # gone: the app already opens on them and "Reset to Defaults" already restores them, so
    # a preset for the same thing was a third door into the same room (issue #51).
    presets = {
        "Broad Somatic": {
            "description": (
                "A wide net for a first pass. Keeps uncertain and conflicting clinical "
                "calls beside the confident ones — CancerVar down to Tier III, CIViC down "
                "to level D, ESCAT tiers I to III — and every ClinVar record except a "
                "benign or likely-benign one, so a risk, protective, drug-response or "
                "unclassified annotation is no longer a reason to set a variant aside. "
                "Sets aside only intergenic and RNA "
                "calls, and drops a variant when every frequency column your file fills in "
                "puts it above 1% of the population — a blank column does not speak for a "
                "variant, and a variant with no frequency anywhere is kept. Use it when you "
                "would rather review a variant than miss one."
            ),
            "file": None,
            "params": copy.deepcopy(SOFT_SOMATIC_PARAMS),
        },
        "Broad Germline": {
            "description": (
                "A wide net for a first pass. Keeps variants of uncertain significance "
                "beside the confident calls, accepts a pathogenic RENOVO prediction at any "
                "confidence, and keeps every ClinVar record except a benign or "
                "likely-benign one, so a risk, protective, drug-response or unclassified "
                "annotation is no longer a reason to set a variant aside. "
                "Sets aside only intergenic and RNA calls, and "
                "drops a variant when every frequency column your file fills in puts it "
                "above 1% of the population — a blank column does not speak for a variant, "
                "and a variant with no frequency anywhere is kept. "
                "Use it when you would rather review a variant than miss one."
            ),
            "file": None,
            "params": copy.deepcopy(SOFT_GERMLINE_PARAMS),
        },
        "Stringent Somatic": {
            "description": (
                "A short list of what you would act on. Keeps only pathogenic and "
                "likely-pathogenic calls with the strongest therapeutic evidence — "
                "CancerVar Tier I and II, CIViC A and B, ESCAT IA to IC — and sets aside "
                "synonymous, splice-region, UTR, flanking, intronic, intergenic and RNA "
                "calls. Drops a variant when every frequency column your file fills in puts "
                "it above 0.5% of the population — a blank column does not speak for a "
                "variant, and a variant with no frequency anywhere is kept."
            ),
            "file": None,
            "params": copy.deepcopy(CLINICAL_SOMATIC_PARAMS),
        },
        "Stringent Germline": {
            "description": (
                "A short list of what you would act on. Keeps only pathogenic and "
                "likely-pathogenic calls, and only high-confidence RENOVO predictions. "
                "Sets aside synonymous, splice-region, UTR, flanking, intronic, intergenic "
                "and RNA calls, and drops a variant when every frequency column your file "
                "fills in puts it above 0.5% of the population — a blank column does not "
                "speak for a variant, and a variant with no frequency anywhere is kept."
            ),
            "file": None,
            "params": copy.deepcopy(CLINICAL_GERMLINE_PARAMS),
        },
    }

    # Display presets
    col1, col2 = st.columns(2)

    with col1:
        st.markdown("**Available Presets:**")

        for preset_name, preset_info in presets.items():
            with st.expander(f"📌 {preset_name}"):
                st.write(preset_info["description"])

                if st.button(
                    f"Load {preset_name}",
                    key=f"load_{preset_name.replace(' ', '_').lower()}",
                ):
                    # Already a private deep copy — the presets are built fresh on every
                    # render, so a widget editing this dict cannot reach back into the
                    # module constant behind it.
                    #
                    # There is no from-a-file branch here any more. Every preset declared
                    # `"file": None`, so `load_preset_from_file` was unreachable — and
                    # being unreachable, it was also the one migration path nobody could
                    # notice going stale (issue #40).
                    #
                    # All four presets are offered on both arms, deliberately: the labels
                    # carry the arm, so "Load Stringent Somatic" from the germline arm is
                    # a one-click route someone can mean. What was missing was the app's
                    # half of that exchange, which `adopt_parameters` now holds.
                    adopt_parameters(
                        preset_info["params"], f"{preset_name} parameters loaded."
                    )

    with col2:
        st.markdown("**Custom Presets:**")

        # Upload custom preset
        uploaded_preset = st.file_uploader(
            "Upload Custom Preset",
            type=["json", "yaml", "yml"],
            help=(
                "Upload a JSON or YAML file with custom parameters. A file written by an "
                "older MAFigate — including a superseded parameter cache — is migrated "
                "onto MAFigate's current parameter set, and anything dropped is named."
            ),
        )

        if uploaded_preset is not None:
            try:
                if uploaded_preset.name.endswith(".json"):
                    document = json.load(uploaded_preset)
                else:
                    document = yaml.safe_load(uploaded_preset)

                if st.button("Load Custom Preset"):
                    _adopt_uploaded_parameters(document, uploaded_preset.name)

            except Exception as e:
                st.error(f"❌ Error loading preset: {str(e)}")

        # Export current parameters
        st.markdown("**Export Current Parameters:**")

        col_json, col_yaml = st.columns(2)

        with col_json:
            st.download_button(
                label="📄 Download as JSON",
                data=export_text(st.session_state.filter_params, "json"),
                file_name="custom_filter_params.json",
                mime="application/json",
                use_container_width=True,
            )

        with col_yaml:
            st.download_button(
                label="📄 Download as YAML",
                data=export_text(st.session_state.filter_params, "yaml"),
                file_name="custom_filter_params.yaml",
                mime="text/yaml",
                use_container_width=True,
            )


#: Where an upload's migration report waits for the render that can show it. See
#: :func:`_adopt_uploaded_parameters`.
MIGRATION_REPORT_KEY = "last_migration_report"


def _adopt_uploaded_parameters(document, filename):
    """Migrate an uploaded parameter file into the session, and report what it cost.

    The report is not decoration. Every one of the migration's per-key rules is a *silent
    widening* if it goes wrong, and several of them legitimately discard something the user
    chose — the other arm's gene list, a germline ESCAT selection, a key this app no longer
    has a parameter for. A user who is told the file was "loaded" and then sees a report
    they do not recognise has nothing to go on. "Unknown keys dropped **and named**" is an
    acceptance criterion, not a nicety.

    Which is why the report is *stored* rather than written here. ``st.rerun`` raises
    immediately and the run's output is discarded, so anything drawn on the way to it is
    thrown away — the messages would flash at best. The rerun is still wanted, because the
    widgets have to be rebuilt around the new dict. So the report is parked in session
    state and :func:`show_migration_report` draws it on the run that follows.
    """
    result = migrate_params(document)
    st.session_state[MIGRATION_REPORT_KEY] = {
        "filename": filename,
        "notes": list(result.notes),
        "dropped": list(result.dropped),
    }
    if result.migrated:
        confirmation = (
            f"{filename} was written by an older MAFigate and has been migrated onto "
            "MAFigate's current parameter set."
        )
    else:
        confirmation = f"Parameters loaded from {filename}."
    # The headline — and the rerun with it — belongs to `adopt_parameters`, because an
    # uploaded file states its own arm and `migrate_params` honours it. It only ever named
    # `sample_type` when the value was *not* an arm; a file that legitimately said the
    # other arm moved the user with the report saying nothing, since a changed key is not
    # a dropped one.
    adopt_parameters(result.params, confirmation)


def show_migration_report():
    """Draw the report an upload parked for us, once, on the run after the rerun.

    Every dropped key is listed, not only the ones a rule had something to say about. Some
    are dropped by a rule that needs no explanation — the other arm's gene list when the
    file already carries ``filter_genes``, a saved panel name — and those are exactly the
    ones that would otherwise disappear without a line anywhere.
    """
    report = st.session_state.pop(MIGRATION_REPORT_KEY, None)
    if not report:
        return

    # The headline moved to `show_parameter_notice`, drawn immediately above this,
    # so that one sentence can carry the arm clause every replacement route shares. What is
    # left here is the detail only an upload has.
    for note in report["notes"]:
        st.info(note)

    dropped = report.get("dropped") or []
    if dropped:
        st.info(
            "Keys from the file that are not parameters here, and were dropped: "
            + ", ".join(f"`{key}`" for key in dropped)
            + "."
        )


def safe_index(options_list, value, default):
    """Safely get the index of a value in a list, return index of default if not found.

    ``default`` is required. It used to default to ``"All"``, from when every option list
    began with that — a fallback that now names nothing, and would silently become index
    ``0`` rather than the caller's intent (issue #36).
    """
    try:
        return options_list.index(value)
    except ValueError:
        try:
            return options_list.index(default)
        except ValueError:
            return 0


def show_basic_filters_tab(sample_type):
    """Show basic filtering options including gene filters."""

    st.subheader("📊 Basic Filtering Parameters")

    col1, col2 = st.columns(2)

    with col1:
        # It was labelled "Minimum Depth (DP)" and this gate has never read `DP`: the
        # vendored filter's clause is `(t_alt_count + t_ref_count) >= coverage`, on both arms.
        # Nor is `DP` that sum under another name, which is what made the label plausible
        # enough to survive — issue #127 measured 322,913 rows across 157 real MAFs and the two
        # numbers differ on 72.7% of them, `DP` being the larger on 57.3%, because `DP` is the
        # depth across every sample and this counts the tumour's reads alone. At this widget's
        # own default of 50 the two answers disagree on 2,117 of those rows. So the label named
        # a column the user could see in the table, holding a number this setting does not use.
        st.session_state.filter_params["min_depth"] = st.number_input(
            label_of("min_depth"),
            min_value=0,
            max_value=1000,
            value=st.session_state.filter_params.get("min_depth", 50),
            help=(
                "Minimum reads at the variant in the tumour — its alt and ref counts added "
                "together — for it to be kept. This is not the MAF's `DP` column, which "
                "counts every sample and is usually higher."
            ),
        )

        # An *exclude* list, which is what the pipeline's parameter of this name means:
        # the classifications named here are dropped and everything else is kept. It read
        # as an include list until issue #36, which was divergence #1 — and an include
        # list cannot be repaired by listing more values, because this page runs before
        # any MAF is loaded, so its vocabulary can never cover what a file actually
        # carries. 211 reference rows hold five classifications the list below has never
        # heard of, two of them minted by the pipeline itself; excluding keeps them.
        #
        # Empty is a real value here, and it means exclude nothing. No back-fill.
        st.session_state.filter_params["filter_variant_classification"] = st.multiselect(
            label_of("exclude_classifications"),
            VARIANT_CLASSIFICATIONS,
            default=filter_terms(
                st.session_state.filter_params.get("filter_variant_classification"),
                VARIANT_CLASSIFICATIONS,
            ),
            help=(
                "Variants with these functional classifications are dropped. "
                "Anything not listed is kept, including classifications this list does "
                "not name. Leave empty to exclude nothing."
            ),
        )

    with col2:
        if sample_type == "somatic":
            st.session_state.filter_params["vaf_threshold"] = st.number_input(
                label_of("vaf_threshold"),
                min_value=0.0,
                max_value=1.0,
                value=st.session_state.filter_params.get("vaf_threshold", 0.05),
                step=0.0001,
                format="%.4f",
                help="Minimum variant allele frequency for somatic variants",
            )
        else:
            st.session_state.filter_params["vaf_threshold_germline"] = st.number_input(
                label_of("vaf_threshold_germline"),
                min_value=0.0,
                max_value=1.0,
                # 0.2, not 0.3: this fallback was the third value in a three-way
                # disagreement — nextflow.config and both germline presets say 0.2.
                # Pinned to the contract by tests/test_param_contract.py.
                value=st.session_state.filter_params.get("vaf_threshold_germline", 0.2),
                step=0.0001,
                format="%.4f",
                help="Minimum variant allele frequency for germline variants",
            )

        # Pathogenic variant retention. Written under the contract's key and polarity,
        # which is the *only* spelling the filter reads once a dict carries both — see
        # _pathogenic_retention_on.
        retain = st.checkbox(
            f"🔬 {label_of('skip_pathogenic')}",
            value=_pathogenic_retention_on(st.session_state.filter_params),
            help="Automatically keep pathogenic variants even if they don't meet other filter criteria (CancerVar Tier I/II, InterVar Pathogenic/Likely pathogenic, ClinVar pathogenic, CIViC A/B evidence)",
        )
        st.session_state.filter_params["skip_pathogenic"] = not retain
        # One truth per concept in the live dict. Leaving the old spelling behind is what
        # made this control dead: the filter prefers `skip_pathogenic` when both are
        # present, so a stale `keep_pathogenic` would sit in every export disagreeing with
        # the setting that was actually applied. Uploaded files are converted on the way
        # in (issue #40); this is what keeps a *preset* loaded mid-session from putting the
        # old spelling back.
        st.session_state.filter_params.pop("keep_pathogenic", None)

    _show_gene_filters(sample_type)


def _show_gene_filters(sample_type):
    """The gene controls: a panel dropdown, a paste box, and an upload.

    One arm-agnostic section where there used to be two near-identical copies, because the
    two copies were the reason a fix could land on one arm and not the other.

    Everything the user does here resolves to exactly **one** filter parameter —
    ``filter_genes``, a list of symbols — via ``filters.gene_lists.parse_gene_list``,
    which is the same tokeniser the filter seam runs. Nothing else about the gene
    selection is a parameter: the panel *name* is UI state (see
    :func:`gene_panel_state_key`), and the per-arm ``somatic_genes``/``germline_genes``
    strings this page used to write are gone. Sharing the tokeniser is what makes the
    label honest — it used to say "comma-separated" while every gene file in the project
    is one symbol per line, and the parser believed the label.
    """
    st.markdown("---")
    st.subheader(f"🧬 {label_of('gene_selection')}")

    _adopt_incoming_gene_list(sample_type)

    col_panel, col_input = st.columns([1, 2])

    with col_panel:
        # No `index=`: the widget's key is the single source of truth for the choice, and
        # _adopt_incoming_gene_list is what puts an incoming list *into* that key rather
        # than fighting it with an `index=`. See gene_panel_state_key for the cost of both.
        panel = st.selectbox(
            "Gene Set",
            options=list(GENE_LISTS.keys()),
            help="Choose a predefined panel, or 'Custom' to paste or upload your own list",
            key=gene_panel_state_key(sample_type),
        )

    with col_input:
        if panel == "Custom":
            selection = _custom_gene_selection(sample_type)
        elif panel in PREDEFINED_GENE_SETS:
            selection = _panel_gene_selection(panel)
        else:  # "All"
            selection = GeneSelection()
            st.info("🌐 Using all genes (no gene filtering)")

        # A list, always — never a joined string. A single-element list joined back into a
        # scalar is what makes the adapter write one character per line.
        st.session_state.filter_params["filter_genes"] = list(selection.symbols)
        # Record what the widgets now represent, so the next run can tell a list *we*
        # produced from one that arrived from somewhere else. See
        # _adopt_incoming_gene_list.
        st.session_state[_gene_adoption_stamp_key(sample_type)] = selection.symbols

        for message in selection.messages():
            st.warning(message)


def _gene_adoption_stamp_key(sample_type):
    """Where we record the gene list the widgets currently represent, per arm."""
    return f"gene_widgets_show_{sample_type}"


#: Distinguishes "no stamp yet" from "stamped with an empty gene list". A tuple can never
#: equal this, so the first render always adopts.
_NEVER_STAMPED = "__never_stamped__"


def _adopt_incoming_gene_list(sample_type):
    """Put a gene list that arrived from outside the page *into* the gene widgets.

    Without this the page silently discards one. The widgets own their state by key, so a
    ``filter_genes`` restored from the ``~/.mafigate`` cache or from an imported JSON/YAML
    reaches a dropdown that has never heard of it, and the assignment at the end of
    :func:`_show_gene_filters` then overwrites the real list with ``[]``. Measured before
    the fix: seeding ``filter_genes = ["TP53", "BRCA1"]`` and rendering the tab once left
    ``filter_genes == []`` and the page reading "Using all genes" — a silent widening, which
    is the one direction issue #28 says must never be silent.

    The naive repair — seed only when the key is absent — fixes the first render and leaves
    the same bug for every *later* import in the same session. So instead the page records
    what its widgets currently represent (the stamp) and adopts whenever the incoming list
    differs from it:

    * a list the page itself just produced equals the stamp, so nothing is adopted and a
      half-typed box is never yanked out from under the typist;
    * a list from anywhere else differs, so it is adopted and shown.

    A list that matches a named panel exactly is adopted *as that panel*, which is what
    makes a saved panel choice survive a round trip through a parameter file that only ever
    stored its symbols.
    """
    params = st.session_state.filter_params
    incoming = parse_gene_list(
        params.get("filter_genes")
        if params.get("filter_genes") is not None
        else params.get(f"{sample_type}_genes", "")
    ).symbols

    stamp_key = _gene_adoption_stamp_key(sample_type)
    if st.session_state.get(stamp_key, _NEVER_STAMPED) == incoming:
        return

    st.session_state[stamp_key] = incoming
    panel = _panel_matching(incoming)
    st.session_state[gene_panel_state_key(sample_type)] = panel
    if panel == "Custom":
        st.session_state[gene_text_state_key(sample_type)] = "\n".join(incoming)


def _panel_matching(symbols):
    """The panel option that denotes exactly ``symbols`` — ``All``, a panel name, or ``Custom``.

    Compared case-insensitively as sets, because that is the granularity at which the
    vendored clause treats two lists as the same restriction.
    """
    if not symbols:
        return "All"
    wanted = {symbol.upper() for symbol in symbols}
    for name, genes in PREDEFINED_GENE_SETS.items():
        if wanted == {str(gene).upper() for gene in genes}:
            return name
    return "Custom"


def _custom_gene_selection(sample_type):
    """A pasted and/or uploaded gene list, tokenised.

    **The text box is the single source of truth, and an upload *fills* it.** So there is no
    precedence rule to learn and nothing hidden: whatever is about to be filtered on is on
    screen, editable, in one place. An upload that instead overrode the box would mean the
    page showed one list and filtered another.

    The box is rendered above the uploader but *evaluated* after it, via a container, so a
    freshly uploaded file can be written into the box's state before the box is drawn.
    Filling is guarded by the file's identity rather than done on every rerun — otherwise
    the upload would overwrite each subsequent edit the user made to the box it filled.
    """
    typed_key = gene_text_state_key(sample_type)
    box = st.container()

    upload = st.file_uploader(
        "Upload a gene list to fill the box",
        type=list(GENE_LIST_EXTENSIONS),
        key=gene_upload_state_key(sample_type),
        help=f"A list file ({', '.join(GENE_LIST_EXTENSIONS)}) with one symbol per line",
    )
    if upload is not None:
        filled_key = f"{typed_key}__filled_from"
        if st.session_state.get(filled_key) != upload.file_id:
            # errors="replace" rather than a refusal: a stray byte in one line must not
            # cost the user the other 499 symbols, and the mangled token is then reported
            # as absent from the MAF rather than silently dropped.
            st.session_state[typed_key] = upload.getvalue().decode(
                "utf-8", errors="replace"
            )
            st.session_state[filled_key] = upload.file_id
        st.caption(f"📄 Filled from `{upload.name}` — edit it in the box above.")

    with box:
        st.text_area(
            "Genes",
            key=typed_key,
            help=(
                "One symbol per line, or separated by commas, semicolons or spaces. "
                "A pasted column heading is dropped. Matching ignores case. "
                "Leave empty for all genes."
            ),
        )

    return parse_gene_list(st.session_state.get(typed_key, ""))


def _panel_gene_selection(panel):
    """A named panel, resolved to its symbols through the one tokeniser."""
    symbols = panel_symbols(panel)
    if not symbols:
        st.warning(
            f"⚠️ The {panel} gene set has no genes defined yet "
            "(placeholder); no gene filtering will be applied."
        )
    else:
        st.info(f"📋 Using {panel} gene set ({len(symbols)} genes)")
        with st.expander(f"📝 View {panel} genes"):
            st.write(", ".join(symbols))
    return GeneSelection(symbols=symbols)


#: The ClinVar **classifications** control's help text, written once and read by both arms.
#:
#: ClinVar is the one guideline source drawn on each arm, so each of its tooltips is an entry
#: in the table below that exists twice — and the note on ``_GUIDELINE_CONTROLS`` says what a
#: silent divergence between two near-identical blocks has already cost this page.
#:
#: **There are two ClinVar tooltips since issue #103**, because there are two ClinVar controls:
#: the pathogenicity calls here, and everything ClinVar asserts that is not one in
#: :data:`_CLINVAR_OTHER_HELP` below. One multiselect holding the whole vocabulary buried the
#: calls a clinician reaches for first. The split follows the classes the report already
#: names — #98's boundary, which is the source's — so neither tooltip explains a grouping this
#: app invented.
#:
#: What needs saying is about what the user sees rather than about how matching works. A
#: ClinVar entry can hold more than one call at once, and each is judged separately — so a
#: variant reported as both pathogenic and likely pathogenic is kept by either choice on its
#: own, and there is no combined entry to look for.
#:
#: The rest is about *spelling*, which is unusual for a tooltip and deliberate: the vocabulary
#: holds pairs whose two entries are one ClinVar call written two ways, and a clinician who
#: does not know that will unselect one of each pair and silently lose every such variant in
#: files written the other way. There are two kinds — the conflicting-classification rename of
#: 2023 (issue #88) and the underscore-prefixed modifiers (issue #99) — and they differ in
#: provenance, not in what the user must do about them. **The two kinds now sit in different
#: controls**: the rename is here, named, because there is one of it and the names are what a
#: clinician recognises; every underscored pair is in the other control, so the rule about
#: them is stated there and not here, where it would describe terms this list does not hold.
#:
#: Neither tooltip spells a raw annotation value. #88 took ``Pathogenic/Likely_pathogenic`` out
#: of this one for exactly that reason — #79's bounded exception to the
#: no-implementation-vocabulary rule covers MAF *column names*, so a user can match them
#: against their own header row, and a ClinVar keep-term is not one. So the rename pair is
#: described rather than quoted, and the underscore pairs are described by their shape, which
#: is what a user comparing them in the list beside the tooltip can actually see.
_CLINVAR_HELP = (
    "ClinVar pathogenicity calls that keep a variant. A variant can carry more than one "
    "ClinVar call at once, and each is judged on its own — so a variant reported as both "
    "pathogenic and likely pathogenic is kept by either choice alone. Conflicting "
    "classifications and conflicting interpretations are ClinVar's newer and older names "
    "for one call, renamed in 2023: keep both selected unless you know which release your "
    "file was annotated against. The three uncertain-significance sub-tiers are a submitting "
    "laboratory's own subdivision of that call, not a scale ClinVar defines, so select them "
    "alongside it rather than instead of it."
)

#: The **other assertions** control's help text — the second ClinVar tooltip (issue #103).
#:
#: One thing distinguishes this control from the one beside it, and the tooltip leads with it:
#: these are not pathogenicity calls, so selecting one keeps a variant that may carry no
#: pathogenicity call at all. That is the only reason to open this control, and a label
#: reading *Other* does not convey it.
#:
#: The underscore rule lives here rather than in :data:`_CLINVAR_HELP` because every
#: underscored pair the vocabulary offers is in this list. Stated as a rule rather than by
#: naming members, so a further such term cannot make the text wrong — the same choice #99
#: made, now applied where the members actually are.
#:
#: It names the five classes in words rather than glossing each term. Which class a term is in
#: is ``components/clinical_summary.py``'s to say, and a per-term gloss here would be a second
#: copy of it — the drift the Help page was rewritten to stop. The Help page carries the list
#: itself, with the three record-state terms glossed from sources this repo holds.
_CLINVAR_OTHER_HELP = (
    "ClinVar assertions that are not pathogenicity calls: a disease risk, a drug response, "
    "an association with a trait, a protective effect, and the records where ClinVar holds "
    "no classification for this variant on its own. Selecting one keeps a variant on that "
    "assertion alone, so a variant carrying no pathogenicity call can still reach the "
    "report — which is the whole of what this control adds. The terms beginning with an "
    "underscore are the same terms as their unprefixed twins, in the form some files spell "
    "them; keep both of each pair selected, or a variant will be reported or missed "
    "according to how its file was written."
)

#: The guideline controls, per arm: parameter key, label, options, and help text. Each is
#: a *keep* list — a value listed is kept, an empty list places no constraint at all.
#:
#: The labels are read from :data:`config.param_labels.PARAM_LABELS` rather than written
#: here, because the results view names these same settings in its filter recap (issue
#: #137) and a second hand-written set of names is how the two would come to disagree.
#: ``tests/test_run_recap.py`` fails if a control's label stops coming from that table.
#:
#: A table rather than four near-identical blocks of widget code. The two blocks it
#: replaces had drifted into each other's shape well enough to hide that the germline one
#: rendered an ESCAT control the germline filter has no argument for.
_GUIDELINE_CONTROLS = {
    "somatic": [
        (
            "filter_cancervar",
            label_of("cancervar_keep"),
            CANCERVAR_OPTIONS,
            "CancerVar tiers that keep a variant in the report.",
        ),
        (
            "filter_civic",
            label_of("civic_keep"),
            CIVIC_OPTIONS,
            "CIViC evidence levels that keep a variant. Applied only when the MAF "
            "carries CIViC annotation.",
        ),
        (
            "filter_clinvar",
            label_of("clinvar_pathogenicity"),
            CLINVAR_PATHOGENICITY_TERMS,
            _CLINVAR_HELP,
        ),
        (
            "filter_clinvar",
            label_of("clinvar_other"),
            CLINVAR_OTHER_ASSERTION_TERMS,
            _CLINVAR_OTHER_HELP,
        ),
        (
            "filter_escat",
            label_of("escat_keep"),
            ESCAT_OPTIONS,
            # "IA=strongest evidence to V=case reports" stood here. The Help page glossed
            # the same scale as "V: Not actionable" and the Pathogenicity Overview maps V to
            # Unknown, so the app told the user three different things about one level and
            # nothing in this repository could say which was right (issue #79). #89 sourced
            # the definitions, and they live in `ESCAT_DEFINITIONS` for the help page to
            # render — deliberately not repeated here. Eight clinical definitions do not fit
            # a multiselect tooltip, and a second copy of them is how the contradiction got
            # in. This states the ordering, reads the strongest level off that same constant
            # rather than spelling it out again, and says where the definitions are.
            f"ESCAT actionability levels that keep a variant, from {ESCAT_STRONGEST} — the "
            "strongest evidence — downwards. The Help page defines each level.",
        ),
    ],
    "germline": [
        (
            "filter_intervar",
            label_of("intervar_keep"),
            INTERVAR_OPTIONS,
            "InterVar ACMG/AMP classifications that keep a variant.",
        ),
        (
            "filter_renovo",
            label_of("renovo_keep"),
            RENOVO_OPTIONS,
            "RENOVO pathogenicity classes that keep a variant.",
        ),
        (
            "filter_clinvar",
            label_of("clinvar_pathogenicity"),
            CLINVAR_PATHOGENICITY_TERMS,
            _CLINVAR_HELP,
        ),
        (
            "filter_clinvar",
            label_of("clinvar_other"),
            CLINVAR_OTHER_ASSERTION_TERMS,
            _CLINVAR_OTHER_HELP,
        ),
    ],
}


def show_clinical_filters_tab(sample_type):
    """Show clinical classification filtering options.

    Every control here contributes one clause to an **OR**: a variant is kept if any
    source keeps it. So emptying one source removes its clause and narrows the report,
    and emptying all of them leaves the guideline block admitting nothing — which is
    warned about below rather than prevented.
    """

    st.subheader("🏥 Clinical Classification Filters")
    st.caption(
        "A variant meets the criteria if **any** of these sources keeps it. "
        "An empty selection means that source places no restriction — it drops out of "
        "the comparison rather than widening it."
    )

    controls = _GUIDELINE_CONTROLS[sample_type]
    columns = st.columns(2)

    # A parameter may be drawn by more than one control, so a selection is *accumulated* per
    # key and written once, after every control feeding it has drawn: `filter_clinvar` is two
    # widgets since issue #103, the pathogenicity calls and the other assertions. Assigning
    # inside the loop — which is what stood here while every key had exactly one control —
    # would leave the second widget's value overwriting the first's, and the report would be
    # cut on half the user's selection with nothing on screen to say so.
    #
    # Each widget's own default does the *partition* with no extra code: `filter_terms` keeps
    # only what its options offer, so one saved `filter_clinvar` list seeds both widgets with
    # its own half. Both read `filter_params` before any of this writes to it, which is why
    # the write is after the loop rather than merely outside the `with`.
    selected = {}
    for position, (key, label, options, help_text) in enumerate(controls):
        with columns[position % 2]:
            chosen = st.multiselect(
                label,
                options,
                default=filter_terms(st.session_state.filter_params.get(key), options),
                help=help_text,
            )
        selected.setdefault(key, []).extend(chosen)

    st.session_state.filter_params.update(selected)

    _warn_if_every_guideline_source_is_empty(sample_type)


def _warn_if_every_guideline_source_is_empty(sample_type):
    """Warn — do not block — when no guideline source keeps anything.

    Deleting the catch-all sentinel made this state reachable for the first time, and the
    acceptance criterion is explicit that it is warned at the widget rather than refused.
    Refusing would make the app stricter than the pipeline, on whose own command line the
    state is perfectly expressible.

    What it produces is worth saying out loud, because it is not zero and it is not an
    error: the criteria path empties, and the report falls back to the variants
    pathogenic retention rescues on their own — 388 somatic and 496 germline rows on the
    GERSOM reference. A user who has emptied the last source otherwise sees a report that
    has quietly changed meaning, with nothing to say why.
    """
    sources = GUIDELINE_SOURCES[sample_type]
    if any(st.session_state.filter_params.get(key) for key in sources):
        return

    if not _pathogenic_retention_on(st.session_state.filter_params):
        st.warning(
            "⚠️ No guideline source keeps anything, and pathogenic retention is off. "
            "Nothing can meet the filter criteria, so the report will be empty."
        )
        return

    st.warning(
        "⚠️ Every guideline source is empty, so no variant can meet the filter "
        "criteria. The report falls back to the **pathogenic-rescue** floor — only "
        "variants pathogenic retention keeps on their own will appear."
    )


def show_population_frequency_tab(sample_type):
    """Show population frequency filtering options."""

    st.subheader("🎯 Population Frequency Filters")

    st.session_state.filter_params["max_freq_population"] = st.number_input(
        label_of("max_freq_population"),
        min_value=0.0,
        max_value=1.0,
        value=st.session_state.filter_params.get("max_freq_population", 0.01),
        step=0.001,
        format="%.4f",
        help=(
            "Maximum population frequency allowed (gnomAD, ExAC, etc.). Set it to 1.0 to "
            "switch this filter off entirely and judge every variant on the other filters "
            "alone."
        ),
    )

    st.info(
        "💡 **Population frequency filtering** helps identify rare variants by excluding "
        "common population variants. Lower values (e.g., 0.01 = 1%) are more restrictive. "
        "A variant is kept if **any** frequency column that has a value for it puts it at "
        "or below the threshold. A blank column is *not seen in this panel* rather than "
        "*common*, so it never counts against a variant — but it does not speak for one "
        "either, and a variant every populated column calls common is dropped. A variant "
        "with no frequency anywhere is kept. **ClinVar-pathogenic variants are exempt**, "
        "so a low-penetrance "
        "pathogenic allele is not dropped for being common; the filter reports how many it "
        "removed and how many it spared."
    )


def reset_parameters(sample_type):
    """Reset parameters to this arm's defaults — the settings the app opens on.

    "Defaults" means the contract, not SOFT (issue #36) — otherwise the one button labelled
    *reset* would be the one way to leave the app's own starting point without choosing to.

    This one cannot move the arm — it is handed the arm the page is on — but it goes through
    :func:`adopt_parameters` all the same, and not only for tidiness: its confirmation had
    the defect the arm change was hiding behind, drawn immediately before the ``st.rerun``
    that discards the frame. Pressing *Reset to Defaults* showed the user nothing at all.
    """
    adopt_parameters(
        pipeline_params(sample_type),
        f"Parameters reset to the settings MAFigate opens with for {sample_type} analysis.",
    )


def show_cache_management_tab():
    """Show cache management options."""

    st.subheader("💾 Parameter Cache Management")
    st.markdown(
        "Manage automatic saving and loading of your parameter settings between sessions."
    )

    # Cache information
    cache_info = get_cache_info()

    col1, col2 = st.columns(2)

    with col1:
        st.markdown("**Cache Status:**")
        if cache_info:
            cache_time = datetime.fromisoformat(cache_info["timestamp"]).strftime(
                "%Y-%m-%d %H:%M:%S"
            )
            st.success(f"✅ Cache available")
            st.write(f"📅 **Last saved:** {cache_time}")
            st.write(f"🔢 **File size:** {cache_info['file_size']} bytes")
            st.write(f"🏷️ **App version:** {cache_info['app_version']}")
            # The stamp, shown beside the app version because the two are not the same
            # thing and reading one as the other is what made the old caches unmigratable.
            st.write(f"🧾 **Parameter format:** {cache_info['schema_version']}")
        else:
            st.warning("⚠️ No cache found")
            st.write("Parameters will be cached automatically when you modify them.")

    with col2:
        st.markdown("**Cache Actions:**")

        # Manual save button
        if st.button(
            "💾 Save Current Parameters",
            help="Manually save current parameters to cache",
        ):
            if save_parameters_to_cache(st.session_state.filter_params):
                # Parked, not drawn: the rerun below discards this frame, so the old
                # inline `st.success` here reached nobody.
                park_confirmation("Parameters saved to your cache.")
                st.rerun()

        # Load from cache button (if cache exists)
        if cache_info and st.button(
            "📥 Reload from Cache", help="Reload parameters from cache"
        ):
            cached_params = load_parameters_from_cache()
            if cached_params:
                # A cache is auto-saved on every change, so it can hold whichever arm was
                # last configured — including from a previous session.
                adopt_parameters(
                    cached_params, "Parameters reloaded from your saved cache."
                )
            else:
                st.error("❌ Failed to load from cache")

        # Clear cache button
        if cache_info and st.button("🗑️ Clear Cache", help="Delete cached parameters"):
            if clear_parameters_cache():
                park_confirmation("Your saved cache has been cleared.")
                st.rerun()
            else:
                st.error("❌ Failed to clear cache")

    # Auto-save status
    st.markdown("---")
    st.markdown("**Auto-Save Feature:**")
    st.info(
        "🔄 **Auto-save is enabled** - Your parameters are automatically saved to cache "
        "whenever you make changes. They are restored when you restart the application, "
        "**as long as they were written in this version's parameter format**. A cache "
        "from an older format is set aside rather than restored, and the app opens at its "
        "own default settings instead."
    )

    # Cache location info
    with st.expander("🔧 Technical Details"):
        st.write(f"**Cache location:** `{PARAMS_CACHE_FILE}`")
        st.write(
            f"**Cache format:** JSON — parameters inside an envelope stamped "
            f"`{SCHEMA_VERSION_KEY}: {PARAM_SCHEMA_VERSION}`, which is the format's own "
            "version and moves independently of the app's."
        )
        st.write("**Privacy:** Cache is stored locally on your computer only")
        st.markdown(
            "💡 **Note:** The cache file is created in your home directory under `.mafigate/` "
            "to persist settings between application sessions."
        )
        st.markdown(
            f"🧾 **Superseded caches** are moved to a file ending `{SUPERSEDED_CACHE_SUFFIX}` "
            "beside this one rather than deleted. Upload one under *Presets → Upload "
            "Custom Preset* to migrate it onto the current parameter set."
        )
