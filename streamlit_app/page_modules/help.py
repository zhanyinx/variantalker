"""
Help page for MAFigate application with comprehensive column information.
"""

import streamlit as st
from config.columns import (
    get_column_description,
    create_column_info_table,
    resolve_visible_columns,
    COLUMN_SOURCE_HEADER,
    REQUIRED_COLUMNS,
)
from config.columns import COLUMN_DESCRIPTIONS
from config.vocabularies import (
    CANCERVAR_OPTIONS,
    CIVIC_DEFINITIONS,
    CIVIC_SOURCE,
    CIVIC_SOURCE_URL,
    CLINVAR_OPTIONS,
    CLINVAR_OTHER_ASSERTION_TERMS,
    CLINVAR_PATHOGENICITY_TERMS,
    ESCAT_DEFINITIONS,
    ESCAT_SOURCE,
    ESCAT_SOURCE_DOI,
    ESCAT_STRONGEST,
    INTERVAR_OPTIONS,
    RENOVO_OPTIONS,
    VARIANT_CLASSIFICATIONS,
)
from filters.absent_columns import PATHOGENIC_INPUTS, REQUIRED_INPUTS
from filters.variant_filters import FREQUENCY_COLUMNS

#: The columns whose absence is escalated, for the missing-column section below. Derived, not
#: listed: it is exactly the pathogenic-retention set intersected with the set that gets
#: filled, which is the same rule ``plan_fills`` applies — so the page cannot come to describe
#: a warning the app does not raise, or miss one it does.
_ESCALATED = sorted(
    {
        column
        for arm, columns in PATHOGENIC_INPUTS.items()
        for column in columns
        if column in REQUIRED_INPUTS[arm]
    }
)


#: The columns this page can honestly call display-only: the app expects them, no filter reads
#: them on any path, and the file is not refused without them.
#:
#: Derived rather than typed out (issue #127). Five names were spelled into the missing-column
#: section and into the reference table's reasoning, and a sixth added to ``REQUIRED_COLUMNS``
#: would have joined the app's behaviour without joining this page's account of it — which is
#: the shape of drift this repo keeps finding. The subtraction is the same one
#: ``validate_required_columns`` applies to decide what its note may name, over the union of
#: both arms because this prose is not drawn per arm.
#: The parentheses are load-bearing: ``-`` binds tighter than ``|``, so without them this reads
#: as ``filtering | somatic | (germline - core - inputs)`` and the somatic list keeps its filter
#: inputs.
_DISPLAY_ONLY = sorted(
    (
        set(REQUIRED_COLUMNS["filtering"])
        | set(REQUIRED_COLUMNS["somatic"])
        | set(REQUIRED_COLUMNS["germline"])
    )
    - set(REQUIRED_COLUMNS["core"])
    - {column for columns in REQUIRED_INPUTS.values() for column in columns}
)


#: The columns whose value decides how a variant is read clinically, as a shortlist for the
#: reader who does not want to search fifty-odd rows. The names only: each description is
#: read from :data:`~config.columns.COLUMN_DESCRIPTIONS` at render time.
#:
#: Five names arrived with issue #219, and every one is a column the variant detail panel began
#: drawing when issue #204 landed: ``ClinVar_VCF_CLNREVSTAT`` as ClinVar's star hierarchy,
#: ``RENOVO_pls`` inline behind the class, ``ESCAT_TISSUE`` and ``ESCAT_CANCER`` as the tumour
#: type inside the tier, and ``am_pathogenicity`` as the AlphaMissense band. Until #204 the
#: heading's claim was true of the nine columns it held; afterwards it was false of five it did
#: not. Each qualifier sits directly after the annotation it qualifies, which is where the panel
#: draws it.
#:
#: Two of the five are **not** rows of the reference table above — ``resolve_visible_columns``
#: carries neither ``ClinVar_VCF_CLNREVSTAT`` nor ``RENOVO_pls`` — so this shortlist is the only
#: place on the page a reader can look up the stars they just saw on a badge. Which columns the
#: app *shows* is behaviour map #199 does not touch; describing one it already draws is not.
#:
#: ``CIViC_Evidence_Level`` stays even though issue #202 took its badge off the panel. The claim
#: here is about reading a variant clinically, not about what the panel draws: the level still
#: decides somatic pathogenic retention, and ``CIViC_Variant_URL`` still draws in the panel's
#: external links.
_KEY_CLINICAL_COLUMNS = (
    "CancerVar",
    "ESCAT",
    "ESCAT_TISSUE",
    "ESCAT_CANCER",
    "CIViC_Evidence_Level",
    "ClinVar_VCF_CLNSIG",
    "ClinVar_VCF_CLNREVSTAT",
    "InterVar",
    "RENOVO_Class",
    "RENOVO_pls",
    "am_pathogenicity",
    "tumor_f",
    "gnomAD_exome_AF",
    "gnomAD_genome_AF",
)


#: The heading over the leftover-columns block on the Required Columns tab.
#:
#: A constant because ``tests/test_help_claims.py`` uses it as the sentinel that ends the two
#: per-arm blocks, and a heading that is a sentinel in one file and a literal in another is a
#: rename away from a guard that reads the wrong region — or, since the sentinel is found with
#: ``next()``, from a ``StopIteration`` reported as an error rather than as a failed claim. It read
#: "🔵 Clinical Annotation Columns" until issue #219, which found the block below it was neither
#: only clinical nor only annotations.
_OTHER_COLUMNS_HEADING = "🔵 Other columns MAFigate reads or shows"


def _columns(names) -> str:
    """Column names for the help prose, backticked."""
    return ", ".join(f"`{name}`" for name in names)


def _frequency_groups() -> str:
    """The frequency columns the filter reads, grouped by source, as markdown bullets.

    Derived from :data:`~filters.variant_filters.FREQUENCY_COLUMNS` rather than typed out. The
    same list was hand-written here until issue #219, which deleted the *other* copy of it — the
    derived one, in the block above — as a duplicate, and would have left the page's only
    frequency list the hand-written one. That is the drift issue #126 fixed in the opposite
    direction: the copy it deleted had omitted ``Freq_esp6500siv2_all`` and named ``MAX_AF``,
    which the filter does not read.

    A column matching no known source is listed rather than dropped, under a label that claims
    nothing about it — a silent omission here is exactly what #126 found, and *"Older panels"*
    would be a guess about a column nobody has seen.
    """
    labels = (
        ("gnomAD exomes", lambda name: name.startswith("gnomAD_exome")),
        ("gnomAD genomes", lambda name: name.startswith("gnomAD_genome")),
        ("Older panels", lambda name: name.startswith("Freq_")),
    )
    remaining = list(FREQUENCY_COLUMNS)
    lines = []
    for label, matches in labels:
        named = [name for name in remaining if matches(name)]
        if not named:
            continue
        remaining = [name for name in remaining if name not in named]
        lines.append(f"    - {label}: {_columns(named)}")
    if remaining:
        lines.append(f"    - Also read: {_columns(remaining)}")
    return "\n".join(lines)


def show_help_page():
    """Display the help page with column information and requirements."""

    st.title("❓ Help & Documentation")

    # Two lines stood between this title and the nav row, and both said what something
    # already on screen says. `Comprehensive guide for **MAFigate v2.0.0**` went at issue
    # #71, which moved the version to the About dialog alone. The breadcrumb that replaced
    # it — "📍 You are here: **Help & Documentation** …" — named the page a third time, after
    # the sidebar's nav radio and the title directly above it, which is the same redundancy
    # #61 removed from Home rather than reworded.

    # Navigation buttons to go back to other pages. The fourth column held the caption
    # "*Navigate to other pages*", labelling three buttons that label themselves.
    col_nav1, col_nav2, col_nav3 = st.columns([1, 1, 1])
    with col_nav1:
        if st.button("🏠 Home", key="nav_home_from_help"):
            st.session_state.current_page = "home"
            st.rerun()
    with col_nav2:
        if st.button("⚙️ Configure", key="nav_config_from_help"):
            st.session_state.current_page = "parameter_config"
            st.rerun()
    with col_nav3:
        if st.button("📊 Data & Analysis", key="nav_data_from_help"):
            st.session_state.current_page = "data_loading"
            st.rerun()

    st.markdown("---")

    # A `help_tab_focus` session key was read here and branched on two values —
    # `filter_options` and `required_columns` — each drawing an `st.info` that named the tab
    # the reader should click, then clearing the request. Read, branches and key are deleted
    # together at issue #167, because **nothing in the app ever wrote it**: that read and its
    # two clears were the only three mentions anywhere in the app, so both sentences were
    # unreachable from every route onto this page.
    #
    # It was never the mechanism its name promised, either. `st.tabs` cannot be selected from
    # Python — measured at issue #161, which wanted to aim a sidebar button at a named Help tab
    # and found it is not a destination — so the ceiling on this branch was always a sentence
    # telling the reader which of five tabs to press. Dead code shaped like a feature is how
    # #161 came to believe the mechanism existed and had to measure to learn otherwise, which
    # is the cost this deletion is paying off.
    #
    # Giving it a writer was the alternative, and it is out of scope rather than merely
    # unargued: a route from the parameter page into Help decides that *what does this setting
    # do?* is answered on this page rather than in the tooltip beside the control, and that is
    # the Help-versus-tooltip question map #157 rules past its destination. #161 had just
    # removed a redundant second door onto this page; adding one back to carry a hint needs an
    # argument nothing in the app was making.
    #
    # The Required Columns hint carried a sentence corrected at issue #127. That correction is
    # not lost with it: `show_required_columns_tab` below states the same thing — the columns
    # MAFigate refuses a file without, and what each other absence costs — in the tab the hint
    # was pointing at.

    # Create tabs for different help sections
    tab_names = [
        "📋 Column Information",
        "📝 Required Columns",
        "🔍 Filter Options",
        "📖 User Guide",
        "❓ FAQ",
    ]

    tab1, tab2, tab3, tab4, tab5 = st.tabs(tab_names)

    with tab1:
        show_column_information_tab()

    with tab2:
        show_required_columns_tab()

    with tab3:
        show_filter_options_tab()

    with tab4:
        show_user_guide_tab()

    with tab5:
        show_faq_tab()


def show_column_information_tab():
    """Show comprehensive column information."""

    st.subheader("📋 Column Information")
    # "detailed descriptions for **all** columns used in MAFigate analysis" was not true of
    # this table, and the omission was not a small one. `create_column_info_table` is built
    # from `resolve_visible_columns` — the columns MAFigate *shows* — so the two columns the
    # frequency filter reads first, `gnomAD_exome_AF` and `gnomAD_genome_AF`, are absent from
    # it while their `_raw` variants are present. A reader searching this table for the
    # column the tab beside it recommends found nothing. The table is left as the resolver's
    # set, which is what makes it the columns you can actually put on screen; the sentence now
    # says so, and says where the filter's own inputs are covered.
    st.markdown(
        """
    Every column MAFigate can show you, with what it holds. Searchable by name or by
    description, and filterable by category.

    The columns the *filters* read are not all shown — population frequency in particular is
    read from columns that never reach the table. Those are listed under **Required Columns**
    and under **Filter Options**.
    """
    )

    # Search functionality
    search_term = st.text_input(
        "🔍 Search columns",
        placeholder="Enter column name or keyword...",
        help="Search for specific columns by name or description keywords",
    )

    # Category filter.
    #
    # No `except ImportError` around this any more, and no second rendering of the same
    # table behind it. The fallback it guarded — a 60-line `show_basic_column_info` under
    # "⚠️ Full table view requires pandas" — was unreachable twice over: this module imported
    # pandas at the top, so a missing pandas failed the page before any of it ran, and
    # `config.columns` imports pandas too. The top-level `import pandas as pd` went with it,
    # having no other user in this file.
    df = create_column_info_table()
    categories = ["All"] + sorted(df["Category"].unique().tolist())
    selected_category = st.selectbox("Filter by category", categories)

    # Apply filters. `regex=False` because this is a search box, not a pattern box: with the
    # default, a user typing `(` into it raised `re.error` out of `str.contains`, which
    # `route_to_page` caught and turned into "Error loading page 'help'" on the Home page.
    filtered_df = df.copy()

    if search_term:
        mask = filtered_df["Column Name"].str.contains(
            search_term, case=False, na=False, regex=False
        ) | filtered_df["Description"].str.contains(
            search_term, case=False, na=False, regex=False
        )
        filtered_df = filtered_df[mask]

    if selected_category != "All":
        filtered_df = filtered_df[filtered_df["Category"] == selected_category]

    # Display results
    st.markdown(f"**Showing {len(filtered_df)} of {len(df)} columns**")

    if len(filtered_df) > 0:
        # Create an interactive table
        st.dataframe(
            filtered_df,
            use_container_width=True,
            hide_index=True,
            column_config={
                "Column Name": st.column_config.TextColumn(
                    "Column Name", width="medium"
                ),
                "Category": st.column_config.TextColumn("Category", width="small"),
                "Sample Types": st.column_config.TextColumn(
                    "Sample Types", width="small"
                ),
                # Keyed by the constant rather than a literal — see COLUMN_SOURCE_HEADER
                # for why. Widened from "medium": the answers are sentences now, not
                # two-word statuses (issue #124).
                COLUMN_SOURCE_HEADER: st.column_config.TextColumn(
                    COLUMN_SOURCE_HEADER, width="large"
                ),
                "Description": st.column_config.TextColumn(
                    "Description", width="large"
                ),
            },
        )
    else:
        st.info("No columns match your search criteria.")

    # Quick reference for key columns. The curation is this page's — which nine of the
    # fifty-odd columns carry the clinical verdict — but the *descriptions* are read from
    # `COLUMN_DESCRIPTIONS`, which is what the searchable table above renders. They used to
    # be hand-written here, a second wording of nine rows the same screen already showed
    # from the constant, free to drift from it. The shortlist is the value; the prose is not.
    with st.expander("🔑 Key Clinical Columns Quick Reference"):
        for col in _KEY_CLINICAL_COLUMNS:
            st.markdown(f"• **{col}**: {get_column_description(col)}")


def show_required_columns_tab():
    """What MAFigate refuses a file without, and what it does when anything else is absent.

    Named for the tab, which keeps the word "Required" because that is what a user looking for
    it will scan for. What the tab *says* stopped depending on the analysis type at issue #127:
    requirement does not, only the annotations each arm's filter reads do.
    """

    st.subheader("📝 Required Columns")
    # "MAFigate requires specific columns depending on your analysis type" was the same
    # overclaim the reference table made, one tab over, and issue #127 settled it in both
    # places: requirement does not depend on the analysis type. Eight columns are required —
    # the file is refused without them, on either arm — and nothing else here is. What does
    # depend on the arm is which columns the filter *reads*, and an absent one of those is
    # filled rather than demanded.
    st.markdown(
        f"""
    {len(REQUIRED_COLUMNS["core"])} columns are genuinely required: without them MAFigate
    cannot open your file at all, on either arm. Everything else on this page is a column
    MAFigate *uses* — some filled with a stand-in value if your file does not carry them, some
    simply not shown — and the sections below say which is which, and what each absence costs
    you.
    """
    )

    # Quick troubleshooting section
    st.markdown("### 🚨 MAF File Loading Issues?")

    # Four of this block's nine lines were wrong, and two of them in the direction that
    # costs the user their data. See the notes on each below.
    with st.expander("🔧 Common Problems & Solutions", expanded=True):
        st.markdown(
            """
        **❌ File won't load or shows errors:**

        🔍 **Check these common issues:**

        1. **File Format**: Ensure your file is tab-separated (.maf, .tsv, .txt)
        2. **Column Names**: Verify column names match standard MAF format exactly
        3. **File size**: files arrive through the browser, where the limit is **200 MB**.
           For a larger MAF, ask whoever set MAFigate up to raise it
        4. **Encoding**: Use UTF-8 encoding (most text editors support this)

        **✅ Quick Fixes:**
        - Check that the first line contains column headers (not data)
        - Ensure there are no completely empty rows in your data
        - Comment lines starting with `#` are handled for you — you do not need to strip them
        """
        )

    # What went from the block above, and why:
    #
    # **"Keep files under 500MB for optimal performance"** — the third of three copies of a
    # number nothing in this repo sets. Streamlit's default `server.maxUploadSize` is 200 MB
    # and there is no `.streamlit/config.toml`, no flag in `run_mafigate.sh` and none in the
    # `Makefile`, so the app *refused* a 300 MB MAF on a screen that had just invited it. The
    # copy is corrected rather than the cap raised: a 500 MB upload has a memory cost this app
    # has never paid — the session holds several copies of the frame — and `README.md` already
    # documents both the real limit and `--server.maxUploadSize` for whoever runs the app. The
    # flag stays out of this page deliberately; a clinician does not launch the process.
    #
    # **"Files with many # comment lines at the top may need manual cleanup"**, and the Quick
    # Fix repeating it, told the user to do by hand what the app does for them. `utils.read_maf`
    # inherits the pipeline's rule — count the leading `#` lines, read the header after them —
    # and `read_maf_without_comment_lines` retries with every comment line removed if that
    # fails. So the line is now the reassurance it should always have been.
    #
    # **"Try opening your file in Excel and re-saving as 'Text (Tab delimited)'"** is deleted
    # outright as the one Quick Fix that can destroy the file it claims to repair. Excel
    # truncates silently past 1,048,576 rows, and it converts gene symbols to dates — `SEPT9`
    # and `MARCH1` are the standard examples. Advising a round-trip through it for clinical
    # variant data is worse than saying nothing.

    # Core required columns, read from the list `validate_required_columns` actually refuses on
    # rather than transcribed beside it. The eight names and their order matched when this was
    # checked; the point is that they cannot come apart later. It is the same reasoning as the
    # derivation behind `_ESCALATED` and the missing-column lists below.
    st.markdown("### 🔴 Core Required Columns")
    st.markdown(
        "Without these, MAFigate cannot open your file at all — it refuses it and says "
        "which are missing:"
    )

    for col in REQUIRED_COLUMNS["core"]:
        st.markdown(f"• **{col}**: {get_column_description(col)}")

    # Not "analysis-specific requirements", and no longer two hand-written lists (issue #127).
    # None of these is required, and the old lists claimed an arm for columns that do not have
    # one: `Tumor_Seq_Allele1`/`2` sat under Somatic and `n_alt_count`/`n_ref_count` under
    # Germline, while the pipeline emits all four on **both** arms — which is exactly what the
    # reference table beside this one reports as "Both", so the two surfaces disagreed on the
    # same screen about the same columns. `tumor_f` and the tumour read counts were in the
    # somatic list on the same false footing: both arms' filters read them.
    #
    # What is genuinely per-arm is which annotations *this arm's filter* reads and the other's
    # does not, and that is derivable rather than curatable — so it is derived, and a column
    # cannot drift into the wrong arm here. On today's vendored source: somatic `CancerVar`,
    # `ESCAT`, `CIViC_Evidence_Level`; germline `InterVar`, `RENOVO_Class`.
    _reads = {
        arm: set(REQUIRED_INPUTS[arm]) | set(PATHOGENIC_INPUTS[arm])
        for arm in ("somatic", "germline")
    }
    col1, col2 = st.columns(2)

    for column, arm, heading in (
        (col1, "somatic", "### 🟡 Somatic Analysis"),
        (col2, "germline", "### 🟢 Germline Analysis"),
    ):
        with column:
            st.markdown(heading)
            st.markdown(
                f"**Annotations only a {arm} report reads — none of them required, and an "
                "absent one is filled with a stand-in value:**"
            )
            other = "germline" if arm == "somatic" else "somatic"
            for col in sorted(_reads[arm] - _reads[other]):
                if col in COLUMN_DESCRIPTIONS:
                    st.markdown(f"• **{col}**: {get_column_description(col)}")

    # What is left once the two lists above and the frequency list below are accounted for
    # (issue #219). Four buckets stood here under the heading "Clinical Annotation Columns", and
    # 10 of their 12 entries were a second copy of a list already on this tab: `ESCAT`,
    # `CancerVar` and `CIViC_Evidence_Level` are rendered five lines above by the per-arm
    # derivation, and all seven `FREQUENCY_COLUMNS` are listed again fifty lines below, grouped
    # by source. Both remaining bucket *names* were false as well — "Essential" over columns the
    # block directly above says are none of them required, and "Functional prediction" over three
    # columns whose own descriptions call them an amino-acid change annotation and HGVS
    # nomenclature. And the heading was false by omission: `InterVar` and `RENOVO_Class` are
    # clinical annotations this app reads and neither was in it.
    #
    # So this is a deletion of the duplicated and the false, not a rewording, and what survives
    # is the five entries no other list on this tab carries. The two groups are **derived** — a
    # column is a filter input or it is not — so `CIViC_Evidence_Rating` cannot come to be
    # described as unread on a day some filter starts reading it.
    st.markdown(f"### {_OTHER_COLUMNS_HEADING}")
    st.markdown(
        "Beyond the per-arm annotations above, and the population-frequency columns listed "
        "further down:"
    )

    _filter_inputs = {
        column
        for source in (REQUIRED_INPUTS, PATHOGENIC_INPUTS)
        for columns in source.values()
        for column in columns
    }
    _other_columns = (
        "ClinVar_VCF_CLNSIG",
        "CIViC_Evidence_Rating",
        "AAChange.refGene",
        "cDNA_Change",
        "Protein_Change",
    )

    for label, wanted in (
        ("Read by both arms' filters", True),
        ("Shown, never filtered on", False),
    ):
        columns = [
            col
            for col in _other_columns
            if (col in _filter_inputs) is wanted and col in COLUMN_DESCRIPTIONS
        ]
        if not columns:
            continue
        with st.expander(f"📋 {label}"):
            for col in columns:
                st.markdown(f"• **{col}**: {get_column_description(col)}")

    # Missing column handling. The two column lists below are *rendered from the derivation*
    # rather than typed out: they are the same sets `filters/absent_columns.py` reads out of
    # the vendored source, and a help page that transcribed them would be exactly the drifted
    # copy this project keeps finding — `config/columns.py` records the last one, a
    # hand-maintained column list in this very file.
    #
    # The *prose* around them was the drifted copy, and issue #150 found it: two claims this
    # block was still making after #136 deleted both from the warnings themselves. It said a
    # MAF missing a filter input "would normally be refused outright" — MAFigate never refuses
    # on an absent column, so that could only have meant *in the pipeline*, the comparison
    # decision 2 of the map retired — and it called pathogenic retention "unconditional", which
    # over-claims: the rescue clears the filter criteria, but `apply_filters` masks it with the
    # app's own population-frequency term, exempted only by the narrower ClinVar-only
    # `pathogenic_exemption`. Derived column lists do not protect hand-written sentences beside
    # them, which is the lesson worth leaving here.
    st.markdown("### ⚠️ Handling Missing Columns")
    st.info(
        f"""
    **What happens if columns are missing?**

    No filter is ever silently switched off. MAFigate fills a missing filter input with a
    stand-in value so you still get a report, and tells you it did — which means **that report
    is not a complete result**.
    - **Core Required columns**: Analysis will fail with error message (Hugo_Symbol, Chromosome, etc.)
    - **Filter-input columns** — somatic: {_columns(REQUIRED_INPUTS["somatic"])}; germline:
      {_columns(REQUIRED_INPUTS["germline"])} — filled neutrally, with a warning naming the
      column and saying whether the fill **added** rows or **removed** them
    - **{_columns(_ESCALATED)}**: the same fill, plus an **escalated** warning — these also
      feed the pathogenic-retention rule, so their absence means variants are *missing*
      from your report rather than merely re-cut
    - **`CIViC_Evidence_Level`**: genuinely optional. MAFigate skips the CIViC clause when the
      column is absent, and your report is unaffected — this one costs you nothing
    - **`Hugo_Symbol`**: a gene-list restriction cannot be applied without it, so the
      restriction is dropped and you see *more* genes than you asked for, with a warning
    - **Display columns** — {_columns(_DISPLAY_ONLY)}: nothing filters on these, so they are
      in neither the table nor anything you download, and every other column is unaffected.
      `DP` in particular is *not* what the depth setting reads — that counts the tumour's own
      alt and ref reads. `Tumor_Seq_Allele2` costs one thing more: a note, and any annotation
      column you invent, is stored against the variant's gene, position **and both alleles**,
      so without it MAFigate cannot tell two variants at one position apart and offers you
      neither field rather than risk showing your writing against the wrong variant
    - **Population frequency columns**: frequency filtering is bypassed if no frequency data available

    **Population frequency columns**: MAFigate reads every one of these that your file
    carries — not one per pair, and it does not rank them:
{_frequency_groups()}

    **Tip**: Open your MAF file first to see which columns are available in your dataset.
    """
    )

    # Return to data loading button
    st.markdown("---")
    col_return1, col_return2 = st.columns([2, 1])
    with col_return1:
        st.success(
            "✅ **Ready to try again?** Choose your MAF file in the sidebar — the chooser "
            "is there on every page."
        )
    with col_return2:
        if st.button(
            "🔙 Back to Data Loading", key="return_to_data_loading", type="primary"
        ):
            st.session_state.current_page = "data_loading"
            st.rerun()


def show_filter_options_tab():
    """Show all available filter options and their meanings."""

    st.subheader("🔍 Filter Options")
    st.markdown(
        """
    MAFigate provides multiple filtering options to refine your variant analysis.
    Understanding these options helps you configure appropriate filters for your research.

    The clinical annotation filters below are **keep** lists combined with OR: a variant
    meets the criteria if any one source keeps it. Leaving a source empty means it places
    no restriction — it drops out of the comparison rather than widening it. Emptying
    *every* source leaves only the variants pathogenic retention rescues on their own,
    and the parameter page warns you when you have done so.
    """
    )

    # Variant classification filters
    st.markdown("### 🧬 Variant Classification Filters")
    st.markdown(
        "This control **excludes** by functional impact: a classification you select is "
        "dropped, and everything else is kept — including classifications this list does "
        "not name. MAFigate starts with `Silent`, `IGR` and `RNA` excluded."
    )

    with st.expander("📋 Classifications you can exclude"):
        for i, variant_class in enumerate(VARIANT_CLASSIFICATIONS, 1):
            st.markdown(f"{i}. **{variant_class}**")

    # Clinical annotation filters. One word changed from "database" (issue #219), on the
    # sentence above and on this comment: three of the six sources these filters read are not
    # databases — `RENOVO` is a predictor this pipeline ran (issue #201), `InterVar` and
    # `CancerVar` are guideline verdicts computed here (issue #187).
    #
    # Map #199 rule 5 leaves filtering behaviour and its documentation alone, and this respects
    # it: nothing here changes what any filter reads, which sources it offers or how their
    # combination is described. What it removes is a category *word* the map falsified — and
    # leaving it would have had the FAQ two tabs over carefully avoiding the phrase this tab
    # leads with, which is the shape of self-contradiction issue #219 exists to end.
    col1, col2 = st.columns(2)

    with col1:
        st.markdown("### 🏥 Somatic Filters")

        with st.expander("CancerVar Classifications"):
            st.markdown("**AI-based cancer variant interpretation:**")
            for option in CANCERVAR_OPTIONS:
                st.markdown(f"• **{option}**: Cancer-specific pathogenicity tier")

        # One gloss per level, read out of `CIVIC_DEFINITIONS` and cited — the ESCAT block
        # below's shape, for the same reason (issue #109). What stood here was the last
        # hand-written vocabulary on this tab, and #79 let it stand as roughly right while
        # nothing else read the scale. Something else does now: the Pathogenicity Overview
        # grades its CI position on these levels, and the gloss it was checked against called
        # `A` "Strong clinical significance" — which names the *other* axis, the one CIViC's
        # levels are not on, and is how a level came to be drawn as a pathogenicity call.
        with st.expander("CIViC Evidence Levels"):
            st.markdown(
                "**How strong the CIViC evidence behind an assertion about this variant is.** "
                "Strongest first, and the scale grades the evidence, not the variant:"
            )
            for level, meaning in CIVIC_DEFINITIONS.items():
                st.markdown(f"• **Level {level} — {meaning.name}**: {meaning.evidence}")
            st.caption(f"Source: {CIVIC_SOURCE} {CIVIC_SOURCE_URL}")

    with col2:
        st.markdown("### 🧪 Germline Filters")

        with st.expander("InterVar Classifications"):
            st.markdown("**ACMG/AMP guidelines-based interpretation:**")
            for option in INTERVAR_OPTIONS:
                st.markdown(f"• **{option}**: ACMG/AMP classification")

        with st.expander("RENOVO Classifications"):
            st.markdown("**Germline pathogenicity assessment:**")
            for option in RENOVO_OPTIONS:
                st.markdown(f"• **{option}**: Pathogenicity prediction")

    # Shared filters
    st.markdown("### 🔗 Shared Filters (Both Sample Types)")

    # Every term the control offers, not the five this expander used to explain. The old
    # version walked a hand-written dict of five meanings and gated each on membership of
    # `CLINVAR_OPTIONS` — so it presented five of the values a user can actually pick as
    # though they were the list, and said nothing about the rest, among them
    # `Conflicting_classifications_of_pathogenicity`, `drug_response` and `not_provided`.
    # Walking the vocabulary and annotating what we can state plainly reverses that: a term
    # without a gloss is still listed as selectable.
    #
    # The counts that used to sit in this note ("five of the eleven", "the other six") are
    # gone rather than corrected. They were true when written and false two tickets later,
    # #88 having taken the two composite terms out of `CLINVAR_OPTIONS` — and a number
    # restated beside the list it counts is exactly the drift this page was rewritten to
    # stop. The list is walked; nothing here needs to know how long it is.
    #
    # The captions carry the one thing walking the vocabulary cannot say for itself: some of
    # these entries are pairs, two spellings of one ClinVar call — the 2023 rename (#88) and
    # the underscore-prefixed modifiers (#99). Said in a caption rather than as glosses,
    # because a gloss here is a *clinical* claim and this page has already been wrong that way
    # once — the ESCAT levels it used to describe from memory (#79). A spelling fact is not a
    # clinical claim, which is why it can be stated plainly.
    #
    # The underscore caption reads its members *out of the vocabulary* and renders nothing when
    # there are none, for the reason the note above gives about counts: a membership list
    # restated beside the list it describes is the drift this page was rewritten to stop, and
    # #103 may well add terms. The rename caption names its two by hand because they are one
    # pair whose names are the recognisable thing about them, and no rule derives them.
    #
    # #103 did add terms, and it grew the underscore caption's membership without a line
    # changing here, which is the whole point of that derivation. The first draft of this note
    # said how many, twice, and got both numbers wrong — beside a paragraph that had just
    # warned against restating a count. They are gone rather than corrected, for that reason.
    #
    # What #103 *did* need was the split into two sections, because the page walks the same
    # vocabulary the controls do and there are two controls now: the pathogenicity calls and
    # the assertions that are not pathogenicity calls, #98's class boundary. The section
    # captions describe that boundary and nothing beneath it — which class each term belongs
    # to is `components/clinical_summary.py`'s, and drawing a per-term class here would be a
    # second copy of it, the drift this note has already warned about twice.
    #
    # The clinical glosses stay at five. The terms #103 added are the ones this repository has
    # least authority over — the institute's table classes them, which is what the app says
    # about them, and it is said where the class is rendered rather than restated here as
    # prose. Three record-state terms are glossed on top of the five, from quotes this repo
    # already holds; the dict below says which and why.
    with st.expander("ClinVar Classifications"):
        # Two kinds of gloss, and the difference is where they come from.
        #
        # The first five are the ACMG/AMP calls, glossed since before this map. The three
        # after them are glossed by issue #103 for a different reason: their *names* describe
        # the state of a ClinVar record rather than a finding about the variant, so a reader
        # cannot infer them the way `Likely_benign` can be inferred. `-` is the extreme case —
        # a bare dash in a list, meaningless on sight.
        #
        # Each of the three is quoted from a source this repository holds: the institute's
        # term table and ClinVar's own definition, both cited in
        # `components/clinical_summary.py` where #98 recorded them. That is the whole reason
        # these three and no more. `no_classifications_from_unflagged_records` stays unglossed
        # because nothing here defines it, and #79 deleted invented clinical prose rather than
        # write it from memory — the same judgement, made again rather than quietly dropped.
        clinvar_meanings = {
            "Pathogenic": "Established pathogenic variant",
            "Likely_pathogenic": "Likely to be pathogenic",
            "Uncertain_significance": "Unknown clinical significance",
            "Likely_benign": "Likely to be benign",
            "Benign": "Established benign variant",
            "other": (
                "ClinVar has a record, but none of its terms fitted the submission — what it "
                "means is in a free-text explanation your MAF does not carry"
            ),
            "-": (
                "Not submitted as a classification: the variant was recorded only as part of "
                "a haplotype, so ClinVar has no call for it on its own"
            ),
            "no_classification_for_the_single_variant": (
                "The same as the dash above, spelled the way an annotated MAF writes it"
            ),
        }

        def _draw(terms):
            for classification in terms:
                gloss = clinvar_meanings.get(classification)
                st.markdown(
                    f"• **{classification}**: {gloss}"
                    if gloss
                    else f"• **{classification}**"
                )

        # Two sections, in the order the two controls are drawn, so a reader who came here
        # from the parameter page finds this list in the shape they just left. The heading
        # words are the controls' own.
        st.markdown("**Clinical significance from ClinVar:**")
        _draw(CLINVAR_PATHOGENICITY_TERMS)
        st.markdown("**Other ClinVar assertions:**")
        _draw(CLINVAR_OTHER_ASSERTION_TERMS)
        st.caption(
            "The second group are not pathogenicity calls — they are what ClinVar records "
            "about a variant when a pathogenicity call is not what it has: a disease risk, a "
            "drug response, an association with a trait, a protective effect, or a record "
            "holding no classification for this variant on its own. Selecting one keeps a "
            "variant on that assertion alone, so a variant with no pathogenicity call can "
            "still reach your report."
        )
        st.caption(
            "A ClinVar entry often carries several of these at once. Your selection is "
            "matched against the individual terms of such a value, not against the whole "
            "string."
        )
        st.caption(
            "**Conflicting classifications** and **conflicting interpretations** are one "
            "ClinVar call under two names: it was renamed in 2023, so which one appears in "
            "your file depends on when the file was annotated. Keep both selected unless "
            "you know which release yours was annotated against."
        )
        underscored = [term for term in CLINVAR_OPTIONS if term.startswith("_")]
        if underscored:
            named = ", ".join(f"**{term}**" for term in underscored)
            example = underscored[0]
            st.caption(
                f"The terms beginning with an underscore — {named} — are not separate "
                "classifications. A ClinVar entry carrying a call *and* a modifier is written "
                "with the two joined together, and files differ over whether the modifier "
                f"keeps a leading underscore, so one MAF spells it `{example.lstrip('_')}` and "
                f"another `{example}`. Keeping both of each pair is what makes your selection "
                "mean the same thing whichever way your file was written; keeping one can miss "
                "every such variant in the file."
            )

    # One gloss per level, read out of `ESCAT_DEFINITIONS` and cited. What stood here was a
    # hand-written meaning per level, wrong in at least two of eight: `IIIA` and `IIIB` were
    # described as "Resistance mechanism with clinical/preclinical evidence", and ESCAT has
    # no resistance tier at all — resistance levels are OncoKB's R1/R2, which this app does
    # not read. `V` was "Not actionable", which the app contradicted twice over: the
    # parameter page's own tooltip glossed the same scale as "IA=strongest evidence to
    # V=case reports", and the Pathogenicity Overview maps tier III to Uncertain
    # Significance and V to Unknown. Three surfaces, three stories, and nothing in this
    # repository was an authority on ESCAT to settle it — so #79 deleted the invented prose
    # rather than rewrite it from memory, and #89 put back a sourced version. The scale's
    # identity and direction, the half #79 could state without a source, survive above.
    #
    # The definitions live in `config/vocabularies.py`, not here, because two surfaces need
    # them and a page is where the copy that drifted lived. This one renders them; the
    # parameter tooltip names the same constant for the direction and glosses nothing.
    with st.expander("ESCAT Evidence Levels"):
        st.markdown(
            "**ESCAT** is the ESMO scale for how actionable a target is, running from "
            f"`{ESCAT_STRONGEST}` — the strongest evidence that acting on it benefits the "
            "patient — downwards."
        )
        group_drawn = None
        for level, meaning in ESCAT_DEFINITIONS.items():
            if meaning.group != group_drawn:
                st.markdown(f"**{meaning.group}**")
                group_drawn = meaning.group
            st.markdown(f"• **`{level}`** — {meaning.evidence} *{meaning.implication}*")
        st.caption(
            f"Definitions from {ESCAT_SOURCE} {ESCAT_SOURCE_DOI}. ESCAT also defines "
            "`IV` (preclinical evidence only) and `X` (no evidence of actionability); "
            "neither is listed above, because neither is a level this filter offers."
        )

    # Population frequency filters
    st.markdown("### 📊 Population Frequency Filters")

    with st.expander("🧬 gnomAD Database", expanded=True):
        st.markdown(
            """
        **gnomAD (Genome Aggregation Database)** is the largest population reference this
        filter reads:

        **gnomAD columns MAFigate reads:**
        - **`gnomAD_exome_AF`** - Allele frequency from exome sequencing, over the genotypes
          that passed gnomAD's own quality filters
        - **`gnomAD_genome_AF`** - The same, from genome sequencing
        - **`gnomAD_exome_AF_raw`** - Exome allele frequency before those quality filters, so
          over every genotype gnomAD called
        - **`gnomAD_genome_AF_raw`** - The same, from genome sequencing

        All four are read the same way, and that is a decision rather than an oversight: the
        rule below takes the **lowest** frequency any column reports, so a pre-QC value can
        only ever keep a variant that the adjusted value would also have kept, and never drop
        one the adjusted value calls rare. Where your file carries only the `_raw` column —
        which is what a MAF that has already been through the VarianTalker pipeline looks
        like — it is the only gnomAD frequency there is to read.

        **How MAFigate uses gnomAD:**
        - Checks every gnomAD column your file carries
        - Keeps a variant if ANY database that has a value for it shows frequency ≤ threshold
        - A blank column is *not seen in this panel*, not *rare*: it never counts against a
          variant, and it does not vouch for one either, so a variant every populated
          column calls common is dropped
        - A variant no database has a frequency for is kept
        - **Exempts ClinVar-pathogenic variants**, so a common low-penetrance pathogenic
          allele is never dropped for its frequency alone
        - A threshold of 1.0 switches this filter off entirely, leaving every variant to be
          judged on the other filters alone

        **Database Coverage:**
        - **Exomes**: ~125,748 individuals, better for coding regions
        - **Genomes**: ~15,708 individuals, better for non-coding regions
        - **Populations**: African, American, Ashkenazi Jewish, East Asian, European, South Asian
        """
        )

    with st.expander("📚 Legacy Population Databases"):
        st.markdown(
            """
        **Also supported for backward compatibility:**
        - **ExAC** (`Freq_ExAC_ALL`) - Predecessor to gnomAD, now integrated
        - **1000 Genomes** (`Freq_1000g2015aug_all`) - International reference panel
        - **ESP6500** (`Freq_esp6500siv2_all`) - NHLBI Exome Sequencing Project
        
        **Note**: These are older panels with smaller sample sizes, so a variant they have
        never seen may still be common. MAFigate does not rank them below gnomAD — it reads
        every frequency column your file carries and takes the lowest value any of them
        reports.
        """
        )

    st.info(
        """
    **💡 Frequency Filtering Tips:**
    - **0.01 (1%)**: Removes common variants, keeps rare ones — what the Broad presets use
    - **0.005 (0.5%)**: Tighter, for a report you would act on — what the Stringent presets use
    - **0.001 (0.1%)**: Very stringent filtering for extremely rare variants
    - **1.0**: Switches frequency filtering off entirely
    - **Missing data**: a variant no database has a frequency for is kept; a blank in one
      database does not rescue a variant the others call common
    """
    )

    # `st.markdown("### 🔍 Other Filter Categories")` closed this tab: a heading with nothing
    # whatsoever under it, rendering as an empty section at the foot of the page. Whatever it
    # was going to introduce was never written.


def show_user_guide_tab():
    """Show comprehensive user guide."""

    st.subheader("📖 User Guide")
    st.markdown("Step-by-step guide to using MAFigate effectively.")

    # Workflow overview
    st.markdown("### 🔄 Analysis Workflow")
    st.markdown(
        """
    Follow this recommended workflow for best results:
    
    1. **📋 Prepare your data** - Ensure MAF file has required columns
    2. **⚙️ Configure parameters** - Set up filtering criteria
    3. **📁 Load data** - Choose your MAF file in the sidebar
    4. **🔍 Apply filters** - Run analysis with your parameters
    5. **📊 Review results** - Examine passed and failed variants
    6. **💾 Export data** - Download results for further analysis
    """
    )

    # Detailed steps
    with st.expander("1️⃣ Data Preparation"):
        st.markdown(
            """
        **Before you open your MAF file:**
        
        ✅ **Verify file format**: Tab-separated values (.maf, .tsv, .txt)
        ✅ **Check required columns**: See the "Required Columns" tab
        ✅ **Validate data quality**: Ensure no major formatting issues
        ✅ **File size**: up to 200 MB through the browser

        **Common MAF sources:**
        - TCGA/GDC data portal
        - cBioPortal downloads
        - Local variant calling pipelines
        - Clinical sequencing reports
        """
        )

    with st.expander("2️⃣ Parameter Configuration"):
        st.markdown(
            """
        **Parameter Configuration Options:**
        
        🎯 **Preset Parameters**:
        - **Broad Somatic/Germline**: a wide net — keeps uncertain calls for review
        - **Stringent Somatic/Germline**: a short list — only what you would act on
        
        🔧 **Custom Parameters**:
        - Adjust individual filter thresholds
        - Choose a gene set from the **Gene Set** list, or paste your own
        - Modify population frequency cutoffs
        - Configure quality metrics

        🧬 **Gene Set Options**:
        - **All**: No gene filtering (analyze all genes)
        - **MSK_IMPACT / MSK_HEME**: MSK-IMPACT solid-tumour / heme panel genes
        - **FOUNDATION_ONE / FOUNDATION_HEME**: FoundationOne CDx / Heme panel genes
        - **COSMIC**: COSMIC Cancer Gene Census genes
        - **OncoKB**: OncoKB actionable cancer genes
        - **Custom**: Paste or upload your own gene list — one symbol per line, or
          separated by commas, semicolons or spaces. A pasted column heading is
          dropped, matching ignores case, and any entry that cannot be a gene symbol
          is listed back to you rather than silently applied.

        💡 **Tips**:
        - Start with preset parameters
        - Use a targeted panel gene set for focused clinical analysis
        - Use OncoKB for comprehensive actionable variants
        - Test with sample data first
        """
        )

    with st.expander("3️⃣ Data Loading"):
        st.markdown(
            """
        **Loading Your Data:**
        
        📁 **The file chooser is in the sidebar**, so a sample can be opened from any
        page:
        - Drag and drop, or click to browse
        - Tab-separated `.maf`, `.tsv` or `.txt`, UTF-8, with any `#` comment lines
          handled for you

        **After loading, you'll see**:
        - Basic statistics (variants, samples, genes)
        - Column information table
        - Data preview
        """
        )

    # This expander promised a feature the app has never had. "**Quick Parameter
    # Adjustment**: Modify key filters without leaving the page / Real-time parameter
    # validation / Immediate feedback on changes" describes no control that exists: the Load
    # section offers **Re-apply Filters** and a button that *leaves* the page for Configure
    # Parameters, and there is nowhere else a threshold can be touched. "Compare different
    # filtering strategies" named a comparison feature that does not exist either. What is
    # left is what the two buttons actually do.
    with st.expander("4️⃣ Filter Application"):
        st.markdown(
            """
        **Applying Filters:**

        Filters run on your file as soon as it opens, with whatever settings you have.

        🔄 **Changed a setting since?** **Re-apply Filters**, on the Load section, runs
        them again — the page warns you when your settings and your report have drifted
        apart.

        ⚙️ **To change a setting**: **Configure Parameters**, from the sidebar or the
        button beside Re-apply Filters. Your settings are kept between sessions, and can
        be downloaded as JSON or YAML and uploaded again.

        **Filter execution provides**:
        - Passed variants (meet all criteria)
        - Failed variants (with failure reasons)
        - Processing warnings and statistics
        """
        )

    with st.expander("5️⃣ Results Analysis"):
        st.markdown(
            """
        **Understanding Your Results:**
        
        📊 **Summary Statistics**:
        - Pass/fail rates
        - Top genes affected
        - Filter performance metrics
        
        📋 **Interactive Tables**:
        - Sort by any column, and filter in the box under each heading — number columns
          take a comparison, such as `>50` or `<0.01`
        - Add columns with **Add columns**, or show every one your file carries
        - Click a variant for the full record, and to write a note

        🔍 **Quality Assessment**:
        - Review failed variants for patterns
        - Adjust filters if needed
        - Validate biological relevance
        """
        )

    # "**TSV**: Tab-separated for downstream tools" was the same false claim #55 took out of
    # both READMEs and #71 out of the About dialog: every download a user can reach writes
    # CSV, plus the summary as text. The one exporter that offered TSV and Excel had no
    # caller anywhere in the app, and #83 has since deleted it outright.
    #
    # Two of these bullets say what they say because #83 decided it. The narrow CSV is **the
    # report's own columns**, deliberately *not* "the columns on screen": the grid carries a
    # *Show all columns* checkbox and an *Add columns* multiselect, so what is on screen is
    # the user's to change while the file is always the resolver's list — "shown" was the
    # wording #83 removed for being a lie, and this page must not put it back. And the
    # filenames lead with the loaded MAF's name rather than a timestamp, because four files
    # called `passed_variants_all.csv` are four files a clinician cannot tell apart.
    with st.expander("6️⃣ Data Export"):
        st.markdown(
            """
        **Exporting Results:**

        💾 **What you get**:
        - **CSV** for the variants — opens in Excel, R and Python
        - **Summary** as a text file: how the report was reached, your settings, and the
          genes it landed on

        📁 **Export Options**:
        - Passed variants, or failed variants with the criterion that decided them
        - The report's own columns, or every column your file carries — each button says
          how many
        - Filenames lead with the name of the MAF you loaded

        **Best Practices**:
        - Keep the summary alongside the table — it records the settings the report used
        - Export both passed and failed for review
        - Notes and any columns you added live for the session only, so download the
          table to keep them
        """
        )

    # Best practices
    st.markdown("### 💡 Best Practices")
    st.success(
        """
    **🎯 For Optimal Results:**
    - Always review failed variants to understand filter impact
    - Start with inclusive filters, then make them more stringent
    - Check column availability before configuring complex filters
    - Export both passed and failed variants for complete analysis
    
    **⚠️ Common Issues:**
    - Missing clinical annotation columns → filled with a stand-in value and reported; the
      report is not a complete result, and for `CancerVar`, `ClinVar_VCF_CLNSIG` or
      `InterVar` it is also **smaller** than it should be
    - Very large files → consider pre-filtering, or narrowing to a gene panel
    - Unusual MAF formats → Check column names match standard format
    """
    )


def show_faq_tab():
    """Show frequently asked questions."""

    st.subheader("❓ Frequently Asked Questions")

    # File format questions
    with st.expander("📁 File Format & Data Questions"):
        st.markdown(
            """
        **Q: What file formats does MAFigate accept?**
        A: MAFigate accepts standard MAF (Mutation Annotation Format) files with extensions .maf, .tsv, or .txt. Files should be tab-separated.
        
        **Q: My MAF file has different column names. Will it work?**
        A: MAFigate expects standard MAF column names. If your file uses different names, you may need to rename columns to match the standard format.
        
        **Q: How large can my MAF file be?**
        A: Up to 200 MB, which is the limit on files arriving through the browser. For a larger MAF, ask whoever set MAFigate up to raise it. Bear in mind that the whole file is held in memory while you work on it.

        **Q: Can I use VCF files instead of MAF?**
        A: No, MAFigate specifically works with MAF format files. You'll need to convert VCF to MAF using tools like Funcotator or VEP.
        """
        )

    with st.expander("🔍 Filtering & Analysis Questions"):
        st.markdown(
            """
        **Q: What's the difference between the Broad and Stringent presets?**
        A: Broad casts a wide net — it keeps uncertain and conflicting calls so you can review them. Stringent keeps only well-supported pathogenic calls, for a report you would act on.
        
        **Q: Why are some variants failing filters?**
        A: Variants fail filters if they don't meet the configured criteria (e.g., low VAF, high population frequency, benign classification).
        
        **Q: Can I save my filter configurations?**
        A: They are saved already — MAFigate reopens with the settings you last used, and says so at the top of the page. To keep a particular set, or move it to another machine, use **Download JSON** or **Download YAML** on the Configure Parameters page; **Upload Custom Preset**, on the same page, reads one back.

        **Q: What happens if clinical annotation columns are missing?**
        A: No filter is ever silently switched off. MAFigate fills the missing column with a stand-in value, filters your file anyway, and tells you it did — which means the report is not a complete result, and for some columns it is smaller than it should be. The **Required Columns** tab sets out which columns those are and what each one's absence costs.
        """
        )

    # The interpretation answer was wrong twice over, and issue #219 replaced it rather than
    # trimming it. It said "refer to the Column Information tab for detailed descriptions of each
    # clinical **database** and their classification systems":
    #
    # * *Database* is false of most of them. RENOVO is an in-silico predictor this pipeline ran
    #   (issue #201) — passing *computed here* is exactly the leg a third-party database fails —
    #   `InterVar` and `CancerVar` are guideline verdicts computed here (issue #187), and
    #   AlphaMissense is a predictor (issue #203). ClinVar and CIViC are the only two the word
    #   fitted. The replacement claims no category at all, which is the one framing that cannot
    #   be false of a member: `MAFigate.py` reached the same wording independently when it
    #   deleted this claim from the About dialog.
    # * The Column Information tab holds **no** classification systems. It is one line per
    #   column, from `COLUMN_DESCRIPTIONS`; every scale — ClinVar's terms, ESCAT's levels,
    #   CIViC's evidence levels, CancerVar, InterVar and RENOVO's classes — is rendered by
    #   `show_filter_options_tab`. So the one FAQ that asks the panel's question sent the reader
    #   to the wrong tab, and had done since before this map.
    #
    # The panel is named because it is the surface that answers the question, and this page never
    # mentioned it — no badge, no hover, no panel anywhere in the file. One sentence, not a
    # description: a second copy of the panel's vocabulary here is the drift this map keeps
    # finding, and issue #217 already decided the row explains itself through its badges rather
    # than through a heading. The route is the two-step one the app actually offers (issue #160):
    # select a row, then press the button, which opens a dialog.
    with st.expander("📊 Results & Export Questions"):
        st.markdown(
            """
        **Q: Why do I see both passed and failed variants?**
        A: This design helps you understand filter impact and review variants that didn't meet criteria, which might still be biologically relevant.
        
        **Q: How do I interpret the clinical annotations?**
        A: Select a variant in the results table and press **View details**. The panel draws each annotation once and carries the scale it is on. Each scale is set out in full under **Filter Options**; the **Column Information** tab says what each column holds.
        
        **Q: Can I re-filter my results without re-uploading?**
        A: Yes! Once data is loaded, you can modify parameters and re-apply filters without uploading the file again.
        
        **Q: What's included in the summary report?**
        A: The summary includes filter parameters, pass/fail statistics, top genes, and analysis metadata.
        """
        )

    # The page's most consequential false claim, on an app whose input is patient data.
    # "MAFigate runs in your browser. Data is processed locally" describes something this is
    # not: MAFigate is a server application that you *reach* through a browser. Your file is
    # uploaded to the Python process, which writes it to a temporary file to read it
    # (`utils._spilled_to_disk`, unlinked after the read) and then holds the frame in memory
    # for the session. Whether "locally" is true depends entirely on where that process runs,
    # which this app cannot know — the launchers bind loopback (issue #182), but nothing stops
    # a colleague serving it from a shared machine, and the answer below is written for that
    # reader. The true half of the old answer is the half worth keeping: nothing is sent to any
    # outside service, because the app queries nothing.
    with st.expander("🔧 Technical Questions"):
        st.markdown(
            """
        **Q: Where does my data go?**
        A: Nowhere outside MAFigate. It reads the annotations already in your file and queries no external service, and your notes are never written to disk.

        Your file is uploaded to the computer running MAFigate, which holds it for as long as your session. If that is your own machine, your data does not leave it. If a colleague runs MAFigate on a shared server for you, your file goes to that server — ask them how it is set up before opening patient data.

        **Q: Why is the application slow with large files?**
        A: Large files require more processing time. Consider pre-filtering your data or using a more powerful computer.
        
        **Q: Can I use MAFigate for non-human data?**
        A: MAFigate is designed for human genomic data. Non-human data may work but clinical annotations won't be relevant.
        
        **Q: Are there any browser requirements?**
        A: MAFigate works best in modern browsers (Chrome, Firefox, Safari, Edge). Ensure JavaScript is enabled.
        """
        )

    with st.expander("📚 Resources & Support"):
        st.markdown(
            """
        **Q: Where can I learn more about MAF format?**
        A: See the [GDC MAF documentation](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/) for detailed MAF format specifications.
        
        **Q: How do I cite MAFigate in publications?**
        A: Give the application name, the version you ran — the ⋮ menu's *About* entry has it — and your institution.
        
        **Q: Where can I report bugs or request features?**
        A: Contact the development team through your institutional channels or create an issue in the project repository.
        
        **Q: Is there documentation for developers?**
        A: Yes, see the project README and inline code documentation for development details.
        """
        )

    # Quick help summary
    st.markdown("### 🚀 Quick Start Checklist")
    st.info(
        """
    **New User? Follow this checklist:**

    ☐ Review the Column Information tab to understand available annotations
    ☐ Try different preset parameters (Broad vs Stringent)
    ☐ Open your own MAF file and check which columns are available
    ☐ Configure filters appropriate for your research question
    ☐ Review both passed and failed variants to understand filter impact
    ☐ Export results with summary report for documentation
    
    **Need help?** Most interface elements have helpful tooltips when you hover over them!
    """
    )
