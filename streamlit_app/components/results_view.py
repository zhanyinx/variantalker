"""The Results section: the passed and failed tabs, and the summary beside them.

The orchestrator of what a user sees after a filter run. It draws the three tabs,
hands each frame to :func:`~components.variant_table.create_data_table`, offers the
per-tab downloads, and calls the summary dashboard. It decides layout and nothing
else — the grid, the charts and the clinical derivation each live in their own module
and are called from here.

A page section rather than a page: ``page_modules/data_loading.py`` calls
:func:`create_enhanced_data_table` and owns the surrounding page.

``create_export_buttons`` sat here with no caller, left standing by issue #76 so that
issue #83 could decide it beside the four downloads it would replace. #83 deleted it:
it serialised the frame it was handed and never called
:func:`~components.variant_table.with_user_columns`, so adopting it as the export
surface would have shipped the bug issue #67 had just fixed — a download that drops
the notes the user typed — and it offered no column choice either. Results download as
CSV alone; ``openpyxl`` left ``requirements.txt`` with it, and ``convert_df``, its only
reader, left ``utils/main_utils.py``.

Every download this section offers is drawn under the thing it exports: the CSVs
beneath each variant table, the run report inside the Summary tab. The ``💾 Download
Results`` block that used to follow the tabs is gone — its two variant buttons served
the same bytes as the ``all columns`` buttons above them, differing only in the name
the file landed under.

Two CSVs beneath a table, or one. Since issue #92 the narrow download delivers the
columns the grid is showing rather than a fixed list, so when the user ticks *Show all
columns* the two coincide and only the ``all columns`` button is drawn — the same rule
that removed the ``💾 Download Results`` block, applied to a state rather than a block.

The Summary tab's ``Column Analysis`` block went at issue #90, and the reason is worth
keeping: it answered *what is my file missing?* against
:data:`~config.columns.PRIORITY_COLUMNS`, which is a presentation **order** — one flat
list of somatic-flavoured names — and not :func:`~config.columns.resolve_visible_columns`,
which asks the vendored ``compute_keep`` what this arm emits for this file. So it reported
absence where the pipeline had made a **choice**: a MAF carrying no CIViC annotation is not
incomplete, ``compute_keep`` drops those columns from its own output list, and the block
drew six red ``❌`` at it regardless. On the germline arm it drew nine, ``CancerVar``,
``cosmic``, ``Freq_ExAC_ALL`` and the CIViC block being columns that arm never emits.
**On both of this repo's reference MAFs the app's verdict is that nothing is missing, and
every marker the block drew there was a false alarm.** Its denominator was rigged the other
way: ``Clinical_Summary`` heads that list and is derived here at filter time, so it counted
as present whatever the file said. Nothing replaced it — the resolver reports a real
absence by name and per arm, ``filters/absent_columns.py`` escalates when a neutral fill
would move the verdict, and a fresh line here would prejudge issue #93, which is deciding
how many times that one sentence should be said. ``tests/test_results_summary.py`` holds
all of this.
"""

from pathlib import PurePosixPath

import streamlit as st
import pandas as pd

from config.columns import prioritize_columns
from config.constants import TABLE_HEIGHT

from .charts import render_summary_dashboard
from .clinical_summary import format_clinical_summary_display
from .variant_table import (
    create_data_table,
    render_fallback_table,
    shown_columns,
    with_user_columns,
)


#: Extensions a MAF arrives under, stripped from the download prefix longest-first so
#: ``.maf.gz`` comes off whole. A closed list rather than repeated ``Path.stem`` calls,
#: which would eat a sample named ``TCGA.02.0001`` down to ``TCGA``.
_MAF_SUFFIXES = (".maf.gz", ".txt.gz", ".tsv.gz", ".maf", ".txt", ".tsv", ".csv", ".gz")


def _download_prefix() -> str:
    """The loaded file's name, stripped of its extension, to lead every download's name.

    A clinician works several samples in a sitting, and four files all called
    ``passed_variants_all.csv`` are four files they cannot tell apart in a Downloads
    folder — the browser resolves that with ``(1)``, which names nothing. The app has
    recorded the source file's name since issue #58 for the sidebar's status block;
    nothing else had thought to read it.
    """
    name = PurePosixPath(str(st.session_state.get("maf_source_name") or "")).name
    for suffix in _MAF_SUFFIXES:
        if name.lower().endswith(suffix):
            return name[: -len(suffix)]
    # No file name recorded — a defensive branch, since every load path stamps one.
    return name or "variants"


def _distinct_count(frames, column: str) -> "int | str":
    """How many distinct values ``column`` holds across ``frames``, or ``N/A`` if none has it.

    ``N/A`` rather than ``0``, which is what the load page's Data Overview printed before
    issue #65 moved these counts here. Zero is a claim about the file — that it contains no
    samples, no genes — and a MAF missing ``Tumor_Sample_Barcode`` still came off somebody's
    sequencer. What is true is that this file does not say.

    Several frames rather than one, and the column is taken *before* the concat. The callers
    hold the file as two frames, and joining them whole to count one column would rebuild
    the 427-column concatenation issue #90 deleted from this tab — over a frame this counts
    two columns of.
    """
    present = [frame[column] for frame in frames if column in frame.columns]
    if not present:
        return "N/A"
    return pd.concat(present).nunique()


def render_file_facts(*frames: pd.DataFrame) -> None:
    """Samples and Genes for a whole MAF, as two metric tiles.

    The two facts the load page's Data Overview held that no other surface reports, rehoused
    when that block was deleted (issue #65). Its third tile, Total Variants, is not here: the
    sidebar names it on every page, which is what made it a duplicate rather than a loss.

    One definition, because two states need it and they are reached by different routes — the
    summary tab, where the report exists and the file arrives as its passed and failed halves,
    and the no-report branch of ``show_results_section``, where the filters failed and the
    whole MAF is the only frame there is. That second call site is the reason this is a
    function: without it these counts would exist only on a tab a refused MAF cannot reach.

    A four-column row with the last two left empty, so the tiles line up under the summary's
    Total Variants and Passed Filters rather than stretching to half the page each.
    """
    samples_col, genes_col, _, _ = st.columns(4)

    with samples_col:
        st.metric("Samples", _distinct_count(frames, "Tumor_Sample_Barcode"))

    with genes_col:
        st.metric("Genes", _distinct_count(frames, "Hugo_Symbol"))


def render_filter_recap(
    recap, stale: bool = False, on_reapply=None, attribution=None
) -> None:
    """What the filters were, above the report they produced, and the way back to them.

    The results view said *"passed all applied filters"* in three places and named none of
    them (issue #137); to find out what had run, a user had to leave the report for
    Configure Parameters and read it off the controls — which show the settings as they are
    **now**, not the ones the report was built from.

    **A recap and a route, never a second set of controls.** Configure Parameters keeps
    sole ownership of the widgets. Duplicating even one of them here would put two writers
    on ``filter_params`` a scroll apart, and this surface exists because the two can already
    disagree.

    **Every setting live on this arm, including the ones that cut nothing** (the dev's
    call). *"Genes: every gene — no gene restriction"* is a line about a filter that did
    not fire, and it is the most useful line here: the mismatch this surface is for —
    *"if this is not in line with what the user thinks"* — is most often a gene panel the
    user believes they set. Listing only what fired would hide exactly that.

    **An expander, not a fourth tab and not always-on.** The section already carries three
    tabs under a segmented control, so a fourth tab would be a fourth thing to choose
    between before reading a report; collapsed, this costs one line. The way-back button is
    *inside* it rather than beside its title because Streamlit's expander has no slot in its
    header — so the expander is the "show me" press and the button within it the "change
    it" press, which is the two-step the request describes.

    ``stale`` is the state ``data_params_hash`` exists to detect: the parameters have moved
    since this report was produced. The recap is then the only surface that says *what*
    changed is not describable from here — so it says plainly which of the two it is
    describing, and offers the re-run rather than leaving the user to find it on another
    section. The sidebar and the load section already warn that the gap exists; this one is
    beside the affected table and names the cut it belongs to.

    ``attribution`` is issue #147's half: what the cut *did*, under what it *was*. It arrives
    the same way and on the same clock — a :class:`~filters.attribution.ReportAttribution`
    snapshotted by the filter run — because a count re-derived here would describe the
    controls as they now stand and could contradict the settings listed directly above it,
    which is the one outcome that would be worse than saying nothing. The dev's call was
    this expander rather than the Summary tab, so the two halves cannot come apart.

    The counts **overlap** and the block says so. That is the framing chosen over the
    partitioning one, and the reason is measured: on real MAFs 70–95% of excluded rows fall
    outside more than one setting, so "the one setting responsible" puts almost everything in
    a single bucket, while these numbers name the two settings actually doing the work. See
    :class:`~filters.attribution.ReportAttribution`.
    """
    if recap is None:
        return

    with st.expander("⚙️ Filters applied to this report"):
        if stale:
            st.warning(
                "⚠️ Your settings have changed since this report was produced. "
                "What follows is what produced it, not what is set now."
            )
            if on_reapply is not None and st.button(
                # Word for word the load section's button, because it is the same act.
                "🔍 Re-apply Filters",
                key="recap_reapply",
            ):
                on_reapply()

        # The arm leads, and is not one line among the others: it decides which of them
        # exist at all — the somatic guideline sources are not arguments the germline
        # filter takes.
        st.markdown(f"**{recap.arm.title()} analysis**")

        if recap.note:
            st.warning(recap.note)

        st.markdown(
            "\n".join(f"- **{line.label}** — {line.value}" for line in recap.lines)
        )

        _render_attribution(attribution)

        # The route the request asked for. A page is enough here, unlike issue #58's
        # button: Configure Parameters has no sections to land wrongly between, and its own
        # `📊 Go to Data Analysis` is the return leg.
        if st.button("⚙️ Change these filters", key="recap_change_filters"):
            st.session_state.current_page = "parameter_config"
            st.rerun()


def _render_attribution(attribution) -> None:
    """How many variants each setting kept out, beneath the settings themselves.

    Drawn only where there is something to attribute. Two silences are deliberate:

    * **no attribution recorded** — a report produced before this existed, or by a call site
      that had no frame to measure. The settings list above is still the whole truth about
      the cut, so nothing is missing; a "not available" line would report on MAFigate rather
      than on the report.
    * **nothing was left out** — every variant is in the report, so there is no exclusion to
      apportion, and a table of zeroes beneath a complete report reads as a fault.

    The overlap sentence is not optional. These counts sum to more than the total, and a
    reader who adds them and overshoots has been told something that looks wrong — the same
    reason :func:`~filters.variant_filters._frequency_note` states its exemption count even
    when it is zero.
    """
    if attribution is None or not attribution.excluded_by:
        return

    st.markdown("---")
    st.markdown(
        f"**What this cut did** — {attribution.rows:,} variants read, "
        f"{attribution.in_report:,} in the report, {attribution.left_out:,} not."
    )
    st.caption(
        f"A variant can fall outside more than one setting, so these add up to more "
        f"than {attribution.left_out:,}."
    )
    st.markdown(
        "\n".join(
            f"- **{label}** — {count:,}"
            for label, count in attribution.excluded_by
        )
    )


def _render_variant_downloads(
    data: pd.DataFrame, dataset: str, grid_columns: list = None
) -> None:
    """The CSVs offered beneath a variant table: what the grid is showing, or every column.

    One function for both tabs, which had it written out twice — and issue #62 found what
    two verbatim copies of a results-page block do given time. They drift.

    **The narrow download follows the grid** (issue #92). ``grid_columns`` is the list
    ``create_data_table`` resolved for this tab on this render, so a column the user adds
    through *Add columns* is in the file offered directly beneath it. Before, this button
    served ``resolve_visible_columns``'s fixed list whatever the grid showed, so a user
    could search for a column, add it, see it, press this, and not find it — issue #60's
    shape one layer out, fixed there for the columns a user *invents* and not for the
    pipeline ones they *add*.

    Following it is safe in the direction that matters because **both grid controls are
    add-only**: *Add columns* appends from what is not already shown and *Show all
    columns* replaces the set with everything, so there is no control that hides one.
    This file is therefore always the report's own columns or a superset — never less
    than the fixed list used to deliver. That is why a download tracking a checkbox is
    not the reproducibility hazard it sounds like: the failure it would risk is
    unreachable.

    **Where the two coincide, one button is drawn.** Tick *Show all columns* and the
    grid's set becomes the whole frame, so the pair would serve one dataset twice —
    exactly the defect issue #83 spent a ticket removing, and one its guard could not
    see: the grid frame is ``prioritize_columns``-reordered and this one is not, so the
    payloads would have differed by column *order* alone and
    ``test_no_two_downloads_hand_over_the_same_bytes`` would have passed over them. The
    test here is set equality, so it also collapses when a user reaches the same place by
    adding every remaining column by hand.

    The label says **shown**, which issue #83's test used to refuse. That ban was right
    for the app it was written against — the line it drew was that a label may state what
    it delivers, never what the screen shows, *while the screen is user-controlled and
    the download is not*. The second half of that condition is what has changed. What
    holds the word true now is not a ban but
    ``test_the_shown_download_carries_exactly_the_grids_columns``, which is the stronger
    guard: it fails when the payload and the grid disagree, rather than when a word
    appears.

    Order is the pipeline's, not the screen's. ``grid_columns`` arrives in resolver order
    with additions appended; ``_PINNED_LEFT`` moves ``Notes`` and the clinical columns to
    the left edge of the *grid* only, downstream of the resolver since issue #60, so what
    the two share is the column set and not the layout.

    Both payloads come from one ``with_user_columns`` call. Each tab used to make two,
    which cost a full derived-column build per render for no difference in output.
    """
    exported = with_user_columns(data)
    # No grid to follow: create_data_table raised and the error fallback drew the table
    # instead. The report's own columns are the right answer then, being what the grid
    # would have opened with.
    columns = [
        col
        for col in (grid_columns if grid_columns is not None else shown_columns(exported))
        if col in exported.columns
    ]
    prefix = _download_prefix()

    if set(columns) == set(exported.columns):
        st.download_button(
            label=f"📥 Download all {len(exported.columns)} columns (CSV)",
            data=exported.to_csv(index=False),
            file_name=f"{prefix}_{dataset}_all.csv",
            mime="text/csv",
            key=f"download_{dataset}_all",
        )
        return

    shown_col, all_col = st.columns(2)
    with shown_col:
        st.download_button(
            label=f"📥 Download the {len(columns)} columns shown (CSV)",
            data=exported[columns].to_csv(index=False),
            file_name=f"{prefix}_{dataset}_{len(columns)}col.csv",
            mime="text/csv",
            key=f"download_{dataset}_shown",
        )
    with all_col:
        st.download_button(
            label=f"📥 Download all {len(exported.columns)} columns (CSV)",
            data=exported.to_csv(index=False),
            file_name=f"{prefix}_{dataset}_all.csv",
            mime="text/csv",
            key=f"download_{dataset}_all",
        )


def create_enhanced_data_table(
    passed_data: pd.DataFrame,
    failed_data: pd.DataFrame,
    height: int = TABLE_HEIGHT,
    summary_report: str = None,
    recap=None,
    stale: bool = False,
    on_reapply=None,
    attribution=None,
) -> None:
    """Create enhanced results view with passed/failed variants.

    ``summary_report`` is the run report as text, offered for download inside the Summary
    tab. It arrives already built rather than being fetched here: the report is
    ``page_modules.data_loading``'s, assembled out of session state, and this package must
    not import a page to draw one of its tabs. Passing ``None`` draws no button.

    ``recap``, ``stale`` and ``on_reapply`` arrive the same way and for the same reason
    (issue #137): the recap is a snapshot the *filter run* took, whether it is stale is a
    comparison against a stamp the page owns, and re-applying is the page's own function.
    Handed in rather than reached for, so this package still imports no page — the pattern
    ``summary_report`` set. ``attribution`` is the fourth of that set (issue #147) and
    travels with ``recap``, being the other half of one snapshot.
    """

    st.subheader("🔍 Variant Analysis Results")

    # Above the tabs rather than inside one: it describes all three of them.
    render_filter_recap(
        recap, stale=stale, on_reapply=on_reapply, attribution=attribution
    )

    # Create tabs for different views
    tab1, tab2, tab3 = st.tabs(["✅ Passed Filters", "❌ Failed Filters", "📊 Summary"])

    with tab1:
        # Heading only when there are rows: the empty case has its own warning below.
        if len(passed_data) > 0:
            st.markdown(f"**{len(passed_data):,} variants passed all applied filters**")

            # Clinical summary is already computed at filter time; just prioritize columns
            display_data = prioritize_columns(passed_data)

            # Create the table with error handling
            try:
                grid_columns = create_data_table(display_data, height, "passed_variants")
            except Exception as e:
                grid_columns = None
                render_fallback_table(
                    display_data, height, e, "passed variants table"
                )

            _render_variant_downloads(passed_data, "passed", grid_columns)

            # Quick stats
            col1, col2, col3 = st.columns(3)
            with col1:
                if "Hugo_Symbol" in display_data.columns:
                    st.write("**Top 5 Genes:**")
                    top_genes = display_data["Hugo_Symbol"].value_counts().head(5)
                    for gene, count in top_genes.items():
                        st.write(f"• {gene}: {count} variants")

            with col2:
                if "Variant_Classification" in display_data.columns:
                    st.write("**Variant Types:**")
                    var_types = (
                        display_data["Variant_Classification"].value_counts().head(5)
                    )
                    for vtype, count in var_types.items():
                        st.write(f"• {vtype}: {count}")

            with col3:
                if "Clinical_Summary" in display_data.columns:
                    st.write("**Clinical Significance:**")
                    clinical_summary = (
                        display_data["Clinical_Summary"].value_counts().head(5)
                    )
                    formatted_summary = format_clinical_summary_display(
                        clinical_summary
                    )
                    for sig, count in formatted_summary.items():
                        st.write(f"• {sig}: {count}")
        else:
            st.warning("⚠️ No variants passed the applied filters.")

    with tab2:
        if len(failed_data) > 0:
            st.markdown(f"**{len(failed_data):,} variants failed one or more filters**")

            # Clinical summary is already computed at filter time; just prioritize columns
            display_data = prioritize_columns(failed_data)

            # Create the table with error handling
            try:
                grid_columns = create_data_table(display_data, height, "failed_variants")
            except Exception as e:
                grid_columns = None
                render_fallback_table(
                    display_data, height, e, "failed variants table"
                )

            _render_variant_downloads(failed_data, "failed", grid_columns)

            # Quick stats
            col1, col2, col3 = st.columns(3)
            with col1:
                if "Hugo_Symbol" in display_data.columns:
                    st.write("**Top 5 Genes (Failed):**")
                    top_genes = display_data["Hugo_Symbol"].value_counts().head(5)
                    for gene, count in top_genes.items():
                        st.write(f"• {gene}: {count} variants")

            with col2:
                if "Variant_Classification" in display_data.columns:
                    st.write("**Variant Types (Failed):**")
                    var_types = (
                        display_data["Variant_Classification"].value_counts().head(5)
                    )
                    for vtype, count in var_types.items():
                        st.write(f"• {vtype}: {count}")

            with col3:
                if "Clinical_Summary" in display_data.columns:
                    st.write("**Clinical Significance (Failed):**")
                    clinical_summary = (
                        display_data["Clinical_Summary"].value_counts().head(5)
                    )
                    formatted_summary = format_clinical_summary_display(
                        clinical_summary
                    )
                    for sig, count in formatted_summary.items():
                        st.write(f"• {sig}: {count}")
        else:
            st.info("ℹ️ All variants passed the applied filters.")

    with tab3:
        st.markdown("**Overall Analysis Summary:**")

        # Create metrics
        col1, col2, col3, col4 = st.columns(4)
        total_variants = len(passed_data) + len(failed_data)

        with col1:
            st.metric("Total Variants", total_variants)

        with col2:
            st.metric("Passed Filters", len(passed_data))

        with col3:
            st.metric("Failed Filters", len(failed_data))

        with col4:
            if total_variants > 0:
                pass_rate = (len(passed_data) / total_variants) * 100
                st.metric("Pass Rate", f"{pass_rate:.1f}%")
            else:
                st.metric("Pass Rate", "N/A")

        # Rehoused from the load page's Data Overview (issue #65) — see `render_file_facts`
        # for what moved and what did not. Both frames, because these count over the whole
        # file rather than over either verdict.
        render_file_facts(passed_data, failed_data)

        # Interactive visualization dashboard
        if total_variants > 0:
            render_summary_dashboard(passed_data, failed_data)

        # The run report, under the tab that shows the run. It is not a third variant
        # table: it is the only artifact anywhere that records the parameters a cut was
        # made with, deliberately so since issue #28 removed the on-screen echoes.
        if summary_report is not None:
            st.markdown("---")
            st.download_button(
                label="📊 Download this report (TXT)",
                data=summary_report,
                file_name=f"{_download_prefix()}_report.txt",
                mime="text/plain",
                key="download_summary",
            )
