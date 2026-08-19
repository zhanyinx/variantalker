"""The summary dashboard: six Plotly charts over the passed and failed frames.

One entry point, :func:`render_summary_dashboard`, which lays out the charts and
draws a caption in place of any whose column the MAF does not carry.

Every chart routes its layout through :func:`_compact_layout`; ask that for a taller
chart rather than re-setting the layout it just wrote. Each also carries a
Plotly-absent fallback, because ``plotly`` is declared in ``requirements.txt`` and is
not guaranteed to be installed.

:data:`CLINICAL_COLORS` and :data:`CHROMOSOME_ORDER` are pure and could live in
``config/``, and deliberately do not. ``CHROMOSOME_ORDER`` has one reader and would be
a module-crossing import for a single ``for``. ``CLINICAL_COLORS`` is a live question
rather than a settled fact — it is doing two jobs at once, giving the clinical axis its
severity order while colouring nothing, next to charts that colour by pass/fail — and
the map holds *what the dashboard's colours mean* open. Rehousing it now would look
like an answer.
"""

import streamlit as st
import pandas as pd
from typing import Optional

from .clinical_summary import strip_summary_glyph

try:
    import plotly.express as px
    import plotly.graph_objects as go

    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False


# ---------------------------------------------------------------------------
# Color constants
# ---------------------------------------------------------------------------

#: The clinical axis's category order — and, for now, nothing else.
#:
#: Keys must be ``Clinical_Summary`` labels with their glyph stripped, because that is what
#: the chart matches against; a key that matches nothing silently drops its category to the
#: unrecognised tail of the axis. Issue #98 added the six classes below ``Benign`` and
#: renamed ``Not Provided`` to ``No Classification``, following
#: ``clinical_summary.CLINICAL_HIERARCHY``, whose order this reproduces — a test holds the
#: two together, since a class added there and forgotten here reorders a clinician's axis
#: without failing anything.
#:
#: The ``Unknown`` *category* is gone from this table because no label reads that word: the
#: class is still keyed ``Unknown`` in ``CLINICAL_HIERARCHY``, but it renders as
#: ``Unrecognised Annotation``, and what this table matches against is the rendered label.
#:
#: The colours remain **unread** — this chart colours by pass/fail — which is the open
#: question the module docstring points at. New entries are given plausible values rather
#: than none so that answering it does not start from a half-filled table.
CLINICAL_COLORS = {
    "Pathogenic": "#d62728",
    "Pathogenic (low penetrance)": "#e377a2",
    "Likely Pathogenic": "#ff7f0e",
    "Likely Pathogenic (low penetrance)": "#ffbb78",
    "Uncertain Significance": "#f0c420",
    "Likely Benign": "#2ca02c",
    "Benign": "#1f77b4",
    "Disease Risk": "#8c564b",
    "Drug Response": "#9467bd",
    "Association or Trait": "#17becf",
    "Protective": "#98df8a",
    "No Classification": "#7f7f7f",
    "Unrecognised Annotation": "#bbbbbb",
    "No Clinical Data": "#c7c7c7",
    "Analysis Error": "#999999",
}

# Natural sort order for chromosomes
CHROMOSOME_ORDER = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _safe_column(df: pd.DataFrame, col: str) -> Optional[pd.Series]:
    """Return column Series if it exists and has non-null values, else None."""
    if col in df.columns:
        series = df[col].dropna()
        if len(series) > 0:
            return df[col]
    return None


def _compact_layout(fig, height: int = 350, top: int = 35) -> None:
    """Apply compact layout settings to a Plotly figure.

    The margins are a floor, not a budget: ``automargin`` lets each axis grow to
    fit the text it draws, which no fixed number can do for labels as different
    as a chromosome number and a gene symbol. Ask here for a taller chart or a
    deeper top margin rather than re-setting the layout this just wrote.
    """
    fig.update_layout(
        margin=dict(l=10, r=10, t=top, b=10),
        height=height,
        font=dict(size=11),
        # Keep category labels that look like numbers as categories. Streamlit
        # registers its own Plotly template and makes it the default, and unlike
        # Plotly's own it leaves autotypenumbers at "convert types" — under which
        # chromosome labels "1".."22" become a number line, so most tick labels
        # are dropped, absent chromosomes leave gaps, and X/Y/MT are placed by
        # coercion rather than named.
        autotypenumbers="strict",
    )
    # No-ops on the donut charts, which have no cartesian axes.
    fig.update_xaxes(automargin=True)
    fig.update_yaxes(automargin=True)


# ---------------------------------------------------------------------------
# Summary Dashboard Charts
# ---------------------------------------------------------------------------


def _plot_vaf_distribution(
    passed_data: pd.DataFrame, failed_data: pd.DataFrame
) -> None:
    """Overlapping histograms of VAF (tumor_f) for passed vs failed."""
    if not PLOTLY_AVAILABLE:
        all_vaf = pd.concat(
            [passed_data.get("tumor_f", pd.Series(dtype=float)),
             failed_data.get("tumor_f", pd.Series(dtype=float))]
        ).dropna()
        all_vaf = pd.to_numeric(all_vaf, errors="coerce").dropna()
        if len(all_vaf) > 0:
            st.bar_chart(all_vaf.value_counts(bins=20).sort_index())
        return

    fig = go.Figure()

    for df, name, color in [
        (passed_data, "Passed", "#2ca02c"),
        (failed_data, "Failed", "#d62728"),
    ]:
        col = _safe_column(df, "tumor_f")
        if col is not None:
            vals = pd.to_numeric(col, errors="coerce").dropna()
            if len(vals) > 0:
                fig.add_trace(
                    go.Histogram(
                        x=vals,
                        name=name,
                        marker_color=color,
                        opacity=0.6,
                        nbinsx=30,
                    )
                )

    fig.update_layout(
        title="VAF Distribution",
        xaxis_title="Variant Allele Frequency",
        yaxis_title="Count",
        barmode="overlay",
        legend=dict(orientation="h", yanchor="bottom", y=1.02),
    )
    # The legend sits above the plot, so the top margin holds it and the title.
    _compact_layout(fig, top=60)
    st.plotly_chart(fig, use_container_width=True)


def _plot_variant_classification(all_data: pd.DataFrame) -> None:
    """Donut chart of Variant_Classification (top 8 + Other)."""
    counts = all_data["Variant_Classification"].value_counts()
    if len(counts) > 8:
        top = counts.head(8)
        other = pd.Series({"Other": counts.iloc[8:].sum()})
        counts = pd.concat([top, other])

    if not PLOTLY_AVAILABLE:
        st.bar_chart(counts)
        return

    fig = px.pie(
        names=counts.index,
        values=counts.values,
        title="Variant Classification",
        hole=0.4,
    )
    fig.update_traces(textposition="inside", textinfo="percent+label")
    _compact_layout(fig)
    st.plotly_chart(fig, use_container_width=True)


def _plot_chromosome_distribution(all_data: pd.DataFrame) -> None:
    """Vertical bar chart of variant counts per chromosome.

    Vertical, unlike Top Genes below, because the chromosome is what a reader
    expects to find along the bottom of a per-chromosome distribution — and
    because the 25 labels of ``CHROMOSOME_ORDER``, one or two characters wide,
    fit along the bottom of a half-width column, which fifteen gene symbols
    would not.
    """
    chrom_col = all_data["Chromosome"].astype(str).str.replace("chr", "", regex=False)
    counts = chrom_col.value_counts()

    # Natural sort
    order = [c for c in CHROMOSOME_ORDER if c in counts.index]
    remaining = [c for c in counts.index if c not in order]
    order += sorted(remaining)
    counts = counts.reindex(order).dropna()

    if not PLOTLY_AVAILABLE:
        st.bar_chart(counts)
        return

    fig = px.bar(
        x=counts.index,
        y=counts.values,
        title="Chromosome Distribution",
        labels={"x": "Chromosome", "y": "Variant Count"},
    )
    # The bars are already in chromosome order; keep them that way rather than
    # letting the axis sort "10" before "9" or push X and Y around.
    fig.update_xaxes(categoryorder="trace")
    _compact_layout(fig, height=400)
    st.plotly_chart(fig, use_container_width=True)


def _plot_top_genes(all_data: pd.DataFrame, top_n: int = 15) -> None:
    """Horizontal bar chart of top mutated genes."""
    counts = all_data["Hugo_Symbol"].value_counts().head(top_n)

    if not PLOTLY_AVAILABLE:
        st.bar_chart(counts)
        return

    fig = px.bar(
        x=counts.values,
        y=counts.index,
        orientation="h",
        title=f"Top {top_n} Mutated Genes",
        labels={"x": "Variant Count", "y": "Gene"},
        color=counts.values,
        color_continuous_scale="Reds",
    )
    fig.update_layout(
        yaxis=dict(autorange="reversed"),
        coloraxis_showscale=False,
    )
    _compact_layout(fig, height=400)
    st.plotly_chart(fig, use_container_width=True)


def _plot_clinical_significance(
    passed_data: pd.DataFrame, failed_data: pd.DataFrame
) -> None:
    """Grouped bar chart of clinical significance: passed vs failed."""
    rows = []
    for df, label in [(passed_data, "Passed"), (failed_data, "Failed")]:
        col = _safe_column(df, "Clinical_Summary")
        if col is not None:
            for val, count in col.value_counts().items():
                # Named for the reader, not for this function: Plotly puts these
                # column names in the tooltip, so they are rendered strings.
                rows.append(
                    {
                        # Stripped by the module that writes the labels. This function used
                        # to hold its own list of nine emoji prefixes, and that copy quietly
                        # decided the axis: a label whose glyph was missing from it kept the
                        # glyph, stopped matching CLINICAL_COLORS, and dropped to the
                        # unrecognised tail of the order instead of its clinical place (#98).
                        "Clinical Significance": strip_summary_glyph(val),
                        "Status": label,
                        "Count": count,
                    }
                )

    if not rows:
        st.caption("No clinical significance data available.")
        return

    chart_df = pd.DataFrame(rows)

    if not PLOTLY_AVAILABLE:
        st.dataframe(chart_df)
        return

    # CLINICAL_COLORS already lists the categories in clinical order, most severe
    # first. Left to itself the axis follows value_counts(), which is frequency
    # order and reads as nonsense to someone scanning for severity — Pathogenic
    # landing between Likely Benign and Benign. Anything unrecognised keeps a
    # stable place at the end rather than being dropped.
    present = set(chart_df["Clinical Significance"])
    in_clinical_order = [name for name in CLINICAL_COLORS if name in present]
    in_clinical_order += sorted(present - set(CLINICAL_COLORS))

    fig = px.bar(
        chart_df,
        x="Clinical Significance",
        y="Count",
        color="Status",
        barmode="group",
        title="Clinical Significance (Passed vs Failed)",
        color_discrete_map={"Passed": "#2ca02c", "Failed": "#d62728"},
        category_orders={"Clinical Significance": in_clinical_order},
    )
    # Every bar is named and the title says what they are, so the axis needs no
    # title of its own. Cleared here rather than through ``labels``, which sets
    # the tooltip's field name too and would leave an orphaned "=" in it.
    fig.update_xaxes(title_text=None)
    _compact_layout(fig, height=380)
    st.plotly_chart(fig, use_container_width=True)


def _plot_mutation_types(all_data: pd.DataFrame) -> None:
    """Donut chart of Variant_Type (SNP, INS, DEL, etc.)."""
    counts = all_data["Variant_Type"].value_counts()

    if not PLOTLY_AVAILABLE:
        st.bar_chart(counts)
        return

    fig = px.pie(
        names=counts.index,
        values=counts.values,
        title="Mutation Type Spectrum",
        hole=0.4,
    )
    fig.update_traces(textposition="inside", textinfo="percent+label")
    _compact_layout(fig)
    st.plotly_chart(fig, use_container_width=True)


# ---------------------------------------------------------------------------
# Dashboard orchestrator
# ---------------------------------------------------------------------------


def render_summary_dashboard(
    passed_data: pd.DataFrame, failed_data: pd.DataFrame
) -> None:
    """Render the full graphical summary dashboard."""
    all_data = pd.concat([passed_data, failed_data], ignore_index=True)
    if len(all_data) == 0:
        st.info("No variant data to visualize.")
        return

    # Row 1: VAF + Variant Classification
    col1, col2 = st.columns(2)
    with col1:
        if _safe_column(all_data, "tumor_f") is not None:
            _plot_vaf_distribution(passed_data, failed_data)
        else:
            st.caption("VAF distribution not available: `tumor_f` column missing.")
    with col2:
        if _safe_column(all_data, "Variant_Classification") is not None:
            _plot_variant_classification(all_data)
        else:
            st.caption("Variant classification not available.")

    # Row 2: Chromosome + Top Genes
    col1, col2 = st.columns(2)
    with col1:
        if _safe_column(all_data, "Chromosome") is not None:
            _plot_chromosome_distribution(all_data)
        else:
            st.caption("Chromosome distribution not available.")
    with col2:
        if _safe_column(all_data, "Hugo_Symbol") is not None:
            _plot_top_genes(all_data)
        else:
            st.caption("Top genes not available: `Hugo_Symbol` column missing.")

    # Row 3: Clinical Significance (full width)
    _plot_clinical_significance(passed_data, failed_data)

    # Row 4: Mutation types
    if _safe_column(all_data, "Variant_Type") is not None:
        col1, col2 = st.columns(2)
        with col1:
            _plot_mutation_types(all_data)
    else:
        st.caption("Mutation type spectrum not available.")
