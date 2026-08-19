"""
Home page for MAFigate application.

The page is deliberately short: the two buttons that start the work, and the explanatory
text folded into collapsed expanders behind them (issue #61). The app name, version and
tagline are drawn once per page by ``render_header``, so this page draws no title of its
own.
"""

import streamlit as st
from config.constants import APP_NAME


def show_home_page():
    """Display the home page."""

    col1, col2 = st.columns(2)

    with col1:
        if st.button("⚙️ Configure Parameters", use_container_width=True):
            st.session_state.current_page = "parameter_config"
            st.rerun()

    with col2:
        if st.button("📊 Load & Analyze Data", use_container_width=True):
            st.session_state.current_page = "data_loading"
            st.rerun()

    with st.expander(f"ℹ️ What {APP_NAME} does"):
        st.markdown(
            f"""
            {APP_NAME} reads an annotated MAF file — `.maf`, `.tsv` or `.txt` — and filters it
            down to the variants you need to look at.

            Choose the analysis type that matches your sample: **somatic** for tumor variants,
            **germline** for hereditary ones. Each has its own thresholds and its own clinical
            classifications.

            Filters cover sequencing depth and variant allele frequency, variant classification,
            gene, population frequency, and the clinical classifications for the analysis type
            you chose. Pathogenic variants can be kept regardless of the other thresholds.

            Nothing is thrown away: the variants that pass and the variants that fail are shown
            separately, and both can be downloaded along with a summary of the run.
            """
        )

    # Only the sources the app actually reads are named here, and only in prose. The catalogue
    # this replaced listed fifteen tools, several of which appear nowhere in the app.
    with st.expander("ℹ️ Annotation sources"):
        st.markdown(
            f"""
            {APP_NAME} does not annotate your file. It reads annotations that are already in it.

            For somatic samples it reads **CancerVar** (pathogenicity tiers), **CIViC**
            (clinical evidence) and **ESCAT** (actionability). For germline samples it reads
            **InterVar** (ACMG classification) and **RENOVO** (predicted pathogenicity).
            **ClinVar** is read for both, as are the population frequency columns your file
            carries — gnomAD, ExAC and 1000 Genomes among them.

            *Help & Documentation* lists the columns {APP_NAME} reads, and what happens when
            your file is missing one.
            """
        )

    with st.expander("ℹ️ Quick start"):
        st.markdown(
            """
            1. **Configure Parameters** — choose somatic or germline, then set your thresholds
               or start from a preset.
            2. **Load & Analyze Data** — choose your MAF file in the sidebar, then apply
               the filters.
            3. Review what passed and what failed, then download the tables you need.

            Every setting has a tooltip explaining what it does; *Help & Documentation* has the
            longer version.
            """
        )
