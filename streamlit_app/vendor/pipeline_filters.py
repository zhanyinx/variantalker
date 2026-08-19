"""Vendored verbatim from ``bin/filter_variants.py`` — DO NOT EDIT BY HAND.

Regenerate with ``python streamlit_app/vendor/_sync.py``; the drift guard in
``streamlit_app/tests/test_vendor_drift.py`` compares this file to ``bin/`` by
``ast.dump()`` equality, so any hand-edit here that changes behaviour turns the
suite red.

The first five functions are copies of the pipeline's own, compared symbol for
symbol. ``compute_keep`` is not: the pipeline has no such function, only eight
statements inside ``main()``. Five of them are lifted into it verbatim and
compared statement by statement; its **first** statement is deliberately not the
pipeline's, and its docstring says why.

Only the import header is local: ``bin/filter_variants.py`` cannot be imported by
the app, because it pulls in ``psycopg`` via ``database.database_utils``.
"""

import os

import pandas as pd

from .pipeline_utils import CLINVAR_PATHO, KEEP, has_clinvar_term


def common_filters(
    maf: pd.DataFrame, coverage: float, variant_classification_filter: list
) -> pd.DataFrame:
    """Set of common filters."""
    # filter variants
    filter_variant_classifications = variant_classification_filter
    return (~maf["Variant_Classification"].isin(filter_variant_classifications)) & (
        (maf["t_alt_count"] + maf["t_ref_count"]) >= coverage
    )


def has_element_from_list(s: str, my_list: list):
    """Check if any of the element in my_list is in the string s."""
    if pd.notna(s):
        for element in my_list:
            if element in s:
                return True
    return False


def check_civic_column_exists(maf: pd.DataFrame):
    """Check if CIViC_Evidence_Level column exists in the DataFrame."""
    civic_exists = 'CIViC_Evidence_Level' in maf.columns
    if not civic_exists:
        print("Warning: CIViC_Evidence_Level column not found in MAF file. CIViC filtering will be skipped.")
    return civic_exists


def somatic_filters(
    maf: pd.DataFrame,
    vaf: float,
    somatic_genes: str,
    cancervar_keep: list,
    civic_keep: list,
    escat_keep: list,
    clinvar_keep: list,
    skip_civic: bool,
    skip_pathogenic: bool = False,
):

    # Check if CIViC columns exist before using them
    civic_column_exists = check_civic_column_exists(maf)
    
    if skip_civic or not civic_column_exists:
        filter_guidelines = (
            maf["CancerVar"].isin(cancervar_keep)
            | maf["ClinVar_VCF_CLNSIG"].apply(
                lambda x: has_clinvar_term(x, clinvar_keep)
            )
            | maf["ESCAT"].isin(escat_keep)
        )
    else:
        filter_guidelines = (
            maf["CancerVar"].isin(cancervar_keep)
            | maf["ClinVar_VCF_CLNSIG"].apply(
                lambda x: has_clinvar_term(x, clinvar_keep)
            )
            | maf["ESCAT"].isin(escat_keep)
            | maf["CIViC_Evidence_Level"].apply(
                lambda x: has_element_from_list(x, civic_keep)
            )
        )

    # filter on variant allele frequency
    filter_vaf = maf["tumor_f"] > vaf

    # filter on genes (default no filters)
    filter_genes = maf["tumor_f"] > -1

    # filter list if defined
    if somatic_genes != "null":
        if os.path.exists(somatic_genes):
            genes = pd.read_csv(somatic_genes, header=None)
            filter_genes = (
                maf["Hugo_Symbol"].str.upper().isin(genes[0].str.upper().values)
            )
        else:
            Warning(f"{somatic_genes} file does not exist. No filters applied")

    # keep all pathogenic variants (unless skip_pathogenic is True)
    if skip_pathogenic:
        # If pathogenic filter is disabled, use empty filter (no pathogenic variants retained)
        filter_patho = pd.Series([False] * len(maf), index=maf.index)
    elif skip_civic or not civic_column_exists:
        filter_patho = (
            maf["CancerVar"].isin(["Tier_II_potential", "Tier_I_strong"])
            | maf["ClinVar_VCF_CLNSIG"].apply(
                lambda x: has_clinvar_term(x, CLINVAR_PATHO)
            )
        )
    else:
        filter_patho = (
            (maf["CancerVar"].isin(["Tier_II_potential", "Tier_I_strong"]))
            | maf["ClinVar_VCF_CLNSIG"].apply(
                lambda x: has_clinvar_term(x, CLINVAR_PATHO)
            )
            | maf["CIViC_Evidence_Level"].apply(
                lambda x: has_element_from_list(x, ["A", "B"])
            )
        )

    return (
        filter_guidelines & filter_vaf & filter_genes,
        filter_patho,
    )


def germline_filters(
    maf: pd.DataFrame,
    vaf: float,
    germline_genes: str,
    intervar_keep: list,
    renovo_keep: list,
    clinvar_keep: list,
    skip_pathogenic: bool = False,
):

    filter_guidelines = (
        (maf["InterVar"].isin(intervar_keep))
        | maf["ClinVar_VCF_CLNSIG"].apply(
            lambda x: has_clinvar_term(x, clinvar_keep)
        )
        | (maf["RENOVO_Class"].isin(renovo_keep))
    )

    # filter on variant allele frequency
    filter_vaf = maf["tumor_f"] > vaf

    # filter on genes (default no filters)
    filter_genes = maf["tumor_f"] > -1

    # filter list if defined
    if germline_genes != "null":
        if os.path.exists(germline_genes):
            genes = pd.read_csv(germline_genes, header=None)
            filter_genes = (
                maf["Hugo_Symbol"].str.upper().isin(genes[0].str.upper().values)
            )
        else:
            Warning(f"{germline_genes} file does not exist. No filters applied")

    # keep all pathogenic variants (unless skip_pathogenic is True)
    if skip_pathogenic:
        # If pathogenic filter is disabled, use empty filter (no pathogenic variants retained)
        filter_patho = pd.Series([False] * len(maf), index=maf.index)
    else:
        filter_patho = (
            maf["InterVar"].isin(["Pathogenic", "Likely pathogenic"])
            | maf["ClinVar_VCF_CLNSIG"].apply(
                lambda x: has_clinvar_term(x, CLINVAR_PATHO)
            )
        )

    return (
        filter_guidelines & filter_vaf & filter_genes,
        filter_patho,
    )


def compute_keep(args, out):
    """The pipeline's output column list — `keep` in `bin/filter_variants.py:main()`.

    The statements after the first are lifted verbatim out of `main()`, so they read
    the pipeline's own names. Supply both as call-site shims: `args` is anything with
    `sample_type` and `skip_civic` attributes, and `out` is the frame whose `.columns`
    decide the CIViC and gnomAD branches. Returns the ordered column list.

    THE FIRST STATEMENT IS NOT THE PIPELINE'S. `main()` writes `keep = KEEP`, which
    aliases the module-level list from `pipeline_utils` and then mutates it in place.
    A pipeline process runs `main()` once and exits, so it never sees the consequence;
    the app keeps this module loaded for the whole session, where a germline call
    leaves the shared list three entries short and the next somatic call silently
    returns germline columns. See `tests/test_vendor_compute_keep.py`.

    The copy is taken here, outside the verbatim block, rather than by editing the
    vendored statements — so the drift guard still compares unedited copies.
    """
    keep = list(KEEP)

    if args.sample_type == "germline":
        keep.remove("Tumor_Sample_Barcode")
        keep.remove("cosmic")
        keep.remove("Freq_ExAC_ALL")
        idx = keep.index("CancerVar")
        keep.remove("CancerVar")
        keep.insert(idx, "InterVar")
        keep.append("RENOVO_Class")
        keep.append("RENOVO_pls")
        keep = [item for item in keep if not item.startswith("CIViC")]

    civic_columns_exist = any(col in out.columns for col in ["CIViC_Evidence_Level", "CIViC_Evidence_Rating", "CIViC_Entity_Disease", "CIViC_Variant_URL", "CIViC_Entity_URL", "CIViC_Entity_Status"])

    if args.skip_civic or not civic_columns_exist:
        keep = [item for item in keep if not item.startswith("CIViC")]
        if not civic_columns_exist and not args.skip_civic:
            print("Warning: CIViC columns not found in MAF file. CIViC columns will be automatically excluded from output.")

    if "gnomAD_exome_AF_raw" in out.columns.values:
        keep.append("gnomAD_exome_AF_raw")
    elif "gnomAD_exome_AF" in out.columns.values:
        keep.append("gnomAD_exome_AF")

    if "gnomAD_genome_AF_raw" in out.columns.values:
        keep.append("gnomAD_genome_AF_raw")
    elif "gnomAD_genome_AF" in out.columns.values:
        keep.append("gnomAD_genome_AF")

    return keep
