#!/usr/bin/env python

import argparse
import os

from numpy import append
import pandas as pd

from utils import *
from database.database_utils import *


def _parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-m",
        "--maf",
        type=str,
        default=None,
        required=True,
        help="Input maf file",
    )
    parser.add_argument(
        "-fc",
        "--filter_cancervar",
        type=str,
        default="Tier_II_potential,Tier_I_strong",
        help="Cancervar filters, available: Tier_II_potential,Tier_I_strong,Tier_III_Uncertain,Tier_IV_benign",
    )
    parser.add_argument(
        "-fci",
        "--filter_civic",
        type=str,
        default="A,B,C",
        help="Civic evidence level filter. Available values: A,B,C,D,E",
    )
    parser.add_argument(
        "-fi",
        "--filter_intervar",
        type=str,
        default="Pathogenic,Likely pathogenic",
        help="Intervar filters, available: Pathogenic,Likely pathogenic,Uncertain significance,Likely benign,Benign",
    )
    parser.add_argument(
        "-fr",
        "--filter_renovo",
        type=str,
        default="LP Pathogenic,IP Pathogenic,HP Pathogenic",
        help="Intervar filters, available: LP Pathogenic,IP Pathogenic,HP Pathogenic,LP Benign,IP Benign,HP Benign",
    )
    parser.add_argument(
        "-fcl",
        "--filter_clinvar",
        type=str,
        default="Pathogenic, Likely pathogenic",
        help="Clinvar filters, available: Pathogenic, Likely_pathogenic, Uncertain_significance, Conflicting_classifications_of_pathogenicity, Likely_benign, Benign, drug_response, not_provided, , no_classifications_from_unflagged_records, risk_factor, Affects, no_classification_for_the_single_variant, association, other, Uncertain_risk_allele, Likely_risk_allele, _low_penetrance, protective, association_not_found, Established_risk_allele, confers_sensitivity",
    )
    parser.add_argument(
        "-fe",
        "--filter_escat",
        type=str,
        default="IA, IB, IC, IIA, IIB",
        help="Escat filters, available: IA, IB, IC, IIA, IIB, IIIA, IIIB, V",
    )
    parser.add_argument(
        "-st",
        "--sample_type",
        type=str,
        default="somatic",
        help="set to germline if sample is germline",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="output.maf",
        help="Output file name, default output.txt",
    )
    parser.add_argument(
        "-md",
        "--min_depth",
        type=int,
        default=50,
        required=False,
        help="Minimum coverage to keep the variant. Default 50.",
    )
    parser.add_argument(
        "-fvc",
        "--filter_variant_classification",
        type=str,
        default="Silent,Intron,3'UTR,5'UTR,IGR,5'Flank,3'Flank,RNA",
        help="Available options can be found in Variant Classification here: https://gatk.broadinstitute.org/hc/en-us/articles/360035531732-Funcotator-Annotation-Specifications",
    )
    parser.add_argument(
        "-vt",
        "--vaf_threshold",
        type=float,
        default=0.01,
        help="VAF threshold for tumor, default 0.01",
    )
    parser.add_argument(
        "-vtg",
        "--vaf_threshold_germline",
        type=float,
        default=0.2,
        help="VAF threshold for normal, default 0.2",
    )
    parser.add_argument(
        "-fgs",
        "--filter_genes_somatic",
        type=str,
        default="null",
        help="file with list of Hugo_Symbol genes to be kept.",
    )
    parser.add_argument(
        "-fgg",
        "--filter_genes_germline",
        type=str,
        default="null",
        help="file with list of Hugo_Symbol genes to be kept.",
    )
    parser.add_argument(
        "-c",
        "--config",
        type=str,
        default=None,
        required=True,
        help="Path to cancervar or intervar config file.",
    )
    parser.add_argument(
        "-f",
        "--funcotator",
        type=str,
        default=None,
        required=True,
        help="Path to funcotator database folder.",
    )
    parser.add_argument(
        "-t",
        "--technology",
        type=str,
        default="illumina",
        help="Technology used for sequencing.",
    )
    parser.add_argument(
        "-av",
        "--annovar_version",
        type=str,
        default="v0",
        help="Annovar version",
    )
    parser.add_argument(
        "-gb",
        "--genome_build",
        type=str,
        default="GRCh38",
        help="Genome build",
    )
    parser.add_argument(
        "-ms",
        "--maf_skip",
        type=str,
        required=False,
        default=None,
        help="Maf file with variants that skipped annotation because found in db",
    )
    parser.add_argument(
        "-sc",
        "--skip_civic",
        action="store_true",
        help="Skip CIViC columns",
    )
    parser.add_argument(
        "-sp",
        "--skip_pathogenic",
        action="store_true",
        help="Skip automatic retention of pathogenic variants",
    )
    parser.add_argument(
        "-pid",
        "--projectid",
        type=str,
        default="other",
        help="The project id corresponding to the analysed sample",
    )
    args = parser.parse_args()
    return args


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



def main():
    """Merge cancervar and with the corresponding maf file."""
    # Parse input
    args = _parse_args()

    keep = KEEP


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

    if not os.path.exists(args.maf):
        raise ValueError(f"Maf file {args.maf} does not exist!")

    out = read_maf(args.maf)

    # Check if CIViC columns exist in the MAF file
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

    # get variants from skipped maf file: variants from db
    if args.maf_skip and os.path.exists(args.maf_skip) and os.path.getsize(args.maf_skip) != 0:
        maf_skip = pd.read_csv(args.maf_skip, sep="\t")
        out = outer_concat_preserve_order(out, maf_skip, fill_value="")
    
    # add project id
    out["project_id"] = args.projectid

    variant_classification_filter = args.filter_variant_classification.split(",")
    out["filter_common"] = common_filters(
        out,
        coverage=args.min_depth,
        variant_classification_filter=variant_classification_filter,
    )

    if args.sample_type not in ["somatic", "germline"]:
        raise ValueError(
            f"sample_type must be somatic or germline; Provided {args.sample_type}"
        )

    if args.sample_type == "somatic":
        cancervar_keep = [x.strip() for x in args.filter_cancervar.split(",") if x.strip()]
        civic_keep = [x.strip() for x in args.filter_civic.split(",") if x.strip()]
        escat_keep = [x.strip() for x in args.filter_escat.split(",") if x.strip()]
        clinvar_keep = [x.strip() for x in args.filter_clinvar.split(",") if x.strip()]
        out["filter_specific"], filter_patho = somatic_filters(
            out,
            vaf=args.vaf_threshold,
            somatic_genes=args.filter_genes_somatic,
            cancervar_keep=cancervar_keep,
            civic_keep=civic_keep,
            escat_keep=escat_keep,
            clinvar_keep=clinvar_keep,
            skip_civic=args.skip_civic,
            skip_pathogenic=args.skip_pathogenic,
        )

    if args.sample_type == "germline":
        intervar_keep = [x.strip() for x in args.filter_intervar.split(",") if x.strip()]
        clinvar_keep = [x.strip() for x in args.filter_clinvar.split(",") if x.strip()]
        renovo_keep = [x.strip() for x in args.filter_renovo.split(",") if x.strip()]
        out["filter_specific"], filter_patho = germline_filters(
            out,
            vaf=args.vaf_threshold_germline,
            germline_genes=args.filter_genes_germline,
            intervar_keep=intervar_keep,
            clinvar_keep=clinvar_keep,
            renovo_keep=renovo_keep,
            skip_pathogenic=args.skip_pathogenic,
        )

    out["filter"] = "NOPASS"
    out.loc[filter_patho, "filter"] = "PASS"
    out.loc[out[["filter_common", "filter_specific"]].all(axis=1), "filter"] = "PASS"
    out = out.drop(["filter_common", "filter_specific"], axis=1)

    writeheader(args.maf, args.output)
    out.to_csv(args.output, sep="\t", index=False, mode="a")

    if len(out[out["filter"] == "NOPASS"]):
        out_nopass = out.loc[
            out["filter"] == "NOPASS", ~(out == "__UNKNOWN__").all()
        ]  # remove unknown columns
        out_nopass = out_nopass[keep]
        out_nopass.to_csv(
            f"filtered.{args.output}.nopass.tsv",
            sep="\t",
            index=False,
        )

    if len(out[out["filter"] == "PASS"]):
        # filtering columns
        out = out.loc[
            out["filter"] == "PASS", ~(out == "__UNKNOWN__").all()
        ]  # remove unknown columns
        out = out[keep]
        out.to_csv(f"filtered.{args.output}.pass.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
