#!/usr/bin/env python

import argparse
import os

from numpy import append
import pandas as pd

from utils import *


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


def somatic_filters(
    maf: pd.DataFrame,
    vaf: float,
    somatic_genes: str,
    cancervar_keep: list,
    civic_keep: list,
):
    """Set of somatic specific filters."""
    clinvar_exclude = CLINVAR_EXCLUDE
    escat_exclude = [
        "IIIA",
        "IIIB",
        "IIIC",
        ".",
        "V",
    ]

    filter_guidelines = (
        (maf["CancerVar"].isin(cancervar_keep))
        | (
            ~maf["ClinVar_VCF_CLNSIG"].isin(clinvar_exclude)
            & (~maf["ClinVar_VCF_CLNSIG"].isna())
        )
        | (~(maf["ESCAT"].isin(escat_exclude)))
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

    # keep all pathogenetic variants
    filter_patho = (
        (maf["CancerVar"].isin(["Tier_II_potential", "Tier_I_strong"]))
        | (maf["ClinVar_VCF_CLNSIG"].isin(CLINVAR_PATHO))
        | maf["CIViC_Evidence_Level"].apply(
            lambda x: has_element_from_list(x, ["A", "B"])
        )
    )

    return (
        filter_guidelines & filter_vaf & filter_genes,
        filter_patho & filter_vaf & filter_genes,
    )


def germline_filters(
    maf: pd.DataFrame,
    vaf: float,
    germline_genes: str,
    intervar_keep: list,
    renovo_keep: list,
):
    """Set of somatic specific filters."""
    clinvar_exclude = CLINVAR_EXCLUDE

    filter_guidelines = (
        (maf["InterVar"].isin(intervar_keep))
        | (
            ~maf["ClinVar_VCF_CLNSIG"].isin(clinvar_exclude)
            & (~maf["ClinVar_VCF_CLNSIG"].isna())
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

    # keep all pathogenetic variants
    filter_patho = (maf["InterVar"].isin(["Pathogenic", "Likely pathogenic"])) | (
        maf["ClinVar_VCF_CLNSIG"].isin(CLINVAR_PATHO)
    )

    return (
        filter_guidelines & filter_vaf & filter_genes,
        filter_patho & filter_vaf & filter_genes,
    )


def main():
    """Merge cancervar and with the corresponding maf file."""
    # Parse input
    args = _parse_args()

    if not os.path.exists(args.maf):
        raise ValueError(f"Maf file {args.maf} does not exist!")

    out = read_maf(args.maf)
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
        cancervar_keep = args.filter_cancervar.split(",")
        civic_keep = args.filter_civic.upper().split(",")
        out["filter_specific"], filter_patho = somatic_filters(
            out,
            somatic_genes=args.filter_genes_somatic,
            cancervar_keep=cancervar_keep,
            civic_keep=civic_keep,
            vaf=args.vaf_threshold,
        )

    if args.sample_type == "germline":
        intervar_keep = args.filter_intervar.split(",")
        renovo_keep = args.filter_renovo.split(",")
        out["filter_specific"], filter_patho = germline_filters(
            out,
            germline_genes=args.filter_genes_germline,
            intervar_keep=intervar_keep,
            vaf=args.vaf_threshold_germline,
            renovo_keep=renovo_keep,
        )

    out["filter"] = "NOPASS"
    out.loc[filter_patho, "filter"] = "PASS"
    out.loc[out[["filter_common", "filter_specific"]].all(axis=1), "filter"] = "PASS"
    out = out.drop(["filter_common", "filter_specific"], axis=1)

    writeheader(args.maf, args.output)
    overlapping_columns = list(set(out.columns) & set(COLUMNS_TO_REMOVE))
    out = out.drop(columns=overlapping_columns)
    out.to_csv(args.output, sep="\t", index=False, mode="a")

    keep = [
        "Tumor_Sample_Barcode",
        "Matched_Norm_Sample_Barcode",
        "project_id",
        "Hugo_Symbol",
        "Annotation_Transcript",
        "Chromosome",
        "Start_Position",
        "End_Position",
        "Variant_Classification",
        "Variant_Type",
        "Reference_Allele",
        "Tumor_Seq_Allele1",
        "Tumor_Seq_Allele2",
        "AAChange.refGene",
        "cDNA_Change",
        "Codon_Change",
        "Protein_Change",
        "Transcript_Exon",
        "tumor_f",
        "DP",
        "t_alt_count",
        "t_ref_count",
        "n_alt_count",
        "n_ref_count",
        "ClinVar_VCF_CLNSIG",
        "CancerVar",
        "ESCAT",
        "ESCAT_TISSUE",
        "ESCAT_CANCER",
        "CIViC_Evidence_Level",
        "CIViC_Evidence_Rating",
        "CIViC_Entity_Disease",
        "CIViC_Variant_URL",
        "CIViC_Entity_URL",
        "CIViC_Entity_Status",
        "am_class",
        "am_pathogenicity",
        "Otherinfo",
        "tumor_tissue",
        "cosmic",
        "Freq_ExAC_ALL",
        "Freq_esp6500siv2_all",
        "Freq_1000g2015aug_all",
        "gnomAD_exome_AF",
    ]

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
