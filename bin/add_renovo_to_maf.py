#!/usr/bin/env python

import argparse
import os

from numpy import append
import pandas as pd

CLINVAR_EXCLUDE = [
        "Affects",
        "Affects|association",
        "Affects|risk_factor",
        "Benign",
        "Benign/Likely_benign",
        "Benign/Likely_benign|association",
        "Benign/Likely_benign|drug_response",
        "Benign/Likely_benign|drug_response|other",
        "Benign/Likely_benign|other",
        "Benign/Likely_benign|other|risk_factor",
        "Benign/Likely_benign|risk_factor",
        "Benign|association",
        "Benign|association|confers_sensitivity",
        "Benign|confers_sensitivity",
        "Benign|drug_response",
        "Benign|other",
        "Benign|protective",
        "Benign|risk_factor",
        "Likely_benign",
        "Likely_benign|drug_response|other",
        "Likely_benign|other",
        "Likely_benign|risk_factor",
        "association_not_found",
        "protective",
        "protective|risk_factor",
    ]

def _parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-m",
        "--maf",
        type=str,
        default=None,
        required=True,
        help="Maf file",
    )
    parser.add_argument(
        "-fm",
        "--filtered_maf",
        type=str,
        default=None,
        help="Filtered maf file, default maf file",
    )
    parser.add_argument(
        "-r",
        "--renovo",
        type=str,
        default=None,
        required=True,
        help="Renovo output file",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="output.maf",
        help="Output file name, default output.txt",
    )
    args = parser.parse_args()
    return args


def read_maf(file: str) -> pd.DataFrame:
    # count comment lines
    c = 0
    with open(file) as f:
        for line in f:
            if line.startswith("#"):
                c = c + 1

    maf = pd.read_csv(file, header=c, sep="\t", low_memory=False)
    return maf


def writeheader(file: str, outfile: str):
    """Write header from input maf file into outfile."""
    out = open(outfile, "w")
    with open(file, "r") as f:
        for line in f:
            if line.startswith("#"):
                out.write(line)

    out.close()


def main():
    """Merge cancervar and with the corresponding maf file."""
    # Parse input
    args = _parse_args()

    if not os.path.exists(args.maf):
        raise ValueError(f"Maf file {args.maf} does not exist!")

    if not os.path.exists(args.renovo):
        raise ValueError(f"Renovo output file {args.renovo} does not exist!")

    maf = read_maf(args.maf)
    renovo = pd.read_csv(args.renovo, sep="\t", low_memory=False)
    renovo = renovo[["Chr", "Start", "Ref", "Alt", "RENOVO_Class", "PL_score"]]
    renovo.columns = ["Chr", "Start", "Ref", "Alt", "RENOVO_Class", "RENOVO_pls"]

    writeheader(args.maf, args.output)

    # merge table based on mutation position, reference and alternative
    out = pd.merge(
        maf,
        renovo,
        right_on=["Chr", "Start", "Ref", "Alt"],
        left_on=[
            "Chromosome",
            "Start_Position",
            "Reference_Allele",
            "Tumor_Seq_Allele2",
        ],
    )
    out = out.drop(["Chr", "Start", "Ref", "Alt"], axis=1)
    out = out.drop_duplicates()
    out.to_csv(args.output, sep="\t", index=False, mode="a")

    filtered_maf = read_maf(args.filtered_maf)
    # if filtered file is not empty
    if "Chromosome" in filtered_maf.columns:
        # merge table based on mutation position, reference and alternative
        out = pd.merge(
            filtered_maf,
            renovo,
            right_on=["Chr", "Start", "Ref", "Alt"],
            left_on=[
                "Chromosome",
                "Start_Position",
                "Reference_Allele",
                "Tumor_Seq_Allele2",
            ],
        )
        out = out.drop(["Chr", "Start", "Ref", "Alt"], axis=1)
        intervar_keep = ["Pathogenic", "Likely pathogenic"]
        clinvar_exclude = CLINVAR_EXCLUDE
        renovo_keep = ["LP Pathogenic", "IP Pathogenic", "HP Pathogenic"]
        out = out[
            (out["InterVar"].isin(intervar_keep))
            | (~out["ClinVar_VCF_CLNSIG"].isin(clinvar_exclude))
            | (out["RENOVO_Class"].isin(renovo_keep))
        ]

        out.to_csv(f"filtered.{args.output}.tsv", sep="\t", index=False)
    else:
        with open(f"filtered.{args.output}.tsv", "a") as f:
            f.write(
                f"EMPTY: No mutations passing filters! Double check VAF to make sure!"
            )


if __name__ == "__main__":
    main()
