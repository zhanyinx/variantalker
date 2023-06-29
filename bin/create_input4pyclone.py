#!/usr/bin/env python
# Calculate TMB given a maf file

import os
import sys
import argparse

import numpy as np
import pandas as pd
import pyranges

from utils import read_maf


def _parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-as",
        "--ascat",
        type=str,
        required=True,
        help="ascat allele count file from sarek.",
    )
    parser.add_argument(
        "-ac",
        "--ascat_cellularity",
        type=str,
        required=False,
        help="ascat cellularity count file from sarek.",
    )
    parser.add_argument(
        "-c",
        "--cellularity",
        type=float,
        default=-1,
        help="Cellularity of the sample, default 1",
    )
    parser.add_argument(
        "-m",
        "--maf",
        type=str,
        default=None,
        required=True,
        help="Maf file",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="output.maf",
        help="Output file name, default output.txt",
    )
    parser.add_argument(
        "-s",
        "--sex",
        type=str,
        default="XX",
        help="XX or XY",
    )

    args = parser.parse_args()
    return args


def main():
    """Create input file for pyclone given ascat copy number and maf file."""
    # Parse input
    args = _parse_args()

    if not os.path.exists(args.maf):
        raise ValueError(f"Maf file {args.maf} does not exist!")

    if not os.path.exists(args.ascat):
        raise ValueError(f"Ascat file {args.ascat} does not exist!")

    # read and process annotated maf file
    maf = read_maf(args.maf)
    sampleid = maf["Tumor_Sample_Barcode"].unique()[0]

    maf = maf[
        [
            "Chromosome",
            "Start_Position",
            "End_Position",
            "t_ref_count",
            "t_alt_count",
            "Genome_Change",
        ]
    ]
    maf.columns = [
        "Chromosome",
        "Start",
        "End",
        "ref_counts",
        "alt_counts",
        "mutation_id",
    ]
    maf["normal_cn"] = 2
    maf = pyranges.PyRanges(maf)

    # read and process ascat cnv file
    try:
        cnv_data = pd.read_csv(args.ascat, sep="\t")
    except pd.io.common.EmptyDataError:
        # write output for pyclone
        maf["sample_id"] = "sampleid"
        maf["ref_counts"] = maf["t_ref_count"]
        maf["alt_counts"] = maf["t_alt_count"]
        maf["normal_cn"] = 2
        maf["minor_cn"] = 0
        maf["major_cn"] = 2
        if args.cellularity != -1:
            maf["tumour_content"] = args.cellularity
        else:
            maf["tumor_content"] = 1
        maf["tumor_f"] = maf["tumor_f"]
        maf = maf[
            [
                "mutation_id",
                "sample_id",
                "ref_counts",
                "alt_counts",
                "normal_cn",
                "major_cn",
                "minor_cn",
                "tumour_content",
            ]
        ]
        maf.to_csv(args.output, index=False, sep="\t", mode="a", header=None)
        sys.exit(0)

    if not any(["chr" in str(x) for x in cnv_data["chr"].values]):
        cnv_data["chr"] = "chr" + cnv_data["chr"].astype(str)
    cnv_data.columns = ["Chromosome", "Start", "End", "major_cn", "minor_cn"]
    cnv_data = pyranges.PyRanges(cnv_data)

    # joining
    joint = maf.join(cnv_data, suffix="_cnv", how="left").df
    joint = joint.drop([x for x in joint.columns if "_cnv" in x], axis=1)
    joint["sample_id"] = sampleid

    # add cellularity: if not provided, use that from ascat
    if args.cellularity != -1:
        joint["tumour_content"] = args.cellularity
    else:
        if not os.path.exists(args.ascat_cellularity):
            raise ValueError(f"Ascat file {args.ascat_cellularity} does not exist!")
        fraction = pd.read_csv(args.ascat_cellularity, sep="\t")
        joint["tumour_content"] = fraction["AberrantCellFraction"].values[0]

    # add sex
    if args.sex == "XX":
        joint["normal_cn"] = 2
    elif args.sex == "XY":
        joint["normal_cn"] = 2
        if joint["Chromosome"].isin(["chrX", "chrY"]).any():
            joint.loc[joint["Chromosome"].isin(["chrX", "chrY"]), "normal_cn"] = 1
    else:
        raise ValueError(f"sex must be XX or XY, provided {args.sex}.")

    joint = joint.drop(["Chromosome", "Start", "End"], axis=1)

    joint = joint[
        [
            "mutation_id",
            "sample_id",
            "ref_counts",
            "alt_counts",
            "normal_cn",
            "major_cn",
            "minor_cn",
            "tumour_content",
        ]
    ]

    # no overlapping regions are set to -1, since we want major and minor = 1 when no cnv is present, we can just take the module
    joint["major_cn"] = np.abs(joint["major_cn"])
    joint["minor_cn"] = np.abs(joint["minor_cn"])
    joint.to_csv(args.output, index=False, sep="\t", mode="a", header=None)


if __name__ == "__main__":
    main()
