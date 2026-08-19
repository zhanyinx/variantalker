#!/usr/bin/env python

import argparse
import os

from numpy import append
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
        help="Maf file",
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
    write_metadata2file(f"Germline filters arguments: {args}", args.output)

    out = merge_with_fallback(left=maf, right = renovo, left_on = [
                "Chromosome",
                "Start_Position",
                "Reference_Allele",
                "Tumor_Seq_Allele2"
            ],
            right_on=["Chr", "Start", "Ref", "Alt"], cols_to_use = renovo.columns)
    
    out = out.drop(["Chr", "Start", "Ref", "Alt"], axis=1, errors="ignore")
    out = out.drop_duplicates()
    out.to_csv(args.output, sep="\t", index=False, mode="a")


if __name__ == "__main__":
    main()
