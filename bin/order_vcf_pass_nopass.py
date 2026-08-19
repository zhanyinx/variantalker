#!/usr/bin/env python

import argparse
import os
import pandas as pd

from utils import *


def _parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-v",
        "--vcf",
        type=str,
        required=True,
        help="Input vcf file",
    )
    parser.add_argument(
        "-n",
        "--nopass_file",
        type=str,
        required=True,
        help="Input nopass file",
    )
    parser.add_argument(
        "-p",
        "--pass_file",
        type=str,
        required=True,
        help="Input pass file",
    )
    parser.add_argument(
        "-vo",
        "--vcf_out",
        type=str,
        default="output.vcf",
        help="Output vcf file name, default output.vcf",
    )
    parser.add_argument(
        "-no",
        "--nopass_out",
        type=str,
        default="output.nopass.tsv",
        help="Output nopass tsv file name, default output.nopass.tsv",
    )
    parser.add_argument(
        "-po",
        "--pass_out",
        type=str,
        default="output.pass.tsv",
        help="Output pass tsv file name, default output.pass.tsv",
    )
    args = parser.parse_args()
    return args

def main():
    """Merge cancervar and with the corresponding maf file."""
    # Parse input
    args = _parse_args()
    
    if not os.path.exists(args.vcf):
        raise ValueError(f"Vcf file {args.vcf} does not exist!")
    if not os.path.exists(args.pass_file):
        raise ValueError(f"Pass file {args.pass_file} does not exist!")
    if not os.path.exists(args.nopass_file):
        raise ValueError(f"Nopass file {args.nopass_file} does not exist!")

    #VCF
    out = read_vcf(args.vcf)

    header_lines = read_header(args.vcf, "##")
    contig_order = extract_contig_order(header_lines)

    out = sort_vcf(out, contig_order)
    with open(args.vcf_out, "w") as f:
        for line in header_lines:
            f.write(line + "\n")
    out.to_csv(args.vcf_out, sep="\t", index=False, mode="a")

    #NOPASS
    try:
        out = pd.read_csv(args.nopass_file, sep="\t")
        out = sort_maf(out, contig_order)
    except pd.errors.EmptyDataError:
        out = pd.DataFrame()
    out.to_csv(args.nopass_out, sep="\t", index=False)

    #PASS
    try:
        out = pd.read_csv(args.pass_file, sep="\t")
        out = sort_maf(out, contig_order)
    except pd.errors.EmptyDataError:
        out = pd.DataFrame()
    out.to_csv(args.pass_out, sep="\t", index=False)

if __name__ == "__main__":
    main()
