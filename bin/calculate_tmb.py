#!/usr/bin/env python
# Calculate TMB given a maf file

import os
import argparse

from utils import *


def _parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-a",
        "--alt_allele",
        type=int,
        default=3,
        help="Minimum number of alt alleles to call a variant, default 3",
    )
    parser.add_argument(
        "-d",
        "--depth",
        type=int,
        default=25,
        help="Depth threshold for TMB calculation, default=25",
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
        "-n",
        "--nmd_scores",
        type=str,
        default="/data/nmd_final_grange.tsv",
        help="NMD scores file. Default scores generated from NMDetectiveA, threshold >=0.25",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="output.maf",
        help="Output file name, default output.txt",
    )
    parser.add_argument(
        "-t",
        "--target_file",
        type=str,
        required=True,
        help="Target region used to calculate TMB (Default filters taken from (https://jitc.bmj.com/content/8/1/e000147.long))",
    )
    parser.add_argument(
        "-v",
        "--vaf_threshold",
        type=float,
        default=0.05,
        help="VAF threshold for tmb calculation, default 0.05",
    )

    args = parser.parse_args()
    return args


def main():
    """Calculate TMB given a maf file."""
    # Parse input
    args = _parse_args()

    if not os.path.exists(args.maf):
        raise ValueError(f"Maf file {args.maf} does not exist!")

    if not os.path.exists(args.target_file):
        raise ValueError(f"Target file {args.target_file} does not exist!")

    maf = read_maf(args.maf)
    target_region = pd.read_csv(args.target_file, sep="\t", header=None)
    size = np.sum(target_region[2] - target_region[1])

    # filter maf
    out = filter_maf4tmb(
        maf,
        depth=args.depth,
        alt_allele=args.alt_allele,
        vaf_threshold=args.vaf_threshold,
    )

    # Calculate TMB
    outfile = open(args.output, "a")
    for vaf in np.arange(0.05, 0.3, 0.05):
        out1 = out[(out["tumor_f"] > vaf)]
        outfile.write(f"TMB for VAF >= {round(vaf,2)}: {len(out1) / size * 1e6}\n")

    tmb_indel = len(out[out["Variant_Type"].isin(["INS", "DEL"])])
    outfile.write(f"TMB indel: {tmb_indel}\n")

    # NMD escapees
    nmd_scores = pd.read_csv(args.nmd_scores, sep="\t")
    gr_nmd = pyranges.PyRanges(df=nmd_scores)

    tmb_nmd_escapees_indel_fs = calculate_nmd(
        out=out,
        gr_nmd=gr_nmd,
        filter_variant_type=["INS", "DEL"],
        filter_variant_classification=["Frame_Shift_Ins", "Frame_Shift_Del"],
    )

    tmb_nmd_escapees = calculate_nmd(
        out=out,
        gr_nmd=gr_nmd,
        filter_variant_type=["INS", "DEL", "Nonsense_Mutation"],
        filter_variant_classification=None,
    )

    outfile.write(f"INDEL-fs NMD escape: {tmb_nmd_escapees_indel_fs}\n")
    outfile.write(f"NMD escape: {tmb_nmd_escapees}\n")
    outfile.close()


if __name__ == "__main__":
    main()
