#!/usr/bin/env python
# Calculate mutational signature
import argparse
import os

import pandas as pd

from SigProfilerMatrixGenerator import install as genInstall
from SigProfilerAssignment import Analyzer as Analyze


def _parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        default=None,
        required=True,
        help="Input folder with vcf file",
    )
    parser.add_argument(
        "-g",
        "--genome",
        type=str,
        default="GRCh38",
        required=False,
        help="Genome build used for calling variants in the vcf file.",
    )
    parser.add_argument(
        "-c",
        "--cosmic_version",
        type=float,
        default=3.4,
        required=False,
        help="Version of cosmic SBS database.",
    )
    parser.add_argument(
        "-cg",
        "--cosmic_group",
        type=str,
        required=False,
        help="File with names to collapse cosmic sbs into groups. E.g. APOBEC group",
    )
    parser.add_argument(
        "-e",
        "--exome",
        action="store_true",
        help="If defined, activate exome analysis",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="output",
        help="File output where to save group-activity values.",
    )

    args = parser.parse_args()
    return args


def main():
    """Extract mutational signatures following https://github.com/AlexandrovLab/SigProfilerAssignment."""
    # Parse input
    args = _parse_args()

    if not os.path.exists(args.input):
        raise ValueError(f"Input folder {args.input} does not exist!")

    if not args.genome in ["GRCh38", "GRCh37"]:
        raise ValueError(
            f"The genome must be GRCh37, GRCh38 but provided {args.genome}"
        )

    # genInstall.install(args.genome)
    Analyze.cosmic_fit(
        samples=args.input,
        output="out",
        input_type="vcf",
        context_type="96",
        collapse_to_SBS96=True,
        cosmic_version=args.cosmic_version,
        exome=False,
        genome_build=args.genome,
        signature_database=None,
        exclude_signature_subgroups=None,
        export_probabilities=False,
        export_probabilities_per_mutation=False,
        make_plots=False,
        sample_reconstruction_plots=False,
        verbose=False,
    )

    # read activity
    activity = pd.read_csv(
        f"out/Assignment_Solution/Activities/Assignment_Solution_Activities.txt",
        sep="\t",
    )

    # create grupping where each SBS is a group
    cosmic_group = pd.DataFrame(
        {"group": activity.columns[1:], "sbs": activity.columns[1:]}
    )

    # if defined custom grouping, overwrite single SBS grouping
    if args.cosmic_group != "null":
        cosmic_group = pd.read_csv(args.cosmic_group)
        cosmic_group.columns = ["group", "sbs"]

    # calculate the sum of activity of each group
    activity = pd.melt(activity)
    activity = pd.merge(activity, cosmic_group, left_on="variable", right_on="sbs")
    activity = activity.groupby("group")["value"].sum()
    activity.to_csv(args.output, sep=":", header=None)


if __name__ == "__main__":
    main()
