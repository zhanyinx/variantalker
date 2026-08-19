#!/usr/bin/env python

import argparse
import os
import hashlib

from numpy import append

from utils import *


def _parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-m",
        "--maf",
        type=str,
        nargs="+",
        default=None,
        required=True,
        help="Input maf file(s). When multiple are given they are union-concatenated "
        "preserving the column order of the first file (new columns appended).",
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
    
    keep = KEEP

    maf_files = args.maf if isinstance(args.maf, list) else [args.maf]
    for maf_file in maf_files:
        if not os.path.exists(maf_file):
            raise ValueError(f"Maf file {maf_file} does not exist!")

    # The first file is the reference for header (# comment lines) and column order.
    first_maf = maf_files[0]

    dfs = []
    all_columns = []

    for maf_file in maf_files:
        df = read_maf(maf_file)
        dfs.append(df)
        for col in df.columns:
            if col not in all_columns:
                all_columns.append(col)

    aligned = []
    for df in dfs:
        missing = [
            c for c in all_columns
            if c not in df.columns
        ]
        if missing:
            df = df.assign(
                **{c: "" for c in missing}
            )
        aligned.append(
            df[all_columns]
        )

    out = pd.concat(
        aligned,
        ignore_index=True,
        copy=False
    )

    # Read header and extract contig order from the reference file
    header_lines = read_header(first_maf, "#")
    contig_order = extract_contig_order(header_lines)

    writeheader(first_maf, args.output)
    overlapping_columns = list(set(out.columns) & set(COLUMNS_TO_REMOVE))
    columns_to_remove = np.array([False] * len(out.columns))
    columns_to_remove[(out == "__UNKNOWN__").all()] = True
    # Remove duplicate columns (same content, different name)
    seen = {}
    for i, col in enumerate(out.columns):
        h = hashlib.md5(
            pd.util.hash_pandas_object(
                out[col],
                index=False
            ).values.tobytes()
        ).hexdigest()

        if h in seen:
            prev = seen[h]
            if out[col].equals(out[prev]):
                columns_to_remove[i] = True
        else:
            seen[h] = col
    columns_to_remove[out.columns.isin(overlapping_columns)] = True
    columns_to_remove[out.columns.isin(keep)] = False
    out = out.loc[:, ~columns_to_remove]

    # Final MAF ordering step
    out = sort_maf(out, contig_order)

    out.to_csv(args.output, sep="\t", index=False, mode="a")

if __name__ == "__main__":
    main()
