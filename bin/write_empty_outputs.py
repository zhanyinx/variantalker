#!/usr/bin/env python

"""Write header-only output files for a sample that has no variants left after
standardization (empty input VCF, or no PASS/post-norm variants).

This keeps an empty/variant-free sample visible as "processed, zero variants"
instead of crashing the whole run. The column set mirrors the filtered TSVs
produced by filter_variants.py so downstream consumers (e.g. MAFigate) can read
the empty tables with the expected header.
"""

import argparse

from utils import KEEP


def _parse_args():
    parser = argparse.ArgumentParser(
        description="Write header-only outputs for a sample with no variants."
    )
    parser.add_argument(
        "--sample_type", required=True, choices=["somatic", "germline"]
    )
    parser.add_argument(
        "--patient", required=True, help="Patient/sample id used for file names."
    )
    parser.add_argument(
        "--skip_civic", action="store_true", help="Drop CIViC columns."
    )
    return parser.parse_args()


def keep_columns(sample_type, skip_civic):
    """Reproduce the keep-column selection of filter_variants.py.

    Operates on a copy of KEEP (KEEP is a shared module-level list).
    """
    keep = list(KEEP)

    if sample_type == "germline":
        for col in ("Tumor_Sample_Barcode", "cosmic", "Freq_ExAC_ALL"):
            if col in keep:
                keep.remove(col)
        if "CancerVar" in keep:
            idx = keep.index("CancerVar")
            keep.remove("CancerVar")
            keep.insert(idx, "InterVar")
        keep.append("RENOVO_Class")
        keep.append("RENOVO_pls")
        keep = [item for item in keep if not item.startswith("CIViC")]

    if skip_civic:
        keep = [item for item in keep if not item.startswith("CIViC")]

    return keep


def main():
    args = _parse_args()
    keep = keep_columns(args.sample_type, args.skip_civic)
    header = "\t".join(keep) + "\n"

    for suffix in (".maf", ".pass.tsv", ".nopass.tsv"):
        with open(f"{args.patient}{suffix}", "w", encoding="utf-8") as f:
            f.write(header)


if __name__ == "__main__":
    main()
