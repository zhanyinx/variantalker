#!/usr/bin/env python

import os
from numpy import append
import argparse

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
        "-c",
        "--cancervar",
        type=str,
        default=None,
        required=True,
        help="Cancervar file",
    )
    parser.add_argument(
        "-cc",
        "--config",
        type=str,
        default=None,
        required=True,
        help="Cancervar configuration file",
    )
    parser.add_argument(
        "-d",
        "--date",
        type=str,
        default=None,
        required=False,
        help="Analysis date.",
    )
    parser.add_argument(
        "-e",
        "--escat",
        type=str,
        default="/data/escat_tiering.csv",
        required=False,
        help="File containig escat scores. Default: /data/escat_tiering.csv",
    )
    parser.add_argument(
        "-g",
        "--germline",
        action="store_true",
        help="If set, germline mode.",
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
        "--tissue",
        type=str,
        default=None,
        required=False,
        help="Tissues available: lung, breast, colorectal, prostate, stomach,pancreatic, liver, other",
    )
    parser.add_argument(
        "-sc",
        "--skip_civic",
        action="store_true",
        help="Skip CIViC columns",
    )
    args = parser.parse_args()
    return args


def main():
    """Merge cancervar and with the corresponding maf file."""
    # Parse input
    args = _parse_args()

    if not os.path.exists(args.maf):
        raise ValueError(f"Maf file {args.maf} does not exist!")

    if not os.path.exists(args.cancervar):
        raise ValueError(f"Cancervar file {args.cancervar} does not exist!")

    maf = read_maf(args.maf)

    cancervar_cancervar = pd.read_csv(
        args.cancervar, sep="\t", low_memory=False
    )  # file with cancervar info
    cancervar = pd.read_csv(
        args.cancervar.replace(".cancervar", ".grl_p").replace(".intervar", ".grl_p"),
        sep="\t",
        low_memory=False,
    )  # file with all pred scores

    # remove unknown records
    maf = maf[maf["Chromosome"] != "__UNKNOWN__"]
    maf.Start_Position = pd.to_numeric(maf.Start_Position)

    # convert xx into chrxx
    if not any(["chr" in str(x) for x in cancervar["Chr"].values]):
        cancervar["Chr"] = "chr" + cancervar["Chr"].astype(str)
    if not any(["chr" in str(x) for x in cancervar_cancervar["#Chr"].values]):
        cancervar_cancervar["#Chr"] = "chr" + cancervar_cancervar["#Chr"].astype(str)
    if not any(["chr" in str(x) for x in maf["Chromosome"].values]):
        maf["Chromosome"] = "chr" + maf["Chromosome"].astype(str)

    # Remove duplicated rows in cancervar
    cancervar_cancervar = cancervar_cancervar.drop_duplicates(subset=["#Chr", "Start", "Ref", "Alt"], keep="first")
    cancervar = cancervar.drop_duplicates(subset=["Chr", "Start", "Ref", "Alt"], keep="first")

    writeheader(args.maf, args.output)
    write_annovar_db_from_cancervar(args.config, args.output)
    write_metadata2file(f"Somatic and germline filters arguments: {args}", args.output)

    # merge table based on mutation position, reference and alternative

    # only not repeated columns
    cols_to_use = cancervar_cancervar.columns.difference(maf.columns)

    # MNVs are not handled properly in funcotator, so we need to manually correct them with MNVTAG
    if "MNVTAG" in maf.columns:
        if len(maf.loc[~maf['MNVTAG'].isna(), "MNVTAG"]): # if there are MNVs
            tmp = maf.loc[~maf['MNVTAG'].isna(), "MNVTAG"].str.extract(r'_(?P<ref>[ACGT]+)->(?P<alt>[ACGT]+)$')
            maf.loc[~maf['MNVTAG'].isna(), ["Reference_Allele",
                        "Tumor_Seq_Allele2"]] = tmp.values

    out = merge_with_fallback(left=maf, right = cancervar_cancervar, left_on = [
                "Chromosome",
                "Start_Position",
                "Reference_Allele",
                "Tumor_Seq_Allele2"
            ],
            right_on=["#Chr", "Start", "Ref", "Alt"], 
            cols_to_use = cols_to_use
            ).drop(["#Chr","Start","Ref","Alt"], axis=1, errors="ignore")

    cols_to_use = cancervar.columns.difference(out.columns)
    out = merge_with_fallback(left=out, right = cancervar,left_on=[
                "Chromosome",
                "Start_Position",
                "Reference_Allele",
                "Tumor_Seq_Allele2",
                #"Ref.Gene",
            ],
            #right_on=["Chr", "Start", "Ref", "Alt", "Gene"], 
            right_on=["Chr", "Start", "Ref", "Alt"],
            cols_to_use = cols_to_use
            ).drop(["Chr", "Start", "Ref", "Alt"], axis=1, errors="ignore")

    # add project id and tumor tissue
    out["tumor_tissue"] = args.tissue

    # read and format escat scores
    escat = pd.read_csv(args.escat).astype(str)
    escat = escat.apply(lambda col: col.str.upper().str.replace("-", " "))
    escat = escat.replace(np.nan, 'NAN')
    escat.to_csv("escat_formatted_test_2.csv", index=False)

    # keep only SNP and INDEL
    if args.germline:
        escat = escat[
            (escat["info"] != "SOMATIC")
            & (escat["mutation_type"].isin(["SNP", "DEL", "INS"]))
        ]
    else:
        escat = escat[
            (escat["info"] != "GERMLINE")
            & (escat["mutation_type"].isin(["SNP", "DEL", "INS"]))
        ]

    # assign escat score to maf
    if args.germline:
        out = assign_escat(out, escat)
    else:
        out = assign_escat(out, escat, args.tissue.upper())

    # reformat CancerVar/Intervar
    if args.germline:
        out["InterVar"] = out[" InterVar: InterVar and Evidence "].str.extract(
            r"InterVar: (.*) PVS"
        )
    else:
        out["CancerVar"] = out[" CancerVar: CancerVar and Evidence "].str.extract(
            r"[\d]#(Tier_[\w\W]+) E"
        )

    # remove duplicated info
    if args.germline:
        duplicated_cols = ["Otherinfo1"]
    else:
        duplicated_cols = ["Otherinfo1", "Ref.Gene", "Gene"]
    out.drop(duplicated_cols, axis=1, inplace=True, errors="ignore")

    # rename cosmic column name
    if not args.germline:
        cosmic = [x for x in out.columns.values if "cosmic" in x]
        assert len(cosmic) <= 1
        if "cosmic" != cosmic[0]:
            out["cosmic"] = out[cosmic[0]]

    # fill in vaf values for iontorrent and alissa
    if out["t_ref_count"].isnull().values.any():
        if "FRO" in out.columns:
            out["t_ref_count"] = out["FRO"]
        elif "DP" in out.columns and "ALTCOUNT" in out.columns:
            out["t_ref_count"] = out["DP"] - out["ALTCOUNT"]
        else:
            Warning("t_ref_count column is empty, probably malformatted VCF input file")

    if out["t_alt_count"].isnull().values.any():
        if "FAO" in out.columns:
            out["t_alt_count"] = out["FAO"]
        elif "DP" in out.columns and "ALTCOUNT" in out.columns:
            out["t_alt_count"] = out["ALTCOUNT"]
        else:
            Warning("t_alt_count column is empty, probably malformatted VCF input file")

    if out["tumor_f"].isnull().values.any():
        if "AF" in out.columns:
            out["tumor_f"] = out["AF"]
        else:
            Warning("tumor_f column is empty, probably malformatted VCF input file")
    elif (
        not out["t_alt_count"].isnull().values.any()
        and not out["t_ref_count"].isnull().values.any()
    ):
        # recalculate tumor frequency in case of dragen
        out["tumor_f"] = out["t_alt_count"] / (out["t_alt_count"] + out["t_ref_count"])

    # recalculate depth
    if "DP" not in out.columns or out["DP"].isnull().values.any():
        out["DP"] = out["t_alt_count"] + out["t_ref_count"]

    # remove variants without support
    out = out[out["DP"] > 0]

    # drop duplicates
    out = out.drop_duplicates()
    
    if not args.germline and not args.skip_civic:
        if not out.empty:
            # split civic into multiple lines
            civic_splitted = out["CIVIC"].apply(split_civic).apply(pd.Series)
            if len(civic_splitted.columns) > 1:
                out[CIVIC_COLUMNS] = civic_splitted
            else:
                for colname in CIVIC_COLUMNS:
                    out[colname] = None
        else:
            # if out is empty, just add the columns
            for colname in CIVIC_COLUMNS:
                out[colname] = pd.Series(dtype="object")

    # write to file
    out.to_csv(args.output, sep="\t", index=False, mode="a")


if __name__ == "__main__":
    main()
