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
        help="If set, germline mode: InterVar instead of CancerVar.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="output.maf",
        help="Output file name, default output.txt",
    )
    parser.add_argument(
        "-p",
        "--projectid",
        type=str,
        default="other",
        help="The project id corresponding to the analysed sample",
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
        "-md",
        "--min_depth",
        type=str,
        default=50,
        required=False,
        help="Minimum coverage to keep the variant. Default 50.",
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

    writeheader(args.maf, args.output)
    write_annovar_db_from_cancervar(args.config, args.output)

    # merge table based on mutation position, reference and alternative

    # only not repeated columns
    cols_to_use = cancervar_cancervar.columns.difference(maf.columns)
    out = pd.merge(
        maf,
        cancervar_cancervar[cols_to_use],
        right_on=["#Chr", "Start", "Ref", "Alt"],
        left_on=[
            "Chromosome",
            "Start_Position",
            "Reference_Allele",
            "Tumor_Seq_Allele2",
        ],
    ).drop(
        [
            "#Chr",
            "Start",
            "Ref",
            "Alt",
        ],
        axis=1,
    )

    # only not repeated columns
    cols_to_use = cancervar.columns.difference(out.columns)
    out = pd.merge(
        out,
        cancervar[cols_to_use],
        left_on=[
            "Chromosome",
            "Start_Position",
            "Reference_Allele",
            "Tumor_Seq_Allele2",
            "Ref.Gene",
        ],
        right_on=["Chr", "Start", "Ref", "Alt", "Gene"],
    ).drop(["Chr", "Start", "Ref", "Alt"], axis=1)

    # recalculate tumor frequency when not present
    out["tumor_f"] = out["t_alt_count"] / (out["t_alt_count"] + out["t_ref_count"])

    # add project id and tumor tissue
    out["project_id"] = args.projectid
    out["tumor_tissue"] = args.tissue

    # read and format escat scores
    escat = pd.read_csv(args.escat)
    for column in escat:
        escat[column] = escat[column].astype(str)
        escat[column] = escat[column].str.upper().replace("-", " ")

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

    # add date
    if args.date:
        out["date"] = args.date

    # remove duplicated info
    if args.germline:
        duplicated_cols = ["Otherinfo1"]
    else:
        duplicated_cols = ["Otherinfo1", "Ref.Gene", "Gene"]
    out.drop(duplicated_cols, axis=1, inplace=True)

    # rename cosmic column name
    if not args.germline:
        cosmic = [x for x in out.columns.values if "cosmic" in x]
        assert len(cosmic) <= 1
        if "cosmic" != cosmic[0]:
            out["cosmic"] = out[cosmic[0]]

        # TODO remove this hack for reply dashboard
        if "cosmic95" not in cosmic:
            out["cosmic95"] = out[cosmic[0]]
        if cosmic[0] != "cosmic95":
            out.drop(cosmic[0], inplace=True, axis=1)

    # fill in vaf values for iontorrent
    if out["tumor_f"].isnull().values.any():
        if "AF" in out.columns:
            out["tumor_f"] = out["AF"]
        else:
            Warning("tumor_f column is empty, probably malformatted VCF input file")

    if out["t_ref_count"].isnull().values.any():
        if "RO" in out.columns:
            out["t_ref_count"] = out["FRO"]
        else:
            Warning("t_ref_count column is empty, probably malformatted VCF input file")

    if out["t_alt_count"].isnull().values.any():
        if "AO" in out.columns:
            out["t_alt_count"] = out["FAO"]
        else:
            Warning("t_alt_count column is empty, probably malformatted VCF input file")

    # save output file
    out = out.drop_duplicates()
    out.to_csv(args.output, sep="\t", index=False, mode="a")

    # filter variants
    filter_variant_classifications = [
        "Silent",
        "Intron",
        "3'UTR",
        "5'UTR",
        "IGR",
        "5'Flank",
        "3'Flank",
        "RNA",
    ]
    out = out[~out["Variant_Classification"].isin(filter_variant_classifications)]

    # filter on coverage
    out = out[(out["t_alt_count"] + out["t_ref_count"]) >= args.min_depth]

    # filter cancervar/intervar clinvar escat
    cancervar_keep = ["Tier_II_potential", "Tier_I_strong"]
    clinvar_exclude = CLINVAR_EXCLUDE
    escat_exclude = [
        "IIIA",
        "IIIB",
        "IIIC",
        ".",
        "V",
    ]

    if not args.germline:
        out = out[
            (out["CancerVar"].isin(cancervar_keep))
            | (~out["ClinVar_VCF_CLNSIG"].isin(clinvar_exclude))
            | (~(out["ESCAT"].isin(escat_exclude)))
        ]

    if args.germline:
        out = out[(out["tumor_f"] > args.vaf_threshold_germline)]
    else:
        out = out[(out["tumor_f"] > args.vaf_threshold)]

    if len(out):
        # filtering columns
        keep = [
            "Tumor_Sample_Barcode",
            "Matched_Norm_Sample_Barcode",
            "project_id",
            "Hugo_Symbol",
            "HGNC_RefSeq_IDs",
            "Chromosome",
            "Start_Position",
            "End_Position",
            "Variant_Classification",
            "Variant_Type",
            "Reference_Allele",
            "Tumor_Seq_Allele1",
            "Tumor_Seq_Allele2",
            "cDNA_Change",
            "Codon_Change",
            "Protein_Change",
            "Transcript_Exon",
            "tumor_f",
            "t_alt_count",
            "t_ref_count",
            "n_alt_count",
            "n_ref_count",
            "ClinVar_VCF_CLNSIG",
            "CancerVar",
            "ESCAT",
            "ESCAT_TISSUE",
            "ESCAT_CANCER",
            "Otherinfo",
            "tumor_tissue",
            "cosmic95",  # TODO cosmic update when change version cosmic
            "Freq_ExAC_ALL",
            "Freq_esp6500siv2_all",
            "Freq_1000g2015aug_all",
            "gnomAD_exome_AF",
        ]

        if args.germline:
            keep.remove("Tumor_Sample_Barcode")
            keep.remove("HGNC_RefSeq_IDs")
            keep.remove("cosmic95")
            keep.remove("Freq_ExAC_ALL")
            idx = keep.index("CancerVar")
            keep.remove("CancerVar")
            keep.insert(idx, "InterVar")

        out = out.loc[:, ~(out == "__UNKNOWN__").all()]  # remove unknown columns
        out = out[keep]
        out.to_csv(f"filtered.{args.output}.tsv", sep="\t", index=False)
    else:
        with open(f"filtered.{args.output}.tsv", "a") as f:
            f.write(
                f"EMPTY! No mutations passing filters! Double check VAF to make sure!"
            )


if __name__ == "__main__":
    main()
