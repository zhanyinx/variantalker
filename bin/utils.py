# Python utility functions

import numpy as np
import pandas as pd
import pyranges

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
    "",
]


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


def filter_maf4tmb(
    maf: pd.DataFrame, depth: int = 25, alt_allele: int = 3, vaf_threshold: float = 0.05
) -> pd.DataFrame:
    """Filter dataframe containing maf data for tmb calculation according to https://jitc.bmj.com/content/8/1/e000147.long"""
    # Calculate TMB
    filter_variant_classifications = ["Silent", "Intron", "3'UTR", "5'UTR"]
    maf = maf[~maf["Variant_Classification"].isin(filter_variant_classifications)]
    maf["tumor_f"] = maf["t_alt_count"] / (maf["t_alt_count"] + maf["t_ref_count"])

    # calculate TMB according to https://jitc.bmj.com/content/8/1/e000147.long
    out = maf[
        ((maf["t_alt_count"] + maf["t_ref_count"]) >= depth)
        & (maf["t_alt_count"] >= alt_allele)
        & ((maf["tumor_f"] >= vaf_threshold))
    ]

    return out


def write_annovar_db_from_cancervar(file: str, outfile: str):
    """Get annovar db from cancervar init file"""
    out = open(outfile, "a")
    with open(file, "r") as f:
        for line in f:
            if line.startswith("database_names"):
                out.write(f"## Annovar {line}")
    out.close()


# def cbind_or_merge(left: pd.DataFrame, right: pd.DataFrame, left_on: list, right_on: list, bind: bool=False) -> pd.DataFrame:
#     """Merge or bind columns between data frames if bind is true and number of rows is the same."""

#     if bind:
#         assert len(left.index) == len(right.index)
#         result = pd.concat([left.reset_index(drop=True), right.reset_index(drop=True)], axis=1)


def assign_escat(
    maf: pd.DataFrame, escat: pd.DataFrame, tissue: str = None
) -> pd.DataFrame:
    """Assign escat score to maf."""
    maf["ESCAT"] = "."
    maf["ESCAT_TISSUE"] = "."
    maf["ESCAT_CANCER"] = "."

    # run over rows
    for idx, row in maf.iterrows():
        gene = row["Hugo_Symbol"]
        variant_type = row["Variant_Type"]

        # if protein change in missing, it's within intron/UTR/splice sites
        try:
            protein_change = row["Protein_Change"].replace("p.", "")
        except:
            continue

        if tissue is not None and tissue != "other":
            # check that all values matches
            subset = escat[
                (escat.mutation_type == variant_type)
                & (escat.gene == gene)
                & (
                    (escat.mutation == protein_change)
                    | (escat.mutation.isin(["INSERTION", "DELETION", "MUTATION"]))
                )
                & (escat.tissue == tissue)
            ].copy()
        else:
            subset = escat[
                (escat.mutation_type == variant_type)
                & (escat.gene == gene)
                & (
                    (escat.mutation == protein_change)
                    | (escat.mutation.isin(["INSERTION", "DELETION", "MUTATION"]))
                )
            ].copy()

        # if there is at least one row matching all the values
        if len(subset):
            # if it is a generatic mutation, check if it is specific to an exon
            if any(subset["transcript_exon"] != "NAN"):
                exon_number = int(
                    float(
                        subset.loc[
                            subset["transcript_exon"] != "NAN", "transcript_exon"
                        ]
                    )
                )
                if exon_number == int(float(row["Transcript_Exon"])):

                    best_score_idx = list(
                        subset.loc[subset["transcript_exon"] != "NAN", "ESCAT"].values
                    )

                    best_score_idx = best_score_idx.index(min(best_score_idx))
                    maf.loc[idx, "ESCAT"] = subset.loc[
                        subset["transcript_exon"] != "NAN", "ESCAT"
                    ].values[best_score_idx]

                    maf.loc[idx, "ESCAT_TISSUE"] = subset.loc[
                        subset["transcript_exon"] != "NAN", "tissue"
                    ].values[best_score_idx]

                    maf.loc[idx, "ESCAT_TISSUE"] = subset.loc[
                        subset["transcript_exon"] != "NAN", "cancer"
                    ].values[best_score_idx]

                    maf.loc[idx, "ESCAT"] = str(
                        list(
                            subset.loc[
                                subset["transcript_exon"] != "NAN", "ESCAT"
                            ].values
                        )
                    )
                    continue
            else:
                best_score_idx = list(subset["ESCAT"].values)
                best_score_idx = best_score_idx.index(min(best_score_idx))
                maf.loc[idx, "ESCAT"] = subset["ESCAT"].values[best_score_idx]
                maf.loc[idx, "ESCAT_TISSUE"] = subset["tissue"].values[best_score_idx]
                maf.loc[idx, "ESCAT_TISSUE"] = subset["cancer"].values[best_score_idx]

                continue
    return maf


def calculate_nmd(out, gr_nmd, filter_variant_type, filter_variant_classification):
    if filter_variant_classification:
        maf_gr = out[
            (out["Variant_Type"].isin(filter_variant_type))
            & (out["Variant_Classification"].isin(filter_variant_classification))
        ][["Chromosome", "Start_Position", "End_Position"]]
    else:
        maf_gr = out[(out["Variant_Type"].isin(filter_variant_type))][
            ["Chromosome", "Start_Position", "End_Position"]
        ]
    if len(maf_gr):
        maf_gr.columns = ["Chromosome", "Start", "End"]
        maf_gr = pyranges.PyRanges(df=maf_gr)
        tmb_nmd_escapees = np.sum(
            pyranges.count_overlaps(
                grs={"nmd": gr_nmd}, how="containment", features=maf_gr
            ).nmd
            != 0
        )
        return tmb_nmd_escapees

    return 0.0
