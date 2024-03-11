# Python utility functions

import numpy as np
import pandas as pd

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
    "not_provided",
]

CLINVAR_PATHO = [
    "Likely_pathogenic",
    "Likely_pathogenic,_low_penetrance",
    "Likely_pathogenic/Likely_risk_allele",
    "Likely_pathogenic/Pathogenic,_low_penetrance",
    "Likely_pathogenic|Affects",
    "Likely_pathogenic|association",
    "Likely_pathogenic|drug_response",
    "Likely_pathogenic|other",
    "Likely_pathogenic|risk_factor",
    "Pathogenic",
    "Pathogenic/Likely_pathogenic",
    "Pathogenic/Likely_pathogenic/Likely_risk_allele",
    "Pathogenic/Likely_pathogenic/Pathogenic,_low_penetrance",
    "Pathogenic/Likely_pathogenic/Pathogenic,_low_penetrance|other",
    "Pathogenic/Likely_pathogenic|drug_response",
    "Pathogenic/Likely_pathogenic|other",
    "Pathogenic/Likely_pathogenic|risk_factor",
    "Pathogenic/Likely_risk_allele",
    "Pathogenic/Likely_risk_allele|risk_factor",
    "Pathogenic|Affects",
    "Pathogenic|association",
    "Pathogenic|association|protective",
    "Pathogenic|confers_sensitivity",
    "Pathogenic|drug_response",
    "Pathogenic|drug_response|other",
    "Pathogenic|other",
    "Pathogenic|protective",
    "Pathogenic|risk_factor",
    "Conflicting_interpretations_of_pathogenicity",
    "Conflicting_interpretations_of_pathogenicity|Affects",
    "Conflicting_interpretations_of_pathogenicity|association",
    "Conflicting_interpretations_of_pathogenicity|association|other",
    "Conflicting_interpretations_of_pathogenicity|association|risk_factor",
    "Conflicting_interpretations_of_pathogenicity|drug_response",
    "Conflicting_interpretations_of_pathogenicity|drug_response|other",
    "Conflicting_interpretations_of_pathogenicity|other",
    "Conflicting_interpretations_of_pathogenicity|other|risk_factor",
    "Conflicting_interpretations_of_pathogenicity|protective",
    "Conflicting_interpretations_of_pathogenicity|risk_factor",
    "drug_response",
    "drug_response|other",
    "drug_response|risk_factor",
]

CIVIC_COLUMNS = [
    "CIViC_Allele",
    "CIViC_Consequence",
    "CIViC_SYMBOL",
    "CIViC_Entrez_Gene_ID",
    "CIViC_Feature_type",
    "CIViC_Feature",
    "CIViC_HGVSc",
    "CIViC_HGVSp",
    "CIViC_Variant_Name",
    "CIViC_Variant_ID",
    "CIViC_Variant_Aliases",
    "CIViC_Variant_URL",
    "CIViC_Molecular_Profile_Name",
    "CIViC_Molecular_Profile_ID",
    "CIViC_Molecular_Profile_Aliases",
    "CIViC_Molecular_Profile_URL",
    "CIViC_HGVS",
    "CIViC_Allele_Registry_ID",
    "CIViC_ClinVar_IDs",
    "CIViC_Molecular_Profile_Score",
    "CIViC_Entity_Type",
    "CIViC_Entity_ID",
    "CIViC_Entity_URL",
    "CIViC_Entity_Source",
    "CIViC_Entity_Variant_Origin",
    "CIViC_Entity_Status",
    "CIViC_Entity_Significance",
    "CIViC_Entity_Direction",
    "CIViC_Entity_Disease",
    "CIViC_Entity_Therapies",
    "CIViC_Entity_Therapy_Interaction_Type",
    "CIViC_Evidence_Phenotypes",
    "CIViC_Evidence_Level",
    "CIViC_Evidence_Rating",
    "CIViC_Assertion_ACMG_Codes",
    "CIViC_Assertion_AMP_Category",
    "CIViC_Assertion_NCCN_Guideline",
    "CIVIC_Assertion_Regulatory_Approval",
    "CIVIC_Assertion_FDA_Companion_Test",
]


def split_civic(row):
    """Split civic values column into multiple columns values."""
    try:
        parts = row.strip("[] ").split(",")
        values = [part.split("|") for part in parts]
        transformed_list = [
            [inner[i] for inner in values] for i in range(len(values[0]))
        ]
        return transformed_list
    except AttributeError:
        return ""


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


def write_metadata2file(text: str, outfile: str):
    """Append text to file."""
    out = open(outfile, "a")
    out.write(f"## Extra metadata: {text} \n")
    out.close()


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
