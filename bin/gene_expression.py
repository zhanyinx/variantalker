#!/usr/bin/env python
"""Extract the three key gene expression values from a count table.

Usage:
python gene_expression.py --fname_in <FNAME> --fname_out <FNAME>
"""

import argparse

import pandas as pd

GEP_WEIGHTS = {
    "CCL5": 0.011787376,
    "CD27": 0.073748858,
    "CD274": 0.044608585,
    "CD276": -0.022120323,
    "CD8A": 0.028755114,
    "CMKLR1": 0.153902619,
    "CXCL9": 0.073274920,
    "CXCR6": 0.009972831,
    "HLA.DQA1": 0.019394848,
    "HLA.DRB1": 0.057739716,
    "HLA.E": 0.070167658,
    "IDO1": 0.058018516,
    "LAG3": 0.121585139,
    "NKG7": 0.078158887,
    "PDCD1LG2 (PD-LG2)": -0.003252444,
    "PSMB10": 0.033727574,
    "STAT1": 0.255206918,
    "TIGIT": 0.082720169,
}
EXPRESSSION_GENES = ["CXCL9", "CD8A", "CD274"]
ENSEMBL_IDS = {
    "CCL5": "ENSG00000271503",
    "CD27": "ENSG00000139193",
    "CD274": "ENSG00000120217",
    "CD276": "ENSG00000103855",
    "CD8A": "ENSG00000153563",
    "CMKLR1": "ENSG00000174600",
    "CXCL9": "ENSG00000138755",
    "CXCR6": "ENSG00000172215",
    "HLA.DQA1": "ENSG00000196735",
    "HLA.DRB1": "ENSG00000196126",
    "HLA.E": "ENSG00000204592",
    "IDO1": "ENSG00000131203",
    "LAG3": "ENSG00000089692",
    "NKG7": "ENSG00000105374",
    "PDCD1LG2 (PD-LG2)": "ENSG00000197646",
    "PSMB10": "ENSG00000205220",
    "STAT1": "ENSG00000115415",
    "TIGIT": "ENSG00000181847",
}
ENSEMBL_KEYS = {v: k for k, v in ENSEMBL_IDS.items()}


def _parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--fname_in", type=str)
    parser.add_argument("-o", "--fname_out", type=str)
    return parser.parse_args()


def main():
    args = _parse_args()
    df = pd.read_csv(args.fname_in, delimiter="\t")
    df[["Name", "Transcript"]] = df["Name"].str.split(".", expand=True)

    # Raw gene signature
    df_genes = df.loc[
        df["Name"].isin([ENSEMBL_IDS[x] for x in EXPRESSSION_GENES]), ["Name", "TPM"]
    ]
    df_genes["Genes"] = df_genes["Name"].apply(ENSEMBL_KEYS.get)

    # Inflamed GEP signature
    df_gep = df.loc[
        df["Name"].isin([ENSEMBL_IDS[x] for x in GEP_WEIGHTS.keys()]), ["Name", "TPM"]
    ]
    df_gep["Genes"] = df_gep["Name"].apply(ENSEMBL_KEYS.get)
    df_gep["Weights"] = df_gep["Genes"].apply(GEP_WEIGHTS.get)
    df_gep["Norm_TPM"] = df_gep["TPM"] * df_gep["Weights"]
    gep_value = df_gep["Norm_TPM"].sum()

    # Merge and save
    df_genes = pd.concat(
        [
            df_genes,
            pd.DataFrame(
                [["T cell inflamed GEP signature", gep_value, pd.NA]],
                columns=["Name", "TPM", "Genes"],
            ),
        ],
    )
    df_genes.to_csv(args.fname_out, index=False)


if __name__ == "__main__":
    main()
