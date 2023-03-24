import base64
import io
import numpy as np
import pandas as pd
import streamlit as st


# columns to keep by default

KEEP = [
    "Hugo_Symbol",
    "ClinVar_VCF_CLNSIG",
    "CancerVar",
    "InterVar",
    "RENOVO_Class",
    "ESCAT",
    "tumor_f",
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
    "t_alt_count",
    "t_ref_count",
    "n_alt_count",
    "n_ref_count",
    "ESCAT_TISSUE",
    "ESCAT_CANCER",
    "RENOVO_pls",
    "Otherinfo",
    "tumor_tissue",
    "COSMIC_total_alterations_in_gene",
    "cosmic",
    "Freq_ExAC_ALL",
    "Freq_esp6500siv2_all",
    "Freq_1000g2015aug_all",
    "gnomAD_exome_AF",
    "dbSNP_ID",
    "Tumor_Sample_Barcode",
    "Matched_Norm_Sample_Barcode",
    "project_id",
    "date",
]


def download_csv(
    df: pd.DataFrame, download_filename: str, download_link_text: str
) -> str:
    """Generates link to download csv of DataFrame.
    Args:
        df: DataFrame to download.
        download_filename: filename and extension of file. e.g. myplot.pdf
        download_link_text: Text to display for download link.
    """
    csv = df.to_csv(index=False).encode()
    b64 = base64.b64encode(csv).decode()
    href = f'<a href="data:file/csv;base64,{b64}" download="{download_filename}" target="_blank">{download_link_text}</a>'
    return href


def read_maf(file: io.StringIO) -> pd.DataFrame:
    # count comment lines
    c = 0
    for line in file:
        if line.startswith("#"):
            c = c + 1

    file.seek(0)  # reset the file cursor to the beginning
    maf = pd.read_csv(file, header=c, sep="\t", low_memory=False)
    return maf


@st.cache_data
def convert_df(df):
    # IMPORTANT: Cache the conversion to prevent computation on every rerun
    return df.to_csv(index=False).encode("utf-8")
