import base64
import numpy as np
import pandas as pd
import streamlit as st


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


def read_in_chunks(file_object, chunk_size=10):
    """Lazy function (generator) to read a file piece by piece.
    Default chunk size: 1k."""
    while True:
        data = file_object.read(chunk_size)
        if not data:
            break
        yield data


def read_maf(file: str) -> pd.DataFrame:
    # count comment lines
    c = 0
    with open(file.read(), "r") as f:
        print(f)

    maf = pd.read_csv(file, header=c, sep="\t", low_memory=False)
    return


@st.cache
def convert_df(df):
    # IMPORTANT: Cache the conversion to prevent computation on every rerun
    return df.to_csv(index=False).encode("utf-8")
