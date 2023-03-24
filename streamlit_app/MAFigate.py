import base64
import io

import numpy as np
import pandas as pd
import streamlit as st

from st_aggrid import AgGrid, GridOptionsBuilder, GridUpdateMode, DataReturnMode
from utils import KEEP, download_csv, convert_df, read_maf


def main():
    """
    This is the main function.
    """

    st.title("MAFigate app")
    st.markdown(
        "This application enables users to view MAF files, apply filtering, and save their results. With this tool, users can easily navigate through MAF files and perform customized filtering to extract relevant information. The saved results can be conveniently accessed and shared for further analysis or reporting purposes."
    )

    # Upload a list of files
    uploaded_file = st.sidebar.file_uploader(
        "Please upload your data file(s)",
        type=["maf", "tsv"],
        accept_multiple_files=False,
    )

    if uploaded_file is not None:
        content = uploaded_file.read().decode("utf-8")
        file_buffer = io.StringIO(content)
        data = read_maf(file_buffer)
        data[" CancerVar: CancerVar and Evidence "] = data[
            " CancerVar: CancerVar and Evidence "
        ].str.extract(r"[\d]#(Tier_[\w\W]+) E")
        out = data.copy()

        keep = np.array([x for x in KEEP if x in data.columns])

        columns2keep = st.sidebar.multiselect(
            "Choose extra columns to keep", data.columns.drop(keep)
        )
        vaf = float(st.sidebar.text_input("VAF threshold", "0.05"))
        out = out.loc[(out["tumor_f"] > vaf), np.append(keep, columns2keep)]
        genes2keep = st.sidebar.multiselect(
            "Choose list of interesting genes",
            data["Hugo_Symbol"].unique(),
            default=None,
        )

        clinvar_options = data["ClinVar_VCF_CLNSIG"].dropna()
        clinvar_keep = st.sidebar.multiselect(
            "Clinvar class to keep",
            data["ClinVar_VCF_CLNSIG"].unique(),
            default=clinvar_options[~clinvar_options.str.contains("enign")].unique(),
        )

        cancervar_keep = st.sidebar.multiselect(
            "Cancervar class to keep",
            data["CancerVar"].unique(),
            default=["Tier_I_strong", "Tier_II_potential", "Tier_III_Uncertain"],
        )

        out = out[
            (
                (out["CancerVar"].isin(cancervar_keep))
                | (out["ClinVar_VCF_CLNSIG"].isin(clinvar_keep))
            )
        ]

        if genes2keep != []:
            out = out[out["Hugo_Symbol"].isin(genes2keep)]

        gb = GridOptionsBuilder.from_dataframe(out)
        gb.configure_default_column(
            enablePivot=True, enableValue=True, enableRowGroup=True
        )
        gb.configure_selection(selection_mode="multiple", use_checkbox=True)
        gb.configure_side_bar()
        gridoptions = gb.build()

        response = AgGrid(
            out,
            height=800,
            width=1500,
            gridOptions=gridoptions,
            enable_enterprise_modules=True,
            update_mode=GridUpdateMode.MODEL_CHANGED,
            data_return_mode=DataReturnMode.FILTERED_AND_SORTED,
            fit_columns_on_grid_load=False,
            header_checkbox_selection_filtered_only=True,
            use_checkbox=True,
        )

        v = response["selected_rows"]
        if v:
            st.write("Selected rows")
            st.dataframe(v)
            dfs = pd.DataFrame(v)
            csv = convert_df(dfs)

            st.download_button(
                label="Download data as CSV",
                data=csv,
                file_name="selected.csv",
                mime="text/csv",
            )


if __name__ == "__main__":
    main()
