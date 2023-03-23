from utils import *
import streamlit as st

from st_aggrid import AgGrid, GridOptionsBuilder, GridUpdateMode, DataReturnMode


def main():
    """
    This is the main function.
    """

    st.title("Navigate maf file")

    # Upload a list of files
    uploaded_file = st.sidebar.file_uploader(
        "Please upload your data file(s)",
        type=["maf", "tsv"],
        accept_multiple_files=False,
    )

    if uploaded_file is not None:
        data = pd.read_csv(uploaded_file, sep="\t", low_memory=False)
        data[" CancerVar: CancerVar and Evidence "] = data[
            " CancerVar: CancerVar and Evidence "
        ].str.extract(r"[\d]#(Tier_[\w\W]+) E")
        out = data.copy()

        keep = [
            "Tumor_Sample_Barcode",
            "Matched_Norm_Sample_Barcode",
            "project_id",
            "date",
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
            "InterVar",
            "RENOVO_Class",
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
        ]

        keep = [x for x in keep if x in data.columns]

        columns2keep = st.sidebar.multiselect(
            "Choose the columns to keep", data.columns, default=keep
        )

        clinvar_keep = st.multiselect(
            "Clinvar class to keep",
            data["ClinVar_VCF_CLNSIG"].unique(),
        )

        cancervar_keep = st.multiselect(
            "Cancervar class to keep",
            data[" CancerVar: CancerVar and Evidence "].unique(),
        )

        out = out[
            (out[" CancerVar: CancerVar and Evidence "].isin(cancervar_keep))
            | (out["ClinVar_VCF_CLNSIG"].isin(clinvar_keep))
        ]

        vaf = float(st.text_input("VAF threshold", "0"))

        out = out.loc[(out["tumor_f"] > vaf), columns2keep]

        gb = GridOptionsBuilder.from_dataframe(out)
        gb.configure_default_column(
            enablePivot=True, enableValue=True, enableRowGroup=True
        )
        gb.configure_selection(selection_mode="multiple", use_checkbox=True)
        gb.configure_side_bar()
        gridoptions = gb.build()

        response = AgGrid(
            out,
            height=200,
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
