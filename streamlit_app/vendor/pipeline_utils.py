# Python utility functions

import os
import numpy as np
import re
import pandas as pd

KEEP = [
        "Tumor_Sample_Barcode",
        "Matched_Norm_Sample_Barcode",
        "project_id",
        "Hugo_Symbol",
        "Annotation_Transcript",
        "Chromosome",
        "Start_Position",
        "End_Position",
        "Variant_Classification",
        "Variant_Type",
        "Reference_Allele",
        "Tumor_Seq_Allele1",
        "Tumor_Seq_Allele2",
        "AAChange.refGene",
        "cDNA_Change",
        "Codon_Change",
        "Protein_Change",
        "Transcript_Exon",
        "tumor_f",
        "DP",
        "t_alt_count",
        "t_ref_count",
        "n_alt_count",
        "n_ref_count",
        "variantalker_naive",
        "ClinVar_VCF_CLNSIG",
        "CancerVar",
        "ESCAT",
        "ESCAT_TISSUE",
        "ESCAT_CANCER",
        "CIViC_Evidence_Level",
        "CIViC_Evidence_Rating",
        "CIViC_Entity_Disease",
        "CIViC_Variant_URL",
        "CIViC_Entity_URL",
        "CIViC_Entity_Status",
        "am_class",
        "am_pathogenicity",
        "Otherinfo",
        "tumor_tissue",
        "cosmic",
        "Freq_ExAC_ALL",
        "Freq_esp6500siv2_all",
        "Freq_1000g2015aug_all",
        "filter"
    ]

COLUMNS_TO_REMOVE = [
    "BAM_File",
    "CADD_phred",
    "CADD_raw",
    "CCLE_ONCOMAP_overlapping_mutations",
    "CCLE_ONCOMAP_total_mutations_in_gene",
    "CGC_Chr",
    "CLNALLELEID",
    "CLNDISDB",
    "CLNDN",
    "CLNREVSTAT",
    "CLNSIG",
    "CosmicFusion_fusion_id",
    "Diag_list",
    "End",
    "Freq_gnomAD_genome_ALL",
    "Freq_gnomAD_genome_POPs",
    "Func.refGene",
    "Gene.ensGene",
    "HGNC_Ensembl_ID(supplied_by_Ensembl)",
    "HGNC_Entrez_Gene_ID(supplied_by_NCBI)",
    "HGNC_UCSC_ID(supplied_by_UCSC)",
    "MUTSIG_Published_Results",
    "Match_Norm_Validation_Allele1",
    "Match_Norm_Validation_Allele2",
    "Matched_Norm_Sample_UUID",
    "Mutation_Status",
    "Oreganno_Build",
    "Score",
    "Sequence_Source",
    "Sequencer",
    "Sequencing_Phase",
    "TCGAscape_Amplification_Peaks",
    "TCGAscape_Deletion_Peaks",
    "Tumor_Sample_UUID",
    "Tumor_Validation_Allele1",
    "Tumor_Validation_Allele2",
    "Tumor_type",
    "Tumorscape_Amplification_Peaks",
    "Tumorscape_Deletion_Peaks",
    "UniProt_AApos",
    "UniProt_Experimental_Info",
    "UniProt_Natural_Variations",
    "UniProt_Region",
    "UniProt_Site",
    "Unnamed: 136",
    "Validation_Method",
    "Validation_Status",
    "Verification_Status",
    "avsnp150",
    "dbSNP_RS",
    "dbSNP_Val_Status",
    "phastCons20way_mammalian",
    "strand",
]

CLINVAR_PATHO = [
    "Established_risk_allele",
    "Likely_pathogenic",
    "Likely_pathogenic,_low_penetrance",
    "Likely_pathogenic/Likely_risk_allele",
    "Likely_pathogenic/Pathogenic,_low_penetrance",
    "Likely_pathogenic|Affects",
    "Likely_pathogenic|association",
    "Likely_pathogenic|drug_response",
    "Likely_pathogenic|other",
    "Likely_pathogenic|risk_factor",
    "Likely_risk_allele",
    "Pathogenic",
    "Pathogenic/Likely_pathogenic",
    "Pathogenic/Likely_pathogenic|association",
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
    with open(file, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                c = c + 1

    maf = pd.read_csv(file, header=c, sep="\t", low_memory=False)
    return maf


def writeheader(file: str, outfile: str):
    """Write header from input maf file into outfile."""
    out = open(outfile, "w", encoding="utf-8")
    with open(file, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                out.write(line)

    out.close()


def outer_concat_preserve_order(df1, df2, fill_value=""):
    """
    Outer joins two DataFrames, keeping all columns, filling missing ones with a
    specified value, and preserving the column order of df1 (df2's new columns
    are added at the end).

    Parameters:
        df1 (pd.DataFrame): The base DataFrame (column order is preserved from here).
        df2 (pd.DataFrame): The second DataFrame to merge.
        fill_value (str): Value to fill for missing columns (default is empty string "").

    Returns:
        pd.DataFrame: A new DataFrame with the outer union of df1 and df2.
    """
    # Combine all unique columns, preserving order of df1 first
    all_columns = list(df1.columns) + [
        col for col in df2.columns if col not in df1.columns
    ]

    # Copy to avoid modifying originals
    df1_extended = df1.copy()
    df2_extended = df2.copy()

    # Add missing columns and fill with specified value
    for col in all_columns:
        if col not in df1_extended.columns:
            df1_extended[col] = fill_value
        if col not in df2_extended.columns:
            df2_extended[col] = fill_value

    # Reorder to match all_columns
    df1_extended = df1_extended[all_columns]
    df2_extended = df2_extended[all_columns]

    # Concatenate
    combined = pd.concat([df1_extended, df2_extended], ignore_index=True)
    return combined


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
    out = open(outfile, "a", encoding="utf-8")
    with open(file, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("database_names"):
                out.write(f"## Annovar {line}")
    out.close()


def write_metadata2file(text: str, outfile: str):
    """Append text to file."""
    out = open(outfile, "a", encoding="utf-8")
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

        # if protein change is missing, it's within intron/UTR/splice sites
        try:
            protein_change = row["Protein_Change"].replace("p.", "")
        except:
            continue

        if tissue is not None and tissue != "OTHER":
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

        # no matches
        if len(subset) == 0:
            continue

        # rows with exon-specific annotation
        subset_exon = subset.loc[subset["transcript_exon"] != "NAN"].copy()

        if len(subset_exon) > 0:
            try:
                transcript_exon = int(float(row["Transcript_Exon"]))
            except:
                continue

            subset_exon = subset_exon[
                subset_exon["transcript_exon"].astype(float).astype(int) == transcript_exon
            ].copy()

            if len(subset_exon) == 0:
                continue

            best_score_idx = list(subset_exon["ESCAT"].values)
            best_score_idx = best_score_idx.index(min(best_score_idx))

            maf.loc[idx, "ESCAT"] = subset_exon["ESCAT"].values[best_score_idx]
            maf.loc[idx, "ESCAT_TISSUE"] = subset_exon["tissue"].values[best_score_idx]
            maf.loc[idx, "ESCAT_CANCER"] = subset_exon["cancer"].values[best_score_idx]
            continue

        # generic (non exon-specific) match
        best_score_idx = list(subset["ESCAT"].values)
        best_score_idx = best_score_idx.index(min(best_score_idx))

        maf.loc[idx, "ESCAT"] = subset["ESCAT"].values[best_score_idx]
        maf.loc[idx, "ESCAT_TISSUE"] = subset["tissue"].values[best_score_idx]
        maf.loc[idx, "ESCAT_CANCER"] = subset["cancer"].values[best_score_idx]

    return maf

def parse_annovar_fs_piece(piece):
    """
    Parse a single ANNOVAR-like AAChange annotation string.

    This function extracts structured variant annotation information from a 
    colon-separated ANNOVAR-style annotation entry. The expected format is:

        GENE:TRANSCRIPT:EXON:cDNA_CHANGE:PROTEIN_CHANGE:...

    Example:
        TP53:NM_000546_exon5:c.524_527del:p.Arg175Hisfs*10

    Parameters
    ----------
    piece : str
        A single annotation entry from an AAChange column 
        (e.g., AAChange.refGene, AAChange.ensGene, etc.).

    Returns
    -------
    dict or None
        A dictionary containing:
            - gene (str): Gene symbol.
            - transcript (str): Transcript ID without version suffix.
            - exon (str or None): Exon number extracted from exon field.
            - transcript_position (str or None): cDNA position range (e.g., "524_527").
            - cdna_change (str): Full cDNA change string (e.g., "c.524_527del").
            - protein_change (str): Protein change without trailing stop codon annotation.
        Returns None if parsing fails or the format is invalid.

    Notes
    -----
    - Transcript version suffixes (e.g., "_v1") are removed.
    - Stop codon annotations (e.g., "*10") are stripped from protein change.
    - Exon number is extracted using regex from strings like "exon5".
    """
    parts = piece.split(":")
    if len(parts) < 5:
        return None

    gene = parts[0]
    transcript_full = parts[1]
    exon_part = parts[2]
    cdna_change = parts[3]
    protein_change_raw = parts[4]

    transcript = transcript_full.split("_")[0]

    exon_match = re.search(r'exon(\d+)', exon_part)
    exon = exon_match.group(1) if exon_match else None

    pos_match = re.search(r'c\.(\d+_\d+)', cdna_change)
    transcript_position = pos_match.group(1) if pos_match else None

    protein_change = re.sub(r'\*.*', '', protein_change_raw)

    return {
        "gene": gene,
        "transcript": transcript,
        "exon": exon,
        "transcript_position": transcript_position,
        "cdna_change": cdna_change,
        "protein_change": protein_change
    }


def maf_to_hgvs_genomic(chrom, start, end, ref, alt):
    """
    Convert MAF-like variant information into HGVS genomic notation.

    This function generates a genomic-level HGVS string (g.) based on
    chromosome, position, reference allele, and alternate allele.
    It also classifies the variant type.

    Parameters
    ----------
    chrom : str or int
        Chromosome identifier (with or without 'chr' prefix).
    start : int
        Start genomic position (1-based).
    end : int
        End genomic position (1-based).
    ref : str
        Reference allele.
    alt : str
        Alternate allele.

    Returns
    -------
    tuple
        (hgvs_string, variant_type)

        hgvs_string : str or None
            HGVS genomic notation string (e.g., "g.17:7674221C>T").
        variant_type : str or None
            One of:
                - "snv"     : Single nucleotide variant
                - "del"     : Deletion
                - "ins"     : Insertion
                - "delins"  : Replacement or complex event
            Returns (None, None) if required coordinates are missing.

    Variant classification logic
    ----------------------------
    - SNV: single base substitution.
    - Deletion: alt allele is "-" or empty.
    - Insertion: ref allele is "-" or empty.
    - Delins: length(ref) != length(alt).
    - MNV: same length >1 substitution (treated as delins).
    """
    if pd.isna(chrom) or pd.isna(start) or pd.isna(end):
        return None, None

    chrom_clean = str(chrom).replace("chr", "")
    ref = str(ref)
    alt = str(alt)

    # --- SNV ---
    if len(ref) == 1 and len(alt) == 1:
        hgvs = f"g.{chrom_clean}:{start}{ref}>{alt}"
        return hgvs, "snv"

    # --- Deletion ---
    if alt in {"-", ""}:
        if start == end:
            hgvs = f"g.{chrom_clean}:{start}del"
        else:
            hgvs = f"g.{chrom_clean}:{start}_{end}del"
        return hgvs, "del"

    # --- Insertion ---
    if ref in {"-", ""}:
        hgvs = f"g.{chrom_clean}:{start}_{int(start)+1}ins{alt}"
        return hgvs, "ins"

    # --- Delins ---
    if len(ref) != len(alt):
        hgvs = f"g.{chrom_clean}:{start}_{end}delins{alt}"
        return hgvs, "delins"

    # --- MNV ---
    hgvs = f"g.{chrom_clean}:{start}_{end}{ref}>{alt}"
    return hgvs, "delins"


def fix_funcotator_fs_annotation(df):
    """
    Harmonize and correct frameshift variant annotations from Funcotator output.

    This function:
    1. Identifies the correct AAChange annotation matching the 
       Annotation_Transcript column.
    2. Parses the selected ANNOVAR-like annotation string.
    3. Generates a genomic HGVS string.
    4. Updates standardized variant annotation columns.
    5. Clears fields that cannot be reliably reconstructed.

    Parameters
    ----------
    df : pandas.DataFrame
        Input dataframe containing at least the following columns:
            - Annotation_Transcript
            - Chromosome
            - Start_Position
            - End_Position
            - Reference_Allele
            - Tumor_Seq_Allele2
            - AAChange.ensGene (optional)
            - AAChange.knownGene (optional)
            - AAChange.refGene (optional)

    Returns
    -------
    pandas.DataFrame
        The updated dataframe with standardized annotation fields:
            - Genome_Change
            - Transcript_Exon
            - cDNA_Change
            - Protein_Change
            - Transcript_Position (cleared)
            - Codon_Change (cleared)
            - Other_Transcripts (cleared)

    Notes
    -----
    - Only the annotation matching the primary transcript is used.
    - If no matching transcript is found, the row is skipped.
    - Certain fields are explicitly cleared because no 
      verifiable reconstruction is possible.
    - Designed specifically for frameshift (FS) variant correction workflow.
    """
    aa_columns = [
        "AAChange.ensGene",
        "AAChange.knownGene",
        "AAChange.refGene"
    ]

    for idx, row in df.iterrows():
        transcript_id = str(row.get("Annotation_Transcript", "")).split("_")[0]
    
        # --- NEW: check if all AAChange columns are empty / '.' ---
        aa_values = [
            row.get(col) for col in aa_columns
        ]
    
        if all(
            pd.isna(v) or str(v).strip() in {"", "."}
            for v in aa_values
        ):
            print(f"For variant {row.get('Hugo_Symbol')} : {row.get('Chromosome')} : {row.get('Start_Position')} AAChange columns empty — skipping annotation fix")
            continue
    
        # --- Search matching transcript ---
        match_piece = None
        for col in aa_columns:
            value = row.get(col)
    
            if pd.isna(value) or str(value).strip() in {"", "."}:
                continue
    
            match_piece = next(
                (p for p in str(value).split(",") if transcript_id in p),
                None
            )
            if match_piece:
                break
    
        if not match_piece:
            print(f"Error, transcript {transcript_id} not found in AAChange columns")
            continue

        parsed = parse_annovar_fs_piece(match_piece)
        if not parsed:
            print("Error, parsing failed")
            continue

        genome_change = maf_to_hgvs_genomic(
            row.get("Chromosome"),
            row.get("Start_Position"),
            row.get("End_Position"),
            row.get("Reference_Allele"),
            row.get("Tumor_Seq_Allele2")
        )

        ref = row.get("Reference_Allele") or ""
        alt = row.get("Tumor_Seq_Allele2") or ""
        
        ref_len = 0 if ref == "-" else len(ref)
        alt_len = 0 if alt == "-" else len(alt)
        
        if (ref_len - alt_len) % 3 == 0:
            variant_classification = "In_Frame_Indel"
        else:
            variant_classification = "Frame_Shift_Indel"
        # print(f"ref len = {ref}, alt len = {alt} so {variant_classification}")
        
        string_cols = [
            "Genome_Change",
            "cDNA_Change",
            "Protein_Change",
            "Variant_Classification",
            "Transcript_Position",
            "Codon_Change",
            "Other_Transcripts"
        ]

        for col in string_cols:
            if col in df.columns:
                df[col] = df[col].fillna("").astype(str)

        df.at[idx, "Genome_Change"] = genome_change[0]
        df.at[idx, "Transcript_Exon"] = int(parsed["exon"]) if parsed["exon"] else np.nan
        df.at[idx, "cDNA_Change"] = parsed["cdna_change"]
        df.at[idx, "Protein_Change"] = parsed["protein_change"]
        df.at[idx, "Variant_Classification"] = variant_classification
        df.at[idx, "Transcript_Position"] = ""
        df.at[idx, "Codon_Change"] = ""
        df.at[idx, "Other_Transcripts"] = ""
        
    return df

def merge_with_fallback(
    left: pd.DataFrame,
    right: pd.DataFrame,
    cols_to_use: list,
    left_on: list,
    right_on: list,
) -> pd.DataFrame:
    """
    Merges two dataframes with a primary exact merge and a positional fallback.
 
    Parameters
    ----------
    left : pd.DataFrame
        Left (base) dataframe, e.g. a MAF file.
    right : pd.DataFrame
        Right (annotation) dataframe, e.g. CancerVar, ClinVar, etc.
    cols_to_use : list
        Columns from `right` to include in the merge.
    left_on : list of 5 elements
        Columns in `left` for the exact match, in this order:
        [chromosome, start_position, reference_allele, alt_allele_2, gene]
        e.g. ["Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2", "Ref.Gene"]
    right_on : list of 5 elements
        Columns in `right` for the exact match, in this order:
        [chromosome, start_position, reference_allele, alt_allele, gene]
        e.g. ["Chr", "Start", "Ref", "Alt", "Gene"]
 
    Returns
    -------
    pd.DataFrame
        Merged dataframe with annotations from `right`.
    """
 
    # Unpack expected 4 columns (chrom, start, ref, alt)
    left_chrom_col, left_start_col, left_ref_col, left_alt2_col = left_on[:4]
    right_chrom_col, right_start_col, right_ref_col, right_alt2_col = right_on[:4]
 
    # Remove rows containing 0 in any of the join key columns in left (if numeric) to avoid false matches
    right = right.loc[~right[right_on].eq('0').any(axis=1)].copy()
    right = right.loc[~right[right_on].eq(0).any(axis=1)].copy()

    # --- Primary merge (exact match) ---
    merged = (
        left.merge(
            right[cols_to_use],
            how="left",
            left_on=left_on,
            right_on=right_on,
            indicator=True,
        )
        .drop(columns=list(right_on), errors="ignore")
    )
 
    matched = merged.loc[merged["_merge"].eq("both")].drop(columns="_merge")
    unmatched_left = merged.loc[merged["_merge"].eq("left_only")].drop(columns="_merge")
 
    # annotation columns from right that are carried into result (excluding join keys)
    annotation_cols = [c for c in cols_to_use if c not in right_on]
 
    # If no unmatched left, nothing to attach in fallback (we can't find "closest on left")
    if unmatched_left.empty:
        return matched.reset_index(drop=True)
 
    # --- Find RIGHT rows that did NOT merge in the primary exact join ---
    # Do right->left merge on exact keys and pick right_only
    right_mark = right[cols_to_use].merge(
        left[left_on].drop_duplicates(),
        how="left",
        left_on=right_on,
        right_on=left_on,
        indicator=True,
    )
    
    right_unmatched = right_mark.loc[right_mark["_merge"].eq("left_only"), cols_to_use].copy()

 
    # If nothing unmatched on right, return primary merged result
    if right_unmatched.empty:
        return pd.concat([matched, unmatched_left], ignore_index=True)
 
    # --- Fallback: loop over right_unmatched, find closest in unmatched_left ---
    # To avoid reusing the same unmatched_left row multiple times (optional),
    # we can drop used rows as we go. If you WANT reuse, remove the "drop" lines.
    pool_left = unmatched_left.copy()
    fallback_rows = []
 
    for _, r in right_unmatched.iterrows():
        if pool_left.empty:
            break
 
        # Prefer same chromosome, otherwise use all remaining unmatched left rows
        same_chr = pool_left[pool_left[left_chrom_col] == r[right_chrom_col]]
        search_pool = same_chr if not same_chr.empty else pool_left
 
        # closest by absolute distance on start positions
        closest_idx = (search_pool[left_start_col] - r[right_start_col]).abs().idxmin()
        lrow = pool_left.loc[closest_idx].copy()
 
        # overwrite variant identity with right Start/Ref/Alt
        orig_ref = lrow[left_ref_col]
        lrow[left_start_col] = r[right_start_col]
        lrow[left_ref_col]   = r[right_ref_col]
        lrow[left_alt2_col]  = r[right_alt2_col]
 
        # update Tumor_Seq_Allele1 if present
        if "Tumor_Seq_Allele1" in lrow.index:
            lrow["Tumor_Seq_Allele1"] = (
                r[right_ref_col]
                if lrow["Tumor_Seq_Allele1"] == orig_ref
                else r[right_alt2_col]
            )
 
        # attach annotations
        for col in annotation_cols:
            lrow[col] = r[col]
 
        fallback_rows.append(lrow)
 
        # OPTIONAL: prevent reusing the same left row for multiple right rows
        pool_left = pool_left.drop(index=closest_idx)
 
    fallback_df = pd.DataFrame(fallback_rows)

    # fix funcotator annotation
    fallback_df = fix_funcotator_fs_annotation(fallback_df)

    # Return: matched (primary hits) + fallback rows + remaining unmatched_left (still unannotated)
    # If you do NOT want to keep the leftover unmatched-left rows, remove the last concat.
    out = pd.concat([matched, fallback_df, pool_left], ignore_index=True)
 
    return out

def read_header(file, string):
    """
    Reads the header of a MAF/VCF file (lines starting with '#')
    and returns it as a list of strings.
    """
    header_lines = []
    with open(file, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(string):
                header_lines.append(line.rstrip("\n"))
            else:
                # stop as soon as the first non-header line is encountered
                break
    return header_lines

def extract_contig_order(header_lines):
    """
    Extracts the contig order from a MAF/VCF header.
    Preserves the original order in which they appear.
    
    Works with both:
    '## contig=<ID=chr1,...>'  (space after ##)
    '##contig=<ID=chr1,...>'   (no space)
    """
    contig_order = []
    # allow optional whitespace after ##
    contig_pattern = re.compile(r'##\s*contig=<ID=([^,>]+)')  

    for line in header_lines:
        match = contig_pattern.match(line)
        if match:
            contig_order.append(match.group(1))
    return contig_order

def sort_maf(df: pd.DataFrame, contig_order) -> pd.DataFrame:
    """
    Sorts a MAF dataframe by:
    1) provided contig order
    2) ascending Start_Position
    3) ascending End_Position
    """
    chrom_order = {c: i for i, c in enumerate(contig_order)}

    df["__chr_order"] = (
        df["Chromosome"]
        .astype(str)
        .map(chrom_order)
        .fillna(len(chrom_order))
        .astype(int)
    )

    for col in ("Start_Position", "End_Position"):
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df.sort_values(
        ["__chr_order", "Start_Position", "End_Position"],
        kind="mergesort",
        inplace=True,
    )

    df.drop(columns="__chr_order", inplace=True)
    df.reset_index(drop=True, inplace=True)

    return df

def read_vcf(file: str) -> pd.DataFrame:
    c = 0
    with open(file, encoding="utf-8") as f:
        for line in f:
            if line.startswith("##"):
                c += 1
            else:
                break

    vcf = pd.read_csv(
        file,
        header=c,
        sep="\t",
        low_memory=False,
        dtype=str
    )

    return vcf

def sort_vcf(df: pd.DataFrame, contig_order) -> pd.DataFrame:
    """
    Sorts a VCF dataframe by:
    1) provided contig order (#CHROM)
    2) ascending POS
    3) REF and ALT for deterministic ordering
    """

    df = df.copy()

    # ensure correct data types
    df["#CHROM"] = df["#CHROM"].astype(str)
    df["POS"] = pd.to_numeric(df["POS"], errors="coerce")

    # create ordered categorical chromosome column
    df["#CHROM"] = pd.Categorical(
        df["#CHROM"],
        categories=contig_order,
        ordered=True
    )

    # stable sort (important for reproducibility)
    df_sorted = (
        df.sort_values(
            by=["#CHROM", "POS", "REF", "ALT"],
            kind="mergesort"
        )
        .reset_index(drop=True)
    )

    return df_sorted

def has_clinvar_term(value, clinvar_keep):
    """
    Return True if at least one parsed ClinVar term is present in clinvar_keep.
    Separators handled: | / ; ,
    """

    if pd.isna(value):
        return False

    if not clinvar_keep:
        return False

    separators_regex = r"[|/;,]"

    parsed_terms = [
        term.strip()
        for term in re.split(separators_regex, str(value))
        if term.strip()
    ]

    clinvar_keep_clean = {
        term.strip()
        for term in clinvar_keep
        if term.strip()
    }

    return any(term in clinvar_keep_clean for term in parsed_terms)
