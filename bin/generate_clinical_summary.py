#!/usr/bin/env python

import argparse
import os
import sys
import pandas as pd

from numpy import append
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
        help="Input Maf file",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="output.maf",
        help="Output file name, default output.maf",
    )
    args = parser.parse_args()
    return args

# Clinvar tags description: https://www.ncbi.nlm.nih.gov/clinvar/docs/clinsig/

CLINICAL_VALUE_MAPPING = {
    # -------------------
    # CancerVar
    # -------------------
    "Tier_I_strong": 1,
    "Tier_II_potential": 2,
    "Tier_III_Uncertain": 3,
    "Tier_IV_benign": 5,

    # -------------------
    # InterVar
    # -------------------
    "Pathogenic": 1,
    "Likely pathogenic": 2,
    "Uncertain significance": 3,
    "Likely benign": 4,
    "Benign": 5,

    # -------------------
    # ClinVar
    # -------------------
    "Pathogenic": 1,
    "Pathogenic/Likely_pathogenic": 1,
    "Likely_pathogenic": 2,
    "Uncertain_significance": 3,
    "Conflicting_classifications_of_pathogenicity": 3,
    "Likely_benign": 4,
    "Benign/Likely_benign": 4,
    "Benign": 5,
    # "_low_penetrance": 1 stood here, and is gone (issue #296). It was an attempt at the
    # same problem main_classification now solves: ClinVar writes low penetrance as a
    # modifier after a comma, never as a call of its own, so the entry could only have
    # fired on a whole cell spelled exactly "_low_penetrance". Measured over 159,573
    # annotated rows (issue #295) it fired on 0, under every candidate rule. A pathogenic
    # call always precedes the modifier, so the split picks the call and the modifier is
    # never what reaches this lookup. Note this dict cannot represent the distinction
    # anyway -- it has no Pathogenic,_low_penetrance key where the app has one ranked just
    # below Pathogenic, and seven tiers where the app has thirteen classes. Reading the
    # call and dropping the qualifier is what the tiers can say.
    "Established_risk_allele": 1,
    "Likely_risk_allele": 2,
    "Uncertain_risk_allele": 3,

    # ClinVar non-patogenic clinical states
    "drug_response": 6,
    "not_provided": 6,
    "": 6,

    # -------------------
    # No clear classification
    # -------------------
    "no_classifications_from_unflagged_records": 7,
    "no_classification_for_the_single_variant": 7,
    "association": 7,
    "other": 7,
    "risk_factor": 7,
    "Affects": 7,
    "protective": 7,
    "association_not_found": 7,
    "confers_sensitivity": 7,

    # -------------------
    # RENOVO
    # -------------------
    "HP Pathogenic": 1,
    "IP Pathogenic": 2,
    "LP Pathogenic": 3,
    "LP Benign": 3,      # borderline / conservative
    "IP Benign": 4,
    "HP Benign": 5,

    # -------------------
    # CIViC
    # -------------------
    "A": 1,
    "B": 2,
    "C": 3,
    "D": 3,
    "E": 5,
}

NUMERIC_TO_CLINICAL = {
    1: "Pathogenic",
    2: "Likely_Pathogenic",
    3: "Uncertain_Significance",
    4: "Likely_Benign",
    5: "Benign",
    6: "Not_Provided",
    7: "Unknown",
}

CLINICAL_COLUMNS = {
    "CancerVar": 1,
    "InterVar": 2,
    "ClinVar_VCF_CLNSIG": 3,
    "RENOVO_Class": 4,
    "CIViC_Evidence_Level": 5,
}



CLINVAR_SEPARATORS = ["|", "/", ";", ","]


def main_classification(value):
    """The one classification a ClinVar cell is read as.

    A cell can hold several calls joined by ``|``, ``/``, ``;`` or ``,``; this picks the
    first that mentions pathogenicity, else the first piece.

    **The piece it picks is read again when it is not a key of CLINICAL_VALUE_MAPPING.**
    The separators are tried in a fixed order, so splitting on the first one present can
    leave a piece that still holds another: ``Pathogenic/Pathogenic,_low_penetrance|other``
    split on ``|`` yields ``Pathogenic/Pathogenic,_low_penetrance``, which was then looked
    up whole, matched nothing, and defaulted to tier 7 -- ``Unknown``, this file's *worst*
    tier. Because the verdict is ``min()``, that failure never surfaced: a second source's
    weaker call outranked it on the same row. Measured over the real corpus (issue #295)
    it was 8 rows in 8 files reading ``Uncertain_Significance`` where the app reads
    ``Pathogenic``.

    **The CLINICAL_VALUE_MAPPING check is what keeps the repair to the rows it is for.**
    Recursing unconditionally instead was measured too, and is wrong: it re-reads the
    piece ``Benign/Likely_benign`` left by ``Benign/Likely_benign|other``, moving 67 rows
    across 66 files from ``Likely_Benign`` to ``Benign`` -- away from the app *and* away
    from what this script itself recorded in those files. Asking only of a piece the split
    has already produced leaves a whole cell splitting exactly as it did before.

    Ported from the app's ``streamlit_app/components/clinical_summary.py``, which arrived
    at this rule as issue #285. Its ``_ATOMIC_TERMS_WITH_SEPARATOR`` guard is deliberately
    *not* ported: it returns terms like ``Pathogenic,_low_penetrance`` whole, and no such
    key exists here, so the term would map to tier 7 where the comma split gives 1 --
    a regression on paper, though one worth 0 rows of the corpus either way.

    The recursion cannot run away: a separator is removed at each step, so the value
    strictly shortens.

    One definition, one reading. This split used to be written out twice in this file,
    once here and once in ``map_to_numeric``. Two copies of one rule is how the app came
    to hold two readings of a ClinVar cell (issues #99, #104), so they are one function
    now rather than two that happen to agree.
    """
    for sep in CLINVAR_SEPARATORS:
        if sep in value:
            parts = value.split(sep)
            pathogenic_parts = [p.strip() for p in parts if "pathogenic" in p.lower()]
            piece = pathogenic_parts[0] if pathogenic_parts else parts[0].strip()
            # Read again only if the piece cannot be named as it stands -- see above.
            if piece in CLINICAL_VALUE_MAPPING:
                return piece
            return main_classification(piece)
    return value


def generate_clinical_summary(row):
    """Generate clinical significance summary for a variant row."""

    numeric_values = []

    for col, priority in CLINICAL_COLUMNS.items():
        if (
            col in row
            and pd.notna(row[col])
            and str(row[col]).strip() not in ["", "nan", "None", "."]
        ):
            value = str(row[col]).strip()

            # Special handling for ClinVar multiple labels
            if col == "ClinVar_VCF_CLNSIG":
                value = main_classification(value)

            mapped_value = CLINICAL_VALUE_MAPPING.get(value, 7)  # default = Unknown

            numeric_values.append(mapped_value)

    if not numeric_values:
        return "🔍 No Clinical Data"

    # numeric comparison
    best_tier = min(numeric_values)

    return NUMERIC_TO_CLINICAL.get(best_tier, "Unknown")



def map_to_numeric(value, column=None):
    if pd.isna(value):
        return np.nan
    
    value = str(value).strip()
    if value in ["", "nan", "None", "."]:
        return np.nan

    if column == "ClinVar_VCF_CLNSIG":
        value = main_classification(value)

    return CLINICAL_VALUE_MAPPING.get(value, 7)  # default Unknown = 7


def compute_final_numeric(row, columns):
    values = [row[f"{col}_NUM"] for col in columns]
    values = [v for v in values if not pd.isna(v)]
    
    if not values:
        return np.nan
    
    return min(values)


def verify_min(row, columns):
    values = [row[f"{col}_NUM"] for col in columns]
    values = [v for v in values if not pd.isna(v)]
    
    if not values:
        return True
    
    return row["FINAL_TIER_NUM"] == min(values)


def main():
    """Generate clinical significance summary for variants in a MAF file."""
    # Parse input
    args = _parse_args()

    if not os.path.exists(args.maf):
        raise ValueError(f"Maf file {args.maf} does not exist!")

    maf = read_maf(args.maf)

    # Apply clinical summary generation to each row
    maf["variantalker_naive"] = maf.apply(generate_clinical_summary, axis=1)
    maf_check = pd.DataFrame()

    CLINICAL_COLUMNS_PRESENT = [col for col in CLINICAL_COLUMNS if col in maf.columns]
    for col in CLINICAL_COLUMNS_PRESENT:
        maf_check[f"{col}_NUM"] = pd.to_numeric(maf[col].apply(lambda x: map_to_numeric(x, col)), errors='coerce')

    maf_check["FINAL_TIER_NUM"] = maf_check.apply(
        lambda row: compute_final_numeric(row, CLINICAL_COLUMNS_PRESENT), axis=1
    )
    maf_check["CHECK_IS_MIN"] = maf_check.apply(
        lambda row: verify_min(row, CLINICAL_COLUMNS_PRESENT), axis=1
    )
    errors = maf_check[maf_check["CHECK_IS_MIN"] == False]
    print(f"Total number of errors in final tier verification: {len(errors)}")

    if not errors.empty:
        print("WARNING: VARIANTALKER_NAIVE:Some rows have a final tier that does not match the minimum of the individual columns", file=sys.stderr)
        sys.exit(1)

    # Write output MAF with new clinical summary column
    writeheader(args.maf, args.output)
    write_metadata2file(f"Clinical summary generation arguments: {args}", args.output)
    maf.to_csv(args.output, sep="\t", index=False, header=True, mode="a")

if __name__ == "__main__":
    main()