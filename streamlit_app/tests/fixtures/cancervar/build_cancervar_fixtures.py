"""Build the two CancerVar evidence fixtures, and record their provenance.

WHY THIS EXISTS AT ALL, GIVEN THE FILES WERE ALREADY HAND-WRITTEN. These two MAFs were
authored by hand under issue #189 and they are genuinely constructed — round coordinates
or public hotspot positions, invented read counts, no barcode or patient identifier
anywhere. What they lacked was any *recorded* claim to that effect, and
the public export refuses to export a MAF that is not recorded as
``"provenance": "generator-constructed"`` in a ``MANIFEST.json`` beside it.

That refusal is not bureaucracy. Five of the seven parity fixtures this repo used to ship
passed every content pattern in the export gate while carrying real variant calls, because
the bytes that must not travel are, by content, indistinguishable from the bytes that must.
Only a claim about where a file came from separates them.

So rather than assert the claim over files nothing can reproduce, the files are emitted
from the table below. The claim then means what it says, and ``--check`` can hold it: the
committed bytes are the ones this script writes, and the manifest's sha256 is the one the
export gate will recompute. The emitted files are byte-identical to the hand-written
originals — this is a graduation of the provenance, not a rewrite of the fixtures.

The row-by-row rationale — which vector, why that tier, how many real rows are in each
state — is in ``README.md`` and is not repeated here.

Usage::

    python streamlit_app/tests/fixtures/cancervar/build_cancervar_fixtures.py
    python streamlit_app/tests/fixtures/cancervar/build_cancervar_fixtures.py --check
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent

#: The one provenance value the public export accepts for a MAF.
CONSTRUCTED = "generator-constructed"

#: The padded spelling, on 121 of the 124 real files that carry the column. The leading
#: and trailing spaces are in the **name and every value**, and they are load-bearing:
#: ``components/variant_detail._evidence_column`` matches by substring and so survives them
#: by luck rather than by design. The column sits on an interior position in both files so
#: no whitespace-trimming tool can quietly repair it.
PADDED = " CancerVar: CancerVar and Evidence "

#: The R-mangled spelling, on the other 3. A file in this state has **no** ``CancerVar``
#: column at all, which is where 15 of the 109 real files carrying an evidence vector sit.
MANGLED = "CancerVar..CancerVar.and.Evidence"

EVIDENCE_COLUMNS = [
    "Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Reference_Allele",
    "Tumor_Seq_Allele2", "Variant_Classification", "Variant_Type", "Protein_Change",
    "CancerVar", PADDED, "t_alt_count", "t_ref_count",
]

EVIDENCE_ROWS = [
    ["PRDM16", "1", "3000001", "3000001", "C", "T", "Silent", "SNP", "p.S148S",
     "Tier_IV_benign",
     " CancerVar: 2#Tier_IV_benign EVS=[0, 0, 0, 0, 0, 0, 0, -1, 2, 0, 1, 0] ", "12", "88"],
    ["SNRNP200", "2", "96000001", "96000001", "G", "A", "Silent", "SNP", "p.R427R",
     "Tier_IV_benign",
     " CancerVar: 0#Tier_IV_benign EVS=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] ", "9", "91"],
    ["AGRN", "1", "1000001", "1000001", "G", "A", "Missense_Mutation", "SNP", "p.A1514T",
     "Tier_IV_benign",
     " CancerVar: -2#Tier_IV_benign EVS=[0, 0, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0] ", "7", "93"],
    ["MPL", "1", "43000001", "43000001", "A", "T", "Missense_Mutation", "SNP", "p.E576V",
     "Tier_III_Uncertain",
     " CancerVar: 3#Tier_III_Uncertain EVS=[0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0] ", "21", "79"],
    ["PIK3CA", "3", "179000001", "179000001", "A", "G", "Missense_Mutation", "SNP",
     "p.M1043V", "Tier_II_potential",
     " CancerVar: 8#Tier_II_potential EVS=[1, 0, 2, 1, 0, 0, 0, 1, 2, -1, 1, 1] ", "34", "66"],
    ["PIK3CA", "3", "179000002", "179000002", "A", "G", "Missense_Mutation", "SNP",
     "p.E545G", "Tier_I_strong",
     " CancerVar: 13#Tier_I_strong EVS=[2, 0, 2, 1, 0, 0, 1, 1, 2, 1, 1, 2] ", "40", "60"],
    ["KRAS", "12", "25000001", "25000001", "C", "T", "Missense_Mutation", "SNP", "p.G12D",
     "Tier_IV_benign", " CancerVar: 1#Tier_IV_benign ", "18", "82"],
    ["BRAF", "7", "140000001", "140000001", "A", "T", "Missense_Mutation", "SNP",
     "p.V600E", ".", ".", "25", "75"],
]

DOTTED_COLUMNS = [
    "Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Reference_Allele",
    "Tumor_Seq_Allele2", "Variant_Classification", "Variant_Type", "Protein_Change",
    MANGLED, "ClinVar_VCF_CLNSIG", "t_alt_count", "t_ref_count",
]

DOTTED_ROWS = [
    ["TP53", "17", "7577120", "7577120", "C", "T", "Missense_Mutation", "SNP", "p.T304A",
     " CancerVar: 5#Tier_III_Uncertain EVS=[1, 0, 1, 1, 0, 0, 0, 0, 1, -1, 1, 1] ",
     "Uncertain_significance", "30", "70"],
    ["PIK3CA", "3", "179000002", "179000002", "A", "G", "Missense_Mutation", "SNP",
     "p.E545G", " CancerVar: 13#Tier_I_strong EVS=[2, 0, 2, 1, 0, 0, 1, 1, 2, 1, 1, 2] ",
     "Pathogenic", "44", "56"],
    ["KRAS", "12", "25398284", "25398284", "C", "T", "Missense_Mutation", "SNP", "p.G12D",
     " CancerVar: 1#Tier_IV_benign ", "Uncertain_significance", "25", "75"],
    ["EGFR", "7", "55249071", "55249071", "C", "T", "Missense_Mutation", "SNP", "p.T790M",
     ".", "Uncertain_significance", "18", "82"],
    ["BRAF", "7", "140453136", "140453136", "A", "T", "Missense_Mutation", "SNP",
     "p.V600E",
     " CancerVar: 0#Tier_IV_benign EVS=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] ", "Benign",
     "20", "80"],
]

FIXTURES = {
    "somatic_cancervar_evidence.maf": (EVIDENCE_COLUMNS, EVIDENCE_ROWS),
    "somatic_cancervar_dotted_no_column.maf": (DOTTED_COLUMNS, DOTTED_ROWS),
}


def render(columns: list[str], rows: list[list[str]]) -> str:
    return "\n".join(["\t".join(columns)] + ["\t".join(row) for row in rows]) + "\n"


def build(out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    fixtures = {}
    for name, (columns, rows) in FIXTURES.items():
        text = render(columns, rows)
        (out_dir / name).write_text(text)
        fixtures[name] = {
            "rows": len(rows),
            "columns": len(columns),
            "bytes": len(text.encode()),
            "provenance": CONSTRUCTED,
            "sha256": hashlib.sha256(text.encode()).hexdigest(),
            "kind": "cancervar evidence",
        }
    (out_dir / "MANIFEST.json").write_text(
        json.dumps({"fixtures": dict(sorted(fixtures.items()))}, indent=1) + "\n"
    )


def check(committed: Path) -> int:
    differing = []
    for name, (columns, rows) in FIXTURES.items():
        path = committed / name
        if not path.exists() or path.read_text() != render(columns, rows):
            differing.append(name)
    manifest = committed / "MANIFEST.json"
    if not manifest.exists():
        differing.append("MANIFEST.json")
    else:
        recorded = json.loads(manifest.read_text())["fixtures"]
        for name in FIXTURES:
            entry = recorded.get(name, {})
            actual = hashlib.sha256((committed / name).read_bytes()).hexdigest()
            if entry.get("provenance") != CONSTRUCTED or entry.get("sha256") != actual:
                differing.append(f"MANIFEST.json:{name}")
    if differing:
        print("the committed fixtures do not match what this script produces:")
        for name in differing:
            print(f"  {name}")
        return 1
    print(f"{len(FIXTURES)} fixture(s) and their manifest match")
    return 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", type=Path, default=HERE)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args(argv)
    if args.check:
        return check(args.out)
    build(args.out)
    print(f"wrote {len(FIXTURES)} fixture(s) and MANIFEST.json to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
