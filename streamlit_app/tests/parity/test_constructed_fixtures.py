"""The properties the constructed fixture set has to keep having (issue #246).

``test_parity.py`` asserts the two sides agree and ``test_mutation_coverage.py`` asserts
the fixtures would notice if they stopped. Neither can see the *shape* of the files, and
shape is what the swap from a real subsample to a constructed set was most likely to lose
quietly: a builder that blanks a column, derives one annotation from another, or lets both
arms' marker columns into one file breaks something real while every verdict still agrees.

So each property below is asserted here, against the committed bytes rather than against
the generator's intentions — a generator that stopped honouring one of these would
otherwise be caught by nothing, since it writes the files the rest of the suite reads.

Cheap and pipeline-free on purpose: this reads the MAFs as text and needs no ``bin/``, so
it runs in a packaged tree and in CI where the parity harness cannot.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

import pytest

STREAMLIT_APP = Path(__file__).resolve().parent.parent.parent
if str(STREAMLIT_APP) not in sys.path:
    sys.path.insert(0, str(STREAMLIT_APP))

FIXTURE_DIR = STREAMLIT_APP / "tests" / "fixtures" / "parity"
GENERATOR = FIXTURE_DIR / "build_parity_fixtures.py"

#: Every MAF in the set, and nothing else.
FIXTURES = [
    "somatic_reference.maf",
    "germline_reference.maf",
    "somatic_synthetic.maf",
    "germline_synthetic.maf",
    "somatic_civic.maf",
    "somatic_gnomad_genome.maf",
    "somatic_dot_numeric.maf",
]

#: The per-arm shape the replaced subsample had: ``(rows, columns, comment lines)``.
REFERENCE_SHAPE = {
    "somatic_reference.maf": (82, 427, 111),
    "germline_reference.maf": (94, 380, 111),
}

#: Arm detection keys on column *presence* and never on values, so a file carrying both
#: marker sets is ambiguous — the exact flaw that got a hand-written sample deleted for
#: "inventing shapes no real MAF has".
SOMATIC_MARKERS = ("CancerVar",)
GERMLINE_MARKERS = ("InterVar", "RENOVO_Class", "RENOVO_pls")

ARM = {
    "somatic_reference.maf": "somatic",
    "germline_reference.maf": "germline",
    "somatic_synthetic.maf": "somatic",
    "germline_synthetic.maf": "germline",
    "somatic_civic.maf": "somatic",
    "somatic_gnomad_genome.maf": "somatic",
    "somatic_dot_numeric.maf": "somatic",
}

GENES_SOMATIC = [
    "APC", "CCDC107", "CUL3", "DROSHA", "ERBB3", "HLA-A",
    "KMT2C", "MAP3K4", "MUC6", "PDS5B", "PMS2", "SCAF4",
]
GENES_GERMLINE = [
    "APOB", "AXIN2", "CDA", "DPYD", "EZH2", "FH",
    "KLK3", "MSH2", "NOTCH3", "PIK3CA", "ROS1", "SRSF2",
]
#: The committed mangling, symbol for symbol. Irregular rather than a rule applied to the
#: symbols, which is the point of keeping it verbatim.
GENES_SOMATIC_MIXED_CASE = [
    "Apc", "ccdc107", "Cul3", "drosha", "Erbb3", "hla-a",
    "Kmt2C", "map3k4", "Muc6", "pds5b", "Pms2", "scaf4",
]

BLANKS = ("", ".")

#: The one provenance value the public export accepts for a MAF. Spelled out here
#: rather than imported, and so is the copy in each generator: the export tool sits outside
#: the app's import path, and a constant that travels by copy is what makes a disagreement
#: between them fail the export loudly instead of passing something through.
CONSTRUCTED = "generator-constructed"


def load(name: str) -> tuple[int, list[str], list[list[str]]]:
    """``(comment line count, header, rows)`` for one fixture, read as text."""
    lines = (FIXTURE_DIR / name).read_text().split("\n")
    index = 0
    while index < len(lines) and lines[index].startswith("#"):
        index += 1
    header = lines[index].split("\t")
    rows = [line.split("\t") for line in lines[index + 1:] if line.strip()]
    return index, header, rows


def column(header: list[str], rows: list[list[str]], wanted: str) -> list[str]:
    at = header.index(wanted)
    return [row[at] if at < len(row) else "" for row in rows]


# ---------------------------------------------------------------------------
# Shape
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", sorted(REFERENCE_SHAPE))
def test_the_reference_files_keep_the_shape_the_subsample_had(name):
    """82 x 427 and 94 x 380, each behind 111 leading comment lines.

    The comment block is part of what the set exercises rather than decoration: ``read_maf``
    counts ``#`` lines instead of assuming a fixed offset, and ``writeheader`` copies the
    whole block into the pipeline's output. A three-line header would leave both untested.
    """
    comments, header, rows = load(name)
    assert (len(rows), len(header), comments) == REFERENCE_SHAPE[name]


@pytest.mark.parametrize("name", FIXTURES)
def test_no_column_is_present_but_entirely_blank(name):
    """A column being present has to stay distinguishable from its being absent.

    ``frequency_mask`` treats the two differently — a blank cannot sink a row, while an
    absent column is not in the mask's list at all — and ``config/presets.py`` records the
    measurement where that decided 17 real variants. The first constructed build left 387
    of 427 columns blank, which is not a smaller fixture but a shape in which that
    distinction cannot be witnessed and cannot fail.
    """
    _, header, rows = load(name)
    blank = [
        header[at]
        for at in range(len(header))
        if all(
            (row[at] if at < len(row) else "") in BLANKS
            for row in rows
        )
    ]
    assert not blank, (
        f"{name}: {len(blank)} column(s) carry no value on any row — {blank[:8]}"
    )


@pytest.mark.parametrize("name", sorted(REFERENCE_SHAPE))
def test_both_blanking_conventions_survive(name):
    """Some columns blank as ``.`` and others as ``""``, as the annotator wrote them.

    A builder that picked one would erase the distinction, and it is not cosmetic: ``.``
    leaves an otherwise-numeric column as ``object`` and makes the pipeline raise, while
    empty reads as NaN and makes it drop the row. ``somatic_dot_numeric.maf`` exists
    because those are different outcomes.
    """
    _, header, rows = load(name)
    dotted = sum(
        1 for at in range(len(header))
        if any((row[at] if at < len(row) else "") == "." for row in rows)
    )
    empty = sum(
        1 for at in range(len(header))
        if any((row[at] if at < len(row) else "") == "" for row in rows)
    )
    assert dotted > 0 and empty > 0, (
        f"{name}: {dotted} column(s) blank as '.' and {empty} as '' — one convention has "
        "been lost"
    )


# ---------------------------------------------------------------------------
# Which columns are there, and which deliberately are not
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", FIXTURES)
def test_each_file_carries_exactly_one_arms_marker_columns(name):
    """Arm detection keys on presence, so a file carrying both marker sets is ambiguous."""
    _, header, _ = load(name)
    present = set(header)
    somatic = [m for m in SOMATIC_MARKERS if m in present]
    germline = [m for m in GERMLINE_MARKERS if m in present]
    if ARM[name] == "somatic":
        assert somatic and not germline, f"{name}: germline markers present — {germline}"
    else:
        assert germline and not somatic, f"{name}: somatic markers present — {somatic}"


#: The two column-presence probes, and the columns each one alone may carry. Every other
#: fixture must lack them or the probe proves nothing.
PROBE_COLUMNS = {
    "somatic_civic.maf": {
        "CIViC_Evidence_Level", "CIViC_Evidence_Rating", "CIViC_Entity_Disease",
        "CIViC_Variant_URL", "CIViC_Entity_URL", "CIViC_Entity_Status",
    },
    "somatic_gnomad_genome.maf": {"gnomAD_genome_AF", "gnomAD_genome_AF_raw"},
}


@pytest.mark.parametrize("name", FIXTURES)
def test_the_absent_column_profile_is_preserved(name):
    """No cross-arm classifier column, and each probe's columns nowhere but its own file.

    ``somatic_civic.maf`` and ``somatic_gnomad_genome.maf`` are column-*presence* probes:
    they take a branch by existing, so they only probe anything if the columns they add are
    absent from the other six. Asserted over all seven rather than the two reference files,
    because the file that would break this is one of the other five.
    """
    _, header, _ = load(name)
    present = set(header)
    for probe, columns in PROBE_COLUMNS.items():
        if name == probe:
            # The other direction, and the half that keeps this from going vacuous: a
            # probe whose columns went missing would satisfy every absence below.
            assert columns <= present, (
                f"{name} is the probe for {sorted(columns - present)}, which it no longer "
                "carries — the branch is taken by nothing"
            )
        else:
            assert not columns & present, (
                f"{name} carries {sorted(columns & present)}, so {probe} no longer probes "
                "anything by carrying them"
            )
    if ARM[name] == "somatic":
        assert "InterVar" not in present
        assert not [c for c in header if c.startswith("RENOVO_")]
    else:
        assert "CancerVar" not in present


def test_the_two_alphamissense_columns_disagree_somewhere():
    """They are different annotations, and the app reads exactly one of them.

    ``components.variant_detail`` reads ``am_pathogenicity``; ``AlphaMissense_score`` is
    dbNSFP's and must not be read. A builder deriving one from the other destroys that
    witness in silence — every guard about which column is read still passes, because the
    two agree everywhere.
    """
    _, header, rows = load("germline_reference.maf")
    funcotator = column(header, rows, "am_pathogenicity")
    dbnsfp = column(header, rows, "AlphaMissense_score")
    disagreeing = [
        (a, b)
        for a, b in zip(funcotator, dbnsfp)
        if a not in BLANKS and b not in BLANKS and a != b
    ]
    assert disagreeing, (
        "no germline row has both AlphaMissense columns populated and disagreeing, so "
        "nothing here would notice the app reading the dbNSFP column instead"
    )


def test_no_keep_column_is_uniformly_the_unknown_placeholder():
    """``filter_variants.py:452`` drops any column whose whole value set is ``__UNKNOWN__``.

    It narrows the NOPASS frame with ``~(out == "__UNKNOWN__").all()``, so a KEEP column
    pseudonymised to that one constant takes the run down with
    ``KeyError: ['Matched_Norm_Sample_Barcode'] not in index`` before any verdict. The
    germline tumour barcode is ``__UNKNOWN__`` on every row deliberately — #15's finding
    depends on it — and germline ``keep`` removes exactly that column, which is why it is
    the one column allowed to be uniform.
    """
    from vendor.pipeline_filters import KEEP

    for name in ("somatic_reference.maf", "germline_reference.maf"):
        _, header, rows = load(name)
        # A zero-count sweep, so what it sweeps over is asserted first: with no KEEP column
        # in the file this would pass by having nothing to look at.
        checked = [c for c in header if c in KEEP and c != "Tumor_Sample_Barcode"]
        assert len(checked) > 30 and rows, (
            f"{name}: only {len(checked)} KEEP column(s) and {len(rows)} row(s) to check, "
            "so this sweep is not looking at the file it thinks it is"
        )
        uniform = [
            header[at]
            for at in range(len(header))
            if header[at] in checked
            and all((row[at] if at < len(row) else "") == "__UNKNOWN__" for row in rows)
        ]
        assert not uniform, f"{name}: KEEP column(s) uniformly __UNKNOWN__ — {uniform}"


# ---------------------------------------------------------------------------
# The gene lists, and the manifest's provenance claim
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "name,expected",
    [
        ("genes_somatic.txt", GENES_SOMATIC),
        ("genes_germline.txt", GENES_GERMLINE),
        ("genes_somatic_mixed_case.txt", GENES_SOMATIC_MIXED_CASE),
    ],
)
def test_the_gene_lists_carry_over_verbatim(name, expected):
    """Twelve HGNC symbols each, and the mixed-case file's exact irregular mangling.

    Not patient data — "APC is in this panel" says nothing about any individual — and the
    mangling is the fixture rather than a rule that could be re-derived from the symbols.
    """
    assert (FIXTURE_DIR / name).read_text().split() == expected


@pytest.mark.parametrize("name", FIXTURES)
def test_the_manifest_records_provenance_and_the_committed_bytes(name):
    """What the public export refuses an export without.

    Content patterns provably cannot separate a constructed MAF from a real one — five of
    the seven files this set replaced passed every one of them while carrying real calls —
    so the export asserts a claim about where the file came from instead. The checksum
    matters as much as the flag: without it, renaming a real MAF over a constructed
    fixture's name defeats the whole assertion in one move.
    """
    manifest = json.loads((FIXTURE_DIR / "MANIFEST.json").read_text())
    entry = manifest["fixtures"][name]
    assert entry["provenance"] == CONSTRUCTED
    assert entry["sha256"] == hashlib.sha256((FIXTURE_DIR / name).read_bytes()).hexdigest()


@pytest.mark.parametrize("name", FIXTURES)
def test_the_manifest_names_every_row(name):
    """The per-row oracle is addressed positionally, so an unnamed row is a hole in it.

    ``rows_named[N]`` is row N's constructed-cell name, and the suite quotes it to say
    which shape moved. It went wrong once and silently: ``record_manifest`` read the MAF
    with ``pd.read_csv(comment="#")``, which strips everything after a ``#`` **anywhere on
    a line**, and this set carries ``#`` inside cell values on purpose — CancerVar writes
    its tier into the evidence string as ``-2#Tier_IV_benign``. Three rows lost their
    ``Otherinfo`` and the manifest recorded ``"nan"`` for them, while every other test
    stayed green. Issue #187 had already recorded the same option
    truncating ``CancerVar`` to ``1``.
    """
    manifest = json.loads((FIXTURE_DIR / "MANIFEST.json").read_text())
    named = manifest["fixtures"][name]["rows_named"]
    _, _, rows = load(name)
    assert len(named) == len(rows)
    assert not [n for n in named if n in ("nan", "")], (
        f"{name}: unnamed row(s) at {[i for i, n in enumerate(named) if n in ('nan', '')]}"
    )
    assert len(set(named)) == len(named), (
        f"{name}: duplicate cell names, so the oracle cannot address a row by one"
    )


# ---------------------------------------------------------------------------
# The generator is the artifact, not just a record of one
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_the_generator_reproduces_the_committed_set():
    """``--check`` rebuilds into a temporary directory and compares, file by file.

    The claim that makes this set reproducible in a public checkout rather than merely
    present in it, and the one the generator it replaced could not make: that one needed
    the institutional reference drive, so the committed files were the only artifact.

    Marked ``slow`` because it runs the whole build; everything else in this module reads
    bytes. It needs no ``bin/`` — the pipeline is only involved in ``record_manifest.py``.
    """
    result = subprocess.run(
        [sys.executable, str(GENERATOR), "--check"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
