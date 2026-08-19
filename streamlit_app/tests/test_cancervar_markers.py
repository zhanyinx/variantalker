"""The markers behind a CBP criterion (issue #198): resolution, grouping, and the disclosure.

Four groups of guards:

* **Resolution** — the indices in a MAF's ``Therap_list`` / ``Diag_list`` / ``Prog_list``
  cell against the vendored marker table. Driven against the **real** vendored file, at
  real indices, so a re-vendored table of a different vintage turns these red rather than
  silently renaming drugs. Every expected value here was read out of the file, not invented.
* **The tumour-type split** — the one thing this module computes rather than reads, and the
  reason the disclosure exists in the grouped shape it does: over the real corpus only
  10.6% of therapy markers match the sample's own tissue.
* **Refusals** — an index out of range, an index whose ``Evidence_type`` disagrees with the
  criterion citing it, and a table that is not there. Each is a *named* state (#187), not a
  silent empty list, and the cited count survives all three.
* **The disclosure** — asserted against the HTML, because the whole table is one
  ``st.markdown`` and #188 established that element counts cannot see inside it.

Each guard here was made to fail before being trusted, per this repo's standing rule.
"""

import csv
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from components.cancervar_markers import (  # noqa: E402
    CRITERION_COLUMNS,
    MARKERS_PATH,
    _load,
    resolve_markers,
)
from components.cbp_evidence import (  # noqa: E402
    _cbp_table_html,
    _CBP_CRITERIA,
    _marker_disclosure,
)

# ---------------------------------------------------------------------------
# Real indices, read out of vendor/cancervar_markers.txt
# ---------------------------------------------------------------------------
#
# Stated as constants with their contents spelled out, so that a mismatch names the row it
# expected rather than only the assertion that failed. `test_the_fixture_indices_still_hold
# _what_these_tests_assume` checks every one of them against the file in one place, which
# is what stops the rest of this module going quietly vacuous after a re-vendor.

#: ERBB2 · Breast · Therapeutic · level A. Its `PMIDs` cell is the literal `FDA` — real,
#: and not a PubMed id, which is why the parser filters on digits.
ERBB2_BREAST_A = 3911

#: ABCB1 · Lung · Therapeutic · level D · response `Resistance`. Off-tissue for a breast
#: sample, and an *adverse* response: the case where naming the drug alone inverts the claim.
ABCB1_LUNG_D_RESISTANCE = 1

#: ABCC3 · Breast · Therapeutic · level D · `Resistance`. In-tissue for a breast sample and
#: still adverse — so "in this tumour type" must not read as "an option".
ABCC3_BREAST_D_RESISTANCE = 8

#: ASXL1 · Bone_Marrow · Prognostic · level C. Drug is the literal `N/A` and the response
#: `not applicable` — the table's own way of saying "no drug", which must not reach a screen.
ASXL1_PROGNOSTIC_NO_DRUG = 841

#: ABL1 · Blood;Bone_Marrow · Diagnostic · level A. Drug empty, response `Positive`. The
#: whole of CBP2's content on a real row, and why #198 let CBP2 in with no drug promised.
ABL1_DIAGNOSTIC = 15

#: BRAF · Therapeutic · level C, whose `PMIDs` mixes seven ids with conference abstracts.
BRAF_MIXED_PMIDS = 1443

EXPECTED = {
    ERBB2_BREAST_A: ("ERBB2", "Therapeutic", "A", "Ado-Trastuzumab Emtansine", "Breast"),
    ABCB1_LUNG_D_RESISTANCE: ("ABCB1", "Therapeutic", "D", "Paclitaxel", "Lung"),
    ABCC3_BREAST_D_RESISTANCE: (
        "ABCC3", "Therapeutic", "D", "Monomethyl Auristatin E,Paclitaxel", "Breast",
    ),
    ASXL1_PROGNOSTIC_NO_DRUG: ("ASXL1", "Prognostic", "C", "N/A", "Bone_Marrow"),
    ABL1_DIAGNOSTIC: ("ABL1", "Diagnostic", "A", "", "Blood;Bone_Marrow"),
    BRAF_MIXED_PMIDS: ("BRAF", "Therapeutic", "C", "Dabrafenib", "Lung;Thyroid"),
}


def _row(get_map):
    """A `row.get`-alike over a plain dict, which is all `resolve_markers` needs."""
    return get_map.get


# ---------------------------------------------------------------------------
# The vendored table, and the assumptions these tests rest on
# ---------------------------------------------------------------------------


def test_the_vendored_marker_table_is_present_and_readable():
    """Without it the app can name no marker at all — and the .dmg/.exe carries only this."""
    assert MARKERS_PATH.is_file(), (
        f"{MARKERS_PATH} is missing. It is vendored precisely because the packaged app "
        f"has no resources/CancerVar/ to read.\nFix: python vendor/_sync.py"
    )
    table = _load(str(MARKERS_PATH))
    assert len(table) > 10_000, f"marker table has only {len(table)} rows"


@pytest.mark.parametrize("index", sorted(EXPECTED))
def test_the_fixture_indices_still_hold_what_these_tests_assume(index):
    """Every index this module asserts on, checked against the file in one place.

    The rest of the module reads drugs and levels through `resolve_markers`; if the table
    were re-vendored at a different vintage those assertions would fail with no hint that
    the *indices* moved rather than the code. This names that explicitly — and it is also
    the guard that would catch the off-by-one a stripped header would cause.
    """
    gene, evidence_type, level, drug, cancervar_type = EXPECTED[index]
    row = _load(str(MARKERS_PATH))[index]

    assert row is not None, f"row {index} of the marker table is malformed"
    assert (row[0], row[9], row[10], row[4], row[8]) == (
        gene, evidence_type, level, drug, cancervar_type
    ), (
        f"marker table row {index} is no longer the row these tests were written against.\n"
        f"  expected {(gene, evidence_type, level, drug, cancervar_type)}\n"
        f"  found    {(row[0], row[9], row[10], row[4], row[8])}\n"
        "If the table was re-vendored, the indices in the MAFs on disk moved with it."
    )


def test_row_zero_is_the_header_which_is_what_makes_indices_zero_based():
    """`CancerVar.py:193` reads the file with no header skip, so index 0 is the header.

    Stated here as well as in `test_vendor_drift.py` because *this* module is the one that
    would resolve every index one row off if it stopped being true.
    """
    assert _load(str(MARKERS_PATH))[0][9] == "Evidence_type"


# ---------------------------------------------------------------------------
# Resolution
# ---------------------------------------------------------------------------


def test_a_real_therap_list_cell_resolves_to_its_markers():
    """The leading comma real cells carry is not an index, and duplicates collapse.

    `CancerVar.py:1092` prepends each hit to a string already ending in a comma and then
    de-duplicates through a `set`, so `,3911,1,3911` is what a cell can genuinely look like.
    """
    sets = resolve_markers(_row({"Therap_list": f",{ERBB2_BREAST_A},{ABCB1_LUNG_D_RESISTANCE},"
                                                f"{ERBB2_BREAST_A}"}))

    assert set(sets) == {"CBP1"}
    cbp1 = sets["CBP1"]
    assert cbp1.cited == 2, "the leading comma is not an index and the repeat is not a marker"
    assert cbp1.unresolved == 0
    assert not cbp1.table_missing
    assert [m.drug for m in cbp1.markers] == ["Ado-Trastuzumab Emtansine", "Paclitaxel"], (
        "markers must be ordered strongest evidence level first (A before D)"
    )
    assert cbp1.best_level == "A"


def test_a_cell_with_an_empty_field_in_the_middle_still_resolves():
    """`926,,925` is a real shape — 541 `Prog_list` cells and more `Therap_list` ones.

    `CancerVar.py:1092` builds the string by prepending `str(pos) + ","` and then
    de-duplicates through a `set`, so the one empty field the trailing comma produces can
    land anywhere once `",".join` puts the set back in arbitrary order. Splitting on `,` and
    calling `int` on the pieces would raise on these; scanning for digit runs does not.
    """
    cbp1 = resolve_markers(
        _row({"Therap_list": f"{ERBB2_BREAST_A},,{ABCB1_LUNG_D_RESISTANCE}"})
    )["CBP1"]

    assert cbp1.cited == 2, "an empty field in the middle is not a marker and not an error"
    assert cbp1.unresolved == 0


def test_a_marker_cell_holding_a_chromosome_would_not_be_read_as_an_index():
    """Issue #194 found that an absent column can read as the *chromosome*.

    That hazard does not reach these three columns — measured: 0 of 112,372 rows have a
    marker cell equal to the row's `Chromosome`, because CancerVar *writes* these columns
    rather than resolving them positionally over a dict whose column 0 is `Chr`. But a bare
    `7` is indistinguishable from index 7, so this pins the reason rather than the luck: if
    a future file ever did carry one, it would resolve to a real marker with a real drug.
    Recorded so the next person measures instead of assuming, and so the shape of the
    refusal is visible if one is ever needed.
    """
    # A chromosome-shaped cell resolves, and that is exactly why the measurement matters.
    cbp1 = resolve_markers(_row({"Therap_list": "7"}))["CBP1"]
    assert cbp1.cited == 1, (
        "a bare chromosome would be read as marker index 7 — the measurement that says no "
        "real file carries one is what makes this safe, not the parser"
    )


def test_each_column_resolves_only_its_own_criterion():
    """All three columns at once, each landing on its own code with its own noun."""
    sets = resolve_markers(
        _row({
            "Therap_list": f",{ERBB2_BREAST_A}",
            "Diag_list": f",{ABL1_DIAGNOSTIC}",
            "Prog_list": f",{ASXL1_PROGNOSTIC_NO_DRUG}",
        })
    )

    assert set(sets) == {"CBP1", "CBP2", "CBP3"}
    assert [sets[code].noun for code in ("CBP1", "CBP2", "CBP3")] == [
        "therapy marker", "diagnostic marker", "prognostic marker"
    ]


def test_the_tables_own_no_drug_spellings_do_not_reach_a_screen():
    """`N/A` with `not applicable`, and an empty drug — 1,238 of 10,442 rows between them.

    Handled in this module rather than by widening `config.missing_values`: `N/A` is not a
    MAF sentinel, and `-` is deliberately excluded there because it is a real value in some
    MAF columns. This is the marker *table's* convention, not a MAF's.
    """
    sets = resolve_markers(
        _row({"Prog_list": f",{ASXL1_PROGNOSTIC_NO_DRUG}", "Diag_list": f",{ABL1_DIAGNOSTIC}"})
    )

    prognostic = sets["CBP3"].markers[0]
    assert prognostic.drug == "", "the literal 'N/A' reached the marker as a drug name"
    assert prognostic.response == "", "'not applicable' reached the marker as a response"

    diagnostic = sets["CBP2"].markers[0]
    assert diagnostic.drug == ""
    assert diagnostic.response == "Positive", "a real response must survive"


def test_only_numeric_pmids_survive():
    """Real `PMIDs` cells hold `FDA`, bare `;`, and conference abstracts.

    A renderer that linked these would produce `pubmed.ncbi.nlm.nih.gov/FDA/`.
    """
    sets = resolve_markers(
        _row({"Therap_list": f",{ERBB2_BREAST_A},{BRAF_MIXED_PMIDS}"})
    )
    by_index = {m.index: m for m in sets["CBP1"].markers}

    assert by_index[ERBB2_BREAST_A].pmids == (), "the literal 'FDA' is not a PubMed id"
    braf = by_index[BRAF_MIXED_PMIDS].pmids
    assert braf, "BRAF's real PubMed ids were dropped along with the abstracts"
    assert all(p.isdigit() for p in braf), f"non-numeric ids survived: {braf}"


def test_an_adverse_response_is_flagged():
    """~11k of ~101k resolved rows say Resistance / Resistant / no benefit / Poor Outcome.

    Naming a drug without its response inverts the claim on every one of them.
    """
    sets = resolve_markers(_row({"Therap_list": f",{ABCB1_LUNG_D_RESISTANCE},{ERBB2_BREAST_A}"}))
    by_index = {m.index: m for m in sets["CBP1"].markers}

    assert by_index[ABCB1_LUNG_D_RESISTANCE].is_adverse, "'Resistance' is not a reason to treat"
    assert not by_index[ERBB2_BREAST_A].is_adverse, "'Responsive' must not read as adverse"


# ---------------------------------------------------------------------------
# The tumour-type split
# ---------------------------------------------------------------------------


def test_markers_are_split_by_the_samples_own_tumour_type():
    """CancerVar's own test, replayed: the tissue string as a case-insensitive regex.

    The list is not tissue-filtered by construction — `CancerVar.py:1092` appends before the
    tumour-type test at `:1104`, which only demotes the score — so over the real corpus only
    10.6% of therapy markers match the sample's tissue and 80.1% of cells have none that do.
    """
    cell = f",{ERBB2_BREAST_A},{ABCB1_LUNG_D_RESISTANCE}"
    cbp1 = resolve_markers(_row({"Therap_list": cell}), tissue="Breast")["CBP1"]

    assert cbp1.tissue_known
    assert [m.index for m in cbp1.in_tissue] == [ERBB2_BREAST_A]
    assert [m.index for m in cbp1.other_tissue] == [ABCB1_LUNG_D_RESISTANCE]


def test_an_in_tissue_marker_can_still_be_adverse():
    """The two axes are independent, and conflating them would be the dangerous reading."""
    cbp1 = resolve_markers(
        _row({"Therap_list": f",{ABCC3_BREAST_D_RESISTANCE}"}), tissue="Breast"
    )["CBP1"]

    marker = cbp1.in_tissue[0]
    assert marker.in_tissue and marker.is_adverse, (
        "a resistance marker in this sample's own tumour type must be both"
    )


def test_no_tumour_type_means_unknown_and_not_off_tissue():
    """15 of 109 CancerVar-bearing files carry no `tumor_tissue` column.

    There "not in this tumour type" is a claim the app cannot make, so `in_tissue` is
    `None` — a third state. Reporting `False` would print "none in ..." as a finding.
    """
    cbp1 = resolve_markers(_row({"Therap_list": f",{ERBB2_BREAST_A}"}), tissue="")["CBP1"]

    assert cbp1.markers[0].in_tissue is None
    assert not cbp1.tissue_known
    assert cbp1.in_tissue == (), "an unknown tissue must not populate the in-tissue group"
    assert len(cbp1.other_tissue) == 1, "with no tissue every marker falls in the flat group"


def test_a_tumour_type_that_is_not_a_valid_regex_is_unknown_rather_than_off_tissue():
    """The tissue is used as a *pattern*, so an unbalanced bracket is possible.

    Declining to guess beats reporting every marker as off-tissue, which would read as a
    finding about the variant rather than a limitation of the comparison.
    """
    cbp1 = resolve_markers(_row({"Therap_list": f",{ERBB2_BREAST_A}"}), tissue="Breast[")["CBP1"]

    assert cbp1.markers[0].in_tissue is None
    assert not cbp1.tissue_known


# ---------------------------------------------------------------------------
# Refusals — each one a named state, with the cited count preserved
# ---------------------------------------------------------------------------


def test_an_index_past_the_end_of_the_table_is_unresolved_not_dropped():
    """The count comes from the MAF cell, so it is still knowable."""
    cbp1 = resolve_markers(_row({"Therap_list": f",{ERBB2_BREAST_A},99999999"}))["CBP1"]

    assert cbp1.cited == 2
    assert len(cbp1.markers) == 1
    assert cbp1.unresolved == 1


def test_an_index_whose_evidence_type_disagrees_with_the_criterion_is_refused():
    """This is what a mismatched table vintage looks like, and why it is detectable.

    A diagnostic row cited by `Therap_list` is not a therapy marker. Issue #198 measured 0
    such rows across 112,328 real indices — the check is what makes that measurement mean
    the vendored table is the right vintage, rather than merely that the indices were in
    range.
    """
    cbp1 = resolve_markers(_row({"Therap_list": f",{ABL1_DIAGNOSTIC}"}))["CBP1"]

    assert cbp1.cited == 1
    assert cbp1.markers == ()
    assert cbp1.unresolved == 1
    assert not cbp1.table_missing, "the table is present; it disagrees with the index"


def test_a_missing_table_is_a_different_state_from_a_missing_marker(tmp_path):
    """"No table" and "no markers" get different sentences on screen, so they are distinct."""
    cbp1 = resolve_markers(
        _row({"Therap_list": ",3911,1"}), path=tmp_path / "not_here.txt"
    )["CBP1"]

    assert cbp1.table_missing
    assert cbp1.cited == 2, "the count is read from the MAF, so it survives a missing table"
    assert cbp1.markers == ()


def test_a_criterion_with_no_cell_is_absent_rather_than_empty():
    """Nothing to disclose and something-unnameable-to-disclose are not the same."""
    assert resolve_markers(_row({})) == {}
    assert resolve_markers(_row({"Therap_list": "."})) == {}
    assert resolve_markers(_row({"Therap_list": ""})) == {}


def test_a_malformed_row_does_not_re_base_the_indices_after_it(tmp_path):
    """A short line is kept as a hole, never removed.

    Dropping it would shift every later index by one — the same silent corruption a
    re-vendored table of the wrong vintage causes, arrived at from inside the loader.
    """
    table = tmp_path / "markers.txt"
    header = ["0Gene"] + [f"c{i}" for i in range(1, 11)]
    good = ["GENE", "MUT", "x", "Cancer", "Drug", "sensitive", "A", "1", "Breast",
            "Therapeutic", "A"]
    with table.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(header)
        writer.writerow(["TRUNCATED", "row"])  # index 1, too short
        writer.writerow(good)  # index 2

    cbp1 = resolve_markers(_row({"Therap_list": ",1,2"}), path=table)["CBP1"]

    assert cbp1.unresolved == 1, "the truncated row must be refused, not skipped over"
    assert [m.index for m in cbp1.markers] == [2], (
        "index 2 must still be index 2 — dropping the malformed row would make it index 1"
    )


# ---------------------------------------------------------------------------
# The criterion table
# ---------------------------------------------------------------------------


def test_every_criterion_column_names_a_criterion_the_cbp_table_has():
    """A code here that `_CBP_CRITERIA` does not carry would attach a disclosure to nothing."""
    codes = {c["code"] for c in _CBP_CRITERIA}
    for code, column, evidence_type, noun in CRITERION_COLUMNS:
        assert code in codes, f"{column} claims to sit behind {code}, which is not a criterion"


def test_the_three_marker_columns_are_the_tumour_type_dependent_ones():
    """CBP1/CBP2/CBP3 are the criteria CancerVar backs with a marker list, and only those.

    CBP12 is also tumour-type dependent but is scored from `check_Pubs`, which writes no
    `*_list` column — so there is nothing to disclose there and no column to read.
    """
    assert [code for code, *_ in CRITERION_COLUMNS] == ["CBP1", "CBP2", "CBP3"]
    assert [column for _, column, *_ in CRITERION_COLUMNS] == [
        "Therap_list", "Diag_list", "Prog_list"
    ]


# ---------------------------------------------------------------------------
# The disclosure
# ---------------------------------------------------------------------------


def _disclosure(cell, tissue="Breast", score=1):
    sets = resolve_markers(_row({"Therap_list": cell}), tissue=tissue)
    return _marker_disclosure(sets["CBP1"], tissue, score)


def test_the_disclosure_is_closed_and_states_the_count_before_it_is_opened():
    """Closed is what keeps the section the height #188 measured (median 2 rows).

    And because the tail runs to 143 markers, the summary must let the reader decide
    whether to open it — which is also why the list itself is not truncated.
    """
    html = _disclosure(f",{ERBB2_BREAST_A},{ABCB1_LUNG_D_RESISTANCE}")

    assert "<details" in html and "<summary" in html
    assert " open" not in html.split("<summary")[0], "the disclosure must start closed"
    assert "2 therapy markers" in html
    assert "best level A" in html


def test_the_summary_names_the_samples_tumour_type_and_how_many_matched():
    """The fact that decides whether the rest is worth opening."""
    matched = _disclosure(f",{ERBB2_BREAST_A}", tissue="Breast")
    assert "1 in Breast" in matched

    unmatched = _disclosure(f",{ABCB1_LUNG_D_RESISTANCE}", tissue="Breast")
    assert "none in Breast" in unmatched, (
        "80.1% of populated cells have no in-tissue marker; that must be said, not implied"
    )


def test_the_groups_are_named_and_the_empty_one_is_still_drawn():
    """An in-tissue group of none is the answer on four rows in five, so it is named."""
    html = _disclosure(f",{ABCB1_LUNG_D_RESISTANCE}", tissue="Breast")

    assert "In this sample&#x27;s tumour type (Breast)" in html or (
        "In this sample's tumour type (Breast)" in html
    )
    assert "Other tumour types" in html
    assert "none" in html


def test_the_disclosure_says_the_list_was_never_tissue_filtered():
    """Otherwise a reader takes the off-tissue group for a mistake rather than the design."""
    html = _disclosure(f",{ERBB2_BREAST_A}")
    assert "whatever the tumour type" in html
    assert "only the score is" in html


def test_every_marker_line_carries_its_response():
    """A drug named without its response inverts the claim on the adverse ones."""
    html = _disclosure(f",{ABCB1_LUNG_D_RESISTANCE}")

    assert "Paclitaxel" in html
    assert "Resistance" in html, "the response must never be dropped"
    assert "#b3261e" in html, "an adverse response must be visually distinguished"


def test_a_non_adverse_response_is_not_coloured_as_one():
    html = _disclosure(f",{ERBB2_BREAST_A}")
    assert "Responsive" in html
    assert "#b3261e" not in html


def test_pmids_are_linked_and_the_tail_is_counted():
    """Up to three linked, the rest as a count — p90 is 23 PMIDs and the max is 129."""
    html = _disclosure(f",{BRAF_MIXED_PMIDS}")

    assert "pubmed.ncbi.nlm.nih.gov" in html
    assert html.count("pubmed.ncbi.nlm.nih.gov") <= 3, "only the first three are linked"
    assert "+" in html, "the remaining ids must be counted, not dropped silently"


def test_the_list_is_not_truncated():
    """Groups reach 143 markers, and the reader opened the disclosure on purpose.

    A cap would be a dead end: static HTML has no "show all" to offer, so the alternative
    to a long list is a list the reader cannot finish.
    """
    table = _load(str(MARKERS_PATH))
    many = [i for i, row in enumerate(table) if row and row[9] == "Therapeutic"][:60]
    html = _disclosure("," + ",".join(str(i) for i in many))

    assert html.count("<li") == 60, (
        f"expected all 60 markers rendered, found {html.count('<li')}"
    )


def test_a_missing_table_says_so_rather_than_falling_silent(tmp_path):
    """#187's rule: absence is named per state."""
    sets = resolve_markers(_row({"Therap_list": ",3911,1"}), path=tmp_path / "gone.txt")
    html = _marker_disclosure(sets["CBP1"], "Breast", 1)

    assert "2 therapy markers" in html
    assert "cannot be named" in html
    assert "unavailable" in html


def test_indices_that_all_fail_to_resolve_are_distinguished_from_a_missing_table():
    """A table that disagrees with the indices is a different fact from no table."""
    sets = resolve_markers(_row({"Therap_list": f",{ABL1_DIAGNOSTIC}"}), tissue="Breast")
    html = _marker_disclosure(sets["CBP1"], "Breast", 1)

    assert "did not resolve" in html or "resolved to a" in html
    assert "unavailable" not in html, "the table is present; only the index disagrees"


def test_a_score_of_zero_behind_a_populated_list_is_explained():
    """#185 §4: `level` resets per transcript while the list accumulates over all of them.

    45 `Therap_list` and 103 `Prog_list` cells in the corpus look like this. Unexplained, a
    populated list under a zero score reads as a contradiction.
    """
    explained = _disclosure(f",{ERBB2_BREAST_A}", score=0)
    assert "scored this criterion 0" in explained
    assert "last transcript" in explained

    not_explained = _disclosure(f",{ERBB2_BREAST_A}", score=1)
    assert "last transcript" not in not_explained, (
        "the note must not appear where the score and the list agree"
    )


def test_marker_text_is_html_escaped(tmp_path):
    """Everything on a marker line is text from a data file, rendered with unsafe_allow_html.

    The vendored table happens to carry no `<`, `>` or `&` today — measured — so this is
    driven from a fixture table rather than pretending the real one proves it.
    """
    table = tmp_path / "markers.txt"
    with table.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["0Gene"] + [f"c{i}" for i in range(1, 11)])
        writer.writerow([
            "GENE", "MUT", "x", "<b>cancer</b>", "Drug & <script>alert(1)</script>",
            "sensitive", "A", "1", "Breast", "Therapeutic", "A",
        ])

    sets = resolve_markers(_row({"Therap_list": ",1"}), tissue="Breast", path=table)
    html = _marker_disclosure(sets["CBP1"], "Breast", 1)

    assert "<script>" not in html, "a marker's text reached the page as markup"
    assert "&lt;script&gt;" in html
    assert "Drug &amp; " in html


# ---------------------------------------------------------------------------
# The disclosure inside the table
# ---------------------------------------------------------------------------


def _fired_row(score=1, code="CBP1"):
    criterion = next(c for c in _CBP_CRITERIA if c["code"] == code)
    return [(criterion, score)]


def test_the_disclosure_lands_in_the_criterions_own_cell():
    """Third column, not a fourth: `table-layout:fixed` splits 14/22/64 (#188)."""
    markers = resolve_markers(_row({"Therap_list": f",{ERBB2_BREAST_A}"}), tissue="Breast")
    html = _cbp_table_html(_fired_row(), markers=markers, tissue="Breast")

    assert "<details" in html
    assert html.count("<col ") == 3, "the disclosure must not have added a column"
    # The disclosure sits inside the last <td> of the row, after the sentence it qualifies.
    last_cell = html.rsplit("<td", 1)[1]
    assert "<details" in last_cell
    assert "therapy marker at" in last_cell.split("<details")[0], (
        "the sentence must still come before the disclosure"
    )


def test_a_criterion_with_no_markers_gets_no_disclosure():
    html = _cbp_table_html(_fired_row(), markers={}, tissue="Breast")
    assert "<details" not in html


def test_the_table_still_renders_with_no_markers_argument_at_all():
    """The markers are a detail of a criterion, never a precondition for showing it."""
    html = _cbp_table_html(_fired_row())
    assert "<table" in html
    assert "<details" not in html


def test_the_scored_zero_rows_get_their_disclosures_too():
    """The rows where the score understates what CancerVar found (#185 §4)."""
    markers = resolve_markers(_row({"Therap_list": f",{ERBB2_BREAST_A}"}), tissue="Breast")
    html = _cbp_table_html(_fired_row(score=0), muted=True, markers=markers, tissue="Breast")

    assert "<details" in html, (
        "a populated list behind a zero score is exactly what must not be hidden"
    )
    assert "last transcript" in html


# ---------------------------------------------------------------------------
# The wiring, driven through render_cbp_evidence
# ---------------------------------------------------------------------------
#
# The three tests above assert `_cbp_table_html` *can* draw a disclosure. They say nothing
# about whether `render_cbp_evidence` hands it the markers — and mutation testing proved
# they do not: dropping `markers=markers` from either call site left the whole suite green.
# These are the guards on the wiring itself.


def _render(value, markers=None, tissue=""):
    """Everything `render_cbp_evidence` draws, as (kind, text) pairs.

    Expanders are entered and labelled so a guard can tell "behind the *scored zero*
    expander" from "in the fired table" — a twin of the harness in `test_cbp_evidence.py`,
    kept separate so neither file's edits can quietly weaken the other's.
    """
    from components import cbp_evidence

    drawn = []

    class _Expander:
        def __init__(self, label):
            self.label = label

        def __enter__(self):
            drawn.append(("expander", self.label))
            return self

        def __exit__(self, *exc):
            return False

    class _Fake:
        def markdown(self, text, *a, **k):
            drawn.append(("markdown", str(text)))

        def caption(self, text, *a, **k):
            drawn.append(("caption", str(text)))

        def info(self, text, *a, **k):
            drawn.append(("info", str(text)))

        def expander(self, label, *a, **k):
            return _Expander(label)

    original = cbp_evidence.st
    cbp_evidence.st = _Fake()
    try:
        cbp_evidence.render_cbp_evidence(value, markers=markers, tissue=tissue)
    finally:
        cbp_evidence.st = original
    return drawn


#: A real evidence string with CBP1 = +1 (so CBP1 fires) and CBP2 = 0 (so CBP2 lands in the
#: *scored zero* expander) — which lets one render exercise both call sites.
FIRED_CBP1_ZERO_CBP2 = " CancerVar: 8#Tier_II_potential EVS=[1, 0, 2, 1, 0, 0, 0, 1, 2, -1, 1, 1] "


def _after_expander(drawn, needle="Scored zero"):
    """Everything drawn from the named expander onwards."""
    for position, (kind, text) in enumerate(drawn):
        if kind == "expander" and needle in text:
            return drawn[position:]
    return []


def test_render_cbp_evidence_passes_the_markers_into_the_fired_table():
    """The wiring, not the renderer: CBP1 fires here, so its disclosure must be drawn."""
    markers = resolve_markers(_row({"Therap_list": f",{ERBB2_BREAST_A}"}), tissue="Breast")
    drawn = _render(FIRED_CBP1_ZERO_CBP2, markers=markers, tissue="Breast")

    before = drawn[: len(drawn) - len(_after_expander(drawn)) or None]
    assert any("<details" in text for kind, text in before if kind == "markdown"), (
        "render_cbp_evidence drew the fired table without CBP1's marker disclosure"
    )
    assert any("Ado-Trastuzumab Emtansine" in text for _, text in drawn)


def test_render_cbp_evidence_passes_the_markers_into_the_scored_zero_expander():
    """CBP2 scores 0 in this vector, so its markers must appear behind that expander.

    This is the #185 §4 case as the user meets it: a populated list under a zero score.
    """
    markers = resolve_markers(_row({"Diag_list": f",{ABL1_DIAGNOSTIC}"}), tissue="Breast")
    drawn = _render(FIRED_CBP1_ZERO_CBP2, markers=markers, tissue="Breast")

    behind = _after_expander(drawn)
    assert behind, "the scored-zero expander was not drawn at all"
    assert any("<details" in text for kind, text in behind if kind == "markdown"), (
        "the scored-zero table was drawn without its marker disclosures"
    )
    assert any("last transcript" in text for _, text in behind), (
        "a populated list behind a zero score must be explained where it is shown"
    )


def test_render_cbp_evidence_names_the_tissue_it_was_given():
    """The tissue must reach the summary line, not merely the resolver."""
    markers = resolve_markers(_row({"Therap_list": f",{ABCB1_LUNG_D_RESISTANCE}"}), tissue="Breast")
    drawn = _render(FIRED_CBP1_ZERO_CBP2, markers=markers, tissue="Breast")

    assert any("none in Breast" in text for _, text in drawn), (
        "the sample's tumour type did not reach the disclosure summary"
    )


def test_render_cbp_evidence_still_works_with_no_markers():
    """The section predates the markers and must not depend on them."""
    drawn = _render(FIRED_CBP1_ZERO_CBP2)

    assert any(kind == "markdown" and "<table" in text for kind, text in drawn)
    assert not any("<details" in text for _, text in drawn)
