"""The predictor scores behind PP3/BP4 and CBP10, and the section that draws them (issue #190).

Four groups, failing for four different reasons:

* **The readings** — pure, so plain value assertions. The interesting ones are the *absences*:
  the two classifiers treat an unreadable input oppositely, and a row that said "no value" and
  stopped there would be the same sentence for two opposite facts.
* **The arithmetic** — ``_beyond_what_this_file_holds``, which is the one claim the section makes
  about the relationship between the file and the verdict. It is a necessary condition on the
  counts, never a re-derivation; the guards pin both that it fires and that it does not fire on a
  file that can account for its verdict.
* **The section** — what a clinician sees, driven through a fake ``st`` and asserted against the
  HTML, because each table is one ``st.markdown`` and element counts cannot see inside it (#188).
* **The panel** — that each table attaches under its own arm's evidence section and neither
  appears on the other arm, driven through ``render_variant_detail_panel`` over the real fixture
  MAFs so the arm gate is exercised as bytes on disk.

Each guard here was made to fail before being trusted, per this repo's standing rule; the harness
is issue #190's mutation script.
"""

import os
import sys

import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from components.predictor_context import (  # noqa: E402
    CBP10_INDEX,
    _ABSENT,
    _BLANK,
    _beyond_what_this_file_holds,
    _cbp10_readings,
    _looks_synonymous,
    _points_from_absence,
    _pp3_bp4_readings,
    _tally,
    render_cbp10_inputs,
    render_pp3_bp4_inputs,
)

CANCERVAR_FIXTURES = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "fixtures", "cancervar"
)

#: A real somatic evidence string, padding included — the same one ``test_cbp_evidence`` uses.
#: ``EVS[9]``, CBP10, is the tenth entry: ``-1``.
REAL_SOMATIC = " CancerVar: 8#Tier_II_potential EVS=[1, 0, 2, 1, 0, 0, 0, 1, 2, -1, 1, 1] "

#: The same vector with CBP10 at ``+1`` and at ``+2``, for the accountability guards.
SOMATIC_CBP10_PLUS_ONE = "CancerVar: 8#Tier_II_potential EVS=[1, 0, 2, 1, 0, 0, 0, 1, 2, 1, 1, 1]"
SOMATIC_CBP10_PLUS_TWO = "CancerVar: 9#Tier_II_potential EVS=[1, 0, 2, 1, 0, 0, 0, 1, 2, 2, 1, 1]"

#: A germline evidence string in InterVar's own shape. BP4 is ``BP`` index 3, PP3 is ``PP`` index 2.
GERMLINE_BP4 = (
    "InterVar: Likely benign PVS1=0 PS=[0, 0, 0, 0] PM=[0, 0, 0, 0, 0, 0] "
    "PP=[0, 0, 0, 0, 0] BA1=0 BS=[0, 0, 0, 0] BP=[0, 0, 0, 1, 0, 0, 0]"
)
GERMLINE_PP3 = (
    "InterVar: Uncertain significance PVS1=0 PS=[0, 0, 0, 0] PM=[0, 0, 0, 0, 0, 0] "
    "PP=[0, 0, 1, 0, 0] BA1=0 BS=[0, 0, 0, 0] BP=[0, 0, 0, 0, 0, 0, 0]"
)
GERMLINE_NEITHER = (
    "InterVar: Uncertain significance PVS1=0 PS=[0, 0, 0, 0] PM=[0, 0, 0, 0, 0, 0] "
    "PP=[0, 0, 0, 0, 0] BA1=0 BS=[0, 0, 0, 0] BP=[0, 0, 0, 0, 0, 0, 0]"
)


def _reading(readings, name):
    for reading in readings:
        if reading.name == name:
            return reading
    raise AssertionError(f"no row named {name!r} in {[r.name for r in readings]}")


# ---------------------------------------------------------------------------
# The readings: CancerVar's five
# ---------------------------------------------------------------------------


def test_the_five_rows_are_the_five_columns_check_prep_reads():
    """Not four, not six, and in ``check_PreP``'s own order."""
    readings = _cbp10_readings(pd.Series(dtype=object), frozenset())
    assert [r.name for r in readings] == [
        "SIFT",
        "PolyPhen-2 HDIV",
        "FATHMM",
        "MutationAssessor",
        "GERP++ RS",
    ]
    assert [r.column for r in readings] == [
        "SIFT_score",
        "Polyphen2_HDIV_pred",
        "FATHMM_pred",
        "MutationAssessor_pred",
        "GERP++_RS",
    ]


@pytest.mark.parametrize(
    ("score", "side", "phrase"),
    [
        ("0.03", "damaging", "below CancerVar's 0.05"),
        ("0.05", "benign", "at or above CancerVar's 0.05"),  # the boundary is CancerVar's `>=`
        ("0.9", "benign", "at or above CancerVar's 0.05"),
    ],
)
def test_sift_is_read_on_cancervars_side_of_its_own_boundary(score, side, phrase):
    """``0.05`` itself is tolerated, because ``check_PreP`` tests ``>= sift_cutoff``."""
    reading = _reading(_cbp10_readings(pd.Series({"SIFT_score": score}), frozenset()), "SIFT")
    assert reading.side == side
    assert phrase in reading.says


@pytest.mark.parametrize(
    ("value", "side"),
    [("4.71", "damaging"), ("2", "damaging"), ("2.0", "damaging"), ("1.99", "benign")],
)
def test_cancervar_calls_gerp_conserved_at_exactly_two(value, side):
    """CancerVar uses ``>=``; InterVar uses ``>``. This is the somatic side of that split."""
    reading = _reading(_cbp10_readings(pd.Series({"GERP++_RS": value}), frozenset()), "GERP++ RS")
    assert reading.side == side


@pytest.mark.parametrize(
    ("column", "name", "letter", "side"),
    [
        ("Polyphen2_HDIV_pred", "PolyPhen-2 HDIV", "D", "damaging"),
        ("Polyphen2_HDIV_pred", "PolyPhen-2 HDIV", "P", "damaging"),
        ("Polyphen2_HDIV_pred", "PolyPhen-2 HDIV", "B", "benign"),
        ("FATHMM_pred", "FATHMM", "D", "damaging"),
        ("FATHMM_pred", "FATHMM", "T", "benign"),
        ("MutationAssessor_pred", "MutationAssessor", "H", "damaging"),
        ("MutationAssessor_pred", "MutationAssessor", "M", "damaging"),
        ("MutationAssessor_pred", "MutationAssessor", "L", "benign"),
        ("MutationAssessor_pred", "MutationAssessor", "N", "benign"),
    ],
)
def test_every_letter_in_the_corpus_is_read_the_way_cancervar_reads_it(column, name, letter, side):
    """The whole measured vocabulary — 60,417 real ``_pred`` cells hold nothing else."""
    reading = _reading(_cbp10_readings(pd.Series({column: letter}), frozenset()), name)
    assert reading.side == side
    assert reading.value == letter


def test_polyphens_p_is_possibly_damaging_and_not_a_polymorphism():
    """Issue #186's finding: ``P`` means different things in different tools' columns.

    PolyPhen-2's ``P`` is *possibly damaging*; MutationTaster's is *polymorphism*, the opposite.
    This section only draws PolyPhen-2's, and it must not borrow the other's word.
    """
    reading = _reading(
        _cbp10_readings(pd.Series({"Polyphen2_HDIV_pred": "P"}), frozenset()), "PolyPhen-2 HDIV"
    )
    assert "possibly damaging" in reading.says
    assert "polymorphism" not in reading.says.lower()


def test_a_letter_cancervar_does_not_recognise_counts_for_neither_side():
    """``check_PreP`` compares with ``==``, so an unknown letter increments nothing at all."""
    reading = _reading(
        _cbp10_readings(pd.Series({"FATHMM_pred": "Z"}), frozenset()), "FATHMM"
    )
    assert reading.side == "neither"
    assert reading.value == "Z"
    assert "recognises" in reading.says


def test_a_multi_valued_pred_cell_is_not_split():
    """CancerVar reads the whole cell, so the app must too.

    Splitting on ``;`` here would report a ``D`` that CancerVar never saw: its ``==`` against
    ``'D;D'`` is false, so the cell counts for neither side. No such cell exists in the corpus,
    which is exactly why the behaviour has to be pinned rather than left to be discovered.
    """
    reading = _reading(
        _cbp10_readings(pd.Series({"Polyphen2_HDIV_pred": "D;D"}), frozenset()),
        "PolyPhen-2 HDIV",
    )
    assert reading.side == "neither"
    assert reading.value == "D;D"


def test_an_absent_column_and_an_empty_cell_are_told_apart():
    """A file-level gap and a variant-level one are different facts and read differently."""
    absent = _reading(_cbp10_readings(pd.Series({"SIFT_score": "0.1"}), frozenset()), "FATHMM")
    blank = _reading(
        _cbp10_readings(pd.Series({"FATHMM_pred": "."}), frozenset()), "FATHMM"
    )
    assert absent.value == _ABSENT
    assert absent.side == "absent"
    assert blank.value == _BLANK
    assert blank.side == "blank"
    assert not absent.readable and not blank.readable


def test_the_polyphen_score_is_never_substituted_for_the_prediction():
    """The measured refusal: 35 of the 52 files missing ``_pred`` hold a chromosome in ``_score``.

    A row that quietly reached for the near-neighbour would print a number on two thirds of the
    files that needed it — and on 35 of them the number is the chromosome.
    """
    row = pd.Series({"Polyphen2_HDIV_score": "0.997", "SIFT_score": "0.01"})
    reading = _reading(_cbp10_readings(row, frozenset()), "PolyPhen-2 HDIV")
    assert reading.value == _ABSENT
    assert "0.997" not in reading.says
    assert reading.column == "Polyphen2_HDIV_pred"


def test_a_contaminated_column_is_drawn_but_is_not_a_reading():
    """Issue #194's rule, carried into this table: show the value, give it no vote."""
    row = pd.Series({"Chromosome": "chr7", "SIFT_score": "7"})
    reading = _reading(_cbp10_readings(row, frozenset({"SIFT_score"})), "SIFT")
    assert reading.value == "⚠ 7"
    assert reading.side == "chromosome"
    assert not reading.readable
    assert "holds the chromosome" in reading.says


def test_a_contaminated_cell_on_a_sex_chromosome_is_drawn_too():
    """The #194 bug that mutation-testing found: ``float('X')`` refuses, and the row went silent."""
    row = pd.Series({"Chromosome": "chrX", "GERP++_RS": "X"})
    reading = _reading(_cbp10_readings(row, frozenset({"GERP++_RS"})), "GERP++ RS")
    assert reading.value == "⚠ X"
    assert reading.side == "chromosome"


def test_a_contaminated_column_absent_from_this_row_is_just_absent():
    """A verdict naming a column this row does not have must not invent a warned row."""
    reading = _reading(
        _cbp10_readings(pd.Series({"SIFT_score": "0.1"}), frozenset({"FATHMM_pred"})), "FATHMM"
    )
    assert reading.value == _ABSENT
    assert reading.side == "absent"


def test_the_dot_mangled_gerp_spelling_is_read_rather_than_called_absent():
    """``GERP.._RS`` on 2 of 167 real MAFs. Read under the file's own spelling of the column."""
    reading = _reading(
        _cbp10_readings(pd.Series({"GERP.._RS": "4.71"}), frozenset()), "GERP++ RS"
    )
    assert reading.value == "4.71"
    assert reading.side == "damaging"


# ---------------------------------------------------------------------------
# The readings: InterVar's four
# ---------------------------------------------------------------------------


def test_the_four_rows_are_the_columns_check_pp3_reads():
    readings = _pp3_bp4_readings(pd.Series(dtype=object), frozenset())
    assert [r.column for r in readings] == [
        "MetaSVM_score",
        "GERP++_RS",
        "dbscSNV_ADA_SCORE",
        "dbscSNV_RF_SCORE",
    ]


@pytest.mark.parametrize(
    ("value", "side", "criterion"),
    [("0.42", "damaging", "PP3"), ("-1.2", "benign", "BP4")],
)
def test_metasvm_names_the_criterion_its_value_pushed_toward(value, side, criterion):
    reading = _reading(
        _pp3_bp4_readings(pd.Series({"MetaSVM_score": value}), frozenset()), "MetaSVM"
    )
    assert reading.side == side
    assert criterion in reading.says


def test_metasvm_at_exactly_zero_is_a_point_for_neither():
    """``check_PP3`` tests ``> 0`` and ``check_BP4`` tests ``< 0``, so ``0`` satisfies neither."""
    reading = _reading(
        _pp3_bp4_readings(pd.Series({"MetaSVM_score": "0"}), frozenset()), "MetaSVM"
    )
    assert reading.side == "neither"
    assert "neither" in reading.says


def test_an_unreadable_metasvm_is_said_to_fire_pp3_rather_than_left_blank():
    """The session's central germline finding, and it is counter-intuitive enough to need saying.

    ``check_PP3``'s ``except ValueError`` branch loops over ``["synon", "coding-synon"]`` and sets
    its point whenever the consequence does *not* contain the token. No ANNOVAR consequence
    contains ``coding-synon``, so the second iteration fires unconditionally: an absent MetaSVM
    score is one of PP3's three points. It happened on 1,278 of the 2,377 real rows where PP3 fired.
    """
    reading = _reading(
        _pp3_bp4_readings(pd.Series({"GERP++_RS": "1.0"}), frozenset()), "MetaSVM"
    )
    assert reading.value == _ABSENT
    assert "PP3" in reading.says
    assert not reading.readable


def test_an_unreadable_metasvm_on_a_synonymous_variant_names_bp4_as_well():
    reading = _reading(
        _pp3_bp4_readings(
            pd.Series({"ExonicFunc.refGene": "synonymous_SNV"}), frozenset()
        ),
        "MetaSVM",
    )
    assert "PP3" in reading.says
    assert "BP4" in reading.says


def test_a_nonsynonymous_consequence_does_not_count_as_synonymous():
    """``check_BP4``'s veto: it looks for ``synon`` and then rules out ``nonsynon``."""
    assert _looks_synonymous(pd.Series({"ExonicFunc.refGene": "synonymous_SNV"}))
    assert not _looks_synonymous(pd.Series({"ExonicFunc.refGene": "nonsynonymous_SNV"}))
    assert not _looks_synonymous(pd.Series({"ExonicFunc.refGene": "."}))
    assert not _looks_synonymous(pd.Series(dtype=object))


def test_intervar_calls_gerp_conserved_only_above_two():
    """The germline side of the split: ``check_PP3`` tests ``>``, so ``2`` is *not* conserved."""
    at_two = _reading(
        _pp3_bp4_readings(pd.Series({"GERP++_RS": "2"}), frozenset()), "GERP++ RS"
    )
    above = _reading(
        _pp3_bp4_readings(pd.Series({"GERP++_RS": "2.01"}), frozenset()), "GERP++ RS"
    )
    assert at_two.side == "benign"
    assert "BP4" in at_two.says
    assert above.side == "damaging"
    assert "PP3" in above.says


def test_an_unreadable_gerp_is_said_to_give_bp4_its_point():
    reading = _reading(_pp3_bp4_readings(pd.Series(dtype=object), frozenset()), "GERP++ RS")
    assert reading.value == _ABSENT
    assert "BP4" in reading.says


@pytest.mark.parametrize(("value", "side"), [("0.71", "damaging"), ("0.6", "benign")])
def test_the_splice_scores_are_read_at_intervars_own_cutoff(value, side):
    reading = _reading(
        _pp3_bp4_readings(pd.Series({"dbscSNV_ADA_SCORE": value}), frozenset()), "dbscSNV ADA"
    )
    assert reading.side == side


def test_an_unreadable_splice_score_is_said_to_give_bp4_its_third_point():
    reading = _reading(_pp3_bp4_readings(pd.Series(dtype=object), frozenset()), "dbscSNV RF")
    assert "BP4" in reading.says
    assert "third point" in reading.says


# ---------------------------------------------------------------------------
# The arithmetic
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("cbp10", "readable", "flagged"),
    [
        (2, 5, False),
        (2, 4, True),   # +2 needs all five damaging
        (1, 3, False),
        (1, 2, True),   # +1 needs more than two
        (-1, 3, False),
        (-1, 2, True),  # -1 needs more than two benign
        (0, 0, False),  # zero needs nothing, so a file with nothing can still hold it
    ],
)
def test_the_accountability_test_is_the_counts_cancervar_grades_on(cbp10, readable, flagged):
    """Verified against a faithful transcription of ``check_PreP`` over 112,299 real rows:
    3,584 flagged, 0 false alarms, 0 misses."""
    assert _beyond_what_this_file_holds(cbp10, readable) is flagged


def test_a_value_cancervar_has_never_emitted_is_not_flagged():
    """``_CBP10_CALLS_NEEDED`` has four keys; a fifth value must not raise or accuse."""
    assert _beyond_what_this_file_holds(7, 0) is False


def test_the_tally_counts_what_the_classifier_could_read_and_names_the_rest():
    """The three counts are distinct on purpose: 2 damaging, 1 benign, 2 unreadable.

    A row with one benign and one unreadable would let ``benign`` be counted as *unreadable* and
    still read 1 — which is exactly the mutation that survived the first pass of issue #190's
    mutation sweep. Every count here has to be its own number.
    """
    readings = _cbp10_readings(
        pd.Series(
            {
                "SIFT_score": "0.01",
                "Polyphen2_HDIV_pred": "B",
                "FATHMM_pred": ".",
                "GERP++_RS": "4.7",
            }
        ),
        frozenset(),
    )
    text = _tally(readings)
    assert "**2** damaging" in text
    assert "**1** benign" in text
    assert "**2** that CancerVar could not read here" in text


def test_a_contaminated_pred_column_is_drawn_but_is_not_a_reading():
    """The letter rows need the same #194 rule as the numeric ones.

    Measured at 0 of 94 somatic and 0 of 57 germline files, so this is prophylaxis — and the
    prophylaxis has to be guarded or it is not there: the first pass of the mutation harness
    removed the check from ``_letter_reading`` and nothing failed, because every contamination
    guard here used a numeric column.
    """
    row = pd.Series({"Chromosome": "chr9", "FATHMM_pred": "9"})
    reading = _reading(_cbp10_readings(row, frozenset({"FATHMM_pred"})), "FATHMM")
    assert reading.value == "⚠ 9"
    assert reading.side == "chromosome"
    assert not reading.readable
    assert "holds the chromosome" in reading.says


def test_bp4s_points_from_absence_treat_the_two_splice_scores_as_one():
    """``check_BP4``'s dbscSNV branch is one ``try``, so one absent column costs the pair."""
    row = pd.Series({"MetaSVM_score": "-1.0", "GERP++_RS": "1.0", "dbscSNV_ADA_SCORE": "0.1"})
    readings = _pp3_bp4_readings(row, frozenset())
    assert _points_from_absence(readings, row) == 1  # only the RF column is missing

    everything_absent = pd.Series({"ExonicFunc.refGene": "synonymous_SNV"})
    assert _points_from_absence(
        _pp3_bp4_readings(everything_absent, frozenset()), everything_absent
    ) == 3


def test_a_nonsynonymous_row_does_not_get_bp4s_first_point_from_absence():
    row = pd.Series({"ExonicFunc.refGene": "nonsynonymous_SNV"})
    assert _points_from_absence(_pp3_bp4_readings(row, frozenset()), row) == 2


# ---------------------------------------------------------------------------
# The sections
# ---------------------------------------------------------------------------


def _drawn(monkeypatch, render, *args):
    """Everything one section writes, as one string."""
    from unittest.mock import MagicMock

    from components import predictor_context

    written = []
    fake = MagicMock()
    fake.markdown.side_effect = lambda text, *a, **k: written.append(str(text))
    fake.caption.side_effect = lambda text, *a, **k: written.append(str(text))
    fake.info.side_effect = lambda text, *a, **k: written.append(str(text))
    monkeypatch.setattr(predictor_context, "st", fake)
    render(*args)
    return "\n".join(written)


def test_the_somatic_section_names_the_criterion_and_attributes_the_cutoffs(monkeypatch):
    row = pd.Series(
        {
            "SIFT_score": "0.01",
            "Polyphen2_HDIV_pred": "D",
            "FATHMM_pred": "D",
            "MutationAssessor_pred": "H",
            "GERP++_RS": "4.7",
        }
    )
    text = _drawn(monkeypatch, render_cbp10_inputs, row, REAL_SOMATIC, frozenset())
    assert "CBP10" in text
    assert "CancerVar recorded **CBP10 = -1**" in text
    assert "cutoffs above are its own, not the guideline's" in text
    assert "0.05" in text  # the SIFT threshold, stated
    assert "PolyPhen-2 HDIV" in text


def test_the_somatic_section_says_when_the_file_cannot_account_for_the_verdict(monkeypatch):
    """The 15 thin files: two readable inputs and a CBP10 of ``+1``, which needs three calls."""
    row = pd.Series({"SIFT_score": "0.01", "GERP++_RS": "4.7"})
    text = _drawn(monkeypatch, render_cbp10_inputs, row, SOMATIC_CBP10_PLUS_ONE, frozenset())
    assert "Only 2 of the five columns" in text
    assert "cannot account for its answer" in text
    assert "still CancerVar's" in text


def test_a_file_with_every_input_gets_no_accountability_warning(monkeypatch):
    row = pd.Series(
        {
            "SIFT_score": "0.01",
            "Polyphen2_HDIV_pred": "D",
            "FATHMM_pred": "D",
            "MutationAssessor_pred": "H",
            "GERP++_RS": "4.7",
        }
    )
    text = _drawn(monkeypatch, render_cbp10_inputs, row, SOMATIC_CBP10_PLUS_TWO, frozenset())
    assert "cannot account" not in text


def test_the_somatic_section_draws_nothing_when_the_evidence_string_does_not_parse(monkeypatch):
    """``render_cbp_evidence`` has already said so directly above; twice is noise."""
    row = pd.Series({"SIFT_score": "0.01"})
    assert _drawn(monkeypatch, render_cbp10_inputs, row, "1#Tier_IV_benign", frozenset()) == ""
    assert _drawn(monkeypatch, render_cbp10_inputs, row, ".", frozenset()) == ""


def test_cbp10_is_read_from_the_tenth_entry_not_the_ninth(monkeypatch):
    """The ticket was charted calling this CBP9. ``EVS[9]`` is zero-based, so it is CBP10.

    The two are different criteria with different value sets — CBP9 is the somatic-database
    criterion and never goes negative — so reading the wrong index would put a somatic-database
    score under five predictor rows.
    """
    assert CBP10_INDEX == 9
    row = pd.Series({"SIFT_score": "0.01"})
    # CBP9 (index 8) is 2 here and CBP10 (index 9) is -1.
    text = _drawn(monkeypatch, render_cbp10_inputs, row, REAL_SOMATIC, frozenset())
    assert "CBP10 = -1" in text
    assert "CBP10 = +2" not in text


def test_the_germline_section_names_which_criterion_intervar_recorded(monkeypatch):
    row = pd.Series({"MetaSVM_score": "-1.2", "GERP++_RS": "1.0",
                     "dbscSNV_ADA_SCORE": "0.1", "dbscSNV_RF_SCORE": "0.1"})
    text = _drawn(monkeypatch, render_pp3_bp4_inputs, row, GERMLINE_BP4, frozenset())
    assert "InterVar recorded **BP4**" in text
    assert "takes PP3 on two of them, BP4 only on all three" in text
    assert "cutoffs above are InterVar's own" in text


def test_the_germline_section_says_when_bp4s_points_came_from_absent_columns(monkeypatch):
    """63,222 of 97,223 real BP4 firings. Silence here reads as a gap in the reasoning."""
    row = pd.Series({"ExonicFunc.refGene": "synonymous_SNV"})
    text = _drawn(monkeypatch, render_pp3_bp4_inputs, row, GERMLINE_BP4, frozenset())
    assert "3 of BP4's three points came from a column this file does not carry" in text
    assert "part of why BP4 fired rather than a gap in the reasoning" in text


def test_the_germline_section_says_when_pp3s_point_came_from_an_absent_metasvm(monkeypatch):
    row = pd.Series({"GERP++_RS": "4.7", "dbscSNV_ADA_SCORE": "0.9", "dbscSNV_RF_SCORE": "0.9"})
    text = _drawn(monkeypatch, render_pp3_bp4_inputs, row, GERMLINE_PP3, frozenset())
    assert "One of PP3's three points came from this file having no readable MetaSVM score" in text


def test_the_germline_section_draws_on_a_variant_where_neither_criterion_fired(monkeypatch):
    """The dev's decision: the numbers explain a criterion that did not fire too."""
    row = pd.Series({"MetaSVM_score": "-1.2", "GERP++_RS": "1.0",
                     "dbscSNV_ADA_SCORE": "0.1", "dbscSNV_RF_SCORE": "0.1"})
    text = _drawn(monkeypatch, render_pp3_bp4_inputs, row, GERMLINE_NEITHER, frozenset())
    assert "PP3 and BP4" in text
    assert "neither PP3 nor BP4" in text
    assert "-1.2" in text


def test_the_germline_section_draws_nothing_when_the_evidence_string_does_not_parse(monkeypatch):
    row = pd.Series({"MetaSVM_score": "-1.2"})
    assert _drawn(monkeypatch, render_pp3_bp4_inputs, row, ".", frozenset()) == ""


def test_neither_table_colours_a_row_by_which_side_it_fell_on(monkeypatch):
    """#188's lesson and #191's risk: five coloured rows would read as five votes.

    Asserted on the HTML, because the only way to tell is to look inside the one ``st.markdown``.
    The muted grey of an unreadable row and the amber of a contaminated one are allowed; the
    evidence palettes are not.
    """
    from components.cbp_evidence import _SCORE_COLORS

    row = pd.Series(
        {
            "SIFT_score": "0.01",
            "Polyphen2_HDIV_pred": "D",
            "FATHMM_pred": "T",
            "MutationAssessor_pred": "H",
            "GERP++_RS": "4.7",
        }
    )
    text = _drawn(monkeypatch, render_cbp10_inputs, row, REAL_SOMATIC, frozenset())
    for color in _SCORE_COLORS.values():
        assert color not in text, f"an evidence-strength colour reached this table: {color}"
    assert "background:" not in text


# ---------------------------------------------------------------------------
# The panel: arm gating
# ---------------------------------------------------------------------------


def _panel(monkeypatch, row):
    """Every string the whole detail panel draws for one variant.

    All three sections' ``st`` are replaced, and the third is why this helper exists rather than
    the one in ``test_cbp_evidence``: a module drawing through the real Streamlit contributes
    nothing to the transcript, which is indistinguishable from its section not being drawn.
    """
    from unittest.mock import MagicMock

    from components import acmg_evidence, cbp_evidence, predictor_context, variant_detail

    drawn = []

    class _Ctx:
        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

    fake = MagicMock()
    fake.columns.side_effect = lambda n, *a, **k: [
        _Ctx() for _ in range(n if isinstance(n, int) else len(n))
    ]
    fake.expander.side_effect = lambda *a, **k: _Ctx()
    fake.markdown.side_effect = lambda text, *a, **k: drawn.append(str(text))
    fake.caption.side_effect = lambda text, *a, **k: drawn.append(str(text))
    fake.info.side_effect = lambda text, *a, **k: drawn.append(str(text))
    fake.metric.side_effect = lambda label, value, *a, **k: drawn.append(f"{label}: {value}")
    for module in (variant_detail, cbp_evidence, acmg_evidence, predictor_context):
        monkeypatch.setattr(module, "st", fake)

    variant_detail.render_variant_detail_panel(row)
    return "\n".join(drawn)


def _fixture(name):
    return pd.read_csv(os.path.join(CANCERVAR_FIXTURES, name), sep="\t", low_memory=False)


def test_a_somatic_row_gets_the_cbp10_table_and_not_the_germline_one(monkeypatch):
    """Driven off a real fixture MAF, so the padded evidence-column spelling is the gate."""
    frame = _fixture("somatic_cancervar_evidence.maf")
    text = _panel(monkeypatch, frame.iloc[5])
    assert "The predictors behind CBP10 (CancerVar)" in text
    assert "The predictors behind PP3 and BP4 (InterVar)" not in text


def test_a_germline_row_gets_the_pp3_bp4_table_and_not_the_somatic_one(monkeypatch):
    germline = pd.Series(
        {
            "Hugo_Symbol": "BRCA2",
            "Chromosome": "chr13",
            "Start_Position": 32339132,
            "InterVar": "Likely benign",
            "InterVar: InterVar and Evidence": GERMLINE_BP4,
            "MetaSVM_score": "-1.1",
            "GERP++_RS": "1.0",
            "dbscSNV_ADA_SCORE": "0.02",
            "dbscSNV_RF_SCORE": "0.03",
        }
    )
    text = _panel(monkeypatch, germline)
    assert "The predictors behind PP3 and BP4 (InterVar)" in text
    assert "The predictors behind CBP10 (CancerVar)" not in text


def test_a_row_with_neither_evidence_column_gets_neither_table(monkeypatch):
    """Both tables are gated on the file's own evidence column, per issue #187's measurement."""
    row = pd.Series(
        {
            "Hugo_Symbol": "EGFR",
            "Chromosome": "chr7",
            "Start_Position": 55191822,
            "SIFT_score": "0.01",
            "MetaSVM_score": "0.5",
        }
    )
    text = _panel(monkeypatch, row)
    assert "The predictors behind" not in text
