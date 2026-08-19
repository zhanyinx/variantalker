"""The app's one answer to *this row does not say* — and the one place that must not ask it.

Issue #131 found four sets of sentinel strings in ``components/``, differing from each other
by members that turned out to be accidents rather than intent. ``config/missing_values.py``
replaces them with two predicates, and these are the claims that keep the replacement honest:

* the two predicates differ on **exactly** the annotator's sentinels, asserted as a derived
  relation rather than as two hand-written lists, so a member added to one is visible here;
* the display predicate catches every spelling the app's own loader can produce, read back
  *through that loader* rather than off a memorised list — the shape ``filters/gene_lists``
  uses for the same reason;
* ``-`` is a real value in every column, which extends #127's decision for the two allele
  columns to the whole app on the strength of the corpus measurement;
* and the note key asks the *other* predicate, which is the difference that would silently
  withdraw notes from a third of real files if someone folded the two together.

Filed as a unit module: nothing here has a pipeline counterpart. What the *filter* does about
a missing annotation is ``test_absent_columns.py``'s, and the two must not be conflated —
``NEUTRAL_FILLS`` answers a different question per column and its values are never rendered.
"""

from __future__ import annotations

import io
import os
import sys

import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from config.missing_values import (  # noqa: E402
    _BLANK,
    _WROTE_NOTHING,
    is_blank,
    says_nothing,
    says_nothing_over,
    shown,
)

#: What an annotator writes when it could not assign a value, with the tool that writes it.
#: Every one is measured in the corpus of 303 byte-distinct real MAFs (issue #131); the
#: module docstring carries the counts.
ANNOTATOR_SENTINELS = [
    (".", "ANNOVAR / dbNSFP: no annotation, in 422 columns"),
    ("__UNKNOWN__", "Funcotator: value not determined, in 87 columns"),
    ("UNKNOWN", "ANNOVAR: no transcript consequence"),
    ("unknown", "ANNOVAR: exonic function not assigned"),
    ("Unknown", "ANNOVAR: no gene assigned, in 136 of 303 files"),
    ("NONE", "ANNOVAR: Ref.Gene unassigned"),
]

#: Values that mean something, in the columns that spell them that way. ``-`` is the one
#: that matters: it is the absent side of an indel in six allele columns and the minus
#: strand in ``Transcript_Strand``, and no column in the corpus uses it for *no value*.
REAL_VALUES = ["-", "0", "0.0", "A", "TP53", "chr7", "p.R273H", "Benign", "false", "NaN_x"]


class TestTheTwoPredicatesDifferOnlyWhereTheyMeanTo:
    """The relation between them, not two transcribed lists.

    A test that spelled both sets out would pass whatever either one became, since the same
    hand writes both; asserting the *difference* means a member added to one shows up here.
    """

    def test_blank_implies_says_nothing(self):
        """Every empty cell is also a cell the interface should not draw."""
        for value in _BLANK:
            assert is_blank(value), value
            assert says_nothing(value), (
                f"{value!r} is blank but the display predicate would render it"
            )

    def test_the_blank_spellings_are_what_pythons_own_missing_objects_render_as(self):
        """Derived from the objects, so the set cannot quietly stop describing them.

        This is the claim that earns ``_BLANK`` its *exact* match: each member is one fixed
        rendering, not a spelling anyone chose. A transcribed list would go on passing after
        pandas changed what ``pd.NA`` prints as.
        """
        for missing in (float("nan"), None, pd.NA, pd.NaT):
            assert str(missing) in _BLANK or is_blank(missing), (
                f"str({missing!r}) == {str(missing)!r} is not read as an empty cell"
            )

    def test_the_gap_between_them_is_the_annotator_sentinels(self):
        """What ``says_nothing`` catches and ``is_blank`` does not, stated as the set."""
        assert _WROTE_NOTHING - {value.lower() for value in _BLANK} == {
            ".", "unknown", "__unknown__",
        }

    @pytest.mark.parametrize("value,who", ANNOTATOR_SENTINELS, ids=lambda v: str(v)[:12])
    def test_a_sentinel_says_nothing_but_is_not_blank(self, value, who):
        """A tool that could not assign a value has still said something about the row.

        Which is why these two answers cannot be one: the cell is not empty, so an identity
        built from it is still an identity, but there is nothing here worth a clinician's
        screen. Both halves are asserted, because a change that made these *blank* would pass
        a test that only checked the display half.
        """
        assert says_nothing(value), f"{value!r} would reach a clinician's screen ({who})"
        assert not is_blank(value), f"{value!r} is not an empty cell ({who})"


class TestWhatTheLoaderCanActuallyProduce:
    """Read back through the app's own reader, not off a remembered list.

    ``filters/gene_lists`` pins its NA-token list the same way and for the same reason: the
    set that matters is whatever *this* pandas calls missing, and it is not ours to memorise.
    The premise underneath is that the app has no nullable dtypes anywhere, so a missing cell
    always renders ``"nan"`` — if that ever stops being true, this is what says so.
    """

    #: Every token pandas treats as missing by default, plus the empty field.
    NA_SPELLINGS = [
        "", "nan", "NaN", "NA", "N/A", "n/a", "None", "null", "NULL", "<NA>", "-NaN", "#N/A",
    ]

    @pytest.mark.parametrize("spelling", NA_SPELLINGS)
    def test_a_missing_cell_says_nothing_however_the_file_spelled_it(self, spelling):
        """What the file wrote is not what the check sees — pandas got there first."""
        frame = pd.read_csv(
            io.StringIO(f"Hugo_Symbol\tCol\nTP53\t{spelling}\nBRCA1\treal\n"),
            sep="\t", low_memory=False,
        )
        value = frame["Col"].iloc[0]

        assert pd.isna(value), (
            f"{spelling!r} is no longer a missing value to this pandas; the module docstring's"
            " account of the loader is what needs revisiting, not this assertion"
        )
        assert is_blank(value)
        assert says_nothing(value)

    @pytest.mark.parametrize("value", REAL_VALUES)
    def test_a_real_value_survives_the_loader_and_is_drawn(self, value):
        """The failure direction issue #14 names: a value read as missing leaves no trace."""
        frame = pd.read_csv(
            io.StringIO(f"Hugo_Symbol\tCol\nTP53\t{value}\nBRCA1\tother\n"),
            sep="\t", low_memory=False,
        )
        read_back = frame["Col"].iloc[0]

        assert not says_nothing(read_back), f"{value!r} would disappear from the interface"
        assert not is_blank(read_back)

    def test_the_dash_is_a_real_value_in_the_columns_that_spell_it(self):
        """#127 decided this for two allele columns; #131 measured it for every column.

        204,251 values across the corpus, and every column the app reads means something by
        them: the absent side of an indel in the six allele columns, the minus strand in
        ``Transcript_Strand``. So neither predicate carries it.
        """
        for column in (
            "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
            "Match_Norm_Seq_Allele1", "Match_Norm_Seq_Allele2", "Transcript_Strand",
        ):
            frame = pd.read_csv(
                io.StringIO(f"Hugo_Symbol\t{column}\nTP53\t-\n"), sep="\t", low_memory=False,
            )
            value = frame[column].iloc[0]

            assert not says_nothing(value), f"an indel or a strand would vanish from {column}"


class TestTheVectorisedTwinCannotDrift:
    """``says_nothing_over`` exists for speed and must not become a second opinion."""

    def test_it_agrees_with_the_row_wise_predicate_value_for_value(self):
        values = (
            list(_BLANK)
            + [value for value, _ in ANNOTATOR_SENTINELS]
            + REAL_VALUES
            + [None, float("nan"), 0, 1.5]
        )
        series = pd.Series(values, dtype=object)

        expected = [says_nothing(v) for v in values]

        assert list(says_nothing_over(series)) == expected

    def test_it_returns_a_mask_aligned_with_its_input(self):
        """Used inside ``Series.where``, so a reset index would misplace every fallback."""
        series = pd.Series(["TP53", "__UNKNOWN__", "BRCA1"], index=[7, 42, 99])

        mask = says_nothing_over(series)

        assert list(mask.index) == [7, 42, 99]
        assert list(mask) == [False, True, False]


class TestShownDrawsADashRatherThanAValueNobodyWrote:
    """For the fields drawn unconditionally, where the choice is *what* and not *whether*."""

    def test_an_empty_cell_in_a_present_column_renders_a_dash(self):
        """The hole the old ``row.get(column, "—")`` default could not defend against.

        A default only fires when the *column* is absent. A column that is present and empty
        went straight through it, so the panel drew the string ``nan`` at a clinician.
        """
        row = pd.Series({"Hugo_Symbol": float("nan")})

        assert shown(row.get("Hugo_Symbol")) == "—"
        assert shown(row.get("Hugo_Symbol", "—")) == "—"

    def test_a_real_value_is_rendered_stripped(self):
        assert shown("  TP53  ") == "TP53"
        assert shown("-") == "-"

    def test_the_dash_is_the_callers_to_choose(self):
        assert shown(float("nan"), dash="not reported") == "not reported"
