"""A predictor column holding the chromosome never reaches the ClinGen scale (issue #194).

Three claims, in three groups.

**The detector** (:mod:`config.contaminated_columns`) answers a whole-column question and gets
both spellings right. It has to: ``normalise_chromosome_spelling`` has already prefixed the
chromosome column by the time it runs, so the score cell holds ``22`` while the chromosome holds
``chr22``, and a column pandas typed as float renders the same cell as ``22.0``. It must also
*not* fire on the measured false-positive case — a legitimate score of exactly ``1.0`` on a
``chr1`` variant, which is 1,622 real ``SIFT_score`` cells and 581 real ``Polyphen2_HDIV_score``
cells across the corpus.

**The panel** draws the value with a warning, gives it no strength, and leaves it out of the
tally. Each of those three is asserted separately, because the failure that matters is a
chromosome being *scored*, and a guard that only checked the caption would pass a table that
warned and voted at the same time.

**The seam** between them: the load page takes the verdict after the chromosome spelling is
settled, not before. Asserted by ordering the two calls, because taking it first is not an
error — it is a wrong answer on the 9 bare-spelled files in the corpus, and nothing raises.

Why the panel cannot answer this itself, which is what shapes the seam: contamination is a fact
about a *column*, and the panel holds one row. Issue #194 measured the per-cell alternative and
rejected it — see the detector's module docstring for the numbers.
"""

from __future__ import annotations

import os
import sys
from unittest.mock import MagicMock

import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from config.contaminated_columns import (  # noqa: E402
    NOTHING_CONTAMINATED,
    PREDICTOR_COLUMNS,
    SESSION_KEY,
    contaminated_columns,
    holds_the_chromosome,
)
from tests.fakes import FakeSessionState  # noqa: E402


# ---------------------------------------------------------------------------
# The detector
# ---------------------------------------------------------------------------


def _frame(chromosomes, **columns):
    """A MAF-shaped frame: a prefixed ``Chromosome`` column and whatever else is asked for."""
    return pd.DataFrame({"Chromosome": chromosomes, **columns})


def test_a_column_of_bare_chromosomes_is_caught():
    """The shape on disk: `CADD_phred` holding `22` beside a `Chromosome` of `chr22`."""
    frame = _frame(
        ["chr22", "chr1", "chr7", "chrX"],
        CADD_phred=["22", "1", "7", "X"],
    )
    assert contaminated_columns(frame) == {"CADD_phred"}


def test_a_float_typed_column_is_caught_too():
    """`22.0` and `chr22` are the same cell; only a numeric comparison sees it.

    This is not hypothetical tidiness. ``read_maf`` hands the column to ``read_csv``, which
    types a column of digits as float64 — so the string comparison alone reads `22.0` against
    `22`, finds them different, and calls a chromosome a CADD score.
    """
    frame = _frame(
        ["chr22", "chr1", "chr7"],
        CADD_phred=[22.0, 1.0, 7.0],
    )
    assert contaminated_columns(frame) == {"CADD_phred"}


def test_the_sex_chromosomes_are_caught_by_the_string_comparison():
    """`X` and `Y` parse as neither score nor number, so only the string comparison sees them.

    Not a corner: 1,066 of the contaminated `CADD_phred` cells in the corpus hold a letter
    rather than a digit, and 546 `phastCons20way_mammalian` cells hold the literal `X`. A
    numeric rule alone reads a column of an X-heavy region as clean, and every one of its cells
    then reaches `reference_scales.parse_score`, whose `float()` refuses them — so they read as
    *missing* and the row silently vanishes instead of being warned about. The failure is
    quieter than the numeric one, not smaller.
    """
    frame = _frame(
        ["chrX", "chrY", "chrX", "chrX"],
        CADD_phred=["X", "Y", "X", "X"],
    )
    assert contaminated_columns(frame) == {"CADD_phred"}


def test_a_real_score_of_one_on_chr1_does_not_fire():
    """The measured false-positive case, at a rate well above the corpus's worst file.

    ``SIFT_score`` is legitimately ``1.0`` on 1,622 real cells and ``Polyphen2_HDIV_score`` on
    581, and where those land on a ``chr1`` variant a per-cell rule flags a correct number.
    Three of these twenty rows collide — 15%, which is the highest fraction any clean file in
    the corpus reaches (one Polyphen2 file, 0.1509). The column-level rule must be quiet here.
    """
    chromosomes = ["chr1"] * 3 + [f"chr{i}" for i in range(2, 19)]
    scores = [1.0, 1.0, 1.0] + [0.03 * i for i in range(2, 19)]
    frame = _frame(chromosomes, SIFT_score=scores)

    assert contaminated_columns(frame) == set()


def test_the_threshold_is_a_majority_not_a_trace():
    """Just over half is contamination; just under is coincidence."""
    chromosomes = [f"chr{i}" for i in range(1, 11)]
    # 6 of 10 cells equal their chromosome.
    contaminated = [str(i) for i in range(1, 7)] + [0.5] * 4
    # 5 of 10 — a tie is not a majority.
    tied = [str(i) for i in range(1, 6)] + [0.5] * 5

    assert holds_the_chromosome(pd.Series(contaminated), pd.Series(chromosomes)) is True
    assert holds_the_chromosome(pd.Series(tied), pd.Series(chromosomes)) is False


def test_missing_cells_are_not_counted_either_way():
    """The fraction is of the cells that *say* something — `.` is a different question.

    A column of two real chromosomes and eighteen dots is contaminated: everything it holds is
    a chromosome. Counting the dots in the denominator would call it clean at 10%, which is how
    a mostly-empty contaminated column would slip through.
    """
    frame = _frame(
        ["chr3", "chr9"] + [f"chr{i}" for i in range(1, 19)],
        CADD_phred=["3", "9"] + ["."] * 18,
    )
    assert contaminated_columns(frame) == {"CADD_phred"}


def test_a_column_that_says_nothing_at_all_is_not_contaminated():
    """An empty column is a missing column's shape, and not this module's business."""
    frame = _frame(["chr1", "chr2"], CADD_phred=[".", "."])
    assert contaminated_columns(frame) == set()


def test_without_a_chromosome_column_there_is_no_verdict():
    """Nothing to compare against, and `validate_required_columns` refuses the file anyway."""
    frame = pd.DataFrame({"CADD_phred": ["22", "1"]})
    assert contaminated_columns(frame) == frozenset()


def test_the_check_is_not_gated_on_the_arm():
    """The same column is caught whichever classifier's columns sit beside it (issue #194 Q3).

    ``Polyphen2_HDIV_score`` is positional on somatic and name-matched on germline, so an
    arm-gated check would skip it on a germline file. This rule measures instead of predicting,
    so provenance moving upstream cannot make it go quiet — the reason is in the detector's
    module docstring, and this is the claim itself.
    """
    cells = {"Chromosome": ["chr22", "chr4"], "Polyphen2_HDIV_score": ["22", "4"]}

    somatic = pd.DataFrame({**cells, "CancerVar": ["Tier_III_Uncertain"] * 2})
    germline = pd.DataFrame({**cells, "InterVar": ["Uncertain significance"] * 2})

    assert contaminated_columns(somatic) == {"Polyphen2_HDIV_score"}
    assert contaminated_columns(germline) == {"Polyphen2_HDIV_score"}


def test_every_column_the_clingen_table_scores_is_asked_about():
    """The five rows on screen are all in scope, so none can be the one that is not checked.

    ``CLINGEN_SVI_TABLE`` is the table the panel draws; :data:`PREDICTOR_COLUMNS` is what
    the detector measures. If a row were added to the first and not the second, the new row
    would score a chromosome in silence — which is precisely the state ``CADD_phred`` was in
    before this ticket.

    The table moved to :mod:`components.reference_scales` in issue #191 and was renamed from
    ``_PREDICTOR_TABLE_CONFIG``, because since issue #190 the panel draws three predictor tables
    and only this one is ClinGen's.
    """
    from components.reference_scales import CLINGEN_SVI_TABLE

    drawn = {col for _, col, _, _, _ in CLINGEN_SVI_TABLE if col is not None}
    assert drawn <= set(PREDICTOR_COLUMNS), sorted(drawn - set(PREDICTOR_COLUMNS))


def test_a_column_the_file_spells_its_own_way_is_still_measured():
    """``GERP.._RS`` on 2 of 167 real MAFs (issue #190).

    The detector looks the column up under the file's own spelling and reports it under the
    canonical one, so every reader keeps asking the question it already asks. Without this, a
    mangled column is measured by nothing and the score-context table would draw a chromosome as a
    conservation score with no warning at all.
    """
    import pandas as pd

    maf = pd.DataFrame(
        {
            "Chromosome": ["chr1", "chr7", "chr22", "chrX"],
            "GERP.._RS": ["1", "7", "22", "X"],
        }
    )
    assert contaminated_columns(maf) == {"GERP++_RS"}


def test_a_mangled_column_that_is_clean_is_not_flagged():
    """The resolver widens what is measured, not what is warned about."""
    import pandas as pd

    maf = pd.DataFrame(
        {
            "Chromosome": ["chr1", "chr7", "chr22", "chrX"],
            "GERP.._RS": ["4.71", "-2.3", "0.5", "5.9"],
        }
    )
    assert contaminated_columns(maf) == frozenset()


def test_every_column_the_score_context_tables_draw_is_asked_about():
    """The same rule for issue #190's two tables, which is now the larger of the two surfaces.

    ``components/predictor_context.py`` draws nine columns — CancerVar's five and InterVar's four,
    overlapping on ``GERP++_RS`` — and states a call for each. A column it drew that the detector
    did not measure would print ``chr7`` as PolyPhen-2's prediction, which is #194's failure in a
    new place.

    Asked of the rows the tables actually build rather than of a config literal, because the two
    tables are assembled by function calls rather than driven off one table; a literal would be a
    second list to keep in step.
    """
    import pandas as pd

    from components.predictor_context import _cbp10_readings, _pp3_bp4_readings

    empty = pd.Series(dtype=object)
    drawn = {r.column for r in _cbp10_readings(empty, frozenset())}
    drawn |= {r.column for r in _pp3_bp4_readings(empty, frozenset())}
    assert len(drawn) == 8, sorted(drawn)  # nine rows, GERP++_RS drawn on both arms
    assert drawn <= set(PREDICTOR_COLUMNS), sorted(drawn - set(PREDICTOR_COLUMNS))


# ---------------------------------------------------------------------------
# The panel
# ---------------------------------------------------------------------------


def _predictor_table(monkeypatch, row, contaminated=()):
    """Every string the predictor section draws, for one variant.

    Returns the HTML table and the captions together, in the order drawn: the value and the
    strength cells are inside the table's markdown while the tally and the warning are captions,
    and each of the three claims below reads a different one of those.

    Drives :func:`components.reference_scales.render_reference_scale` since issue #191 moved the
    table there. The contaminated set is now an *argument* rather than something the section
    reads out of the session — the panel still reads it from there and hands it down — so it is
    passed both ways here: as the argument the function takes, and in the fake session, so a
    section that started reading the session directly would not quietly pass on the argument
    alone.
    """
    from components import reference_scales

    drawn = []

    fake = MagicMock()
    fake.session_state = FakeSessionState(**{SESSION_KEY: frozenset(contaminated)})
    fake.markdown.side_effect = lambda text, *a, **k: drawn.append(str(text))
    fake.caption.side_effect = lambda text, *a, **k: drawn.append(str(text))
    monkeypatch.setattr(reference_scales, "st", fake)

    reference_scales.render_reference_scale(row, "InterVar", frozenset(contaminated))
    return drawn


#: One variant with a clean REVEL and a CADD of 22.0 — the value a `chr22` variant's
#: contaminated `CADD_phred` holds, and one that classifies as BP4 Supporting (<=22.7) if the
#: panel is allowed to score it.
def _row(**overrides):
    cells = {
        "Hugo_Symbol": "TP53",
        "Chromosome": "chr22",
        "REVEL_score": 0.412,
        "CADD_phred": 22.0,
    }
    cells.update(overrides)
    return pd.Series(cells)


def test_a_contaminated_score_is_shown_but_not_scored(monkeypatch):
    """The value reaches the screen; the ClinGen verdict does not."""
    drawn = _predictor_table(monkeypatch, _row(), contaminated=["CADD_phred"])
    table = "\n".join(t for t in drawn if "<table" in t)

    assert "⚠ 22" in table, "the value the file holds is not on screen"
    # The row is present, and its two strength cells hold the em dash rather than a strength.
    assert "CADD" in table
    assert "Supporting" not in table.split("CADD", 1)[1], (
        "a chromosome was given a place on the ClinGen scale"
    )


def test_a_contaminated_cell_on_chrx_is_warned_about_too(monkeypatch):
    """The letter `X` in a score column gets the same row and the same warning as `22`.

    The order the panel does two things in, asserted by its effect. `reference_scales.parse_score`
    refuses `X` — rightly, it is not a score — so a panel that parsed before it checked the
    verdict would drop this row silently. Every row of a contaminated file holds a chromosome,
    so warning on the autosomes and going quiet on X and Y would put the panel at its least
    trustworthy where the reader had least reason to suspect it.
    """
    row = _row(Chromosome="chrX", CADD_phred="X")
    drawn = _predictor_table(monkeypatch, row, contaminated=["CADD_phred"])
    table = "\n".join(t for t in drawn if "<table" in t)

    assert "⚠ X" in table, table
    assert [t for t in drawn if t.startswith("⚠")], "no warning was drawn for a chrX cell"


def test_a_contaminated_column_absent_from_this_row_draws_nothing(monkeypatch):
    """A verdict about a column this variant has no cell in is not a row to warn about."""
    row = _row()
    del row["CADD_phred"]
    drawn = _predictor_table(monkeypatch, row, contaminated=["CADD_phred"])

    assert not [t for t in drawn if t.startswith("⚠")], drawn
    assert "CADD" not in "\n".join(t for t in drawn if "<table" in t)


def test_a_contaminated_score_is_out_of_the_tally(monkeypatch):
    """The denominator counts what scored, not what was drawn.

    With REVEL alone scoring, the caption reads `1/1`. Counting the warned CADD row would make
    it `1/2` — a larger evidence base than actually spoke.
    """
    drawn = _predictor_table(monkeypatch, _row(), contaminated=["CADD_phred"])
    tally = [t for t in drawn if "PP3 (pathogenic)" in t]

    assert tally, "the tally caption stopped being drawn"
    assert "/1 predictors" in tally[0], tally[0]


def test_the_warning_names_the_row_and_the_column(monkeypatch):
    """Both, so the reader can find the row on screen and check it against their own file."""
    drawn = _predictor_table(monkeypatch, _row(), contaminated=["CADD_phred"])
    warnings = [t for t in drawn if t.startswith("⚠")]

    assert len(warnings) == 1, drawn
    assert "CADD" in warnings[0]
    assert "CADD_phred" in warnings[0]
    assert "chromosome" in warnings[0]


def test_a_clean_file_is_scored_exactly_as_before(monkeypatch):
    """No verdict, no change: the same CADD value is scored and counted.

    The pair with the first two tests. Without this, suppressing every CADD row unconditionally
    would pass both of them.
    """
    drawn = _predictor_table(monkeypatch, _row(), contaminated=[])
    table = "\n".join(t for t in drawn if "<table" in t)
    tally = [t for t in drawn if "PP3 (pathogenic)" in t]

    assert "⚠" not in table
    assert "Supporting" in table, "a clean CADD of 22.0 is BP4 Supporting and should say so"
    assert "/2 predictors" in tally[0], tally[0]


def test_a_panel_with_no_verdict_stashed_draws_what_the_file_says(monkeypatch):
    """A render reached without one of the load paths behaves as it did before issue #194.

    The panel cannot measure a column from one row, so when it has not been told it draws the
    file's own values — the same branch a clean file takes.

    Since issue #191 the section takes the verdict as an argument and the *panel* is what reads
    it out of the session, so this drives both halves: it asks ``variant_detail`` for the verdict
    with no key stashed, and hands whatever it gets to the section — which is exactly the wiring
    ``render_variant_detail_panel`` performs.
    """
    from components import reference_scales, variant_detail

    drawn = []
    fake = MagicMock()
    fake.session_state = FakeSessionState()  # no key at all
    fake.markdown.side_effect = lambda text, *a, **k: drawn.append(str(text))
    fake.caption.side_effect = lambda text, *a, **k: drawn.append(str(text))
    monkeypatch.setattr(variant_detail, "st", fake)
    monkeypatch.setattr(reference_scales, "st", fake)

    verdict = variant_detail._untrustworthy_columns()
    reference_scales.render_reference_scale(_row(), "InterVar", verdict)

    assert verdict is NOTHING_CONTAMINATED
    assert not [t for t in drawn if t.startswith("⚠")]


def test_when_nothing_scores_no_tally_is_drawn(monkeypatch):
    """`0/0 predictors` is not a statement; the warning is.

    The one-row case: a file whose only drawable predictor is contaminated. The table still
    draws — the reader sees what their file holds — and the tally, which would have no base at
    all, is skipped rather than printed empty.
    """
    row = _row(REVEL_score=".")
    drawn = _predictor_table(monkeypatch, row, contaminated=["CADD_phred"])

    assert [t for t in drawn if t.startswith("⚠")], "the warning must still be drawn"
    assert not [t for t in drawn if "PP3 (pathogenic)" in t], (
        "a tally over an empty base was drawn"
    )


# ---------------------------------------------------------------------------
# The seam
# ---------------------------------------------------------------------------


def test_the_verdict_is_taken_after_the_chromosome_spelling_is_settled():
    """Order, not presence — and the two are not interchangeable.

    The verdict is a comparison against ``Chromosome``. Taken before
    ``normalise_chromosome_spelling``, it compares a bare-spelled file's ``22`` against its own
    bare ``22`` and is *still* right by luck — but a prefixed file's ``chr22`` against ``22``
    only agrees because :func:`holds_the_chromosome` strips the prefix. The order is what makes
    the answer a fact rather than a coincidence of which spelling arrived, and the corpus holds
    9 bare-spelled files to be wrong about.

    Read from the source rather than by running the page, because the load tail needs a
    Streamlit session and this claim is about a sequence.
    """
    import ast
    import inspect
    from pathlib import Path

    source = Path(inspect.getsourcefile(_load_page())).read_text(encoding="utf-8")
    tree = ast.parse(source)

    calls = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        name = getattr(node.func, "id", None) or getattr(node.func, "attr", None)
        if name in ("normalise_chromosome_spelling", "contaminated_columns"):
            calls.append((node.lineno, name))

    ordered = [name for _, name in sorted(calls)]
    assert ordered.count("normalise_chromosome_spelling") == 1, ordered
    assert ordered.count("contaminated_columns") == 1, ordered
    assert ordered.index("normalise_chromosome_spelling") < ordered.index(
        "contaminated_columns"
    ), ordered


def test_the_load_page_stashes_the_verdict_under_the_shared_key():
    """One key, spelled once — the writer and the reader import the same name.

    Without this, the page could stash the verdict under a name the panel does not read, and
    every panel test above would still pass because they set the key themselves.
    """
    import ast
    import inspect
    from pathlib import Path

    page = Path(inspect.getsourcefile(_load_page())).read_text(encoding="utf-8")
    panel = Path(
        inspect.getsourcefile(__import__("components.variant_detail", fromlist=["x"]))
    ).read_text(encoding="utf-8")

    for source, who in ((page, "the load page"), (panel, "the panel")):
        tree = ast.parse(source)
        imported = {
            alias.asname or alias.name
            for node in ast.walk(tree)
            if isinstance(node, ast.ImportFrom)
            and node.module == "config.contaminated_columns"
            for alias in node.names
        }
        assert imported, f"{who} does not import the shared key"
        assert "SESSION_KEY" in {
            alias.name
            for node in ast.walk(tree)
            if isinstance(node, ast.ImportFrom)
            and node.module == "config.contaminated_columns"
            for alias in node.names
        }, f"{who} spells the session key itself instead of importing it"

    # And the literal is not also spelled inline anywhere it could drift.
    assert SESSION_KEY not in page.replace(f'SESSION_KEY = "{SESSION_KEY}"', "")
    assert SESSION_KEY not in panel


def _load_page():
    """``page_modules.data_loading``, or a skip if Streamlit is not installed."""
    try:
        import page_modules.data_loading as page
    except ImportError as exc:  # pragma: no cover - streamlit-free CI job
        pytest.skip(f"the load page needs Streamlit: {exc}")
    return page
