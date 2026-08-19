"""How many times the missing-column sentence is said, and in which states (issue #93).

The resolver reports what a MAF does not carry, and that report used to be drawn where it
was made: ``st.warning`` inside :func:`components.variant_table._visible_columns`, on every
call. The results section called it four times — two grids, two downloads — so one sentence
rendered four times verbatim on one screen. Issue #92 took the downloads off the resolver
and made them follow the grid, which brought it to two.

Two was the harder half, because the accidental-looking explanation was gone: two is one
per tab, which is that function's documented design working as intended. What settled it
was a state the ticket had not measured. **Ticking *Show all columns* rendered the sentence
zero times** — that branch of ``create_data_table`` skipped ``shown_columns`` and so skipped
the resolver, meaning the app went silent about a MAF missing three pipeline columns exactly
when the user asked to see everything. So the count was never two; it was nought, one or two
depending on a checkbox that has nothing to do with which columns the file carries.

What is asserted here
---------------------
The sentence is said **once per render, at page level**, in every state where it is true:

* once with the grid in its default state, and once with *Show all columns* ticked on
  either or both tabs — the same number, because the fact does not depend on the view;
* once for a MAF missing only *output* columns, where nothing else in the app speaks. The
  filter's own warning covers filter *inputs*, so a MAF missing ``am_class`` and ``cosmic``
  gets nothing from it, and this sentence is the whole of what the user is told;
* never for a complete MAF, and never on the Load section, which draws no table.

And the wiring, which is the part that could fail silently: the section must draw none
itself, and the page must drain what the section collected. A collector with no drain is
the silence ``_visible_columns``'s docstring has forbidden since it was written — it would
pass any test that only counted what the page shows, since nought is also what a complete
MAF shows.

Every assertion below was mutation-checked; the mutations are named in the test that
catches them, because a green test over a wired-up page proves nothing on its own (the
lesson of issues #67, #83 and #90).
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

STREAMLIT_APP = Path(__file__).resolve().parent.parent
if str(STREAMLIT_APP) not in sys.path:
    sys.path.insert(0, str(STREAMLIT_APP))

REFERENCE_MAF = STREAMLIT_APP / "tests" / "fixtures" / "parity" / "somatic_reference.maf"

#: A filter input and two output-only columns. The mix is the point: ``ESCAT`` is read by
#: the filter, so its absence is also reported by ``filters/absent_columns.py`` at filter
#: time, while ``am_class`` and ``cosmic`` are emitted by the pipeline and read by nothing,
#: so the resolver is the only thing that can speak for them.
A_FILTER_INPUT = "ESCAT"
OUTPUT_ONLY = ["am_class", "cosmic"]

#: The clause that identifies this sentence wherever it is drawn. Deliberately not the
#: whole literal: the count and the column names vary per file, and a test that pinned the
#: full string would have to be rewritten to change a word rather than to change a claim.
#: The wording is owned by ``config/columns.py``; what this module owns is how often it is
#: said and from where.
SENTENCE = "does not carry"

_PAGE = """
import os, sys
sys.path.insert(0, {app!r})

import streamlit as st
from config.pipeline_params import pipeline_params
from page_modules import data_loading
from utils import read_maf

st.session_state.filter_params = pipeline_params("somatic")
st.session_state.maf_source_name = os.path.basename({maf!r})
st.session_state.maf_data = read_maf({maf!r})
data_loading.apply_filters_to_data(show_messages=False)
st.session_state["data_page_section"] = {section!r}

{seed}

data_loading.show_data_loading_page()
"""

#: The results section alone, with no page around it to drain what it collects. Used to
#: assert the section draws the sentence nowhere itself — the half of "page level" that
#: counting the page's own warnings cannot see.
_SECTION_ONLY = """
import os, sys
sys.path.insert(0, {app!r})

import streamlit as st
from config.pipeline_params import pipeline_params
from page_modules import data_loading
from utils import read_maf

st.session_state.filter_params = pipeline_params("somatic")
st.session_state.maf_source_name = os.path.basename({maf!r})
st.session_state.maf_data = read_maf({maf!r})
data_loading.apply_filters_to_data(show_messages=False)

{seed}

data_loading.show_results_section()
"""

#: Ticking *Show all columns* before the grid is drawn. Seeded as widget state rather than
#: driven through ``AppTest``, for the reason issue #92's tests record: an interaction plus
#: a re-run draws the results section twice in one run and trips the duplicate-widget-key
#: error the page's own ``except`` swallows, so the run ends having reported nothing.
SHOW_ALL_PASSED = 'st.session_state["show_all_passed_variants"] = True\n'
SHOW_ALL_FAILED = 'st.session_state["show_all_failed_variants"] = True\n'


def _maf_without(tmp_path: Path, columns) -> Path:
    """The somatic reference with ``columns`` dropped, as a MAF on disk."""
    from utils import read_maf

    maf = read_maf(str(REFERENCE_MAF))
    present = [col for col in columns if col in maf.columns]
    assert present == list(columns), (
        f"the reference MAF does not carry {sorted(set(columns) - set(present))}, so "
        "dropping them proves nothing"
    )
    path = tmp_path / "short.maf"
    maf.drop(columns=present).to_csv(path, sep="\t", index=False)
    return path


def _render(script: str, maf: Path, seed: str = "", section: str = "results"):
    pytest.importorskip("streamlit", reason="streamlit not installed")
    from streamlit.testing.v1 import AppTest

    app = AppTest.from_string(
        script.format(app=str(STREAMLIT_APP), maf=str(maf), seed=seed, section=section)
    )
    app.run(timeout=300)
    assert not app.exception, [str(e.value) for e in app.exception]
    # The page catches its own render errors and turns them into `st.error`, so an
    # exception-free run is not by itself a rendered results table.
    assert not app.error, [element.value for element in app.error]
    return app


def _copies(app) -> list:
    """Every rendering of the sentence, however many containers deep."""
    return [w.value for w in app.warning if SENTENCE in w.value]


@pytest.fixture
def incomplete_maf(tmp_path):
    """A MAF missing one filter input and two output-only columns."""
    return _maf_without(tmp_path, [A_FILTER_INPUT] + OUTPUT_ONLY)


class TestHowManyTimesItIsSaid:
    """Once per render, whatever the grid is set to."""

    def test_it_is_said_once(self, incomplete_maf):
        """Mutation: restore the ``st.warning`` in ``_visible_columns`` → **three** copies.

        Three rather than two, which is worth stating because two is the number the ticket
        was written against: the grids would draw their two *and* the page would still draw
        the drained one, the two mechanisms adding rather than one replacing the other.

        Mutation: return the collection from ``drain_missing_column_reports`` without
        deduplicating → two copies, the grids' identical sentences both reaching the slot.
        """
        copies = _copies(_render(_PAGE, incomplete_maf))
        assert len(copies) == 1, f"expected one copy, got {len(copies)}: {copies}"
        for column in [A_FILTER_INPUT] + OUTPUT_ONLY:
            assert column in copies[0]

    @pytest.mark.parametrize(
        "seed,state",
        [
            (SHOW_ALL_PASSED, "passed only"),
            (SHOW_ALL_FAILED, "failed only"),
            (SHOW_ALL_PASSED + SHOW_ALL_FAILED, "both tabs"),
        ],
    )
    def test_show_all_columns_does_not_silence_it(self, incomplete_maf, seed, state):
        """The hole issue #93 found: measured at **zero** copies before this change.

        Mutation: move ``shown_columns`` back inside the ``else`` branch of
        ``create_data_table`` → nought copies with both tabs ticked, one with either.

        Showing every column the file *has* cannot show one it has not, so the sentence is
        exactly as true here as in the default view. It is worth stating in those terms
        because the opposite reading is available and wrong: that *show all* means nothing
        is being withheld, so there is nothing left to say.
        """
        copies = _copies(_render(_PAGE, incomplete_maf, seed=seed))
        assert len(copies) == 1, (
            f"with Show all columns ticked on {state}, expected one copy, "
            f"got {len(copies)}: {copies}"
        )

    def test_a_complete_maf_says_nothing(self, tmp_path):
        """The reference MAF carries every column the somatic arm emits."""
        assert _copies(_render(_PAGE, REFERENCE_MAF)) == []

    @pytest.mark.parametrize("section", ["filter_run", "results"])
    def test_a_leftover_from_a_previous_render_is_not_redrawn(
        self, incomplete_maf, section
    ):
        """What a previous render collected is cleared, not said again.

        The failure this guards is the sentence outliving the screen that produced it: the
        user reads it under the table, presses Load Data, and is told again about a table
        no longer on screen — or stays on Results and reads last render's copy beside this
        render's. So the leftover is planted rather than produced, and it carries text no
        resolver would ever emit, which is what lets the assertion tell a survivor from a
        fresh report on the ``results`` leg.

        Mutation: drop ``clear_missing_column_reports`` from the top of
        ``show_data_loading_page`` → the planted sentence is drawn on both legs. Mutation:
        make ``drain_missing_column_reports`` read without emptying → the same, one render
        later.

        Planted rather than driven through two renders because ``AppTest`` cannot re-run
        this page at all: reconstructing widget state calls ``SegmentedControl.value``,
        which indexes its *formatted* options with an unformatted one and raises
        ``ValueError`` on any page carrying the section switch. Nothing else in this suite
        re-runs this page, so the limit was not known before.
        """
        stale = "⚠️ This MAF does not carry 99 column(s) from a previous render"
        app = _render(
            _PAGE,
            incomplete_maf,
            seed=f'st.session_state["missing_column_reports"] = [{stale!r}]\n',
            section=section,
        )
        assert stale not in [w.value for w in app.warning], (
            "a report collected before this render reached the screen"
        )
        assert len(_copies(app)) == (1 if section == "results" else 0)


class TestWhereItIsSaid:
    """Page level: the section collects, the page speaks."""

    def test_the_results_section_draws_none_itself(self, incomplete_maf):
        """Rendered alone, the section draws the sentence nowhere.

        This is the half that counting the page's warnings cannot see. Both surviving
        copies used to be inside tabs — one of them the Failed tab, which a user may never
        open — and a page-level count of one is equally consistent with the sentence
        having simply moved into the first tab.

        Mutation: restore the ``st.warning`` in ``_visible_columns`` → two copies here.
        The page test catches that one too, at three; what this adds is *where*, which a
        count alone cannot say — it fails identically whether the page draws three or the
        tabs do.
        """
        app = _render(_SECTION_ONLY, incomplete_maf)
        assert _copies(app) == []

    def test_the_section_still_collects_what_it_does_not_draw(self, incomplete_maf):
        """Drawing nothing is only right if something was collected.

        Without this, ``test_the_results_section_draws_none_itself`` passes just as well
        against a resolver that reports nothing at all — which is precisely the silent
        export loss ``_visible_columns`` has forbidden since it was written.

        Mutation: drop the ``setdefault(...).append(...)`` in ``_visible_columns`` → the
        collection is empty here. The page tests catch that mutation as well, since a
        report never collected is never drawn; what this one adds is the *count*, which is
        how a grid that quietly stopped reporting would be told from the dedup doing its
        job. Both of them saying nothing looks identical from the page.
        """
        from components.variant_table import _MISSING_COLUMN_REPORTS

        app = _render(_SECTION_ONLY, incomplete_maf)
        collected = app.session_state[_MISSING_COLUMN_REPORTS]
        assert len(collected) == 2, (
            "both grids should report, so that dedup is what collapses them rather than "
            f"one of them staying quiet; collected {collected}"
        )
        assert len(set(collected)) == 1, collected

    def test_the_page_drains_what_it_showed(self, incomplete_maf):
        """Nothing is left behind for the next render to redraw.

        Mutation: make ``drain_missing_column_reports`` read without emptying → the
        collection survives the render, and a second render of a page whose section drew
        no grid would repeat the sentence.
        """
        from components.variant_table import _MISSING_COLUMN_REPORTS

        app = _render(_PAGE, incomplete_maf)
        assert app.session_state[_MISSING_COLUMN_REPORTS] == []


class TestWhatOnlyThisSentenceSays:
    """Why it could not simply be dropped in favour of the filter's account."""

    def test_output_only_columns_are_reported_by_nothing_else(self, tmp_path):
        """A MAF missing ``am_class`` and ``cosmic`` gets one warning, and it is this one.

        These are emitted by the pipeline and read by no filter, so
        ``filters/absent_columns.py`` has nothing to say about them — which is what makes
        this sentence load-bearing rather than a third copy of the filter's. The two
        messages were deliberately left separate (issue #93) on exactly this evidence.

        Mutation: delete the drain from ``_show_stashed_banners`` → nought warnings on a
        MAF that is genuinely missing two columns, with no other surface covering it.
        """
        maf = _maf_without(tmp_path, OUTPUT_ONLY)
        app = _render(_PAGE, maf)

        assert len(_copies(app)) == 1
        assert len(app.warning) == 1, (
            "the filter reads neither column, so this is the whole of what the app says: "
            f"{[w.value for w in app.warning]}"
        )
        for column in OUTPUT_ONLY:
            assert column in _copies(app)[0]

    def test_a_missing_filter_input_is_still_reported_by_the_filter(
        self, incomplete_maf
    ):
        """Three adjacent messages, each on its own clock, kept apart on purpose.

        The filter's answers *your verdict may be wrong*; the resolver's answers *your
        table and your export are narrower*. They overlap on ``ESCAT`` and on the clause
        about the table and the export, and that overlap was weighed and accepted rather
        than merged, because a merged message would have to be true of both a filter run
        and a render — two different clocks.

        The third is the filter recap's (issue #137), and it is a third clock rather than a
        third copy: the first two are drawn once, by the render that follows the run, and
        are gone from the screen after it — a user who scrolls back to this report an hour
        later has no record that a column was filled at all. The recap's is carried with
        the report for as long as the report exists. It is also the only one of the three
        the user has to **ask** for: it is inside the collapsed ``Filters applied to this
        report`` expander, which ``AppTest`` counts whether or not a reader would have
        opened it.

        The count is what this asserts, so it fails when a fourth voice appears — which is
        how it caught this one.
        """
        app = _render(_PAGE, incomplete_maf)
        warnings_shown = [w.value for w in app.warning]

        assert len(warnings_shown) == 3, warnings_shown
        filter_account = [
            w for w in warnings_shown if "Missing filter-input column(s)" in w
        ]
        recap_account = [w for w in warnings_shown if "to filter it" in w]
        assert len(filter_account) == 1, warnings_shown
        assert len(recap_account) == 1, warnings_shown
        assert A_FILTER_INPUT in filter_account[0]
        assert A_FILTER_INPUT in recap_account[0]
        # Neither the filter nor its recap says anything about the output-only pair; only
        # the resolver can, which is what keeps its sentence load-bearing.
        for column in OUTPUT_ONLY:
            assert column not in filter_account[0]
            assert column not in recap_account[0]


class TestTheSentenceItself:
    """The wording has to be true in every state it is now said in."""

    def test_it_puts_the_absence_on_the_file(self, incomplete_maf):
        """Not "left out of the table", which reads as MAFigate doing the leaving.

        It is said with *Show all columns* ticked now, and in that state the table is
        showing every column the file has. Nothing is being left out; the columns were
        never there.

        Mutation: restore the old literal → this fails on the clause, and six others fail
        on the count, because :data:`SENTENCE` identifies the sentence by the very phrase
        that changed. That coupling is deliberate: a reworded sentence should be looked at
        rather than silently kept. This is the one that says *which* words, and its message
        is the one worth reading when they all go red at once.
        """
        copy = _copies(_render(_PAGE, incomplete_maf, seed=SHOW_ALL_PASSED))[0]
        assert "left out of the table" not in copy
        assert "This MAF does not carry" in copy
        assert "in neither the table nor the export" in copy


# --- The *load-time* note, which is a different clock (issue #127) ---------------------
#
# Everything above is the resolver's account, made at render time about what the table and the
# export cannot carry. `validate_required_columns` makes a second one at **load** time, about
# the display columns the arm expects, and #93 measured the two as naming disjoint sets — so
# neither is a duplicate of the other and both are guarded.
#
# Nothing guarded this one at all, which is how it came to promise a consequence that was false
# of every column it can name: "some display columns and **charts** will be blank". No chart
# reads `DP`, the allele columns or the normal read counts. The one column whose absence does
# blank a chart is `tumor_f` — the VAF histogram's input — and the note subtracts the arm's
# filter inputs, so that is precisely the column it never names.


class _RecordingStreamlit:
    """Just enough ``st`` for :func:`validate_required_columns`, recording what it said."""

    def __init__(self, sample_type):
        self.warnings: list[str] = []
        self.errors: list[str] = []
        self.session_state = type(
            "_S", (), {"filter_params": {"sample_type": sample_type}}
        )()

    def warning(self, text, *args, **kwargs):
        self.warnings.append(text)

    def error(self, text, *args, **kwargs):
        self.errors.append(text)


def _load_note(monkeypatch, sample_type, drop):
    """What the load-time note says about a reference MAF with ``drop`` taken out of it."""
    from utils import read_maf
    from page_modules import data_loading

    frame = read_maf(REFERENCE_MAF).drop(columns=list(drop))
    recorder = _RecordingStreamlit(sample_type)
    monkeypatch.setattr(data_loading, "st", recorder)

    assert data_loading.validate_required_columns(frame) is True, (
        "the MAF was refused, so there is no display-column note to read"
    )
    return [note for note in recorder.warnings if "Not in this MAF" in note]


@pytest.mark.parametrize(
    "sample_type, dropped",
    [
        ("somatic", ("DP", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2")),
        ("germline", ("DP", "n_alt_count", "n_ref_count")),
    ],
)
def test_the_load_note_names_the_columns_it_can_speak_for(monkeypatch, sample_type, dropped):
    """One note, naming every absent display column, on either arm.

    The columns are asserted by name rather than derived from the same expression the code
    uses, which would restate it: these are the six the ticket is about, and if the lists
    behind them change, the sentence a user reads should be looked at.
    """
    notes = _load_note(monkeypatch, sample_type, dropped)

    assert len(notes) == 1, f"expected one note, got {len(notes)}"
    for column in dropped:
        assert column in notes[0], f"{column} is absent and the note does not name it"


@pytest.mark.parametrize("sample_type", ["somatic", "germline"])
def test_the_load_note_says_the_columns_are_absent_not_empty(monkeypatch, sample_type):
    """"Empty in the table" was the first correction here, and it was false the same way.

    ``resolve_visible_columns`` ends ``[col for col in ordered if col in present]``, so a column
    the file does not carry is **dropped**: it is not an empty column, it is not a column. The
    old sentence promised a blank chart, the first replacement promised a blank cell, and only
    the third says what happens.

    Asserted against the resolver rather than against the sentence, so it is the behaviour that
    holds the copy: the columns the note names must be absent from what the resolver returns for
    a file lacking them.
    """
    from config.columns import MissingColumnsWarning, resolve_visible_columns

    dropped = (
        ("DP", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2")
        if sample_type == "somatic"
        else ("DP", "n_alt_count", "n_ref_count")
    )
    note = _load_note(monkeypatch, sample_type, dropped)[0]

    assert "empty in the table" not in note, (
        "the note says these columns will be empty, and the resolver drops them instead"
    )
    assert "neither the table nor" in note

    present = [
        column for column in resolve_visible_columns(sample_type) if column not in dropped
    ]
    with pytest.warns(MissingColumnsWarning):
        shown = resolve_visible_columns(sample_type, available_columns=present)
    still_there = sorted(column for column in dropped if column in shown)
    assert not still_there, (
        f"{still_there} survive into the resolver's list for a file that does not carry them, "
        "so the note is wrong to say they are in neither the table nor the download"
    )


def test_the_load_note_says_when_notes_go_with_the_column(monkeypatch):
    """The one absence here that costs more than a blank column, and only that one.

    Without ``Tumor_Seq_Allele2`` two variants at one position cannot be told apart, so
    :func:`components.variant_table._variant_key` returns ``None`` and the dialog offers no
    note field. A note saying only that the column will be "empty in the table" would be
    describing half of what happened, and the user meeting the missing field has nothing else
    to consult.

    The negative case is the half that keeps this honest: on a MAF missing the *other* display
    columns the sentence must not appear, or it becomes a warning about notes on every
    incomplete file.
    """
    with_allele = _load_note(monkeypatch, "somatic", ("DP", "Tumor_Seq_Allele2"))[0]
    assert "neither the note field nor" in with_allele, (
        "the alternate allele is absent, so notes are unavailable, and the note does not say so"
    )
    # Both go, because the dialog returns before either — a sentence naming only the note field
    # would describe half of what the user meets.
    assert "annotation columns" in with_allele

    without = _load_note(monkeypatch, "somatic", ("DP", "Tumor_Seq_Allele1"))[0]
    assert "note" not in without.lower(), (
        f"notes are unaffected by these absences and the note mentions them: {without!r}"
    )


@pytest.mark.parametrize("sample_type", ["somatic", "germline"])
def test_the_load_note_does_not_promise_a_blank_chart(monkeypatch, sample_type):
    """The false half of the old sentence, and the fact that makes it false.

    Two assertions rather than one, because dropping the word "chart" is not the same as the
    claim being true: the second requires that none of the columns this note can name appears in
    ``components/charts.py`` as a **string the module actually evaluates**. Broader than a
    subscript, because that module reads columns three ways — ``df["x"]``, ``data.get("x")`` and
    its own ``_safe_column(df, "x")``.

    Collected from the parsed AST rather than by grepping the text, which review caught as the
    weaker instrument: that file already carries comments quoting column names, so a comment
    added to explain *why no chart reads DP* would have failed the guard it was explaining.

    Mutation: put "charts will be blank" back → the first fails. Point a chart at ``DP`` →
    the second fails, and the note would then need to say so again.
    """
    dropped = (
        ("DP", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2")
        if sample_type == "somatic"
        else ("DP", "n_alt_count", "n_ref_count")
    )
    note = _load_note(monkeypatch, sample_type, dropped)[0]

    assert "chart" not in note.lower(), (
        "the note promises a blank chart for columns no chart reads"
    )

    import ast

    charts = ast.parse((STREAMLIT_APP / "components" / "charts.py").read_text())
    literals = {
        node.value
        for node in ast.walk(charts)
        if isinstance(node, ast.Constant) and isinstance(node.value, str)
    }
    drawn = sorted(column for column in dropped if column in literals)
    assert not drawn, (
        f"{drawn} now appear in components/charts.py, so an absence of them may well blank a "
        "chart — and the load note, which says only that the table and the export go empty, "
        "has stopped being the whole truth"
    )
