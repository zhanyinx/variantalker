"""The single output-column resolver, `config.columns.resolve_visible_columns`.

Three separately drifted copies of the pipeline's column list used to answer "which
columns does the app show?", and they disagreed:

* `config/columns.py:get_filtered_columns` — a hand-maintained 44-entry transcription
  of `bin/utils.py:KEEP` that no display code ever called;
* `_DEFAULT_VISIBLE_COLUMNS`, then in `components/ui_components.py` and now in neither
  that file nor the tree — the 25 columns that actually drove the grid and the "shown
  columns" export. The grid is `components/variant_table.py` since issue #76;
* a `KEEP` in `utils/main_utils.py`, already deleted by issue #32.

This module pins the replacement. The resolver builds the pipeline's list by calling
the vendored `compute_keep`, then appends the app's own extras, so the app's list is a
deliberate **superset** with the pipeline's list as an exact prefix.

The prefix is asserted element-by-element, never by length. The drifted lists were once
measured at 40 against the pipeline's 40 while differing by a substitution, so a length
check is exactly the check that misses this bug.
"""

import sys
import warnings
from pathlib import Path

import pytest

STREAMLIT_APP = Path(__file__).resolve().parents[1]
if str(STREAMLIT_APP) not in sys.path:
    sys.path.insert(0, str(STREAMLIT_APP))

from config.columns import (  # noqa: E402
    APP_EXTRA_COLUMNS,
    COLUMN_DESCRIPTIONS,
    COLUMN_SOURCE_HEADER,
    COLUMNS_MAFIGATE_ADDS,
    COLUMNS_THE_MAF_CARRIES_OUTSIDE_THE_KEEP_LIST,
    OPTIONAL_COLUMNS,
    PIPELINE_COLUMNS_THE_APP_REPLACES,
    REQUIRED_COLUMNS,
    MissingColumnsWarning,
    _COLUMN_SOURCE_PROSE,
    _MISSING_DESCRIPTION,
    _classify_column_source,
    create_column_info_table,
    get_column_description,
    pipeline_output_columns,
    resolve_visible_columns,
)
from vendor.pipeline_utils import KEEP  # noqa: E402


def shown_pipeline_columns(arm: str, skip_civic: bool = False, available=None) -> list:
    """The pipeline's columns the app actually shows: its list, less the replaced ones.

    **Independent in one half and not the other, and the docstring has to say which.**
    The *order and membership* of the pipeline's list is built here from
    `pipeline_output_columns`, so a resolver that reordered or dropped anything fails the
    tests below. Which names are excused is read from
    `PIPELINE_COLUMNS_THE_APP_REPLACES` — the same constant the resolver reads — so on that
    half this is `X` against `X`, and **no test in this module can fail on a name being
    added to that constant**.

    That is deliberate rather than overlooked, and something else holds it: the parity
    baseline records `pipeline_columns_not_shown` by name, so adding one fails
    `tests/parity/test_parity.py::test_baseline_matches` on 32 cases until someone
    regenerates the baseline deliberately. The unit test says the resolver obeys the list;
    the committed baseline is what makes the list itself hard to grow.

    An earlier draft of this docstring claimed the whole thing was independent. It was the
    shape issue #77 found and issue #100 found again in a key that compared itself with
    itself — claimed *about* a guard this time rather than built into one, which is the
    version that survives a green test run.
    """
    return [
        col
        for col in pipeline_output_columns(arm, skip_civic, available)
        if col not in PIPELINE_COLUMNS_THE_APP_REPLACES
    ]

#: The gnomAD spellings a real MAF carries and `compute_keep` appends.
GNOMAD_RAW = ["gnomAD_exome_AF_raw", "gnomAD_genome_AF_raw"]

#: gnomAD **v4** column names. `_DEFAULT_VISIBLE_COLUMNS` listed both, and both are
#: fiction here: neither appears in a reference MAF on either genome build, and the
#: vendored computation only ever appends the four older `gnomAD_{exome,genome}_AF[_raw]`
#: spellings — so even a MAF that carried them could not get them into the list.
GNOMAD_V4_NAMES = ["gnomAD_exome_AF_grpmax", "gnomad41_genome_AF"]

#: Columns the app derives for display, which no MAF supplies.
DERIVED_COLUMNS = ["Clinical_Summary", "Pathogenicity_Overview", "Notes", "Sample_Name"]


def maf_columns(arm: str = "somatic") -> list:
    """The columns a filtered frame carries: the pipeline's `KEEP` plus MAF extras."""
    columns = list(KEEP) + GNOMAD_RAW + ["dbSNP_ID"]
    if arm == "germline":
        columns += ["InterVar", "RENOVO_Class", "RENOVO_pls"]
    return columns


@pytest.fixture(autouse=True)
def keep_is_pristine():
    """The shared `KEEP` must survive every resolver call unmutated.

    `compute_keep` copies before it mutates; the resolver must not undo that by
    handing the list on to something that mutates in place.
    """
    before = list(KEEP)
    yield
    assert KEEP == before, "vendor.pipeline_utils.KEEP was mutated by a resolver call"


# ---------------------------------------------------------------------------
# The prefix property — the whole point of the ticket
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("arm", ["somatic", "germline"])
@pytest.mark.parametrize("skip_civic", [False, True])
def test_app_list_opens_with_the_pipeline_list_as_an_exact_prefix(arm, skip_civic):
    """Element-by-element, not by length. See the module docstring.

    Since issue #117 the prefix is the pipeline's list **less**
    :data:`PIPELINE_COLUMNS_THE_APP_REPLACES`. That is a weaker claim than the one this
    test made before, so the strength has to come back somewhere: it comes back in
    :func:`test_the_only_pipeline_columns_missing_from_the_view_are_the_named_ones`,
    which pins the difference to exactly that list. The two together say what this one
    used to say alone — every pipeline column is shown, in order, except the ones a
    ticket has argued for by name.
    """
    available = maf_columns(arm) + DERIVED_COLUMNS
    pipeline = shown_pipeline_columns(arm, skip_civic, available)
    app = resolve_visible_columns(arm, skip_civic, available)

    assert app[: len(pipeline)] == pipeline, (
        f"{arm}/skip_civic={skip_civic}: the app's list does not open with the "
        f"pipeline's.\n  pipeline: {pipeline}\n  app:      {app[: len(pipeline)]}"
    )
    assert len(app) > len(pipeline), "the app contributes no extras at all"


@pytest.mark.parametrize("arm", ["somatic", "germline"])
def test_pipeline_prefix_is_the_vendored_computation_verbatim(arm):
    """`pipeline_output_columns` must be `compute_keep`, not a transcription of it."""
    from types import SimpleNamespace

    import pandas as pd

    from vendor.pipeline_filters import compute_keep

    available = maf_columns(arm)
    frame = pd.DataFrame({name: pd.Series(dtype="object") for name in available})
    expected = compute_keep(SimpleNamespace(sample_type=arm, skip_civic=False), frame)

    assert pipeline_output_columns(arm, False, available) == expected


@pytest.mark.parametrize("arm", ["somatic", "germline"])
def test_extras_follow_the_pipeline_columns_rather_than_mixing_in(arm):
    """No app extra may appear before the last pipeline column."""
    available = maf_columns(arm) + DERIVED_COLUMNS
    pipeline = pipeline_output_columns(arm, False, available)
    app = resolve_visible_columns(arm, False, available)

    for extra in app[len(pipeline) :]:
        assert extra not in pipeline, f"{extra} is both a pipeline column and an extra"
    assert set(app[: len(pipeline)]).isdisjoint(app[len(pipeline) :])


# ---------------------------------------------------------------------------
# What the tail may and may not contain
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", GNOMAD_V4_NAMES)
def test_the_gnomad_v4_names_are_not_in_the_app_tail(name):
    """They match no reference MAF, and `compute_keep` could not append them if they did."""
    assert name not in APP_EXTRA_COLUMNS


@pytest.mark.parametrize("name", GNOMAD_V4_NAMES)
def test_a_gnomad_v4_name_present_in_the_data_still_does_not_appear(name):
    """Even handed the column, the resolver must not adopt it.

    The pipeline half is `compute_keep`, which hardcodes the older spellings; the app
    half is a fixed list. Neither reacts to a v4 name, and this pins that.
    """
    available = maf_columns("somatic") + DERIVED_COLUMNS + [name]
    assert name not in resolve_visible_columns("somatic", False, available)


def test_verdict_and_reason_are_the_last_two_columns():
    """The two columns the app's own filter writes, appended after everything else."""
    from filters.variant_filters import MAFIGATE_FILTER, MAFIGATE_REASON

    available = maf_columns("somatic") + DERIVED_COLUMNS + [MAFIGATE_FILTER, MAFIGATE_REASON]
    app = resolve_visible_columns("somatic", False, available)

    assert app[-2:] == [MAFIGATE_FILTER, MAFIGATE_REASON]


def test_the_extras_are_named_by_one_list_the_resolver_actually_uses():
    """`APP_EXTRA_COLUMNS` must be the tail, not documentation of it."""
    available = maf_columns("somatic") + DERIVED_COLUMNS + APP_EXTRA_COLUMNS
    pipeline = shown_pipeline_columns("somatic", False, available)
    app = resolve_visible_columns("somatic", False, available)

    assert app[len(pipeline) :] == APP_EXTRA_COLUMNS


# ---------------------------------------------------------------------------
# The one subtraction: issue #117
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("arm", ["somatic", "germline"])
@pytest.mark.parametrize("skip_civic", [False, True])
def test_the_only_pipeline_columns_missing_from_the_view_are_the_named_ones(arm, skip_civic):
    """The app may show less than the pipeline emits **only** by that named list.

    This is the guard that keeps the weakened prefix property honest. Written as an
    equality on the difference rather than as "the replaced column is absent", because
    the claim worth holding is not that one name is gone — it is that *nothing else* is,
    and an absence test cannot say that.

    It compares against the constant, so it holds the resolver to the list and **not the
    list to anything** — see :func:`shown_pipeline_columns` for what does that, and why
    this module cannot.
    """
    available = maf_columns(arm) + DERIVED_COLUMNS
    emitted = pipeline_output_columns(arm, skip_civic, available)
    app = resolve_visible_columns(arm, skip_civic, available)

    dropped = [col for col in emitted if col not in app]
    assert dropped == [c for c in PIPELINE_COLUMNS_THE_APP_REPLACES if c in emitted], (
        f"{arm}/skip_civic={skip_civic}: the view drops pipeline columns nothing has "
        f"argued for. Dropped: {dropped}"
    )


@pytest.mark.parametrize("column", PIPELINE_COLUMNS_THE_APP_REPLACES)
def test_the_replaced_column_is_still_what_the_pipeline_emits(column):
    """Taking a column off the screen must not change the app's account of the pipeline.

    `pipeline_output_columns` feeds the drift guard, the missing-column account and
    `circle_sources`. If the subtraction reached it, the app would stop being able to say
    what the pipeline's own output is — and the drift guard would go quiet about exactly
    the column this ticket touched.

    It doubles as the check that a name here is not a typo: a name the pipeline emits on
    neither arm subtracts nothing, and would sit in this list looking like a decision.
    *At least one* arm rather than both, because that is the claim — a future entry may
    well be somatic-only, and requiring both would refuse it for the wrong reason.
    """
    emitting = [arm for arm in ("somatic", "germline") if column in pipeline_output_columns(arm)]

    assert emitting, (
        f"{column} is in no arm's pipeline output, so replacing it is a no-op and the "
        "name is probably a typo"
    )


def test_a_maf_without_the_replaced_column_is_not_warned_about_it():
    """Its absence is not a loss to a table that was never going to show it.

    Measured on the real corpus at issue #117: 74 of the 176 annotated MAFs carry no
    `variantalker_naive`, and every one of them drew a warning naming it as a column
    missing from "neither the table nor the export" — of a column the app has now decided
    is in neither by choice.
    """
    available = [
        c
        for c in maf_columns("somatic")
        if c not in PIPELINE_COLUMNS_THE_APP_REPLACES and c != "Otherinfo"
    ]

    with pytest.warns(MissingColumnsWarning) as caught:
        resolve_visible_columns("somatic", False, available)

    said = str(caught[0].message)
    assert "Otherinfo" in said, "the warning stopped reporting a genuinely absent column"
    for column in PIPELINE_COLUMNS_THE_APP_REPLACES:
        assert column not in said, f"{column} is named as missing from a view it is not in"


@pytest.mark.parametrize("column", PIPELINE_COLUMNS_THE_APP_REPLACES)
def test_a_column_taken_off_the_screen_can_still_be_looked_up(column):
    """The glossary is what makes leaving it out honest rather than silent.

    It is one tick of *Show all columns* away, so a reader can still meet it; this table
    is the only thing the app says about it now. The description must be a real one — the
    `get_column_description` fallback is a string, so a test that only checks for truth
    passes over "Description not available for variantalker_naive".
    """
    table = create_column_info_table()
    row = table[table["Column Name"] == column]

    assert not row.empty, f"{column} is in neither the view nor the column reference"
    assert "not available" not in row.iloc[0]["Description"]
    assert row.iloc[0]["Category"] == "Clinical Significance", (
        "a column the app replaces is one it answers itself, so it belongs beside the "
        "sources it summarises rather than in the chain's 'Other' fallthrough"
    )

    # Which arms the table says it belongs to, derived from which arms actually emit it
    # rather than written down — this row is built from the pipeline's list, not the
    # resolver's, and that substitution is exactly the thing to get wrong.
    emitting = {arm for arm in ("somatic", "germline") if column in pipeline_output_columns(arm)}
    expected = {
        frozenset({"somatic", "germline"}): "Both",
        frozenset({"somatic"}): "Somatic",
        frozenset({"germline"}): "Germline",
    }[frozenset(emitting)]
    assert row.iloc[0]["Sample Types"] == expected


def test_the_two_clinical_verdicts_are_both_described_and_say_they_are_the_same_question():
    """What the app says to tell them apart, which is the whole of issue #117's copy.

    Both columns had **no** entry at all before it: the column reference rendered
    "Description not available for Clinical_Summary" for the column the table pins first.
    Asserted through `get_column_description`, the accessor the page calls, rather than
    against the dict — the page renders what the accessor returns.
    """
    summary = get_column_description("Clinical_Summary")
    pipeline_copy = get_column_description("variantalker_naive")

    for text in (summary, pipeline_copy):
        assert "not available" not in text
    assert "Clinical_Summary" in pipeline_copy, (
        "the pipeline column's entry must name the column it duplicates — that pairing is "
        "the only thing on any surface that says they answer one question"
    )
    assert summary != pipeline_copy


def test_no_circle_source_is_a_column_the_app_takes_off_the_screen():
    """The Pathogenicity Overview may only summarise a column the reader can go and check.

    `circle_sources` draws a circle for a source whose column the arm emits and this MAF
    carries — and its docstring promised that was "the same test `resolve_visible_columns`
    applies", so the summarised column is always in the table beside it. Since issue #117
    it is *not* the same test: the resolver applies it and then subtracts. The promise
    survives only because no circle source is in that list, which is a fact about the list.

    So it is asserted rather than assumed. Putting a source in
    `PIPELINE_COLUMNS_THE_APP_REPLACES` would give the strip a circle summarising a column
    the user cannot open — issue #95's defect exactly, arriving from the other direction:
    that ticket stopped the overview naming sources the file lacks, and this stops it
    naming one the app has hidden.
    """
    from components.clinical_summary import CLINICAL_CIRCLE_SOURCES

    summarised = {column for _abbreviation, _display, column in CLINICAL_CIRCLE_SOURCES}
    hidden = summarised & set(PIPELINE_COLUMNS_THE_APP_REPLACES)

    assert not hidden, (
        f"{sorted(hidden)} would be summarised by a circle over a column the table no "
        "longer shows, so the reader has nothing to check it against"
    )


def test_every_column_the_reference_table_lists_has_a_real_description():
    """The guard issue #120 exists for, and the reason the fallback may stay a string.

    `get_column_description` returns a plausible sentence for a name nobody described, so
    the column reference rendered "Description not available for Pathogenicity_Overview"
    about the column the grid pins **second** — and six rows of fifty-seven said it before
    anyone measured. The fallback is kept (a `KeyError` from an accessor the help page
    calls with its own curated shortlists would take out the whole page over one missing
    sentence), so this is what makes a missing entry loud instead.

    Built from the table rather than from a list of names, which is the whole point: the
    table is derived from the resolver, so a new `APP_EXTRA_COLUMNS` entry joins it by
    itself and fails here until someone describes it. A hand-listed set would have to be
    updated by the person adding the column, and that is precisely what did not happen.

    Compared against the exact fallback rather than sniffing for "not available", because a
    real description may legitimately contain that phrase — one describing a column that is
    blank when an annotation is absent very nearly does.

    This subsumes a second guard that was written beside it and then deleted: one looping
    `APP_EXTRA_COLUMNS` to make the point that the app's *own* columns have no explanation
    anywhere else, so this table is the only place they can be looked up. Every one of those
    names is a row of this table, so it could not fail on its own — issue #77's shape, a
    guard that only ever agrees with another. The argument was the value, and it is written
    here instead.
    """
    table = create_column_info_table()
    undescribed = sorted(
        row["Column Name"]
        for _, row in table.iterrows()
        if row["Description"] == _MISSING_DESCRIPTION.format(column=row["Column Name"])
    )

    assert not undescribed, (
        f"the column reference tells the user {undescribed} is undescribed, in a table "
        "whose own introduction promises every column MAFigate can show them"
    )


def test_every_app_extra_column_is_declared_on_one_side_or_the_other():
    """The guard issue #124 exists for, and the reason its two halves are lists.

    `COLUMNS_MAFIGATE_ADDS` and `COLUMNS_THE_MAF_CARRIES_OUTSIDE_THE_KEEP_LIST` must
    partition `APP_EXTRA_COLUMNS`. Neither is inferred as the complement of the other,
    because inferring it is exactly what would let a seventh app extra be classified by
    accident — and being classified by accident is the whole bug: the six columns the app
    derives are in neither `REQUIRED_COLUMNS` nor `OPTIONAL_COLUMNS`, so they fell off the
    end of the old chain and the column reference called them "Optional".

    So a new entry fails **here**, where the question "does the file carry this, or do we
    make it?" is asked once, rather than rendering a wrong claim on a clinician's screen.
    """
    declared = COLUMNS_MAFIGATE_ADDS + COLUMNS_THE_MAF_CARRIES_OUTSIDE_THE_KEEP_LIST

    assert len(declared) == len(set(declared)), (
        "a column is declared both as one MAFigate adds and one the file carries"
    )
    assert set(declared) == set(APP_EXTRA_COLUMNS), (
        "every app extra must say which kind it is: "
        f"undeclared {sorted(set(APP_EXTRA_COLUMNS) - set(declared))}, "
        f"declared but not an app extra {sorted(set(declared) - set(APP_EXTRA_COLUMNS))}"
    )


def test_every_token_is_worded_and_every_wording_is_reached():
    """The classifier's tokens and their prose must be the same set, in both directions.

    This asserted that every `REQUIRED_COLUMNS` **category** had wording, because the
    classifier returned a category name as its own token. Issue #127 ended that: the status
    now comes from what an absence does, so only `core` is still both a category and a token,
    and the old assertion would pass over four wordings nothing can reach.

    Both directions matter and they fail differently. A token with no wording raises
    `KeyError` out of `column_source_status`, inside the Help page. A wording nothing reaches
    is the quieter failure: it reads as a live answer the table can give, and #126 found
    exactly that shape in a dead `OPTIONAL_COLUMNS` loop returning what its own fallback
    returned.

    **The tokens are read out of the classifier's own source**, and it took two goes to get
    there. Drawing the subjects from the reference table was the first attempt: building that
    table calls `column_source_status`, so an unworded token raised `KeyError` *before* the
    assertion below could run and the "unworded" half was unreachable — a guard passing on the
    strength of an exception raised elsewhere. Feeding it one name per branch was the second,
    and review caught what it still could not see: a *new* branch returning a new token has no
    subject that reaches it, so the case this test exists for stayed invisible.

    Parsing `_classify_column_source` for its `return` literals has neither hole. It enumerates
    the branches rather than sampling them, so a sixth branch is in `tokens` the moment it is
    written, whether or not anyone adds a column that reaches it.
    """
    import ast
    import inspect

    from config import columns as columns_module

    classifier = next(
        node
        for node in ast.walk(ast.parse(inspect.getsource(columns_module)))
        if isinstance(node, ast.FunctionDef) and node.name == "_classify_column_source"
    )
    tokens = {
        node.value.value
        for node in ast.walk(classifier)
        if isinstance(node, ast.Return)
        and isinstance(node.value, ast.Constant)
        and isinstance(node.value.value, str)
    }
    assert len(tokens) >= 3, (
        f"only {sorted(tokens)} was parsed out of _classify_column_source, so this guard has "
        "stopped reading the function and is checking almost nothing"
    )

    unworded = sorted(tokens - set(_COLUMN_SOURCE_PROSE))
    assert not unworded, (
        f"_classify_column_source returns {unworded}, which _COLUMN_SOURCE_PROSE does not "
        "word, so the column reference raises KeyError for every column in it"
    )

    unreached = sorted(set(_COLUMN_SOURCE_PROSE) - tokens)
    assert not unreached, (
        f"_COLUMN_SOURCE_PROSE words {unreached}, which no branch of _classify_column_source "
        "returns — so the table cannot give that answer and the wording claims it can"
    )


def test_no_filter_input_is_called_optional():
    """The claim issue #127 exists for, in the words the user reads.

    "Optional" says *your MAF may or may not carry this, and MAFigate copes either way*. For a
    column the filter reads on every path, coping means **filling it with a stand-in value and
    handing back a report that is not a complete result** — issue #39 measured that at up to
    70% of the somatic report — so "Optional" is the wrong answer and "Filled if absent" is
    the one the table has always given the five clinical annotations.

    The subjects come from the derivation over the vendored source, and the thing under test is
    the classifier: two different sources, which is #124's rule for a guard. Sourcing them from
    a list in `config/columns.py` is what the deleted `clinical` mirror did, and a mirror can
    only be as right as whoever last edited it.

    `core` is allowed and is not an exception being smuggled in: `Variant_Classification` is
    both core and a filter input, and a file missing it is refused before the filter could fill
    anything, so the refusal is the true answer and the stronger one.
    """
    from filters.absent_columns import REQUIRED_INPUTS

    inputs = {column for columns in REQUIRED_INPUTS.values() for column in columns}
    assert inputs, "the derivation found no filter inputs, so this guard checks nothing"

    optional_sentence = _COLUMN_SOURCE_PROSE["optional"]
    table = create_column_info_table()
    misdescribed = sorted(
        row["Column Name"]
        for _, row in table.iterrows()
        if row["Column Name"] in inputs
        and row[COLUMN_SOURCE_HEADER] == optional_sentence
    )

    assert not misdescribed, (
        f"the column reference calls {misdescribed} optional, and the filter reads every one "
        "of them on every path — an absent one is filled with a stand-in value and the report "
        "is not a complete result, which is not what 'optional' tells a clinician"
    )


def test_no_column_mafigate_derives_is_called_optional():
    """The false claim itself, in the words the user reads.

    "Optional" is a statement about the *user's file* — your MAF may or may not carry this,
    and MAFigate copes either way. It was false in both directions for all six columns the
    app derives: a MAF never carries them, and the app always makes them.

    Built from the table, like #120's description guard and for the same reason: the table
    is derived from the resolver, so a new derived column joins it by itself.

    **The set checked is the complement, not `COLUMNS_MAFIGATE_ADDS` itself** — this test
    was written the obvious way first and was #77-vacuous under the exact mutation its own
    docstring advertised. Filtering the table's rows by membership of the list under test
    means deleting `Notes` from that list also deletes it from the filter: it went back to
    the fallthrough, rendered "Optional", and this still passed. So the question asked here
    is the one the user's file answers — *is this an app extra that the MAF does not
    carry?* — sourced from `APP_EXTRA_COLUMNS` less the names declared file-carried, which
    is a different statement from the one being guarded.
    """
    must_not_be_optional = set(APP_EXTRA_COLUMNS) - set(
        COLUMNS_THE_MAF_CARRIES_OUTSIDE_THE_KEEP_LIST
    )
    table = create_column_info_table()
    said_optional = sorted(
        row["Column Name"]
        for _, row in table.iterrows()
        if row["Column Name"] in must_not_be_optional
        and row[COLUMN_SOURCE_HEADER] == _COLUMN_SOURCE_PROSE["optional"]
    )

    assert not said_optional, (
        f"the column reference tells the user {said_optional} may or may not be in their "
        "MAF, about columns no MAF carries and MAFigate makes on every report"
    )


def test_dbsnp_id_is_answered_as_a_column_the_file_carries():
    """The one app extra "Optional" was always right about, asserted rather than argued.

    The app appends `dbSNP_ID` because `compute_keep` does not, not because it makes it —
    the file supplies it — and #120 measured 167 of 290 real MAFs carrying it. An absent
    app extra is dropped from the view in silence, which is what "optional" means at its
    strongest. Nothing bound that to what the table draws until now: the partition guard
    only said `dbSNP_ID` was declared on *some* side.
    """
    table = create_column_info_table()
    row = table[table["Column Name"] == "dbSNP_ID"]

    assert len(row) == 1, "dbSNP_ID is no longer a row of the column reference"
    assert row.iloc[0][COLUMN_SOURCE_HEADER] == _COLUMN_SOURCE_PROSE["optional"]


def test_the_two_optionals_are_two_facts():
    """Why the classifier returns tokens and the prose lives in a separate mapping.

    `get_column_requirement_status` returned "Optional" from two places — a genuine
    `OPTIONAL_COLUMNS` member, and the final fallback for a name in no category at all —
    and two different facts spelled the same way is what let the derived columns hide among
    them. They still *render* the same sentence, deliberately: for a column the file
    supplies, "you do not have to carry this" is true either way. What must not come back
    is the two being indistinguishable in code, and no assertion over the rendered string
    can tell them apart, so this asserts over the tokens.
    """
    assert _classify_column_source("CIViC_Evidence_Level") == "optional"
    assert _classify_column_source("cosmic") == "uncategorised"
    assert _COLUMN_SOURCE_PROSE["optional"] == _COLUMN_SOURCE_PROSE["uncategorised"], (
        "if these ever diverge, this test is the place to decide what each should say"
    )

    # And the fallthrough that carried the bug is closed to the columns it was wrong about:
    # nothing the app derives may reach it, whatever else does. Sourced from the complement
    # for the reason `test_no_column_mafigate_derives_is_called_optional` sets out — looping
    # `COLUMNS_MAFIGATE_ADDS` here asks the list under test to nominate its own members.
    for name in set(APP_EXTRA_COLUMNS) - set(COLUMNS_THE_MAF_CARRIES_OUTSIDE_THE_KEEP_LIST):
        assert _classify_column_source(name) == "derived", (
            f"{name} is an app extra nobody declared the file carries, so the reference "
            "table must not answer for it as though the MAF might supply it"
        )


def test_the_column_reference_header_asks_a_question_both_kinds_of_row_can_answer():
    """The rename, held against the values drawn under it.

    A fifth value under the old "Required" header would have answered a question the header
    did not ask. The header was widened instead, so this asserts the pairing rather than
    either half: the frame carries `COLUMN_SOURCE_HEADER` as its own key — which is what
    `page_modules/help.py` configures the column by — and every cell under it is prose the
    module owns rather than a string built at the call site.
    """
    table = create_column_info_table()

    assert COLUMN_SOURCE_HEADER in table.columns
    assert "Required" not in table.columns

    drawn = set(table[COLUMN_SOURCE_HEADER])
    assert drawn <= set(_COLUMN_SOURCE_PROSE.values()), (
        f"the table draws {sorted(drawn - set(_COLUMN_SOURCE_PROSE.values()))}, which "
        "nothing in config.columns worded"
    )
    # Both kinds of row are actually present, so the pairing is being exercised rather than
    # asserted over a table that happens to contain only one of them.
    assert _COLUMN_SOURCE_PROSE["derived"] in drawn
    assert _COLUMN_SOURCE_PROSE["core"] in drawn


def test_the_verdict_and_reason_entries_name_every_value_they_can_hold():
    """A glossary entry for a coded column is only usable if it decodes the codes.

    Derived from the filter's own constants, not from a list typed out beside them: a fifth
    reason cell, or a renamed verdict, fails here until the description catches up. That is
    the failure this ticket is about arriving one layer in — a column whose *values* nobody
    described reads exactly like a column whose values are self-evident.

    The raw spellings are in the prose deliberately. Issue #79's bounded exception admits
    raw MAF column names to this table so a user can match it against their own header row;
    these are values the app itself puts on that user's screen, and the argument is the
    same one — the glossary is where you look up what you are looking at.
    """
    from filters.variant_filters import (
        MAFIGATE_FILTER,
        MAFIGATE_REASON,
        NOPASS,
        PASS,
        _CELL_REASONS,
    )

    verdict = get_column_description(MAFIGATE_FILTER)
    # Counted, not `in`. `PASS` is a substring of `NOPASS`, so the obvious spelling of this
    # check is half vacuous: an entry that mentioned only NOPASS would satisfy both
    # iterations of `for value in (PASS, NOPASS): assert value in verdict`. The count
    # separates them — two occurrences of the shorter, one of the longer.
    assert verdict.count(PASS) == 2, (
        f"the verdict column's entry does not say what both {PASS} and {NOPASS} mean"
    )
    assert verdict.count(NOPASS) == 1

    reason = get_column_description(MAFIGATE_REASON)
    for value in _CELL_REASONS.values():
        assert value in reason, (
            f"{value} can appear in this column and the glossary does not mention it"
        )


def test_the_overview_entry_points_at_the_key_rather_than_copying_it():
    """Say what the column is; let the key say what it drew (issue #120's second bullet).

    The key beside the grid names every glyph it can draw and every source it drew them
    for, and **both halves are derived** — per arm, per file — precisely because a
    hand-written copy had drawn ⚪ beside a key that did not list it. A glossary that
    transcribed either half would be a third surface that cannot be right for every report
    even on the day it is written, which is how `docs/CLINICAL_SUMMARY_FEATURE.md` came to
    hold three wrong lines.

    So this asserts the *absence* of a copy, against the constants the key is drawn from.
    The entry still has to say something, which the check above holds.
    """
    from components.clinical_summary import (
        CLINICAL_CIRCLE_SOURCES,
        pathogenicity_circle_legend,
    )

    entry = get_column_description("Pathogenicity_Overview")

    copied_glyphs = [glyph for glyph, _label in pathogenicity_circle_legend() if glyph in entry]
    assert not copied_glyphs, (
        f"{copied_glyphs} is transcribed into the glossary, so it can outlive the key it "
        "was copied from — which is exactly how a circle came to be drawn unexplained"
    )

    copied_sources = [
        display
        for _abbreviation, display, _column in CLINICAL_CIRCLE_SOURCES
        if display in entry
    ]
    assert not copied_sources, (
        f"{copied_sources} is named in the glossary, but which sources a report draws "
        "depends on the arm and on the file, so a fixed list here is wrong for some report"
    )

    assert "key" in entry.lower(), (
        "the entry drops the copy without saying where the reader should look instead"
    )


def test_the_two_pinned_clinical_columns_are_filed_together():
    """The category filter must not split what the grid pins side by side.

    Issue #117 moved `Clinical_Summary` out of the chain's "Other" fallthrough and left its
    twin behind, because it was only categorising the columns it was describing. Describing
    `Pathogenicity_Overview` is what brings the other half: the two summarise the same
    sources for the same variant, are pinned first and second, and a user filtering this
    table by Clinical Significance was getting one of them.

    Read off the built table, not off the `if/elif` chain, because what the chain assigns
    and what the page shows are the same thing only if you read the thing the page renders.
    """
    table = create_column_info_table()
    categories = {
        column: table[table["Column Name"] == column].iloc[0]["Category"]
        for column in ("Clinical_Summary", "Pathogenicity_Overview")
    }

    assert categories["Pathogenicity_Overview"] == categories["Clinical_Summary"], (
        f"the two pinned verdict columns are filed apart, as {categories} — a reader "
        "filtering by one category meets one of the pair and not the other"
    )
    assert categories["Clinical_Summary"] == "Clinical Significance"


def test_the_derived_sample_column_is_filed_with_the_barcodes_it_is_derived_from():
    """`Sample_Name` is `Tumor_Sample_Barcode` with a fallback, and was filed under "Other".

    So Sample Information listed the two barcodes and not the column the grid actually
    shows the user — the same fallthrough as the pair above, found by the same measurement.

    Asserted against the barcode's own row rather than against the literal "Sample
    Information", which is the claim the name of this test makes: renaming the category
    should move both or fail, not quietly leave the derived column behind again.
    """
    table = create_column_info_table()
    categories = {
        column: table[table["Column Name"] == column].iloc[0]["Category"]
        for column in ("Sample_Name", "Tumor_Sample_Barcode")
    }

    assert categories["Sample_Name"] == categories["Tumor_Sample_Barcode"], (
        f"the derived sample column is filed apart from the barcode it is derived from, "
        f"as {categories}"
    )


def test_the_replaced_columns_are_not_also_claimed_as_app_extras():
    """A name cannot be both subtracted from the pipeline half and appended in the tail."""
    assert set(PIPELINE_COLUMNS_THE_APP_REPLACES).isdisjoint(APP_EXTRA_COLUMNS)
    for column in PIPELINE_COLUMNS_THE_APP_REPLACES:
        assert column in COLUMN_DESCRIPTIONS


# ---------------------------------------------------------------------------
# Missing columns: intersect and warn, never crash
# ---------------------------------------------------------------------------


def test_a_maf_missing_an_expected_column_warns_rather_than_crashing():
    available = [c for c in maf_columns("somatic") if c != "Otherinfo"]

    with pytest.warns(MissingColumnsWarning, match="Otherinfo"):
        columns = resolve_visible_columns("somatic", False, available)

    assert "Otherinfo" not in columns
    assert columns, "the resolver returned nothing at all"


def test_every_returned_column_is_actually_present():
    """The caller indexes a frame with this list, so an absent name is a KeyError."""
    available = maf_columns("germline")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", MissingColumnsWarning)
        columns = resolve_visible_columns("germline", False, available)

    assert set(columns) <= set(available)


def test_absent_app_extras_are_dropped_without_a_warning():
    """The extras are derived at display time, so their absence is normal, not drift.

    Warning about them would fire on every call made before the display columns are
    computed — which is most of them — and drown the signal the warning exists for.
    """
    available = maf_columns("somatic")  # carries no derived columns

    with warnings.catch_warnings():
        warnings.simplefilter("error", MissingColumnsWarning)
        columns = resolve_visible_columns("somatic", False, available)

    assert "Clinical_Summary" not in columns


def test_no_available_columns_means_assume_everything_is_present():
    """The help pages describe the schema with no MAF in hand."""
    with warnings.catch_warnings():
        warnings.simplefilter("error", MissingColumnsWarning)
        columns = resolve_visible_columns("somatic", False, None)

    assert "Clinical_Summary" in columns
    for name in ("Hugo_Symbol", "CIViC_Evidence_Level", "cosmic"):
        assert name in columns


# ---------------------------------------------------------------------------
# Arm and CIViC behaviour, inherited from the vendored computation
# ---------------------------------------------------------------------------


def test_skip_civic_drops_the_civic_columns():
    available = maf_columns("somatic")
    columns = resolve_visible_columns("somatic", True, available)

    assert not [c for c in columns if c.startswith("CIViC")]


def test_germline_columns_do_not_leak_into_a_later_somatic_call():
    """The aliasing bug `compute_keep` fixes, re-checked through the resolver.

    The resolver is what the app calls, so a fix that only holds one layer down is
    not a fix the app benefits from.
    """
    somatic_available = maf_columns("somatic")
    first = resolve_visible_columns("somatic", False, somatic_available)
    resolve_visible_columns("germline", False, maf_columns("germline"))
    second = resolve_visible_columns("somatic", False, somatic_available)

    assert first == second
    assert "Tumor_Sample_Barcode" in second
    assert "InterVar" not in second


def test_two_germline_calls_in_a_row_do_not_raise():
    """The second call used to raise `ValueError: list.remove(x): x not in list`."""
    available = maf_columns("germline")
    assert resolve_visible_columns("germline", False, available) == (
        resolve_visible_columns("germline", False, available)
    )


def test_germline_substitutes_intervar_for_cancervar():
    columns = resolve_visible_columns("germline", False, maf_columns("germline"))

    assert "InterVar" in columns
    assert "CancerVar" not in columns
