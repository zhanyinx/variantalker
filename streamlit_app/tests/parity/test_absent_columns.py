"""The one deviation from the pipeline the app takes on purpose, recorded as such.

Everything else in ``tests/parity/`` asks *does the app reach the pipeline's verdict*. This
file asks a different question, because on these files there is no pipeline verdict to reach:
the pipeline raises ``KeyError`` on a MAF missing a filter-input column, and issue #39 chose to
have the app fill the column neutrally and report anyway, so that an incomplete MAF is still
usable.

That is **off parity by construction**, and it lives in the parity suite rather than beside the
app-owned tests for one reason: a deviation recorded only where the app tests itself is a
deviation nobody comparing the two tools will ever find. Issue #20's handoff asked for exactly
this — fill-neutral runs excluded from the parity assertion *or* asserted against their own
expectations, never quietly absent from both.

What is asserted, per arm and per column:

* the pipeline really does refuse the file — measured, not assumed, because "the pipeline would
  have raised" is the entire justification for the deviation;
* the app really does report — every row back, labelled;
* the two are recognised as the third state, :attr:`~tests.parity.harness.ParityResult
  .off_parity_by_construction`, rather than as parity or as an unexplained divergence.

The fixtures are derived here rather than committed. A MAF missing a column is a MAF the
generator would have to build from the reference drive, and it carries no information the
committed one does not — dropping a column is a total function of the file already in git.
"""

from __future__ import annotations

import sys
import warnings
from dataclasses import replace
from pathlib import Path

import pytest

STREAMLIT_APP = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(STREAMLIT_APP))

from config.columns import MissingColumnsWarning, pipeline_output_columns  # noqa: E402
from filters.absent_columns import PATHOGENIC_INPUTS, REQUIRED_INPUTS  # noqa: E402
from tests.parity import harness as H  # noqa: E402
from tests.parity.contract import CASES_BY_NAME  # noqa: E402

#: The case each arm's derived fixtures are built from — contract defaults, no gene list, the
#: same parameters every other parity assertion runs under.
BASE_CASES = {
    "somatic": "somatic_defaults",
    "germline": "germline_defaults",
}


def _cases():
    return [
        pytest.param(arm, column, id=f"{arm}-{column}")
        for arm in BASE_CASES
        for column in REQUIRED_INPUTS[arm]
    ]


def _without(case_name: str, column: str, into: Path):
    """A copy of the case's fixture with one column dropped, and a Case pointing at it.

    The MAF's leading comment lines are carried over verbatim: ``read_maf`` counts them to find
    the header row, and both sides use ``read_maf``, so a file that lost them would be read
    differently rather than being the same file minus a column.

    ``Case.maf_path`` joins ``FIXTURE_DIR`` with ``fixture``, and joining an absolute path
    discards the left-hand side — so an absolute path in ``fixture`` is a fixture outside the
    committed directory, with no change to the contract needed.
    """
    case = CASES_BY_NAME[case_name]
    source = case.maf_path
    comments = []
    with source.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("#"):
                break
            comments.append(line)

    frame = H.read_maf(str(source)).drop(columns=[column])
    target = into / f"{source.stem}_no_{column}.maf"
    with target.open("w", encoding="utf-8") as handle:
        handle.writelines(comments)
        frame.to_csv(handle, sep="\t", index=False)

    return replace(case, name=f"{case.name}_no_{column}", fixture=str(target))


def _compare_without(arm: str, column: str, into: Path) -> H.ParityResult:
    """``H.compare`` on a derived fixture, with the *other* missing-column signal let through.

    Every column this file drops is also a column the pipeline's output carries, so
    ``config.columns.resolve_visible_columns`` warns that the table and the export will be one
    column short. That is a true and separate statement — it is about what the user can *see*,
    where this file is about what the filter *read* — and ``pytest.ini`` promotes it to an
    error so that no call site can forget to handle it.

    Downgraded rather than asserted here, because it fires inside the harness rather than in
    anything this file calls directly. That it still fires at all is asserted on its own, by
    ``test_the_column_resolver_says_its_own_half`` below.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", MissingColumnsWarning)
        return H.compare(_without(BASE_CASES[arm], column, into))


@pytest.mark.skipif(not H.PIPELINE_AVAILABLE, reason="needs bin/ for the pipeline side")
@pytest.mark.parametrize("arm,column", _cases())
def test_the_pipeline_refuses_the_file_the_app_fills(arm, column, tmp_path):
    """The premise, measured. Without it the deviation has nothing to be a deviation from.

    Asserted per column rather than once: "the pipeline dies on an incomplete MAF" is a claim
    about seven distinct columns per arm, and a column the pipeline turned out to tolerate
    would be one the app should not be filling at all — it would be inventing a divergence
    rather than working around one.
    """
    result = _compare_without(arm, column, tmp_path)

    assert result.pipeline_error is not None, (
        f"the pipeline reported on a MAF with no {column} column, so the app has no reason "
        "to deviate here: filling it invents a divergence instead of working around one"
    )
    assert "KeyError" in result.pipeline_error


@pytest.mark.skipif(not H.PIPELINE_AVAILABLE, reason="needs bin/ for the pipeline side")
@pytest.mark.parametrize("arm,column", _cases())
def test_a_filled_run_is_recorded_as_the_third_state(arm, column, tmp_path):
    """Not parity, not an unexplained divergence — off parity by construction.

    The harness must be able to say which of the three a case is in. If a filled run read as
    parity, the ratchet would accept any app verdict against a pipeline ``KeyError``; if it
    read as a plain divergence, this deliberate deviation would sit in the baseline looking
    like a bug someone forgot to fix, which is how the ESCAT blind spot survived.
    """
    result = _compare_without(arm, column, tmp_path)

    assert result.app_error is None, (
        f"the app failed on a MAF with no {column} column instead of filling it: "
        f"{result.app_error}"
    )
    assert result.app_filled_columns == [column]
    assert result.off_parity_by_construction
    assert not result.in_parity
    assert not result.errors_in_parity


@pytest.mark.parametrize("arm,column", _cases())
def test_the_column_resolver_says_its_own_half(arm, column):
    """The signal :func:`_compare_without` downgrades is still being produced.

    Two different things are missing on these files and the user needs both said. The filter
    lost an *input*, which this file is about; the table and the export lose an *output*, which
    ``config/columns.py`` is about. Suppressing the second inside the harness is only safe
    while something still asserts it fires, or the suppression quietly becomes the behaviour.

    Needs no pipeline: the resolver is a pure function of the column list.
    """
    pipeline_columns = pipeline_output_columns(sample_type=arm, skip_civic=False)
    if column not in pipeline_columns:
        pytest.skip(f"{column} is a filter input the pipeline does not emit")

    available = [name for name in pipeline_columns if name != column]
    with pytest.warns(MissingColumnsWarning, match=column):
        from config.columns import resolve_visible_columns

        resolve_visible_columns(
            sample_type=arm, skip_civic=False, available_columns=available
        )


@pytest.mark.skipif(not H.PIPELINE_AVAILABLE, reason="needs bin/ for the pipeline side")
@pytest.mark.parametrize("arm", sorted(BASE_CASES))
def test_a_complete_maf_is_not_in_the_third_state(arm):
    """The control. Every committed case must be judged by parity, not excused by this flag.

    Cheap and load-bearing: a bug that made ``filled_columns`` non-empty on a complete MAF
    would silently move every parity case into a state the ratchet does not assert.
    """
    result = H.compare(CASES_BY_NAME[BASE_CASES[arm]])

    assert result.app_filled_columns == []
    assert not result.off_parity_by_construction


@pytest.mark.skipif(not H.PIPELINE_AVAILABLE, reason="needs bin/ for the pipeline side")
@pytest.mark.parametrize("arm", sorted(BASE_CASES))
def test_filling_a_pathogenic_source_costs_the_report_rows(arm, tmp_path):
    """The direction the escalated warning exists to announce, measured against the pipeline.

    On the app's own suite this is asserted against the app's baseline; here it is asserted
    against the number the *pipeline* produced on the intact file, which is the report the user
    was actually sent. That is the comparison the warning is making on their behalf: not "fewer
    rows than another app setting" but "fewer rows than your report".
    """
    intact = H.compare(CASES_BY_NAME[BASE_CASES[arm]])
    assert intact.pipeline_error is None

    for column in PATHOGENIC_INPUTS[arm]:
        if column not in REQUIRED_INPUTS[arm]:
            # CIViC_Evidence_Level: guarded by the pipeline itself, never filled, on parity.
            continue
        filled = _compare_without(arm, column, tmp_path)
        assert filled.app_pass < intact.pipeline_pass, (
            f"filling {column} did not shrink the report below the pipeline's "
            f"{intact.pipeline_pass} rows, so the escalated warning would be announcing "
            "a loss that did not happen"
        )
