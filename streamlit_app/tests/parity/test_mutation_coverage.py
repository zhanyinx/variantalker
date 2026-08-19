"""What the parity fixtures can still catch, asserted by breaking the app on purpose.

The suite's answer to "is this fixture set still worth anything?". ``test_parity.py``
asserts that the app and the pipeline agree; this asserts that the fixtures would
**notice if they stopped agreeing**, one known divergence at a time. See
``mutation_coverage.py`` for how the injection works and why the pipeline side is
untouched by construction.

It exists because the instrument it replaces stopped measuring. Issue #33 closed the last
divergence, so the recorded diff went empty, so no attribution probe had a witness, and
``test_attribution_coverage.py``'s soundness assertion — "every probe explains a diverging
row" — relaxed into its own negation: with nothing diverging it asserted that *no* case
carries an attribution block, and passed. Any fixture set reaching row parity went green
there, including one with no coverage at all. Issue #242 removed it and put this in its
place.

**This module is the safety net for the fixture swap** (issue #246). Its claim is not
"the fixtures contain the right shapes" — shape was measured and found not to be coverage
(the fixture README's "Why 'the column is populated' was not coverage") — but "some case
notices each divergence coming back".

Runtime **~11 s** measured (``tests/parity/`` totals ~31 s for 252 tests). Almost all of
it is the 33 pipeline subprocesses, which run once and are cached per fixture directory;
the ~800 in-process app runs across both sweeps cost so little that the falsifiability
control — a second full measurement, with a divergence's witnesses withheld — is
effectively free, which is why it re-measures rather than deriving its answer from the run
that included them. Marked ``integration`` like the rest of the parity suite and
deliberately **not** ``slow``: ``make test-fast`` deselects that marker, and this is
exactly the check someone wants before changing fixtures.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

STREAMLIT_APP = Path(__file__).resolve().parent.parent.parent
if str(STREAMLIT_APP) not in sys.path:
    sys.path.insert(0, str(STREAMLIT_APP))

from tests.parity import harness as H  # noqa: E402
from tests.parity import mutation_coverage as M  # noqa: E402
from tests.parity.contract import CASES, CASES_BY_NAME  # noqa: E402

pytestmark = [
    pytest.mark.integration,
    pytest.mark.skipif(
        not H.PIPELINE_AVAILABLE,
        reason="pipeline bin/ not present (packaged app or partial checkout)",
    ),
]


@pytest.fixture(scope="module")
def report() -> M.CoverageReport:
    """One measurement of the committed fixture set, shared by every assertion below."""
    return M.measure()


#: The fixture README's baseline table, transcribed: ``case -> (PASS, NOPASS, columns)``.
#: Twelve of its thirteen rows; the thirteenth, ``dot_numeric``, has no verdict counts to
#: transcribe and is asserted separately below.
#:
#: Transcribed from the README rather than read out of ``baseline.json``, and that is the
#: point of this table existing at all: it is an *independent* statement of what the two
#: sides do on these files, so reproducing it is evidence that the instrument is driving
#: the real harness rather than a mock of it. Reading the numbers out of the same file
#: ``test_baseline_matches`` reads would prove only that two tests agree.
#:
#: A fixture swap must move these numbers and the README's table together — a different
#: fixture set has different counts, and neither of them is allowed to change in silence.
README_BASELINE_TABLE = {
    "somatic_defaults": (60, 22, 40),
    "somatic_skip_patho": (54, 28, 40),
    "somatic_genes": (58, 24, 40),
    "somatic_genes_mixed_case": (58, 24, 40),
    "germline_defaults": (59, 35, 39),
    "germline_skip_patho": (51, 43, 39),
    "germline_genes": (57, 37, 39),
    "somatic_synthetic": (7, 11, 40),
    "germline_synthetic": (7, 10, 39),
    "civic_present": (8, 4, 46),
    "civic_skipped": (0, 12, 40),
    "gnomad_genome_present": (0, 6, 41),
}


# ---------------------------------------------------------------------------
# The instrument is measuring the real harness
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("case_name", sorted(README_BASELINE_TABLE))
def test_the_instrument_reproduces_the_fixture_readmes_baseline_table(report, case_name):
    """Every PASS / NOPASS / output-column count the fixture README publishes.

    The validation issue #242 asks for before anything downstream is believed. An
    instrument that mocked either side, or drifted onto different parameters, would still
    happily report full mutation coverage — of a harness nobody runs.
    """
    measured = report.table[case_name]
    expected = README_BASELINE_TABLE[case_name]
    assert (measured["pass"], measured["nopass"], measured["output_columns"]) == expected, (
        f"{case_name}: measured {measured}, README table says "
        f"PASS={expected[0]} NOPASS={expected[1]} cols={expected[2]}. Either the "
        "instrument is not running the real harness, or the fixtures changed and the "
        "README's table and this one both need re-recording."
    )


def test_the_unreadable_fixture_is_the_readmes_thirteenth_row(report):
    """``dot_numeric`` has no PASS count in the README, and must have none here either.

    The row the table records as "pipeline raises ``TypeError``; app refuses". Asserted
    separately rather than folded in as a sentinel triple, because what makes it parity is
    the *shape* of the two non-verdicts, not a count.
    """
    row = report.table["dot_numeric"]
    assert row["pipeline_error"].startswith("TypeError")
    assert row["app_error"] == "UnreadableNumericColumns"
    assert row["app_refused_columns"] == ["t_alt_count", "t_ref_count", "tumor_f"]
    assert row["agree"] is True


def test_no_case_is_off_parity_before_any_mutation(report):
    """A case that already disagrees cannot report whether it noticed a mutation.

    The instrument excludes such a case and names it, so this failing means the coverage
    numbers below were measured over fewer cases than the set contains — read this before
    reading them.
    """
    assert not report.disagreeing_unmutated, (
        f"cases off parity before any mutation was injected: {report.disagreeing_unmutated}"
    )


# ---------------------------------------------------------------------------
# Coverage
# ---------------------------------------------------------------------------


def assert_covered(outcome: M.MutationOutcome) -> None:
    """The suite's verdict on one divergence. Shared with the falsifiability control.

    Factored out so the control below can assert that *this assertion* fails on a report
    with a gap, rather than only that the report records one. A guard nobody has seen fail
    is a guard nobody knows works.
    """
    assert outcome.covered, (
        f"{outcome.mutation}: no case notices this divergence when it is re-injected into "
        f"the app ({outcome.why}). Silent on {sorted(outcome.silent_on)}; not applicable "
        f"to {sorted(outcome.not_applicable)}"
        + (
            f"; the app crashed rather than diverging on {sorted(outcome.caught_by_crash)}"
            if outcome.caught_by_crash
            else ""
        )
        + ". The fixture set has lost the rows that witness it — add a cell to "
        "build_parity_fixtures.py pinning a row's whole path to the verdict."
    )


@pytest.mark.parametrize("mutation_name", sorted(m.name for m in M.MUTATIONS))
def test_every_scoreable_divergence_is_caught(report, mutation_name):
    """13 of 13 on the committed fixtures. One red light per divergence, named."""
    assert_covered(report.outcomes[mutation_name])


def test_the_committed_set_scores_thirteen_of_thirteen(report):
    """The headline, as a literal rather than as an expression over the same data.

    ``score`` is derived from the outcomes, so comparing it to a formula over
    ``len(MUTATIONS)`` could only ever restate itself. The number issue #242 asks for is
    **13/13**, so that is what is written here: dropping a mutation from the table lowers
    the denominator and fails, which is the direction a derived assertion could not see.
    """
    assert report.score == "13/13", f"gaps: {report.gaps}"
    assert len(report.outcomes) == len(M.MUTATIONS) == 13
    assert len({m.name for m in M.MUTATIONS}) == 13, "two mutations share a name"


def test_coverage_is_verdict_level_rather_than_crash_level(report):
    """No divergence may be covered only because the mutated app blew up.

    A crash proves the mutation ran, not that the fixtures discriminate: the same crash
    would be reported by a file with no interesting rows in it at all. The instrument
    already keeps crashes out of :attr:`MutationOutcome.covered`; this says the committed
    set does not need them, so a future fixture set cannot quietly start relying on them.
    """
    crashed = {
        name: sorted(outcome.caught_by_crash)
        for name, outcome in report.outcomes.items()
        if outcome.caught_by_crash
    }
    assert not crashed, (
        f"mutations that make the app raise rather than diverge: {crashed}. Either the "
        "mutation is malformed, or the app has a latent crash on these rows."
    )


# ---------------------------------------------------------------------------
# The denominator
# ---------------------------------------------------------------------------


def test_every_one_of_the_twelve_divergences_is_accounted_for():
    """Scored, unscoreable-by-name, or measured directly — exactly once each.

    The guard on the denominator. A divergence dropped from the mutation table would
    otherwise raise the score by shrinking the map: 12/12 of eleven divergences reads
    exactly like full coverage.
    """
    scored = {m.divergence for m in M.MUTATIONS}
    named = scored | set(M.UNSCOREABLE) | set(M.MEASURED_DIRECTLY)
    assert named == set(M.THE_TWELVE), (
        f"divergences with no account: {sorted(set(M.THE_TWELVE) - named)}; "
        f"accounted for but not in the map: {sorted(named - set(M.THE_TWELVE))}"
    )
    overlap = scored & (set(M.UNSCOREABLE) | set(M.MEASURED_DIRECTLY))
    assert not overlap, f"divergence both scored and excused: {sorted(overlap)}"


def test_the_two_unscoreable_divergences_are_named_and_not_counted(report):
    """#5 and #10 are excused explicitly, and the excuse says why.

    Neither can be scored by mutation — #5 has no MAF shape, and no parity case exercises
    #10 at all, since every case runs at the neutral ``max_freq_population = 1.0``. What
    this pins is that they are *named* and stay outside the tally: an unscoreable
    divergence dropped from the report would leave the reader counting a full house, and
    one counted as caught would inflate the score by two.
    """
    assert set(M.UNSCOREABLE) == {"#5", "#10"}
    assert not {"#5", "#10"} & {m.divergence for m in M.MUTATIONS}
    assert not {"#5", "#10"} & {o.divergence for o in report.outcomes.values()}

    rendered = M.render(report)
    for number in ("#5", "#10"):
        assert M.THE_TWELVE[number] in rendered, f"{number} is not named in the report"
        assert M.UNSCOREABLE[number][:40] in rendered


def test_column_parity_is_measured_directly_and_the_widths_are_exercised(report):
    """#11, by the harness's prefix rule rather than by mutation.

    Both halves matter. Parity says the app's list opens with the pipeline's, position by
    position. The **widths** say the fixture set still reaches every output shape: 39 and
    40 for the two arms, 41 for the ``gnomAD_genome`` column-presence probe, and 46 for
    the CIViC fixture — the only case that has ever emitted all 45 ``KEEP`` entries.
    """
    mismatched = sorted(
        name for name, row in report.column_parity.items() if not row["in_parity"]
    )
    assert not mismatched, f"the app's columns do not open with the pipeline's: {mismatched}"
    assert report.widths == [39, 40, 41, 46], (
        f"output widths exercised: {report.widths}. A set that no longer reaches 46 has "
        "lost the KEEP-completeness the CIViC fixture bought; one that no longer reaches "
        "41 has lost compute_keep's gnomAD_genome branch."
    )


# ---------------------------------------------------------------------------
# The instrument can report a gap
# ---------------------------------------------------------------------------


def _thinnest_witnessed(full: M.CoverageReport) -> tuple[str, list[str]]:
    """The divergence with the fewest witnesses, and all of them. Deterministic by name.

    *All* of them, not one: a control that withheld a single case would stop working the
    day the fixture set improved enough to give every divergence a second witness — it
    would break on good news, which is the wrong way round for a safety net.
    """
    witnessed = sorted(
        (len(outcome.caught_by), name, sorted(outcome.caught_by))
        for name, outcome in full.outcomes.items()
        if outcome.caught_by
    )
    assert witnessed, "nothing is covered at all, so there is no witness to withhold"
    _, name, witnesses = witnessed[0]
    return name, witnesses


def test_withholding_a_divergences_witnesses_makes_the_instrument_report_a_gap(report):
    """Take the witnesses away and the instrument goes 12/13, red, and names which one.

    Issue #242's third acceptance criterion, and the reason to believe the other twelve
    green lights above. The whole sweep is re-measured over the remaining cases — all
    thirteen mutations, really re-run — rather than recomputed from the run that included
    the witnesses, so the number below is a measurement and not an inference about how the
    measurement composes.

    The expected gap set is derived from the full run (a mutation gaps exactly when every
    case that caught it was withheld) rather than hardcoded, so the control keeps its
    meaning if the fixtures change. Today that derivation is a single divergence, and the
    literal it produces is the 12/13 the ticket asks to see.
    """
    mutation_name, witnesses = _thinnest_witnessed(report)
    remaining = [case for case in CASES if case.name not in set(witnesses)]

    without = M.measure(cases=remaining)

    expected_gaps = sorted(
        name
        for name, outcome in report.outcomes.items()
        if not set(outcome.caught_by) - set(witnesses)
    )
    assert mutation_name in expected_gaps
    assert without.gaps == expected_gaps, (
        f"withholding {witnesses} was expected to cost exactly {expected_gaps}, but the "
        f"re-measurement reports {without.gaps} — the per-case measurements are not "
        "independent, or the full run under-reported."
    )
    assert without.score == f"{len(M.MUTATIONS) - len(expected_gaps)}/{len(M.MUTATIONS)}"
    assert without.score == "12/13", (
        f"the control is meant to demonstrate 12/13; it produced {without.score} because "
        f"withholding {witnesses} cost {expected_gaps}"
    )

    with pytest.raises(AssertionError, match="no case notices this divergence"):
        assert_covered(without.outcomes[mutation_name])


def test_the_vendored_symbols_are_restored_after_every_mutation(report):
    """No mutation may outlive the measurement that installed it.

    A leaked mutation would not fail here — it would fail in whichever parity module ran
    next, as a divergence nobody introduced. The restoration is ``finally``-guaranteed;
    this checks the guarantee held over a real run, by identity against the vendored
    module the app imports from. It takes ``report`` so that the full sweep has actually
    happened by the time the identities are compared.
    """
    assert report.outcomes, "no mutation ran, so there is nothing to have restored"

    import filters.variant_filters as V
    from vendor import pipeline_filters as P

    assert V.pipeline_common_filters is P.common_filters
    assert V.pipeline_somatic_filters is P.somatic_filters
    assert V.pipeline_germline_filters is P.germline_filters


def test_a_mutation_is_installed_only_inside_its_context():
    """And the context really installs one — otherwise every mutation would read silent.

    The other direction of the test above, and the cheaper of the two failures to make:
    a context manager that patched nothing would leave the instrument reporting a fixture
    set that catches nothing, which reads identically to a fixture set that is empty.
    """
    import filters.variant_filters as V

    mutation = next(m for m in M.MUTATIONS if m.target == M.TARGET_COMMON)
    before = V.pipeline_common_filters
    with M.mutated(mutation):
        assert V.pipeline_common_filters is not before
    assert V.pipeline_common_filters is before


def test_pointing_the_instrument_elsewhere_puts_the_fixture_directory_back():
    """``--fixture-dir`` must not outlive the measurement that asked for it.

    The flag exists so a *candidate* fixture set can be measured on the same mutations
    (issues #233 and #246), and it works by rebinding a module global in two modules. Left
    rebound inside a pytest session it would hand every later parity assertion a different
    set of MAFs, which is a failure that would surface far from its cause.
    """
    from tests.parity import contract as C

    committed = C.FIXTURE_DIR
    with M.fixtures_from(Path("/nonexistent/fixture/set")):
        assert C.FIXTURE_DIR != committed
        assert H.FIXTURE_DIR != committed
    assert C.FIXTURE_DIR == committed
    assert H.FIXTURE_DIR == committed


def test_the_case_set_the_score_is_measured_over_is_the_whole_contract(report):
    """33 cases, not a subset someone narrowed while debugging."""
    assert sorted(report.cases) == sorted(CASES_BY_NAME)
