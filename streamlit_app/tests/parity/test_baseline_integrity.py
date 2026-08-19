"""What the committed baseline must say, read with no ``bin/`` anywhere in sight.

The pipeline-less half of the parity suite. ``test_parity.py`` and
``test_mutation_coverage.py`` both need the pipeline checkout and skip without it, which
issue #24 identified as the place the parity suite silently disappears; this module runs
off ``baseline.json`` alone, so the claims below hold in a pipeline-less CI job and in a
partial checkout.

It is what remains of ``test_attribution_coverage.py``, which issue #242 removed. That
module asserted two halves of one property — *completeness* (no diverging row is
unattributed) and *soundness* (every probe explains a row that actually diverges) — and
the soundness half had inverted. Issue #33 closed the last divergence, so the recorded
diff went empty, so no probe could have a witness, so the assertion was relaxed to hold
"while divergence remains" and in its absence asserted instead that **no case carries an
attribution block at all**: the negation of what it was written to detect, passing on any
fixture set that reaches row parity, including one with no coverage in it. Coverage is now
measured by mutation (``mutation_coverage.py``), which is a measurement rather than a
restatement of parity.

The completeness half did not invert, and is kept. It is *dormant* on the committed
baseline — with an empty diff there is nothing to leave unattributed — and dormant is a
state this repo has learned not to trust on a docstring's word, so each check below is
paired with a test that runs the same predicate over a hand-made record and requires it to
fire. A guard nobody has seen fail is a guard nobody knows works.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

STREAMLIT_APP = Path(__file__).resolve().parent.parent.parent
if str(STREAMLIT_APP) not in sys.path:
    sys.path.insert(0, str(STREAMLIT_APP))

from tests.parity import harness as H  # noqa: E402

_BASELINE = H.load_baseline()["cases"]


def diverging_cases(baseline: dict) -> list[str]:
    """The cases a baseline records as not in row parity."""
    return sorted(name for name, record in baseline.items() if not record.get("in_parity"))


def unattributed_cases(baseline: dict) -> dict[str, int]:
    """Cases whose diverging rows no probe explains, and how many rows each."""
    return {
        name: (record.get("attribution") or {})["unattributed"]
        for name, record in baseline.items()
        if "unattributed" in (record.get("attribution") or {})
    }


def unknown_attribution_keys(baseline: dict) -> list[str]:
    """Attribution keys no probe emits — a rename that left its old bucket behind."""
    known = set(H.PROBES) | {"unattributed"}
    seen = {
        key
        for record in baseline.values()
        for key in (record.get("attribution") or {})
    }
    return sorted(seen - known)


# ---------------------------------------------------------------------------
# The committed baseline
# ---------------------------------------------------------------------------


def test_the_fixtures_record_no_remaining_row_divergence():
    """Issue #33's destination, read straight off the committed file.

    The assertion ``test_parity.py`` names as the reason its own green is not
    green-by-skipping: with no ``bin/`` the whole of that module skips, and this says the
    recorded result of the last full run was parity on every case.
    """
    diverging = diverging_cases(_BASELINE)
    assert not diverging, (
        f"{len(diverging)} case(s) are no longer in row parity: {diverging}. Row parity "
        "was reached on all 33 by issue #33; if this is a deliberate change re-record the "
        "baseline, otherwise it is a regression."
    )


def test_no_diverging_row_is_unattributed():
    """Completeness: every diverging row is explained by some known divergence.

    An ``unattributed`` entry is not a failure of the app — it means the map does not yet
    name whatever is moving those rows, which is the most informative state the harness can
    report. It goes red here so it is investigated rather than accumulated.

    Dormant while the diff is empty, and it is the *diff* being empty that makes it so —
    :func:`test_the_fixtures_record_no_remaining_row_divergence` is what stops that from
    being an unchecked excuse, and :func:`test_the_unattributed_check_can_fail` is what
    stops the predicate from being unfalsifiable in the meantime.
    """
    unattributed = unattributed_cases(_BASELINE)
    assert not unattributed, (
        f"{unattributed} diverging rows no probe explains. Either the map is missing a "
        "divergence, or PROBES is missing its probe — measure it as a candidate in "
        "reference.py before editing harness.PROBES, so the instrument is not tuned to "
        "the bucket it found."
    )


def test_baseline_covers_every_probe_and_nothing_else():
    """The two vocabularies do not drift apart.

    Every attribution key is either a probe name or ``unattributed``; a leftover key from a
    renamed probe would otherwise sit in the baseline forever, counted by nobody and read
    as coverage by anyone skimming the file.
    """
    unknown = unknown_attribution_keys(_BASELINE)
    assert not unknown, f"baseline.json holds attribution keys no probe emits: {unknown}"


# ---------------------------------------------------------------------------
# The checks above, made to fail
# ---------------------------------------------------------------------------


def test_the_divergence_check_can_fail():
    """The predicate finds a diverging case when there is one to find.

    Cheap, and the only thing standing between "no case diverges" and a predicate that
    returns the empty list whatever it is handed — which is how the module this replaces
    stopped meaning anything.
    """
    assert diverging_cases({"a": {"in_parity": True}, "b": {"in_parity": False}}) == ["b"]


def test_the_unattributed_check_can_fail():
    """An ``unattributed`` bucket is found, and an explained case is not mistaken for one."""
    baseline = {
        "explained": {"attribution": {"#6 germline_escat_unmirrored": 4}},
        "not_explained": {"attribution": {"unattributed": 3}},
    }
    assert unattributed_cases(baseline) == {"not_explained": 3}


def test_the_vocabulary_check_can_fail():
    """A key no probe emits is reported, while a real probe name and ``unattributed`` pass.

    Both halves: a check that flagged everything would be as useless as one that flagged
    nothing, and it would be discovered only by whoever next added an attribution key.
    """
    probe = sorted(H.PROBES)[0]
    baseline = {"a": {"attribution": {probe: 1, "#99 renamed_away": 2, "unattributed": 1}}}
    assert unknown_attribution_keys(baseline) == ["#99 renamed_away"]


@pytest.mark.parametrize(
    "predicate", [diverging_cases, unattributed_cases, unknown_attribution_keys]
)
def test_every_predicate_reads_an_empty_baseline_as_clean(predicate):
    """The other direction: nothing fires on nothing, so a red above means a real record."""
    assert not predicate({})
