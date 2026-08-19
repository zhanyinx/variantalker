"""One answer to "what does an absent parameter mean", on the path that filters.

Issue #289. ``filter_params`` was resolved **twice**, by two resolvers that disagreed, and
only one of them was behind the completion boundary issue #280 added:

* ``config.param_migration.complete_params``, at the top of ``show_parameter_config_page``,
  fills an absent key from **the arm's contract**;
* ``filters.variant_filters._Settings.from_params`` fills an absent key from **its own
  literal** — ``min_depth=0``, ``vaf_threshold=0.0``, every keep-list ``[]`` — and the data
  page, the report and the grids reached it without completing anything first.

Seven of the engine's literals contradict the contract on each arm. On the reference MAFs a
dict missing ``min_depth`` filtered **2 rows wider**, missing ``vaf_threshold`` 1 wider on
the somatic arm and **10 wider** on the germline one, and missing the last guideline source
its arm had **41 rows narrower**. The page and the engine were not two ways of saying the
same thing.

Which resolver wins, and why the engine keeps its literals
---------------------------------------------------------
The app completes; the engine is left alone. Three measurements decided it:

* ``filters/`` imports nothing from ``config/`` — deliberately, so the parity harness and a
  bare ``pytest`` can read the filter without booting a UI. Teaching the engine the app's
  contract inverts that edge.
* the engine's neutral literals are load-bearing for *isolation*: instrumenting one whole
  suite run counted **2083** calls into ``from_params``, **68** of them deliberately partial
  dicts across **37** unit tests, which is how a test of one filter layer says "nothing else
  is switched on". Completing there, or refusing there, breaks exactly those 37.
* refusing an incomplete dict at the engine cannot be right either: **every one of the four
  shipped presets** is incomplete by that measure — none carries ``skip_civic`` and all four
  spell ``keep_pathogenic`` — so the engine would reject a first-class route. Completion
  repairs all four, and is behaviour-preserving: ``+0`` rows on both reference MAFs for all
  four presets and both contracts.

So the boundary belongs where session state enters the filtering path, which is the top of
``show_data_loading_page`` — the mirror of the line the parameters page already runs.

The one key completion cannot supply
------------------------------------
``sample_type`` is the contract key ``complete_params`` deliberately never writes: arm
identity is not a filter setting, and writing it would make completion a path that moves a
user between arms without saying so (issue #280). Nothing else completed it either, and
every reader defaulted it to ``"somatic"`` — so a germline set that had lost its arm filtered
as somatic, silently, keeping **14** rows of ``germline_reference.maf`` where its own contract
keeps 59. That one is refused rather than completed, which is the only honest answer left.

Reachability, stated plainly: no route the app itself writes produces an armless dict today —
the parameters page always has an arm, ``migrate_params`` names a ``sample_type`` that is not
one, and issue #286 retires a cache that does not state an arm *before* it is restored. The
completion half of this, by contrast, was live: ``load_parameters_from_cache`` returns a
stamped cache verbatim by design, so a partial dict could reach the data page without the
parameters page ever rendering.

A unit module in the sense ``test_stale_banners.py`` means it: the pipeline has no session
state, so the *boundary* has no counterpart in ``bin/``. What the boundary must produce —
the contract's own numbers — is checked against ``apply_filters`` itself rather than against
literals typed in here, so this file cannot go green by agreeing with a stale copy.
"""

from __future__ import annotations

import ast
import copy
import sys
from pathlib import Path

import pytest

STREAMLIT_APP = Path(__file__).resolve().parent.parent
if str(STREAMLIT_APP) not in sys.path:
    sys.path.insert(0, str(STREAMLIT_APP))

from config.param_migration import complete_params  # noqa: E402
from config.pipeline_params import pipeline_params, stated_arm  # noqa: E402
from filters.variant_filters import MAFIGATE_FILTER, PASS, apply_filters  # noqa: E402
from page_modules.param_store import missing_contract_keys  # noqa: E402
from utils import read_maf  # noqa: E402

FIXTURES = STREAMLIT_APP / "tests" / "fixtures" / "parity"
GERMLINE_MAF = FIXTURES / "germline_reference.maf"

#: The keys dropped to make a partial set. All three are keys whose engine literal
#: contradicts the contract, and the germline arm is chosen because ``vaf_threshold`` is
#: where the two disagree most — 0.2 against the engine's 0.0.
DROPPED = ("min_depth", "vaf_threshold", "filter_variant_classification")

#: The page module's own file, read as text by the last section.
DATA_LOADING = STREAMLIT_APP / "page_modules" / "data_loading.py"

#: Every contract key, which is every key a defaulting read on the filter path could be
#: inventing. Read off the contract rather than restated, so a key added to it is covered.
CONTRACT_KEYS = frozenset(pipeline_params("somatic")) | frozenset(
    pipeline_params("germline")
)

#: The page, rendered over a seeded session with a chosen file — the sidebar's own hand-over
#: (issue #64), which is all the load path reads. Copied in shape from
#: ``test_data_page_sections.py``; the one difference that matters is that this script must
#: **not** seed ``filter_params`` itself, since the seed is what is under test.
_SCRIPT = """
import os
import sys
sys.path.insert(0, {app!r})

import streamlit as st
from components.sidebar import PENDING_UPLOAD_KEY
from page_modules import data_loading


class UploadStub:
    def __init__(self, path):
        self.name = os.path.basename(path)
        self._path = path

    def getvalue(self):
        with open(self._path, "rb") as handle:
            return handle.read()


st.session_state[PENDING_UPLOAD_KEY] = UploadStub({maf!r})

data_loading.show_data_loading_page()
"""


def _render(params):
    """The data page, rendered once over ``params``, with the germline MAF chosen."""
    pytest.importorskip("streamlit", reason="streamlit not installed")
    from streamlit.testing.v1 import AppTest

    app = AppTest.from_string(
        _SCRIPT.format(app=str(STREAMLIT_APP), maf=str(GERMLINE_MAF))
    )
    app.session_state["filter_params"] = params
    app.run(timeout=180)
    assert not app.exception, [str(e.value) for e in app.exception]
    return app


def _state(app, key, default=None):
    """``AppTest.session_state`` has no ``get``."""
    try:
        return app.session_state[key]
    except (KeyError, AttributeError):
        return default


def _partial_germline_set():
    """The germline contract, minus three keys the engine would answer differently."""
    params = pipeline_params("germline")
    for key in DROPPED:
        del params[key]
    return params


def _passing_rows(params):
    """How many rows of the germline MAF ``params`` keeps, straight through the filter."""
    labelled, _ = apply_filters(read_maf(str(GERMLINE_MAF)), copy.deepcopy(params))
    return int((labelled[MAFIGATE_FILTER] == PASS).sum())


def _messages(app):
    """Every error and warning the render drew, as plain strings."""
    return [element.value for element in app.error] + [
        element.value for element in app.warning
    ]


# ---------------------------------------------------------------------------
# The completion: the page filters with a complete set for its arm
# ---------------------------------------------------------------------------


def test_the_partial_set_and_the_contract_really_do_filter_differently():
    """The premise, first: without it every test below could pass for nothing.

    Both numbers come from ``apply_filters``, so this is the divergence as the engine
    actually resolves it today and not a remembered one.
    """
    partial = _partial_germline_set()
    completed = complete_params(copy.deepcopy(partial))

    assert _passing_rows(partial) != _passing_rows(completed), (
        "dropping "
        f"{list(DROPPED)} no longer changes what the filter keeps, so these guards can "
        "no longer tell a completed set from an uncompleted one — pick keys whose engine "
        "literal still contradicts the contract"
    )


def test_the_data_page_completes_a_partial_set_against_its_own_arms_contract():
    """Every contract key present after the render, and the arm untouched."""
    app = _render(_partial_germline_set())

    resolved = _state(app, "filter_params")
    assert missing_contract_keys(resolved) == (), (
        "the data page filtered with a set that is not complete for its arm, so the "
        f"engine's own literals answered for {missing_contract_keys(resolved)}"
    )
    contract = pipeline_params("germline")
    for key in DROPPED:
        assert resolved[key] == contract[key], (
            f"{key} was completed to {resolved[key]!r}, not to the contract's "
            f"{contract[key]!r}"
        )
    assert resolved["sample_type"] == "germline", (
        "completion moved the user between arms, which is the one thing it may never do"
    )


def test_the_data_page_cuts_the_report_its_contract_asks_for_not_the_engines():
    """The row count, which is what a user would have seen go wrong.

    Asserted against ``apply_filters`` on the completed set rather than a literal, so the
    number cannot go stale — and against the *uncompleted* count too, so a page that
    quietly stopped completing fails here rather than agreeing with itself.
    """
    partial = _partial_germline_set()
    expected = _passing_rows(complete_params(copy.deepcopy(partial)))

    app = _render(partial)

    filtered = _state(app, "filtered_data")
    assert filtered is not None, "the page drew no report at all"
    assert len(filtered) == expected, (
        f"the page kept {len(filtered)} rows where the germline contract keeps "
        f"{expected}; an uncompleted set keeps {_passing_rows(partial)}"
    )


# ---------------------------------------------------------------------------
# The refusal: the one key completion cannot supply
# ---------------------------------------------------------------------------


def test_the_data_page_refuses_a_set_that_does_not_state_an_arm():
    """No arm, no report — and the message says which control settles it.

    The alternative is what stood here: ``params.get("sample_type", "somatic")``, which
    filters a germline MAF on the somatic arm and says nothing. The count that route
    produced is asserted *against*, so a silent fallback reintroduced anywhere on this path
    fails this test rather than passing it quietly.
    """
    armless = pipeline_params("germline")
    del armless["sample_type"]

    app = _render(armless)

    assert _state(app, "filtered_data") is None, (
        "the page filtered a set that does not state an arm, so it picked one: "
        f"{len(_state(app, 'filtered_data'))} rows"
    )
    said = " ".join(_messages(app))
    assert "Sample Type" in said, (
        f"the refusal does not name the control that settles it: {said!r}"
    )


def test_the_data_page_refuses_a_sample_type_that_is_not_an_arm():
    """``"tumour"`` is refused by name, as a refusal and not as a filter error.

    It used to reach the engine, raise ``ValueError`` there and be caught two hundred lines
    below as ``❌ Error applying filters`` — accurate about where it happened and wrong
    about whose fault it was, since it sends the user off to change a filter.
    """
    misspelled = pipeline_params("germline")
    misspelled["sample_type"] = "tumour"

    app = _render(misspelled)

    said = " ".join(_messages(app))
    assert "tumour" in said, f"the refusal does not name the value it refused: {said!r}"
    assert "Error applying filters" not in said, (
        "a parameter this page can read for itself is still being reported as a failure of "
        f"the filter run: {said!r}"
    )
    assert _state(app, "filtered_data") is None, "a report was cut on an arm that is not one"


def test_stated_arm_never_answers_for_a_dict_that_states_nothing():
    """The shared reading, directly: it resolves an arm or it resolves nothing."""
    assert stated_arm({"sample_type": "germline"}) == "germline"
    assert stated_arm({"sample_type": "somatic"}) == "somatic"
    assert stated_arm({}) is None
    assert stated_arm({"sample_type": "tumour"}) is None
    assert stated_arm({"sample_type": None}) is None
    assert stated_arm(None) is None
    assert stated_arm("germline") is None


def test_both_boundaries_read_the_arm_the_same_way():
    """The cache's retirement and the page's refusal agree, because they share one reading.

    Not a restatement of the two tests above: what this pins is that the *same* dicts are
    unresolvable to both, so a cache the app declines to restore can never be a dict the
    page then filters with, or the other way round.
    """
    for params in ({}, {"sample_type": "tumour"}, {"sample_type": None}):
        assert stated_arm(params) is None
        assert missing_contract_keys(params) == ("sample_type",), (
            f"{params!r} is unresolvable to the page but not to the cache's retirement"
        )


# ---------------------------------------------------------------------------
# The shape, so a fourth resolver cannot appear in silence
# ---------------------------------------------------------------------------


def _defaulting_contract_reads(tree):
    """Every ``.get("<contract key>", <default>)`` in ``tree``, with its line.

    From the AST rather than by grep, because the single-argument form is fine — it asks
    whether a key is there — and only the two-argument form *invents* a value. The receiver
    is not inspected: in this module a contract key is only ever asked of the parameters.
    """
    found = []
    for inner in ast.walk(tree):
        if not isinstance(inner, ast.Call):
            continue
        if not (isinstance(inner.func, ast.Attribute) and inner.func.attr == "get"):
            continue
        if len(inner.args) != 2:
            continue
        first = inner.args[0]
        if isinstance(first, ast.Constant) and first.value in CONTRACT_KEYS:
            found.append((first.value, inner.lineno))
    return found


def test_the_data_page_never_answers_for_an_absent_contract_key():
    """Below the completion, a contract key is *indexed* — the issue #280 shape.

    A whole-file rule rather than a list of functions, because every function in this module
    runs inside ``show_data_loading_page``: the loads, the arm notice, the filter run and the
    report are all reached from it, so all of them are downstream of the boundary and none of
    them has anything to gain from a fallback. Three stood here before issue #289 —
    ``validate_required_columns``, ``current_arm`` and the filter run's ``circle_sources``
    call, that last one three hundred lines after the point where the dict is guaranteed to
    carry the keys it read with a default.

    ``KeyError`` is the intended outcome if the boundary is ever bypassed, exactly as it is
    for the parameter page's controls since issue #280 — a raise names the missing key, and a
    default hides it behind a report that looks fine.

    Deliberately *not* app-wide. ``MAFigate.initialize_session_state`` and
    ``parameter_config.adopt_parameters`` read ``sample_type`` with a default and keep it:
    both sit upstream of any completion, where the default seeds a set rather than answering
    for one already in play. The line this guard draws is the boundary, not the spelling.
    """
    tree = ast.parse(DATA_LOADING.read_text())

    invented = _defaulting_contract_reads(tree)
    assert invented == [], (
        f"{DATA_LOADING.name} answers for an absent contract key itself, which is the "
        f"second resolver issue #289 removed: {invented}. Index the key instead — the "
        "completion at the top of the page guarantees it."
    )
