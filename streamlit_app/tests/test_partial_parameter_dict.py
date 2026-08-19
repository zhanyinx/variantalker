"""The parameters page invents an answer for a key nobody set (issue #276).

The reported symptom was *"the parameters page loses your selections"*, and the artifact
behind it was a parameter cache holding ``[]`` for every list parameter while every number
held a plausible value. That asymmetry read as a defect of the multiselect *read-back* —
something about how the list widgets are drawn or read that the number widgets escape.

It is not. Neither kind survived. The difference is only what each one invents when the key
is absent from ``filter_params``:

* a number reads ``.get(key, <literal>)``, so it invents its own hardcoded literal — which
  *looks* like a value the user chose, and for ``max_freq_population`` is ``0.01`` where
  both arms' contract says ``1.0``;
* a list reads ``.get(key)`` with **no fallback**, and ``filter_terms(None, options)`` is
  ``[]`` — so it invents *keep nothing*.

Then ``show_clinical_filters_tab`` writes what the widgets returned back into
``filter_params`` as though the user had chosen it, and the auto-save persists it. The page's
own written key set is a **fixed point**: a cache in this shape is re-written identically on
every later render, so it never heals and the app opens on nothing for good.

What makes a partial dict reach the page at all: ``load_parameters_from_cache`` returns the
cached parameters **verbatim** once the schema stamp matches (``param_store.py:95-101``). It
does not complete them against the arm's contract — unlike the upload route, whose
``migrate_params`` is documented to always produce *a complete parameter set for one arm*. So
any cache missing a key, whatever wrote it, is adopted with that key missing.

Measured on Streamlit 1.47.0. The three tests below all fail before the fix.

On the instrument
-----------------
Two ``AppTest`` results were discarded on the way here, and are recorded so nobody spends the
session re-finding them. ``unselect`` on a multiselect, and ``set_value([])``, can both report
the previous value as if the interaction were ignored — but the *same* probe behaves the same
way on a **keyed** multiselect, which a browser clears without trouble. So that is the
instrument, not the app, and no claim about the read-back is made from it. What is asserted
here needs no interaction at all: one render of a partial dict is enough.
"""

from __future__ import annotations

from contextlib import contextmanager
from pathlib import Path
from unittest import mock

import pytest

STREAMLIT_APP = Path(__file__).resolve().parents[1]

#: The page, seeded with a caller-supplied ``filter_params`` and nothing else. The seed is
#: installed directly rather than through the cache loader so the test names the state it is
#: testing, instead of depending on which arm this developer last configured.
_PAGE_SCRIPT = """
import sys
sys.path.insert(0, {app!r})
import streamlit as st
from page_modules import parameter_config

if "filter_params" not in st.session_state:
    st.session_state.filter_params = {seed!r}
parameter_config.show_parameter_config_page()
"""

#: What the dev's ``~/.mafigate/cached_parameters.json`` held on 2026-08-19T16:57:59, under
#: the current stamp (``schema_version`` 1, ``app_version`` 1.0.0). One render of
#: ``{{"sample_type": "germline"}}`` reproduces it exactly — including the *absence* of
#: ``filter_escat``, which the germline arm never draws and so never invents.
DEV_CACHE = {
    "sample_type": "germline",
    "min_depth": 50,
    "filter_variant_classification": [],
    "vaf_threshold_germline": 0.2,
    "skip_pathogenic": False,
    "filter_genes": [],
    "filter_intervar": [],
    "filter_renovo": [],
    "filter_clinvar": [],
    "max_freq_population": 0.01,
}


@contextmanager
def _page(seed):
    """One render of the whole parameters page, with the cache stubbed both ways.

    ``load_parameters_from_cache`` is stubbed for the reason ``test_parameter_adoption``
    gives — an unstubbed render consults ``~/.mafigate`` and lets the developer's own last
    session decide what is under test — and ``save_parameters_to_cache`` is stubbed because
    this file drives the auto-save deliberately: unstubbed, these tests would write the very
    poisoned cache they are about onto the machine running them.
    """
    pytest.importorskip("streamlit", reason="streamlit not installed")
    from streamlit.testing.v1 import AppTest

    from page_modules import param_store, parameter_config

    written = []

    with mock.patch.object(
        parameter_config, "load_parameters_from_cache", lambda: None
    ), mock.patch.object(
        parameter_config,
        "save_parameters_to_cache",
        lambda params: written.append(dict(params)) or True,
    ), mock.patch.object(
        parameter_config, "get_cache_info", lambda: None
    ), mock.patch.object(
        parameter_config, "clear_parameters_cache", lambda: True
    ), mock.patch.object(param_store, "discard_stale_cache", lambda: None):
        app = AppTest.from_string(_PAGE_SCRIPT.format(app=str(STREAMLIT_APP), seed=seed))
        app.run(timeout=60)
        assert not app.exception, [str(e.value) for e in app.exception]
        yield app, written


def _contract(arm):
    from config.pipeline_params import pipeline_params

    return pipeline_params(arm)


def test_a_dict_holding_only_the_arm_is_not_filled_with_empty_guideline_lists():
    """The reported loss, with no interaction in it at all.

    A user who never touches a control still ends the render with every guideline source
    emptied — which is not a narrower report but a *different* one: the criteria path empties
    and the report falls back to the pathogenic-rescue floor, which
    ``_warn_if_every_guideline_source_is_empty`` describes and does not prevent.
    """
    with _page({"sample_type": "germline"}) as (app, _written):
        params = app.session_state["filter_params"]
        contract = _contract("germline")

        emptied = [
            key
            for key in ("filter_clinvar", "filter_intervar", "filter_renovo")
            if contract.get(key) and not params.get(key)
        ]
        assert not emptied, (
            f"the page invented an empty keep-list for {emptied} from a dict that simply "
            f"did not carry them; the contract keeps "
            f"{ {k: contract[k] for k in emptied} }"
        )


def test_an_absent_number_is_not_invented_against_the_contract():
    """The numbers did not survive either — they were re-invented from literals.

    ``max_freq_population`` is the one that shows it, because the widget's own fallback
    (``0.01``) and the contract (``1.0``) disagree. That disagreement is what pinned the
    artifact's provenance: a cache holding ``0.01`` cannot have come from either arm's
    contract, so the value on screen was never a setting anybody chose.
    """
    with _page({"sample_type": "germline"}) as (app, _written):
        params = app.session_state["filter_params"]
        contract = _contract("germline")

        assert params["max_freq_population"] == contract["max_freq_population"], (
            "the page invented a population-frequency threshold from its own literal "
            f"({params['max_freq_population']}) rather than the contract's "
            f"({contract['max_freq_population']}), so a number the user never set reads "
            "on screen as one they did"
        )


def test_the_invented_emptiness_is_not_persisted_to_the_cache():
    """The step that makes one bad render permanent.

    ``filter_params.update(selected)`` writes the widgets' invented values back as though
    they were chosen, and the auto-save puts them on disk. Because the page's written key set
    is a fixed point, the file it writes here is exactly what it will read at the next boot —
    so the loss stops being a render and becomes the app's opening state.

    Whether a *repair* should also refuse to persist a wholesale collapse is deliberately not
    decided here; that belongs to the ticket this one blocks.
    """
    with _page({"sample_type": "germline"}) as (app, written):
        assert written, "the page did not auto-save, so this test is not measuring the write"
        document = written[-1]

        assert document != DEV_CACHE, (
            "one render of a dict holding nothing but the arm wrote the dev's poisoned "
            "cache byte-for-byte — every list emptied, every number an invented literal, "
            "and `filter_escat` absent because the germline arm never draws it"
        )


def test_guideline_lists_that_are_present_and_empty_stay_empty():
    """The other half of the fix, and the half a symptom-clearing repair would get wrong.

    **An absent key and a present-but-empty key are different states, and only the first is
    a defect.** An empty keep-list is legal and expressible — it is what the pipeline's own
    ``--filter_cancervar ""`` denotes — which is why issue #36 *deleted* the empty→``["All"]``
    backfill and why the page warns about an all-empty guideline block instead of refusing
    it. A completion that re-filled these from the contract would clear the reported symptom
    while overruling a choice the user is entitled to make, and would take the app's
    behaviour further from the pipeline's rather than closer.

    So this pins the distinction from both sides at once: the empties survive the render *and
    the write*, and the warning that describes what they do still draws.
    """
    from config.pipeline_params import pipeline_params

    emptied = pipeline_params("germline")
    for key in ("filter_clinvar", "filter_intervar", "filter_renovo"):
        emptied[key] = []

    with _page(emptied) as (app, written):
        params = app.session_state["filter_params"]
        for key in ("filter_clinvar", "filter_intervar", "filter_renovo"):
            assert params[key] == [], (
                f"{key} was chosen empty and came back as {params[key]!r}; completing an "
                "*absent* key from the contract must not re-fill a present one"
            )

        assert written, "the page did not auto-save, so the write is not being measured"
        for key in ("filter_clinvar", "filter_intervar", "filter_renovo"):
            assert written[-1][key] == [], (
                f"{key} was chosen empty and {written[-1][key]!r} was written to the cache"
            )

        warnings = [str(warning.value) for warning in app.warning]
        assert any("Every guideline source is empty" in text for text in warnings), (
            "the all-empty warning did not draw, so a user who has emptied every source is "
            f"not told the report has collapsed to the pathogenic-rescue floor; saw {warnings}"
        )
