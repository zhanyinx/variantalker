"""The app opens at parity, and the guideline sentinel is gone (issue #36).

Two changes, one theme: the controls stop lying about what they do.

**The classification control is an exclude list**, because that is what the pipeline's
``filter_var_classification`` means. The include-list complement was rejected on
measurement, not on taste: 211 reference rows carry five classifications the app's
hardcoded vocabulary has never heard of, two of them minted by the pipeline itself, and
the parameter page runs *before any MAF is loaded* — so the vocabulary can never be
data-driven and an include list is permanently off parity by whatever it has not heard
of. An exclude list keeps those rows, as the pipeline does.

**The catch-all sentinel is deleted**, from the classification key and from all six
guideline keys. Reading ``All`` as "no restriction" was measured to be *identical to the
bug it was meant to avoid*: three of the sources match a listed value in 100% of
reference rows, so one saturated source made the OR true for every row — 68,364 somatic
rows against a parity baseline of 411, a 167-fold widening reached by selecting the
option that reads as "don't filter".

What replaces it is nothing at all: an empty multiselect *is* the pipeline's empty
keep-list, and ``isin([])`` drops that source out of the union. Neutral, not narrowing,
and not saturating.

What is asserted here, and where
--------------------------------
The unreachability of the widening is asserted three ways, because one alone would be
weak:

1. **Structurally** — the sentinel is in no option list, so it cannot be selected.
2. **Behaviourally, on the value** — passing the literal string ``"All"`` as a keep-list
   (which a stale ``~/.mafigate`` cache still can) produces the *same verdict on every
   row* as an empty list. There is no branch left for it to trigger.
3. **Behaviourally, on the bound** — no reachable guideline selection makes the report
   larger than the contract's, and the all-empty state bottoms out at the
   pathogenic-rescue floor rather than at the row count.

Numbers are quoted per *cell*, never as a bare total: the somatic reference is
legitimately both 408 and 411 depending on whether the criteria path or the union is
meant, and conflating the two is how a real discrepancy survived several tickets
(issue #28).

The widget tests drive the real page through ``streamlit.testing.v1.AppTest``, following
``test_gene_lists.py``. They are what proves the *control* changed, as against the values
behind it.

"The app's defaults" means the app, not a constant (issue #77)
--------------------------------------------------------------
This file used to answer "does the app open at the contract?" by comparing two constants —
``DEFAULT_PARAMS`` against ``PIPELINE_PARAMS["somatic"]`` — while no app module read the
first of them. Since ``DEFAULT_PARAMS`` was itself ``pipeline_params("somatic")``, that
comparison was ``deepcopy(X) == X``: unfailable, and blind to the only thing it claimed to
watch, because the app was free to seed a session from anything at all without moving
either side of it.

So every default asserted here is now read off something running: ``MAFigate.py`` booted,
for the seed a fresh session gets, and the parameter page rendered, for both arms and for
each control. Values the *contract* must hold are not asserted here at all — they belong to
``test_param_contract.py``, which checks them against ``nextflow.config``, which is where
they come from. The division is what stops this file drifting back into a second, quieter
statement of the same defaults.
"""

from __future__ import annotations

import copy
import os
import sys
from pathlib import Path
from typing import NamedTuple
from unittest import mock

import pytest

STREAMLIT_APP = Path(__file__).resolve().parent.parent
if str(STREAMLIT_APP) not in sys.path:
    sys.path.insert(0, str(STREAMLIT_APP))

from config.presets import (  # noqa: E402
    CLINICAL_GERMLINE_PARAMS,
    CLINICAL_SOMATIC_PARAMS,
    SOFT_GERMLINE_PARAMS,
    SOFT_SOMATIC_PARAMS,
)
from config.vocabularies import (  # noqa: E402
    CANCERVAR_OPTIONS,
    CIVIC_OPTIONS,
    CLINVAR_OPTIONS,
    CLINVAR_OTHER_ASSERTION_TERMS,
    CLINVAR_PATHOGENICITY_TERMS,
    ESCAT_OPTIONS,
    INTERVAR_OPTIONS,
    RENOVO_OPTIONS,
    VARIANT_CLASSIFICATIONS,
)
from config.pipeline_params import ARMS, PIPELINE_PARAMS, pipeline_params  # noqa: E402
from filters.variant_filters import (  # noqa: E402
    MAFIGATE_FILTER,
    MAFIGATE_REASON,
    PASS,
    REASON_RESCUE,
    apply_filters,
)
from vendor.pipeline_utils import has_clinvar_term, read_maf  # noqa: E402

FIXTURE_DIR = STREAMLIT_APP / "tests" / "fixtures" / "parity"

#: The catch-all this ticket deletes. Written once so a test cannot look for a different
#: spelling than the one the widgets used to offer.
SENTINEL = "All"

#: Every option list a filter multiselect is drawn from, by the parameter it feeds. The
#: classification list plus the six guideline sources — which is the whole of "that key
#: and all six guideline keys" in the acceptance criteria, enumerated so a seventh list
#: cannot be added without this file being edited.
OPTION_LISTS = {
    "filter_variant_classification": VARIANT_CLASSIFICATIONS,
    "filter_cancervar": CANCERVAR_OPTIONS,
    "filter_civic": CIVIC_OPTIONS,
    "filter_clinvar": CLINVAR_OPTIONS,
    "filter_intervar": INTERVAR_OPTIONS,
    "filter_renovo": RENOVO_OPTIONS,
    "filter_escat": ESCAT_OPTIONS,
}

#: The guideline sources each arm's vendored filter ORs together. Imported from the page
#: rather than restated: a third copy would be exactly the drift the house rule exists to
#: prevent, and a test that restates a relation cannot catch the page getting it wrong.
#: ``test_param_contract.test_the_pages_guideline_sources_follow_the_pipeline_filter_signatures``
#: is what checks *that* copy against the vendored signatures.
from config.vocabularies import GUIDELINE_SOURCES  # noqa: E402

#: One reference fixture per arm, and the row count it carries. The counts are asserted
#: rather than assumed so that "the report did not become the whole file" is measured
#: against a number this file knows, not against whatever the frame happens to be.
REFERENCE_FIXTURES = {
    "somatic": ("somatic_reference.maf", 82),
    "germline": ("germline_reference.maf", 94),
}


@pytest.fixture(scope="module")
def reference_mafs():
    """Both reference fixtures, loaded once by the pipeline's own reader."""
    return {
        arm: read_maf(str(FIXTURE_DIR / name))
        for arm, (name, _) in REFERENCE_FIXTURES.items()
    }


def _verdicts(maf, arm, **overrides):
    """``(verdict series, diagnostics)`` for one arm at the contract plus overrides."""
    params = pipeline_params(arm)
    params.update(overrides)
    labelled, diagnostics = apply_filters(maf, params)
    return labelled[MAFIGATE_FILTER], diagnostics


# ---------------------------------------------------------------------------
# The app opens at parity
# ---------------------------------------------------------------------------


def test_a_fresh_session_opens_at_the_contract():
    """The criterion "the pipeline parameter set becomes the app default", read off the app.

    This replaces ``test_the_apps_default_is_the_pipelines_own_parameter_set``, which
    asserted ``DEFAULT_PARAMS == PIPELINE_PARAMS["somatic"]`` over a constant *defined* as
    ``pipeline_params("somatic")`` — that is, ``deepcopy(X) == X``. It could not fail, and
    the failure it announced — "the app no longer opens at the pipeline's defaults" — was
    one it had no way to see: no app module read that constant, so the app could have been
    seeded from anything at all and the assertion would have gone on passing. Issue #77
    deleted the constant and moved the question to where it can be answered.

    The seed is read off a booted ``MAFigate.py``, which is the module that decides it.
    ``test_switching_arms_opens_that_arm_at_parity`` below covers the *parameter page*'s
    own initialisation, on both arms; neither it nor anything else touched app startup.

    The cache is stubbed for the reason ``_render_full_page`` gives, and here it is not
    merely isolation: ``MAFigate.py`` consults ``~/.mafigate`` *before* falling back to the
    contract, so the fallback is the branch under test and a developer's own cache would
    take the run down the other one.
    """
    pytest.importorskip("streamlit", reason="streamlit not installed")
    from streamlit.testing.v1 import AppTest

    from page_modules import param_store

    # Patched on ``param_store``, not on a ``MAFigate`` module object: the script runs its
    # own ``from page_modules.param_store import ...`` on every run, so it binds whatever
    # that module holds at the moment it executes. ``show_discarded_cache_banner`` is
    # stubbed whole because it is the path to ``discard_stale_cache``, which really does
    # move a file on the machine the suite runs on.
    with mock.patch.object(param_store, "load_parameters_from_cache", lambda: None), \
        mock.patch.object(param_store, "show_discarded_cache_banner", lambda: None), \
        mock.patch.dict(os.environ):
        # Otherwise the boot auto-loads a MAF and this stops being a test of the seed.
        os.environ.pop("MAFIGATE_OPEN_FILE", None)
        app = AppTest.from_file(str(STREAMLIT_APP / "MAFigate.py"), default_timeout=60)
        app.run()

    assert not app.exception, [str(e.value) for e in app.exception]
    assert app.session_state["filter_params"] == PIPELINE_PARAMS["somatic"], (
        "a fresh session does not open at the pipeline's defaults, so a user's first "
        "view is a cut nobody chose"
    )
    assert app.session_state["cache_loaded"] is False, (
        "the stubbed cache was reported as loaded, so this run did not take the branch "
        "under test"
    )

    # Second property, same boot, and deliberately not a second test sharing a fixture: the
    # check below *mutates* the dict this run produced, so a module-scoped booted app would
    # let this test corrupt the seed the assertion above reads — passing or failing on test
    # order. One boot, asked both questions in sequence, has no such ordering.
    #
    # The page mutates `filter_params` in place and the keep-lists are nested, so anything
    # shallower than a deep copy leaves a widget able to `append` to the contract's own
    # list. `test_param_contract.test_the_contract_is_not_mutable_shared_state` proves this
    # of `pipeline_params`; what is checked here is that the app's seed really is such a
    # copy. The test this replaces mutated `DEFAULT_PARAMS` itself — a module constant, in a
    # suite whose `AppTest` runs share one process — to make the same point.
    #
    # Compared against a snapshot rather than against written-out values: a literal `50` here
    # would be this file restating a contract value it has just said it holds none of, and
    # would be what the repair below wrote back if the contract's own default ever moved.
    contract_before = copy.deepcopy(PIPELINE_PARAMS["somatic"])
    app.session_state["filter_params"]["min_depth"] = 12345
    app.session_state["filter_params"]["filter_clinvar"].append("Benign")
    try:
        assert PIPELINE_PARAMS["somatic"] == contract_before, (
            "editing the session's parameters reached the contract itself, so a widget "
            "would redefine the app's default for every session after it"
        )
    finally:
        # A no-op on the happy path, by construction — and the guard says so rather than a
        # comment claiming it: only a real failure above leaves anything to put back, and a
        # polluted contract would make every later test in this process lie about parity.
        if PIPELINE_PARAMS["somatic"] != contract_before:
            PIPELINE_PARAMS["somatic"].clear()
            PIPELINE_PARAMS["somatic"].update(contract_before)


@pytest.mark.parametrize("arm", ARMS)
def test_switching_arms_opens_that_arm_at_parity(arm, reference_mafs):
    """"Without the user doing anything" has to hold on the arm they are actually on.

    A session opens somatic, so before this the germline arm inherited the somatic dict —
    which carries no ``filter_intervar`` and no ``filter_renovo``. Both arrived at
    the germline tab absent, became empty, and were written back: a criteria path of 13
    against the contract's 27 on ``germline_reference.maf``, and *silent*, because the
    all-empty warning cannot fire while the shared ``filter_clinvar`` is still populated.

    Asserted through the real page, because the bug lived in the page's arm handling and
    a unit test on the constants would have passed throughout.
    """
    app = _render_full_page(select_arm=arm)
    params = app.session_state["filter_params"]

    assert params["sample_type"] == arm
    for key in GUIDELINE_SOURCES[arm]:
        assert params[key] == PIPELINE_PARAMS[arm][key], (
            f"after switching to {arm}, {key} is {params.get(key)!r} rather than the "
            f"contract's {PIPELINE_PARAMS[arm][key]!r}"
        )

    # And it filters at parity, quoted on the criteria cell rather than the union.
    _, contract = _verdicts(reference_mafs[arm], arm)
    _, actual = apply_filters(reference_mafs[arm], dict(params))
    assert actual.cells() == contract.cells()


def test_the_pathogenic_retention_checkbox_actually_reaches_the_filter(reference_mafs):
    """It went dead the moment the contract became the default, and nothing said so.

    The checkbox wrote ``keep_pathogenic``; the filter prefers ``skip_pathogenic``
    whenever a dict carries it, and the contract carries it. So unticking "auto-retain
    pathogenic variants" changed nothing at all. Measured on ``somatic_reference.maf``:
    ticked and unticked both gave 18 rescued rows and a union of 20, where the real
    setting empties both rescue cells and moves the union to 19.

    Driven through the widget rather than asserted on the dict, because writing the wrong
    key is exactly the failure and only the widget can commit it.
    """
    app = _render("show_basic_filters_tab")
    checkbox = next(c for c in app.checkbox if "pathogenic" in c.label.lower())
    assert checkbox.value is True, "parity retains pathogenic variants"

    maf = reference_mafs["somatic"]
    _, retained = apply_filters(maf, dict(app.session_state["filter_params"]))

    checkbox.uncheck().run()
    assert not app.exception, [str(e.value) for e in app.exception]
    _, dropped = apply_filters(maf, dict(app.session_state["filter_params"]))

    assert dropped.cells() != retained.cells(), (
        "unticking auto-retain changed no cell of the decomposition — the control is "
        "writing a key the filter does not read"
    )
    assert dropped.rescue_only == 0 and dropped.cells()["both"] == 0, (
        f"pathogenic retention is off, so both rescue cells must be empty: "
        f"{dropped.cells()}"
    )


# `test_the_classification_default_is_the_configs_three_exclusions` stood here, asserting
# `DEFAULT_PARAMS["filter_variant_classification"] == ["Silent", "IGR", "RNA"]`. Issue #77
# deleted it: `test_param_contract.test_variant_classification_carries_the_configs_three_exclusions`
# asserts the same three values for *both* arms, and reads them out of `nextflow.config`
# itself rather than out of a copy the module under test made of the contract. That one takes
# the `config_params` fixture, so it skips where `nextflow.config` is absent; what holds the
# exclusions in a packaged tree is `test_filter_entry_point.py`, which filters at the contract
# against frozen pipeline verdicts and needs no `bin/`.


def test_filtering_at_the_apps_default_reaches_the_contracts_criteria_cell(
    reference_mafs,
):
    """The default reaches the filter, asserted on the cell the change actually moves.

    The union PASS total is the wrong number to check here and checking it would have
    passed vacuously: pathogenic retention re-admits the rows an inverted classification
    list drops, so the two settings union to the same total. Issue #28 names this trap —
    "counts that coincide" — and it is why the assertion quotes the criteria path instead.
    On this fixture set the contract's criteria path is 54 of 82 rows against a union of
    60, so the six rows the rescue contributes are the ones the total would have hidden.

    This used to filter twice and compare the two runs — once through ``pipeline_params``
    and once through ``dict(DEFAULT_PARAMS)``, which was a copy of the same dict. Both
    halves were the contract, so the comparison was between a value and itself; issue #77
    dropped it and kept the measurement, which is the part that was ever load-bearing.
    """
    maf = reference_mafs["somatic"]
    _, from_default = _verdicts(maf, "somatic")

    assert from_default.criteria_path == 54


def test_a_deviating_preset_still_deviates_after_being_converted(reference_mafs):
    """Converting SOFT's classification list must not have converted SOFT into parity.

    The complement was taken within the control's own vocabulary, so a preset keeps its
    intent — but a preset that came out *identical* to the contract would mean the
    conversion had quietly collapsed the deviation, and the presets exist to be
    deviations.
    """
    maf = reference_mafs["somatic"]
    _, at_parity = _verdicts(maf, "somatic")
    _, under_soft = apply_filters(maf, dict(SOFT_SOMATIC_PARAMS))

    assert under_soft.cells() != at_parity.cells(), (
        "SOFT_SOMATIC_PARAMS now decides every cell exactly as the contract does; it has "
        "stopped being a deviation"
    )


# ---------------------------------------------------------------------------
# The sentinel is gone
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("param,options", sorted(OPTION_LISTS.items()))
def test_no_option_list_offers_the_catch_all_sentinel(param, options):
    """It cannot be selected, because it is not on offer.

    The weakest of the three unreachability arguments and deliberately first: it is the
    one a reader can check by eye, and the two behavioural tests below are what make it
    more than a spelling change.
    """
    assert SENTINEL not in options, (
        f"{param}'s option list still offers {SENTINEL!r}. Selecting it reads as "
        "'no restriction' but saturates the guideline union — 68,364 somatic rows "
        "against a parity baseline of 411 (issue #28)."
    )
    assert options, f"{param}'s option list is empty"


@pytest.mark.parametrize("arm", ARMS)
def test_the_sentinel_value_decides_every_row_exactly_as_an_empty_list_does(
    arm, reference_mafs
):
    """The 167-fold widening, asserted unreachable *by value* rather than by option list.

    A stale ``~/.mafigate`` cache or an old exported file can still carry ``["All"]``
    long after the widgets stop offering it, so removing the option is not on its own a
    proof. What makes it safe is that nothing branches on the value any more: issue #33
    deleted the guard that skipped a whole guideline clause when its source read ``All``,
    and the vendored code has no such concept — ``isin(["All"])`` is simply a test no row
    passes.

    Compared per row rather than by total, because equal counts either side of a
    substitution is the specific coincidence this codebase keeps producing.
    """
    maf = reference_mafs[arm]
    sources = GUIDELINE_SOURCES[arm]

    empty, _ = _verdicts(maf, arm, **{key: [] for key in sources})
    sentinel, _ = _verdicts(maf, arm, **{key: [SENTINEL] for key in sources})
    assert sentinel.equals(empty), (
        f"{arm}: setting every guideline source to {SENTINEL!r} decides some row "
        "differently from setting them all empty, so the sentinel is still live"
    )

    # And one source at a time, which is how a user reaches it: the saturation only ever
    # needed *one* source to read as unrestricted.
    for key in sources:
        one_empty, _ = _verdicts(maf, arm, **{key: []})
        one_sentinel, _ = _verdicts(maf, arm, **{key: [SENTINEL]})
        assert one_sentinel.equals(one_empty), f"{arm}: {key} still honours the sentinel"


@pytest.mark.parametrize("arm", ARMS)
def test_no_guideline_selection_can_widen_the_report_past_the_contracts(
    arm, reference_mafs
):
    """Emptying a source is neutral or narrowing — never saturating.

    The bug's signature was a *widening*: one source reading as unrestricted made the OR
    true for every row, so the report became the file. This asserts the property that
    signature violates, per source and for all of them at once, as a subset relation on
    the passing rows rather than as a count.
    """
    maf = reference_mafs[arm]
    rows = REFERENCE_FIXTURES[arm][1]
    assert len(maf) == rows, "fixture changed shape; the bound below is no longer known"

    contract, contract_diagnostics = _verdicts(maf, arm)
    passing = set(contract.index[contract == PASS])
    assert passing, f"{arm}: the contract passes no row, so there is nothing to bound"

    for keys in [(key,) for key in GUIDELINE_SOURCES[arm]] + [GUIDELINE_SOURCES[arm]]:
        widened, diagnostics = _verdicts(maf, arm, **{key: [] for key in keys})
        emptied = set(widened.index[widened == PASS])
        assert emptied <= passing, (
            f"{arm}: emptying {list(keys)} admitted rows the contract rejects "
            f"({sorted(emptied - passing)[:5]}) — an empty keep-list is saturating a "
            "guideline source instead of dropping it out of the union"
        )
        assert diagnostics.passed < rows, (
            f"{arm}: emptying {list(keys)} passed all {rows} rows. This is the shape of "
            "the 167-fold widening, at fixture scale."
        )

    # The contract itself is nowhere near the file, which is what makes the bound above
    # worth having rather than vacuous.
    assert contract_diagnostics.passed < rows


@pytest.mark.parametrize("arm", ARMS)
def test_the_all_empty_state_is_expressible_and_stops_at_the_rescue_floor(
    arm, reference_mafs
):
    """Warned at the widget, not blocked — so the filter must survive it.

    Emptying every guideline source is expressible on the pipeline's own command line,
    so refusing it in the app would be the app being stricter than the pipeline. What it
    produces is the pathogenic-rescue floor: the guideline block admits nothing, the
    criteria path is empty, and every surviving row is there by rescue alone.
    """
    maf = reference_mafs[arm]
    verdicts, diagnostics = _verdicts(
        maf, arm, **{key: [] for key in GUIDELINE_SOURCES[arm]}
    )

    assert diagnostics.criteria_path == 0, (
        f"{arm}: with no guideline source keeping anything, no row can meet the "
        f"criteria — but {diagnostics.criteria_path} did"
    )
    assert diagnostics.rescue_only == diagnostics.passed
    assert 0 < diagnostics.passed < len(maf), (
        f"{arm}: the floor is neither zero nor the whole file; got "
        f"{diagnostics.passed} of {len(maf)}"
    )
    assert diagnostics.cells()["criteria_only"] == 0
    assert diagnostics.cells()["both"] == 0


@pytest.mark.parametrize("arm", ARMS)
def test_every_row_at_the_floor_says_it_was_rescued(arm, reference_mafs):
    """The reason column agrees, so a clinician can see *why* the report collapsed."""
    maf = reference_mafs[arm]
    params = pipeline_params(arm)
    params.update({key: [] for key in GUIDELINE_SOURCES[arm]})
    labelled, _ = apply_filters(maf, params)

    survivors = labelled[labelled[MAFIGATE_FILTER] == PASS]
    assert not survivors.empty
    assert set(survivors[MAFIGATE_REASON]) == {REASON_RESCUE}


# ---------------------------------------------------------------------------
# The presets, converted rather than left reading backwards
# ---------------------------------------------------------------------------

PRESETS = {
    "SOFT_SOMATIC_PARAMS": SOFT_SOMATIC_PARAMS,
    "SOFT_GERMLINE_PARAMS": SOFT_GERMLINE_PARAMS,
    "CLINICAL_SOMATIC_PARAMS": CLINICAL_SOMATIC_PARAMS,
    "CLINICAL_GERMLINE_PARAMS": CLINICAL_GERMLINE_PARAMS,
}


@pytest.mark.parametrize("name,preset", sorted(PRESETS.items()))
def test_a_presets_classification_list_reads_as_an_exclude_list(name, preset):
    """SOFT and CLINICAL survive as deviations — but not as *inverted* ones.

    Every preset used to hold an include list headed by ``Missense_Mutation``. Read as
    the exclude list the parameter now is, that says "drop every missense call", which no
    preset means. The two are told apart by exactly this: an exclude list never names the
    classifications the preset exists to report.
    """
    excluded = preset["filter_variant_classification"]
    assert isinstance(excluded, list)
    for kept in ("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del"):
        assert kept not in excluded, (
            f"{name} excludes {kept}, so its classification list is still the old "
            "include list being read backwards"
        )
    assert set(excluded) <= set(VARIANT_CLASSIFICATIONS), (
        f"{name} excludes a classification the control cannot offer: "
        f"{sorted(set(excluded) - set(VARIANT_CLASSIFICATIONS))}"
    )


def test_clinical_excludes_everything_soft_does_and_more():
    """The one relation between the two converted lists, guarded rather than assumed.

    CLINICAL is the stricter preset — it reported protein-coding and canonical
    splice-site calls only, where SOFT also reported regulatory and intronic ones. Under
    the exclude polarity that inverts into a containment: CLINICAL must exclude
    everything SOFT excludes, and strictly more. The two lists were hand-computed
    complements, so this is what stops one being edited into disagreeing with the other.

    Note what deliberately has *no* guard: neither list is derived from
    ``VARIANT_CLASSIFICATIONS``, so adding a classification to the vocabulary leaves both
    presets keeping it. That is correct under an exclude list and is the whole point of
    the polarity — a value nobody has named is kept, as the pipeline keeps it.
    """
    soft = set(SOFT_SOMATIC_PARAMS["filter_variant_classification"])
    clinical = set(CLINICAL_SOMATIC_PARAMS["filter_variant_classification"])

    assert soft < clinical, (
        f"CLINICAL must exclude a strict superset of SOFT's exclusions; SOFT excludes "
        f"{sorted(soft)} and CLINICAL excludes {sorted(clinical)}"
    )
    for arm_preset, name in (
        (SOFT_GERMLINE_PARAMS, "SOFT_GERMLINE_PARAMS"),
        (CLINICAL_GERMLINE_PARAMS, "CLINICAL_GERMLINE_PARAMS"),
    ):
        counterpart = soft if "SOFT" in name else clinical
        assert set(arm_preset["filter_variant_classification"]) == counterpart, (
            f"{name} excludes different classifications from its somatic twin; the two "
            "arms filter classification identically in the pipeline"
        )


@pytest.mark.parametrize("name,preset", sorted(PRESETS.items()))
def test_a_preset_is_still_a_deliberate_deviation_from_parity(name, preset):
    """They are meant to differ from the contract, and to stay selectable."""
    assert preset["filter_variant_classification"] != ["Silent", "IGR", "RNA"] or (
        preset["max_freq_population"] < 1.0
    ), f"{name} has quietly become the parity preset"
    assert preset["sample_type"] in ARMS


# ---------------------------------------------------------------------------
# The dead CIViC rating options
# ---------------------------------------------------------------------------


def test_the_dead_civic_rating_options_are_gone():
    """Imported but never used, and no parameter behind it.

    ``CIVIC_RATING_OPTIONS`` offered evidence *ratings* 1–5, which nothing in the app or
    the pipeline reads: the pipeline's only CIViC parameter is
    ``filter_civic_evidence_level``. A control that cannot affect the report is the same
    lie about what a control does that the rest of this ticket is fixing.
    """
    import config.vocabularies as vocabularies

    assert not hasattr(vocabularies, "CIVIC_RATING_OPTIONS")

    for module in ("page_modules/data_loading.py", "page_modules/parameter_config.py"):
        source = (STREAMLIT_APP / module).read_text()
        assert "CIVIC_RATING_OPTIONS" not in source, f"{module} still imports it"


# ---------------------------------------------------------------------------
# The ClinVar keep-terms that could not keep anything (issue #88)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("term", CLINVAR_OPTIONS)
def test_every_offered_clinvar_term_can_keep_a_variant(term):
    """The offered term, asked of the matcher that actually decides — not of a regex.

    ``CLINVAR_OPTIONS`` carried ``Pathogenic/Likely_pathogenic`` and
    ``Benign/Likely_benign`` for the app's whole life. Neither could ever keep a row:
    ``has_clinvar_term`` splits the *annotation* on ``[|/;,]`` and asks whether a piece is
    in the keep-list, so a keep-term holding a separator cannot equal any piece of
    anything. The control offered two doors that opened onto nothing.

    The property checked is the weakest one that catches it and the strongest one that is
    true of every honest term: **a variant annotated with exactly this term is kept by
    exactly this term.** A term failing that is unusable by construction, whatever the
    reason — a separator today, some future normalisation tomorrow.

    Written against the vendored matcher on purpose. A test spelling the ``[|/;,]`` rule
    out again would pass while the pipeline's own splitting moved underneath it, which is
    the failure this map keeps finding: the guard restates the claim instead of running it.
    """
    assert has_clinvar_term(term, [term]), (
        f"{term!r} is offered by the ClinVar control but cannot keep even a variant "
        f"annotated exactly {term!r} — the filter splits the annotation on [|/;,] and "
        "matches the pieces, so a keep-term containing a separator matches nothing"
    )


@pytest.mark.parametrize("name,preset", sorted(PRESETS.items()))
def test_a_presets_clinvar_terms_are_all_offered_and_all_live(name, preset):
    """A preset may not ship a term the control cannot offer, or one that keeps nothing.

    Both Broad presets shipped ``Pathogenic/Likely_pathogenic`` until #88 — dead weight in
    a list a clinician reads as *what will be kept*. Dropping it changed no row, because
    both presets also carry the atoms it splits into; that it changed no row is precisely
    why it survived so long, and precisely why the preset was lying for free.

    The subset half is the one that stops the next drift: a term removed from the
    vocabulary but left in a preset is invisible in the UI — ``filter_terms`` silently
    drops it from the widget default — while still sitting in the file a reader trusts.
    """
    terms = preset["filter_clinvar"]
    assert terms, f"{name} keeps no ClinVar term at all"

    unofferable = [t for t in terms if t not in CLINVAR_OPTIONS]
    assert not unofferable, (
        f"{name} carries ClinVar term(s) {unofferable} that the control does not offer; "
        "the widget drops them from its default, so the preset says one thing and the "
        "screen shows another"
    )

    dead = [t for t in terms if not has_clinvar_term(t, [t])]
    assert not dead, (
        f"{name} carries ClinVar term(s) {dead} that cannot keep any variant; a preset "
        "names what it keeps, and these keep nothing at any value on any MAF"
    )


#: Offered ClinVar terms the Clinical Summary cannot name. **Empty since issue #98.**
#:
#: It held ``other`` — the term #88 declined to classify, on the ground that every bucket
#: available was a clinical claim this repository had no authority to make. That ground was
#: the part that turned out to be wrong: an authority does exist. ClinVar's own
#: documentation groups its terms, and the institute keeps a term table that classes each
#: one and states what to do with it, which #98 read as the source for the mapping.
#:
#: ``other`` is now ``No_Classification`` — not *mapped to* ``not_provided``, which is what
#: #88 refused and was right to refuse, but sharing a class named for the class both terms
#: belong to. ClinVar defines ``other`` as "ClinVar does not have the appropriate term for
#: your submission", and the table excludes it because its meaning depends on a free-text
#: explanation the app never receives; both say there is no usable classification, which is
#: what the class asserts and is not a severity.
#:
#: ``_other`` — the spelling #99 put on offer beside it — takes the same class, for #99's own
#: reason: the two are one ClinVar call written two ways, so whatever decides one decides
#: both, and a class that depended on how a file was written is the shape #99 removed.
#: Neither rescues a row by itself, the summary reading the call rather than the modifier.
#:
#: Kept rather than deleted, empty, because it is the device that catches the *next* term:
#: a term added to the control and not to the mapping fails below, and the only way to make
#: that failure go away is to write the reason down here.
_UNMAPPED_BY_DESIGN = set()


@pytest.mark.parametrize("term", [t for t in CLINVAR_OPTIONS if t not in _UNMAPPED_BY_DESIGN])
def test_every_offered_clinvar_term_is_readable_by_the_clinical_summary(term):
    """Keeping a variant is half the job; the report then has to say what it is.

    A ClinVar term reaches two places, and adding one only to the vocabulary leaves the
    second silently wrong. ``CLINICAL_VALUE_MAPPING`` translates an annotation value into a
    class, and a value it does not hold cannot be named — so the report keeps the variant
    and then cannot say what it was kept for.

    That is not hypothetical. #88 put ``Conflicting_interpretations_of_pathogenicity`` on
    offer and into both Broad presets, making 77 rows of the measured data newly reportable
    — and every one of them arrived in the report labelled as having no clinical data,
    because the mapping knew only the post-2023 spelling. Review found it; this is what
    stops the next term doing the same.
    """
    pytest.importorskip("streamlit", reason="streamlit not installed")

    assert term not in _terms_the_summary_cannot_name(), (
        f"a variant annotated {term!r} — a term the ClinVar control offers and a preset "
        "can keep — reaches the report with no class of its own. The value is missing from "
        "CLINICAL_VALUE_MAPPING; map it against its class in the institute's term table"
    )


def _terms_the_summary_cannot_name():
    """Offered ClinVar terms that reach the summary's *unrecognised* label.

    Read the signal carefully. Until issue #98 an unreadable value was **discarded**, so a
    term the mapping did not hold surfaced as ``No Clinical Data`` and that string was how
    these guards detected it. It is no longer: an unreadable value is now kept and rendered
    as its own class, and ``No Clinical Data`` is reserved for a variant no source annotated
    at all. Testing for the old string here would have made both guards below **vacuous** —
    passing for every term, including one genuinely missing from the mapping.
    """
    import pandas as pd

    from components.clinical_summary import _SUMMARY_LABELS, NO_CLINICAL_DATA
    from components.clinical_summary import generate_clinical_summary

    unnameable = set()
    for term in CLINVAR_OPTIONS:
        summary = generate_clinical_summary(pd.Series({"ClinVar_VCF_CLNSIG": term}))
        if summary in (_SUMMARY_LABELS["Unknown"], NO_CLINICAL_DATA):
            unnameable.add(term)
    return unnameable


def test_the_unmapped_clinvar_terms_are_exactly_the_documented_ones():
    """The exception list is pinned, so it cannot quietly become the rule.

    An exception carried in a set beside the test that honours it is one edit away from
    being where unmapped terms go to stop failing. This asserts the set is exactly what is
    unmapped: a term added to it that the summary *can* read fails here, and so does a
    newly-unreadable term that someone added to the exception instead of to the mapping.

    The set is empty since #98, which is the strongest state it has been in — every term the
    control offers is named by a class — and the test is kept for the term after next.
    """
    pytest.importorskip("streamlit", reason="streamlit not installed")

    unnameable = _terms_the_summary_cannot_name()

    assert unnameable == _UNMAPPED_BY_DESIGN, (
        f"the Clinical Summary cannot name {sorted(unnameable)}, but the documented "
        f"exception is {sorted(_UNMAPPED_BY_DESIGN)}. Map the term in "
        "CLINICAL_VALUE_MAPPING against its class in the institute's term table, or — if no "
        "class is honest for it — widen the exception here with the reason written down"
    )


# ---------------------------------------------------------------------------
# The control offers everything the report can name (issue #103)
# ---------------------------------------------------------------------------
#
# The other direction of the contract above. #88 caught a term the control offered and the
# report could not name; #103 caught the mirror image — nineteen terms the report named and
# the control would not offer, so 1,215 rows of a 159,580-row real corpus could be kept by no
# selection available anywhere in the app. One defect with two faces, and it kept recurring
# because only one face had a guard. Both do now.


#: The other five sources' vocabularies, so ClinVar's own mapping keys can be told apart.
#:
#: ``CLINICAL_VALUE_MAPPING`` is one flat dict over several annotation sources and carries no
#: source tag, so *which of its keys are ClinVar's* cannot be read off it. What can be read
#: off it is which keys another source claims — every other source's vocabulary is a constant
#: in ``config/vocabularies.py`` — leaving ClinVar's as the remainder. Derived rather than
#: listed, so a term added to any source's vocabulary moves this set with nobody remembering.
#:
#: It said "six annotation sources", and that was **false when it was written**: ESCAT has never
#: been in that dict, its levels being mapped by group from its own table, so the dict spanned
#: five. CIViC left it at issue #109 for the same reason — its levels grade evidence rather than
#: the variant — which leaves four and is why the number is gone rather than decremented.
#:
#: ``CIVIC_OPTIONS`` stays in the subtraction deliberately, though nothing it names is a key
#: any more: a bare ``A`` arriving back in the mapping must not be counted as a ClinVar term the
#: control fails to offer. Subtracting a vocabulary the dict does not hold costs nothing;
#: forgetting one is how this derivation would start making claims about ClinVar out of another
#: source's values.
_OTHER_SOURCE_VOCABULARIES = (
    set(CANCERVAR_OPTIONS)
    | set(INTERVAR_OPTIONS)
    | set(RENOVO_OPTIONS)
    | set(CIVIC_OPTIONS)
)

#: ClinVar terms the report can name that the control may **not** offer, and why.
#:
#: Every one of them carries a separator in its own name, so the split that matches a keep-term
#: against the pieces of an annotation can never hand it back whole: offering one would be
#: offering a door onto nothing, which is exactly what #88 deleted. They stay in the *mapping*
#: because a real cell spells them that way and the report has to name what it renders — the
#: asymmetry is the point, and it is why this list exists rather than the two constants simply
#: agreeing.
#:
#: Pinned rather than computed so that a fifth such term cannot join the mapping and slip
#: silently past the offering rule; the test below re-derives the property and compares.
_NOT_OFFERABLE_BY_CONSTRUCTION = {
    "Pathogenic/Likely_pathogenic",
    "Benign/Likely_benign",
    "Pathogenic,_low_penetrance",
    "Likely_pathogenic,_low_penetrance",
}


def test_every_clinvar_term_the_report_can_name_is_offered():
    """The reverse of the mapping contract, and the defect issue #103 was opened for.

    A term the report can label is a term a user should be able to ask for. Before #103 the
    app labelled ``⚠️ Disease Risk`` on a variant and offered no way to select it — the
    control's list had stopped at ten classifications while ``ClinVar_VCF_CLNSIG`` carried
    twenty-odd, so the vocabulary was a transcription that stopped early rather than a
    shortlist anyone defended. Measured through the vendored matcher over 176 byte-distinct
    real MAFs: **1,215 of 159,580 ClinVar-annotated rows could be kept by no selection the
    control offered**, 854 of them carrying ``no_classification_for_the_single_variant``
    alone, in 102 of those files, while the app offered ``not_provided`` — the same class.

    The exception is structural and it is checked, not trusted: a term may go unoffered only
    if the matcher cannot hand it back, which is asserted here rather than asserted of a list
    of names. So a ClinVar term added to the mapping either becomes selectable or fails.
    """
    pytest.importorskip("streamlit", reason="streamlit not installed")

    from components.clinical_summary import CLINICAL_VALUE_MAPPING

    clinvar_only = set(CLINICAL_VALUE_MAPPING) - _OTHER_SOURCE_VOCABULARIES

    # Run the matcher rather than restate the [|/;,] rule -- #99's line, and the reason this
    # exception cannot rot as the pipeline's splitting moves.
    unofferable = {term for term in clinvar_only if not has_clinvar_term(term, [term])}
    assert unofferable == _NOT_OFFERABLE_BY_CONSTRUCTION, (
        "the terms the matcher cannot hand back are "
        f"{sorted(unofferable)}, but the documented exception is "
        f"{sorted(_NOT_OFFERABLE_BY_CONSTRUCTION)}. A new one means a term was added to "
        "CLINICAL_VALUE_MAPPING whose own name carries a separator; a missing one means a "
        "term became matchable and should now be offered"
    )

    missing = sorted(clinvar_only - set(CLINVAR_OPTIONS) - unofferable)
    assert not missing, (
        f"the Clinical Summary can name ClinVar term(s) {missing}, and the ClinVar control "
        "does not offer them — so the app labels a variant with a classification the user "
        "cannot filter for. Add them to CLINVAR_PATHOGENICITY_TERMS or "
        "CLINVAR_OTHER_ASSERTION_TERMS according to the class the institute's term table "
        "gives them, or, if the matcher cannot hand the term back, to "
        "_NOT_OFFERABLE_BY_CONSTRUCTION with that reason"
    )

    # The derivation above is blind to the two names ClinVar shares with InterVar, which is
    # harmless only because they are the two least likely terms in the app to go missing --
    # so it is said out loud rather than left as a hole in a subtraction.
    shared = {"Pathogenic", "Benign"}
    assert shared <= set(CLINVAR_OPTIONS), (
        f"{sorted(shared - set(CLINVAR_OPTIONS))} is spelled identically by InterVar, so the "
        "ClinVar-only derivation above cannot see it; the ClinVar control must offer it"
    )


def test_the_two_clinvar_controls_partition_the_vocabulary():
    """Two widgets, one vocabulary — with nothing lost between them and nothing doubled.

    The split is a fact about the screen (issue #103): thirty-one terms in one multiselect
    buried ``Pathogenic``, so the classifications and the other assertions are drawn as two
    controls. Every other reader — validation, the Help page's walk, ``param_migration``, the
    parity harness — sees ``CLINVAR_OPTIONS`` and must not know or care.

    A term in neither list is offered nowhere while every list-shaped guard still passes; a
    term in both is drawn twice and, because the page accumulates the two selections, lands
    in the parameter file twice.
    """
    assert (
        list(CLINVAR_PATHOGENICITY_TERMS) + list(CLINVAR_OTHER_ASSERTION_TERMS)
        == list(CLINVAR_OPTIONS)
    ), (
        "CLINVAR_OPTIONS is no longer the two control lists in order, so a reader of the "
        "whole vocabulary and a user reading the two controls see different things"
    )

    overlap = sorted(set(CLINVAR_PATHOGENICITY_TERMS) & set(CLINVAR_OTHER_ASSERTION_TERMS))
    assert not overlap, f"term(s) {overlap} are drawn by both ClinVar controls"

    duplicated = sorted({t for t in CLINVAR_OPTIONS if CLINVAR_OPTIONS.count(t) > 1})
    assert not duplicated, f"the ClinVar vocabulary repeats {duplicated}"


def test_the_penetrance_piece_is_deliberately_not_offered():
    """``_low_penetrance`` is the one place #99's every-spelling rule is not applied.

    It is not a classification. It is the tail of ``Pathogenic,_low_penetrance``, left behind
    by a split that cuts on the comma in that term's own name — so there is no honest class
    for it and ``CLINICAL_VALUE_MAPPING`` does not name it. Refusing it costs nothing, and
    that is measured rather than assumed: all 8 of its occurrences in the real corpus sit
    beside a real call, so it never stands as a whole cell and no row anywhere is reachable
    only through it.

    The pipeline's own ``--filter_clinvar`` help does advertise it. That help is the
    vocabulary #103 declined to mirror — it also advertises an empty term — so its offering
    this is not an argument, and this test records that judgement where it can fail if
    someone adds the term without revisiting it.
    """
    pytest.importorskip("streamlit", reason="streamlit not installed")

    from components.clinical_summary import CLINICAL_VALUE_MAPPING

    assert "_low_penetrance" not in CLINVAR_OPTIONS, (
        "`_low_penetrance` is offered. It is a penetrance modifier rather than a "
        "classification, and the report has no class for it -- so offering it would keep "
        "variants the summary then labels as an unrecognised annotation"
    )
    assert "_low_penetrance" not in CLINICAL_VALUE_MAPPING, (
        "`_low_penetrance` has been given a class. If that is deliberate then it is now "
        "nameable, and test_every_clinvar_term_the_report_can_name_is_offered requires it "
        "to be offered too -- decide both together, and rewrite this test's reason"
    )


def test_both_broad_presets_keep_the_no_classification_class_whole():
    """A preset that keeps a class keeps the class (issue #103).

    Both Broad presets already kept ``not_provided``, ``other`` and ``_other`` — three
    members of the class #98 named **No Classification**. They omitted the rest only because
    the control did not offer it, which is not a decision. So #103 completed the class rather
    than deciding a fresh one, at a measured +855 rows on the ClinVar clause over 176 real
    MAFs, 854 of them one term.

    Asserted as *the class*, read off the mapping, rather than as a list of four names: a term
    joining the class later is caught, and the test cannot drift from the reason it exists.
    Stringent is deliberately excluded — it keeps two pathogenicity calls and this class is
    not among them.
    """
    pytest.importorskip("streamlit", reason="streamlit not installed")

    from components.clinical_summary import CLINICAL_VALUE_MAPPING

    no_classification = {
        term
        for term in CLINVAR_OPTIONS
        if CLINICAL_VALUE_MAPPING.get(term) == "No_Classification"
    }
    assert no_classification, (
        "no offered ClinVar term is classed No_Classification any more; this test's subject "
        "has gone and its reasoning needs rereading, not its assertion relaxing"
    )

    for name in ("SOFT_SOMATIC_PARAMS", "SOFT_GERMLINE_PARAMS"):
        kept = set(PRESETS[name]["filter_clinvar"])
        assert no_classification <= kept, (
            f"{name} keeps some of the No Classification class and not all of it — missing "
            f"{sorted(no_classification - kept)}. That is the incoherence #103 found in the "
            "vocabulary (not_provided offered, no_classification_for_the_single_variant not) "
            "moved into the preset. Keep the class whole, or record here why this member is "
            "different"
        )


def test_a_broad_preset_keeps_each_clinvar_class_whole_or_not_at_all():
    """Half a class is the shape of every ClinVar defect this map has found.

    #88: a preset kept a call under one spelling. #99: a call was offered under one spelling.
    #103: the vocabulary held three members of a class and not the fourth, and the app offered
    ``not_provided`` while refusing ``no_classification_for_the_single_variant`` — 854 rows of
    the same class. Each time the report a clinician gets depended on which member of one
    class their file happened to spell, and each time nothing failed.

    The invariant that would have caught all three: **a preset keeps a ClinVar class wholly or
    not at all.** Read off ``CLINICAL_VALUE_MAPPING`` rather than off a list of names, so a
    term joining a class later is covered without anyone remembering this test.

    **There is no exception any more, and the set that held one is deleted rather than left
    standing empty.** ``_PARTIALLY_KEPT_BY_DESIGN`` named ``Drug_Response``: both Broad presets
    kept ``drug_response`` without ``confers_sensitivity``, because that term is the one
    placement #98 made with no citable definition and #103 declined to guess at it. #114 keeps
    the class whole, so the exception is gone — and gone as a *symbol*, since an empty set
    beside its own test is one edit from being where inconsistencies go to stop failing, which
    is what its own note said about it. Adding an exception now means arguing for it in this
    docstring.

    Stringent is excluded by design — it keeps two pathogenicity calls as its whole point, and
    ``Likely_Pathogenic`` is a class of one, so it satisfies this trivially and asserting it
    would suggest the test had checked something.
    """
    pytest.importorskip("streamlit", reason="streamlit not installed")

    from components.clinical_summary import CLINICAL_VALUE_MAPPING

    members = {}
    for term in CLINVAR_OPTIONS:
        members.setdefault(CLINICAL_VALUE_MAPPING[term], set()).add(term)

    for name in ("SOFT_SOMATIC_PARAMS", "SOFT_GERMLINE_PARAMS"):
        kept = set(PRESETS[name]["filter_clinvar"])
        partial = {
            klass: sorted(terms - kept)
            for klass, terms in members.items()
            if (terms & kept) and (terms - kept)
        }
        assert not partial, (
            f"{name} keeps these ClinVar classes only in part: "
            f"{ {k: v for k, v in sorted(partial.items())} }. A class kept in part means "
            "whether a variant reaches the report depends on which member of one class its "
            "file spells — the shape of #88, #99 and #103. Keep the class whole, or add the "
            "exception to this test's docstring with the reason it is one"
        )


#: The ClinVar classes a Broad preset holds out, and the whole of the rule (issue #114).
#:
#: **These are class names, not term names**, and the two are one capital letter apart here —
#: the classes are ``Likely_Benign``/``Benign`` and the ClinVar terms are
#: ``Likely_benign``/``Benign``. Spelled as classes because the class is the unit a preset keeps
#: (``test_a_broad_preset_keeps_each_clinvar_class_whole_or_not_at_all``), so a new benign *term*
#: — some future ``Benign_low_penetrance`` mapping into one of these classes — is held out
#: without anyone editing this.
#:
#: **A new benign *class* is not**, and that limit is worth stating rather than discovering: the
#: rule below derives what Broad must keep by subtracting this set from the classes actually in
#: use, so a class this does not name is a class Broad is *required* to keep. If ClinVar's
#: vocabulary ever gains a benign class of its own, the test fails and this set is where the
#: answer goes — which is the intended failure, since whether Broad keeps a new class is exactly
#: the judgement #114 had to be told.
_CLASSES_BROAD_HOLDS_OUT = {"Benign", "Likely_Benign"}


def test_broad_keeps_every_clinvar_class_except_the_benign_ones():
    """What Broad keeps is one line, and #114 is where it became one.

    Until then the keep-list was a set of classes chosen one ticket at a time, and the test
    above could only ask that each be kept *whole* — a class dropped entirely satisfied it, so
    nothing recorded which classes Broad had actually decided on. #103 left four selectable and
    unselected precisely because the decision had not been made; #114 made it, and the answer
    turned out to be a rule rather than four calls: **Broad keeps every ClinVar class except
    🟢 Likely Benign and 🔵 Benign.**

    That is worth pinning in both directions, and the second half is the one that carries the
    meaning. Dropping a class would leave the whole-or-not-at-all test green while narrowing a
    preset described to the user as *a wide net*; adding a benign class would make Broad keep
    everything ClinVar offers, at which point the ClinVar clause stops filtering anything and
    the preset silently becomes a different cut.

    Read off ``CLINICAL_VALUE_MAPPING``, never a list of term names, for the reason #103 and
    #98 both give: the mapping is the record of what the report can *name*, and a guard written
    against names goes stale the moment ClinVar adds a term — which is how the app came to
    offer ``not_provided`` while refusing ``no_classification_for_the_single_variant``, 854 rows
    of one class, for as long as it did.

    Stringent is excluded for the same reason as above: two pathogenicity calls are its whole
    point, and #103 and #114 both left it untouched — nothing yet argues a report you would act
    on should carry a non-pathogenicity assertion.

    **This makes ``test_both_broad_presets_keep_the_no_classification_class_whole`` redundant,
    and it stays anyway** — deliberately, unlike ``_PARTIALLY_KEPT_BY_DESIGN``, which this change
    deleted for being redundant. The difference is what each one holds. That set held an
    *exception*, and an unused exception is a place for the next inconsistency to be filed
    quietly; #103's test holds a *finding* — that the app offered ``not_provided`` while refusing
    ``no_classification_for_the_single_variant``, 854 rows of one class across 102 files — and
    deleting it would delete the measurement that argued for the rule this test now states in
    general. A redundant assertion costs a few milliseconds; a lost reason costs the next ticket.
    """
    pytest.importorskip("streamlit", reason="streamlit not installed")

    from components.clinical_summary import CLINICAL_VALUE_MAPPING

    classes = {CLINICAL_VALUE_MAPPING[term] for term in CLINVAR_OPTIONS}
    unknown = _CLASSES_BROAD_HOLDS_OUT - classes
    assert not unknown, (
        f"{sorted(unknown)} is held out by name but is not a class of any offered ClinVar "
        "term, so this test is guarding a class that no longer exists — re-read the rule "
        "against CLINICAL_VALUE_MAPPING rather than adjusting this constant"
    )
    expected = classes - _CLASSES_BROAD_HOLDS_OUT

    for name in ("SOFT_SOMATIC_PARAMS", "SOFT_GERMLINE_PARAMS"):
        kept = set(PRESETS[name]["filter_clinvar"])
        actual = {CLINICAL_VALUE_MAPPING[term] for term in kept}

        missing = sorted(expected - actual)
        assert not missing, (
            f"{name} no longer keeps the ClinVar class(es) {missing}. Broad keeps everything "
            "ClinVar has not called benign (issue #114) and its own description promises a "
            "wide net, so narrowing it is a decision to argue for beside "
            "`_SOFT_CLINVAR_TERMS`, not a list to edit"
        )

        benign = sorted(actual & _CLASSES_BROAD_HOLDS_OUT)
        assert not benign, (
            f"{name} now keeps the ClinVar class(es) {benign}. Benign is the one line Broad "
            "draws — keeping it leaves the ClinVar clause admitting every term the control "
            "offers, which is not a wider net but no net at all"
        )


class _DoubledSpelling(NamedTuple):
    """One ClinVar call the annotation writes two ways, and real values witnessing it."""

    #: What a clinician would call it, for the failure message.
    call: str
    #: The two offered spellings. Order carries no meaning: neither is the canonical one.
    spellings: tuple
    #: ``ClinVar_VCF_CLNSIG`` strings **copied from real annotated MAFs**, not constructed.
    values: tuple
    #: The spellings the measured corpus was found to *carry* — which is not always both.
    #:
    #: Recorded rather than derived, because it is a finding about the data and the point of
    #: some rows: the post-2023 conflicting spelling occurs in **zero** rows of the corpus, and
    #: that absence is #88's whole result. Each spelling named here must be witnessed by
    #: :attr:`values` and must be load-bearing, so the field cannot overstate the sample and
    #: the sample cannot quietly stop covering a spelling.
    carried: tuple
    #: Why one call has two spellings, since the two causes are unrelated.
    why: str


#: Every ClinVar call this vocabulary offers under two spellings.
#:
#: One table because there is one rule — *a term is offered in every spelling the annotation
#: writes it in* — reached twice by different routes. #88 found the 2023 rename; #99 found the
#: underscore-prefixed modifiers. They differ in provenance and in nothing the app does about
#: them, and holding them apart made the second look like a new principle when it was the
#: first applied again.
#:
#: The values are copied out of real files — 48 byte-distinct annotated MAFs, 16,776
#: ClinVar-annotated rows for the first three rows, and a wider sweep of 176 files and 159,580
#: rows for the four #103 added — because what is under test is a property of the *data* and
#: not one this repository could restate: which spelling reaches a MAF is decided by the
#: annotation.
#:
#: They are deliberately **not** required to cover both spellings of a pair. The
#: post-2023 conflicting spelling occurs in **zero** rows of that corpus, which is not a gap in
#: the sample but #88's whole finding — the offered spelling was the dead one. So the witness
#: clause below asks only that *some* single spelling be measurably insufficient.
#:
#: **This table is the whole set, and keeping it so is the point** — the rule is *every* call
#: with two spellings, so a call added to the vocabulary in a pair and not added here is a
#: silent exception to the one rule three tickets have now applied. #103 added four
#: (``risk_factor``, ``association``, ``protective``, ``confers_sensitivity``) and they are
#: below; ``test_every_offered_underscored_term_is_in_this_table`` is what stops the eighth
#: being forgotten, since nothing else would have failed.
#:
#: What #103 did change is the *preset* clause, and the reason is on the test rather than
#: here: these four were offered and kept by no preset, so the clause became all-or-none.
#: **#114 puts all four into both Broad presets** — Broad now keeps every ClinVar class except
#: the two benign ones — so every spelling tabled here is once again one a preset keeps. The
#: clause stays all-or-none regardless; see the test for why that is the property rather than a
#: concession to how the presets happen to stand.
_DOUBLED_SPELLINGS = (
    _DoubledSpelling(
        call="conflicting classifications",
        spellings=(
            "Conflicting_classifications_of_pathogenicity",
            "Conflicting_interpretations_of_pathogenicity",
        ),
        # Every observed value carries the pre-2023 spelling; the post-2023 one, which the app
        # offered alone until #88, appears nowhere in the corpus.
        values=(
            "Conflicting_interpretations_of_pathogenicity",
            "Conflicting_interpretations_of_pathogenicity|other",
            "Conflicting_interpretations_of_pathogenicity|risk_factor",
        ),
        # One spelling only, and that is the finding: 338 rows carry the pre-2023 name and none
        # carries the post-2023 one the app offered by itself until #88.
        carried=("Conflicting_interpretations_of_pathogenicity",),
        why=(
            "ClinVar renamed the call in 2023 and this repo pins no release, so which name "
            "lands in a MAF is decided by when the file was annotated (issue #88)"
        ),
    ),
    _DoubledSpelling(
        call="other",
        spellings=("other", "_other"),
        # The same ClinVar content, two spellings, in different files: 12 rows spell it with
        # the underscore, 7 without. No file measured carries both spellings of `other`, so on
        # those two files the offered `other` kept nothing at all.
        values=(
            "Conflicting_interpretations_of_pathogenicity|_other",
            "Conflicting_interpretations_of_pathogenicity|other",
            "Benign|_other",
            "Likely_benign|other",
            "Benign/Likely_benign|_other",
        ),
        # Both, in different files, and evenly: 14 rows spell it with the underscore, 14
        # without, no file carrying both. The underscored 14 are 7 rows in each of two files —
        # which are two annotation runs of **one sample**, so this spelling rests on one sample
        # plus the pipeline's own help listing `_low_penetrance`.
        carried=("other", "_other"),
        why=(
            "a modifier following a call keeps its leading underscore in some files and loses "
            "it in others, and the split does not strip it (issue #99)"
        ),
    ),
    _DoubledSpelling(
        call="drug response",
        spellings=("drug_response", "_drug_response"),
        # This one is mixed *within* one file: 19 bare against 2 underscored, so no per-file
        # rule would have covered it.
        values=(
            "Uncertain_significance|_drug_response",
            "Benign|_drug_response",
            "Likely_benign|drug_response|other",
        ),
        # Both, and here within one file: 64 rows bare against 4 underscored.
        carried=("drug_response", "_drug_response"),
        why="the same underscore artefact as `other`, but mixed within one file (issue #99)",
    ),
    # -- The four calls issue #103 put on offer. Same artefact, same rule, fourth through
    # -- seventh application of it. Values and row counts measured over 176 byte-distinct
    # -- real MAFs and 159,580 ClinVar-annotated rows; every one of these eight spellings is
    # -- load-bearing, keeping a real value its partner does not.
    _DoubledSpelling(
        call="risk factor",
        spellings=("risk_factor", "_risk_factor"),
        values=(
            "risk_factor",
            "Conflicting_classifications_of_pathogenicity|risk_factor",
            "Benign|_risk_factor",
            "Conflicting_interpretations_of_pathogenicity|_other|_risk_factor",
        ),
        # 247 rows bare across 102 files, 4 underscored across 2. The largest of the four,
        # and the term the ticket for #103 led with.
        carried=("risk_factor", "_risk_factor"),
        why="the same underscore artefact as `other` (issues #99, #103)",
    ),
    _DoubledSpelling(
        call="association",
        spellings=("association", "_association"),
        values=(
            "association",
            "Likely_benign|Affects|association",
            "Affects|_association",
            "Benign|_association|_confers_sensitivity",
        ),
        # 152 rows bare across 93 files, 4 underscored across 2 -- and `Affects|_association`
        # is the cell #103's ticket named as reachable by neither spelling of either term.
        carried=("association", "_association"),
        why="the same underscore artefact as `other` (issues #99, #103)",
    ),
    _DoubledSpelling(
        call="protective",
        spellings=("protective", "_protective"),
        values=(
            "protective",
            "Benign|protective",
            "Pathogenic|_protective",
        ),
        # 107 rows bare across 92 files, 2 underscored across 2.
        carried=("protective", "_protective"),
        why="the same underscore artefact as `other` (issues #99, #103)",
    ),
    _DoubledSpelling(
        call="confers sensitivity",
        spellings=("confers_sensitivity", "_confers_sensitivity"),
        values=(
            "confers_sensitivity",
            "confers_sensitivity|other",
            "Benign|_association|_confers_sensitivity",
        ),
        # 15 rows bare across 5 files, 2 underscored across 2. The one term of the four that
        # ClinVar's current documentation does not list -- a legacy value, placed in the
        # Drug Response class by the dev's call and flagged as such in the mapping.
        carried=("confers_sensitivity", "_confers_sensitivity"),
        why="the same underscore artefact as `other` (issues #99, #103)",
    ),
)


@pytest.mark.parametrize("pair", _DOUBLED_SPELLINGS, ids=lambda p: p.call)
def test_both_spellings_of_a_call_stay_on_offer(pair):
    """One ClinVar call, two live spellings — and a user cannot re-annotate to suit us.

    Replaces ``test_both_spellings_of_the_conflicting_term_stay_on_offer``, which held this
    property for the 2023 rename alone. The rename is now the first row of
    :data:`_DOUBLED_SPELLINGS`, asserted by this same test: one rule, every instance of it.

    Both halves matter. A spelling off the vocabulary cannot be selected at all, so a MAF
    written that way has those calls kept by nothing. A spelling missing from a Broad preset is
    subtler and was the state #88 found: the preset still reads as keeping the call, while
    whether it *reports* one depends on how the file was written.

    **The preset half is all-or-none, not always-present** — and that is #103's correction to
    this test rather than a loosening of it. The first three rows are all terms a Broad preset
    keeps, so requiring both spellings *in* the preset and requiring the preset to be
    consistent looked like one assertion. #103 then offered four calls that no preset kept,
    because whether a discovery report should carry a protective or risk annotation was the
    institute table's filtering column to answer and that column was not in this repository;
    demanding their presence would have forced a clinical decision through a spelling test.

    **#114 has since made that decision and every row here is now kept by both Broad presets,
    which does not restore the stronger form.** All-or-none is the property this test is *for* —
    a preset that keeps a call keeps every spelling of it — and always-present would make this
    guard fail the day a preset legitimately drops a class, reporting a spelling defect where
    there is a scope change. The two are only the same assertion while the preset happens to
    keep everything tabled, which is the coincidence #103 found and #114 restored.

    A pair in no preset passes; a pair half in one fails, which is the state #88 found and the
    only state this half was ever able to detect.
    """
    for spelling in pair.spellings:
        assert spelling in CLINVAR_OPTIONS, (
            f"{spelling!r} is not offered, so a MAF spelling {pair.call} that way cannot have "
            f"those variants kept at all — {pair.why}"
        )

    for name in ("SOFT_SOMATIC_PARAMS", "SOFT_GERMLINE_PARAMS"):
        kept = PRESETS[name]["filter_clinvar"]
        present = [s for s in pair.spellings if s in kept]
        if not present:
            continue
        missing = [s for s in pair.spellings if s not in kept]
        assert not missing, (
            f"{name} keeps {pair.call} as {present} but not {missing}, so whether the Broad "
            f"preset reports such a variant depends on how its file was written — {pair.why}"
        )


@pytest.mark.parametrize("pair", _DOUBLED_SPELLINGS, ids=lambda p: p.call)
def test_the_pair_reaches_what_one_spelling_alone_misses(pair):
    """The question #88's guard could not ask, asked of the matcher that decides.

    ``test_every_offered_clinvar_term_can_keep_a_variant`` asserts that a term keeps a variant
    annotated *exactly* that term. ``other`` passes it and always did — which is why nothing
    ever reported that on some real files it kept nothing. That guard is not wrong; it answers a
    different question, and this is the missing one: **a call is reached however its file
    spells it.**

    Run through the vendored matcher on real annotation values, never a regex of our own — the
    split is the pipeline's and a test restating it would pass while the pipeline moved.

    The second half is what makes this fail on a mutation, and the obvious weaker form of it
    does not: asking merely that *some* spelling alone be insufficient stays green even after
    the witnesses stop covering a spelling, because an underscored spelling always misses the
    bare-spelled values. So every spelling :attr:`_DoubledSpelling.carried` names must be
    **load-bearing** — it must keep a real value the rest of the pair does not — which fails
    both when a spelling leaves the vocabulary and when the sample stops witnessing it.
    """
    for value in pair.values:
        assert has_clinvar_term(value, list(pair.spellings)), (
            f"a variant annotated {value!r} is not kept by {list(pair.spellings)}, so "
            f"selecting {pair.call} reports variants according to how their file was written"
        )

    assert pair.carried, f"no spelling of {pair.call} is recorded as observed in real data"
    unoffered = sorted(set(pair.carried) - set(pair.spellings))
    assert not unoffered, (
        f"{pair.call} records carrying {unoffered}, which is not among its offered spellings"
    )

    for spelling in pair.carried:
        rest = [other for other in pair.spellings if other != spelling]
        needs_it = [
            value
            for value in pair.values
            if has_clinvar_term(value, [spelling]) and not has_clinvar_term(value, rest)
        ]
        assert needs_it, (
            f"no observed value needs {spelling!r} — every one is kept by {rest} alone. Either "
            f"{spelling!r} is carrying no weight, or these witnesses have stopped covering a "
            "spelling the corpus was measured to hold; re-measure before trusting this test"
        )


def test_every_offered_underscored_term_is_in_this_table():
    """The table has to keep up with the vocabulary, and nothing else made it.

    ``_DOUBLED_SPELLINGS`` says it holds *every* ClinVar call the vocabulary offers under two
    spellings — one table because there is one rule. Nothing enforced that. Both tests above
    are parametrised **over the table**, so a pair added to the vocabulary and forgotten here
    is not a failure, it is two tests that never run: the rule silently stops applying to the
    newest term while every message still claims it is universal.

    That is not hypothetical. #103 offered four new pairs and the first draft of it added
    none of them here, leaving the docstring false and the preset consistency of eight
    spellings unasserted. Review caught it; this is what catches the next one.

    Keyed on the underscore because that is the mechanical half of the rule — the 2023 rename
    is a pair whose names no rule derives, which is why the table names it by hand and why
    this test cannot be the only guard.
    """
    tabled = {spelling for pair in _DOUBLED_SPELLINGS for spelling in pair.spellings}
    underscored = [term for term in CLINVAR_OPTIONS if term.startswith("_")]
    assert underscored, "no underscored term is offered; this test is guarding nothing"

    missing = sorted(term for term in underscored if term not in tabled)
    assert not missing, (
        f"{missing} are offered as underscored spellings and appear in no _DOUBLED_SPELLINGS "
        "row, so the rule that a call is offered in every spelling its files use is asserted "
        "for every term but these. Add a row with real values measured from the corpus"
    )

    # And the other half: a bare partner offered without its underscored twin would leave the
    # pair unrepresentable here at all, which is the #99 defect in its original direction.
    for pair in _DOUBLED_SPELLINGS:
        unoffered = [s for s in pair.spellings if s not in CLINVAR_OPTIONS]
        assert not unoffered, (
            f"{pair.call} is tabled with spellings {unoffered} the vocabulary no longer "
            "offers; a row here must describe the control as it is"
        )


# ---------------------------------------------------------------------------
# The widgets themselves
# ---------------------------------------------------------------------------

_TAB_SCRIPT = """
import sys
sys.path.insert(0, {app!r})
import streamlit as st
from config.pipeline_params import pipeline_params
from page_modules.parameter_config import {function}

if "filter_params" not in st.session_state:
    st.session_state.filter_params = pipeline_params({arm!r})
{seed}
{function}({arm!r})
"""


def _render(function, arm="somatic", seed=""):
    """One parameter tab under ``AppTest``, over an optionally pre-seeded dict.

    A tab rather than the whole page, for the reason ``test_gene_lists.py`` gives: the
    full page consults the ``~/.mafigate`` cache, which decides the arm from whatever the
    developer last configured. ``seed`` runs inside the script so it lands before the
    render, exactly as the cache and an imported file do.
    """
    pytest.importorskip("streamlit", reason="streamlit not installed")
    from streamlit.testing.v1 import AppTest

    app = AppTest.from_string(
        _TAB_SCRIPT.format(
            app=str(STREAMLIT_APP), function=function, arm=arm, seed=seed
        )
    )
    app.run()
    assert not app.exception, [str(e.value) for e in app.exception]
    return app


_FULL_PAGE_SCRIPT = """
import sys
sys.path.insert(0, {app!r})
import streamlit as st
from page_modules import parameter_config

parameter_config.show_parameter_config_page()
"""


def _render_full_page(select_arm=None):
    """The whole parameter page, with the parameter cache stubbed out.

    Stubbed rather than pointed at a temporary ``HOME``, because the cache is read through
    a module-level path constant evaluated at import. Stubbing the cache functions is the
    smaller lie and it isolates what is under test: this file asserts the *default*, and a
    developer's own ``~/.mafigate`` deciding the arm has already cost this project a
    confusing failure.

    ``discard_stale_cache`` is stubbed for a stronger reason than isolation. Since issue
    #40 the page really does set an unstamped cache aside on the filesystem, so an
    unstubbed render would move the developer's own file — a test suite with a side effect
    on the machine it runs on.

    The patches are applied from here rather than by assignment inside the script, and
    that matters: the script runs in *this* process and mutates the shared module, so
    assignments there are permanent for the rest of the session. They used to be, and the
    stubbed ``save_parameters_to_cache`` silently disabled cache writes for every test
    that ran afterwards.
    """
    pytest.importorskip("streamlit", reason="streamlit not installed")
    from streamlit.testing.v1 import AppTest

    from page_modules import param_store, parameter_config

    # Two patch targets, and the split is not arbitrary. The first three are names the page
    # imported, so its own module dict is what a call resolves against and stubbing them on
    # `param_store` would not be seen. `discard_stale_cache` the page never names: it is
    # reached inside `param_store.show_discarded_cache_banner`, where the lookup is that
    # module's global — so that one has to be stubbed there, or the render moves the
    # developer's own `~/.mafigate` after all.
    with mock.patch.object(parameter_config, "load_parameters_from_cache", lambda: None), \
        mock.patch.object(parameter_config, "save_parameters_to_cache", lambda params: True), \
        mock.patch.object(parameter_config, "get_cache_info", lambda: None), \
        mock.patch.object(param_store, "discard_stale_cache", lambda: None):
        app = AppTest.from_string(_FULL_PAGE_SCRIPT.format(app=str(STREAMLIT_APP)))
        app.run(timeout=30)
        assert not app.exception, [str(e.value) for e in app.exception]

        if select_arm is not None:
            app.selectbox[0].select(select_arm).run(timeout=30)
            assert not app.exception, [str(e.value) for e in app.exception]
    return app


def test_the_classification_control_is_an_exclude_style_multiselect():
    """Its label says exclude, and it opens on the config's three values."""
    app = _render("show_basic_filters_tab")

    labels = [widget.label for widget in app.multiselect]
    classification = [label for label in labels if "lassification" in label]
    assert classification, f"no classification multiselect on the tab: {labels}"
    assert any("xclude" in label for label in classification), (
        "the classification control does not say it excludes. It is an exclude list "
        f"now, and a control labelled otherwise is the divergence restated: {labels}"
    )

    assert app.session_state["filter_params"]["filter_variant_classification"] == [
        "Silent",
        "IGR",
        "RNA",
    ]


def test_the_depth_control_does_not_name_a_column_it_does_not_read():
    """It was labelled ``Minimum Depth (DP)`` and the gate has never read ``DP`` (issue #127).

    The vendored clause is ``(t_alt_count + t_ref_count) >= coverage``, on both arms, and
    ``DP`` is not that sum by another name: across 322,913 rows of 157 real MAFs the two
    differ on 72.7% of them, ``DP`` the larger on 57.3%, because ``DP`` counts every sample
    and this counts the tumour's reads. At this control's own default of 50 the two answers
    disagree on 2,117 of those rows.

    Both the label and the help text are checked, because the help is where a reader goes when
    a label is ambiguous, and a label that stopped saying ``DP`` over help that still did would
    be the same claim one hover away. ``DP`` may still be *mentioned* — the help says it is
    **not** what this reads, which is the useful thing to say to someone who has the column on
    screen — so what is asserted is that the two are not equated, by requiring the disclaimer
    wherever the name appears.
    """
    app = _render("show_basic_filters_tab")

    depth = [widget for widget in app.number_input if "epth" in widget.label]
    assert len(depth) == 1, f"expected one depth control, got {[w.label for w in depth]}"
    control = depth[0]

    assert "DP" not in control.label, (
        f"the depth control is labelled {control.label!r}, naming a column the gate does not "
        "read — it sums t_alt_count and t_ref_count"
    )
    if "DP" in (control.help or ""):
        assert "not" in control.help, (
            f"the depth control's help mentions DP without saying it is not what this reads: "
            f"{control.help!r}"
        )


def test_clearing_the_classification_control_excludes_nothing():
    """Empty means empty — it must not be back-filled with a catch-all.

    ``validate_multiselect_params`` used to rewrite every empty selection to ``["All"]``,
    which is what made the sentinel reachable without the user ever choosing it.
    """
    app = _render("show_basic_filters_tab")
    widget = next(w for w in app.multiselect if "lassification" in w.label)
    widget.set_value([]).run()
    assert not app.exception, [str(e.value) for e in app.exception]

    assert app.session_state["filter_params"]["filter_variant_classification"] == []


@pytest.mark.parametrize("arm", ARMS)
def test_the_clinvar_control_explains_its_duplicate_entry(arm):
    """The screen must account for the two conflicting-classification entries.

    Offering both spellings is the right answer to a vocabulary ClinVar renamed under the
    app (see ``CLINVAR_OPTIONS``), but it leaves a clinician reading two near-identical
    lines in one dropdown with no way to tell them apart — which is its own small lie
    about what a control does. The tooltip is where that gets said, so this asserts the
    tooltip says it rather than that the constant is right.

    Checked on both arms because the ClinVar control is the one guideline source drawn on
    each of them, from two separate entries in ``_GUIDELINE_CONTROLS`` — the shape that
    let the germline copy drift into rendering an ESCAT control the germline filter has no
    argument for.
    """
    app = _render("show_clinical_filters_tab", arm=arm)

    clinvar = [widget for widget in app.multiselect if "ClinVar" in widget.label]
    assert clinvar, f"{arm}: no ClinVar control on the clinical filters tab"

    help_text = (clinvar[0].help or "").lower()

    # Both words, because on screen the two entries differ by exactly this one word --
    # "classifications" against "interpretations" -- and a tooltip naming only one of them
    # has not told the reader which two lines it is talking about. A bare "2023" would
    # pass on any sentence mentioning the year.
    for word in ("classifications", "interpretations", "2023"):
        assert word in help_text, (
            f"{arm}: the ClinVar control offers two spellings of the conflicting-"
            f"classification call and its tooltip never says {word!r}. The two entries "
            "read as two different classifications on screen. Tooltip was: "
            f"{clinvar[0].help!r}"
        )

    # The rename is what makes them one call rather than two, so the tooltip has to say a
    # rename happened -- naming both words and the year while implying they are separate
    # findings would be worse than silence.
    assert "renamed" in help_text or "older name" in help_text, (
        f"{arm}: the tooltip names both spellings but never says they are one call under "
        f"an old and a new name. Tooltip was: {clinvar[0].help!r}"
    )


@pytest.mark.parametrize("arm", ARMS)
def test_both_clinvar_controls_write_one_parameter(arm):
    """Two widgets, one ``filter_clinvar`` — and the saved list partitions itself.

    The single most breakable thing issue #103 introduced. Every other guideline key is drawn
    by exactly one control, and the render loop wrote each widget's value straight into
    ``filter_params[key]``; with two controls on one key that assignment silently discards the
    first, and the report is cut on half the user's selection with nothing on screen to say
    so. The loop accumulates instead, and this is what proves it.

    It also pins the seeding, which is the half with no code of its own: each widget's default
    is ``filter_terms`` of the *whole* saved list against its own options, so one stored
    ``filter_clinvar`` — from a preset, a ``~/.mafigate`` cache or an uploaded JSON — arrives
    split across the two controls without anything partitioning it explicitly.

    Both arms, because the ClinVar control is the one guideline source drawn on each of them
    and ``_GUIDELINE_CONTROLS`` holds a separate copy per arm.
    """
    seed = (
        "st.session_state.filter_params['filter_clinvar'] = "
        "['Pathogenic', 'risk_factor', 'no_classification_for_the_single_variant']"
    )
    app = _render("show_clinical_filters_tab", arm=arm, seed=seed)

    controls = [widget for widget in app.multiselect if "ClinVar" in widget.label]
    assert len(controls) == 2, (
        f"{arm}: expected two ClinVar controls, found {[w.label for w in controls]}"
    )

    # Order is load-bearing on screen and in the vocabulary, and one sibling test reads
    # `clinvar[0]` for the rename tooltip -- so it is asserted here rather than assumed there.
    assert controls[0].label == "Keep ClinVar Classifications", controls[0].label
    assert controls[1].label == "Keep Other ClinVar Assertions", controls[1].label

    assert list(controls[0].value) == ["Pathogenic"], (
        f"{arm}: the classifications control was seeded {list(controls[0].value)}; it should "
        "hold only the pathogenicity call from the saved list"
    )
    assert sorted(controls[1].value) == [
        "no_classification_for_the_single_variant",
        "risk_factor",
    ], (
        f"{arm}: the other-assertions control was seeded {sorted(controls[1].value)}; it "
        "should hold exactly the saved terms that are not pathogenicity calls"
    )

    written = app.session_state["filter_params"]["filter_clinvar"]
    assert sorted(written) == [
        "Pathogenic",
        "no_classification_for_the_single_variant",
        "risk_factor",
    ], (
        f"{arm}: the two controls wrote {sorted(written)} into filter_clinvar. Both "
        "selections must survive — a parameter drawn by two widgets is the union of them, "
        "and anything less means the report is cut on part of what the user asked for"
    )

    # And a change in the second control reaches the parameter, which assignment-in-loop
    # would also have produced -- so this is asserted on the *first*, where it would not.
    controls[0].set_value([]).run()
    assert not app.exception, [str(e.value) for e in app.exception]
    assert sorted(app.session_state["filter_params"]["filter_clinvar"]) == [
        "no_classification_for_the_single_variant",
        "risk_factor",
    ], (
        f"{arm}: emptying the classifications control left "
        f"{sorted(app.session_state['filter_params']['filter_clinvar'])}; the other "
        "control's selection must be untouched by it"
    )


@pytest.mark.parametrize("arm", ARMS)
def test_the_other_clinvar_control_says_what_it_is_for(arm):
    """The second control has to account for itself on screen.

    A clinician meeting *Keep Other ClinVar Assertions* has no way to know from the label
    that selecting one of its terms can bring in a variant carrying **no pathogenicity call
    at all** — which is the whole of what the control adds, and the only reason to open it.
    The tooltip is where that gets said, so this reads the tooltip the page drew rather than
    the constant behind it.

    Deliberately not asserting a gloss per term: which class a term belongs to is
    ``components/clinical_summary.py``'s to say, and a second copy of it in a tooltip is the
    drift #79 spent a ticket removing from the Help page.
    """
    app = _render("show_clinical_filters_tab", arm=arm)

    other = [
        widget
        for widget in app.multiselect
        if widget.label == "Keep Other ClinVar Assertions"
    ]
    assert other, f"{arm}: no other-assertions control on the clinical filters tab"

    help_text = (other[0].help or "").lower()
    assert "not pathogenicity calls" in help_text, (
        f"{arm}: the tooltip never says these are not pathogenicity calls, which is the one "
        f"thing distinguishing this control from the one beside it: {other[0].help!r}"
    )
    assert "no pathogenicity call" in help_text, (
        f"{arm}: the tooltip never says a variant carrying no pathogenicity call can reach "
        f"the report through this control, which is what it is for: {other[0].help!r}"
    )


@pytest.mark.parametrize("arm", ARMS)
def test_emptying_every_guideline_source_warns_rather_than_blocks(arm):
    """The now-reachable all-empty state is warned at the widget, not refused.

    It is expressible on the pipeline's own command line, so blocking it would make the
    app stricter than the pipeline. The user is told instead — because the report
    collapsing to the pathogenic-rescue floor is otherwise indistinguishable from a
    filter that simply found nothing.
    """
    seed = "\n".join(
        f"st.session_state.filter_params[{key!r}] = []"
        for key in GUIDELINE_SOURCES[arm]
    )
    app = _render("show_clinical_filters_tab", arm=arm, seed=seed)

    warnings = " ".join(warning.value for warning in app.warning)
    assert warnings, f"{arm}: no warning when every guideline source is empty"
    assert "pathogenic" in warnings.lower(), (
        f"{arm}: the warning does not name the pathogenic-rescue floor, which is what "
        f"the report has collapsed to: {warnings!r}"
    )

    # Warned, not blocked: the parameters are left exactly as chosen.
    for key in GUIDELINE_SOURCES[arm]:
        assert app.session_state["filter_params"][key] == [], (
            f"{arm}: {key} was back-filled rather than left empty"
        )


@pytest.mark.parametrize("arm", ARMS)
def test_a_populated_guideline_selection_is_not_warned_about(arm):
    """The warning fires on the state it names, not on every render."""
    app = _render("show_clinical_filters_tab", arm=arm)
    warnings = " ".join(warning.value for warning in app.warning)
    assert "pathogenic-rescue" not in warnings, (
        f"{arm}: warned about an all-empty selection at the contract's own defaults: "
        f"{warnings!r}"
    )


#: Which tab renders which multiselects, so the assertion below is made only about the
#: controls the tab under test actually draws.
_TAB_CONTROLS = [
    ("show_basic_filters_tab", "somatic", ("filter_variant_classification",)),
    ("show_clinical_filters_tab", "somatic", GUIDELINE_SOURCES["somatic"]),
    ("show_clinical_filters_tab", "germline", GUIDELINE_SOURCES["germline"]),
]


@pytest.mark.parametrize("function,arm,keys", _TAB_CONTROLS)
def test_a_stale_sentinel_in_the_parameters_does_not_break_the_page(function, arm, keys):
    """A ``~/.mafigate`` cache written before this ticket still opens.

    Streamlit raises if a ``multiselect``'s default is not among its options, so simply
    deleting the option would turn every pre-#36 cache entry into a stack trace on the
    parameter page — the render inside ``_render`` asserts no exception, which is half of
    what this test is for. The value is dropped instead, and dropping it is safe rather
    than merely convenient: issue #33 already made it a value nothing matches, so
    ``["All"]`` and ``[]`` decide every row identically.
    """
    seed = "\n".join(
        f"st.session_state.filter_params[{key!r}] = ['All']" for key in OPTION_LISTS
    )
    app = _render(function, arm=arm, seed=seed)

    for key in keys:
        value = app.session_state["filter_params"][key]
        assert SENTINEL not in value, (
            f"{key} kept the stale sentinel {SENTINEL!r} after a render: {value}"
        )


def test_the_sentinel_is_stripped_from_keys_the_current_arm_never_renders():
    """Including the ones no widget on this arm ever sees.

    ``selectable()`` only cleans a value on its way into a widget, so a germline session
    never touches ``filter_cancervar``. Left alone, the page's auto-save would write the
    stale value straight back to the cache and into every parameter file the user
    exported from that session — the sentinel outliving the ticket that deleted it.
    """
    from config.vocabularies import validate_multiselect_params

    stale = {key: [SENTINEL] for key in OPTION_LISTS}
    stale["filter_clinvar"] = ["Pathogenic", SENTINEL]
    stale["filter_renovo"] = SENTINEL  # a hand-edited JSON's bare string

    cleaned = validate_multiselect_params(dict(stale))

    for key in OPTION_LISTS:
        assert SENTINEL not in cleaned[key], f"{key}: {cleaned[key]}"
    assert cleaned["filter_clinvar"] == ["Pathogenic"], "a real term was stripped with it"
    assert cleaned["filter_renovo"] == []


def test_an_empty_selection_is_not_back_filled_by_the_validator():
    """The back-fill is what made the sentinel reachable without anyone choosing it."""
    from config.vocabularies import validate_multiselect_params

    cleaned = validate_multiselect_params({key: [] for key in OPTION_LISTS})
    assert all(cleaned[key] == [] for key in OPTION_LISTS), cleaned


def test_the_germline_arm_offers_no_escat_control():
    """ESCAT is somatic-only, so a germline ESCAT control cannot do anything.

    ``germline_filters()`` takes no ESCAT argument — an app-side germline ESCAT clause
    was the single largest divergence on real data (540 rows, 51% of all attributed
    divergence, issue #28), and issue #33 deleted it. What survived was the *widget*: a
    control the germline arm still drew, which a user could set and which provably could
    not change a single row. Same lie about what a control does as the sentinel, so it
    goes on the same ticket.
    """
    app = _render("show_clinical_filters_tab", arm="germline")
    labels = [widget.label for widget in app.multiselect]
    assert not [label for label in labels if "ESCAT" in label], (
        f"the germline tab still offers an ESCAT control: {labels}"
    )

    somatic = _render("show_clinical_filters_tab", arm="somatic")
    assert [label for label in (w.label for w in somatic.multiselect) if "ESCAT" in label], (
        "ESCAT vanished from the somatic arm too, where it is a real guideline source"
    )
