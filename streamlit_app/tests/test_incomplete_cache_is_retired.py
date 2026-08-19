"""A stamped cache that is not a complete set for its arm is retired (issue #286).

Issue #279 decided what such a file *means*, and the decision is the whole design: the app
never asks whether an empty keep-list was meant. It asks whether the dict is one this app
could have written **without inventing** — and issue #276 established that inventing is
exactly what produces an incomplete one. So:

* **chosen empty** → a *complete* dict with ``[]`` in it, restored untouched;
* **defect-written** → a dict *missing keys its own arm's contract carries*, set aside.

Both sides are pinned below, on both arms, because a repair that cleared the reported
symptom by re-filling empty lists would be issue #36's deleted backfill returning and would
overrule the one user this must not touch.

Why the condition is read on the raw document
---------------------------------------------
Issue #280 completes a partial dict at the top of ``show_parameter_config_page``. Anything
measured *after* that completion is complete by construction, so a completeness test placed
behind it could never fire — it would guard a door while the room was being repaired out
from under it. The test therefore runs in ``load_parameters_from_cache`` and
``discard_stale_cache``, both of which see the file as it is on disk.

Why the load path refuses as well as the banner moving
------------------------------------------------------
The boot order is ``initialize_session_state`` (which restores the cache) and then
``render_header`` (which draws the banner) — ``MAFigate.main``. A read that adopted an
incomplete cache and left the moving to the banner would open the session on the very
parameters it was about to set aside, and the user would be told their settings were
retired while looking at them.

On the instrument
-----------------
``test_a_render_of_the_contract_writes_a_complete_cache`` is the false-positive guard, and
four passing render probes prove nothing on their own — this repo's recorded failure mode
is a guard that never could fail. ``test_the_completeness_check_can_see_a_loss`` is the
live fire: it neuters ``complete_params`` and shows the same harness reporting the loss,
naming the two keys issue #276 measured.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from unittest import mock

import pytest

from config.param_migration import PARAM_SCHEMA_VERSION, SCHEMA_VERSION_KEY
from config.pipeline_params import ARMS, pipeline_params
from config.vocabularies import GUIDELINE_SOURCES

# The artifact and the render harness are issue #276's, imported rather than copied: a
# second transcription of the dev's cache is a second thing to drift.
from tests.test_param_migration import cache_at  # noqa: F401  (a fixture, used by name)
from tests.test_partial_parameter_dict import DEV_CACHE, _page

STREAMLIT_APP = Path(__file__).resolve().parents[1]

#: What the dev's cache is missing, measured: the two keys no widget on the germline arm
#: draws, so the page could not re-invent them and the file never regained them.
DEV_CACHE_MISSING = ("skip_civic", "vaf_threshold")


def _stamped(params, timestamp="2026-08-19T16:57:59", version=PARAM_SCHEMA_VERSION):
    """A cache document as this app writes one."""
    return {
        SCHEMA_VERSION_KEY: version,
        "parameters": params,
        "timestamp": timestamp,
        "app_version": "1.0.0",
    }


def _emptied(arm):
    """A complete parameter set for ``arm`` whose every guideline source is empty.

    The state a user is entitled to choose, and the one a symptom-clearing repair would
    silently overrule. Complete, so the retirement condition never looks at it.
    """
    params = pipeline_params(arm)
    for key in GUIDELINE_SOURCES[arm]:
        params[key] = []
    return params


# ---------------------------------------------------------------------------
# The condition
# ---------------------------------------------------------------------------


def test_the_shape_the_defect_wrote_is_retired(cache_at):
    """The dev's own cache, as it was on disk at 2026-08-19T16:57:59.

    It carries this version's stamp and every number, so nothing about it reads as old.
    What gives it away is structural: the germline arm's contract names ``skip_civic`` and
    ``vaf_threshold`` and the file carries neither, because no widget on that arm draws
    them — so the page could not invent them, and the file could not heal.
    """
    from page_modules.param_store import (
        DISCARDED_INCOMPLETE,
        discard_stale_cache,
        load_parameters_from_cache,
        missing_contract_keys,
    )

    cache_at.write_text(json.dumps(_stamped(DEV_CACHE)))

    assert missing_contract_keys(DEV_CACHE) == DEV_CACHE_MISSING

    discarded = discard_stale_cache()

    assert discarded is not None, (
        "the cache the app itself wrote by inventing values was restored as though the "
        "user had chosen it"
    )
    assert discarded.reason == DISCARDED_INCOMPLETE
    assert discarded.missing == DEV_CACHE_MISSING
    assert not cache_at.exists(), "the retired cache is still there to be read next boot"
    assert load_parameters_from_cache() is None


def test_the_load_path_refuses_it_too_not_just_the_banner(cache_at):
    """The boot restores the cache *before* the header draws the banner.

    So a read that adopted an incomplete cache and left the moving to the banner would open
    the session on the parameters it was in the act of retiring. Asserted without calling
    ``discard_stale_cache`` at all, which is the ordering the app actually runs in.
    """
    from page_modules.param_store import load_parameters_from_cache

    cache_at.write_text(json.dumps(_stamped(DEV_CACHE)))

    assert load_parameters_from_cache() is None
    assert cache_at.exists(), "the read moved the file; that is the banner's job, not its"


@pytest.mark.parametrize("arm", ARMS)
def test_a_complete_cache_with_every_source_emptied_on_purpose_is_left_alone(cache_at, arm):
    """The user this change must not touch, on both arms.

    Emptying every guideline source is a legal, expressible choice — the pipeline's own
    ``--filter_cancervar ""`` — which is why issue #36 deleted the empty→``["All"]``
    backfill and why the page *warns* about the state instead of preventing it. Their dict
    is complete, so the condition never looks at it: nothing is repaired, back-filled or
    overruled, and it comes back exactly as empty as they left it.
    """
    from page_modules.param_store import discard_stale_cache, load_parameters_from_cache

    chosen = _emptied(arm)
    cache_at.write_text(json.dumps(_stamped(chosen)))

    assert discard_stale_cache() is None, (
        f"a complete {arm} dict with every guideline source emptied on purpose was set "
        "aside; the condition read emptiness, which issue #279 decided it must never do"
    )
    assert cache_at.exists()
    restored = load_parameters_from_cache()
    assert restored == chosen
    for key in GUIDELINE_SOURCES[arm]:
        assert restored[key] == [], f"{key} was chosen empty and came back {restored[key]!r}"


def test_a_newer_format_is_left_alone_even_though_it_is_incomplete(cache_at):
    """A format this version has never seen is not ours to judge — completeness included.

    The keys a newer MAFigate's contract names are not this one's, so measuring its file
    against this contract answers a question nobody asked. It is declined by the reader and
    left exactly where it is, because the user is very likely to run that version again.
    """
    from page_modules.param_store import (
        discard_stale_cache,
        load_parameters_from_cache,
        missing_contract_keys,
    )

    params = {"sample_type": "somatic", "min_depth": 7}
    assert missing_contract_keys(params), "this fixture is meant to be incomplete"
    cache_at.write_text(json.dumps(_stamped(params, version=PARAM_SCHEMA_VERSION + 1)))

    assert discard_stale_cache() is None
    assert cache_at.exists()
    assert load_parameters_from_cache() is None


def test_a_cache_that_states_no_arm_is_retired(cache_at):
    """There is no contract such a file could be complete against.

    The app reads an absent ``sample_type`` as somatic everywhere, but it never *writes*
    one that way: the Sample Type control sets the key on every render of the parameters
    page, which is the only thing that writes the cache. So a stamped file without an arm
    did not come from this app's own hands either.
    """
    from page_modules.param_store import discard_stale_cache, missing_contract_keys

    params = dict(pipeline_params("somatic"))
    params.pop("sample_type")
    assert missing_contract_keys(params) == ("sample_type",)

    cache_at.write_text(json.dumps(_stamped(params)))
    discarded = discard_stale_cache()

    assert discarded is not None
    assert "which analysis arm" in discarded.summary(), discarded.summary()


def test_the_condition_reads_the_contract_rather_than_a_restated_key_list():
    """A parameter added to the contract must not silently stop being required.

    The failure this forbids is the quiet one: a new key is wired into
    ``pipeline_params`` and a list of names copied into ``param_store`` goes on describing
    the old contract, so the very defect this ticket retires becomes invisible again for
    whichever parameter is newest.
    """
    from page_modules import param_store

    invented = "filter_something_new"
    contract = pipeline_params("somatic")
    contract[invented] = ["value"]

    with mock.patch.object(param_store, "pipeline_params", lambda arm: dict(contract)):
        params = pipeline_params("somatic")
        assert param_store.missing_contract_keys(params) == (invented,), (
            "a key added to the contract was not required of a cache, so param_store is "
            "measuring against a list of its own rather than against the contract"
        )


# ---------------------------------------------------------------------------
# What the banner says
# ---------------------------------------------------------------------------


def _retirement_banner(cache_at, params):
    from page_modules.param_store import discard_stale_cache

    cache_at.write_text(json.dumps(_stamped(params)))
    discarded = discard_stale_cache()
    assert discarded is not None
    return discarded


def test_the_retirement_banner_does_not_reuse_the_unstamped_wording(cache_at):
    """Three of that sentence's claims are false of this case, and one is harmful.

    It says the file was written *by an older MAFigate* — this one wrote it; that *nothing
    in the file says which version wrote it* — it is stamped; and it tells the user to
    upload it back under *Presets → Upload Custom Preset*, which for a v1-stamped file is a
    verbatim no-op, so following the app's own remedy would restore the damage.
    """
    summary = _retirement_banner(cache_at, DEV_CACHE).summary()

    for false_claim in (
        "by an older MAFigate",
        "nothing in the file says",
        "Upload Custom Preset",
        "discarded, not migrated",
    ):
        assert false_claim not in summary, (
            f"the retirement banner still says {false_claim!r}, which is not true of a "
            f"cache this MAFigate wrote and stamped: {summary}"
        )


def test_the_retirement_banner_names_the_evidence_and_where_the_file_went(cache_at):
    """What the user needs: which parameters were not theirs, and that nothing is lost."""
    discarded = _retirement_banner(cache_at, DEV_CACHE)
    summary = discarded.summary()

    assert "not the ones you chose" in summary, summary
    for key in DEV_CACHE_MISSING:
        assert key in summary, f"{key} was filled in for the user without being named: {summary}"
    assert "germline" in summary, summary
    assert discarded.kept_at in summary, summary
    assert "2026-08-19 16:57" in summary, summary


def test_the_rescue_floor_clause_is_drawn_only_when_every_source_was_empty(cache_at):
    """The backstop for a user who never opened the parameters page.

    ``_warn_if_every_guideline_source_is_empty`` fires only there, so someone who loaded a
    MAF and cut a report while this cache was in effect was never told the report had
    collapsed to the pathogenic-rescue floor. The clause is conditional because a cache can
    be incomplete *and* keep its keep-lists, and a sentence that said this of that file
    would be false.
    """
    collapsed = _retirement_banner(cache_at, DEV_CACHE).summary()
    assert "pathogenic-rescue" in collapsed, collapsed

    kept_something = dict(DEV_CACHE, filter_clinvar=["Pathogenic"])
    Path(cache_at).unlink(missing_ok=True)
    partial = _retirement_banner(cache_at, kept_something).summary()
    assert "pathogenic-rescue" not in partial, (
        "the banner told a user their report had collapsed to the rescue floor while "
        f"their ClinVar keep-list was populated: {partial}"
    )


def test_an_unstamped_cache_still_gets_its_own_wording(cache_at):
    """Issue #40's banner is not replaced, only branched — its case is still real."""
    from page_modules.param_store import DISCARDED_UNSTAMPED, discard_stale_cache

    cache_at.write_text(
        json.dumps({"parameters": {"sample_type": "germline", "min_depth": 50}})
    )

    discarded = discard_stale_cache()

    assert discarded.reason == DISCARDED_UNSTAMPED
    assert "by an older MAFigate" in discarded.summary()
    assert "Upload Custom Preset" in discarded.summary()


# ---------------------------------------------------------------------------
# The false positive, and the live fire that proves the guard can see one
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("arm", ARMS)
def test_a_render_of_the_contract_writes_a_complete_cache(arm):
    """The false-positive guard: what the app writes must never trip its own condition.

    One full render of each arm, through issue #276's harness, and the assertion is on the
    dict the auto-save actually handed the cache — not on session state, because the cache
    is what the next boot reads.
    """
    from page_modules.param_store import missing_contract_keys

    with _page(pipeline_params(arm)) as (_app, written):
        assert written, "the page did not auto-save, so this is not measuring the write"
        assert missing_contract_keys(written[-1]) == (), (
            f"one ordinary render of the {arm} contract wrote a cache this app would then "
            "retire, so the condition would fire on every user on their next boot"
        )


@pytest.mark.parametrize("arm", ARMS)
def test_a_render_of_a_deliberately_emptied_set_writes_a_complete_cache(arm):
    """The same guard for the user whose choice must survive, on both arms.

    Their settings pass through the page and the auto-save on every render, so if that
    round trip dropped a key the condition would retire the very choice it exists to
    protect.
    """
    from page_modules.param_store import missing_contract_keys

    with _page(_emptied(arm)) as (_app, written):
        assert written
        assert missing_contract_keys(written[-1]) == (), (
            f"a render of a complete {arm} set with every guideline source emptied on "
            "purpose wrote an incomplete cache, so that user's choice would be retired"
        )
        for key in GUIDELINE_SOURCES[arm]:
            assert written[-1][key] == [], f"{key} was re-filled behind the user"


@pytest.mark.parametrize(
    "preset",
    [
        "SOFT_SOMATIC_PARAMS",
        "SOFT_GERMLINE_PARAMS",
        "CLINICAL_SOMATIC_PARAMS",
        "CLINICAL_GERMLINE_PARAMS",
    ],
)
def test_a_render_of_a_preset_writes_a_complete_cache(preset):
    """The false positive issue #281 makes real, and the one nothing else here would catch.

    **All four presets lack ``skip_civic``.** So a preset dict is itself incomplete by this
    condition's measure, and if one ever reached the cache unaltered every user who loaded
    a preset would have it retired on their next boot — the condition firing on the app's
    own shipped settings.

    It does not, and this pins why rather than trusting it: the presets arrive through
    ``adopt_parameters``, which reruns before the auto-save, so what is written is always a
    dict that has been through ``complete_params`` at the top of the page. The assertion is
    on the write, because the write is what the next boot reads.
    """
    from config import presets
    from page_modules.param_store import missing_contract_keys

    params = json.loads(json.dumps(getattr(presets, preset)))
    assert "skip_civic" not in params, (
        f"{preset} now carries skip_civic, so this test no longer measures the case issue "
        "#281 found; check the other three before deleting it"
    )

    with _page(params) as (_app, written):
        assert written
        assert missing_contract_keys(written[-1]) == (), (
            f"a render of {preset} wrote a cache this app would retire, so every user who "
            "loaded that preset would lose it on their next boot"
        )


def test_the_completeness_check_can_see_a_loss_through_the_render_harness():
    """Live fire. Without this, the four probes above prove nothing.

    The completion is made to miss exactly one key and the same harness is asked the same
    question about the same seed. It must report that key — which is what makes the passing
    probes above evidence rather than decoration.

    The key is ``skip_civic`` for the reason that made the dev's cache permanent: **no
    widget draws it**, on either arm. So the page renders perfectly happily without it and
    the auto-save writes the gap straight to disk, which is the exact shape this ticket
    retires. Removing the completion altogether would not do: since issue #280 every
    control *indexes* its key, so a page that reaches a widget with a key missing raises
    rather than inventing a value — a stronger guarantee, but not one this harness can be
    calibrated against.
    """
    from config.param_migration import complete_params
    from page_modules import parameter_config
    from page_modules.param_store import missing_contract_keys

    def _completion_that_misses_one(params):
        complete_params(params)
        params.pop("skip_civic", None)
        return params

    with mock.patch.object(
        parameter_config, "complete_params", _completion_that_misses_one
    ):
        with _page({"sample_type": "germline"}) as (_app, written):
            assert written, "the page did not auto-save, so nothing is being measured"
            missing = missing_contract_keys(written[-1])

    assert missing == ("skip_civic",), (
        "a render that dropped a contract key wrote a cache the condition called complete "
        f"(saw {missing!r}); the harness cannot see the loss it is supposed to catch"
    )


@pytest.mark.parametrize("arm", ARMS)
def test_every_contract_key_is_one_the_condition_would_miss(arm):
    """Mutation, key by key: removing any one of them must be seen.

    A condition that noticed only the two keys the dev's cache happened to lack would be a
    restatement of one artifact rather than a reading of the contract, and would go on
    passing while a newer parameter went missing.
    """
    from page_modules.param_store import missing_contract_keys

    contract = pipeline_params(arm)
    assert missing_contract_keys(contract) == ()

    for key in contract:
        mutated = {k: v for k, v in contract.items() if k != key}
        assert missing_contract_keys(mutated) == (key,), (
            f"a {arm} cache missing {key} was not seen as incomplete"
        )


# ---------------------------------------------------------------------------
# The route, end to end
# ---------------------------------------------------------------------------


def test_a_real_boot_retires_the_dev_cache_and_opens_at_the_contract(tmp_path):
    """The whole recovery, through ``MAFigate.main`` rather than through its parts.

    Nothing is stubbed but the cache location: the boot reads the file, refuses it, the
    header moves it and draws the banner, and the session opens on the contract. An earlier
    superseded copy is placed first, so the collision numbering is exercised where it
    actually matters — the dev's ``~/.mafigate`` already holds one from 2026-06-16, and
    overwriting it would take back the promise that moving rather than deleting makes.
    """
    pytest.importorskip("streamlit", reason="streamlit not installed")
    from streamlit.testing.v1 import AppTest

    from page_modules import param_store

    home = tmp_path / ".mafigate"
    home.mkdir()
    cache_file = home / "cached_parameters.json"
    cache_file.write_text(json.dumps(_stamped(DEV_CACHE)))
    earlier = home / "cached_parameters.superseded.json"
    earlier.write_text(json.dumps({"parameters": {"sample_type": "somatic"}}))

    with mock.patch.object(param_store, "CACHE_DIR", home), mock.patch.object(
        param_store, "PARAMS_CACHE_FILE", cache_file
    ), mock.patch.dict(os.environ):
        # Otherwise the boot auto-loads a MAF and this stops being a test of the cache.
        os.environ.pop("MAFIGATE_OPEN_FILE", None)
        app = AppTest.from_file(str(STREAMLIT_APP / "MAFigate.py"), default_timeout=60)
        app.run()

    assert not app.exception, [str(e.value) for e in app.exception]

    assert not cache_file.exists(), "the boot left the retired cache where the next one reads it"
    kept = home / "cached_parameters.superseded-2.json"
    assert kept.exists(), sorted(p.name for p in home.iterdir())
    assert json.loads(kept.read_text())["parameters"] == DEV_CACHE
    assert json.loads(earlier.read_text())["parameters"] == {"sample_type": "somatic"}, (
        "the earlier superseded copy was overwritten"
    )

    assert app.session_state["filter_params"] == pipeline_params("somatic")
    assert app.session_state["filter_params"]["filter_clinvar"] == [
        "Pathogenic",
        "Likely_pathogenic",
    ], "the app opened on a ClinVar keep-list nobody chose"

    warnings = [str(warning.value) for warning in app.warning]
    assert any("not the ones you chose" in text for text in warnings), (
        f"the cache was retired without telling the user; saw {warnings}"
    )
