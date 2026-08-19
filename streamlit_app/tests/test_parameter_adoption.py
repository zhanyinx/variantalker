"""No path changes the arm without saying so (issue #133).

Three controls replace ``filter_params`` *whole* — a preset button, an uploaded custom
preset, and ``📥 Reload from Cache`` — and every parameter set carries its own
``sample_type``. So all three could move a user between the germline and somatic arms, and
all three did it in silence. Reported from a real session: the arm was set to germline and
*"then it automatically switched to somatic"*.

Why silence was the whole distance to the symptom: the arm decides which guideline sources
are read, which thresholds apply, which gene panel is offered and which columns the filter
requires. A germline MAF filtered on the somatic arm draws
``❌ PATHOGENIC RETENTION DEGRADED — CancerVar column not found``, which is true of the
somatic arm and blames the file for the absence of a column it was never supposed to carry.

One rule, three call sites, so this file asserts the rule rather than three fixes:

* each site announces what it loaded, and names the arm when the arm moved;
* the announcement survives ``st.rerun``, which is where the old ones died;
* and no *fourth* site can be added without failing the guard at the bottom.

On the instrument
-----------------
The end-to-end tests drive the real page through ``AppTest``, and one thing about that
instrument had to be measured before any of them could be trusted: **``AppTest`` reports
elements drawn immediately before ``st.rerun``** when they sit inside a container, though a
real user never sees them. A one-element probe settles it — a ``st.success`` inside an
``st.expander`` followed by ``st.rerun`` is reported; the same call at page level is not.
That is why the app's ``✅ … parameters loaded!`` read as working while showing nobody
anything, and why these tests assert the notice is *consumed* from session state as well as
drawn: a consumed key can only have been read on a frame that survived.
"""

from __future__ import annotations

import ast
import inspect
from pathlib import Path
from unittest import mock

import pytest

STREAMLIT_APP = Path(__file__).resolve().parents[1]

_PAGE_SCRIPT = """
import sys
sys.path.insert(0, {app!r})
import streamlit as st
from config.pipeline_params import pipeline_params
from page_modules import parameter_config

if "filter_params" not in st.session_state:
    st.session_state.filter_params = pipeline_params({arm!r})
{seed}
parameter_config.show_parameter_config_page()
"""


def _page(arm="germline", seed="", cached=None):
    """The whole parameter page, with the parameter cache stubbed.

    Stubbed for the reason ``test_app_defaults._render_full_page`` gives: the page consults
    ``~/.mafigate`` before falling back, so an unstubbed render lets whichever arm the
    developer last configured decide what this file is testing. ``discard_stale_cache`` is
    stubbed as well, because an unstubbed render moves the developer's own cache file.

    ``clear_parameters_cache`` is stubbed for a stronger reason than either, learned the
    expensive way: this file drives the *🗑️ Clear Cache* button, and unstubbed that button
    does exactly what it says on the machine running the suite. It deleted this developer's
    ``~/.mafigate/cached_parameters.json`` once, and passed while doing it — then failed on
    the second run, there being nothing left to delete. A test that only works once is
    announcing a side effect.
    """
    pytest.importorskip("streamlit", reason="streamlit not installed")
    from streamlit.testing.v1 import AppTest

    from page_modules import param_store, parameter_config

    info = None if cached is None else {
        "timestamp": "2026-08-15T09:00:00",
        "file_size": 512,
        "app_version": "2.0.0",
        "schema_version": "1",
    }
    with mock.patch.object(
        parameter_config, "load_parameters_from_cache", lambda: cached
    ), mock.patch.object(
        parameter_config, "save_parameters_to_cache", lambda params: True
    ), mock.patch.object(
        parameter_config, "get_cache_info", lambda: info
    ), mock.patch.object(
        parameter_config, "clear_parameters_cache", lambda: True
    ), mock.patch.object(param_store, "discard_stale_cache", lambda: None):
        app = AppTest.from_string(
            _PAGE_SCRIPT.format(app=str(STREAMLIT_APP), arm=arm, seed=seed)
        )
        app.run(timeout=60)
        assert not app.exception, [str(e.value) for e in app.exception]
        yield app


def _run_page(arm="germline", seed="", cached=None):
    """``_page`` as a context manager, so the stubs outlive the clicks made inside it."""
    from contextlib import contextmanager

    return contextmanager(_page)(arm=arm, seed=seed, cached=cached)


def _click(app, label):
    button = [b for b in app.button if b.label == label]
    assert button, f"no button labelled {label!r}: {[b.label for b in app.button]}"
    button[0].click().run(timeout=60)
    assert not app.exception, [str(e.value) for e in app.exception]
    return app


def _arm(app):
    return app.session_state["filter_params"]["sample_type"]


def _drawn(app):
    return " ".join(
        [e.value for e in app.info] + [e.value for e in app.success]
    )


# ---------------------------------------------------------------------------
# The reported symptom, end to end
# ---------------------------------------------------------------------------


def test_a_cross_arm_preset_still_switches_the_arm_and_now_says_so():
    """The switch is kept — and announced.

    Both halves are the decision. All four presets stay offered on both arms, because the
    labels carry the arm and a one-click route from *germline, default* to *somatic,
    stringent* is something a user can mean. What was missing was the app's half of the
    exchange.
    """
    from page_modules import parameter_config

    with _run_page(arm="germline") as app:
        assert _arm(app) == "germline"
        offered = [b.label for b in app.button if b.label.startswith("Load ")]
        assert "Load Broad Somatic" in offered, offered

        _click(app, "Load Broad Somatic")

        assert _arm(app) == "somatic", (
            "the preset no longer moves the arm — if that is deliberate, this file is the "
            "wrong place to decide it, since #133 kept the cross-arm route on purpose"
        )
        drawn = _drawn(app)
        assert "Switched to somatic" in drawn, (
            f"the arm moved with nothing saying so, which is the reported bug: {drawn}"
        )
        assert "Broad Somatic parameters loaded." in drawn, drawn
        assert parameter_config.PARAM_NOTICE_KEY not in app.session_state, (
            "the notice was left in session state, so it was drawn on a frame that did "
            "not survive — or it will be drawn again on the next render"
        )


def test_a_same_arm_preset_says_it_loaded_and_claims_no_switch():
    """The confirmation is restored; the arm clause is only for an arm that moved."""
    with _run_page(arm="germline") as app:
        _click(app, "Load Stringent Germline")

        assert _arm(app) == "germline"
        drawn = _drawn(app)
        assert "Stringent Germline parameters loaded." in drawn, drawn
        assert "Switched to" not in drawn, (
            f"an arm that did not move was announced as a switch: {drawn}"
        )


def test_the_preset_lands_whole_rather_than_re_seeded_from_the_contract():
    """A preset is a complete set for its own arm, so it is adopted as-is.

    The Sample Type control deliberately re-seeds from ``pipeline_params(arm)`` and
    discards edits. Doing that here would discard the very preset the click asked for by
    name — which is checkable, because the Broad presets differ from the contract.

    Keep-lists are compared as sets, and that is not a weakening: the render after the
    adoption rebuilds every multiselect around the new dict and writes each selection back
    in the *control's* option order, so the order in session state after a click is the
    widget's rather than the preset's. What the filter reads is membership — ``isin`` —
    so the set is the claim. A re-seed changes membership, which is what this catches.

    ``keep_pathogenic`` is checked through the page's own reader for the same reason. The
    presets spell pathogenic retention that way and the contract spells it
    ``skip_pathogenic``; the retention checkbox writes the contract's spelling and drops
    the preset's on the render after adoption. That is a deliberate translation predating
    this ticket — what must survive is the setting, not the key.
    """
    from config.presets import SOFT_SOMATIC_PARAMS
    from page_modules.parameter_config import _pathogenic_retention_on

    with _run_page(arm="germline") as app:
        _click(app, "Load Broad Somatic")

        landed = app.session_state["filter_params"]
        assert _pathogenic_retention_on(landed) == _pathogenic_retention_on(
            SOFT_SOMATIC_PARAMS
        ), "the preset's pathogenic-retention setting did not survive adoption"

        for key, value in SOFT_SOMATIC_PARAMS.items():
            if key == "keep_pathogenic":
                continue
            got = landed[key]
            same = (
                sorted(got) == sorted(value)
                if isinstance(value, list) and isinstance(got, list)
                else got == value
            )
            assert same, (
                f"{key} came out as {got!r}, not the preset's {value!r} — the preset was "
                "re-seeded or merged rather than adopted"
            )


# ---------------------------------------------------------------------------
# The sentence itself
# ---------------------------------------------------------------------------


def _notice(switched_to, confirmation="Broad Somatic parameters loaded.", maf=None):
    from unittest.mock import MagicMock

    from page_modules import parameter_config
    from tests.fakes import FakeSessionState

    session = FakeSessionState({
        parameter_config.PARAM_NOTICE_KEY: {
            "confirmation": confirmation,
            "switched_to": switched_to,
        },
        "maf_data": maf,
    })
    with mock.patch.object(parameter_config, "st", MagicMock()) as st_mock:
        st_mock.session_state = session
        parameter_config.show_parameter_notice()
        drawn = [call[0][0] for call in st_mock.info.call_args_list] + [
            call[0][0] for call in st_mock.success.call_args_list
        ]
    return " ".join(drawn), session


def test_the_notice_names_the_open_file_only_when_there_is_one():
    """The parameter page is reachable with nothing loaded.

    *MAFigate now filters this file as somatic* would be a claim about a file that is not
    there, so the subject is dropped rather than the sentence.
    """
    with_file, _ = _notice("somatic", maf=object())
    assert "filters this file as somatic." in with_file, with_file

    without_file, _ = _notice("somatic", maf=None)
    assert "filters as somatic." in without_file, without_file
    assert "this file" not in without_file, without_file


def test_the_notice_does_not_borrow_the_selectbox_sentence():
    """It borrows that route's *vocabulary*, and cannot borrow its promise.

    ``🔄 Switched to {arm}: parameters reset to the settings MAFigate opens with`` is true
    of the Sample Type control, which re-seeds from the contract. It is false of a preset,
    which brings its own parameters — so the two routes share an opening and diverge after
    the colon.
    """
    drawn, _ = _notice("somatic")
    assert "Switched to somatic" in drawn
    assert "opens with" not in drawn, (
        f"the preset route claims it reset to the opening settings, which it did not: "
        f"{drawn}"
    )


def test_the_notice_is_drawn_once():
    """Consumed on the run after the rerun, and not on every run after that."""
    from unittest.mock import MagicMock

    from page_modules import parameter_config

    drawn, session = _notice("somatic")
    assert drawn
    assert parameter_config.PARAM_NOTICE_KEY not in session

    with mock.patch.object(parameter_config, "st", MagicMock()) as st_mock:
        st_mock.session_state = session
        parameter_config.show_parameter_notice()
        assert not st_mock.info.called and not st_mock.success.called, (
            "the notice was drawn a second time"
        )


# ---------------------------------------------------------------------------
# The other two sites
# ---------------------------------------------------------------------------


def _adopt(previous, params):
    """``adopt_parameters`` over a session, returning what it parked."""
    from unittest.mock import MagicMock

    from page_modules import parameter_config
    from tests.fakes import FakeSessionState

    session = FakeSessionState({} if previous is None else {"filter_params": previous})
    with mock.patch.object(parameter_config, "st", MagicMock()) as st_mock:
        st_mock.session_state = session
        parameter_config.adopt_parameters(params, "Parameters loaded.")
    return session[parameter_config.PARAM_NOTICE_KEY]


def test_a_first_adoption_announces_no_switch():
    """There is nothing to have moved away from.

    Found by review. Comparing the incoming arm against a missing one made *every* first
    adoption a switch, which tells a user they were somewhere they have never been. The
    page seeds before any button can be pressed, so this is reached through
    ``_adopt_uploaded_parameters`` rather than from the screen — but the app has no other
    guard that the notice is about a movement.
    """
    assert _adopt(None, {"sample_type": "somatic"})["switched_to"] is None


def test_a_set_that_states_no_arm_is_read_as_somatic():
    """Silence is not neutrality — it is what every other reader calls somatic.

    Also found by review, and it is the one route by which the arm could still change
    without saying so: a cache carrying no ``sample_type`` lands on a germline session,
    the page then reads it as somatic (``MAFigate.initialize_session_state``,
    ``data_loading.validate_required_columns``), and a
    plain "loaded" line would be the only thing said about it.
    """
    assert _adopt({"sample_type": "germline"}, {"min_depth": 30})["switched_to"] == (
        "somatic"
    )
    assert _adopt({"min_depth": 30}, {"sample_type": "somatic"})["switched_to"] is None


def test_an_uploaded_preset_naming_the_other_arm_says_so():
    """A *changed* key is not a dropped one, which is why the report missed this.

    ``migrate_params`` names ``sample_type`` only when the value is not an arm at all. A
    file that legitimately said the other arm moved the user with the report saying
    nothing.
    """
    from unittest.mock import MagicMock

    from page_modules import parameter_config
    from tests.fakes import FakeSessionState

    session = FakeSessionState({"filter_params": {"sample_type": "germline"}})
    with mock.patch.object(parameter_config, "st", MagicMock()) as st_mock:
        st_mock.session_state = session
        parameter_config._adopt_uploaded_parameters(
            {"sample_type": "somatic"}, "my_preset.json"
        )
        assert st_mock.rerun.called

    notice = session[parameter_config.PARAM_NOTICE_KEY]
    assert notice["switched_to"] == "somatic", notice
    assert "my_preset.json" in notice["confirmation"], notice


def test_reloading_the_cache_says_which_arm_it_restored():
    """A cache is auto-saved on every change, so it can hold a previous session's arm.

    Its own ``✅ Parameters reloaded from cache!`` never reached anyone either: it was
    drawn immediately before the rerun that threw the frame away.
    """
    cached = {"sample_type": "somatic", "min_depth": 30}

    with _run_page(arm="germline", cached=cached) as app:
        # The page seeds from the script, not the cache, so the arm is germline until the
        # button is pressed.
        assert _arm(app) == "germline"
        _click(app, "📥 Reload from Cache")

        assert _arm(app) == "somatic"
        drawn = _drawn(app)
        assert "Switched to somatic" in drawn, drawn
        assert "reloaded from your saved cache." in drawn, drawn


def test_saving_and_clearing_the_cache_say_so_where_the_user_can_read_it():
    """Two more confirmations that were drawn into a discarded frame.

    Neither touches the arm, and both are here because the defect underneath this ticket
    was not about the arm: a message written on the way to ``st.rerun`` reaches nobody.
    Found by review, in the same two functions this change was already editing.
    """
    with _run_page(arm="germline", cached={"sample_type": "germline"}) as app:
        _click(app, "💾 Save Current Parameters")
        assert "Parameters saved to your cache." in _drawn(app), _drawn(app)

    with _run_page(arm="germline", cached={"sample_type": "germline"}) as app:
        _click(app, "🗑️ Clear Cache")
        assert "Your saved cache has been cleared." in _drawn(app), _drawn(app)


def test_one_definition_of_how_the_app_says_the_arm_moved():
    """Both routes open with the same words, and there is one place to change them.

    The Sample Type control draws its notice inline — it does not rerun, so its frame is
    the one the user sees — while a replacement parks one. Two mechanisms, deliberately;
    a second *vocabulary* is what must not happen, and a literal copied into both is how
    it would.

    Docstrings are excluded, and that is the distinction the rule rests on everywhere in
    this app: a docstring explaining the words is not the app saying them. Written the
    other way round, this guard failed on the very function whose docstring argues for
    sharing them.
    """
    from page_modules import parameter_config

    source = Path(inspect.getfile(parameter_config)).read_text()
    tree = ast.parse(source)

    prose = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
            first = node.body[0] if node.body else None
            if isinstance(first, ast.Expr) and isinstance(first.value, ast.Constant):
                prose.add(id(first.value))

    holders = set()
    for node in ast.walk(tree):
        if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            continue
        for inner in ast.walk(node):
            if (
                isinstance(inner, ast.Constant)
                and isinstance(inner.value, str)
                and id(inner) not in prose
                and "Switched to" in inner.value
            ):
                holders.add(node.name)

    assert holders == {"switched_to"}, (
        "the words for an arm change are written in more than one place, which is how the "
        f"two routes drift into two vocabularies: {sorted(holders)}"
    )


# ---------------------------------------------------------------------------
# The rule, so a fourth site cannot be added in silence
# ---------------------------------------------------------------------------


def _app_modules():
    """Every module of the app a replacement could hide in.

    The sweep is app-wide rather than one file, because the rule is about the app and not
    about ``parameter_config``: ``data_loading`` and ``MAFigate`` both seed the same dict,
    which a single-file guard could not see. Found by review.
    """
    modules = [STREAMLIT_APP / "MAFigate.py"]
    for package in ("page_modules", "components"):
        modules.extend(sorted((STREAMLIT_APP / package).glob("*.py")))
    return modules


def _wholesale_replacements():
    """Every assignment that replaces ``filter_params`` itself, with where it is.

    Read from the AST rather than by grep, because the app writes single keys constantly
    (``st.session_state.filter_params["min_depth"] = …``) and those are a different act:
    the target there is a subscript *of* the dict, not the dict.
    """
    found = []
    for path in _app_modules():
        tree = ast.parse(path.read_text())
        for node in ast.walk(tree):
            if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                continue
            for inner in ast.walk(node):
                if not isinstance(inner, ast.Assign):
                    continue
                for target in inner.targets:
                    replaces = (
                        isinstance(target, ast.Attribute)
                        and target.attr == "filter_params"
                    ) or (
                        isinstance(target, ast.Subscript)
                        and isinstance(target.slice, ast.Constant)
                        and target.slice.value == "filter_params"
                    )
                    if replaces:
                        found.append((path.name, node.name, inner.lineno))
    return found


def test_only_the_sanctioned_places_replace_the_parameters_wholesale():
    """The rule, held by structure rather than by everyone remembering it.

    Two functions may replace the dict, and each already says what it did:

    * ``adopt_parameters`` — the shared route, which is what announces the arm;
    * ``show_parameter_config_page`` — the session's own seeding, which happens before
      there is a previous arm to move away from, and the Sample Type re-seed, which draws
      its own ``🔄 Switched to`` notice inline because it does not rerun.

    A third function replacing it is another silent path, and that is the bug this ticket
    closed. This guard earned its place immediately: it found ``reset_parameters``, which
    the ticket did not list. That one cannot move the arm — it is handed the arm the page
    is on — but it carried the same defect underneath, a confirmation drawn into the frame
    ``st.rerun`` throws away, so *Reset to Defaults* said nothing either. It goes through
    the helper now.
    """
    #: Every sanctioned writer, with how many writes it may make and why it is allowed —
    #: written out rather than derived, so adding one is a deliberate act with a reason
    #: beside it.
    sanctioned = {
        ("parameter_config.py", "adopt_parameters"): (
            1,
            "the shared route: this is what announces the arm",
        ),
        ("parameter_config.py", "show_parameter_config_page"): (
            5,
            "the cached seed and the contract seed (both before any arm exists to move "
            "away from), the Sample Type re-seed (which announces itself inline, having "
            "no rerun to survive), the `validate_multiselect_params` normalisation "
            "(which repairs the shape of the dict it was given and cannot reach "
            "`sample_type`), and the `complete_params` completion (which fills the keys "
            "the arm's contract names and the dict does not carry, and deliberately skips "
            "`sample_type` for this exact reason — see issue #280)",
        ),
        ("data_loading.py", "show_data_loading_page"): (
            2,
            "seed-if-absent, so there is no previous arm to move away from, and the "
            "`complete_params` completion — sanctioned for the same reason it is on the "
            "parameters page, that it fills the keys the arm's contract names and cannot "
            "reach `sample_type`, so it is the one kind of wholesale write that provably "
            "cannot move the user between arms (issue #289)",
        ),
        ("MAFigate.py", "initialize_session_state"): (
            2,
            "the session's own seeding, from the cache or from the contract, before any "
            "page has run",
        ),
    }

    counted = {}
    for module, function, lineno in _wholesale_replacements():
        counted.setdefault((module, function), []).append(lineno)

    unsanctioned = {site: lines for site, lines in counted.items() if site not in sanctioned}
    assert not unsanctioned, (
        "a parameter set is being adopted somewhere that cannot announce the arm: "
        f"{sorted(unsanctioned)}. Route it through `adopt_parameters`, or add it here "
        "with the reason it cannot move the arm."
    )

    for site, (expected, why) in sanctioned.items():
        lines = counted.get(site, [])
        assert len(lines) == expected, (
            f"{site[0]}:{site[1]} writes `filter_params` {len(lines)} times, not "
            f"{expected}. The sanctioned ones are: {why}. A new one has to say which "
            f"kind it is: lines {sorted(lines)}"
        )
