"""The filter's slot says something about the open file's latest run, or nothing at all.

Everything the page-level slot draws is a statement about one file and one cut of it — the
refusal, the run notes, the chromosome notice. None of them may outlive either. Until issue
#155 two of the three could, by the same mechanism from two directions:

* the **writer** was partial. ``elif notes:`` wrote the stash when a run had something to
  say and left it alone when it did not, so the previous run's notes stayed standing for
  whichever render came next to draw. ``_report_filter_failure`` had the mirror of it — it
  wrote the refusal on every silent failure and cleared it on no success.
* the **shared load tail** cleared neither, so a file replaced with no render of this page
  between the two loads handed its account to its successor.

Both halves are needed and neither subsumes the other, which is the claim
:func:`test_a_refused_file_does_not_inherit_the_last_files_notes` exists to make: a MAF
refused by ``validate_required_columns`` returns *before* the filter runs, so there is no
run to replace anything, and only the tail can clear it.

Why this was not a line in issue #151
-------------------------------------
Because the empty-note run read as the rare branch, and it is the ordinary one. #151 measured
the four shipped presets (two per arm) on the two reference MAFs — 4 runs of 4, each drawing
exactly one note — and the population-frequency report is 97% of every note the filter can
draw. What that configuration is not, is the one the app starts in.

Measured over **158 placeable real MAFs**, each on its own arm under both of that arm's
presets: **116 of 316 runs produce no notes at all — 36.7%**, 72 of the 158 files under at
least one preset and 44 under both. Zero notes is not a missing frequency column, either;
every one of the 44 carries three to five populated frequency columns, and the note fires only
``if removed or exempted``. It correlates with a *small* MAF — median 55 rows against 1,785
for the files that always speak.

The app's opening state is the silent one by construction, which is what
:func:`test_the_apps_opening_state_produces_no_notes_at_all` pins: ``pipeline_params`` ships
``max_freq_population = 1.0``, the set ``MAFigate.py`` falls back to, *Reset to Defaults*
installs and the mismatch notice's arm switch installs — and ``apply_filters`` skips the whole
frequency layer above 1.0, so 310 of those same 316 runs fall silent at that threshold.

Reachability, and why the guards below are shaped as they are
-------------------------------------------------------------
The sequence the ticket asked about — two loads with no render of the data page between them —
does not appear to be reachable through the interface today. What makes it unreachable is a
*pairing*, not a single barrier: every load happens near the top of ``show_data_loading_page``
and the slot is drained near the bottom of the same call, so a load drains in its own run.

Three ``st.rerun`` calls can stand between a stash and that drain, and it is worth naming all
three rather than the obvious one — the page's own jump rerun, which fires only on a jump the
*section* set and so cannot precede the load; ``_reapply_from_results``, which stashes and then
reruns unconditionally; and ``adopt_parameters``, reached from the arm switch's ``then`` hook,
which stashes and reruns from **inside this slot**, ahead of the two drains. In each case the
run that follows renders this same page and drains what was left.

The upload token is not the barrier it first looks like, either: it is stamped only after a
*successful* filter run, so a MAF that loads and then fails to filter is re-read on every
render. That does not open the sequence, because such a re-read ends in a refusal the same run
draws — but it is the reason the argument rests on the drain following the load rather than on
a second load being impossible.

What *is* reachable is a run that stashes and never draws — a page that raises between the two,
which the router reports in place since issue #144, or a widget event landing during the
stashing run — after which the next load draws the last file's account. So this is
precautionary, exactly as issue #149's own clear-before-set is, and reachable by no wider a
route than that one.

So these are unit claims on the three sites rather than a driven sequence: the contract of the
writer, of the load tail, and of the shape that ties them. A behavioural test of a sequence
the UI cannot currently reach would be a test of a fixture, and it would go green for the
wrong reason the moment the page's ordering changed.

A unit module: the pipeline has no session state and no banner slot, so nothing here has a
counterpart in ``bin/``. The last section reads this app's own source, in the manner of the
note-construction tripwire in ``test_filter_notes.py`` — a claim about the *shape* of the
three sites, which no behavioural test can make about a stash that has not been invented yet.
"""

from __future__ import annotations

import ast
from pathlib import Path
from unittest.mock import MagicMock

import pandas as pd
import pytest

from config import presets
from config.pipeline_params import pipeline_params
from filters.notes import INFO, Note
from filters.numeric_columns import UnreadableNumericColumns
from filters.variant_filters import apply_filters
from tests.fakes import FakeSessionState
from utils import read_maf

FIXTURES = Path(__file__).parent / "fixtures" / "parity"

DATA_LOADING_SOURCE = (
    Path(__file__).resolve().parent.parent / "page_modules" / "data_loading.py"
)

#: A note from an earlier run, recognisable in an assertion. Its wording is nobody's copy —
#: what is under test is whether it survives, never what it says.
STALE = Note(INFO, "🌍 a note about the file you were looking at before")


# ---------------------------------------------------------------------------
# The measurement that reprices the branch
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "name, arm",
    [("somatic_reference.maf", "somatic"), ("germline_reference.maf", "germline")],
)
def test_the_apps_opening_state_produces_no_notes_at_all(name, arm):
    """A zero-note run is what the app does on a fresh session, not a corner of one.

    ``pipeline_params(arm)`` is the app's opening state: ``MAFigate.py`` falls back to it when
    no cache is present, the data page seeds it, *Reset to Defaults* installs it and so does
    the mismatch notice's arm switch. It carries ``max_freq_population = 1.0``, which removes
    nothing — and the population-frequency note, the one message issue #151 found the slot
    always drawing, fires only when the filter removed or exempted a row.

    Both files are complete and each is on the arm it was annotated for, so there is nothing
    else to report either. That is the state in which the old ``elif notes:`` neither wrote
    the stash nor cleared it — and it is not a property of these two fixtures: across 158 real
    MAFs, holding everything else and setting only this threshold to 1.0 silences 310 of 316
    runs, because ``apply_filters`` skips the frequency layer entirely above 1.0.
    """
    maf = read_maf(FIXTURES / name)
    params = pipeline_params(arm)

    assert params["max_freq_population"] == 1.0, (
        "the contract no longer opens at 1.0, so this test is measuring something else — "
        "re-derive whether the app's opening state still produces a silent run"
    )

    _, diagnostics = apply_filters(maf, params)

    assert list(diagnostics.notes) == [], (
        f"{name} on {arm} under the app's opening state drew "
        f"{[(n.level, n.text[:60]) for n in diagnostics.notes]}. If this now says something, "
        "the empty-note run has stopped being the default — which changes how often the "
        "clearing this module guards is what stands between two files"
    )


@pytest.mark.parametrize(
    "name, arm, preset",
    [
        ("somatic_reference.maf", "somatic", presets.SOFT_SOMATIC_PARAMS),
        ("somatic_reference.maf", "somatic", presets.CLINICAL_SOMATIC_PARAMS),
        ("germline_reference.maf", "germline", presets.SOFT_GERMLINE_PARAMS),
        ("germline_reference.maf", "germline", presets.CLINICAL_GERMLINE_PARAMS),
    ],
)
def test_a_preset_does_produce_a_note(name, arm, preset):
    """The other half of the measurement, and the reason the first one is surprising.

    Issue #151's "4 of 4" is reproduced here rather than trusted, because the whole point of
    the pair is the *contrast*: change one number in the contract and the opening state starts
    behaving like a preset, at which point the test above passes for a reason that has nothing
    to do with what it claims.

    What the pair does **not** claim is that a preset always speaks. On the real corpus a
    preset run is silent 36.7% of the time; these two files are among the loud ones because
    they are large enough for something to cross the threshold. The two references are the
    right instrument for *this* contrast and the wrong one for that rate, which is why the rate
    is recorded in the module docstring and not asserted here.

    The reference table is re-derived rather than imported from ``test_filter_notes.py``, which
    holds the same four runs. That is not an oversight: ``tests/fakes.py`` exists because
    ``test_app_identity.py`` borrowed a helper *out of* a collected test module and so made that
    module double as a library — the coupling being invisible to ``test_suite_organisation.py``.
    Importing four lines of parametrisation from a sibling would reinstate exactly that.
    """
    maf = read_maf(FIXTURES / name)
    params = dict(preset)
    params.setdefault("sample_type", arm)

    _, diagnostics = apply_filters(maf, params)

    assert diagnostics.notes, (
        f"{name} on {arm} under {params['max_freq_population']} drew nothing, so the "
        "contrast this pair rests on has gone"
    )


# ---------------------------------------------------------------------------
# The writer replaces the slot rather than adding to it
# ---------------------------------------------------------------------------


def _filter(monkeypatch, notes=(), refusal=None, **seed):
    """Run ``apply_filters_to_data`` silently over a stubbed filter, and return the session.

    Stubbed at ``apply_filters`` rather than driven with a real MAF, because what is under
    test is the *bookkeeping* around a run and not the run: the interesting cases are a run
    with no notes and a run that refuses, and reaching either through a real file would make
    each assertion depend on a fixture chosen to produce it. The five collaborators after the
    filter are stubbed for the same reason — they are what the run stores, not what it says.
    """
    from page_modules import data_loading

    state = FakeSessionState(**seed)
    state.setdefault("maf_data", pd.DataFrame({"Hugo_Symbol": ["TP53"]}))
    state.setdefault("filter_params", pipeline_params("somatic"))

    fake_st = MagicMock()
    fake_st.session_state = state
    monkeypatch.setattr(data_loading, "st", fake_st)

    def fake_apply(maf, params):
        if refusal is not None:
            raise refusal
        labelled = maf.copy()
        labelled[data_loading.MAFIGATE_FILTER] = data_loading.PASS
        return labelled, MagicMock(notes=list(notes))

    monkeypatch.setattr(data_loading, "apply_filters", fake_apply)
    monkeypatch.setattr(data_loading, "circle_sources", lambda *a, **k: [])
    monkeypatch.setattr(
        data_loading, "add_clinical_summary_column", lambda frame, sources: frame
    )
    monkeypatch.setattr(data_loading, "params_hash", lambda params: "hash")
    monkeypatch.setattr(data_loading, "describe_run", lambda *a, **k: None)
    monkeypatch.setattr(data_loading, "attribute_report", lambda *a, **k: None)

    produced = data_loading.apply_filters_to_data(show_messages=False)
    return state, produced


def test_a_run_with_notes_stashes_them(monkeypatch):
    """The behaviour that already worked, pinned so the fix cannot be a deletion.

    Mutation: clear unconditionally instead of replacing → the slot goes permanently silent
    and every account of a filled column, an unusable gene list or a frequency cut is lost on
    exactly the path that produces them, which is every load.
    """
    from page_modules.data_loading import _FILTER_NOTES

    fresh = Note(INFO, "🌍 this run's own report")
    state, produced = _filter(monkeypatch, notes=[fresh])

    assert produced is True
    assert state[_FILTER_NOTES] == [fresh]


def test_a_run_with_nothing_to_say_clears_the_last_runs_notes(monkeypatch):
    """**The finding.** Silence replaces the account; it does not leave it standing.

    This is the test that fails on the code before issue #155: ``elif notes:`` skipped both
    the write and the clear, so the seeded note below survived a run that had nothing to do
    with it and was drawn by the next render of the page.

    Mutation: restore ``elif notes:`` → the stale note is still here.
    """
    from page_modules.data_loading import _FILTER_NOTES

    state, produced = _filter(monkeypatch, notes=[], **{_FILTER_NOTES: [STALE]})

    assert produced is True
    assert _FILTER_NOTES not in state, (
        f"a run that produced no notes left {state[_FILTER_NOTES]!r} in the slot, so the "
        "page draws the previous run's account against this one's report"
    )


def test_a_report_clears_the_last_refusal(monkeypatch):
    """The refusal is half of the same account, and a report is the news that it is over.

    Mutation: drop the ``filter_error`` pop from ``_stash_the_runs_account`` → a file that was
    refused and then filtered successfully — which is exactly what switching the arm does for
    a wrong-arm MAF — is reported as refused above the report it produced.
    """
    from page_modules.data_loading import _FILTER_ERROR

    state, produced = _filter(
        monkeypatch, notes=[], **{_FILTER_ERROR: "🛑 the last file could not be filtered"}
    )

    assert produced is True
    assert _FILTER_ERROR not in state


def test_a_refusal_clears_the_notes_it_cannot_be_about(monkeypatch):
    """A run either produced a report to describe or was stopped before it could.

    ``apply_filters`` raises before this run reaches the note stash, so the notes standing in
    the slot belong to some earlier run. Left there they are drawn *underneath* the refusal:
    an account of a cut beside the news that no cut was made.

    Mutation: drop the ``_FILTER_NOTES`` pop from ``_report_filter_failure`` → the stale note
    is drawn under the refusal banner.
    """
    from page_modules.data_loading import _FILTER_ERROR, _FILTER_NOTES

    state, produced = _filter(
        monkeypatch,
        refusal=UnreadableNumericColumns({"tumor_f": {".": 3}}),
        **{_FILTER_NOTES: [STALE]},
    )

    assert produced is False
    assert state[_FILTER_ERROR].startswith("🛑"), (
        "the fixture did not actually refuse — the framing is what tells this case apart "
        "from a crash"
    )
    assert _FILTER_NOTES not in state


# ---------------------------------------------------------------------------
# The shared load tail clears what the last file left
# ---------------------------------------------------------------------------


def _open(monkeypatch, refused=False, **seed):
    """Take a frame through the shared load tail, and return the session it leaves.

    The same seam ``test_chromosome_spelling.py`` drives for the chromosome notice, and for
    the same reason: this is the one door all three load paths funnel through, so a claim made
    here is a claim about every way a file can arrive.

    Its near-twin over there is **deliberately not shared**, which is worth saying because this
    repo hoists duplicated readers into ``tests/fakes.py`` and did exactly that for
    ``note_texts``. The two differ on the one line that matters: this one stubs
    ``normalise_chromosome_spelling`` out, because the re-spelling is irrelevant to a claim
    about stashes and a real one would need a frame shaped for it — while over there that
    function *is* the subject. A shared helper would need a flag deciding whether to disarm the
    thing its other caller exists to test, which is a worse object than two short fixtures.
    """
    from page_modules import data_loading

    state = FakeSessionState(**seed)
    state.setdefault("filter_params", pipeline_params("somatic"))

    fake_st = MagicMock()
    fake_st.session_state = state
    monkeypatch.setattr(data_loading, "st", fake_st)
    monkeypatch.setattr(data_loading, "normalise_chromosome_spelling", lambda data: False)
    monkeypatch.setattr(
        data_loading, "validate_required_columns", lambda data: not refused
    )
    # The filter is stubbed *out*, not stubbed to produce nothing: a run that replaced the
    # stash would answer the question this section is asking, which is what the tail does on
    # its own. The refused case never reaches it at all.
    monkeypatch.setattr(
        data_loading, "apply_filters_to_data", lambda show_messages=True: True
    )

    data_loading._open_the_file_just_read(
        pd.DataFrame({"Hugo_Symbol": ["TP53"]}), "next_sample.maf"
    )
    return state


def test_notes_never_drawn_do_not_survive_the_next_file(monkeypatch):
    """A file replaced with no render of this page between the two loads.

    The chromosome notice two lines above is cleared here for precisely this reason (issue
    #149) and these two were not, though they are drawn by the same render of the same slot.
    #149's comment called that flag "the only stash in this tail that a later render rather
    than a later run consumes" — true of the keys the tail writes, and the reason these were
    passed over.

    Mutation: drop the ``_FILTER_NOTES`` pop from the reset tail → the new file's report is
    introduced by the last file's notes.
    """
    from page_modules.data_loading import _FILTER_NOTES

    state = _open(monkeypatch, **{_FILTER_NOTES: [STALE]})

    assert _FILTER_NOTES not in state


def test_a_refusal_never_drawn_does_not_survive_the_next_file(monkeypatch):
    """The same claim for the other key. Mutation: drop its pop from the reset tail."""
    from page_modules.data_loading import _FILTER_ERROR

    state = _open(monkeypatch, **{_FILTER_ERROR: "🛑 the last file could not be filtered"})

    assert _FILTER_ERROR not in state


def test_a_refused_file_does_not_inherit_the_last_files_notes(monkeypatch):
    """**Why the writer's clearing does not make the tail's redundant.**

    ``validate_required_columns`` returning false ends the load *before* the filter runs, so
    no run happens and there is nothing to replace the slot's contents. If only the writer
    cleared, this is the file that would meet the previous one's notes — drawn beside the
    error saying this MAF could not be opened.

    Mutation: move both pops from the tail into ``_stash_the_runs_account`` alone and this is
    the one test that goes red.
    """
    from page_modules.data_loading import _FILTER_ERROR, _FILTER_NOTES

    state = _open(
        monkeypatch,
        refused=True,
        **{_FILTER_NOTES: [STALE], _FILTER_ERROR: "🛑 an older refusal"},
    )

    assert state["maf_data"] is None, "the fixture did not actually refuse the file"
    assert _FILTER_NOTES not in state
    assert _FILTER_ERROR not in state


# ---------------------------------------------------------------------------
# The shape, so a fourth stash cannot repeat this
# ---------------------------------------------------------------------------


def _function(tree: ast.AST, name: str) -> ast.FunctionDef:
    """The named top-level function of a parsed module."""
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name == name:
            return node
    raise AssertionError(
        f"{name} is no longer a function in data_loading.py, so this guard has stopped "
        "reading what it claims to read"
    )


def _session_keys_taken(function: ast.FunctionDef) -> set[str]:
    """Every ``st.session_state`` key this function *removes*, as source expressions.

    Removal rather than every mention, because taking a key is what identifies a stash: a
    read leaves it standing, and only ``del`` and ``pop`` are how one is consumed or cleared.
    Returned unresolved — the caller turns names into values by importing them, which is the
    lesson issue #54 recorded when the ``GUIDELINE_SOURCES`` check stopped parsing a page's
    source and started importing the object.

    **Both spellings of a delete**, and review is why. Streamlit's ``session_state`` supports
    attribute access, so ``del st.session_state.filter_notes`` and
    ``del st.session_state[_FILTER_NOTES]`` are the same act — and the attribute form is what
    this page used until this change replaced it. A first version of this matched the subscript
    alone, which left the guard blind to the *more idiomatic* of the two: a fourth stash
    drained that way would have satisfied the rule below and left
    :func:`test_the_extractor_finds_the_three_stashes_the_slot_has_today` unmoved, since the
    three it names would still be all it could see. An attribute name is already a literal, so
    it is returned quoted to match the subscript form's source text.
    """
    taken = set()
    for node in ast.walk(function):
        if isinstance(node, ast.Delete):
            for target in node.targets:
                if isinstance(target, ast.Subscript) and _is_session_state(target.value):
                    taken.add(ast.unparse(target.slice))
                if isinstance(target, ast.Attribute) and _is_session_state(target.value):
                    taken.add(repr(target.attr))
        if (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and node.func.attr == "pop"
            and _is_session_state(node.func.value)
            and node.args
        ):
            taken.add(ast.unparse(node.args[0]))
    return taken


def _is_session_state(node: ast.AST) -> bool:
    """Whether an expression is ``st.session_state``."""
    return (
        isinstance(node, ast.Attribute)
        and node.attr == "session_state"
        and isinstance(node.value, ast.Name)
        and node.value.id == "st"
    )


def _resolved(expressions: set[str]) -> dict[str, str]:
    """Each source expression mapped to the key it names, resolved through the module."""
    from page_modules import data_loading

    resolved = {}
    for expression in sorted(expressions):
        if expression.startswith(("'", '"')):
            resolved[expression] = ast.literal_eval(expression)
            continue
        assert hasattr(data_loading, expression), (
            f"the slot takes {expression}, which is not a module-level name of "
            "data_loading — this guard resolves keys by importing them, so a key built at "
            "runtime cannot be checked and must not be introduced without changing this"
        )
        resolved[expression] = getattr(data_loading, expression)
    return resolved


def test_every_stash_the_slot_drains_is_cleared_by_a_change_of_file():
    """The rule, derived from the two functions rather than from a list beside them.

    A stash drawn by ``_show_stashed_banners`` describes the file that was open when it was
    written. If a load can replace that file without the slot having rendered, the tail has to
    clear it — which is the whole of issue #155 and, before it, of issue #149. Both were found
    one key at a time, by a person reading the tail and noticing an omission.

    Derived so the third time is caught by the suite instead: a fourth stash added to the slot
    and not to the tail fails here on the commit that adds it. That is the difference between
    this and #149's comment, which stated the same rule in prose and was read past.
    """
    tree = ast.parse(DATA_LOADING_SOURCE.read_text())

    drained = _resolved(_session_keys_taken(_function(tree, "_show_stashed_banners")))
    cleared = _resolved(_session_keys_taken(_function(tree, "_open_the_file_just_read")))

    assert drained, (
        "no stash was found in the slot at all, so this guard is asserting nothing — the "
        "extractor has stopped matching how the drains are written"
    )

    missed = {
        expression: key
        for expression, key in drained.items()
        if key not in set(cleared.values())
    }
    assert not missed, (
        f"_show_stashed_banners drains {sorted(missed.values())} and "
        "_open_the_file_just_read does not clear it, so a file replaced with no render of "
        "the data page between the two loads hands its account to the next file. Add the pop "
        "to the reset tail beside the others, or say in the tail why this key is a statement "
        "about something other than the open file"
    )


def test_the_extractor_finds_the_three_stashes_the_slot_has_today():
    """Anti-vacuity by name, not only by count.

    A bug in ``_session_keys_taken`` that returned the empty set would make the rule above
    pass on any code at all — the failure mode this suite has hit with a
    ``deepcopy(X) == X`` (issue #77) and with an ``AppTest`` probe (issue #133). Naming the
    three means a fourth arriving is a *deliberate* edit here, which is where its reason gets
    written down.

    Three **direct** drains, which is not everything the slot resolves. ``show_parameter_notice``
    and ``drain_missing_column_reports`` each take something too, through their own module's
    function rather than out of ``session_state`` here — so each owns its own lifecycle and
    neither is a key this file can see or should be clearing. The count is of what the load tail
    is answerable for.

    Naming the set is *not* on its own enough to keep the extractor honest: an extractor that
    went blind to one spelling would leave this set unchanged. That is
    :func:`test_the_extractor_reads_both_spellings_of_a_delete`'s job, and it exists because
    review found exactly that hole here.
    """
    from page_modules.data_loading import (
        _CHROMOSOMES_RESPELLED,
        _FILTER_ERROR,
        _FILTER_NOTES,
    )

    tree = ast.parse(DATA_LOADING_SOURCE.read_text())
    drained = set(
        _resolved(_session_keys_taken(_function(tree, "_show_stashed_banners"))).values()
    )

    assert drained == {_CHROMOSOMES_RESPELLED, _FILTER_ERROR, _FILTER_NOTES}, (
        f"the slot now drains {sorted(drained)}. A new stash here needs a pop in the load "
        "tail and a line in this set; a departed one needs removing from both"
    )


def test_the_rule_can_see_a_stash_the_tail_forgets():
    """The guard above, made to fail on purpose.

    The source is the shape ``data_loading.py`` had before issue #155 — a slot that drains a
    key the tail never clears — kept as a string so that reinstating the defect is not
    required to know the rule works. It drives the extractor itself, not a re-implementation
    of it, for the reason ``test_discarded_frames.py`` gives at the same seam.
    """
    module = ast.parse(
        "import streamlit as st\n"
        "\n"
        "def _open_the_file_just_read(data, name):\n"
        "    st.session_state.pop('chromosomes_respelled', None)\n"
        "\n"
        "def _show_stashed_banners():\n"
        "    if 'filter_notes' in st.session_state:\n"
        "        del st.session_state['filter_notes']\n"
        "    st.session_state.pop('chromosomes_respelled', False)\n"
    )

    drained = _session_keys_taken(_function(module, "_show_stashed_banners"))
    cleared = _session_keys_taken(_function(module, "_open_the_file_just_read"))

    assert drained == {"'filter_notes'", "'chromosomes_respelled'"}
    assert drained - cleared == {"'filter_notes'"}


def test_the_extractor_reads_both_spellings_of_a_delete():
    """``del st.session_state.foo`` counts as much as ``del st.session_state['foo']``.

    Streamlit's ``session_state`` takes either, so the two are one act — and the attribute form
    is the one this page used until issue #155 replaced it, which makes it the likelier spelling
    for whoever adds the fourth stash. The first version of the extractor matched the subscript
    alone: the rule above stayed green on a slot draining a key the tail never cleared, and the
    named-set test stayed green too, because the three keys it names were still all the
    extractor could see. Two guards, one blind spot, no red — found by review, not by a run.

    Both forms are asserted from **one** fixture drained two ways, so a regression in either
    branch is a difference between the two halves of the same set.
    """
    module = ast.parse(
        "import streamlit as st\n"
        "\n"
        "def _show_stashed_banners():\n"
        "    del st.session_state.spelled_as_an_attribute\n"
        "    del st.session_state['spelled_as_a_subscript']\n"
    )

    assert _session_keys_taken(_function(module, "_show_stashed_banners")) == {
        "'spelled_as_an_attribute'",
        "'spelled_as_a_subscript'",
    }
