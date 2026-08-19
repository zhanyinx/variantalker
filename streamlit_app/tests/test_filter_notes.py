"""The three levels of the filter's slot, and the invariants issue #151 settled.

The slot on the data page receives every message the filter produces. Until #151 it drew them
at two levels off one predicate, and the reason that was wrong is measurable rather than
aesthetic: :func:`test_a_complete_maf_on_its_own_arm_never_warns` is the finding that opened
the ticket, and it fails on the code that shipped before it.

What is asserted here, and what is deliberately not
---------------------------------------------------
The **level** of every note, because that is what the change decided. Not the prose — each
sentence is asserted where its wording belongs (``test_absent_columns.py``,
``test_gene_lists.py``, ``test_filter_app_extras.py``), and duplicating it here would make a
rephrasing look like a re-levelling.

The one thing about the text that *is* asserted is the glyph's agreement with its tier, in the
direction that can rot: the two warning tiers own ``❌`` and ``⚠️``, so an info note must not
open with either. Asserted as a prohibition rather than against a list of permitted info
glyphs, because the info tier's glyph marks the *topic* (``🌍``, ``🧬``) and a list of topics is
a hand-maintained table that falls behind the moment a filter gains a note — the exact failure
mode this suite has hit before.
"""

import ast
from pathlib import Path

import pandas as pd
import pytest

from config import presets
from filters.absent_columns import ESCALATION_MARKER, is_escalated, plan_fills
from filters.gene_lists import parse_gene_list
from filters.notes import ERROR, INFO, LEVELS, WARNING, Note
from filters.variant_filters import apply_filters
from utils import read_maf

FIXTURES = Path(__file__).parent / "fixtures" / "parity"

FILTERS_DIR = Path(__file__).resolve().parent.parent / "filters"

#: The reference MAFs, each with the arm it was annotated for and the presets for that arm.
#: Both are complete files — every filter input present — which is what makes them the right
#: instrument for the measurement below.
REFERENCES = {
    "somatic_reference.maf": (
        "somatic",
        (presets.SOFT_SOMATIC_PARAMS, presets.CLINICAL_SOMATIC_PARAMS),
    ),
    "germline_reference.maf": (
        "germline",
        (presets.SOFT_GERMLINE_PARAMS, presets.CLINICAL_GERMLINE_PARAMS),
    ),
}


def _reference_runs():
    """Every reference MAF against every preset for its own arm — four runs."""
    for name, (arm, arm_presets) in REFERENCES.items():
        maf = read_maf(FIXTURES / name)
        for preset in arm_presets:
            params = dict(preset)
            params.setdefault("sample_type", arm)
            yield name, arm, params, apply_filters(maf, params)[1]


# ---------------------------------------------------------------------------
# The finding that opened the ticket
# ---------------------------------------------------------------------------


def test_a_complete_maf_on_its_own_arm_never_warns():
    """A run with nothing wrong with it draws no warning and no error. Only notes.

    **This is the measurement issue #151 turned on** and it fails on the code that preceded
    it: before the third level existed, all four of these runs drew exactly one *yellow* box,
    and it reported a population-frequency filter that had done precisely what was asked. Set
    against issue #136's finding that the escalated warning fires on the file's own arm for 2
    of 173 placeable real MAFs, the slot's ordinary output was a warning about success — so
    yellow had nothing left to mean by the time a user met a real one.

    Both files are complete and each is filtered on the arm it was annotated for, so there is
    no honest warning available to draw. Anything louder than ``INFO`` here is the regression.
    """
    for name, arm, _params, diagnostics in _reference_runs():
        assert diagnostics.notes, (
            f"{name} on {arm} produced no notes at all, so this test is asserting nothing. "
            "The population-frequency report should always be here: it fires whenever the "
            "filter removed or exempted a row, which every shipped preset does on both "
            "references."
        )
        louder = [note for note in diagnostics.notes if note.level != INFO]
        assert not louder, (
            f"{name} filtered on its own arm as {arm} drew {len(louder)} box(es) louder than "
            f"blue: {[(n.level, n.text[:70]) for n in louder]}. This file carries every filter "
            "input and is on the arm it was annotated for, so there is nothing to warn about "
            "— and a warning on a healthy run is what makes the real ones unreadable."
        )


# ---------------------------------------------------------------------------
# The renderer cannot fall behind the levels
# ---------------------------------------------------------------------------


def test_every_level_has_a_renderer():
    """A level added to ``filters.notes`` and not drawn is a ``KeyError`` at render time.

    The renderer is a mapping rather than a chain of ``if``s precisely so this can be checked
    against :data:`~filters.notes.LEVELS` instead of trusting the two to be read side by side.
    Imports Streamlit, unlike everything else in this module, because the mapping is about
    Streamlit calls — there is nothing to check without it.
    """
    from page_modules.data_loading import _NOTE_RENDERERS

    assert set(_NOTE_RENDERERS) == set(LEVELS), (
        "filters.notes.LEVELS and page_modules.data_loading._NOTE_RENDERERS have come apart: "
        f"levels with no renderer {sorted(set(LEVELS) - set(_NOTE_RENDERERS))}, renderers for "
        f"no level {sorted(set(_NOTE_RENDERERS) - set(LEVELS))}. A level with no renderer "
        "raises KeyError on the render that meets it."
    )


def test_the_renderer_resolves_streamlit_at_call_time():
    """Each level is drawn through this module's ``st`` global, not a reference frozen at import.

    Caught in review, and it is the sort of break that leaves the suite **green** while
    disarming it. The renderer began as ``{ERROR: st.error, ...}``, which binds the real
    callables once when the module is imported. Eight sites in this suite stub the page by
    replacing ``data_loading.st``; against a frozen mapping every one of them would be ignored
    and the live Streamlit called instead — silently, because a bare-mode Streamlit call does
    not raise.

    The predicate-based renderer this replaced had the property by accident, resolving
    ``st.warning`` through the global on each call. Asserted so it is kept on purpose: the
    mapping stores method *names*, and only a lookup at call time can see a patched ``st``.
    """
    from unittest.mock import patch

    from page_modules import data_loading

    class Recorder:
        def __init__(self):
            self.calls = []

        def __getattr__(self, name):
            def record(text):
                self.calls.append((name, text))

            return record

    recorder = Recorder()
    with patch.object(data_loading, "st", recorder):
        data_loading._show_filter_notes(
            [Note.error("❌ e"), Note.warning("⚠️ w"), Note.info("🌍 i")]
        )

    assert recorder.calls == [("error", "❌ e"), ("warning", "⚠️ w"), ("info", "🌍 i")], (
        "a stubbed data_loading.st did not receive the notes; the renderer is holding "
        f"references bound at import instead of resolving through the global. Saw: "
        f"{recorder.calls}"
    )


def test_every_note_carries_a_level_the_renderer_knows():
    """Across every note the reference runs can produce, no level outside :data:`LEVELS`."""
    for name, _arm, _params, diagnostics in _reference_runs():
        for note in diagnostics.notes:
            assert note.level in LEVELS, (
                f"{name} produced a note at level {note.level!r}, which is not one of "
                f"{list(LEVELS)}: {note.text[:80]}"
            )


# ---------------------------------------------------------------------------
# The two answers about escalation cannot disagree
# ---------------------------------------------------------------------------


def test_the_escalation_marker_and_the_error_level_agree():
    """``is_escalated`` and ``level == ERROR`` are the same answer, on a MAF that escalates.

    :func:`~filters.absent_columns.is_escalated` survived #151 for the callers that only hold
    the string, and its docstring claims it cannot disagree with the level. That claim is only
    true while the escalated note is the sole builder writing :data:`ESCALATION_MARKER`, so it
    is asserted rather than asserted-in-prose.
    """
    maf = read_maf(FIXTURES / "germline_reference.maf").drop(columns=["InterVar"])
    _, diagnostics = apply_filters(maf, {"sample_type": "germline"})

    escalated = [note for note in diagnostics.notes if is_escalated(note.text)]
    errors = [note for note in diagnostics.notes if note.level == ERROR]

    assert escalated, (
        "dropping InterVar from the germline reference produced no escalated note, so this "
        f"test is vacuous. Notes were: {[(n.level, n.text[:60]) for n in diagnostics.notes]}"
    )
    assert escalated == errors, (
        "is_escalated and the ERROR level disagree about which notes are escalated. "
        f"by marker: {[n.text[:60] for n in escalated]}; by level: "
        f"{[n.text[:60] for n in errors]}. Some builder other than _degraded_note has "
        f"written {ESCALATION_MARKER!r}, or the escalated note is no longer Note.error."
    )


def test_a_fill_is_never_merely_informational():
    """``FillPlan.warnings()`` speaks only when a stand-in was used, so nothing there is INFO.

    A fill means the report was computed from a value the user's file does not contain, which
    is by definition not the report they asked for.
    """
    maf = read_maf(FIXTURES / "germline_reference.maf").drop(columns=["InterVar"])
    plan = plan_fills(maf, "germline")

    notes = plan.warnings()
    assert notes, "dropping InterVar produced no fill notes, so this test is vacuous"
    assert all(note.level in {ERROR, WARNING} for note in notes), (
        "a fill note was drawn at INFO, which says the run did what was asked. It did not — "
        f"a filter input was replaced by a stand-in. Notes: "
        f"{[(n.level, n.text[:60]) for n in notes]}"
    )


# ---------------------------------------------------------------------------
# The question the ticket actually asked
# ---------------------------------------------------------------------------


def _did_not_run_notes():
    """The three notes that mean *a filter you asked for did not apply*.

    Issue #151's body names two of these. There are three, and before #151 they carried two
    different glyphs: the tokeniser's all-invalid note was ``⚠️`` while the other two were
    ``ℹ️``, which is the inconsistency #150 left behind while aligning the other pair.
    """
    complete = read_maf(FIXTURES / "somatic_reference.maf")

    # No frequency column anywhere, so the app's own frequency layer cannot run.
    no_frequency = complete.drop(
        columns=[column for column in complete.columns if "Freq" in column or "gnomAD" in column]
    )
    _, skipped = apply_filters(
        no_frequency, {"sample_type": "somatic", "max_freq_population": 0.01}
    )

    # A gene restriction asked for, on a frame that cannot be restricted by gene.
    _, no_symbol = apply_filters(
        complete.drop(columns=["Hugo_Symbol"]),
        {"sample_type": "somatic", "filter_genes": "TP53, BRCA1"},
    )

    # Every token unusable, so the restriction vanishes before the frame is consulted.
    _, all_invalid = apply_filters(complete, {"sample_type": "somatic", "filter_genes": "123, ---"})

    return {
        "frequency skipped": [n for n in skipped.notes if "frequency columns" in n.text],
        "no Hugo_Symbol": [n for n in no_symbol.notes if "was not applied" in n.text],
        "no usable symbols": [
            n for n in all_invalid.notes if "No gene filter was applied" in n.text
        ],
    }


@pytest.mark.parametrize("which", ["frequency skipped", "no Hugo_Symbol", "no usable symbols"])
def test_a_filter_that_did_not_run_warns(which):
    """All three of them, at WARNING and with ``⚠️`` — the answer #151 settled.

    The report is wider than the one the user asked for, which is what WARNING says. It is not
    escalated, because issue #28's *extra rows are visible, missing rows are not* holds: the
    rows the missing restriction let through are on screen to be judged.

    Parametrised over the three so a fourth cannot be added at a different level in silence —
    and so a failure names *which* of them drifted, rather than reporting that one of three
    did.
    """
    found = _did_not_run_notes()[which]

    assert len(found) == 1, (
        f"expected exactly one {which!r} note, found {len(found)}: "
        f"{[n.text[:70] for n in found]}. The scenario that produces it has changed shape, so "
        "this assertion is no longer about what it says it is about."
    )
    note = found[0]
    assert note.level == WARNING, (
        f"the {which!r} note is drawn at {note.level!r}. A restriction the user asked for did "
        "not apply, so the report is not the one they requested — that is WARNING. INFO would "
        "say the run did what was asked, and issue #28 calls a silently widened report the "
        f"dangerous one to hide: {note.text[:90]}"
    )
    assert note.text.startswith("⚠️"), (
        f"the {which!r} note opens with {note.text[:2]!r} rather than ⚠️. In the warning tiers "
        "the glyph marks severity and must agree with its box — that agreement is what #151 "
        f"bought, and #150's ℹ️ in a yellow box is what it replaced: {note.text[:90]}"
    )


def test_the_three_did_not_run_notes_agree_with_each_other():
    """What #150 was protecting, kept: they are still identical in level and glyph.

    Asserted as one statement about all three rather than left to the per-note cases above,
    because *agreement* is the property being kept — three right answers reached independently
    can still drift apart pairwise, and pairwise is how #150's version came apart.
    :mod:`filters.notes` records what #150 chose and why.
    """
    notes = [found[0] for found in _did_not_run_notes().values()]

    assert len({note.level for note in notes}) == 1, (
        "the notes meaning 'a filter you asked for did not apply' are drawn at different "
        f"levels: {[(n.level, n.text[:50]) for n in notes]}"
    )
    assert len({note.text[:2] for note in notes}) == 1, (
        "the notes meaning 'a filter you asked for did not apply' carry different glyphs: "
        f"{[(n.text[:2], n.text[:50]) for n in notes]}. This is the state #150 left and #151 "
        "repaired — one of them has drifted back."
    )


# ---------------------------------------------------------------------------
# The gene list's own three, where #151 drew a line inside the tokeniser
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "raw, expected, why",
    [
        (
            "Hugo_Symbol\nTP53\nBRCA1",
            INFO,
            "a heading was dropped and the restriction the user asked for still applied",
        ),
        (
            "TP53, 123, ---",
            WARNING,
            "tokens the user typed were discarded, so what they believe they restricted to "
            "is not what was filtered on",
        ),
        (
            "123, ---",
            WARNING,
            "no restriction applied at all and the report is wider than requested",
        ),
    ],
)
def test_the_tokeniser_levels_its_own_messages(raw, expected, why):
    """The partial rejection warns, and that was the dev's call rather than an obvious one.

    An unusable token could never have matched a gene, so on the frame's own terms nothing
    changed and this could have been an INFO note about a run that did what was asked. It
    warns because the token is evidence of *intent*: a mistyped ``TP-53`` is a gene the
    clinician believes they restricted to and did not, and nothing in the table shows that.
    """
    notes = parse_gene_list(raw).messages()

    assert notes, f"{raw!r} produced no messages, so this case asserts nothing"
    assert all(note.level == expected for note in notes), (
        f"{raw!r} should be {expected} — {why} — but produced "
        f"{[(n.level, n.text[:60]) for n in notes]}"
    )


def test_a_requested_gene_absent_from_the_maf_is_informational():
    """The restriction applied; a symbol simply has no rows here. That is not a warning.

    The contrast with the tokeniser's rejections next door is the whole distinction: there the
    symbol was never filtered on, here it was and matched nothing.
    """
    maf = read_maf(FIXTURES / "somatic_reference.maf")
    present = str(maf["Hugo_Symbol"].dropna().iloc[0])

    _, diagnostics = apply_filters(
        maf,
        {"sample_type": "somatic", "filter_genes": f"{present}, ZZZNOTAGENE1"},
    )
    absent = [note for note in diagnostics.notes if "not present in this MAF" in note.text]

    assert len(absent) == 1, f"expected one 'not present' note, got {[n.text[:70] for n in absent]}"
    assert absent[0].level == INFO, (
        f"a requested gene with no variants in this file was drawn at {absent[0].level!r}. "
        "The gene restriction the user asked for did apply, so the report is the one they "
        f"asked for: {absent[0].text[:90]}"
    )


def test_the_frequency_report_is_informational():
    """The note that made the case for a third level existing at all."""
    for name, _arm, params, diagnostics in _reference_runs():
        reports = [
            note
            for note in diagnostics.notes
            if "variant(s) above" in note.text and "were removed" in note.text
        ]
        assert len(reports) == 1, (
            f"{name} drew {len(reports)} population-frequency reports at "
            f"{params['max_freq_population']}, expected one"
        )
        assert reports[0].level == INFO, (
            f"{name}: the population-frequency filter read the columns it needed and removed "
            f"what it was asked to remove, and the note saying so is drawn at "
            f"{reports[0].level!r}. That is a warning about success — the thing #151 measured."
        )


# ---------------------------------------------------------------------------
# The glyph invariant, and the tripwire on new notes
# ---------------------------------------------------------------------------


def _every_note_the_filter_can_draw():
    """As many distinct notes as the scenarios here can provoke, with their levels.

    Not a claim of exhaustiveness — :func:`test_no_note_is_built_without_a_case_here` is what
    holds that line, by counting construction sites in the source.
    """
    complete = read_maf(FIXTURES / "somatic_reference.maf")
    germline = read_maf(FIXTURES / "germline_reference.maf")
    notes: list[Note] = []

    for found in _did_not_run_notes().values():
        notes.extend(found)

    # The fill pair: one escalated column and one ordinary one, on the germline arm.
    _, filled = apply_filters(
        germline.drop(columns=["InterVar", "RENOVO_Class"]), {"sample_type": "germline"}
    )
    notes.extend(filled.notes)

    # The tokeniser's benign case, and a requested gene that matches nothing.
    present = str(complete["Hugo_Symbol"].dropna().iloc[0])
    _, gene = apply_filters(
        complete,
        {
            "sample_type": "somatic",
            "filter_genes": f"Hugo_Symbol\n{present}\nZZZNOTAGENE1",
        },
    )
    notes.extend(gene.notes)

    return notes


def test_an_info_note_never_wears_a_severity_glyph():
    """``❌`` and ``⚠️`` belong to the two warning tiers; the info tier marks topic instead.

    A prohibition rather than a list of permitted info glyphs. The info tier's glyph says what
    the note is *about* — ``🌍`` for population frequency, ``🧬`` for the gene list — and a
    committed list of topics is a table that stops applying the moment a filter gains a note.
    """
    for note in _every_note_the_filter_can_draw():
        if note.level != INFO:
            continue
        assert not note.text.startswith(("❌", "⚠️")), (
            f"an info note opens with a severity glyph: {note.text[:90]}. The box is already "
            "blue, so the glyph says 'warning' while the box says 'note' — which is the "
            "mismatch #151 existed to remove, re-created in the other direction."
        )


def test_a_warning_or_error_glyph_agrees_with_its_box():
    """The agreement #151 bought, across every note the scenarios above can provoke."""
    expected = {ERROR: "❌", WARNING: "⚠️"}
    for note in _every_note_the_filter_can_draw():
        if note.level not in expected:
            continue
        assert note.text.startswith(expected[note.level]), (
            f"a {note.level} note opens with {note.text[:2]!r}, not "
            f"{expected[note.level]!r}: {note.text[:90]}"
        )


#: Where every note in ``filters/`` is built, and at which level — the one hand-maintained
#: table in this module, and the one line here whose job is to *break*.
#:
#: Per level and per file rather than a bare total, which review caught as the weaker shape: a
#: single number cannot see a ``Note.info`` swapped for a ``Note.warning``, nor one note added
#: while another is deleted, and re-levelling by accident is the whole subject of issue #151.
#: There is no derivation available that would replace this — a level is a judgement about what
#: a sentence means, so the point is to force the judgement, not to compute it.
NOTES_BY_LEVEL = {
    # The fill pair.
    ("absent_columns.py", ERROR): 1,
    ("absent_columns.py", WARNING): 1,
    # The tokeniser: a dropped heading, and its two rejections.
    ("gene_lists.py", INFO): 1,
    ("gene_lists.py", WARNING): 2,
    # The frequency skip and the absent Hugo_Symbol; the frequency report and a gene with no rows.
    ("variant_filters.py", WARNING): 2,
    ("variant_filters.py", INFO): 2,
}


def _note_construction_sites():
    """Every ``Note.<level>(...)`` call under ``filters/``, by file and line.

    Parsed rather than grepped, for the reason this suite prefers an AST everywhere: a call
    inside a docstring or a comment is not a construction site, and a grep cannot tell.

    ``rglob``, so a note built in a future subpackage of ``filters/`` is still counted — the
    directory is flat today and a guard that quietly stops covering a new subdirectory is the
    failure this suite has hit before. Keyed on the literal name ``Note``, which is how every
    module here imports it; an aliased import would slip past, and the honest limit is written
    down rather than guarded against, because the alias is not a thing anyone does here.
    """
    sites = []
    for path in sorted(FILTERS_DIR.rglob("*.py")):
        tree = ast.parse(path.read_text())
        for node in ast.walk(tree):
            if (
                isinstance(node, ast.Call)
                and isinstance(node.func, ast.Attribute)
                and node.func.attr in LEVELS
                and isinstance(node.func.value, ast.Name)
                and node.func.value.id == "Note"
            ):
                sites.append((path.name, node.lineno, node.func.attr))
    return sites


def test_every_note_in_filters_is_built_at_the_level_recorded_here():
    """A note added, removed, or re-levelled in ``filters/`` fails this until it is decided here.

    Issue #151 exists because notes had accumulated in a two-level slot over several tickets,
    each addition defensible on its own, and nothing anywhere made the accumulation a decision.
    This is the thing that makes it one.
    """
    counted: dict[tuple[str, str], int] = {}
    for name, _line, level in _note_construction_sites():
        counted[(name, level)] = counted.get((name, level), 0) + 1

    assert counted == NOTES_BY_LEVEL, (
        f"the notes built in filters/ no longer match NOTES_BY_LEVEL.\n"
        f"  found:    {dict(sorted(counted.items()))}\n"
        f"  recorded: {dict(sorted(NOTES_BY_LEVEL.items()))}\n"
        f"  sites:    {_note_construction_sites()}\n"
        "If you added a note, decide its level against the three in filters/notes.py — is the "
        "report missing rows invisibly (ERROR), not the one the user asked for (WARNING), or is "
        "this an account of a run that did what was asked (INFO)? — then add a case for it "
        "above and record it here. If you re-levelled one, this is the decision being asked for."
    )


def test_the_levels_are_ordered_loudest_first():
    """``LEVELS`` is documented as descending severity, and a reader relies on that."""
    assert LEVELS == (
        ERROR,
        WARNING,
        INFO,
    ), f"LEVELS reads {LEVELS}, which is no longer loudest-first as its comment claims"


def test_a_note_survives_the_stash_as_a_plain_tuple():
    """``Note`` is a NamedTuple so the silent-load stash can leave it in session state.

    The load paths filter with ``show_messages=False`` and stash, and the page draws the notes
    on the next rerender — so whatever is stashed has to survive being put in a dict and taken
    out again with no custom handling anywhere.
    """
    note = Note.warning("⚠️ something")
    assert isinstance(note, tuple)
    assert tuple(note) == (WARNING, "⚠️ something")
    assert Note(*tuple(note)) == note


def test_the_reference_frames_are_actually_complete():
    """The premise of the headline test: both references carry every filter input.

    Without this, ``test_a_complete_maf_on_its_own_arm_never_warns`` could pass on a file that
    warns about nothing because it was never asked to.
    """
    for name, (arm, _presets) in REFERENCES.items():
        maf = read_maf(FIXTURES / name)
        plan = plan_fills(maf, arm)
        assert plan.warnings() == (), (
            f"{name} is missing a filter input for the {arm} arm, so it is not the complete "
            f"file the measurement assumes: {[n.text[:70] for n in plan.warnings()]}"
        )
        assert isinstance(maf, pd.DataFrame)
