"""The neutral fill for an absent filter-input column, and the escalation that guards it.

Issue #39. Needs neither ``bin/`` nor the clinical drive, like ``test_filter_entry_point.py``
and ``test_numeric_columns.py`` — this is app-owned behaviour with no pipeline counterpart to
compare against, so it *cannot* be asserted by the parity harness. The pipeline raises
``KeyError`` on every file this module is about.

Four things are asserted here, and they are not equally strong:

1. **The derived contract, against the vendored code as oracle.** For each arm and each
   column, dropping it and calling the vendored functions raises ``KeyError`` if and only if
   the column is in :data:`~filters.absent_columns.REQUIRED_INPUTS`. That is a biconditional
   against the thing being modelled, not a restatement of the list — and it is what makes
   "the app fills exactly the columns the pipeline would have died on" checkable.

2. **The escalated trigger set, exactly.** Including the two boundary members that read as
   bugs: ``RENOVO_Class``, excluded despite being the single largest source loss, and
   ``CIViC_Evidence_Level``, in the pathogenic set and still never escalated. Both are named
   in tests of their own, because a set assertion alone tells the next reader *what* is true
   and not *why*, and the why is the only thing stopping someone "fixing" it.

3. **Containment.** No fill value reaches the returned frame, per column and per arm. Nothing
   else can see this: the parity harness compares verdicts, and a filled ``tumor_f`` of 1.0 in
   the grid produces a perfectly plausible verdict while telling the user their variant had
   100% VAF.

4. **The direction of each fill**, measured on the fixtures rather than asserted from the
   docstring — a guideline fill never grows the report, a quality-gate fill never shrinks it.
   That is the property the warnings are written around, and it is the one issue #20 found
   inverts everyone's intuition.
"""

from __future__ import annotations

import sys
import tracemalloc
from pathlib import Path

import pandas as pd
import pytest

STREAMLIT_APP = Path(__file__).resolve().parent.parent
if str(STREAMLIT_APP) not in sys.path:
    sys.path.insert(0, str(STREAMLIT_APP))

from config.pipeline_params import pipeline_params  # noqa: E402
from filters.absent_columns import (  # noqa: E402
    ADDS,
    ARM_FILTERS,
    NEUTRAL_FILLS,
    PATHOGENIC_INPUTS,
    REMOVES,
    REQUIRED_INPUTS,
    FillPlan,
    derive_required_inputs,
    is_escalated,
    plan_fills,
)
from filters.variant_filters import (  # noqa: E402
    MAFIGATE_FILTER,
    PASS,
    apply_filters,
)
from vendor.pipeline_filters import (  # noqa: E402
    common_filters,
    germline_filters,
    somatic_filters,
)
from tests.fakes import note_texts  # noqa: E402
from vendor.pipeline_utils import read_maf  # noqa: E402

FIXTURE_DIR = STREAMLIT_APP / "tests" / "fixtures" / "parity"

#: One fixture per arm, both real reference extracts. The synthetic fixtures would do for the
#: mechanics, but the *directions* asserted below are claims about real annotation densities,
#: and a hand-built frame can be made to show any direction at all.
FIXTURES = {
    "somatic": "somatic_reference.maf",
    "germline": "germline_reference.maf",
}

ARMS = list(FIXTURES)


def _frame(arm: str) -> pd.DataFrame:
    return read_maf(str(FIXTURE_DIR / FIXTURES[arm]))


def _params(arm: str, **overrides) -> dict:
    """The contract, from the module guarded against ``nextflow.config``. Never a copy."""
    params = pipeline_params(arm)
    params.update(overrides)
    return params


def _passed(labelled: pd.DataFrame) -> int:
    return int((labelled[MAFIGATE_FILTER] == PASS).sum())


def _arm_columns():
    """One parameter per (arm, required column): fourteen."""
    return [
        pytest.param(arm, column, id=f"{arm}-{column}")
        for arm in ARMS
        for column in REQUIRED_INPUTS[arm]
    ]


# ---------------------------------------------------------------------------
# 1. The derived contract, against the vendored code
# ---------------------------------------------------------------------------


def _call_vendored(frame: pd.DataFrame, arm: str) -> None:
    """Both vendored calls the app makes for ``arm``, at the contract's parameters.

    The gene path is ``"null"`` — no restriction — because that is the configuration in which
    ``Hugo_Symbol`` is *not* read, and the whole point of the every-path rule is that a column
    read only under a condition the app controls is not one the app has to fill.
    """
    params = pipeline_params(arm)
    common_filters(frame, params["min_depth"], params["filter_variant_classification"])
    if arm == "somatic":
        somatic_filters(
            frame,
            vaf=params["vaf_threshold"],
            somatic_genes="null",
            cancervar_keep=params["filter_cancervar"],
            civic_keep=params["filter_civic"],
            escat_keep=params["filter_escat"],
            clinvar_keep=params["filter_clinvar"],
            skip_civic=False,
            skip_pathogenic=False,
        )
    else:
        germline_filters(
            frame,
            vaf=params["vaf_threshold"],
            germline_genes="null",
            intervar_keep=params["filter_intervar"],
            renovo_keep=params["filter_renovo"],
            clinvar_keep=params["filter_clinvar"],
            skip_pathogenic=False,
        )


@pytest.mark.parametrize("arm", ARMS)
def test_required_inputs_are_exactly_the_columns_the_vendored_code_dies_without(arm):
    """The biconditional: derived membership *is* "the vendored call raises without it".

    Every column of a real MAF is dropped in turn — 427 of them on the somatic arm — and the
    vendored functions are called on what is left. This is the assertion that makes the AST
    derivation trustworthy: it never consults the parser's answer to decide what to try, so a
    derivation that gained a column, lost one, or invented one fails here rather than shifting
    the app's behaviour quietly.
    """
    frame = _frame(arm)
    required = set(REQUIRED_INPUTS[arm])

    raises = set()
    for column in frame.columns:
        try:
            _call_vendored(frame.drop(columns=[column]), arm)
        except KeyError:
            raises.add(column)

    assert raises == required, (
        "the columns the vendored filters cannot do without have moved away from the derived "
        f"contract: only-in-code {sorted(raises - required)}, "
        f"only-in-derivation {sorted(required - raises)}"
    )


def test_the_every_path_rule_is_what_is_implemented():
    """A branch-local read is not required; a read in both branches, or in a test, is.

    Against a synthetic source rather than the vendored one, because the rule has to be
    checkable on shapes the pipeline does not currently contain — otherwise the test is a
    second copy of today's answer and says nothing about what the derivation would do to a
    pipeline change.
    """
    source = '''
import pandas as pd

def common_filters(maf: pd.DataFrame, coverage, excluded):
    return maf["always"] > coverage

def somatic_filters(maf: pd.DataFrame, flag):
    if maf["in_the_test"].any():
        both = maf["in_both_branches"]
        only_here = maf["only_in_the_if"]
    else:
        both = maf["in_both_branches"]
    if flag:
        guarded = maf["only_in_a_bare_if"]
    return both

def germline_filters(maf: pd.DataFrame, column):
    computed = maf[column]
    return maf["germline_always"]
'''
    derived = derive_required_inputs(source)

    assert derived["somatic"] == ("always", "in_both_branches", "in_the_test")
    assert derived["germline"] == ("always", "germline_always")


def test_a_renamed_column_moves_the_derivation_with_it():
    """The derivation follows the source, rather than remembering an answer.

    The mutation is applied to the real vendored source, so this is the drift check the
    hand-written list could not have: rename the column in the pipeline and the app fills the
    new name, without anyone editing a list here.
    """
    from filters.numeric_columns import vendored_source

    mutated = vendored_source().replace('maf["ESCAT"]', 'maf["ESCAT_V2"]')
    derived = derive_required_inputs(mutated)

    assert "ESCAT_V2" in derived["somatic"]
    assert "ESCAT" not in derived["somatic"]


def test_a_computed_column_key_is_not_claimed():
    """``maf[column]`` is not readable from the source, so it yields nothing.

    The limit is on the record deliberately. Under-collecting means the app does not fill a
    column and the vendored code raises — today's behaviour — while a guess would mean filling
    a column nobody named.
    """
    source = '''
import pandas as pd

def common_filters(maf: pd.DataFrame, column):
    return maf[column] > 0

def somatic_filters(maf: pd.DataFrame, column):
    return maf[column]

def germline_filters(maf: pd.DataFrame, column):
    return maf[column]
'''
    assert derive_required_inputs(source) == {"somatic": (), "germline": ()}


def test_every_required_column_has_a_neutral_fill():
    """The derived set and the hand-written values are held together, not assumed together.

    ``filters/absent_columns`` raises at import if they part company; this states the property
    where a reader looking for it will find it.
    """
    required = {column for columns in REQUIRED_INPUTS.values() for column in columns}
    assert required <= set(NEUTRAL_FILLS)


def test_no_fill_is_recorded_for_a_column_that_never_needs_one():
    """And the other direction, which the import-time check cannot see.

    A stale entry is not harmless: it reads as a claim that the pipeline once needed the column
    and is the kind of thing that gets copied into the next list.
    """
    required = {column for columns in REQUIRED_INPUTS.values() for column in columns}
    assert set(NEUTRAL_FILLS) == required


def test_the_help_pages_schema_categories_have_not_drifted_from_the_derivation():
    """``config/columns.py``'s status answers must describe the MAF the filter reads.

    There was a **mirror** here until issue #127: a hand-written ``REQUIRED_COLUMNS["clinical"]``
    listing the five annotations the filter reads, which this test held against the derivation
    because the module was said to be unable to import it. It can — nothing in
    ``filters/absent_columns.py`` imports ``config/``, and ``_classify_column_source`` now takes
    ``REQUIRED_INPUTS`` function-locally exactly as ``pipeline_output_columns`` takes
    ``compute_keep`` — so the mirror is gone and there is no drift left to check on that half.

    What is still a decision, and therefore still worth a guard, is the **token** a derived
    filter input gets. The failure mode is unchanged and is the one the whole project is about:
    the help page going on calling a column optional after the pipeline started requiring it,
    with nothing saying so. It is now asserted one layer earlier than the rendered sentence, so
    it holds for every surface that reads the classifier rather than only for the table;
    ``test_no_filter_input_is_called_optional`` guards the rendered words.
    """
    from config.columns import OPTIONAL_COLUMNS, REQUIRED_COLUMNS, _classify_column_source

    derived = {column for columns in REQUIRED_INPUTS.values() for column in columns}
    assert derived, "the derivation found no filter inputs, so this guard checks nothing"

    # `Variant_Classification` is the one column that is both core and a filter input, and it is
    # named here rather than carved out by a membership test. Writing the carve-out as
    # ``"core" if column in REQUIRED_COLUMNS["core"] else "filled"`` reads the same list the
    # branch under test reads, so it agreed with the classifier by construction — review caught
    # it. Named, a *second* such column fails this test, which is the right outcome: whether the
    # refusal or the fill is the truer answer for it is a decision, not a lookup.
    both = sorted(derived & set(REQUIRED_COLUMNS["core"]))
    assert both == ["Variant_Classification"], (
        f"{both} are both core and filter inputs. Only Variant_Classification was, and the "
        "reference table answers it as the refusal because a file missing it never reaches the "
        "filter — decide the same question for the newcomer and name it here"
    )
    assert _classify_column_source("Variant_Classification") == "core"

    for column in sorted(derived - set(REQUIRED_COLUMNS["core"])):
        assert _classify_column_source(column) == "filled", (
            f"the filter reads {column} on every path, and the column reference classifies it "
            f"as {_classify_column_source(column)!r} rather than 'filled' — so the help "
            "pages describe a different MAF from the one the filter reads"
        )

    civic = {column for columns in PATHOGENIC_INPUTS.values() for column in columns} - derived
    assert set(OPTIONAL_COLUMNS["civic"]) == civic, (
        "the genuinely-optional category no longer matches the filter inputs the pipeline "
        "guards for itself"
    )


def test_the_arm_map_names_functions_the_vendored_module_defines():
    """The one thing here that *is* written down, guarded rather than trusted."""
    import vendor.pipeline_filters as vendored

    for names in ARM_FILTERS.values():
        for name in names:
            assert callable(getattr(vendored, name, None)), name


# ---------------------------------------------------------------------------
# 2. The escalated trigger set, exactly
# ---------------------------------------------------------------------------


def test_pathogenic_inputs_are_exactly_the_four_columns_that_feed_retention():
    """Asserted as an exact set, per arm, and by name.

    An assertion on the *size* would pass a substitution, which is the failure that put issue
    #35 in this repo's history.
    """
    assert PATHOGENIC_INPUTS == {
        "somatic": ("CIViC_Evidence_Level", "CancerVar", "ClinVar_VCF_CLNSIG"),
        "germline": ("ClinVar_VCF_CLNSIG", "InterVar"),
    }

    union = {column for columns in PATHOGENIC_INPUTS.values() for column in columns}
    assert union == {
        "CancerVar",
        "ClinVar_VCF_CLNSIG",
        "CIViC_Evidence_Level",
        "InterVar",
    }


@pytest.mark.parametrize(
    "arm,expected",
    [
        ("somatic", ("CancerVar", "ClinVar_VCF_CLNSIG")),
        ("germline", ("ClinVar_VCF_CLNSIG", "InterVar")),
    ],
)
def test_the_escalated_trigger_set_is_exactly_this(arm, expected):
    """Which absent columns actually escalate — the pathogenic set that can also be filled."""
    frame = _frame(arm)
    triggered = tuple(
        column
        for column in REQUIRED_INPUTS[arm]
        if plan_fills(frame.drop(columns=[column]), arm).degraded
    )
    assert triggered == expected


def test_renovo_is_not_escalated_although_it_is_the_largest_single_source_loss():
    """The counter-intuition, named so it does not get "fixed".

    Filling ``RENOVO_Class`` costs **−322 germline rows (−39%)** against the 818-row baseline —
    more than any other single fill on either arm, and more than ``InterVar``'s −147. It is
    still not escalated, and the reason is not severity but *kind*: germline pathogenic
    retention is ``InterVar | ClinVar`` and has no RENOVO clause, so filling RENOVO drops one
    term from the guideline OR and leaves the pipeline's unconditional safety net intact
    (``filter_patho`` 496 → 496). ``InterVar`` takes the same net to **90**.

    The escalated warning means *the safety net is gone*. Escalating the largest loss instead
    of the ones that disarm the net would make it mean "a lot of rows moved", which is a
    different and much less actionable statement — and would leave nothing that says the one
    thing the user cannot recover by re-reading their parameters.
    """
    assert "RENOVO_Class" in REQUIRED_INPUTS["germline"]
    assert "RENOVO_Class" not in PATHOGENIC_INPUTS["germline"]
    assert not plan_fills(_frame("germline").drop(columns=["RENOVO_Class"]), "germline").degraded


@pytest.mark.parametrize(
    "arm,column",
    [("somatic", "CancerVar"), ("germline", "InterVar")],
)
def test_nothing_is_escalated_when_the_user_turned_pathogenic_retention_off(arm, column):
    """``skip_pathogenic`` empties the rescue mask, so a fill cannot degrade it.

    The escalated warning says *variants that would have been retained as pathogenic are absent
    from this report*. Under ``skip_pathogenic`` no variant would have been retained — the user
    asked for that — so escalating would name a loss that did not happen, in the app's loudest
    voice, on a setting they chose deliberately. The column is still filled and still reported;
    it just gets the ordinary dropout line, which is the true statement about what happened.

    The fill is not free even here, which is why the warning does not disappear: the guideline
    term is still gone.
    """
    frame = _frame(arm).drop(columns=[column])
    diagnostics = apply_filters(frame, _params(arm, skip_pathogenic=True))[1]

    assert diagnostics.filled_columns == (column,)
    assert diagnostics.degraded_columns == ()
    assert diagnostics.notes
    assert not any(is_escalated(text) for text in note_texts(diagnostics))


def test_civic_is_the_one_genuinely_optional_filter_input():
    """Read by pathogenic retention, never filled, and therefore never escalated.

    ``CIViC_Evidence_Level`` is the single filter input the pipeline itself guards, with
    ``check_civic_column_exists``. So an absent CIViC column is *real pipeline behaviour*: the
    pipeline skips the clause and reports, and the app does the same by doing nothing at all.
    That run stays **on parity**, unlike every other absent filter input, which is why it must
    not pick up an off-parity warning on the way through.

    It is also the only one that actually happens: absent in 100 of 100 reference files, where
    every other filter input is present in 100 of 100.
    """
    assert "CIViC_Evidence_Level" in PATHOGENIC_INPUTS["somatic"]
    assert "CIViC_Evidence_Level" not in REQUIRED_INPUTS["somatic"]
    assert "CIViC_Evidence_Level" not in NEUTRAL_FILLS

    frame = _frame("somatic")
    assert "CIViC_Evidence_Level" not in frame.columns

    plan = plan_fills(frame, "somatic")
    assert plan.filled == ()
    assert plan.warnings() == ()


# ---------------------------------------------------------------------------
# 3. Containment — no fill value reaches the returned frame
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("arm,column", _arm_columns())
def test_the_absent_column_is_still_absent_from_the_returned_frame(arm, column):
    """The frame the caller gets back is the frame it passed in, plus two verdict columns.

    Stated per column because the failure it guards against is per column: a filled
    ``tumor_f`` reaching the grid reads as a genuine 100% VAF, and a filled
    ``Variant_Classification`` reads as a real annotation. Both would be the app inventing
    data and presenting it as the user's own.
    """
    frame = _frame(arm).drop(columns=[column])
    labelled, diagnostics = apply_filters(frame, _params(arm))

    assert column not in labelled.columns
    assert diagnostics.filled_columns == (column,)


@pytest.mark.parametrize("arm,column", _arm_columns())
def test_no_fill_value_appears_anywhere_in_the_returned_frame(arm, column):
    """Not merely "the column is absent" — the *value* is nowhere in the frame either.

    A stricter reading of the same rule, and the one that survives a refactor: if some future
    caller reinstates the column under another name, or writes the fill into a neighbouring
    column, the column-name assertion above passes and this one does not.

    The sentinel and ``inf`` are values no MAF carries, so their absence is the whole check.
    ``tumor_f``'s 1.0 is a value a real VAF *could* take, so the check is that the returned
    frame holds no more 1.0s than the input did.
    """
    original = _frame(arm)
    frame = original.drop(columns=[column])
    labelled, _ = apply_filters(frame, _params(arm))

    fill = NEUTRAL_FILLS[column].value
    if isinstance(fill, float) and pd.isna(fill):
        # NaN cannot be told apart from a real missing value, so counting it proves nothing.
        # The parameterisation still has to assert *something* rather than skip: what is
        # checkable here is that the frame did not gain the column back under any name.
        assert column not in labelled.columns
        assert len(labelled.columns) == len(frame.columns) + 2  # the two verdict columns
        return

    before = sum(int((frame[name] == fill).sum()) for name in frame.columns)
    after = sum(int((labelled[name] == fill).sum()) for name in labelled.columns)
    assert after == before


@pytest.mark.parametrize("arm", ARMS)
def test_the_verdict_is_unchanged_by_a_column_the_filter_does_not_read(arm):
    """A control: dropping a non-filter-input column changes nothing and warns about nothing.

    Without this, every assertion above would pass on an implementation that filled and warned
    about *any* missing column — which would bury the seven that matter under four hundred that
    do not.
    """
    frame = _frame(arm)
    spare = next(
        column
        for column in frame.columns
        if column not in REQUIRED_INPUTS[arm] and column != "Hugo_Symbol"
    )

    baseline = apply_filters(frame, _params(arm))[1]
    reduced = apply_filters(frame.drop(columns=[spare]), _params(arm))[1]

    assert reduced.cells() == baseline.cells()
    assert reduced.filled_columns == ()
    assert reduced.notes == baseline.notes


# ---------------------------------------------------------------------------
# 4. The direction of each fill
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("arm,column", _arm_columns())
def test_each_fill_moves_the_report_the_way_it_says_it_does(arm, column):
    """Guideline fills never grow the report; quality-gate fills never shrink it.

    The property the warnings are built on, measured rather than asserted from the table in the
    module docstring. It is also the one that inverts intuition: the *removing* family is the
    dangerous one, because a row that is not in the report cannot be looked at and questioned.
    """
    frame = _frame(arm)
    baseline = _passed(apply_filters(frame, _params(arm))[0])
    filled = _passed(apply_filters(frame.drop(columns=[column]), _params(arm))[0])

    if NEUTRAL_FILLS[column].direction == REMOVES:
        assert filled <= baseline, f"{column} is recorded as removing rows but added some"
    else:
        assert filled >= baseline, f"{column} is recorded as adding rows but removed some"


@pytest.mark.parametrize("arm", ARMS)
def test_a_degraded_fill_shrinks_the_pathogenic_rescue(arm):
    """What "degraded" means, in the numbers: the rescue path itself loses rows.

    The distinguishing property of the escalated set, and the reason it is a different warning
    rather than a louder one. A guideline-only fill can cost the report rows while leaving the
    rescue exactly where it was; these cannot.
    """
    frame = _frame(arm)
    baseline = apply_filters(frame, _params(arm))[1]
    baseline_rescue = baseline.rescue_only + baseline.both

    for column in REQUIRED_INPUTS[arm]:
        diagnostics = apply_filters(frame.drop(columns=[column]), _params(arm))[1]
        rescue = diagnostics.rescue_only + diagnostics.both
        if column in PATHOGENIC_INPUTS[arm]:
            assert rescue < baseline_rescue, column
        else:
            assert rescue == baseline_rescue, column


# ---------------------------------------------------------------------------
# The warnings
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("arm,column", _arm_columns())
def test_every_fill_is_named_and_declared_incomplete(arm, column):
    """A fill is never silent, and never lets the report pass as a complete one.

    The warning no longer says "off parity" — that word left the interface with issue #51 —
    but the statement it carried has to survive the rename, so this pins the replacement
    rather than dropping the assertion.
    """
    diagnostics = apply_filters(_frame(arm).drop(columns=[column]), _params(arm))[1]

    assert diagnostics.notes, "an absent filter-input column produced no warning at all"
    text = " ".join(note_texts(diagnostics))
    assert column in text
    assert "not a complete result" in text


@pytest.mark.parametrize(
    "arm,column",
    [("somatic", "CancerVar"), ("germline", "InterVar")],
)
def test_the_escalated_warning_says_rows_are_missing(arm, column):
    """The escalated block is distinguishable from a dropout line, and says which way it hurts.

    Not asserted word for word — the wording is display prose and pinning it would make every
    rephrasing a test failure. What is asserted is the part a user acts on: that this is
    pathogenic retention, and that the fill can only have taken rows away.

    It used to read ``"missing" in escalated.lower()``, matching the clause *"rows are
    missing, not added"*. Issue #136 replaced that clause because it named no referent and
    the mismatch notice directly above it prices the **other arm** — on the germline
    reference filtered as somatic, 21 kept against 31 with 14 in both, so seven rows are on
    screen that the correct arm would not return. The direction is now stated as a property
    of the fill, which has one referent and cannot be read against those numbers.
    """
    diagnostics = apply_filters(_frame(arm).drop(columns=[column]), _params(arm))[1]

    assert diagnostics.degraded_columns == (column,)
    escalated = next(
        text for text in note_texts(diagnostics) if "PATHOGENIC RETENTION" in text
    )
    assert column in escalated
    assert "can only take rows out of your report" in escalated
    # The direction, not merely a sentence about it: a note claiming the fill could add rows
    # would pass the clause above while saying the opposite of what issue #20 measured.
    assert "never put them in" in escalated


def test_an_ordinary_dropout_is_one_message_however_many_columns_it_covers():
    """Pooled, deliberately: issue #20 measured eight yellow boxes per arm and called it noise.

    The escalations stay per column — there are at most two, and each names a distinct
    consequence — but the plain fills are one line, so a MAF missing five columns produces one
    message a user finishes reading rather than five they scroll past.
    """
    frame = _frame("somatic").drop(columns=["ESCAT", "Variant_Classification", "tumor_f"])
    diagnostics = apply_filters(frame, _params("somatic"))[1]

    assert len(diagnostics.notes) == 1
    for column in ("ESCAT", "Variant_Classification", "tumor_f"):
        assert column in note_texts(diagnostics)[0]


def test_both_directions_are_stated_when_both_happened():
    """A user whose report both gained and lost rows is told both, not the tidier half."""
    frame = _frame("somatic").drop(columns=["ESCAT", "Variant_Classification"])
    diagnostics = apply_filters(frame, _params("somatic"))[1]

    line = note_texts(diagnostics)[0]
    assert "missing" in line
    # The ADDS clause stopped naming the pipeline with issue #51; the direction it states
    # is what this test is about, so it pins the replacement phrasing.
    assert "otherwise have been dropped" in line


def _escalated_warning(arm: str, column: str) -> str:
    """The escalated note for ``column``, from a real filter run rather than built by hand.

    Calling ``_degraded_note`` directly would assert the string against itself. Going through
    ``apply_filters`` means the escalation has to actually have fired for these to say
    anything, which is the same reasoning ``test_only_the_escalated_warnings_are_marked_escalated``
    rests on.
    """
    diagnostics = apply_filters(_frame(arm).drop(columns=[column]), _params(arm))[1]
    return next(t for t in note_texts(diagnostics) if "PATHOGENIC RETENTION" in t)


# ---------------------------------------------------------------------------
# What a fill warning may not claim about the file (issue #136)
# ---------------------------------------------------------------------------
#
# These cover the copy the whole suite was silent about: every assertion below passed
# vacuously before the wording changed, because nothing had ever read these sentences.


@pytest.mark.parametrize("arm,column", _arm_columns())
def test_no_fill_warning_says_the_file_would_have_been_refused(arm, column):
    """MAFigate does not refuse an absent filter input, so it must not say it would.

    The removed clause read *"a file without this column would normally be refused rather
    than filtered"*. Absence fills and warns, here; only an unreadable **value** refuses,
    in ``filters/numeric_columns.py``. So "normally" could only have meant *in the
    pipeline*, which is the comparison issue #51 took out of the interface.

    Its cost was not abstract. Issue #136 measured that this warning fires on the file's own
    arm for 2 of 173 placeable real MAFs; every other time, the arm is wrong — and a germline
    MAF told its data would be refused sends the user to their annotation pipeline instead of
    the Sample Type control one click away.

    Parametrised over every arm and column rather than the escalating pair, because the
    sentence was shared and both notes made it.
    """
    diagnostics = apply_filters(_frame(arm).drop(columns=[column]), _params(arm))[1]

    text = " ".join(note_texts(diagnostics)).lower()
    assert "refused" not in text
    assert "would normally be" not in text


@pytest.mark.parametrize(
    "arm,column",
    [("somatic", "CancerVar"), ("germline", "InterVar")],
)
def test_the_escalated_warning_does_not_claim_the_rescue_collapses(arm, column):
    """``filter_patho`` is an OR, so filling one member drops one term and the rest survive.

    The removed clause claimed the opposite in as many words — *"does not just drop one term
    from an OR — it collapses the safety net"* — against a vendored mask that is
    ``CancerVar | ClinVar_VCF_CLNSIG`` somatic and ``InterVar | ClinVar_VCF_CLNSIG``
    germline. "Unconditional" went with it: the rescue clears the criteria but is still
    subject to the app's own population-frequency filter, exempted only by the narrower
    ClinVar-only ``pathogenic_exemption``.

    The magnitude is why one sentence could not carry the old claim. ``ClinVar_VCF_CLNSIG``
    escalates on **both** arms and drew the same paragraph at a measured 388 → 334 and
    496 → 475 — a −4% loss described as a collapse.
    """
    escalated = _escalated_warning(arm, column)

    assert "collapse" not in escalated.lower()
    assert "unconditional" not in escalated.lower()


@pytest.mark.parametrize(
    "arm,column",
    [("somatic", "CancerVar"), ("germline", "InterVar")],
)
def test_the_escalated_warning_borrows_the_run_report_vocabulary(arm, column):
    """The warning and the run report describe pathogenic retention in one set of words.

    ``_decomposition_summary`` renders *"kept only by pathogenic retention (these did not
    meet the criteria)"* on the same page. Saying it a second way here would be a fourth
    voice for one mechanism — the failure issue #71 spent a ticket undoing. Asserted on the
    clause rather than the sentence, so the surrounding prose stays free to change.
    """
    escalated = _escalated_warning(arm, column)

    assert "did not meet your filter criteria" in escalated


@pytest.mark.parametrize("arm,column", _arm_columns())
def test_every_fill_warning_says_mafigate_filtered_anyway(arm, column):
    """The replacement statement, in one place and read by both notes.

    The claim the old sentence carried has to survive its removal: a report built on a
    stand-in value must not be read as a whole one. What changed is the subject — MAFigate's
    own behaviour rather than a verdict about the user's file.
    """
    diagnostics = apply_filters(_frame(arm).drop(columns=[column]), _params(arm))[1]

    text = " ".join(note_texts(diagnostics))
    assert "MAFigate filtered your file anyway" in text
    assert "not a complete result" in text


def test_the_dropout_clause_verb_agrees_with_one_column():
    """The pooled note is reached with one column as readily as with five.

    ``RENOVO_Class`` alone is exactly what a somatic MAF filtered on the germline arm draws
    beside the escalated ``InterVar`` note — 110 of 110 somatic-detected files in issue
    #136's corpus — so *"``RENOVO_Class`` ... now match nothing"* was the app's commonest
    wrong-arm rendering, not an edge case.
    """
    frame = _frame("germline").drop(columns=["RENOVO_Class"])
    diagnostics = apply_filters(frame, _params("germline"))[1]

    line = next(t for t in note_texts(diagnostics) if "RENOVO_Class" in t)
    assert "now matches nothing" in line
    assert "now match nothing" not in line


def test_the_dropout_clause_verb_agrees_with_several_columns():
    """The other half of the same rule, so a fix in one direction cannot pass alone.

    Without this, making the verb unconditionally singular passes the test above.

    It needs ``skip_pathogenic``, and finding that out is the sharper half of this pair.
    Each arm has exactly **one** guideline source that does not escalate — ``ESCAT`` somatic,
    ``RENOVO_Class`` germline — and every other ``REMOVES`` column is pulled out of the pooled
    line by the escalation. So on an ordinary run the plural is unreachable and the singular
    is the *only* form this clause ever renders, which is why the disagreement survived: the
    app has never once drawn the reading that would have made it look right.
    """
    frame = _frame("germline").drop(
        columns=["InterVar", "RENOVO_Class", "ClinVar_VCF_CLNSIG"]
    )
    diagnostics = apply_filters(frame, _params("germline", skip_pathogenic=True))[1]

    line = next(t for t in note_texts(diagnostics) if "RENOVO_Class" in t)
    assert "now match nothing" in line
    assert "now matches nothing" not in line


@pytest.mark.parametrize("arm,column", _arm_columns())
def test_only_the_escalated_warnings_are_marked_escalated(arm, column):
    """``is_escalated`` picks out exactly the notes that name a degraded column.

    It no longer describes how the UI renders anything, and that clause was deleted rather than
    updated. Before issue #151 this docstring said *"the app shows escalated warnings as errors
    and everything else as a warning box"*, which was true of a two-level slot and is false of a
    three-level one — and the renderer no longer consults :func:`is_escalated` at all, so a
    restatement here would be a second account of a decision that now lives with each message.

    What survives is the property this module can actually speak for: escalation is a fact about
    the *messages*, and it lines up with ``degraded_columns`` one for one. That the escalated
    note and the ``ERROR`` level agree is asserted in ``tests/test_filter_notes.py``, next to the
    renderer it concerns.
    """
    diagnostics = apply_filters(_frame(arm).drop(columns=[column]), _params(arm))[1]
    escalated = [
        text for text in note_texts(diagnostics) if is_escalated(text)
    ]

    assert len(escalated) == len(diagnostics.degraded_columns)
    for warning in escalated:
        assert any(name in warning for name in diagnostics.degraded_columns)


def test_nothing_tells_the_user_a_filter_was_skipped_for_a_missing_column():
    """The seven "X filter not applied" warnings are gone and must not come back.

    They described a state the app cannot reach. Each guarded a filter clause with
    ``if column in maf.columns`` and skipped it when the column was absent — but the vendored
    filters index every input unguarded, so the app never got as far as the warning; it raised.
    Issue #20 measured all eleven filter inputs and found **11 of 11 raise**, with
    ``CIViC_Evidence_Level`` the sole exception the *pipeline* guards.

    They are worth a standing assertion rather than a one-off deletion, because the phrasing is
    the obvious thing to reach for when adding a column check — and it is now not merely
    unreachable but false: no filter is skipped for a missing column, the column is filled and
    the filter runs.
    """
    forbidden = ("filter not applied", "filters will be disabled", "filter will be skipped")
    searched = [
        STREAMLIT_APP / "filters",
        STREAMLIT_APP / "page_modules",
        STREAMLIT_APP / "components",
        STREAMLIT_APP / "config",
    ]

    offenders = []
    for directory in searched:
        for path in directory.rglob("*.py"):
            text = path.read_text(encoding="utf-8").lower()
            offenders += [
                f"{path.relative_to(STREAMLIT_APP)}: {phrase}"
                for phrase in forbidden
                if phrase in text
            ]

    assert not offenders, (
        "a 'filter skipped for a missing column' claim is back in the app: "
        f"{offenders}. Missing filter inputs are filled and reported, not skipped."
    )


# ---------------------------------------------------------------------------
# The happy path stays cheap
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("arm", ARMS)
def test_a_complete_maf_is_not_copied_and_not_warned_about(arm):
    """The identity that keeps the common case free.

    ``frame_for_masks`` returning the caller's own object is the assertion, not a proxy for it:
    every other guarantee in this module is about the rare frame, and this is the one about the
    frame every real user loads.
    """
    frame = _frame(arm)
    plan = plan_fills(frame, arm)

    assert plan.filled == ()
    assert plan.warnings() == ()
    assert plan.frame_for_masks(frame) is frame


def test_the_fill_step_allocates_nothing_on_a_large_maf():
    """Asserted against a frame big enough that a copy could not hide in the noise.

    400,000 rows of float64 — more than twice the whole 181,566-row reference — so a full-frame
    copy is tens of megabytes and any threshold below that catches it. The measurement is of
    the fill step alone: ``apply_filters`` copies once in ``_label``, which is documented,
    deliberate, and not what this ticket added.
    """
    columns = REQUIRED_INPUTS["somatic"]
    frame = pd.DataFrame({column: [1.0] * 400_000 for column in columns})
    footprint = int(frame.memory_usage(deep=True).sum())

    tracemalloc.start()
    try:
        plan = plan_fills(frame, "somatic")
        masked = plan.frame_for_masks(frame)
        _, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()

    assert masked is frame
    assert peak < footprint / 10, (
        f"the happy path allocated {peak:,} bytes against a {footprint:,}-byte frame, which "
        "is the shape of a copy"
    )


def test_the_unhappy_path_does_copy_and_leaves_the_caller_frame_alone():
    """The other side of the same guarantee: the fill lands on a throwaway, never in place.

    Streamlit hands this seam its session state. Filling in place would put the fill in front
    of the user by the shortest possible route — the very next rerender.
    """
    frame = _frame("somatic").drop(columns=["tumor_f"])
    plan = plan_fills(frame, "somatic")
    masked = plan.frame_for_masks(frame)

    assert masked is not frame
    assert "tumor_f" in masked.columns
    assert "tumor_f" not in frame.columns


def test_an_unknown_arm_is_not_quietly_given_an_empty_contract():
    """A plan for an arm with no derived contract would fill nothing and filter anyway."""
    with pytest.raises(KeyError):
        plan_fills(pd.DataFrame({"a": [1]}), "tumour")


# ---------------------------------------------------------------------------
# Hugo_Symbol — the input the app guards rather than fills
# ---------------------------------------------------------------------------


def test_a_gene_restriction_on_a_maf_with_no_hugo_symbol_widens_rather_than_empties():
    """The one filter input whose neutral answer is *no file*, not a value.

    Any value stood in for a gene symbol matches none of the requested genes, so filling would
    empty the report — the invisible direction. Dropping the restriction instead leaves the
    user with more rows than they asked for, and a warning telling them so.
    """
    frame = _frame("somatic").drop(columns=["Hugo_Symbol"])
    params = _params("somatic", filter_genes=["TP53", "BRCA1"])

    labelled, diagnostics = apply_filters(frame, params)

    assert _passed(labelled) > 0
    assert diagnostics.filled_columns == ()
    assert any("Hugo_Symbol" in text for text in note_texts(diagnostics))


def _gene_warning() -> str:
    """The note a somatic MAF with no ``Hugo_Symbol`` draws against a gene restriction.

    Raises rather than returns ``None`` if it stops being drawn, so the assertions below cannot
    pass by the warning disappearing — the failure mode #136's own copy guards were written
    after, having all passed vacuously for want of anything reading the sentences.
    """
    frame = _frame("somatic").drop(columns=["Hugo_Symbol"])
    params = _params("somatic", filter_genes=["TP53", "BRCA1"])
    diagnostics = apply_filters(frame, params)[1]
    return next(t for t in note_texts(diagnostics) if "Hugo_Symbol" in t)


@pytest.mark.parametrize("forbidden", ["would normally", "refused"])
def test_the_gene_warning_does_not_say_the_file_would_have_been_refused(forbidden):
    """The last copy of the sentence issue #136 deleted, removed by issue #150.

    It closed *"A file without that column would normally be refused rather than filtered"* —
    the same clause, in another module, still shipping because #136 was reading
    ``absent_columns.py`` and this copy was out of its view. MAFigate refuses no file for an
    absent column: here the restriction is dropped, and next door an unreadable *value* is what
    refuses.

    The gene case is the one place the sentence was arguably true, and that is the trap in it —
    it was true of the **pipeline**, whose gene clause reads ``Hugo_Symbol``, and never of
    MAFigate, which writes the gene file and simply does not write one. A warning states what
    MAFigate did.
    """
    assert forbidden not in _gene_warning().lower()


def test_the_gene_warning_still_says_which_way_the_report_moved():
    """Deleting the false clause must not delete the direction, which is the point of the note.

    This is the one member of the missing-column family whose direction runs the *other* way:
    #136 settled that a fill "can only take rows out of your report", and nothing is filled
    here — the restriction is dropped and the report gets **wider**. Issue #28's *extra rows are
    visible, missing rows are not* is why that is said rather than escalated, and why it has to
    be said at all.

    Held against *the report the user asked for*, which is a comparison the user can make, and
    not against what another tool would have done with their file.
    """
    warning = _gene_warning()

    assert "was not applied" in warning
    assert "wider than the one you asked for" in warning


def test_the_gene_warning_warns_without_escalating():
    """``⚠️`` since issue #151, and still not escalated. Both halves matter.

    It was ``ℹ️`` between #150 and #151; :mod:`filters.notes` is where that is explained.

    Not escalated is the unchanged half, and the half worth asserting beside the wording: issue
    #28's *extra rows are visible, missing rows are not* still holds, and the rows a dropped gene
    restriction lets through are on screen to be judged. Escalation is for rows that are **gone**.

    The level itself, and this note's agreement with the frequency-skip note, are asserted in
    ``tests/test_filter_notes.py`` — this module owns the sentence.
    """
    warning = _gene_warning()

    assert warning.startswith("⚠️")
    assert not is_escalated(warning)


def test_hugo_symbol_is_not_filled_because_the_app_is_what_guards_it():
    """Recorded as a property of the derivation, not just of the wiring above.

    ``Hugo_Symbol`` is read inside the vendored gene-file branch, and the app writes that file
    — so unlike ``CIViC_Evidence_Level``, which the *pipeline* guards, this one is guarded by a
    decision the app makes on every call.
    """
    for arm in ARMS:
        assert "Hugo_Symbol" not in REQUIRED_INPUTS[arm]
    assert "Hugo_Symbol" not in NEUTRAL_FILLS


def test_a_maf_with_no_hugo_symbol_and_no_gene_list_is_not_warned_about():
    """No restriction was asked for, so nothing was lost and there is nothing to say."""
    frame = _frame("somatic").drop(columns=["Hugo_Symbol"])
    diagnostics = apply_filters(frame, _params("somatic"))[1]

    assert diagnostics.notes == ()


# ---------------------------------------------------------------------------
# The plan object itself
# ---------------------------------------------------------------------------


def test_the_plan_reports_the_same_columns_it_warns_about():
    """The structured fields and the prose come from one object, and this pins that.

    A caller rendering ``filled_columns`` and a caller rendering ``warnings`` must not be able
    to disagree about which columns were filled — the UI does both, side by side.
    """
    frame = _frame("germline").drop(columns=["InterVar", "ESCAT_NOT_A_COLUMN"], errors="ignore")
    labelled, diagnostics = apply_filters(frame, _params("germline"))

    plan = FillPlan(arm="germline", filled=("InterVar",), degraded=("InterVar",))
    assert diagnostics.filled_columns == plan.filled
    assert diagnostics.degraded_columns == plan.degraded
    assert set(diagnostics.notes) == set(plan.warnings())
    assert labelled is not frame
