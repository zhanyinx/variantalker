"""What the app says when the loaded MAF and the selected arm disagree (issue #135).

The app has always been arm-*correct* and arm-*blind*: every per-arm decision reads the
selected arm, and nothing read the file back to ask whether that setting suits it. The
reported session is a germline MAF filtered on the somatic arm, which draws
``❌ PATHOGENIC RETENTION DEGRADED — CancerVar column not found`` — true of the somatic arm,
and useless, because the file is fine and the setting is wrong.

Detect and warn; never auto-switch
----------------------------------
Settled by the dev before this surface was designed. The file does not get to override a
deliberate choice, and the arm drives presets, thresholds, the gene panel and the
Pathogenicity Overview's sources, so a silent switch would move all of them.

Why the consequence sentence is measured and not written down
--------------------------------------------------------------
A wrong-arm report is not simply a smaller one, and *how* it is wrong depends on the
direction. Over real files the asymmetry is sharp — a somatic MAF on the germline arm gained
nothing in 8 of 8, while a germline MAF on the somatic arm gained variants in 8 of 8, one of
them with **zero** overlap between the two reports (17 kept on germline, 13 on somatic, none
in common) and another where the *wrong* arm's report was the larger of the two, so size is
not even a cue. No single sentence is true both ways, so the sentence is chosen by
:class:`~filters.arm_detection.ArmComparison`'s own numbers for the open file rather than by a
rule about arms that a change to either contract could quietly falsify.

The same measurement supplies the counts the notice quotes, which are the counts the button
delivers: it re-seeds from the other arm's contract, and that contract is what was priced.

What the escalated warning beneath this says, and why nothing here decides it
-----------------------------------------------------------------------------
Issue #136 settled it, and settled it **without** either surface learning about the other.
The disagreement recorded here was real: that warning ended *"Rows are missing, not added"*,
while this module reports 21 kept against germline's 31 with 14 in both on the reference — so
seven rows are on screen that the correct arm would not return. Both sentences were true of
their own counterfactual, the fill's being *the same arm with a real column*, and nothing on
screen said which. It now states the direction as a property of the fill — *a stand-in value
can only take rows out of your report* — which has one referent and cannot be read against
these numbers.

The rest of #136's answer is the same shape. The clause that sent a user to their annotation
pipeline — *"a file without this column would normally be refused rather than filtered"* — was
false of MAFigate on **either** arm, absence being filled rather than refused, so it was
replaced for every run rather than branched on a detected mismatch. ``filters/absent_columns``
therefore still knows nothing about arms, and this notice, drawn above it, remains the whole
of the app's account of *why* the column is missing.
"""

import streamlit as st

from filters.arm_detection import ANNOTATOR_MARKERS, ArmComparison, ArmEvidence


def _columns(names) -> str:
    """Name column names the way the user's own header row spells them.

    Raw MAF column spellings, which the Style section permits since issue #79 in exactly this
    setting: these are the user's file's names, so a description of them would be harder to
    match against the header than the name itself.
    """
    quoted = [f"`{name}`" for name in names]
    if len(quoted) <= 1:
        return "".join(quoted)
    return f"{', '.join(quoted[:-1])} and {quoted[-1]}"


def _absent(names) -> str:
    """The same names under a negation, which needs a different conjunction.

    ``no `InterVar` and `RENOVO_Class``` reads as denying the pair while allowing either
    one, which is the opposite of what the rule found: this file carries **neither**. The
    positive joiner cannot be reused here.
    """
    quoted = [f"`{name}`" for name in names]
    if len(quoted) == 1:
        return f"no {quoted[0]}"
    if len(quoted) == 2:
        return f"neither {quoted[0]} nor {quoted[1]}"
    return f"none of {', '.join(quoted[:-1])} or {quoted[-1]}"


def _variants(count: int) -> str:
    """``1 variant``, ``3 variants`` — the report can hold exactly one."""
    return f"{count} variant" if count == 1 else f"{count} variants"


def _evidence_sentence(evidence: ArmEvidence, detected: str, other: str) -> str:
    """What in this file says which arm annotated it — both halves of the rule.

    The rule is an exclusive-or, so the absence is evidence as much as the presence is, and
    the absent column is the one the escalated warning beside this is complaining about. A
    notice that named only what the file carries would leave the two messages looking like
    they were about different files.

    The absent half is read off :data:`ANNOTATOR_MARKERS` rather than off the file: a
    detected mismatch means the file carries **none** of ``other``'s markers, so asking the
    evidence what it carries for that arm can only ever return nothing to name.
    """
    carried = _columns(evidence.evidence_for(detected))
    return f"It carries {carried}, and {_absent(ANNOTATOR_MARKERS[other])}."


def _consequence(comparison: ArmComparison | None, current: str, detected: str) -> str:
    """What filtering on the wrong arm has done to *this* file.

    Four sentences for four measured states, and the fourth is the one a rule would have got
    wrong: the two arms can keep exactly the same variants, which 7 of 8 sampled small
    somatic files did. Claiming a loss there would name one that did not happen — the
    argument ``plan_fills`` already makes for not escalating under ``skip_pathogenic``.
    """
    if comparison is None:
        # The counterfactual could not be computed — either the run was refused, so there is
        # no report to compare, or the other arm could not filter this file. The mismatch is
        # still real, being read from the header, so the notice keeps the claim that holds
        # in every direction and drops only the numbers. It must not say "the report below":
        # a refused run has no report at all.
        return (
            f"{current.title()} and {detected} settings read different evidence columns, so "
            f"filtering this file as {current} does not give the report it would get on "
            f"{detected}."
        )

    if comparison.identical:
        return (
            f"On this file the two happen to keep the same {_variants(comparison.kept)}, so "
            "your report is unaffected. Anything you change from here will differ, though — "
            "the arm decides which evidence the filter reads."
        )

    if comparison.gained == 0:
        return (
            f"{current.title()} settings only drop variants from this file — they add "
            f"nothing. Your current settings keep {_variants(comparison.kept)}; "
            f"{detected}'s opening settings would keep {comparison.would_keep}."
        )

    overlap = (
        "and **no** variant is in both reports"
        if comparison.in_both == 0
        else f"and only {comparison.in_both} of those variants are in both reports"
    )
    return (
        f"{current.title()} settings do not simply cut this file harder — they keep a "
        f"*different* set of variants. Your current settings keep "
        f"{_variants(comparison.kept)}; {detected}'s opening settings would keep "
        f"{comparison.would_keep}, {overlap}."
    )


def render_mismatch_notice(
    *,
    evidence: ArmEvidence,
    current_arm: str,
    comparison: ArmComparison | None,
    customised: bool,
    on_switch,
) -> None:
    """Say that this file was annotated for the other arm, and offer the way across.

    Args:
        evidence: what the header says, from :func:`~filters.arm_detection.read_arm_evidence`.
        current_arm: the arm MAFigate is set to.
        comparison: what the detected arm would keep, or ``None`` if it could not be priced.
        customised: whether the session's parameters differ from ``current_arm``'s contract.
            The switch re-seeds, so it discards edits — but only where there are any, and on
            the reported session (a cached arm, untouched settings) there are none. Saying
            "your settings will be lost" to a user with nothing to lose is how the next real
            warning gets ignored.
        on_switch: called when the user presses the button. Takes the detected arm.
    """
    detected = evidence.detected
    if detected is None:
        return

    st.warning(
        f"⚠️ **This file was annotated for {detected} analysis, and MAFigate is set to "
        f"{current_arm}.**\n\n"
        f"{_evidence_sentence(evidence, detected, current_arm)}\n\n"
        f"{_consequence(comparison, current_arm, detected)}"
    )

    if customised:
        # Drawn above the button rather than inside its label: the label says what the
        # button does, and this says what it costs, which is a different sentence and only
        # sometimes true.
        st.caption(
            f"Switching replaces your current settings with the ones MAFigate opens with "
            f"for {detected}. Changes you have made on the {current_arm} settings will be "
            "lost."
        )

    if st.button(
        f"🔄 Switch to {detected} and re-filter",
        key="switch_arm_to_detected",
        type="primary",
    ):
        on_switch(detected)


def render_ambiguous_notice(evidence: ArmEvidence, current_arm: str) -> None:
    """Say that the file carries both arms' annotations, so the app cannot place it.

    The *only* unplaceable case this surface speaks about. A file carrying **neither** arm's
    markers is left alone deliberately: the escalated fill warning is already on that screen
    saying the report is not complete, and a second box beside it is the wall issue #93 spent
    a ticket avoiding. This case has no such warning — neither arm fills anything, because
    every column both arms read is present — so without this the app is silent on both arms
    about a file that is genuinely ambiguous.

    No switch button: there is no detected arm to switch to. Pointing at where the arm is
    set is all this can honestly offer.
    """
    if not evidence.ambiguous:
        return

    somatic = _columns(evidence.evidence_for("somatic"))
    germline = _columns(evidence.evidence_for("germline"))
    # "which comes from" cannot be shared by the two halves: one arm's markers are a single
    # column and the other's are a pair, so one verb is wrong whichever is chosen. Naming
    # the source without a verb is true of both.
    st.warning(
        "⚠️ **MAFigate cannot tell which analysis this file was annotated for.**\n\n"
        f"It carries {somatic} from somatic annotation, and {germline} from germline "
        "annotation. A MAF normally carries one set or the other.\n\n"
        f"It is being filtered as {current_arm}. If that is not right, change it under "
        "Configure Parameters."
    )
