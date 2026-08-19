"""Which arm a MAF was *annotated* for, and what filtering it on the other one costs.

The third member of the ``filters/`` family that reads a MAF's header rather than its values,
beside ``absent_columns`` and ``numeric_columns`` — and the one that asks the converse
question. Those two ask *does this file carry what this arm's filter reads?*; this asks *which
arm's annotator wrote this file?* The app has always been arm-correct and arm-blind: every
per-arm decision reads ``filter_params["sample_type"]`` and nothing anywhere read the file back
to check that setting suits it (issue #135).

Which columns the annotator wrote is not which columns the filter reads
----------------------------------------------------------------------
This is the distinction the whole module turns on, and issue #134 measured it rather than
argued it. ``REQUIRED_INPUTS`` differenced gives ``{CancerVar, ESCAT}`` against
``{InterVar, RENOVO_Class}`` — the right shape with one wrong member. **ESCAT is not
somatic-only**: it is in 55 of 57 real germline files, every one that was annotated at all.
The germline *filter* does not read it; the annotation step *writes* it on both arms. A
detector derived straight from the filter contract therefore detects **0 of 57** germline
files, because ESCAT drags every annotated one into "both".

So :data:`ANNOTATOR_MARKERS` is written down rather than derived, and guarded *against* the
derivation instead — see ``tests/test_arm_detection.py``. Each arm's set must stay a non-empty
subset of ``REQUIRED_INPUTS[arm] - REQUIRED_INPUTS[other]``, which is one narrowing, ESCAT
removed. A vendored change that stops the somatic filter requiring ``CancerVar`` then fails a
test rather than quietly changing what arm the app infers; the part the source can justify
stays derived, and the part it cannot — that ESCAT is written on both arms — is written once
with that measurement as its reason.

What the rule scores, on 142 real MAFs with provenance
------------------------------------------------------
===========  ====  =========  ======================  =======
truth          n    correct    called the other arm    unknown
===========  ====  =========  ======================  =======
somatic        85         75                       0       10
germline       57         55                       0        2
===========  ====  =========  ======================  =======

130 of 142 classified, **0 wrong in either direction**, 12 unknown. Counts rather than a
percentage, because the denominator is small enough that a percentage would flatter it. On the
2026 GERSOM cohort — the pipeline output a user is most likely to arrive with — it is 104 for
104 with no ambiguity at all; the unknowns cluster in old material.

**Presence settles it; no value pass is needed.** None of these three columns is ever a header
over an empty column across the corpus (0 of 98, 0 of 59, 0 of 61). ``ESCAT`` is, in 60 of 164
files — the one column that is routinely present-and-blank is the one column that cannot be
trusted anyway, so ``absent_columns``'s presence-only question is the right one here too.

Why an exclusive-or, and not "look for CancerVar first"
-------------------------------------------------------
No real file carries ``CancerVar`` and ``InterVar`` together — 0 of 179. But a hand-merged MAF
could, and this repo's own test fixture once did, which is why :func:`read_arm_evidence`
reports that case as **ambiguous** rather than resolving it by precedence. Guessing there
would be the one way this rule could be confidently wrong.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Mapping, Sequence

import pandas as pd

from filters.variant_filters import MAFIGATE_FILTER, PASS, apply_filters

#: The columns each arm's **annotation step** writes and the other's never does. Not derived
#: from the filter contract — see the module docstring for the measurement that says why, and
#: ``tests/test_arm_detection.py`` for the guard that keeps this a narrowing of it rather than
#: an unrelated list.
ANNOTATOR_MARKERS: Mapping[str, tuple[str, ...]] = {
    "somatic": ("CancerVar",),
    "germline": ("InterVar", "RENOVO_Class"),
}

#: The other arm, for a caller that has one. Written here rather than spelled as a conditional
#: at each call site, so "the two arms" is one fact with one home.
OTHER_ARM: Mapping[str, str] = {"somatic": "germline", "germline": "somatic"}


@dataclass(frozen=True)
class ArmEvidence:
    """What a MAF's header says about the arm it was annotated for.

    Attributes:
        detected: ``"somatic"``, ``"germline"``, or ``None`` where the header cannot place
            the file. ``None`` is a legitimate answer and callers must be able to hold it —
            12 of 142 real files land here.
        ambiguous: whether the file carries **both** arms' markers. A different kind of
            ``None`` from carrying neither: no real file does this, but a hand-merged one
            can, and it is the only unplaceable case the rest of the app is silent about.
            A file carrying neither already draws the escalated fill warning on either arm.
        present: the marker columns this file actually carries, per arm, in the order
            :data:`ANNOTATOR_MARKERS` lists them — the evidence a notice quotes back.
    """

    detected: str | None
    ambiguous: bool
    present: Mapping[str, tuple[str, ...]]

    def disagrees_with(self, arm: str) -> bool:
        """Whether this file was annotated for an arm other than ``arm``.

        False when the header cannot place the file: *cannot tell* is not *disagrees*, and a
        detector that reported one as the other would spend its credibility on the 12 files
        it knows least about.
        """
        return self.detected is not None and self.detected != arm

    def evidence_for(self, arm: str) -> tuple[str, ...]:
        """The marker columns this file carries for ``arm``."""
        return self.present.get(arm, ())


def read_arm_evidence(columns: Iterable[str]) -> ArmEvidence:
    """Read the arm off a MAF's column names.

    Args:
        columns: the MAF's header. Only names are read — never values, for the reason in the
            module docstring: these three columns are never present-and-blank, and the one
            that is cannot be trusted anyway.
    """
    header = set(columns)
    present = {
        arm: tuple(column for column in markers if column in header)
        for arm, markers in ANNOTATOR_MARKERS.items()
    }

    somatic = bool(present["somatic"])
    germline = bool(present["germline"])
    if somatic and not germline:
        detected: str | None = "somatic"
    elif germline and not somatic:
        detected = "germline"
    else:
        detected = None

    return ArmEvidence(detected=detected, ambiguous=somatic and germline, present=present)


@dataclass(frozen=True)
class ArmComparison:
    """What this file's report would be on the other arm, measured rather than predicted.

    Attributes:
        kept: how many variants the report on screen holds. Taken from that report rather
            than recomputed, so the notice cannot quote a number the user cannot see.
        would_keep: how many the other arm's opening settings keep.
        in_both: how many variants are in both reports.
    """

    kept: int
    would_keep: int
    in_both: int

    @property
    def gained(self) -> int:
        """Variants the current settings keep that the other arm's would not."""
        return self.kept - self.in_both

    @property
    def lost(self) -> int:
        """Variants the other arm's settings would keep that the current ones do not."""
        return self.would_keep - self.in_both

    @property
    def identical(self) -> bool:
        """Whether the two arms keep exactly the same variants of this file.

        Not a curiosity: 7 of 8 sampled somatic files with a provenance label kept one
        variant on either arm, so the same-report case is common on small MAFs. A notice
        that named a loss there would name one that did not happen — ``plan_fills``'s own
        argument for not escalating under ``skip_pathogenic``, in a different place.
        """
        return self.gained == 0 and self.lost == 0


def price_other_arm(
    maf: pd.DataFrame, other_params: Mapping, kept_index: Sequence
) -> ArmComparison | None:
    """Filter ``maf`` on the other arm's settings, to count what switching would get.

    The direction of a wrong-arm error is **not** a property of the arms that can be written
    down once. Measured over real files it is sharply asymmetric — a somatic file on the
    germline arm gained nothing in 8 of 8, while a germline file on the somatic arm gained
    variants in 8 of 8, often with almost no overlap (one file: 17 kept on germline, 13 on
    somatic, **zero in common**). So the notice reads the direction off *this* file's numbers
    rather than off a rule that could go stale, which is what this function exists to supply.

    Args:
        maf: the open file, unfiltered.
        other_params: the other arm's opening settings — the contract, because that is
            exactly what the switch button re-seeds, so the count quoted is the count
            delivered.
        kept_index: the index of the report currently on screen.

    Returns:
        The comparison, or ``None`` if the other arm cannot filter this file. A caller that
        gets ``None`` still has a mismatch to report; it just cannot price it. Filtering is
        ~11ms on a real 1,846-row, 388-column MAF, which is why this runs per render rather
        than being cached against a file that may have been replaced.
    """
    try:
        labelled, _ = apply_filters(maf, other_params)
    except Exception:
        # Deliberately broad. This is a *counterfactual* — nothing the user asked for — so a
        # failure to compute it must cost the notice its numbers and nothing else. The
        # mismatch itself is read from the header and cannot fail.
        return None

    would_keep = set(labelled.index[labelled[MAFIGATE_FILTER] == PASS])
    kept = set(kept_index)
    return ArmComparison(
        kept=len(kept),
        would_keep=len(would_keep),
        in_both=len(kept & would_keep),
    )
