"""What the app says when the loaded MAF and the selected arm disagree (issue #135).

The app was arm-*correct* and arm-*blind*: every per-arm decision read the selected arm and
nothing read the file back. A germline MAF filtered on the somatic arm drew
``❌ PATHOGENIC RETENTION DEGRADED — CancerVar column not found`` — true of the somatic arm,
and useless, because the file is fine and the setting is wrong.

What is asserted here
---------------------
* **The detector set stays a narrowing of the filter contract.** Written down rather than
  derived, because ``ESCAT`` is written by the annotator on *both* arms while being a somatic
  filter input — so the derivation detects 0 of 57 real germline files. The guard is issue
  #134's: a non-empty subset of ``REQUIRED_INPUTS[arm] - REQUIRED_INPUTS[other]``, so a
  vendored change that stops the somatic filter requiring ``CancerVar`` fails here rather
  than quietly changing what arm the app infers.
* **Cannot tell is not disagrees.** A file carrying neither arm's markers draws no mismatch
  notice at all — the escalated fill warning is already on that screen — and a file carrying
  *both* draws a different notice with no switch offered, because there is no arm to switch
  to.
* **The consequence sentence is chosen by this file's numbers**, not by a rule about arms.
  Measured over real files the direction is asymmetric (a somatic MAF on the germline arm
  gained nothing in 8 of 8; a germline MAF on the somatic arm gained variants in 8 of 8, once
  with zero overlap), so no single sentence is true both ways.
* **The jump to Results is withheld on a mismatch**, and the upload token is stamped anyway —
  it guards the re-read, not the navigation.
* **The switch re-seeds, re-filters, grants the withheld jump, and says so in issue #133's
  vocabulary** rather than a fifth one.

Every assertion here was mutation-checked; the mutations are named on the test that catches
them, because a green test over a wired-up page proves nothing on its own (issues #67, #83,
#90).
"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path
from unittest.mock import Mock, patch

import pytest

STREAMLIT_APP = Path(__file__).resolve().parent.parent
if str(STREAMLIT_APP) not in sys.path:
    sys.path.insert(0, str(STREAMLIT_APP))

from tests.fakes import FakeSessionState  # noqa: E402

FIXTURES = STREAMLIT_APP / "tests" / "fixtures" / "parity"
GERMLINE_MAF = FIXTURES / "germline_reference.maf"
SOMATIC_MAF = FIXTURES / "somatic_reference.maf"

#: The clause that identifies the mismatch notice wherever it is drawn. Deliberately not the
#: whole literal — the arms and the counts vary — so a reworded sentence does not fail a test
#: about *whether* the app speaks.
#:
#: Not ``"was annotated for"``, which was the first choice and is a **substring of the
#: ambiguous notice's own headline** ("cannot tell which analysis this file was annotated
#: for"). The two clauses have to be unsatisfiable by one another or a build drawing the
#: wrong notice passes both tests — which is exactly what happened, caught here rather than
#: reasoned about.
MISMATCH = "and MAFigate is set to"

#: The ambiguous notice's own clause. See above: disjoint from ``MISMATCH`` by construction,
#: since that notice never claims the app *is set to* a detected arm.
AMBIGUOUS = "cannot tell which analysis"


def _maf_with_columns_removed(tmp_path, columns, source=SOMATIC_MAF, name="short.maf"):
    """``source`` with ``columns`` dropped, written out as a MAF on disk."""
    from utils import read_maf

    maf = read_maf(str(source))
    present = [column for column in columns if column in maf.columns]
    assert present == list(columns), (
        f"{source.name} does not carry {sorted(set(columns) - set(present))}, so dropping "
        "them proves nothing"
    )
    path = tmp_path / name
    maf.drop(columns=present).to_csv(path, sep="\t", index=False)
    return path


def _maf_carrying_both_arms(tmp_path):
    """A hand-merged MAF: somatic markers *and* germline ones.

    No real file does this — 0 of 179 measured — but this repo's own fixture once did, which
    is why the rule is an exclusive-or rather than "look for ``CancerVar`` first".
    """
    from utils import read_maf

    maf = read_maf(str(SOMATIC_MAF))
    path = tmp_path / "both_arms.maf"
    maf.assign(InterVar="Pathogenic", RENOVO_Class="HP").to_csv(
        path, sep="\t", index=False
    )
    return path


class TestTheDetectorSetIsANarrowingOfTheContract(unittest.TestCase):
    """Issue #134's guard: the part the vendored source can justify stays derived.

    Mutation-checked by adding ``ESCAT`` back to the somatic markers — the very member the
    measurement removed — which fails the subset assertion below.
    """

    def test_each_arms_markers_are_that_arms_filter_inputs_alone(self):
        from filters.absent_columns import REQUIRED_INPUTS
        from filters.arm_detection import ANNOTATOR_MARKERS, OTHER_ARM

        for arm, markers in ANNOTATOR_MARKERS.items():
            other = OTHER_ARM[arm]
            only_this_arm = set(REQUIRED_INPUTS[arm]) - set(REQUIRED_INPUTS[other])
            self.assertTrue(
                markers, f"{arm} has no markers, so it can never be detected"
            )
            self.assertTrue(
                set(markers) <= only_this_arm,
                f"{arm} markers {sorted(set(markers) - only_this_arm)} are not columns "
                f"only the {arm} filter reads — see issue #134 before widening this",
            )

    def test_escat_is_deliberately_not_a_marker(self):
        """The one narrowing, pinned with its reason.

        ``ESCAT`` *is* a somatic-only filter input, so the subset guard above admits it. What
        rules it out is a fact no derivation can reach: the annotation step writes it on both
        arms, in 55 of 57 real germline files. Without this the guard would pass a set that
        detects 0 of 57 germline files.
        """
        from filters.absent_columns import REQUIRED_INPUTS
        from filters.arm_detection import ANNOTATOR_MARKERS

        self.assertIn("ESCAT", REQUIRED_INPUTS["somatic"])
        self.assertNotIn("ESCAT", ANNOTATOR_MARKERS["somatic"])

    def test_the_two_arms_markers_do_not_overlap(self):
        from filters.arm_detection import ANNOTATOR_MARKERS

        self.assertEqual(
            set(ANNOTATOR_MARKERS["somatic"]) & set(ANNOTATOR_MARKERS["germline"]),
            set(),
            "a shared marker would make every annotated file ambiguous",
        )


class TestReadingTheArmOffTheHeader(unittest.TestCase):
    """Presence settles it — no value pass, for the reason the module docstring measures."""

    def _evidence(self, path):
        from filters.arm_detection import read_arm_evidence
        from utils import read_maf

        return read_arm_evidence(read_maf(str(path)).columns)

    def test_the_reference_mafs_are_placed_correctly(self):
        self.assertEqual(self._evidence(GERMLINE_MAF).detected, "germline")
        self.assertEqual(self._evidence(SOMATIC_MAF).detected, "somatic")

    def test_a_placed_file_disagrees_only_with_the_other_arm(self):
        germline = self._evidence(GERMLINE_MAF)
        self.assertTrue(germline.disagrees_with("somatic"))
        self.assertFalse(germline.disagrees_with("germline"))

    def test_a_file_carrying_neither_marker_is_unplaced_and_disagrees_with_nothing(self):
        """*Cannot tell* is not *disagrees* — the claim that keeps the 12 quiet.

        Mutation-checked by having ``disagrees_with`` return ``detected != arm``, which
        reports every unplaceable file as a mismatch on both arms.
        """
        from filters.arm_detection import read_arm_evidence
        from utils import read_maf

        maf = read_maf(str(SOMATIC_MAF)).drop(columns=["CancerVar"])
        evidence = read_arm_evidence(maf.columns)

        self.assertIsNone(evidence.detected)
        self.assertFalse(evidence.ambiguous)
        self.assertFalse(evidence.disagrees_with("somatic"))
        self.assertFalse(evidence.disagrees_with("germline"))

    def test_a_file_carrying_both_markers_is_ambiguous_rather_than_resolved(self):
        """Not "look for CancerVar first" — guessing is the one way this could be wrong."""
        from filters.arm_detection import read_arm_evidence
        from utils import read_maf

        maf = read_maf(str(SOMATIC_MAF)).assign(InterVar="Pathogenic")
        evidence = read_arm_evidence(maf.columns)

        self.assertIsNone(evidence.detected)
        self.assertTrue(evidence.ambiguous)

    def test_one_germline_marker_is_enough(self):
        """``RENOVO_Class`` alone places a file; the rule is an *or* over that arm's set."""
        from filters.arm_detection import read_arm_evidence
        from utils import read_maf

        maf = read_maf(str(GERMLINE_MAF)).drop(columns=["InterVar"])
        self.assertEqual(read_arm_evidence(maf.columns).detected, "germline")


class TestPricingTheOtherArm(unittest.TestCase):
    """The counterfactual the notice quotes, and the button then delivers."""

    def _priced(self, path, current):
        from config.pipeline_params import pipeline_params
        from filters.arm_detection import OTHER_ARM, price_other_arm
        from filters.variant_filters import MAFIGATE_FILTER, PASS, apply_filters
        from utils import read_maf

        maf = read_maf(str(path))
        labelled, _ = apply_filters(maf, pipeline_params(current))
        kept = labelled.index[labelled[MAFIGATE_FILTER] == PASS]
        return price_other_arm(maf, pipeline_params(OTHER_ARM[current]), kept)

    def test_a_germline_maf_on_the_somatic_arm_keeps_variants_germline_would_drop(self):
        """The finding that decides the copy: the wrong arm is not a subset.

        This is why no single consequence sentence works, and it is also what issue #136 had
        to answer for the escalated fill warning beside this notice. That warning used to end
        *"Rows are missing, not added"* — with no referent, beside these numbers, on one
        screen. #136 gave the direction a subject rather than a comparison: the **fill** can
        only take rows out. That is true whatever this counterfactual reports, so the two
        surfaces no longer have to agree about a number to stop contradicting each other.
        """
        comparison = self._priced(GERMLINE_MAF, "somatic")

        self.assertGreater(comparison.gained, 0)
        self.assertGreater(comparison.lost, 0)
        self.assertFalse(comparison.identical)

    def test_a_somatic_maf_on_the_germline_arm_only_loses(self):
        comparison = self._priced(SOMATIC_MAF, "germline")

        self.assertEqual(comparison.gained, 0)
        self.assertGreater(comparison.lost, 0)

    def test_the_counts_are_the_reports_own(self):
        from config.pipeline_params import pipeline_params
        from filters.arm_detection import price_other_arm
        from filters.variant_filters import MAFIGATE_FILTER, PASS, apply_filters
        from utils import read_maf

        maf = read_maf(str(GERMLINE_MAF))
        somatic, _ = apply_filters(maf, pipeline_params("somatic"))
        germline, _ = apply_filters(maf, pipeline_params("germline"))
        kept = somatic.index[somatic[MAFIGATE_FILTER] == PASS]

        comparison = price_other_arm(maf, pipeline_params("germline"), kept)

        self.assertEqual(comparison.kept, len(kept))
        self.assertEqual(
            comparison.would_keep, int((germline[MAFIGATE_FILTER] == PASS).sum())
        )

    def test_a_file_the_other_arm_cannot_filter_is_priced_as_none(self):
        """A counterfactual nobody asked for must cost the numbers and nothing else."""
        from filters.arm_detection import price_other_arm
        from utils import read_maf

        maf = read_maf(str(GERMLINE_MAF))
        with patch(
            "filters.arm_detection.apply_filters", side_effect=ValueError("nope")
        ):
            self.assertIsNone(price_other_arm(maf, {}, maf.index[:3]))

    def test_identical_reports_are_neither_a_gain_nor_a_loss(self):
        from filters.arm_detection import ArmComparison

        same = ArmComparison(kept=1, would_keep=1, in_both=1)

        self.assertTrue(same.identical)
        self.assertEqual((same.gained, same.lost), (0, 0))


class TestWhatTheNoticeSays(unittest.TestCase):
    """The claims, read off what renders — issue #79's standard, made structural.

    Mutation-checked by dropping the absent-column half of the evidence sentence, which the
    naming assertions below catch, and by returning the positive joiner from ``_absent``,
    which ``no `InterVar` and `RENOVO_Class``` fails.
    """

    def _sentence(self, path, current):
        from components.arm_notice import _consequence, _evidence_sentence
        from config.pipeline_params import pipeline_params
        from filters.arm_detection import OTHER_ARM, price_other_arm, read_arm_evidence
        from filters.variant_filters import MAFIGATE_FILTER, PASS, apply_filters
        from utils import read_maf

        maf = read_maf(str(path))
        evidence = read_arm_evidence(maf.columns)
        labelled, _ = apply_filters(maf, pipeline_params(current))
        kept = labelled.index[labelled[MAFIGATE_FILTER] == PASS]
        comparison = price_other_arm(maf, pipeline_params(OTHER_ARM[current]), kept)
        detected = evidence.detected
        return (
            _evidence_sentence(evidence, detected, current),
            _consequence(comparison, current, detected),
        )

    def test_the_evidence_names_what_is_there_and_what_is_not(self):
        evidence, _ = self._sentence(GERMLINE_MAF, "somatic")

        self.assertIn("`InterVar`", evidence)
        self.assertIn("`RENOVO_Class`", evidence)
        self.assertIn("`CancerVar`", evidence)

    def test_the_absent_half_reads_as_a_denial_of_both(self):
        """``no `A` and `B``` denies the pair while allowing either; the file carries neither."""
        evidence, _ = self._sentence(SOMATIC_MAF, "germline")

        self.assertIn("neither `InterVar` nor `RENOVO_Class`", evidence)

    def test_a_gaining_direction_says_different_set_and_gives_the_overlap(self):
        _, consequence = self._sentence(GERMLINE_MAF, "somatic")

        self.assertIn("different", consequence)
        self.assertIn("in both reports", consequence)

    def test_a_losing_direction_does_not_claim_extra_variants(self):
        _, consequence = self._sentence(SOMATIC_MAF, "germline")

        self.assertIn("only drop variants", consequence)
        self.assertNotIn("different", consequence)

    def test_an_unpriced_notice_makes_no_claim_about_a_report_that_may_not_exist(self):
        """The refused case reaches this branch, and there is no report "below" it."""
        from components.arm_notice import _consequence

        consequence = _consequence(None, "somatic", "germline")

        self.assertNotIn("report below", consequence)
        self.assertIn("germline", consequence)

    def test_identical_reports_are_not_described_as_a_loss(self):
        """``plan_fills``'s cry-wolf argument, in a different place.

        7 of 8 sampled small somatic files kept one variant on either arm, so this state is
        common rather than a curiosity. Mutation-checked by removing the ``identical`` branch,
        which then reports "only drop variants … they add nothing" over a report that lost
        nothing.
        """
        from components.arm_notice import _consequence
        from filters.arm_detection import ArmComparison

        consequence = _consequence(ArmComparison(1, 1, 1), "germline", "somatic")

        self.assertIn("same 1 variant", consequence)
        self.assertNotIn("variants,", consequence.split("so your report")[0])
        self.assertNotIn("only drop", consequence)

    def test_no_internal_parameter_names_reach_the_copy(self):
        """The Style section's rule. Column names are the #79 exception; this is not."""
        from components.arm_notice import _consequence
        from filters.arm_detection import ArmComparison

        for consequence in (
            self._sentence(GERMLINE_MAF, "somatic")[1],
            self._sentence(SOMATIC_MAF, "germline")[1],
            _consequence(None, "somatic", "germline"),
            _consequence(ArmComparison(1, 1, 1), "germline", "somatic"),
        ):
            self.assertNotIn("sample_type", consequence)
            self.assertNotIn("filter_params", consequence)
            self.assertNotIn("#13", consequence)


class TestTheJumpIsWithheld(unittest.TestCase):
    """A report cut on the wrong arm is not one to drop the user into.

    Driven through ``_load_pending_upload`` rather than asserted on a flag, because the
    withholding and the token stamping sit in the same branch and the interesting failure is
    getting one of them wrong.
    """

    def _load(self, maf_path, arm):
        from components.sidebar import PENDING_UPLOAD_KEY
        from config.pipeline_params import pipeline_params
        from page_modules import data_loading
        from utils import read_maf

        upload = Mock(name="upload")
        upload.name = maf_path.name
        upload.size = 1234

        state = FakeSessionState(filter_params=pipeline_params(arm))
        state[PENDING_UPLOAD_KEY] = upload

        with patch.object(data_loading, "st") as page_st, patch.object(
            data_loading, "read_maf", return_value=read_maf(str(maf_path))
        ):
            page_st.session_state = state
            page_st.spinner.return_value.__enter__ = Mock()
            page_st.spinner.return_value.__exit__ = Mock(return_value=False)
            data_loading._load_pending_upload()

        return state

    def test_a_matching_file_still_opens_its_report(self):
        state = self._load(SOMATIC_MAF, "somatic")

        self.assertTrue(state.get("jump_to_results"))

    def test_a_mismatched_file_lands_on_the_filter_run_section(self):
        """Mutation-checked by dropping the ``file_and_arm_disagree`` guard, which jumps."""
        state = self._load(GERMLINE_MAF, "somatic")

        self.assertFalse(state.get("jump_to_results"))

    def test_the_upload_token_is_stamped_either_way(self):
        """It guards the re-read, not the navigation.

        Mutation-checked by moving the stamp inside the jump branch: the file is then re-read
        and re-filtered on every render, for exactly the files this notice is about.
        """
        state = self._load(GERMLINE_MAF, "somatic")

        self.assertEqual(state.get("last_upload_token"), (GERMLINE_MAF.name, 1234))


class TestTheSwitch(unittest.TestCase):
    """Re-seed, re-filter, grant the withheld jump — and say so in issue #133's words."""

    def _switch(self, arm_before, switch_to, customise=None):
        from config.pipeline_params import pipeline_params
        from page_modules import data_loading, parameter_config
        from utils import read_maf

        params = pipeline_params(arm_before)
        if customise:
            params.update(customise)
        state = FakeSessionState(
            filter_params=params, maf_data=read_maf(str(GERMLINE_MAF))
        )

        # Both modules share one session state, as in the app: the switch runs on the data
        # page and the replacement it routes through lives on the parameter page.
        with patch.object(data_loading, "st") as page_st, patch.object(
            parameter_config, "st"
        ) as param_st, patch.object(
            parameter_config, "park_confirmation"
        ) as parked:
            page_st.session_state = state
            param_st.session_state = state
            page_st.spinner.return_value.__enter__ = Mock()
            page_st.spinner.return_value.__exit__ = Mock(return_value=False)
            param_st.rerun.side_effect = RuntimeError("rerun")
            with self.assertRaises(RuntimeError):
                data_loading._switch_arm_and_refilter(switch_to)

        return state, parked

    def test_the_arm_is_re_seeded_from_that_arms_contract(self):
        """The same call ``price_other_arm`` was handed, so the count quoted is delivered."""
        from config.pipeline_params import pipeline_params

        state, _ = self._switch("somatic", "germline")

        self.assertEqual(state["filter_params"], pipeline_params("germline"))

    def test_edits_on_the_old_arm_do_not_survive(self):
        state, _ = self._switch("somatic", "germline", customise={"min_depth": 999})

        self.assertNotEqual(state["filter_params"].get("min_depth"), 999)

    def test_the_report_is_recut_and_the_jump_granted(self):
        state, _ = self._switch("somatic", "germline")

        self.assertIsNotNone(state.get("filtered_data"))
        self.assertTrue(state.get("jump_to_results"))

    def test_the_confirmation_is_parked_naming_the_arm(self):
        """Parked, not drawn: the ``st.rerun`` beneath it discards the frame (issue #133)."""
        _, parked = self._switch("somatic", "germline")

        parked.assert_called_once()
        self.assertEqual(parked.call_args.kwargs.get("arm"), "germline")
        self.assertIn("germline", parked.call_args.args[0])


class TestPressingTheButtonReachesTheSwitch(unittest.TestCase):
    """The one link ``AppTest`` cannot drive, so it is asserted here instead.

    Clicking the button through the whole app is not available: replaying a run trips
    ``element_tree``'s reconstruction of the section control's ``format_func`` — the harness
    limitation issue #92 recorded, not an app defect. What that leaves untested is the
    shortest and most breakable link in the chain, since ``app-load-check`` proves only that
    the button is *drawn*. A button wired to nothing looks identical from outside.
    """

    def test_the_button_calls_what_it_was_handed(self):
        from components import arm_notice
        from filters.arm_detection import read_arm_evidence
        from utils import read_maf

        evidence = read_arm_evidence(read_maf(str(GERMLINE_MAF)).columns)
        pressed = []

        with patch.object(arm_notice, "st") as notice_st:
            notice_st.button.return_value = True
            arm_notice.render_mismatch_notice(
                evidence=evidence,
                current_arm="somatic",
                comparison=None,
                customised=False,
                on_switch=pressed.append,
            )

        self.assertEqual(pressed, ["germline"])

    def test_an_unpressed_button_switches_nothing(self):
        """Mutation-checked by calling ``on_switch`` unconditionally, which fires on render."""
        from components import arm_notice
        from filters.arm_detection import read_arm_evidence
        from utils import read_maf

        evidence = read_arm_evidence(read_maf(str(GERMLINE_MAF)).columns)
        pressed = []

        with patch.object(arm_notice, "st") as notice_st:
            notice_st.button.return_value = False
            arm_notice.render_mismatch_notice(
                evidence=evidence,
                current_arm="somatic",
                comparison=None,
                customised=False,
                on_switch=pressed.append,
            )

        self.assertEqual(pressed, [])

    def test_the_page_hands_the_button_the_real_switch(self):
        """The other half: the notice calls what it is given, and this is what it is given."""
        from config.pipeline_params import pipeline_params
        from page_modules import data_loading
        from utils import read_maf

        state = FakeSessionState(
            filter_params=pipeline_params("somatic"),
            maf_data=read_maf(str(GERMLINE_MAF)),
            filtered_data=None,
        )

        with patch.object(data_loading, "st") as page_st, patch.object(
            data_loading, "render_mismatch_notice"
        ) as notice:
            page_st.session_state = state
            data_loading._show_arm_notice()

        notice.assert_called_once()
        self.assertIs(
            notice.call_args.kwargs["on_switch"], data_loading._switch_arm_and_refilter
        )


class TestWhatThePageDraws(unittest.TestCase):
    """The wiring, run rather than read — the half a component test cannot see."""

    _PAGE = """
import os, sys
sys.path.insert(0, {app!r})

import streamlit as st
from config.pipeline_params import pipeline_params
from page_modules import data_loading
from utils import read_maf

st.session_state.filter_params = pipeline_params({arm!r})
st.session_state.maf_source_name = os.path.basename({maf!r})
st.session_state.maf_data = read_maf({maf!r})
data_loading.apply_filters_to_data(show_messages=False)

data_loading.show_data_loading_page()
"""

    def _render(self, maf, arm):
        pytest.importorskip("streamlit", reason="streamlit not installed")
        from streamlit.testing.v1 import AppTest

        app = AppTest.from_string(
            self._PAGE.format(app=str(STREAMLIT_APP), maf=str(maf), arm=arm)
        )
        app.run(timeout=300)
        assert not app.exception, [str(e.value) for e in app.exception]
        return app

    def _warnings(self, app, clause):
        return [w.value for w in app.warning if clause in w.value]

    def test_a_mismatched_file_draws_the_notice_once(self):
        app = self._render(GERMLINE_MAF, "somatic")

        self.assertEqual(len(self._warnings(app, MISMATCH)), 1)

    def test_a_matching_file_draws_no_notice(self):
        """Mutation-checked by dropping the ``disagrees_with`` guard in ``_show_arm_notice``."""
        app = self._render(SOMATIC_MAF, "somatic")

        self.assertEqual(self._warnings(app, MISMATCH), [])
        self.assertEqual(self._warnings(app, AMBIGUOUS), [])

    def test_the_notice_offers_the_switch(self):
        app = self._render(GERMLINE_MAF, "somatic")

        labels = [button.label for button in app.button]
        self.assertIn("🔄 Switch to germline and re-filter", labels)

    def test_a_pristine_session_is_not_told_it_will_lose_settings(self):
        """Crying wolf on a user with nothing to lose is how the next real warning dies."""
        app = self._render(GERMLINE_MAF, "somatic")

        self.assertEqual(
            [c.value for c in app.caption if "will be lost" in c.value], []
        )

    def test_a_file_carrying_neither_marker_draws_no_mismatch_notice(self):
        """Silence where the filter already shouts — and the shout is asserted, not assumed.

        Mutation-checked by making ``read_arm_evidence`` fall back to the selected arm, which
        turns every unplaceable file into a mismatch.
        """
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            path = _maf_with_columns_removed(Path(tmp), ["CancerVar"])
            app = self._render(path, "somatic")

            self.assertEqual(self._warnings(app, MISMATCH), [])
            self.assertEqual(self._warnings(app, AMBIGUOUS), [])
            self.assertTrue(
                [e.value for e in app.error if "PATHOGENIC RETENTION" in e.value],
                "the escalated warning is what makes this silence safe; if it stops "
                "firing, this file says nothing at all",
            )

    def test_a_file_carrying_both_markers_says_so_and_offers_no_switch(self):
        """The one unplaceable case nothing else in the app mentions."""
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            path = _maf_carrying_both_arms(Path(tmp))
            app = self._render(path, "somatic")

            self.assertEqual(len(self._warnings(app, AMBIGUOUS)), 1)
            self.assertEqual(self._warnings(app, MISMATCH), [])
            self.assertEqual(
                [b.label for b in app.button if b.label.startswith("🔄 Switch")], []
            )


if __name__ == "__main__":
    unittest.main()
