"""Opening a MAF from the sidebar, and what happens to the last file's notes (issue #64).

Loading used to happen in exactly one place: an ``st.file_uploader`` on the data page's load
section. Opening a different sample therefore meant navigating to that page and finding the right
section first, which the walkthrough in #53 reported as friction.

The app now has **one** chooser and it lives in the sidebar, beside the block that says which
file is open. Two claims here, and they are different kinds of claim:

*Where the chooser is, and what choosing does.* It is drawn plainly when nothing is open —
then it is the only thing to do — and folded into an expander once a file is open, so the
sidebar's first job, saying what you are looking at, keeps the top of the column. Choosing a
file parks it in ``PENDING_UPLOAD_KEY`` and routes to the data page: the page is where a MAF is
read, validated and filtered, and where the account of that has room to be read.

*What happens to notes when the file underneath them changes.* Notes and custom annotation
values are keyed by variant identity, not by file, so without a rule they re-attach to any
matching variant in the next sample — one patient's note surfacing on another's report. They
are cleared, the columns are kept, and the user is told. The tests that matter most here are
the two negative ones: re-opening the *same* file keeps your notes, and a file that could not
be read takes nothing with it.

Whether a note should survive the *session* is issue #67's question, not this module's.

A unit module: nothing here has a pipeline counterpart — the pipeline has no sidebar — and
none of it needs ``bin/``.
"""

from __future__ import annotations

import os
import sys
import unittest
from unittest.mock import Mock, patch

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


class FakeSessionState(dict):
    """Streamlit's session_state: dict access plus attributes. See ``test_components.py``."""

    def __getattr__(self, key):
        try:
            return self[key]
        except KeyError as exc:
            raise AttributeError(key) from exc

    def __setattr__(self, key, value):
        self[key] = value


class TestWhereTheChooserIsDrawn(unittest.TestCase):
    """One chooser, placed by whether there is already a file to keep in sight."""

    def _render(self, mock_st, **session):
        from components.sidebar import render_load_status

        mock_st.session_state = FakeSessionState(**session)
        slot = Mock()
        slot.button.return_value = False
        render_load_status(slot)
        return slot

    @patch("components.sidebar.st")
    def test_the_empty_state_draws_the_chooser_in_the_open(self, mock_st):
        slot = Mock()
        slot.button.return_value = False
        mock_st.session_state = FakeSessionState(
            current_page="parameter_config", maf_data=None
        )
        from components.sidebar import render_load_status

        render_load_status(slot)

        self.assertEqual(slot.file_uploader.call_count, 1)
        self.assertFalse(
            slot.expander.called, "nothing to keep in sight, so nothing to fold behind"
        )

    @patch("components.sidebar.st")
    def test_an_open_file_folds_the_chooser_behind_an_expander(self, mock_st):
        import pandas as pd

        slot = self._render(
            mock_st,
            current_page="parameter_config",
            maf_source_name="sample_42.maf",
            maf_data=pd.DataFrame({"Hugo_Symbol": ["TP53"] * 10}),
            filtered_data=None,
            filter_params={"sample_type": "somatic"},
        )

        self.assertFalse(
            slot.file_uploader.called,
            "the chooser was drawn straight into the column, above the nav, pushing the "
            "answer to *which file* down every page",
        )
        label = str(slot.expander.call_args.args[0])
        self.assertIn("different file", label)
        self.assertEqual(slot.expander.return_value.file_uploader.call_count, 1)


class TestWhatChoosingAFileDoes(unittest.TestCase):
    """The wiring, exercised through the callback the chooser is actually given."""

    def _on_change(self, mock_st, **session):
        """Render the empty state and return the chooser's own ``on_change``."""
        from components.sidebar import render_load_status

        state = FakeSessionState(current_page="help", maf_data=None, **session)
        mock_st.session_state = state
        slot = Mock()
        slot.button.return_value = False
        render_load_status(slot)

        kwargs = slot.file_uploader.call_args.kwargs
        return state, kwargs

    @patch("components.sidebar.st")
    def test_a_chosen_file_is_handed_to_the_data_page(self, mock_st):
        from components.sidebar import PENDING_UPLOAD_KEY, UPLOAD_CHOOSER_KEY

        state, kwargs = self._on_change(mock_st)
        chosen = object()
        state[UPLOAD_CHOOSER_KEY] = chosen

        kwargs["on_change"]()

        self.assertIs(state.get(PENDING_UPLOAD_KEY), chosen)
        self.assertEqual(
            state.get("current_page"),
            "data_loading",
            "choosing a file from another page left the user on that page, with the load "
            "and everything it has to say happening somewhere they are not looking",
        )

    @patch("components.sidebar.st")
    def test_clearing_the_chooser_withdraws_the_file_and_closes_nothing(self, mock_st):
        """The ✕ has to be a way out, without being a way to lose your session.

        A file the page cannot read is offered until it is withdrawn, so the ✕ is the escape
        hatch from a banner that would otherwise redraw on every render. What it must not do
        is close the file you are reading: the open MAF and its notes are not what was cleared.
        """
        from components.sidebar import PENDING_UPLOAD_KEY, UPLOAD_CHOOSER_KEY

        state, kwargs = self._on_change(mock_st)
        state[PENDING_UPLOAD_KEY] = object()
        state["maf_data"] = "the file being read"
        state[UPLOAD_CHOOSER_KEY] = None

        kwargs["on_change"]()

        self.assertIsNone(state.get(PENDING_UPLOAD_KEY))
        self.assertEqual(state.get("maf_data"), "the file being read")
        self.assertEqual(state.get("current_page"), "help")

    @patch("components.sidebar.st")
    def test_the_chooser_is_keyed_so_its_value_can_be_read_back(self, mock_st):
        from components.sidebar import UPLOAD_CHOOSER_KEY

        _, kwargs = self._on_change(mock_st)

        self.assertEqual(kwargs["key"], UPLOAD_CHOOSER_KEY)
        self.assertEqual(kwargs["type"], ["txt", "tsv", "maf"])


class TestNotesWhenTheFileChanges(unittest.TestCase):
    """Notes survive a change of file, and that is the decision rather than the default.

    Issue #64 cleared them here. It keyed on the same fact this class does — ``_variant_key``
    identifies a variant by gene, chromosome, position and ref>alt, so a note written on one
    sample re-attaches to any matching variant in the next — and read it as one sample's note
    leaking onto another sample's report.

    Issue #67 then answered the question #64's docstring had referred to it, and answered it
    the other way: a note is a statement *about the variant*, intended to be read by everyone
    at the institute once there is a server to hold it. So the surfacing is the feature, and
    what protects the patient is the variant dialog's line asking for no patient identifiers,
    not a shorter life for the note.

    This class is therefore a **guard against silent reversal**. Deleting a clearing rule
    leaves nothing behind, and the next reader of ``_variant_key`` will reach for exactly the
    rule #64 reached for. Anyone who re-adds it fails here and has to reopen #67 first.
    """

    def _load_a_different_file(self, **session):
        """Drive the real shared load tail, not a helper that only this test calls.

        #64's version of this class called the clearing function directly, which is what let
        the rule and its tests agree with each other while the app's three load doors were
        free to disagree. Going through ``_open_the_file_just_read`` means the claim is about
        the app: whichever door a file arrives by, what the user wrote is still there.
        """
        import pandas as pd

        with patch("page_modules.data_loading.apply_filters_to_data") as apply_filters, patch(
            "page_modules.data_loading.validate_required_columns"
        ) as validate, patch("components.sidebar.st") as mock_st, patch(
            "page_modules.data_loading.st"
        ) as page_st:
            from page_modules.data_loading import _open_the_file_just_read

            validate.return_value = True
            apply_filters.return_value = True

            state = FakeSessionState(maf_source_name="sample_a.maf", **session)
            page_st.session_state = state
            mock_st.session_state = state

            _open_the_file_just_read(pd.DataFrame({"Hugo_Symbol": ["BRAF"]}), "sample_b.maf")
            return state, page_st, mock_st

    def test_a_different_file_keeps_the_notes(self):
        notes = {"TP53:17:7577120:C>T": "hotspot, class 5, seen before"}
        state, _, _ = self._load_a_different_file(
            variant_notes=dict(notes), custom_annotations={}
        )

        self.assertEqual(state["variant_notes"], notes)

    def test_a_different_file_keeps_what_was_written_in_the_users_own_columns(self):
        """The values, not just the column names — #64 kept the columns and emptied them."""
        annotations = {"Reviewed by": {"TP53:17:7577120:C>T": "MR"}}
        state, _, _ = self._load_a_different_file(
            variant_notes={}, custom_annotations={"Reviewed by": dict(annotations["Reviewed by"])}
        )

        self.assertEqual(state["custom_annotations"], annotations)

    def test_the_new_file_is_still_the_one_the_sidebar_names(self):
        """The rule's removal must not take the name recording with it."""
        state, _, _ = self._load_a_different_file(variant_notes={}, custom_annotations={})

        self.assertEqual(state["maf_source_name"], "sample_b.maf")

    def test_nothing_tells_the_user_their_notes_were_cleared(self):
        """Because nothing clears them. A leftover banner would be worse than none.

        #64's banner named the outgoing file and counted what went. Left in place over a rule
        that no longer runs, it would report a loss that did not happen.
        """
        _, page_st, mock_st = self._load_a_different_file(
            variant_notes={"TP53:17:7577120:C>T": "hotspot"},
            custom_annotations={"Reviewed by": {"KRAS:12:25398284:C>A": "MR"}},
        )

        for streamlit in (page_st, mock_st):
            said = " ".join(str(call) for call in streamlit.info.call_args_list)
            self.assertNotIn("cleared", said)

    def test_no_one_has_re_added_a_clearing_rule_beside_the_notes(self):
        """The reversal this class exists to catch, named at its most likely address.

        A future clearing rule would most plausibly arrive as a function in
        ``components/variant_table.py``, where the notes and ``_variant_key`` live and where
        #64 put it — in ``ui_components.py``, which #76 dissolved into the seven modules
        this one is the notes' half of. The behavioural tests above catch one wired into
        the load tail; this catches one added and wired in somewhere they do not reach.
        """
        from components import variant_table

        named = [name for name in dir(variant_table) if "discard_notes" in name]

        self.assertEqual(named, [], "see issue #67 before re-adding a notes-clearing rule")


class TestAFileThatCannotBeRead(unittest.TestCase):
    """A mis-picked file costs a banner and nothing else.

    The read happens before anything is written, so a ``.txt`` that is not a MAF — or a
    truncated download — leaves the session it arrived into exactly as it was: same open file,
    same report, same notes. And the file is *withdrawn* rather than left on offer, because the
    page runs on every render and a file that raised once will raise identically forever.
    """

    def _fail_to_read(self, mock_st, page_st, read_maf, read_fallback, **session):
        from components.sidebar import PENDING_UPLOAD_KEY
        from page_modules.data_loading import _load_pending_upload

        read_maf.side_effect = ValueError("not a MAF")
        read_fallback.side_effect = ValueError("not a MAF either")

        upload = Mock(name="upload")
        upload.name = "not_a_maf.txt"
        state = FakeSessionState(**session)
        state[PENDING_UPLOAD_KEY] = upload
        # One session state behind both modules, as in the app: the page reads the chooser's
        # key and the notes rule lives with the notes.
        mock_st.session_state = state
        page_st.session_state = state

        _load_pending_upload()
        return state

    @patch("page_modules.data_loading.read_maf_without_comment_lines")
    @patch("page_modules.data_loading.read_maf")
    @patch("components.sidebar.st")
    @patch("page_modules.data_loading.st")
    def test_it_leaves_the_open_file_and_its_notes_alone(
        self, page_st, mock_st, read_maf, read_fallback
    ):
        notes = {"TP53:17:7577120:C>T": "discussed at MTB"}
        state = self._fail_to_read(
            mock_st,
            page_st,
            read_maf,
            read_fallback,
            maf_source_name="sample_a.maf",
            maf_data="the file being read",
            variant_notes=dict(notes),
            custom_annotations={},
        )

        self.assertEqual(state["variant_notes"], notes)
        self.assertEqual(
            state["maf_data"],
            "the file being read",
            "a file the app could not even read closed the file the user was reading",
        )
        self.assertEqual(
            state["maf_source_name"],
            "sample_a.maf",
            "the sidebar would name a file the app never managed to open",
        )
        self.assertTrue(page_st.error.called)

    @patch("page_modules.data_loading.read_maf_without_comment_lines")
    @patch("page_modules.data_loading.read_maf")
    @patch("components.sidebar.st")
    @patch("page_modules.data_loading.st")
    def test_it_is_withdrawn_rather_than_retried_forever(
        self, page_st, mock_st, read_maf, read_fallback
    ):
        from components.sidebar import PENDING_UPLOAD_KEY

        state = self._fail_to_read(
            mock_st, page_st, read_maf, read_fallback, variant_notes={}
        )

        self.assertNotIn(
            PENDING_UPLOAD_KEY,
            state,
            "the file stays on offer, so its banner redraws on every render of the page",
        )


class TestWhatThePreviousFileLeavesBehind(unittest.TestCase):
    """Opening a file drops what the *last* file's filter run left in the session.

    ``failed_data`` and ``data_params_hash`` both outlived the file they described. The stamp is
    the one that bites: the sidebar compares it with the parameters now in force to decide
    whether to call your results *current*, so a new MAF whose own filter run then failed left
    the old run's stamp standing over no report at all.
    """

    @patch("page_modules.data_loading.apply_filters_to_data")
    @patch("page_modules.data_loading.validate_required_columns")
    @patch("components.sidebar.st")
    @patch("page_modules.data_loading.st")
    def test_a_new_file_keeps_none_of_the_last_run(
        self, page_st, mock_st, validate, apply_filters
    ):
        import pandas as pd

        from page_modules.data_loading import _open_the_file_just_read

        validate.return_value = True
        apply_filters.return_value = False  # read and validated, but the filters raised

        state = FakeSessionState(
            maf_source_name="sample_a.maf",
            filtered_data=pd.DataFrame({"Hugo_Symbol": ["TP53"]}),
            failed_data=pd.DataFrame({"Hugo_Symbol": ["KRAS"] * 3}),
            data_params_hash="the stamp of the last cut",
            variant_notes={},
            custom_annotations={},
        )
        page_st.session_state = state
        mock_st.session_state = state

        produced = _open_the_file_just_read(pd.DataFrame({"Hugo_Symbol": ["BRAF"]}), "b.maf")

        self.assertFalse(produced)
        self.assertIsNone(state["filtered_data"])
        self.assertTrue(state["failed_data"].empty)
        self.assertIsNone(state["data_params_hash"])
        self.assertEqual(state["maf_source_name"], "b.maf")


if __name__ == "__main__":
    unittest.main()
