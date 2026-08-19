"""The list of files you have opened: what is written, when, and what the app says about it.

Issue #158. Three claims, and they are different kinds of claim.

*What reaches the disk.* An entry is a name and a time, in a stamped envelope, capped, and
deduplicated. The stamp is issue #40's lesson applied to a second file in ``~/.mafigate``:
an unstamped document there was silently restored as though it were current. The cap and the
name-only rule are the map's, and they are the reason the feature was acceptable at all —
153 of 205 real MAF basenames carry a sample code, so this file names patients, and it
**outlives the MAF it names**.

*When it is written.* Once, in the tail all three load doors share, after the read has
succeeded and after the column check has accepted the file. So the list means *files you have
looked at*: a MAF that could not be read and a MAF refused for a missing column both leave no
trace, while a file whose *filters* refuse is recorded, because it is open and the sidebar
names it.

*What the app promises.* The sentence beside the list has to state the cap and say that a
name here re-opens nothing, in words a clinician reads — the map's Style section — and the
clear control has to clear this list without taking the user's filter settings with it.

Every test redirects the store into a temporary directory. Left unredirected they would
delete and rewrite the developer's own ``~/.mafigate``.

A unit module: the pipeline has no sidebar and no recents list, so nothing here has a
counterpart in ``bin/``.
"""

from __future__ import annotations

import os
import sys
import tempfile
import unittest
from datetime import datetime
from pathlib import Path
from unittest.mock import Mock, patch

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

#: A fabricated sample-code-shaped MAF name, and the reason it is a constant rather than a
#: literal repeated down the file.
#:
#: The shape is the point. Most of what this module asserts is about a *sample code* reaching
#: the disk and reaching the sidebar — that is why the feature needed a decision at all — and a
#: name like ``sample.maf`` would make those assertions read as though the concern were
#: hypothetical. The generic names are still used where the shape is irrelevant (``a.maf``,
#: ``first.maf``, ``sample_3.maf``).
#:
#: These used to be a real case identifier out of the cohort, in nineteen places. That is a
#: clinical identifier in a file that ships publicly, and ``tools/export_public.py``'s scan
#: gate refuses to export it — correctly (issue #247). One constant, so the next person who
#: needs a twentieth occurrence reaches for this instead of pasting a real one, and so nothing
#: has to be de-identified twice.
SAMPLE_NAME = "CASE01ABCD-1N.maf"

#: A second one, for the tests that need two entries to be distinguishable.
OTHER_SAMPLE_NAME = "CASE02EFGH-1T.maf"


class FakeSessionState(dict):
    """Streamlit's session_state: dict access plus attributes. See ``test_components.py``."""

    def __getattr__(self, key):
        try:
            return self[key]
        except KeyError as exc:
            raise AttributeError(key) from exc

    def __setattr__(self, key, value):
        self[key] = value


class RedirectedHistory(unittest.TestCase):
    """A base class whose only job is to keep the suite out of the real home directory.

    ``config.file_history`` writes to ``~/.mafigate`` and :func:`clear_history` unlinks. Both
    module constants are redirected, not just the file: :func:`record_open_file` calls
    ``mkdir`` on the directory, so a test that redirected only the file would still create a
    directory in the developer's home on a machine that had none.
    """

    def setUp(self):
        from config import file_history

        self._tmp = tempfile.TemporaryDirectory()
        self.home = Path(self._tmp.name) / ".mafigate"
        self._patches = [
            patch.object(file_history, "HISTORY_DIR", self.home),
            patch.object(file_history, "HISTORY_FILE", self.home / "opened_files.json"),
        ]
        for entry in self._patches:
            entry.start()
        self.addCleanup(self._tmp.cleanup)
        for entry in self._patches:
            self.addCleanup(entry.stop)

    @property
    def store(self) -> Path:
        from config import file_history

        return file_history.HISTORY_FILE


class TestWhatReachesTheDisk(RedirectedHistory):
    """A name, a time, a stamp, and a cap — and nothing else about the file."""

    def test_a_recorded_file_can_be_read_back(self):
        from config.file_history import NAME_KEY, read_history, record_open_file

        record_open_file(SAMPLE_NAME)

        history = read_history()
        self.assertEqual([entry[NAME_KEY] for entry in history], [SAMPLE_NAME])

    def test_the_document_carries_the_format_stamp(self):
        """Issue #40's lesson, applied to the second file in ``~/.mafigate``.

        Without a stamp this file is what the parameter cache was: something a later MAFigate
        cannot tell apart from a file five formats old, leaving it no option but to discard
        every one of them.
        """
        import json

        from config.file_history import (
            HISTORY_SCHEMA_VERSION,
            SCHEMA_VERSION_KEY,
            record_open_file,
        )

        record_open_file("sample.maf")

        document = json.loads(self.store.read_text())
        self.assertEqual(document[SCHEMA_VERSION_KEY], HISTORY_SCHEMA_VERSION)

    def test_an_entry_holds_the_name_and_the_time_and_nothing_else(self):
        """The map's names-only rule, asserted as an exact key set.

        A variant count or an arm added here would be a fact about a patient's file living on
        in a file that outlives it, which is more at rest than the map weighed and accepted.
        Widening this is a decision, so it fails here first.
        """
        from config.file_history import NAME_KEY, OPENED_KEY, read_history, record_open_file

        record_open_file("sample.maf")

        (entry,) = read_history()
        self.assertEqual(set(entry), {NAME_KEY, OPENED_KEY})

    def test_only_the_basename_is_stored_even_when_the_caller_has_a_path(self):
        """The OS "Open With" door does have a real path. It is still not stored.

        A path is a route back to the file, which this list deliberately is not, and a second
        statement about where a patient's data sits on this disk.

        The path is fabricated, and it has to be one: this test used to pass a real cloud-drive
        path under a real home directory, and both are things the publication scan gate refuses
        to export. What the assertion needs is a path with *directories* in it, which any
        fabricated one supplies.
        """
        from config.file_history import NAME_KEY, read_history, record_open_file

        record_open_file(f"/data/clinical/cohort_share/{SAMPLE_NAME}")

        (entry,) = read_history()
        self.assertEqual(entry[NAME_KEY], SAMPLE_NAME)
        self.assertNotIn("cohort_share", self.store.read_text())

    def test_the_newest_file_is_first(self):
        from config.file_history import NAME_KEY, read_history, record_open_file

        record_open_file("first.maf")
        record_open_file("second.maf")

        self.assertEqual(
            [entry[NAME_KEY] for entry in read_history()], ["second.maf", "first.maf"]
        )

    def test_re_opening_a_file_moves_its_entry_rather_than_adding_a_second(self):
        """A recents list answers *which files*, not *how often*."""
        from config.file_history import NAME_KEY, read_history, record_open_file

        record_open_file("a.maf")
        record_open_file("b.maf")
        record_open_file("a.maf")

        self.assertEqual([entry[NAME_KEY] for entry in read_history()], ["a.maf", "b.maf"])

    def test_the_list_is_capped(self):
        """The cap is what keeps this a recents list rather than an archive of patient names."""
        from config.file_history import HISTORY_LIMIT, NAME_KEY, read_history, record_open_file

        for index in range(HISTORY_LIMIT + 5):
            record_open_file(f"sample_{index}.maf")

        history = read_history()
        self.assertEqual(len(history), HISTORY_LIMIT)
        self.assertEqual(history[0][NAME_KEY], f"sample_{HISTORY_LIMIT + 4}.maf")
        self.assertNotIn("sample_0.maf", self.store.read_text())

    def test_the_cap_holds_on_the_disk_and_not_only_on_the_way_out(self):
        """Counted in the stored document, because :func:`read_history` slices too.

        Measured: with the write's own cap removed and every other assertion above still
        green, the file held eleven names — the ten it read back plus one. A cap enforced
        only on the read is not a cap on what is at rest, and what is at rest is the whole
        reason this list has a cap.
        """
        import json

        from config.file_history import FILES_KEY, HISTORY_LIMIT, record_open_file

        for index in range(HISTORY_LIMIT + 5):
            record_open_file(f"sample_{index}.maf")

        document = json.loads(self.store.read_text())
        self.assertEqual(len(document[FILES_KEY]), HISTORY_LIMIT)

    def test_the_time_is_the_time_the_file_was_opened(self):
        from config.file_history import OPENED_KEY, read_history, record_open_file

        record_open_file("sample.maf", now=datetime(2026, 8, 16, 9, 20, 5))

        (entry,) = read_history()
        self.assertEqual(entry[OPENED_KEY], "2026-08-16T09:20:05")


class TestAStoreThisVersionCannotRead(RedirectedHistory):
    """Every unreadable shape reads as *no history*, in silence.

    Silence is the difference from the parameter cache, and it is deliberate. A cache is the
    user's only copy of settings they configured, so it is moved aside and announced. A
    recents list is rebuilt by opening files, so there is nothing to rescue — and a banner
    about the internals of a recents list would be the app talking about itself.
    """

    def test_no_file_is_an_empty_history(self):
        from config.file_history import read_history

        self.assertEqual(read_history(), [])

    def test_an_unstamped_document_is_not_read(self):
        import json

        from config.file_history import FILES_KEY, read_history

        self.home.mkdir(parents=True)
        self.store.write_text(json.dumps({FILES_KEY: [{"name": "old.maf"}]}))

        self.assertEqual(read_history(), [])

    def test_a_document_from_a_newer_format_is_not_read(self):
        import json

        from config.file_history import (
            FILES_KEY,
            HISTORY_SCHEMA_VERSION,
            SCHEMA_VERSION_KEY,
            read_history,
        )

        self.home.mkdir(parents=True)
        self.store.write_text(
            json.dumps(
                {
                    SCHEMA_VERSION_KEY: HISTORY_SCHEMA_VERSION + 1,
                    FILES_KEY: [{"name": "future.maf"}],
                }
            )
        )

        self.assertEqual(read_history(), [])

    def test_a_file_that_does_not_parse_costs_the_list_and_not_the_session(self):
        from config.file_history import read_history

        self.home.mkdir(parents=True)
        self.store.write_text("{ this was a truncated write")

        self.assertEqual(read_history(), [])

    def test_an_entry_with_no_name_is_dropped_rather_than_drawn_blank(self):
        """Checked per entry, because the stamp says who wrote the file, not that every row
        in it survived the write."""
        import json

        from config.file_history import (
            FILES_KEY,
            HISTORY_SCHEMA_VERSION,
            NAME_KEY,
            SCHEMA_VERSION_KEY,
            read_history,
        )

        self.home.mkdir(parents=True)
        self.store.write_text(
            json.dumps(
                {
                    SCHEMA_VERSION_KEY: HISTORY_SCHEMA_VERSION,
                    FILES_KEY: [{"opened": "2026-08-16T09:20:05"}, {NAME_KEY: "real.maf"}],
                }
            )
        )

        self.assertEqual([entry[NAME_KEY] for entry in read_history()], ["real.maf"])

    def test_a_home_that_cannot_be_written_does_not_raise(self):
        """The user keeps the file they just opened; they lose only the recents list."""
        from config import file_history

        nowhere = Path("/proc/mafigate-nowhere")
        with patch.object(file_history, "HISTORY_DIR", nowhere), patch.object(
            file_history, "HISTORY_FILE", nowhere / "opened_files.json"
        ):
            file_history.record_open_file("sample.maf")


class TestClearing(RedirectedHistory):
    def test_clearing_removes_the_names_from_the_disk(self):
        """Deleted, not moved aside. A user clearing this list is asking for these names to
        stop being on their computer, so keeping a superseded copy beside it — which is right
        for a parameter cache — would take back what the control promises."""
        from config.file_history import clear_history, read_history, record_open_file

        record_open_file(SAMPLE_NAME)
        clear_history()

        self.assertEqual(read_history(), [])
        self.assertFalse(any(self.home.glob("opened_files*")))

    def test_clearing_an_empty_history_is_not_an_error(self):
        from config.file_history import clear_history

        self.assertTrue(clear_history())


class TestHowTheTimeReads(unittest.TestCase):
    """Relative for today, absolute beyond it — and nothing at all when it cannot be read."""

    def test_today_is_relative(self):
        from config.file_history import format_opened

        self.assertEqual(
            format_opened("2026-08-17T14:02:00", now=datetime(2026, 8, 17, 18, 0, 0)),
            "today 14:02",
        )

    def test_an_earlier_day_is_absolute(self):
        from config.file_history import format_opened

        self.assertEqual(
            format_opened("2026-08-16T09:20:00", now=datetime(2026, 8, 17, 18, 0, 0)),
            "16 Aug 09:20",
        )

    def test_a_stamp_that_cannot_be_read_leaves_the_name_standing_alone(self):
        """Rather than "Unknown" beside it, which is the app reporting on its own storage."""
        from config.file_history import format_opened

        for stamp in ("", "yesterday", None):
            self.assertEqual(format_opened(stamp), "")


class TestWhenAnEntryIsWritten(RedirectedHistory):
    """Once, in the tail all three load doors share — and only for a file that opened.

    Driven through ``_open_the_file_just_read`` rather than by calling the store, for the
    reason ``test_file_chooser.py`` gives about the same function: the claim is about the app,
    so it must hold whichever door the file arrived by.
    """

    def _open(self, name="sample_b.maf", accepted=True, filters_pass=True):
        import pandas as pd

        with patch("page_modules.data_loading.apply_filters_to_data") as apply_filters, patch(
            "page_modules.data_loading.validate_required_columns"
        ) as validate, patch("page_modules.data_loading.st") as page_st:
            from page_modules.data_loading import _open_the_file_just_read

            validate.return_value = accepted
            apply_filters.return_value = filters_pass
            page_st.session_state = FakeSessionState()

            opened = _open_the_file_just_read(pd.DataFrame({"Hugo_Symbol": ["BRAF"]}), name)
            return opened, page_st.session_state

    def test_a_file_that_opened_is_recorded(self):
        from config.file_history import NAME_KEY, read_history

        self._open(SAMPLE_NAME)

        self.assertEqual([entry[NAME_KEY] for entry in read_history()], [SAMPLE_NAME])

    def test_a_file_refused_for_a_missing_column_is_not_recorded(self):
        """The refusal has just unloaded it. Naming it in a list of files you have looked at
        would be the app's two surfaces disagreeing about whether it opened."""
        from config.file_history import read_history

        opened, state = self._open("incomplete.maf", accepted=False)

        self.assertFalse(opened)
        self.assertIsNone(state["maf_data"])
        self.assertEqual(read_history(), [])

    def test_a_file_whose_filters_refuse_is_still_recorded(self):
        """It is open, and the sidebar names it — the status block's middle state exists for
        exactly this. So it is a file the user has looked at."""
        from config.file_history import NAME_KEY, read_history

        opened, _ = self._open("awkward.maf", filters_pass=False)

        self.assertFalse(opened)
        self.assertEqual([entry[NAME_KEY] for entry in read_history()], ["awkward.maf"])

    def test_the_history_is_written_where_the_name_is_recorded(self):
        """One writer, at the one point three doors converge.

        A second call site added to a page would be a door that recorded twice, or one that
        recorded a file the shared tail then refused. Read statically, because the failure is
        the existence of the second call rather than anything it does.
        """
        import ast

        app = Path(__file__).resolve().parent.parent
        callers = []
        for path in app.rglob("*.py"):
            if "vendor" in path.parts or "tests" in path.parts or "build" in path.parts:
                continue
            for node in ast.walk(ast.parse(path.read_text())):
                if (
                    isinstance(node, ast.Call)
                    and isinstance(node.func, ast.Name)
                    and node.func.id == "record_open_file"
                ):
                    callers.append(f"{path.relative_to(app).as_posix()}:{node.lineno}")

        self.assertEqual(
            len(callers),
            1,
            f"the history has one writer, in the shared load tail; found {callers}",
        )
        self.assertTrue(callers[0].startswith("page_modules/data_loading.py"), callers)


class TestWhereTheListIsDrawn(RedirectedHistory):
    """Inside the chooser, so the sidebar column gains nothing unfolded when a file is open."""

    def _render(self, mock_st, **session):
        from components.sidebar import render_load_status

        state = FakeSessionState(**session)
        mock_st.session_state = state
        slot = Mock()
        slot.button.return_value = False
        slot.expander.return_value.button.return_value = False
        render_load_status(slot)
        return state, slot

    def _open_file_session(self):
        import pandas as pd

        return dict(
            current_page="parameter_config",
            maf_source_name="sample_42.maf",
            maf_data=pd.DataFrame({"Hugo_Symbol": ["TP53"] * 10}),
            filtered_data=None,
            filter_params={"sample_type": "somatic"},
        )

    @staticmethod
    def _said(box) -> str:
        return " ".join(
            str(call) for call in box.caption.call_args_list + box.markdown.call_args_list
        )

    @patch("components.sidebar.st")
    def test_an_open_file_puts_the_list_behind_the_chooser_expander(self, mock_st):
        """The column already holds four things. A fifth unfolded one is the question the map
        left in the fog, and this ticket does not get to answer it by taking the space."""
        from config.file_history import record_open_file

        record_open_file(SAMPLE_NAME)
        _, slot = self._render(mock_st, **self._open_file_session())

        self.assertIn(SAMPLE_NAME, self._said(slot.expander.return_value))
        self.assertNotIn(SAMPLE_NAME, self._said(slot))

    @patch("components.sidebar.st")
    def test_with_no_file_open_the_list_is_on_screen(self, mock_st):
        """The one state where it is unfolded is the state where it is most worth reading, and
        where the column is at its shortest."""
        from config.file_history import record_open_file

        record_open_file(SAMPLE_NAME)
        _, slot = self._render(mock_st, current_page="home", maf_data=None)

        self.assertIn(SAMPLE_NAME, self._said(slot))
        self.assertFalse(slot.expander.called)

    @patch("components.sidebar.st")
    def test_an_empty_history_draws_nothing(self, mock_st):
        """A first-time user has no history to explain, and a paragraph about what the app
        would keep, in a column whose job is to say what is open, introduces a feature to
        someone who has not used it."""
        _, slot = self._render(mock_st, current_page="home", maf_data=None)

        self.assertNotIn("Recently opened", self._said(slot))
        self.assertFalse(slot.button.called)

    @patch("components.sidebar.st")
    def test_the_time_is_drawn_beside_the_name(self, mock_st):
        from config.file_history import record_open_file

        record_open_file(SAMPLE_NAME, now=datetime.now())
        _, slot = self._render(mock_st, current_page="home", maf_data=None)

        self.assertIn("today", self._said(slot))

    @patch("components.sidebar.st")
    def test_nothing_in_the_list_is_clickable(self, mock_st):
        """The map ruled out re-opening from the history on **both** routes, so the list has
        exactly one control and it is the clear button. An entry that became a button would be
        the ruled-out feature arriving by the back door."""
        from config.file_history import record_open_file

        record_open_file(SAMPLE_NAME)
        record_open_file(OTHER_SAMPLE_NAME)
        _, slot = self._render(mock_st, current_page="home", maf_data=None)

        labels = [str(call.args[0]) if call.args else "" for call in slot.button.call_args_list]
        for label in labels:
            self.assertNotIn(".maf", label, f"an entry became a control: {label}")
        self.assertEqual(len(labels), 1, labels)
        self.assertFalse(slot.link_button.called)
        self.assertFalse(slot.download_button.called)


class TestWhatTheAppSaysAboutTheList(RedirectedHistory):
    """The sentence the dev made a condition of keeping the names at all."""

    def _sentence(self, mock_st) -> str:
        from components.sidebar import render_load_status
        from config.file_history import record_open_file

        record_open_file(SAMPLE_NAME)
        mock_st.session_state = FakeSessionState(current_page="home", maf_data=None)
        slot = Mock()
        slot.button.return_value = False
        render_load_status(slot)
        return " ".join(str(call) for call in slot.caption.call_args_list)

    @patch("components.sidebar.st")
    def test_the_cap_is_stated_and_derived(self, mock_st):
        """Derived rather than transcribed: a cap raised in the store and left unsaid here is
        the app promising to keep less than it keeps."""
        from config.file_history import HISTORY_LIMIT

        self.assertIn(str(HISTORY_LIMIT), self._sentence(mock_st))

    @patch("components.sidebar.st")
    def test_it_says_a_name_here_re_opens_nothing(self, mock_st):
        """Otherwise the list reads as a set of shortcuts that do not work."""
        self.assertIn("re-open", self._sentence(mock_st))

    @patch("components.sidebar.st")
    def test_it_says_a_name_outlives_the_file(self, mock_st):
        """The one thing about this list that looking at it cannot tell you, and the whole
        reason the map required a sentence: delete the MAF and the name is still here."""
        self.assertIn("gone", self._sentence(mock_st))

    @patch("components.sidebar.st")
    def test_it_is_written_for_a_clinician(self, mock_st):
        """The map's Style section: no implementation vocabulary reaches the interface. The
        exception it grants is MAF *column* names, which this surface has none of."""
        said = self._sentence(mock_st).lower()

        for word in ("json", "schema", "cache", "opened_files", ".mafigate", "stamp"):
            self.assertNotIn(word, said)


class TestTheClearControl(RedirectedHistory):
    """It clears the names, says so, and leaves the user's settings alone."""

    def _clear_control(self, mock_st, state):
        from components.sidebar import render_load_status

        mock_st.session_state = state
        slot = Mock()
        slot.button.return_value = False
        render_load_status(slot)
        return slot

    @patch("components.sidebar.st")
    def test_pressing_it_forgets_the_names(self, mock_st):
        from config.file_history import read_history, record_open_file

        record_open_file(SAMPLE_NAME)
        state = FakeSessionState(current_page="home", maf_data=None)
        slot = self._clear_control(mock_st, state)

        slot.button.call_args.kwargs["on_click"]()

        self.assertEqual(read_history(), [])

    @patch("components.sidebar.st")
    def test_it_does_not_touch_the_saved_filter_settings(self, mock_st):
        """Two different things a user might want gone. The parameter page has its own control
        for the settings, and clearing file names must not reconfigure anybody."""
        from page_modules import param_store
        from config.file_history import record_open_file

        record_open_file(SAMPLE_NAME)
        state = FakeSessionState(current_page="home", maf_data=None)
        slot = self._clear_control(mock_st, state)

        with patch.object(param_store, "clear_parameters_cache") as clearing_the_settings:
            slot.button.call_args.kwargs["on_click"]()

        self.assertFalse(clearing_the_settings.called)
        self.assertIn("filter settings", str(slot.button.call_args.kwargs["help"]))

    @patch("components.sidebar.st")
    def test_the_confirmation_is_drawn_on_the_render_that_follows(self, mock_st):
        """A callback, so the render below it reads an already-empty list — and the sentence
        cannot be drawn from inside a callback, where nothing is rendering yet.

        This is the shape issue #140 arrived at for the parameter page's confirmations, and
        ``test_discarded_frames.py`` guards the failure it avoids: a message drawn in a frame
        that is about to be thrown away reaches nobody.
        """
        from components.sidebar import HISTORY_CLEARED
        from config.file_history import record_open_file

        record_open_file(SAMPLE_NAME)
        state = FakeSessionState(current_page="home", maf_data=None)
        slot = self._clear_control(mock_st, state)
        slot.button.call_args.kwargs["on_click"]()
        self.assertTrue(state[HISTORY_CLEARED])

        after = self._clear_control(mock_st, state)
        said = " ".join(str(call) for call in after.caption.call_args_list)

        self.assertIn("has been cleared", said)
        self.assertNotIn(SAMPLE_NAME, said)
        self.assertNotIn(HISTORY_CLEARED, state, "the sentence would be drawn again")


class TestTheSuiteCannotWriteToYourHome(unittest.TestCase):
    """``tests/conftest.py``'s redirect, made load-bearing.

    The classes above redirect the store themselves, so they would pass with the fixture
    deleted. Every *other* module in the suite would not: several drive
    ``_open_the_file_just_read`` for real, and before the fixture existed one full run left
    seven fixture names in the developer's own ``~/.mafigate/opened_files.json``. Nothing
    failed when that happened, which is why this is asserted rather than trusted.
    """

    def test_the_store_is_redirected_for_every_test(self):
        """Both constants, not just the file.

        ``record_open_file`` calls ``mkdir`` on the *directory* before it opens the file, so a
        fixture that redirected the file alone would still reach into a real home — creating
        ``~/.mafigate`` on a machine that had none. Checking only the file left exactly that
        half-redirect passing.
        """
        from config import file_history

        for name in ("HISTORY_DIR", "HISTORY_FILE"):
            self.assertFalse(
                Path(getattr(file_history, name)).is_relative_to(Path.home()),
                f"file_history.{name} points into a real home directory during a test run; "
                "see tests/conftest.py",
            )


if __name__ == "__main__":
    unittest.main()
