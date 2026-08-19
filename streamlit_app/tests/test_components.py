"""
Test UI components and user interface elements.
"""

import unittest
import re
import sys
import os
from pathlib import Path
from unittest.mock import MagicMock, Mock, patch

import pandas as pd

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Defined in `tests/fakes.py`, and re-exported here because it was declared in this module
# first: every reference below still reads `FakeSessionState`, and a second reader appeared
# (`tests/test_app_identity.py`) that should not have to import a collected test module to
# get at it.
from tests.fakes import FakeSessionState  # noqa: E402,F401


class TestUIComponents(unittest.TestCase):
    """Test UI component functions."""

    def setUp(self):
        """Set up test fixtures."""
        # Mock streamlit to avoid import issues in testing
        self.st_mock = Mock()

    @patch("components.sidebar.st")
    def test_sidebar_navigation_creation(self, mock_st):
        """Test that sidebar navigation maps the selected radio label to a page key."""
        try:
            from components.sidebar import create_sidebar_navigation
        except ImportError:
            # If components.sidebar doesn't exist or has streamlit dependencies
            self.skipTest("UI components require Streamlit environment")

        # Default current page, and the label the user selected in the radio
        mock_st.session_state = FakeSessionState(current_page="home", maf_data=None)
        mock_st.sidebar.button.return_value = False
        mock_st.sidebar.radio.return_value = "📊 Load & Analyze Data"

        result = create_sidebar_navigation()

        # The function returns the page key corresponding to the selected label
        self.assertEqual(result, "data_loading")

    @patch("components.sidebar.st")
    def test_the_navigation_draws_no_button_of_its_own(self, mock_st):
        """The nav radio is the whole of the navigation below the status block (issue #161).

        A `❓ Need Help? Click here` button was drawn here, under a rule of its own, and it
        went to `help` — the radio's own fourth entry, two elements above it. The one button
        left in this column belongs to `render_load_status`, which is drawn from a different
        call, so *this* function drawing a button at all is the defect returning.

        `tests/test_sidebar_doors.py` holds the general rule by reading the source. This is
        the same claim made by running the code, which is what catches a button drawn from
        somewhere this file's AST sweep is not looking.
        """
        try:
            from components.sidebar import create_sidebar_navigation
        except ImportError:
            self.skipTest("UI components require Streamlit environment")

        mock_st.session_state = FakeSessionState(current_page="help", maf_data=None)
        mock_st.sidebar.radio.return_value = "❓ Help & Documentation"

        create_sidebar_navigation()

        mock_st.sidebar.button.assert_not_called()


class TestLoadStatus(unittest.TestCase):
    """The sidebar block that says what is open and offers the way back (issue #58).

    Leaving the results for the parameters keeps the loaded file; these assert that the
    interface now says so, in all three states — including the middle one, a file open
    with no results behind it, which must not read as no file at all.
    """

    def _render(self, mock_st, **session):
        """Render the block into a stand-in for the sidebar's reserved slot."""
        from components.sidebar import render_load_status

        mock_st.session_state = FakeSessionState(**session)
        slot = Mock()
        slot.button.return_value = False
        render_load_status(slot)

        written = [str(c.args[0]) for c in slot.markdown.call_args_list]
        written += [str(c.args[0]) for c in slot.caption.call_args_list]
        labels = [str(c.args[0]) for c in slot.button.call_args_list]
        return written, labels

    @patch("components.sidebar.st")
    def test_no_file_says_so_and_offers_no_route_of_its_own(self, mock_st):
        """The empty state says so, and offers the chooser rather than a way to reach one.

        The route button this state used to carry existed to take the user to the page that
        held the uploader. Since issue #64 the uploader is the next thing in this column, so
        the button would be the longer way round to what is already on screen. Where the
        chooser is drawn is ``test_file_chooser.py``'s claim; this asserts only that nothing
        else is offered beside it.
        """
        written, labels = self._render(
            mock_st, current_page="parameter_config", maf_data=None
        )

        self.assertTrue(any("No file open" in line for line in written))
        self.assertEqual(labels, [])

    @patch("components.sidebar.st")
    def test_loaded_and_filtered_names_the_file_and_both_counts(self, mock_st):
        written, labels = self._render(
            mock_st,
            current_page="parameter_config",
            maf_source_name="sample_42.maf",
            maf_data=pd.DataFrame({"Hugo_Symbol": ["TP53"] * 2718}),
            filtered_data=pd.DataFrame({"Hugo_Symbol": ["TP53"] * 412}),
            filter_params={"sample_type": "somatic"},
        )
        blob = " ".join(written)

        self.assertIn("sample_42.maf", blob)
        self.assertIn("2,718", blob)
        self.assertIn("412", blob)
        self.assertIn("Somatic", blob)
        # The way back is offered, and it is named in terms of the results
        self.assertEqual(len(labels), 1)
        self.assertIn("results", labels[0].lower())

    @patch("components.sidebar.st")
    def test_loaded_but_unfiltered_does_not_read_as_no_file(self, mock_st):
        """The middle state: `maf_data` set, `filtered_data` still None."""
        written, labels = self._render(
            mock_st,
            current_page="parameter_config",
            maf_source_name="sample_42.maf",
            maf_data=pd.DataFrame({"Hugo_Symbol": ["TP53"] * 2718}),
            filtered_data=None,
            filter_params={"sample_type": "germline"},
        )
        blob = " ".join(written)

        self.assertNotIn("No file open", blob)
        self.assertIn("sample_42.maf", blob)
        self.assertIn("2,718", blob)
        self.assertIn("Germline", blob)
        # Still a route back, even with nothing filtered yet
        self.assertEqual(len(labels), 1)

    def _press_the_way_back(self, mock_st, **session):
        """Render the block, press its button, and report what the press asked for.

        ``page_modules.data_loading``'s own ``st`` is patched alongside this module's,
        because the press reaches into that module to request its Results section and
        would otherwise write to the real session state.
        """
        from components.sidebar import render_load_status

        state = FakeSessionState(**session)
        mock_st.session_state = state
        slot = Mock()
        slot.button.return_value = True

        with patch("page_modules.data_loading.st") as mock_page_st:
            mock_page_st.session_state = state
            render_load_status(slot)

        return state

    @patch("components.sidebar.st")
    def test_the_way_back_to_results_asks_for_the_results_section(self, mock_st):
        """Naming the page is not enough since issue #59.

        The data page keeps a fixed section order and remembers the section you left it
        on, so a button labelled *back to your results* pressed after a visit to the
        filter-run section would land there instead (it was called Load Data until issue
        #65). Two features each right alone; the label is what lies.
        """
        state = self._press_the_way_back(
            mock_st,
            current_page="parameter_config",
            maf_source_name="sample_42.maf",
            maf_data=pd.DataFrame({"Hugo_Symbol": ["TP53"] * 2718}),
            filtered_data=pd.DataFrame({"Hugo_Symbol": ["TP53"] * 412}),
            filter_params={"sample_type": "somatic"},
        )

        self.assertEqual(state.get("current_page"), "data_loading")
        self.assertTrue(
            state.get("jump_to_results"),
            "the way back to your results did not ask for the Results section",
        )

    @patch("components.sidebar.st")
    def test_the_way_back_to_an_unfiltered_file_asks_for_no_section(self, mock_st):
        """The middle state has no report to open, so it must not request one."""
        state = self._press_the_way_back(
            mock_st,
            current_page="parameter_config",
            maf_source_name="sample_42.maf",
            maf_data=pd.DataFrame({"Hugo_Symbol": ["TP53"] * 2718}),
            filtered_data=None,
            filter_params={"sample_type": "germline"},
        )

        self.assertEqual(state.get("current_page"), "data_loading")
        self.assertFalse(state.get("jump_to_results"))

    @patch("components.sidebar.st")
    def test_no_way_back_button_when_already_on_the_data_page(self, mock_st):
        """A button that navigates where you already are is noise, not a route."""
        written, labels = self._render(
            mock_st,
            current_page="data_loading",
            maf_source_name="sample_42.maf",
            maf_data=pd.DataFrame({"Hugo_Symbol": ["TP53"] * 10}),
            filtered_data=pd.DataFrame({"Hugo_Symbol": ["TP53"] * 3}),
            filter_params={"sample_type": "somatic"},
        )

        self.assertEqual(labels, [])
        self.assertIn("sample_42.maf", " ".join(written))

    @patch("components.sidebar.st")
    def test_file_without_a_recorded_name_still_reads_as_open(self, mock_st):
        written, _ = self._render(
            mock_st,
            current_page="home",
            maf_data=pd.DataFrame({"Hugo_Symbol": ["TP53"] * 5}),
            filtered_data=None,
        )
        blob = " ".join(written)

        self.assertNotIn("No file open", blob)
        self.assertIn("Your MAF file", blob)

    @patch("components.sidebar.st")
    def test_unfiltered_copy_does_not_claim_the_filters_never_ran(self, mock_st):
        """Both load paths filter at once, and a refusal leaves no results behind.

        So the usual way into this state is the filters running and *failing* — copy that
        says they have not run tells the clinician the opposite of what happened.
        """
        written, _ = self._render(
            mock_st,
            current_page="parameter_config",
            maf_source_name="sample_42.maf",
            maf_data=pd.DataFrame({"Hugo_Symbol": ["TP53"] * 2718}),
            filtered_data=None,
            filter_error="🛑 This MAF cannot be filtered: …",
            filter_params={"sample_type": "somatic"},
        )
        blob = " ".join(written)

        self.assertNotIn("have not run", blob)
        self.assertNotIn("No file open", blob)
        self.assertIn("No results for this file yet", blob)

    @patch("components.sidebar.st")
    def test_results_are_not_called_current_once_the_settings_have_moved_on(self, mock_st):
        """The data page warns that the settings changed; this block must not disagree."""
        from utils.main_utils import params_hash

        written, _ = self._render(
            mock_st,
            current_page="parameter_config",
            maf_source_name="sample_42.maf",
            maf_data=pd.DataFrame({"Hugo_Symbol": ["TP53"] * 2718}),
            filtered_data=pd.DataFrame({"Hugo_Symbol": ["TP53"] * 412}),
            filter_params={"sample_type": "somatic", "min_depth": 30},
            data_params_hash=params_hash({"sample_type": "somatic", "min_depth": 20}),
        )
        blob = " ".join(written)

        self.assertIn("412", blob)
        self.assertNotIn("current filters", blob)
        self.assertIn("last filter run", blob)
        self.assertIn("Settings have changed since", blob)

    @patch("components.sidebar.st")
    def test_results_are_current_while_the_settings_still_match(self, mock_st):
        from utils.main_utils import params_hash

        params = {"sample_type": "somatic", "min_depth": 20}
        written, _ = self._render(
            mock_st,
            current_page="parameter_config",
            maf_source_name="sample_42.maf",
            maf_data=pd.DataFrame({"Hugo_Symbol": ["TP53"] * 2718}),
            filtered_data=pd.DataFrame({"Hugo_Symbol": ["TP53"] * 412}),
            filter_params=params,
            data_params_hash=params_hash(params),
        )
        blob = " ".join(written)

        self.assertIn("passed your current filters", blob)
        self.assertNotIn("Settings have changed", blob)

    @patch("components.sidebar.st")
    def test_status_slot_is_filled_after_the_page_not_in_call_order(self, mock_st):
        """The sidebar is built before the page that opens the file.

        So the block must be drawn into a reserved slot at the end of the render, or it
        announces "No file open" on the very render that opens the file.
        """
        from components.sidebar import (
            LOAD_STATUS_SLOT,
            create_sidebar_navigation,
            render_into_status_slot,
        )

        mock_st.session_state = FakeSessionState(current_page="home", maf_data=None)
        mock_st.sidebar.button.return_value = False
        mock_st.sidebar.radio.return_value = "🏠 Home"

        create_sidebar_navigation()

        # Space claimed, nothing drawn into it yet
        slot = mock_st.session_state[LOAD_STATUS_SLOT]
        self.assertIs(slot, mock_st.sidebar.empty.return_value)
        self.assertFalse(slot.container.called)

        # The page has now run and opened a file; the slot reflects it
        mock_st.session_state["maf_data"] = pd.DataFrame({"Hugo_Symbol": ["TP53"] * 7})
        mock_st.session_state["maf_source_name"] = "opened_late.maf"
        render_into_status_slot()

        body = slot.container.return_value
        written = [str(c.args[0]) for c in body.markdown.call_args_list]
        written += [str(c.args[0]) for c in body.caption.call_args_list]
        blob = " ".join(written)
        self.assertIn("opened_late.maf", blob)
        self.assertNotIn("No file open", blob)

    @patch("components.sidebar.st")
    def test_status_slot_render_is_a_no_op_without_a_sidebar(self, mock_st):
        from components.sidebar import render_into_status_slot

        mock_st.session_state = FakeSessionState()
        render_into_status_slot()  # must not raise

        self.assertFalse(mock_st.sidebar.markdown.called)

    @patch("components.sidebar.st")
    def test_copy_carries_no_implementation_vocabulary(self, mock_st):
        """The map's Style section: this block is read by clinicians."""
        written, labels = self._render(
            mock_st,
            current_page="home",
            maf_source_name="sample_42.maf",
            maf_data=pd.DataFrame({"Hugo_Symbol": ["TP53"] * 8}),
            filtered_data=pd.DataFrame({"Hugo_Symbol": ["TP53"] * 4}),
            filter_params={"sample_type": "somatic"},
        )
        blob = " ".join(written + labels).lower()

        for jargon in ("filtered_data", "maf_data", "session", "row", "dataframe"):
            self.assertNotIn(jargon, blob)

    def test_parameter_validation_functions(self):
        """Test parameter validation utilities."""
        # Test numeric parameter validation
        self.assertTrue(self._validate_numeric_parameter(10, 0, 100))
        self.assertFalse(self._validate_numeric_parameter(-5, 0, 100))
        self.assertFalse(self._validate_numeric_parameter(150, 0, 100))

    def test_file_validation_functions(self):
        """Test file validation utilities."""
        # Test file extension validation
        self.assertTrue(self._validate_file_extension("test.maf", [".maf", ".tsv"]))
        self.assertTrue(self._validate_file_extension("test.tsv", [".maf", ".tsv"]))
        self.assertFalse(self._validate_file_extension("test.txt", [".maf", ".tsv"]))

    def _validate_numeric_parameter(self, value, min_val, max_val):
        """Helper function to validate numeric parameters."""
        try:
            num_value = float(value)
            return min_val <= num_value <= max_val
        except (ValueError, TypeError):
            return False

    def _validate_file_extension(self, filename, allowed_extensions):
        """Helper function to validate file extensions."""
        if not filename:
            return False

        filename_lower = filename.lower()
        return any(filename_lower.endswith(ext.lower()) for ext in allowed_extensions)


class TestSessionState(unittest.TestCase):
    """Test session state management."""

    def test_session_state_initialization(self):
        """Test that session state can be properly initialized."""
        # Mock session state
        mock_session_state = {}

        # Test initialization logic
        default_values = {
            "current_page": "home",
            "maf_data": None,
            "filter_params": {},
            "last_filter_timestamp": None,
        }

        for key, default_value in default_values.items():
            if key not in mock_session_state:
                mock_session_state[key] = default_value

        # Verify initialization
        for key, expected_value in default_values.items():
            self.assertIn(key, mock_session_state)
            self.assertEqual(mock_session_state[key], expected_value)

    def test_parameter_persistence(self):
        """Test parameter persistence logic."""
        # Mock parameter storage and retrieval
        params = {
            "sample_type": "somatic",
            "min_depth": 20,
            "vaf_threshold": 0.1,
            "filter_variant_classification": ["Missense_Mutation"],
        }

        # Test parameter serialization/deserialization
        import json

        # Should be able to serialize to JSON
        serialized = json.dumps(params)
        self.assertIsInstance(serialized, str)

        # Should be able to deserialize back
        deserialized = json.loads(serialized)
        self.assertEqual(params, deserialized)


class TestErrorHandling(unittest.TestCase):
    """Test error handling in UI components."""

    def test_file_upload_error_handling(self):
        """Test file upload error scenarios."""
        # Test empty file handling
        self.assertFalse(self._is_valid_upload(None))
        self.assertFalse(self._is_valid_upload(""))

        # Test invalid file types
        self.assertFalse(self._is_valid_upload_type("test.pdf"))
        self.assertTrue(self._is_valid_upload_type("test.maf"))

    def test_parameter_validation_errors(self):
        """Test parameter validation error cases."""
        # Test invalid depth values
        invalid_depths = [-1, "abc", None, ""]
        for depth in invalid_depths:
            self.assertFalse(self._is_valid_depth(depth))

        # Test valid depth values
        valid_depths = [0, 10, 50, "20", "100"]
        for depth in valid_depths:
            self.assertTrue(self._is_valid_depth(depth))

    def test_data_processing_error_handling(self):
        """Test data processing error scenarios."""
        # Test empty data handling
        import pandas as pd

        empty_df = pd.DataFrame()
        self.assertTrue(empty_df.empty)

        # Test malformed data handling
        malformed_data = {"invalid": "structure"}
        self.assertFalse(self._is_valid_maf_structure(malformed_data))

    def _is_valid_upload(self, file):
        """Helper function to validate file uploads."""
        return file is not None and file != ""

    def _is_valid_upload_type(self, filename):
        """Helper function to validate upload file types."""
        if not filename:
            return False
        allowed_extensions = [".maf", ".tsv", ".txt"]
        return any(filename.lower().endswith(ext) for ext in allowed_extensions)

    def _is_valid_depth(self, depth):
        """Helper function to validate depth parameters."""
        if depth is None or depth == "":
            return False
        try:
            num_depth = float(depth)
            return num_depth >= 0
        except (ValueError, TypeError):
            return False

    def _is_valid_maf_structure(self, data):
        """Helper function to validate MAF data structure."""
        import pandas as pd

        return isinstance(data, pd.DataFrame) and not data.empty


class TestAccessibility(unittest.TestCase):
    """Test accessibility features of UI components."""

    def test_form_labels(self):
        """Test that form elements have proper labels."""
        # Verify that form elements would have descriptive labels
        form_elements = {
            "min_depth": "Minimum Read Depth",
            "vaf_threshold": "Variant Allele Frequency Threshold",
            "sample_type": "Sample Type Selection",
        }

        for element_id, expected_label in form_elements.items():
            self.assertIsInstance(expected_label, str)
            self.assertTrue(len(expected_label) > 0)

    def test_help_text_availability(self):
        """Test that help text is available for complex parameters."""
        help_texts = {
            "max_freq_population": "Maximum population frequency in public databases",
            "filter_variant_classification": "Types of variants to include in analysis",
            "clinical_databases": "Clinical significance databases for variant annotation",
        }

        for parameter, help_text in help_texts.items():
            self.assertIsInstance(help_text, str)
            self.assertTrue(len(help_text) > 10)  # Should be descriptive


class TestChartLayout(unittest.TestCase):
    """The dashboard charts leave room for the text they draw.

    Whether a label is legible is a rendering question, and these tests do not
    pretend to answer it — that was checked by eye against a real MAF. What they
    hold is the two layout settings that answer it, which are invisible in the
    app until someone looks at a chart, and so would otherwise only be missed by
    another walkthrough:

    * ``automargin``, without which a fixed margin clipped the axis titles and
      cut gene symbols in half;
    * ``autotypenumbers="strict"``, without which chromosome labels are read as
      numbers and the axis becomes a number line.

    They also hold the shape of the fix: a chart asks ``_compact_layout`` for the
    height it wants, rather than re-setting the layout the helper just wrote.

    **Where these run, said out loud.** They need ``plotly``, and so they still
    skip in the ``vendor drift`` CI job, which is deliberately narrow (pytest and
    pandas, per its own header comment). Saying so is the point: issue #24 was not
    a skip, it was a skip nobody had written down. Whoever widens CI to the full
    suite inherits six live tests — five when this note was first written, and the
    clinical-order test has joined them since.

    They also used to skip on the developer's own machine, which is the half of
    this note that has changed. ``plotly`` was declared in ``requirements.txt``
    and installed in no interpreter on it, so these six ran nowhere at all and
    every Summary chart fell back to ``st.bar_chart``. Issue #162 installed it
    into the interpreter that runs the app and made ``setup.sh`` verify the
    install rather than assume it, and all six passed first time under plotly
    6.9.0 — a major version above the ``>=5.15.0`` floor they were written
    against. So this class is now a live guard wherever the app itself can run,
    and still not a gate anything crosses on the way to merge.
    """

    def setUp(self):
        try:
            import plotly  # noqa: F401
        except ImportError:
            self.skipTest("plotly not installed")

        # Chromosome labels that look like numbers, plus one that does not.
        self.frame = pd.DataFrame(
            {
                "Chromosome": ["chr1"] * 3 + ["chr2"] * 2 + ["chr17", "chrX"],
                "Hugo_Symbol": ["TP53", "TP53", "BRCA2", "EGFR", "MTOR", "PIK3CA", "MUTYH"],
                "Variant_Classification": ["Missense_Mutation"] * 5 + ["Silent"] * 2,
                "Variant_Type": ["SNP"] * 6 + ["DEL"],
                "tumor_f": [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7],
                "Clinical_Summary": [clinical_summary._SUMMARY_LABELS["Pathogenic"]] * 4
                + [clinical_summary._SUMMARY_LABELS["Benign"]] * 3,
            }
        )

    def _figure_from(self, plot_name, *args):
        """Call a chart function and return the figure it hands to Streamlit."""
        from components import charts

        captured = []
        with patch.object(charts, "st") as mock_st:
            mock_st.plotly_chart.side_effect = lambda fig, *a, **k: captured.append(fig)
            getattr(charts, plot_name)(*args)

        self.assertEqual(len(captured), 1, f"{plot_name} drew no chart")
        return captured[0]

    #: The dashboard, as one table: the chart's function, the height it is
    #: entitled to ask ``_compact_layout`` for, whether it compares passed
    #: against failed (and so takes two frames), and whether it has axes at all
    #: — the donuts do not. One row per chart, so a chart cannot be in one of
    #: these tests and quietly absent from another.
    CHARTS = [
        # function name, height, splits passed/failed, has axes
        ("_plot_vaf_distribution", 350, True, True),
        ("_plot_variant_classification", 350, False, False),
        ("_plot_chromosome_distribution", 400, False, True),
        ("_plot_top_genes", 400, False, True),
        ("_plot_clinical_significance", 380, True, True),
        ("_plot_mutation_types", 350, False, False),
    ]

    def _figure_for(self, chart):
        """Build the figure for one CHARTS row, feeding it the frames it takes."""
        plot_name, _, splits, _ = chart
        args = (self.frame.iloc[:4], self.frame.iloc[4:]) if splits else (self.frame,)
        return self._figure_from(plot_name, *args)

    def test_every_chart_keeps_the_height_it_asked_for(self):
        """Each chart ends at the height it is entitled to.

        This one pins the heights rather than catching the old shape: a chart
        that re-sets its height after ``_compact_layout`` lands on the same
        number, so the difference is invisible from the figure. What it holds is
        that changing a height is a deliberate act.
        """
        for chart in self.CHARTS:
            with self.subTest(chart=chart[0]):
                self.assertEqual(self._figure_for(chart).layout.height, chart[1])

    def test_every_chart_keeps_numeric_looking_labels_as_categories(self):
        """Streamlit's template would otherwise turn '1'..'22' into a number line."""
        for chart in self.CHARTS:
            with self.subTest(chart=chart[0]):
                fig = self._figure_for(chart)
                self.assertEqual(fig.layout.autotypenumbers, "strict")

    def test_the_axis_charts_let_their_axes_claim_the_room_they_need(self):
        """automargin on both axes of every chart that has axes."""
        for chart in [c for c in self.CHARTS if c[3]]:
            with self.subTest(chart=chart[0]):
                fig = self._figure_for(chart)
                self.assertTrue(fig.layout.xaxis.automargin, "x axis cannot grow")
                self.assertTrue(fig.layout.yaxis.automargin, "y axis cannot grow")

    def test_the_chromosome_chart_names_every_chromosome_it_plots(self):
        """Along the bottom, in chromosome order, including the non-numeric ones."""
        fig = self._figure_from("_plot_chromosome_distribution", self.frame)
        self.assertEqual(list(fig.data[0].x), ["1", "2", "17", "X"])
        self.assertEqual(fig.layout.xaxis.title.text, "Chromosome")
        self.assertEqual(fig.layout.yaxis.title.text, "Variant Count")

    def test_the_clinical_chart_says_nothing_to_the_reader_in_its_own_words(self):
        """Neither the axis nor the tooltip carries this function's column name.

        The tooltip is the half that is easy to miss: Plotly builds its field
        names from the frame's columns, so a column named for the code leaks into
        rendered text even when the axis title is blank. Blanking the title
        *through* ``labels`` is worse still — it empties the tooltip's field name
        too and leaves a dangling "=" on that line.
        """
        fig = self._figure_from(
            "_plot_clinical_significance", self.frame.iloc[:4], self.frame.iloc[4:]
        )
        self.assertIn(fig.layout.xaxis.title.text, (None, ""))

        for trace in fig.data:
            self.assertNotIn("Category", trace.hovertemplate)
            self.assertNotIn("<br>=", trace.hovertemplate, "orphaned = in the tooltip")
            self.assertIn("Clinical Significance=", trace.hovertemplate)

    def test_the_clinical_chart_runs_in_clinical_order_not_frequency_order(self):
        """Severity is the order a clinician scans for; counts are not.

        Left to ``value_counts()`` the axis came out in frequency order, which
        put Pathogenic between Likely Benign and Benign. ``CLINICAL_COLORS``
        already holds the right order, so that is the one to follow.
        """
        from components.charts import CLINICAL_COLORS

        every_category = pd.DataFrame(
            # Descending counts, deliberately the reverse of clinical order, so
            # frequency order and clinical order cannot agree by accident.
            {
                "Clinical_Summary": [
                    name
                    for i, name in enumerate(reversed(list(CLINICAL_COLORS)))
                    for _ in range(i + 1)
                ]
            }
        )
        fig = self._figure_from(
            "_plot_clinical_significance", every_category, every_category
        )

        drawn = [c for c in fig.layout.xaxis.categoryarray or () if c]
        self.assertEqual(drawn, list(CLINICAL_COLORS))


# ---------------------------------------------------------------------------
# Where the user's own writing sits in the table (issue #60)
# ---------------------------------------------------------------------------
#
# `Notes` is free text the user typed into the variant dialog. Unpinned it sat at
# APP_EXTRA_COLUMNS' third slot, which on the somatic arm is past roughly forty pipeline
# columns — the one column they wrote themselves, further right than every column they
# did not.
#
# The fix is pinning, not reordering, and the distinction is the whole point:
# `config/columns.py` requires the resolver's output to open with the pipeline's list as
# an exact prefix, so nothing can be moved up it. Pinning is a viewport decision applied
# after the resolver, so the table can lead with `Notes` while the column list, and
# therefore the export, is untouched. These tests pin both halves of that — that the
# lead is right, and that the list underneath it did not move.

import pytest

pytest.importorskip("streamlit", reason="UI components require Streamlit")

from components import charts, clinical_summary, variant_table  # noqa: E402
from config.columns import APP_EXTRA_COLUMNS, resolve_visible_columns  # noqa: E402
from config.vocabularies import CIVIC_OPTIONS  # noqa: E402


@pytest.fixture
def annotations(monkeypatch):
    """Set `st.session_state.custom_annotations` to a real dict and return a setter.

    A `Mock` session state would answer `.get("custom_annotations", {})` with a `Mock`,
    which sorts into nonsense rather than failing — so the double here is a real dict.
    `FakeSessionState` is that dict and was already in `tests/fakes.py` for this exact
    reason; this fixture and the two below each declared their own copy of it until #127.
    """
    session = FakeSessionState()
    monkeypatch.setattr(variant_table.st, "session_state", session, raising=False)

    def _set(*names):
        session["custom_annotations"] = {name: {} for name in names}

    return _set


def test_notes_is_buried_in_the_column_list_it_is_pinned_out_of(annotations):
    """The premise, asserted — otherwise the pin is protecting against nothing.

    If a later change moves `Notes` forward in the resolver's output this fails, and
    whoever made that change gets to decide whether pinning is still wanted. That is the
    conversation worth forcing; a silently redundant pin is not.
    """
    annotations()
    columns = resolve_visible_columns("somatic")

    assert columns.index("Notes") > 30, (
        "Notes has moved forward in the resolver's output — re-read issue #60 before "
        "assuming the pin is still what makes it findable"
    )


def test_notes_leads_the_table_next_to_the_clinical_verdict(annotations):
    """What the user sees first: the app's verdict, then their own writing."""
    annotations()
    columns = resolve_visible_columns("somatic")

    assert variant_table._leading_columns(columns) == [
        "Clinical_Summary",
        "Pathogenicity_Overview",
        "Notes",
    ]


def test_a_users_own_annotation_columns_lead_too_and_follow_notes(annotations):
    """`custom_annotations` invents columns by the same mechanism, so it gets the same
    treatment — behind `Notes`, since `Notes` is the one every session has."""
    annotations("Review_Status", "Actionable")
    columns = resolve_visible_columns("somatic") + ["Actionable", "Review_Status"]

    assert variant_table._leading_columns(columns) == [
        "Clinical_Summary",
        "Pathogenicity_Overview",
        "Notes",
        "Actionable",
        "Review_Status",
    ]


def test_a_lead_column_the_frame_does_not_carry_is_not_conjured(annotations):
    """`_leading_columns` is filtered by presence: callers index a frame with it.

    An annotation column the user created *after* this frame was built, or `Notes`
    before `_add_derived_columns` has run, must drop out rather than raise a `KeyError`
    deep in AgGrid.
    """
    annotations("Actionable")

    assert variant_table._leading_columns(["Hugo_Symbol", "Clinical_Summary"]) == [
        "Clinical_Summary"
    ]


def test_the_fallback_table_shows_every_column_the_grid_would(annotations):
    """`_display_order` reorders for the AgGrid-less path. It must not also *select*.

    A permutation, checked as a multiset: dropping a column here would lose it from the
    fallback table only, which is exactly the kind of divergence nobody would notice
    until the machine without st_aggrid is the one in front of a user.
    """
    annotations("Actionable")
    columns = resolve_visible_columns("somatic") + ["Actionable"]

    ordered = variant_table._display_order(columns)

    assert sorted(ordered) == sorted(columns)
    assert ordered[:3] == ["Clinical_Summary", "Pathogenicity_Overview", "Notes"]


def test_pinning_did_not_move_notes_up_the_resolvers_output(annotations):
    """The contract `config/columns.py` documents, restated from this side of it.

    `test_column_resolver.py` owns the prefix rule itself. What this adds is that the
    display fix left it alone: the extras are still the tail, in their own order, with
    `Notes` third among them and nowhere near the front of the list that the export
    writes out.
    """
    annotations()
    columns = resolve_visible_columns("somatic")

    assert columns[-len(APP_EXTRA_COLUMNS) :] == APP_EXTRA_COLUMNS
    assert APP_EXTRA_COLUMNS[2] == "Notes"


# --- What a download carries (issue #67) -------------------------------------------
#
# Nothing stores a note: `variant_notes` and `custom_annotations` live in session state
# and no code writes either to disk. That was decided rather than merely observed — the
# real store is institute-wide and hosted, and building a per-laptop one first would be
# the wrong reach and a migration to drain later. The consequence is that a download is
# the *only* way a note leaves the session, which promotes these from tidiness to the
# app's one unaffordable silent failure.


@pytest.fixture
def noted_frame(monkeypatch):
    """A two-variant frame plus a note and an invented annotation column on the first.

    Session state is a real dict for the reason the `annotations` fixture gives: a `Mock`
    answers `.get(...)` with a `Mock`, which flows onward instead of failing. It also
    answers to attribute access, because the code under test reaches these stores both
    ways — `st.session_state.variant_notes` when building the columns,
    `st.session_state.get("custom_annotations", {})` when listing them — and a double
    that served only one would pass while the app raised. That is `FakeSessionState`'s whole
    job, so this uses it rather than redeclaring it.
    """
    session = FakeSessionState()
    monkeypatch.setattr(variant_table.st, "session_state", session, raising=False)

    frame = pd.DataFrame(
        {
            "Hugo_Symbol": ["TP53", "BRCA1"],
            "Chromosome": ["17", "17"],
            "Start_Position": [7577120, 41246481],
            "Reference_Allele": ["C", "G"],
            "Tumor_Seq_Allele2": ["T", "A"],
        }
    )
    key = variant_table._variant_key(frame.iloc[0])
    session["variant_notes"] = {key: "hotspot, discussed at board"}
    session["custom_annotations"] = {"Actionable": {key: "yes"}}
    return frame


def test_a_download_off_the_stored_frame_carries_what_the_user_typed(noted_frame):
    """`with_user_columns` is what the *Download Results* buttons were missing.

    They serialised `st.session_state.filtered_data` directly. `Notes` is built onto a
    copy at display time, so that frame has never carried it and those two buttons
    exported every pipeline column and none of the user's own — the export being, at the
    same time, the only place a note can survive.
    """
    exported = variant_table.with_user_columns(noted_frame)

    assert exported.loc[0, "Notes"] == "hotspot, discussed at board"
    assert exported.loc[0, "Actionable"] == "yes"
    assert exported.loc[1, "Notes"] == ""


def test_the_stored_frame_is_left_alone(noted_frame):
    """The derived columns go on a copy, and this pins that they still do.

    It is the property that made the bug invisible: every call site *looked* like it was
    annotating the frame everyone reads. A `with_user_columns` that mutated its argument
    would fix the downloads by accident and put user-authored columns into the frame the
    filters and the parity harness read, which is a different and worse change.
    """
    variant_table.with_user_columns(noted_frame)

    assert "Notes" not in noted_frame.columns
    assert "Actionable" not in noted_frame.columns


def test_the_shown_columns_download_carries_the_columns_shown(noted_frame):
    """The grid and the "shown columns" download read one list, so they cannot disagree.

    They used to compute it twice: the grid appended the user's invented columns after
    the resolver, the download called the resolver alone. `Notes` survived that because
    it is in `APP_EXTRA_COLUMNS`; a column the user named did not, so it was on screen
    and absent from the CSV offered directly beneath it. The resolver cannot know these
    names — it is streamlit-free and they live in session state — so the append is the
    app's, and one place is the only number of places it can live.

    This holds the *seam*. Since issue #92 the download no longer calls it: it takes the
    list ``create_data_table`` returned, so it follows the grid after the user adds a
    column rather than only agreeing with what the grid opened with. That wiring is
    ``test_download_contents.py``'s, which drives the multiselect and compares the CSV's
    header against the grid's own list — the same division of labour as issue #67, where
    both seams existed and the page did not use them.
    """
    exported = variant_table.with_user_columns(noted_frame)

    columns = variant_table.shown_columns(exported)

    assert "Notes" in columns
    assert "Actionable" in columns
    # The resolver's order is the pipeline's, and the user's own columns come after it —
    # pinning is what brings them to the left edge of the grid, not this list.
    assert columns.index("Actionable") > columns.index("Notes")


def test_an_annotation_column_the_frame_does_not_carry_is_not_asked_for(noted_frame):
    """A column invented after this frame was built must drop out, not raise.

    Callers index a frame with this list. The dialog can create a column at any time,
    including from a rerun that happens after the download's frame was assembled.
    """
    exported = variant_table.with_user_columns(noted_frame)

    # Created after the frame was assembled — the order is the whole test.
    variant_table.st.session_state["custom_annotations"]["Invented_Later"] = {}

    assert "Invented_Later" not in variant_table.shown_columns(exported)


# ---------------------------------------------------------------------------
# What happens when the grid cannot be drawn (issues #73, #76)
# ---------------------------------------------------------------------------
#
# The handler that answers this had never run. It referred four times to `cleaned_data`,
# a local of `create_data_table` several hundred lines away, so *every* AgGrid failure
# raised `NameError` from inside the handler for the original error — losing the table
# and the reason for it together. Nothing here exercised the path, and `app-load-check`
# renders the real grid successfully, so neither instrument could see it; it was found by
# `flake8 --select=F821` and fixed with the split that gave the grid its own module.
#
# So these tests make AgGrid throw, which is the one thing no other test does.


class _GridFailure(Exception):
    """Stands in for what AgGrid raises on data it cannot serialise."""


@pytest.fixture
def broken_grid(monkeypatch):
    """Wire `variant_table` so the AgGrid call raises, and capture what Streamlit draws.

    `GridOptionsBuilder` is doubled too, so this holds in a checkout where `st_aggrid`
    was never installed — the failure being tested is the app's, not the library's.
    """
    fake_st = MagicMock()
    monkeypatch.setattr(variant_table, "st", fake_st)
    monkeypatch.setattr(variant_table, "GridOptionsBuilder", MagicMock(), raising=False)

    def _explode(*args, **kwargs):
        raise _GridFailure("could not serialise column 'Otherinfo'")

    monkeypatch.setattr(variant_table, "AgGrid", _explode, raising=False)
    return fake_st


@pytest.fixture
def two_variants():
    return pd.DataFrame(
        {
            "Clinical_Summary": [
                clinical_summary._SUMMARY_LABELS["Pathogenic"],
                clinical_summary._SUMMARY_LABELS["Benign"],
            ],
            "Hugo_Symbol": ["TP53", "BRCA2"],
            "Chromosome": ["17", "13"],
            "Start_Position": [7577120, 32900000],
        }
    )


def test_a_grid_that_cannot_be_drawn_falls_back_to_a_table(broken_grid, two_variants):
    """The fallback renders, which for four references' worth of history it did not.

    This is the regression guard: repair the name and this passes, revert it and the
    `NameError` fires before anything is drawn, so `st.dataframe` is never reached.
    """
    variant_table._render_aggrid_with_detail(
        two_variants, two_variants, 500, "passed_variants", list(two_variants.columns)
    )

    assert broken_grid.dataframe.called, (
        "the grid failed and nothing was drawn in its place — the fallback is not "
        "falling back (see issue #73)"
    )


def test_the_fallback_shows_the_frame_the_grid_was_handed(broken_grid, two_variants):
    """`display_data`, not `full_data`.

    The user chose these columns; answering a failure to show forty of them by showing
    four hundred is not a fallback, it is a different screen. `full_data` here carries a
    column the display frame does not, so passing the wrong one is visible.
    """
    display = two_variants[["Clinical_Summary", "Hugo_Symbol"]]

    variant_table._render_aggrid_with_detail(
        display, two_variants, 500, "passed_variants", list(display.columns)
    )

    drawn = broken_grid.dataframe.call_args.args[0]
    assert list(drawn.columns) == ["Clinical_Summary", "Hugo_Symbol"]


def test_the_fallback_says_what_actually_broke(broken_grid, two_variants):
    """The original error survives the handler rather than being replaced by its own."""
    variant_table._render_aggrid_with_detail(
        two_variants, two_variants, 500, "passed_variants", list(two_variants.columns)
    )

    said = " ".join(str(call) for call in broken_grid.error.call_args_list)
    assert "could not serialise column 'Otherinfo'" in said


def test_nothing_after_the_failed_render_is_attempted(broken_grid, two_variants):
    """The handler returns, because there is no grid response to read.

    Falling through would reach `grid_response`, which the failed call never bound —
    the original defect in a second form. The caption below the grid explains *its*
    filter boxes, so drawing it under a table that has none would be false as well.
    """
    variant_table._render_aggrid_with_detail(
        two_variants, two_variants, 500, "passed_variants", list(two_variants.columns)
    )

    assert not broken_grid.caption.called


def test_a_failure_opening_a_variant_is_not_reported_as_a_broken_table(
    monkeypatch, two_variants
):
    """The guard ends where the render ends (issue #73's second question).

    This `try` used to reach as far as `_show_variant_dialog`, so a failure while opening
    one variant's details was reported as "Error displaying table" and replaced a grid
    that had drawn perfectly well with a fallback for it. Catching broadly is right for
    the AgGrid call — there is no narrower exception to name — and wrong for everything
    the user does afterwards.
    """
    fake_st = MagicMock()
    fake_st.button.return_value = True
    monkeypatch.setattr(variant_table, "st", fake_st)
    monkeypatch.setattr(variant_table, "GridOptionsBuilder", MagicMock(), raising=False)

    grid_response = MagicMock()
    grid_response.selected_rows = [{"Hugo_Symbol": "TP53", "Chromosome": "17"}]
    monkeypatch.setattr(
        variant_table, "AgGrid", lambda *a, **k: grid_response, raising=False
    )

    def _dialog_explodes(row):
        raise _GridFailure("the detail panel raised")

    monkeypatch.setattr(variant_table, "_show_variant_dialog", _dialog_explodes)

    with pytest.raises(_GridFailure):
        variant_table._render_aggrid_with_detail(
            two_variants, two_variants, 500, "passed_variants", list(two_variants.columns)
        )

    assert not fake_st.dataframe.called, (
        "a dialog failure drew the table's fallback — the grid rendered fine, and "
        "replacing it misreports which part of the screen broke"
    )


# --- The Pathogenicity Overview circles and the key beside them (issue #100) ------------
#
# Nothing covered `generate_pathogenicity_circles` or its key, which is how two glyphs came
# to be drawn beside a key that did not list them: `💊` for ClinVar `drug_response`, and `⚪`
# for anything the mapping could not classify — reached by ESCAT `V` on every somatic MAF
# carrying a `TP53` SNP in a head-and-neck sample, since the only `TP53` row in
# `resources/escat_tiering.csv` is that tier's wildcard. The user met an unexplained white
# circle immediately beside the `⬜` the same key defines, as "No data".
#
# The claims here are read off what *renders*, not off the constants alone, which is #79's
# line: the key was a hand-written HTML string, and a guard that reads a module cannot see
# what the module draws. `broken_grid` serves for this — the key is written before the
# `try`, so a grid that cannot draw still draws its key.


def _key_text(fake_st):
    """Everything the grid wrote above the table, as one string."""
    return " ".join(str(call) for call in fake_st.markdown.call_args_list)


@pytest.fixture
def escat_rows():
    """One row per level `ESCAT_DEFINITIONS` holds, and two the scale does not define."""
    from config.vocabularies import ESCAT_DEFINITIONS

    levels = list(ESCAT_DEFINITIONS) + ["IV", "X"]
    return pd.DataFrame({"ESCAT": levels, "Hugo_Symbol": ["TP53"] * len(levels)})


@pytest.fixture
def civic_rows():
    """Cells in the shape the annotation writes them, plus the two it writes for *nothing*.

    A list per row rather than a letter per row, because that is the format: `['B', 'C', 'D']`
    is one variant's three CIViC evidence items, and reading it as one value is issue #109.
    """
    cells = [cell for cell, _levels in _REAL_CIVIC_CELLS] + ["[]", ""]
    return pd.DataFrame(
        {"CIViC_Evidence_Level": cells, "Hugo_Symbol": ["KRAS"] * len(cells)}
    )


def _sources_for(columns, arm="somatic"):
    """What the overview would draw for a frame with these columns, as the app derives it.

    #95 settled the question this file's `_escat_circle` was written against: the circles
    follow the arm *and* the file, so a frame carrying only `ESCAT` draws one position. The
    tests below seed this rather than the six-source catalogue, because the catalogue is no
    longer what any report draws.
    """
    return clinical_summary.circle_sources(arm, list(columns))


def _escat_circle(level):
    """The circle drawn at the ES position for a row whose only annotation is this level.

    The position is *looked up*, not assumed: which positions exist depends on the arm and
    the file since #95, so it is not this test's to hardcode.
    """
    row = pd.Series({"ESCAT": level})
    sources = _sources_for(["ESCAT"])
    position = [entry[0] for entry in sources].index("ES")
    # The glyph *list*, not the joined string: a glyph can be two code points (issue #98),
    # so indexing the join stopped meaning "the nth source".
    return clinical_summary.pathogenicity_circle_glyphs(row, sources)[position]


def test_every_escat_level_the_scale_defines_has_a_circle():
    """The mapping is keyed on ESCAT's own groups, so it cannot silently miss a level.

    This is the guard the string-prefix version could not have: `V` matched none of its
    three branches and fell through, and nothing said so.
    """
    from config.vocabularies import ESCAT_DEFINITIONS

    groups = {meaning.group for meaning in ESCAT_DEFINITIONS.values()}
    unmapped = groups - set(clinical_summary._ESCAT_GROUP_CIRCLES)

    assert not unmapped, (
        f"ESCAT groups with no circle: {sorted(unmapped)} — a level in this group would "
        "draw the unreadable glyph while the help page defines it"
    )
    assert set(clinical_summary._ESCAT_GROUP_CIRCLES) == groups, (
        "the circle map names a group `ESCAT_DEFINITIONS` does not, so it describes a "
        "level the app cannot receive"
    )


def test_a_level_draws_the_circle_its_own_group_earns():
    """Derived from `ESCAT_DEFINITIONS`, not from the shape of the level's name."""
    from config.vocabularies import ESCAT_DEFINITIONS

    for level, meaning in ESCAT_DEFINITIONS.items():
        expected = clinical_summary._CLASS_GLYPHS[
            clinical_summary._ESCAT_GROUP_CIRCLES[meaning.group]
        ]
        assert _escat_circle(level) == expected, (
            f"{level} is in ESCAT's '{meaning.group}' group and drew something else"
        )


def test_the_combination_development_level_is_not_drawn_as_a_blank():
    """`V` is the defect this ticket exists for.

    Its definition says a matched drug produces objective responses without improving
    outcome — a real target under combination development. It drew `⚪`, a glyph the key did
    not list, immediately beside the `⬜` the key calls "No data". Both would have reported
    that the app has nothing to say about it.
    """
    drawn = _escat_circle("V")

    assert drawn == clinical_summary._CLASS_GLYPHS["Drug_Response"], (
        f"ESCAT V drew {drawn!r}; its definition is about a drug that acts, so it draws "
        "the drug glyph"
    )
    assert drawn not in (
        clinical_summary._UNREADABLE_CIRCLE,
        clinical_summary._ABSENT_CIRCLE[0],
    ), "ESCAT V is being reported as absent or unreadable, which is what #100 fixed"


def test_a_level_the_scale_does_not_define_is_not_drawn_as_the_strongest():
    """`IV` used to take the `Pathogenic` branch on `"IV".startswith("I")`.

    ESCAT's `IV` is preclinical evidence only — its *weakest* actionable level — so the old
    mapping drew it indistinguishably from `IA`. It was unreachable because
    `resources/escat_tiering.csv` assigns none, which is a fact about a data file rather
    than a check this function made; now the function makes the check.
    """
    for undefined in ("IV", "X", "['IA']"):
        drawn = _escat_circle(undefined)
        assert drawn != clinical_summary._CLASS_GLYPHS["Pathogenic"], (
            f"{undefined!r} drew the strongest circle, which the scale does not support"
        )
        assert drawn == clinical_summary._UNREADABLE_CIRCLE


def test_the_seven_reachable_levels_still_draw_what_they_drew():
    """This was a repair, not a re-grading.

    Measured on four real annotated MAFs before the change: `IA`/`IB`/`IC` red, `IIA`/`IIB`
    orange, `IIIA`/`IIIB` yellow. The group derivation reproduces every one of them and
    moves only `V`, so a future edit to `_ESCAT_GROUP_CIRCLES` that quietly re-grades the
    scale fails here rather than on a clinician's screen.
    """
    expected = {
        "IA": "🔴",
        "IB": "🔴",
        "IC": "🔴",
        "IIA": "🟠",
        "IIB": "🟠",
        "IIIA": "🟡",
        "IIIB": "🟡",
    }
    assert {level: _escat_circle(level) for level in expected} == expected


def test_the_key_names_every_glyph_the_circles_can_draw(broken_grid, escat_rows):
    """The claim the old key could not make, and did not keep.

    Every glyph is collected by *drawing* — over every ESCAT level plus a row for each
    entry in `CLINICAL_VALUE_MAPPING` and an unmapped ClinVar term — and then looked for in
    what the grid rendered. The transcribed key named six of the eight.
    """
    rows = [pd.Series({"ESCAT": level}) for level in list(escat_rows["ESCAT"])]
    rows += [
        pd.Series({"ClinVar_VCF_CLNSIG": value})
        for value in clinical_summary.CLINICAL_VALUE_MAPPING
    ]
    # CIViC is read on its own scale since #109, so its glyphs are not reachable through the
    # mapping above — a level whose circle the key does not explain would otherwise pass here.
    rows += [
        pd.Series({"CIViC_Evidence_Level": cell})
        for cell in [level for level in clinical_summary._CIVIC_LEVEL_CLASSES] + ["[]", "?"]
    ]
    rows += [pd.Series({"ClinVar_VCF_CLNSIG": "risk_factor"}), pd.Series({"Hugo_Symbol": "TP53"})]

    # The glyph *list*, never the joined string: iterating the join yields U+FE0F as if it
    # were a glyph, so a two-code-point glyph would be collected as two broken halves (#98).
    drawn = set()
    for row in rows:
        sources = _sources_for(row.index)
        drawn.update(clinical_summary.pathogenicity_circle_glyphs(row, sources))
        if not sources:
            # A row with no source column draws no circles at all since #95, so the ⬜ this
            # test needs to see collected has to come from a row that does draw a position.
            drawn.update(
                clinical_summary.pathogenicity_circle_glyphs(row, _sources_for(["ESCAT"]))
            )

    broken_grid.session_state = {
        "overview_sources": _sources_for(escat_rows.columns)
    }
    variant_table._render_aggrid_with_detail(
        escat_rows, escat_rows, 500, "passed_variants", ["Pathogenicity_Overview", "ESCAT"]
    )
    key = _key_text(broken_grid)

    undocumented = {glyph for glyph in drawn if glyph not in key}
    assert not undocumented, (
        f"the circles draw {sorted(undocumented)} and the key does not explain them — "
        "this is exactly how ⚪ and 💊 reached a clinician's screen unexplained (#100)"
    )


def test_the_key_explains_no_glyph_the_circles_cannot_draw(broken_grid, escat_rows):
    """The other direction, so deriving the key cannot document a state that never occurs.

    #88 and #89 both drew this line: an option that can keep nothing, or a level the
    annotation cannot carry, is not a gap to fill.

    Reachability is read off the two *mappings*, never off the glyph table — the legend is
    derived from that table, so comparing the legend against it compares the legend with
    itself. Mutation-checking caught exactly that: an invented `🟣 Resistant` entry passed,
    because adding it to the legend added it to the map as well (#77's `deepcopy(X) == X`,
    in a new place). #104 moved the table under both columns, which changes which constant
    must not be read here and not the reason.
    """
    broken_grid.session_state = {
        "overview_sources": _sources_for(escat_rows.columns)
    }
    variant_table._render_aggrid_with_detail(
        escat_rows, escat_rows, 500, "passed_variants", ["Pathogenicity_Overview", "ESCAT"]
    )
    key = _key_text(broken_grid)

    reachable = (
        set(clinical_summary.CLINICAL_VALUE_MAPPING.values())
        | set(clinical_summary._ESCAT_GROUP_CIRCLES.values())
        # Three mappings now, CIViC having left the shared one at #109. Named here even though
        # its classes are all reachable through ClinVar too: a mapping that draws and is not
        # read here is how this derivation would start comparing the legend with a subset of
        # what the circles can say.
        | set(clinical_summary._CIVIC_LEVEL_CLASSES.values())
    )
    for name, (glyph, label) in clinical_summary._CIRCLE_LEGEND.items():
        assert name in reachable, (
            f"the key explains {glyph} {label!r}, and no annotation value maps to "
            f"{name!r} — it documents a circle the app cannot draw"
        )
        assert label in key, f"the key drew {glyph} without the words that explain it"

    absent_glyph, absent_label = clinical_summary._ABSENT_CIRCLE
    assert absent_label in key, (
        f"the key does not explain {absent_glyph}, the circle drawn for a source with "
        "nothing to say — the one glyph that is not a classification"
    )


def test_the_key_says_which_axis_the_escat_circle_is_on(broken_grid, escat_rows):
    """One of six sources is graded on a different scale from the other five.

    ESCAT ranks how actionable a *target* is; the rest classify how pathogenic a *variant*
    is, and on real files they disagree — 34 of 48 ESCAT-annotated rows on one measured MAF
    draw 🔴 or 🟠 here while every other populated source on the row says benign. The
    circle stays; the key has to say what it means.

    Since #95 the note follows the ES *position*: a report that draws no ESCAT circle
    should not carry a sentence explaining one, which
    :func:`test_the_escat_note_is_not_drawn_without_an_escat_circle` holds from the other
    side.
    """
    broken_grid.session_state = {
        "overview_sources": _sources_for(escat_rows.columns)
    }
    variant_table._render_aggrid_with_detail(
        escat_rows, escat_rows, 500, "passed_variants", ["Pathogenicity_Overview", "ESCAT"]
    )

    assert clinical_summary.ESCAT_CIRCLE_NOTE in _key_text(broken_grid)


def test_the_escat_note_is_not_drawn_without_an_escat_circle(broken_grid, escat_rows):
    """A germline MAF with no ESCAT column draws no ES position, so the note has no subject.

    The note is a sentence about one of the circles. Drawn beside a key that does not list
    ES, it would explain a position the reader cannot find — the shape #95 removed from the
    Order line, arriving one line lower.
    """
    broken_grid.session_state = {
        "overview_sources": _sources_for(["ClinVar_VCF_CLNSIG"], arm="germline")
    }
    variant_table._render_aggrid_with_detail(
        escat_rows, escat_rows, 500, "passed_variants", ["Pathogenicity_Overview", "ESCAT"]
    )

    key = _key_text(broken_grid)
    assert "Pathogenicity Overview" in key, "the key was not drawn at all"
    assert clinical_summary.ESCAT_CIRCLE_NOTE not in key


# ---------------------------------------------------------------------------
# The CIViC position: the format it reads, and the axis it is graded on (issue #109)
# ---------------------------------------------------------------------------
#
# Every value used below is copied out of real annotated MAF or this repo's `somatic_civic.maf`
# fixture, because the defect was a *format* nobody had looked at: the annotation writes one
# entry per CIViC evidence item, so a cell holds a stringified list and the app was reading the
# first characters of a list literal. 133 of 140 CIViC-annotated rows across 34 real MAFs drew
# the unreadable circle, and a test built on invented single letters would have passed
# throughout — the seven rows that worked are exactly the bare-value ones.

#: Cells lifted verbatim from real files, with the levels each one mentions.
#:
#: Verbatim, and the count beside each is how many rows of the measured corpus carry that exact
#: cell — because two of these were **invented** on the first pass, in a table whose comment
#: claimed they were real: `"['B', 'C', 'D']"` and `"['B', 'B', 'B', 'C', 'C', 'D']"` are
#: plausible, are what an abbreviated example looks like, and occur **0 times** in 297 MAFs.
#: `/code-review` caught it. The repeated levels are the point of the real ones — a cell holds
#: one entry per evidence item, so a variant with nine items whose levels are three distinct
#: values carries all nine — and an example that quietly tidies that away is testing the format
#: this ticket wishes the column had.
_REAL_CIVIC_CELLS = [
    ("['B', 'C', 'D', 'D', 'D']", ["B", "C", "D"]),  # 3 rows
    ("['A']", ["A"]),  # 1 row
    ("['A', 'C']", ["A", "C"]),  # 7 rows
    ("['D', 'E']", ["D", "E"]),  # 1 row
    ("['A', 'B', 'B', 'C', 'D', 'D', 'D']", ["A", "B", "C", "D"]),  # 5 rows
    ("['B', 'B', 'C', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D']", ["B", "C", "D"]),  # 6
    ("A", ["A"]),  # 6 rows
    ("C", ["C"]),  # 2 rows
]


def _civic_circle(value):
    """The circle drawn at the CI position for a row whose only annotation is this cell."""
    row = pd.Series({"CIViC_Evidence_Level": value})
    sources = _sources_for(["CIViC_Evidence_Level"])
    position = [entry[0] for entry in sources].index("CI")
    return clinical_summary.pathogenicity_circle_glyphs(row, sources)[position]


@pytest.mark.parametrize("cell,levels", _REAL_CIVIC_CELLS)
def test_a_civic_cell_is_read_as_every_level_it_mentions(cell, levels):
    """The format the column is actually written in, bare values and lists alike."""
    assert clinical_summary.civic_levels(cell) == levels


@pytest.mark.parametrize("cell,levels", _REAL_CIVIC_CELLS)
def test_the_civic_reading_is_the_vendored_filters_own_rule(cell, levels):
    """Run the matcher rather than restate the rule — #99's line, and the point of the fix.

    ``has_element_from_list`` is what the pipeline asks of this column, and it is why CIViC
    *filtering* worked on the very files where the circle failed. The report's reading is that
    rule applied per level, so the two agree by construction; this asserts it against the frozen
    matcher instead of against a second description of it, so a change to either side fails.

    A unit test with a vendored function as its oracle, filed in the unit suite rather than under
    ``guard``: ``tests/README.md``'s own table puts ``test_numeric_columns.py`` here on exactly
    that footing — *"the refusal, asserted as a biconditional with the vendored functions as
    oracle"*. What the guards hold is ``vendor/`` against ``bin/``; what this holds is a
    component against ``vendor/``.

    The levels are ``CIVIC_OPTIONS``, never a literal list: a table of copied levels is how a
    parametrised guard silently stops covering the thing it names, which is a defect this repo
    has met before.
    """
    from config.vocabularies import CIVIC_OPTIONS
    from vendor.pipeline_filters import has_element_from_list

    assert [
        level for level in CIVIC_OPTIONS if has_element_from_list(cell, [level])
    ] == levels


@pytest.mark.parametrize(
    "cell,glyph",
    [
        ("['A', 'C']", "🔴"),
        ("['B', 'C', 'D']", "🟠"),
        ("['C', 'D', 'E']", "🟡"),
        ("['D', 'E']", "🟡"),
    ],
)
def test_the_strongest_level_on_the_row_decides_the_civic_circle(cell, glyph):
    """One circle, several evidence items — 110 of 140 measured rows mention two or more.

    The strongest is the reading the filter already takes (a row is kept when *any* level it
    mentions is selected) and the reading the summary beside it takes across sources, so a
    third answer here would be the app disagreeing with itself about one cell.
    """
    assert _civic_circle(cell) == glyph


@pytest.mark.parametrize("level", CIVIC_OPTIONS)
def test_no_civic_level_draws_a_benign_circle(level):
    """The falsity this ticket removes, held from the outside.

    Parametrised over the vocabulary rather than over a copy of it, so a level added to the
    control is covered here without anyone remembering — a table of copied values is how a
    parametrised guard stops applying to what its name claims.

    ``E`` used to map to ``Benign``, so a variant whose CIViC evidence was merely *indirect*
    drew 🔵 beside sources that had said nothing of the kind. No CIViC level asserts that a
    variant is harmless, so none of them may draw a benign glyph — asserted over the glyphs
    rather than over the mapping, because the mapping is what a change would edit.
    """
    benign_glyphs = {
        clinical_summary._CLASS_GLYPHS["Benign"],
        clinical_summary._CLASS_GLYPHS["Likely_Benign"],
    }
    assert _civic_circle(level) not in benign_glyphs, (
        f"CIViC level {level} draws a benign circle; CIViC's levels grade the evidence, and "
        "weak evidence is not a benign call"
    )


def test_a_civic_cell_holding_no_evidence_items_draws_no_data():
    """``[]`` is CIViC saying *no evidence*, which is ⬜'s meaning and not ⚪'s.

    Both are drawn for a cell mentioning no level, and the difference matters on screen: ⬜ says
    this source has no call for this variant, ⚪ says the app could not read what it said.
    """
    assert _civic_circle("[]") == clinical_summary._ABSENT_CIRCLE[0]
    assert _civic_circle("[]") != clinical_summary._UNREADABLE_CIRCLE


def test_a_civic_value_this_app_cannot_read_is_still_unreadable():
    """The other half of that pair, so the empty-list case cannot be widened into a catch-all.

    A cell mentioning no level and spelled like nothing the annotation writes has not been
    read, and saying "No data" about it would be the mistake this ticket removes in reverse.

    **This pins the symbol, not the glyph.** Issue #109's third bullet asked whether a value the
    app failed to parse should be reported with a glyph at all, and left it where it already sat:
    issue #104 owns what an unclassifiable value draws, and is open. So the assertion is against
    ``_UNREADABLE_CIRCLE`` rather than against ``⚪``, and #104 can change what that is — what
    this holds is only that an unreadable CIViC cell is not quietly reclassified as *no data*,
    which is this ticket's own distinction and not #104's to inherit.
    """
    assert _civic_circle("Level unknown") == clinical_summary._UNREADABLE_CIRCLE


def test_every_civic_level_the_control_offers_is_one_the_report_can_draw():
    """#88's rule, applied to CIViC: a level a user may select is one the report can name.

    Both directions, so a level cannot be offered without a circle *or* mapped without being
    offered — the pair of defects #88 and #103 each spent a ticket on, in one assertion.
    """
    from config.vocabularies import CIVIC_DEFINITIONS, CIVIC_OPTIONS

    assert set(CIVIC_OPTIONS) == set(clinical_summary._CIVIC_LEVEL_CLASSES)
    assert set(CIVIC_OPTIONS) == set(CIVIC_DEFINITIONS), (
        "the control offers a level the help page cannot define, or defines one it does not "
        "offer"
    )


def test_the_civic_levels_run_strongest_first():
    """The claim the help page makes about the order it renders them in.

    Checkable rather than trusted: the definitions are rendered in ``CIVIC_DEFINITIONS`` order
    and described as strongest first, so walking that order must never move *up* the clinical
    hierarchy. ``ESCAT_OPTIONS`` was in its paper's order by accident and asserted nowhere
    (issue #89); this is the same guard before the same accident.
    """
    from config.vocabularies import CIVIC_DEFINITIONS, CIVIC_OPTIONS

    assert list(CIVIC_DEFINITIONS) == CIVIC_OPTIONS
    ranks = [
        clinical_summary.CLINICAL_HIERARCHY.index(clinical_summary._CIVIC_LEVEL_CLASSES[level])
        for level in CIVIC_DEFINITIONS
    ]
    assert ranks == sorted(ranks), (
        f"the levels are rendered as strongest first and map to {ranks} in the hierarchy — a "
        "later level reads as stronger than an earlier one"
    )


def test_the_summary_and_the_circle_read_a_civic_cell_the_same_way():
    """One read, two readers — the seam the two columns came apart on.

    ``Clinical_Summary`` had the same broken reading of this column and the same repair, so a
    variant whose only annotation is CIViC must be labelled with the class its circle draws.
    Measured over 140 CIViC-annotated rows of real MAF, this moves **no** label — every one of
    them carries a stronger call from another source — so the state under test here is real but
    unobserved in that corpus, which is exactly why it needs a test rather than a measurement.

    **It asks the circle rather than recomputing what the circle should say.** The first version
    looked up ``_CIVIC_LEVEL_CLASSES[levels[0]]`` on both sides of the comparison, so it held the
    summary against a mapping and not against the other column: a divergence introduced in the
    CI branch of ``pathogenicity_circle_glyphs`` — the very place this repair lives — would have
    left it green. `/code-review` caught that, and it is #77's ``deepcopy(X) == X`` in a third
    place: a cross-column agreement test that never asks the second column.
    """
    for cell, _levels in _REAL_CIVIC_CELLS:
        row = pd.Series({"CIViC_Evidence_Level": cell})
        label = clinical_summary.generate_clinical_summary(row)
        summarised_as = clinical_summary.strip_summary_glyph(label)
        drawn = _civic_circle(cell)

        # The two columns speak different vocabularies by design (#104 owns whether they should),
        # so what is compared is the class underneath: the glyph the strip drew for this cell must
        # be the glyph of the class the label names.
        assert drawn == clinical_summary._CLASS_GLYPHS[
            next(
                name
                for name, rendered in clinical_summary._SUMMARY_LABELS.items()
                if clinical_summary.strip_summary_glyph(rendered) == summarised_as
            )
        ], (
            f"a variant annotated only {cell!r} is summarised as {label!r} while its own circle "
            f"draws {drawn} — the two columns are reading one cell two ways again"
        )


def test_a_civic_cell_with_no_evidence_items_is_not_an_annotation():
    """``[]`` reaches the summary too, and it is not evidence that CIViC said something.

    Left to the general path it would be an unreadable *value*, so a variant no source
    annotated would be labelled ``Unrecognised Annotation`` on the strength of an empty list.
    """
    row = pd.Series({"CIViC_Evidence_Level": "[]"})

    assert clinical_summary.generate_clinical_summary(row) == clinical_summary.NO_CLINICAL_DATA


def test_the_key_says_which_axis_the_civic_circle_is_on(broken_grid, civic_rows):
    """The second position on the strip that is not read on the pathogenicity ladder.

    ESCAT grades a target and CIViC grades *evidence*, so four of six positions classify the
    variant and two do not. Both clauses are drawn, and each only where its position is.
    """
    broken_grid.session_state = {"overview_sources": _sources_for(civic_rows.columns)}
    variant_table._render_aggrid_with_detail(
        civic_rows, civic_rows, 500, "passed_variants",
        ["Pathogenicity_Overview", "CIViC_Evidence_Level"],
    )

    assert clinical_summary.CIVIC_CIRCLE_NOTE in _key_text(broken_grid)


def test_switching_civic_filtering_off_takes_the_position_and_its_clause():
    """The case where the column is *present* and the position still must not exist.

    A user can switch CIViC filtering off on a MAF that carries CIViC, and ``compute_keep``
    then drops those columns from the report — so #95's ``skip_civic`` is the one route to a
    drawn-column-less position that is not "this file lacks the column". The clause follows the
    sources rather than the columns, so it goes with the circle; asserted here because every
    other test of the pair reaches the no-position state by taking the column away.
    """
    columns = ["ClinVar_VCF_CLNSIG", "CancerVar", "CIViC_Evidence_Level", "ESCAT"]
    sources = clinical_summary.circle_sources("somatic", columns, skip_civic=True)

    assert [entry[0] for entry in sources] == ["CV", "CA", "ES"]
    assert clinical_summary.CIVIC_CIRCLE_NOTE not in clinical_summary.circle_axis_notes(sources)
    assert clinical_summary.ESCAT_CIRCLE_NOTE in clinical_summary.circle_axis_notes(sources)


def test_the_civic_note_is_not_drawn_without_a_civic_circle(broken_grid, civic_rows):
    """A MAF with no CIViC column draws no CI position, so the clause has no subject.

    #95's line, one line lower and for the second source: a sentence about a circle the reader
    cannot find explains nothing.
    """
    broken_grid.session_state = {
        "overview_sources": _sources_for(["ClinVar_VCF_CLNSIG"], arm="germline")
    }
    variant_table._render_aggrid_with_detail(
        civic_rows, civic_rows, 500, "passed_variants",
        ["Pathogenicity_Overview", "CIViC_Evidence_Level"],
    )

    key = _key_text(broken_grid)
    assert "Pathogenicity Overview" in key, "the key was not drawn at all"
    assert clinical_summary.CIVIC_CIRCLE_NOTE not in key


def test_no_key_is_drawn_for_a_column_that_is_not_shown(broken_grid, escat_rows):
    """The key belongs to `Pathogenicity_Overview`, and follows it off the screen.

    The sources are seeded, so this fails if the column check is removed rather than
    passing because nothing was carried for the key to name.
    """
    broken_grid.session_state = {
        "overview_sources": _sources_for(escat_rows.columns)
    }
    variant_table._render_aggrid_with_detail(
        escat_rows, escat_rows, 500, "passed_variants", ["ESCAT"]
    )

    assert "Pathogenicity Overview" not in _key_text(broken_grid)


# ---------------------------------------------------------------------------
# What a ClinVar term that is not a pathogenicity call is reported as (issue #98)
# ---------------------------------------------------------------------------
#
# The ticket was written about `other` and the measurement moved it: across 62 byte-distinct
# real annotated MAFs (147,152 rows) `other` never stands alone — every occurrence is a
# modifier on a cell that also carries a real call, which the split resolves — so it was
# reported as "No Clinical Data" on **zero** rows. What *was* reported that way is five
# sibling terms the control does not even offer.
#
# The figures, kept together here so no two docstrings can drift apart: 121 rows change
# summary; 31 of those are the `Not Provided` -> `No Classification` rename; 89 leave "No
# Clinical Data" (82 carrying a ClinVar value, 7 a `CancerVar` value of `1`); and 76 of the
# 89 carry one of the six terms below — risk_factor 45, protective 13, association 12,
# Affects 3, confers_sensitivity 3, other 0.

_NON_PATHOGENICITY_TERMS = [
    "risk_factor",
    "association",
    "protective",
    "Affects",
    "confers_sensitivity",
    "other",
]

#: The ACMG/AMP ladder — the classes at or above `Benign` in the hierarchy. Nothing in the
#: list above may be reported as one of these.
#:
#: Sliced from `CLINICAL_HIERARCHY` rather than typed out, so it cannot become the second
#: copy of the class list that this change removed from production code.
_SEVERITY_LABELS = {
    clinical_summary._SUMMARY_LABELS[name]
    for name in clinical_summary.CLINICAL_HIERARCHY[
        : clinical_summary.CLINICAL_HIERARCHY.index("Benign") + 1
    ]
}


def _summary_of(**annotations):
    return clinical_summary.generate_clinical_summary(pd.Series(annotations))


@pytest.mark.parametrize("term", _NON_PATHOGENICITY_TERMS)
def test_a_clinvar_assertion_is_classified_rather_than_called_absent(term):
    """The defect this ticket exists for, held per term.

    Measured before the fix: 76 rows carried one of these and were summarised as having no
    clinical data — said of a variant whose only distinguishing feature is the ClinVar
    assertion it carries. `protective` is the sharpest, being a claim in the opposite
    direction from the one the reader would infer from silence.

    Both halves are asserted, and the second is the one that bites. Since the fall-through
    was repaired, deleting a term from the mapping no longer produces "No Clinical Data" —
    it produces the *unrecognised* label, so a test checking only for the old sentence would
    pass over exactly the regression it was written to catch. Mutation-checking found that:
    removing `protective` from the mapping left the first assertion green.
    """
    summary = _summary_of(ClinVar_VCF_CLNSIG=term)

    assert summary != clinical_summary.NO_CLINICAL_DATA, (
        f"a variant ClinVar calls {term!r} is reported as carrying no clinical data"
    )
    assert summary != clinical_summary._SUMMARY_LABELS["Unknown"], (
        f"{term!r} has no class of its own, so the report can only say it did not "
        "understand the annotation — true, but the term table says what this one means"
    )


@pytest.mark.parametrize("term", _NON_PATHOGENICITY_TERMS)
def test_no_non_pathogenicity_term_is_reported_as_a_severity(term):
    """The other half, and the reason #88 refused to fix this by picking a bucket.

    Repairing the sentence by mapping these onto the ACMG ladder would trade a false
    *absence* for a false *classification* — `protective` as Benign, `risk_factor` as
    Pathogenic — which is a clinical claim about a patient's variant that no source made.
    ClinVar groups these terms separately from the ACMG/AMP five, and the institute's term
    table gives each its own class; both keep them off the ladder.

    Asserted on what the summary *renders*, not on the mapping entry. Reading the dict would
    pass over a term whose entry the split makes unreachable — which is not hypothetical,
    it is what the low-penetrance test below caught in this very change.
    """
    summary = _summary_of(ClinVar_VCF_CLNSIG=term)

    assert summary not in _SEVERITY_LABELS, (
        f"a variant ClinVar calls {term!r} is reported as {summary!r}, a point on the "
        "pathogenicity ladder. No source says that; it is the app inventing a severity"
    )


def test_terms_from_different_classes_are_told_apart():
    """A shared bucket would have merged claims that point opposite ways.

    The first design for this fix was one 'other ClinVar assertion' class. The term table
    refused it: `protective` decreases risk, `risk_factor` increases it, and `association`
    is a GWAS finding rather than either. One label for all three would have been true of
    none of them.
    """
    labels = {
        term: _summary_of(ClinVar_VCF_CLNSIG=term)
        for term in ("protective", "risk_factor", "association")
    }

    assert len(set(labels.values())) == 3, (
        f"these three terms share a label: {labels}. They are different classes in the "
        "source that assigns them"
    )


def test_other_is_classed_as_having_no_classification():
    """The ticket's title question, answered from the term table rather than from taste.

    ClinVar defines `other` as "ClinVar does not have the appropriate term for your
    submission", and the institute's table excludes it because its meaning depends on a
    free-text explanation this app never receives. Both say: there is a record, and no
    usable classification in it.
    """
    assert _summary_of(ClinVar_VCF_CLNSIG="other") == clinical_summary._SUMMARY_LABELS[
        "No_Classification"
    ]
    assert _summary_of(ClinVar_VCF_CLNSIG="other") == _summary_of(
        ClinVar_VCF_CLNSIG="not_provided"
    )


def test_no_clinical_data_is_said_only_when_nothing_annotated_the_variant():
    """The sentence has to keep meaning something, or the repair is just a rename.

    A row no source annotated still says it; a row with an annotation nothing can read does
    not, and says so in its own words instead.
    """
    nothing = _summary_of(Hugo_Symbol="TP53")
    unreadable = _summary_of(ClinVar_VCF_CLNSIG="a_term_clinvar_has_never_issued")

    assert nothing == clinical_summary.NO_CLINICAL_DATA
    assert unreadable != clinical_summary.NO_CLINICAL_DATA
    assert unreadable == clinical_summary._SUMMARY_LABELS["Unknown"]


def test_an_unreadable_value_in_any_source_is_not_called_absent():
    """The repair is not ClinVar's alone, because the discard never was.

    All five columns went through the same drop. Measured on the corpus this moves 7 rows
    whose `CancerVar` cell is `1` — an annotation artefact, and "unrecognised" is true of it
    where "no clinical data" was not.
    """
    assert _summary_of(CancerVar="1") != clinical_summary.NO_CLINICAL_DATA


def test_a_low_penetrance_call_is_not_reported_as_a_fully_penetrant_one():
    """The term table's one explicit warning, and the split was breaking it.

    `Pathogenic,_low_penetrance` carries a comma, which is one of the separators the
    multi-value split uses — so the cell was cut down to `Pathogenic` before the mapping saw
    it, asserting the equivalence the table says must not be asserted. Caught by running the
    mapping rather than reading it: the entries for these two terms were unreachable, and a
    test that only asserted the dict would have passed over that.
    """
    for term, plain in (
        ("Pathogenic,_low_penetrance", "Pathogenic"),
        ("Likely_pathogenic,_low_penetrance", "Likely_pathogenic"),
    ):
        summary = _summary_of(ClinVar_VCF_CLNSIG=term)
        assert summary != _summary_of(ClinVar_VCF_CLNSIG=plain), (
            f"{term!r} is reported exactly as {plain!r}, so the penetrance qualifier the "
            "term carries is dropped between the file and the screen"
        )
        assert "low penetrance" in summary


def test_a_genuine_composite_is_still_split():
    """The counterpart, pinning how narrow the fix above is.

    Skipping the split for *any* value the mapping knows would also route
    `Benign/Likely_benign` through its own entry — a defensible reading that moves 2,088
    rows from Benign to Likely Benign on the measured corpus. That is a clinical
    reclassification this ticket did not set out to make, and this test fails if a later
    change makes it in passing.
    """
    assert _summary_of(ClinVar_VCF_CLNSIG="Benign/Likely_benign") == _summary_of(
        ClinVar_VCF_CLNSIG="Benign"
    )


def test_every_class_has_a_label_a_rank_and_a_colour():
    """Three lists that must agree, and did not: the else-branch classes were in none of them.

    The renderer, the sort order and the chart axis each held their own copy of the label
    set. A class missing from the second sorted to the bottom; missing from the third, it
    lost its place in the axis order silently. Both are invisible failures, so they are
    asserted rather than hoped for.
    """
    for name in clinical_summary.CLINICAL_HIERARCHY:
        assert name in clinical_summary._SUMMARY_LABELS, f"{name} renders as nothing"

        label = clinical_summary._SUMMARY_LABELS[name]
        assert label in clinical_summary._SORT_ORDER, f"{label} has no sort rank"
        assert clinical_summary.strip_summary_glyph(label) in charts.CLINICAL_COLORS, (
            f"{label} is not a category of the clinical chart, so it drops to the "
            "unrecognised tail of that axis instead of its clinical place"
        )


#: Cells measured in real annotated MAFs whose pieces do not all carry the same call.
#:
#: A composite is where the two columns could still part company after #104: they agree about
#: a glyph by construction now, but each has to pick the *same piece* first, and until this
#: ticket they looked for a different substring while doing it.
_MEASURED_COMPOSITES = [
    "Benign/Likely_benign",
    "Uncertain_significance|_drug_response",
    "Affects|_association",
    "Benign|_other",
    "Pathogenic/Likely_pathogenic|other",
    "Benign;Pathogenic,_low_penetrance",
]


@pytest.mark.parametrize(
    "value",
    sorted(clinical_summary.CLINICAL_VALUE_MAPPING)
    + _MEASURED_COMPOSITES
    + ["a_value_no_source_issues"],
)
def test_a_class_draws_the_same_glyph_in_both_clinical_columns(value):
    """Issue #104's answer, held where it can fail: on what the two columns *render*.

    `Clinical_Summary` and `Pathogenicity_Overview` are pinned side by side, and they used to
    spell two of the thirteen classes differently — ✅ Benign against 🔵 on **248,034 of
    330,189 measured rows**, and ❓ No Classification against ⚪ on 427, the same words under
    two glyphs. Both now read `_CLASS_GLYPHS`.

    **Asserting that the two dicts agree would prove nothing**, since both are derived from
    that one table — the shape #100 was caught by, where the key was compared against a map
    built out of the key. So this runs the two columns instead and compares their output: a
    row carrying one value, summarised by one function and drawn by another.

    The last case is a value in no source's vocabulary, which is the pair the ticket was
    opened for — ⚪ Unrecognised Annotation beside ⚪ — and the composites are cells measured
    in real MAFs, which is where picking a *different piece* would show up.
    """
    row = pd.Series({"ClinVar_VCF_CLNSIG": value})
    sources = _sources_for(["ClinVar_VCF_CLNSIG"])

    summary = clinical_summary.generate_clinical_summary(row)
    circle = clinical_summary.pathogenicity_circle_glyphs(row, sources)[0]

    assert summary.partition(" ")[0] == circle, (
        f"{value!r} is drawn {circle} by the overview and labelled {summary!r} by the "
        "summary beside it — one class, two glyphs, on the same row"
    )


def test_the_shared_circle_is_not_explained_as_the_class_it_only_half_covers(
    broken_grid, escat_rows
):
    """⚪ is two claims, so the words beside it may not be either one of them alone.

    The dev's call on #104 was that an unreadable value keeps sharing `No_Classification`'s
    circle rather than earning an eleventh glyph: measured across 183 real annotated MAFs, ⚪
    is drawn at 1,108 positions and only 121 of them are unreadable — and after #98 every one
    of those 121 is a value the app cannot **parse** rather than a clinical term it declines
    to classify, so a glyph of its own would have been a permanent key row for other tickets'
    unfinished parsing. The two have never shared a cell: 0 rows draw ⚪ from both causes.

    What sharing costs is the wording, which is what this holds. "No classification" is a
    claim about what the *source* recorded and is false of a value that was merely unreadable
    — as "Not provided" before it was false of both. The premise is asserted too, so that
    splitting the glyph later fails this test rather than quietly emptying it.
    """
    assert (
        clinical_summary._CLASS_GLYPHS["Unknown"]
        == clinical_summary._CLASS_GLYPHS["No_Classification"]
    ), "the two no longer share a circle, so this test's reason for existing has changed"

    shared = clinical_summary._KEY_WORDS["No_Classification"]
    for name in ("No_Classification", "Unknown"):
        assert shared.casefold() != clinical_summary._SUMMARY_WORDS[name].casefold(), (
            f"the key explains ⚪ as {shared!r}, which is exactly what the app calls "
            f"{name!r} — one of the two cases that circle stands for, asserted over the other"
        )

    broken_grid.session_state = {"overview_sources": _sources_for(escat_rows.columns)}
    variant_table._render_aggrid_with_detail(
        escat_rows, escat_rows, 500, "passed_variants", ["Pathogenicity_Overview", "ESCAT"]
    )
    assert shared in _key_text(broken_grid)


def test_every_class_a_value_can_map_to_has_a_glyph():
    """The lookup in `pathogenicity_circle_glyphs` is total, and this is what makes it so.

    It used to end in `.get(mapped, _UNREADABLE_CIRCLE)`, which was load-bearing only because
    `Unknown` had no glyph. #104 gave it one, and a default that can no longer fire would
    have gone on quietly catching a class someone forgot to name — drawing it as unreadable
    rather than failing. So the default is gone and the totality is asserted instead.

    Every branch feeds the same table: ESCAT maps through its own groups and CIViC through
    its own levels — the two sources graded on their own scale (#100, #109) — and the other
    four through `CLINICAL_VALUE_MAPPING`. Every value any of them can produce has to be a
    class this app can draw.
    """
    produced = (
        set(clinical_summary.CLINICAL_VALUE_MAPPING.values())
        | set(clinical_summary._ESCAT_GROUP_CIRCLES.values())
        | set(clinical_summary._CIVIC_LEVEL_CLASSES.values())
        | {"Unknown"}
    )

    unnamed = produced - set(clinical_summary._CLASS_GLYPHS)
    assert not unnamed, (
        f"{sorted(unnamed)} can come out of a mapping and has no glyph, so a real annotation "
        "value would raise while the overview column is being derived"
    )


def test_no_source_spells_a_call_that_would_have_split_the_two_readings():
    """Why dropping `main_classification`'s marker parameter was safe, asserted not argued.

    The summary looked for `"pathogenic"` in a piece and the circles for `"athogenic"`, and
    #104 removed the difference. The argument for its being inert is structural — the second
    is a substring of the first and the piece is lowercased before either is sought — so the
    two can only ever part over a value carrying `athogenic` **without** `pathogenic`. That
    holds only while nobody adds such a term, which is what this checks: every call any of
    the six sources can spell, plus the ESCAT levels the circles read instead.

    Measured as well as reasoned: none of the 125 distinct values across 183 real annotated
    MAFs is one either. This holds the vocabularies, which is the half that can change here.
    """
    spellings = set(clinical_summary.CLINICAL_VALUE_MAPPING) | set(
        clinical_summary.ESCAT_DEFINITIONS
    )

    for value in spellings:
        piece = value.lower()
        assert not ("athogenic" in piece and "pathogenic" not in piece), (
            f"{value!r} contains 'athogenic' without 'pathogenic', so the marker the two "
            "columns once used would pick different pieces of a composite carrying it — "
            "the difference #104 removed as inert is no longer inert"
        )


def test_the_feature_doc_lists_the_classes_the_app_actually_renders():
    """`docs/CLINICAL_SUMMARY_FEATURE.md` names all thirteen classes, and nothing held it.

    It is the third place the glyph-and-words pairs are written, after the summary column and
    the key — and the only one a change to :data:`_CLASS_GLYPHS` does not reach. That is the
    shape this map keeps finding: #100's key named six of eight glyphs, #61's Home page
    catalogued annotation sources the app had stopped reading, and #79 corrected nineteen
    claims on a page nobody had checked against the code. Three of these lines were wrong the
    moment #104 changed a glyph, and the suite would have stayed green.

    Only the hierarchy section is read. The rest of the file describes each source's own
    vocabulary, which is `CLINICAL_VALUE_MAPPING`'s business and a different guard's.
    """
    doc = (
        Path(__file__).resolve().parent.parent / "docs" / "CLINICAL_SUMMARY_FEATURE.md"
    ).read_text(encoding="utf-8")

    section = doc.split("### Priority Hierarchy")[1].split("### Data Sources")[0]
    listed = re.findall(r"^(?:\d+\.|-) \*\*(.+?)\*\*", section, flags=re.MULTILINE)

    rendered = set(clinical_summary._SUMMARY_LABELS.values()) | {
        clinical_summary.NO_CLINICAL_DATA
    }
    for label in listed:
        assert label in rendered, (
            f"the feature doc calls a class {label!r} and the app renders no such label — "
            "a reader matching the doc against their screen finds neither in the other"
        )

    for name in clinical_summary.CLINICAL_HIERARCHY:
        assert clinical_summary._SUMMARY_LABELS[name] in listed, (
            f"{name} is a class the app can put a variant in and the doc's hierarchy does "
            "not list it"
        )

    # The section also quotes what the key calls the shared circle, which is a fourth copy of
    # `_KEY_WORDS["No_Classification"]` and the one thing here not read off a list item.
    shared = clinical_summary._KEY_WORDS["No_Classification"]
    assert shared in section, (
        f"the doc explains the shared circle without using the key's own words ({shared!r}), "
        "so a reader comparing the two finds the app calling it something else"
    )


def test_the_analysis_error_glyph_stands_for_nothing_else():
    """❔ named two things in one column until #104, which is the collision in miniature.

    It was `Unrecognised Annotation`'s glyph *and* the label shown when the derivation itself
    raised — so one column said ❔ for a value it could not read and ❔ for a frame that had
    thrown. Moving the unreadable case onto the circles' ⚪ leaves this glyph naming one
    thing, and this fails if a later change points a class back at it.
    """
    error_glyph = clinical_summary.ANALYSIS_ERROR.partition(" ")[0]

    assert error_glyph not in set(clinical_summary._CLASS_GLYPHS.values()), (
        f"{error_glyph} labels a failed derivation and now also a class, so the column says "
        "the same thing about a value it read and a frame it could not process"
    )


def test_the_chart_axis_runs_in_the_same_order_as_the_hierarchy():
    """One severity order, not two, now that there are thirteen classes rather than eight."""
    from_hierarchy = [
        clinical_summary.strip_summary_glyph(clinical_summary._SUMMARY_LABELS[name])
        for name in clinical_summary.CLINICAL_HIERARCHY
    ]
    axis = [name for name in charts.CLINICAL_COLORS if name in set(from_hierarchy)]

    assert axis == from_hierarchy, (
        "CLINICAL_COLORS orders the classes differently from CLINICAL_HIERARCHY, so the "
        "chart and the table disagree about which class is more severe"
    )


def test_one_glyph_per_source_survives_a_two_code_point_glyph():
    """The hazard `pathogenicity_circle_glyphs` exists for, which nothing else catches.

    `⚠️` and `🛡️` are a symbol plus U+FE0F. The moment one of those is drawn, the *joined*
    strip stops having one character per source, so indexing or iterating it silently
    reports the wrong source — and every existing circle test happens to miss it, because
    they set one annotation and leave the other positions on the single-code-point ⬜.

    So this row draws a two-code-point glyph at the ClinVar position and then asks about a
    later one. Found by mutation: re-joining and re-splitting the list left the suite green.
    """
    row = pd.Series({"ClinVar_VCF_CLNSIG": "protective", "ESCAT": "IA"})
    sources = _sources_for(["ClinVar_VCF_CLNSIG", "ESCAT"])
    positions = [entry[0] for entry in sources]

    glyphs = clinical_summary.pathogenicity_circle_glyphs(row, sources)

    assert len(glyphs) == len(positions), (
        f"{len(glyphs)} glyphs for {len(positions)} sources — the strip no longer lines up "
        "with the sources it is drawn from"
    )
    assert glyphs[positions.index("CV")] == clinical_summary._CLASS_GLYPHS["Protective"]
    assert glyphs[positions.index("ES")] == clinical_summary._CLASS_GLYPHS["Pathogenic"]

    # And the joined form is exactly these glyphs, in order — the column's own value.
    assert clinical_summary.generate_pathogenicity_circles(row, sources) == "".join(glyphs)


@pytest.mark.parametrize(
    "label",
    ["⚠️ Disease Risk", "🛡️ Protective", "🔴 Pathogenic", "🔗 Association or Trait"],
)
def test_a_glyph_of_any_length_is_stripped_cleanly(label):
    """`⚠️` and `🛡️` are two code points, and both old strippers assumed one.

    One sliced `label[2:]`, which left a leading space; the other matched a hand-written
    tuple of nine emoji, and a glyph missing from it was not stripped at all — which on the
    chart meant the category stopped matching `CLINICAL_COLORS`. Neither failure raises.
    """
    stripped = clinical_summary.strip_summary_glyph(label)

    assert stripped == stripped.strip(), f"{label!r} stripped to {stripped!r}, with padding"
    assert stripped[0].isalpha(), f"{label!r} stripped to {stripped!r}, glyph still attached"
    assert stripped in charts.CLINICAL_COLORS


# --- What a note is stored against (issue #127) --------------------------------------
#
# Four of `_variant_key`'s five components are core columns, so a MAF missing one of those is
# refused before any of this runs. `Tumor_Seq_Allele2` is the fifth and nothing requires it —
# the reference table says so in as many words — and its absence used to be spelled into the
# key as the empty string, which is a hole rather than an identity. Issue #67 makes a note
# institute-wide and about the variant, so the failure to avoid is a note read against the
# wrong variant, and these guard that direction rather than the convenience of having a key.


def _keyed_session(monkeypatch, notes=None):
    """A real-dict session state, wired into `variant_table`. See `noted_frame` on why."""
    session = FakeSessionState()
    session["variant_notes"] = dict(notes or {})
    session["custom_annotations"] = {}
    monkeypatch.setattr(variant_table.st, "session_state", session, raising=False)
    return session


@pytest.mark.parametrize(
    "alt, why",
    [
        (None, "the column is absent from the file"),
        ("", "the column is there and the cell is empty"),
        (float("nan"), "the cell is NaN, which reads back as the string 'nan'"),
    ],
)
def test_a_row_that_does_not_say_its_alt_allele_has_no_note_key(alt, why):
    """No identity, no key — and the check is on the *value*, not on the frame's columns.

    `str(row.get(name, ""))` renders an absent column and a missing value identically, so a
    guard over `frame.columns` would catch only the first of these three. Only one of them is
    visible from the header row, and all three leave the same hole.
    """
    row = pd.Series(
        {
            "Hugo_Symbol": "TP53",
            "Chromosome": "17",
            "Start_Position": 7577120,
            "Reference_Allele": "C",
        }
    )
    if alt is not None:
        row["Tumor_Seq_Allele2"] = alt

    assert variant_table._variant_key(row) is None, (
        f"a key was built where {why}, so it identifies a position rather than a variant"
    )


@pytest.mark.parametrize(
    "component", ["Hugo_Symbol", "Chromosome", "Start_Position", "Reference_Allele"]
)
def test_a_blank_core_value_holes_the_key_too(component):
    """The four components that are core columns, which is not the same as core *values*.

    A MAF without one of these columns is refused before any of this runs, so it is tempting to
    treat them as guaranteed — and the temptation is the bug: `_variant_key` reads values, and a
    present column can hold nothing. Measured at 0 blank values across 326,204 rows of 165 real
    MAFs, so this is belt-and-braces; what makes it worth asserting is that the first version of
    the dialog's message named `Tumor_Seq_Allele2` for every hole, and would have sent a user
    with a blank `Hugo_Symbol` to check a column that was fine.
    """
    row = pd.Series(
        {
            "Hugo_Symbol": "TP53",
            "Chromosome": "17",
            "Start_Position": 7577120,
            "Reference_Allele": "C",
            "Tumor_Seq_Allele2": "T",
        }
    )
    row[component] = ""

    assert variant_table._variant_key(row) is None
    assert variant_table._missing_key_components(row) == (component,), (
        "the missing component must be named, not merely detected — the dialog quotes it"
    )


@pytest.mark.parametrize("allele", ["-", "."])
def test_a_value_a_maf_legitimately_writes_is_not_read_as_missing(allele):
    """`-` is how a MAF spells the absent side of an indel, and `.` the pipeline's own blank.

    Both are real values in these two columns, so treating them as holes would refuse a note
    on variants that are perfectly well identified — every insertion, for a start. Asserted
    rather than left to the reader of `config.missing_values.is_blank`, because "looks like
    nothing" is exactly the intuition that would add them — and because `.` is now read as
    missing by that module's *other* predicate, which is the one everything a user reads asks
    (issue #131). This is the assertion that keeps the note key on the right one.
    """
    row = pd.Series(
        {
            "Hugo_Symbol": "TP53",
            "Chromosome": "17",
            "Start_Position": 7577120,
            "Reference_Allele": "C",
            "Tumor_Seq_Allele2": allele,
        }
    )

    assert variant_table._variant_key(row) == f"TP53:17:7577120:C>{allele}"


def test_two_variants_at_one_position_do_not_share_one_note(monkeypatch):
    """The misattachment itself, in the column the user reads.

    Two different variants at one position, in a file that does not record the alternate
    allele. Both used to key to `TP53:17:7577120:C>`, so a note written against either was
    shown against both — on this file and on every later file carrying that position. The
    seeded note is under exactly the key the old code built, which is what makes this fail
    before the fix rather than merely pass after it.
    """
    _keyed_session(monkeypatch, {"TP53:17:7577120:C>": "hotspot, discussed at board"})
    frame = pd.DataFrame(
        {
            "Hugo_Symbol": ["TP53", "TP53"],
            "Chromosome": ["17", "17"],
            "Start_Position": [7577120, 7577120],
            "Reference_Allele": ["C", "C"],
            # The alt allele is what tells these two apart, and the file does not carry it.
            "Protein_Change": ["p.R273H", "p.R273C"],
        }
    )

    noted = variant_table._add_derived_columns(frame.copy())

    assert list(noted["Notes"]) == ["", ""], (
        "a note stored against an incomplete identity surfaced against variants it was not "
        "written about"
    )


def test_a_note_still_reaches_the_variant_it_was_written_about(monkeypatch):
    """The paired positive, so the guard above cannot pass by emptying every note.

    Same shape, same store, one column added — and now the note lands on one row and not the
    other.
    """
    _keyed_session(monkeypatch, {"TP53:17:7577120:C>T": "hotspot, discussed at board"})
    frame = pd.DataFrame(
        {
            "Hugo_Symbol": ["TP53", "TP53"],
            "Chromosome": ["17", "17"],
            "Start_Position": [7577120, 7577120],
            "Reference_Allele": ["C", "C"],
            "Tumor_Seq_Allele2": ["T", "A"],
        }
    )

    noted = variant_table._add_derived_columns(frame.copy())

    assert list(noted["Notes"]) == ["hotspot, discussed at board", ""]


@pytest.fixture
def dialog(monkeypatch):
    """`_show_variant_dialog` outside a Streamlit run, recording the widgets it drew.

    Called through `__wrapped__` because `@st.dialog` refuses to run without a script
    context. The detail panel above the notes block is stubbed out: it renders the whole
    variant and is not what this is asking about.
    """
    drawn = {"info": [], "text_area": [], "button": []}

    class _Context:
        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

    class _St:
        session_state = None  # replaced per test by `_keyed_session`

        @staticmethod
        def markdown(*args, **kwargs):
            pass

        @staticmethod
        def caption(*args, **kwargs):
            pass

        @staticmethod
        def info(text, *args, **kwargs):
            drawn["info"].append(text)

        @staticmethod
        def text_area(label, *args, **kwargs):
            drawn["text_area"].append(label)
            return ""

        @staticmethod
        def text_input(label, *args, **kwargs):
            return ""

        @staticmethod
        def expander(*args, **kwargs):
            return _Context()

        @staticmethod
        def button(label, *args, **kwargs):
            drawn["button"].append(label)
            return False

    monkeypatch.setattr(variant_table, "st", _St, raising=False)
    monkeypatch.setattr(variant_table, "render_variant_detail_panel", lambda row: None)
    _keyed_session(monkeypatch)
    return variant_table._show_variant_dialog.__wrapped__, drawn


@pytest.mark.parametrize("holed", ["Tumor_Seq_Allele2", "Hugo_Symbol"])
def test_the_dialog_says_why_rather_than_offer_a_note_it_cannot_store(dialog, holed):
    """A note field over a keyless row is the app inviting writing it will silently lose.

    The caption directly above it has just promised the note appears on any file carrying this
    variant, and a key with a hole in it cannot keep that promise. The widget keys were a
    second reason: they are built from the key too, so every such row shared one text area.

    Both components are driven because the message used to name `Tumor_Seq_Allele2` whatever the
    actual hole was — review caught it — so a blank `Hugo_Symbol` sent the reader to a column
    that was fine.

    The annotation widgets are asserted absent alongside the note field: the guard returns
    before them, an invented column is keyed exactly as a note is, and offering to create one
    here would offer a column the user cannot then fill for this variant.
    """
    show, drawn = dialog
    row = pd.Series(
        {
            "Hugo_Symbol": "TP53",
            "Chromosome": "17",
            "Start_Position": 7577120,
            "Reference_Allele": "C",
            "Tumor_Seq_Allele2": "T",
        }
    )
    row[holed] = ""
    show(row)

    assert not drawn["text_area"], "a note field was offered for a row with no identity"
    assert not drawn["button"], "and a Save button that would have stored it somewhere wrong"
    assert len(drawn["info"]) == 1
    assert holed in drawn["info"][0], (
        f"the message must name {holed}, the component actually missing, or it sends the "
        f"reader to the wrong column: {drawn['info'][0]!r}"
    )
    assert "annotation" in drawn["info"][0], (
        "the annotation widgets are withheld too, and the message does not say so"
    )


def test_the_dialog_offers_the_note_field_when_the_row_is_identified(dialog):
    """The paired positive: the guard above must not have turned notes off for everyone."""
    show, drawn = dialog
    show(
        pd.Series(
            {
                "Hugo_Symbol": "TP53",
                "Chromosome": "17",
                "Start_Position": 7577120,
                "Reference_Allele": "C",
                "Tumor_Seq_Allele2": "T",
            }
        )
    )

    assert drawn["text_area"] == ["Notes"]
    assert not drawn["info"]


# ---------------------------------------------------------------------------
# What the detail panel does with a value that says nothing (issue #131)
# ---------------------------------------------------------------------------


def _panel_text(monkeypatch, row):
    """Every string the detail panel draws, for one variant.

    The panel writes through `st.markdown`, `st.caption` and `st.metric`, and a claim about
    what a clinician sees has to read all three: the fields these guards are about are drawn
    by different ones, and a guard that watched only `markdown` would pass over a metric that
    had silently stopped being drawn.
    """
    from components import variant_detail

    drawn = []

    class _Column:
        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

    fake = MagicMock()
    fake.columns.side_effect = lambda n, *a, **k: [_Column() for _ in range(n if isinstance(n, int) else len(n))]
    fake.markdown.side_effect = lambda text, *a, **k: drawn.append(str(text))
    fake.caption.side_effect = lambda text, *a, **k: drawn.append(str(text))
    fake.metric.side_effect = lambda label, value, *a, **k: drawn.append(f"{label}: {value}")
    monkeypatch.setattr(variant_detail, "st", fake)
    monkeypatch.setattr(variant_detail, "render_acmg_evidence", lambda *a, **k: None)

    variant_detail.render_variant_detail_panel(row)
    return drawn


def _variant(**overrides):
    row = {
        "Hugo_Symbol": "TP53",
        "Variant_Classification": "Missense_Mutation",
        "Variant_Type": "SNP",
        "Chromosome": "17",
        "Start_Position": 7577120,
        "End_Position": 7577120,
        "Reference_Allele": "C",
        "Tumor_Seq_Allele2": "T",
    }
    row.update(overrides)
    return pd.Series(row)


def test_a_protein_change_the_file_carries_is_not_lost_to_a_dead_fallback(monkeypatch):
    """`Protein_Change or AAChange.refGene` never fell through, because NaN is truthy.

    So the fallback column was unreachable whenever `Protein_Change` was present and empty,
    which is the only state it exists for. Measured over 303 real MAFs: 6,267 rows carry a
    protein change in `AAChange.refGene` that the panel did not draw.
    """
    drawn = _panel_text(
        monkeypatch,
        _variant(Protein_Change=float("nan"), **{"AAChange.refGene": "p.R273H"}),
    )

    assert any("p.R273H" in text for text in drawn), (
        "the fallback column carries a protein change and the panel drew nothing"
    )


@pytest.mark.parametrize("sentinel", [".", "UNKNOWN"])
def test_the_revived_fallback_does_not_draw_the_annotators_sentinel(monkeypatch, sentinel):
    """The other half of the same fix, and the reason the two had to land together.

    `AAChange.refGene` spells "no annotation" as `.` on 130,523 rows and `UNKNOWN` on 788 of
    the corpus. Reviving the fallback against the old check — which knew neither — would have
    put those on screen as if they were findings.
    """
    drawn = _panel_text(
        monkeypatch,
        _variant(Protein_Change=float("nan"), **{"AAChange.refGene": sentinel}),
    )

    assert not any("Protein Change" in text for text in drawn), (
        f"{sentinel!r} is the annotator saying it had nothing, and it reached the panel"
    )


def test_a_read_count_of_zero_is_a_measurement_and_is_drawn(monkeypatch):
    """`if t_alt and t_ref` discarded a real zero, and zero is common.

    `t_ref_count` is 0 on 67,797 rows across 157 of the 303 MAFs measured — every variant
    called at 100% VAF — and the Alt/Ref caption was absent on all of them.
    """
    drawn = _panel_text(monkeypatch, _variant(t_alt_count=61, t_ref_count=0))

    assert any("61/0" in text for text in drawn), (
        "a variant with no reference reads reported no read counts at all"
    )


def test_a_vaf_of_zero_is_drawn_rather_than_dropped(monkeypatch):
    """The same falsy-zero reading, one field up."""
    drawn = _panel_text(monkeypatch, _variant(tumor_f=0.0))

    assert any(text.startswith("VAF: ") for text in drawn)


@pytest.mark.parametrize(
    "value, why",
    [
        (float("nan"), "the column is present and the cell is empty"),
        ("Unknown", "ANNOVAR could not assign a gene"),
        (".", "the annotator had no value for this row"),
    ],
)
def test_a_field_with_nothing_in_it_draws_a_dash_rather_than_the_word(monkeypatch, value, why):
    """`row.get(column, "—")` only ever defended against an absent *column*.

    A column that is present and says nothing went straight through the default, so the panel
    drew the string `nan` — or the annotator's own sentinel — at a clinician, in the one field
    that names the gene.
    """
    drawn = _panel_text(monkeypatch, _variant(Hugo_Symbol=value))

    gene = next(text for text in drawn if text.startswith("**Gene:**"))
    assert gene == "**Gene:** `—`", why


def test_a_gene_the_annotator_could_not_name_still_carries_a_note(monkeypatch):
    """The line between the two predicates, asserted where crossing it would cost something.

    `Hugo_Symbol` is ANNOVAR's `Unknown` on 2,605 rows across 136 of the 303 MAFs measured. The
    panel above must not *draw* that word, and the note key must still *accept* it: the row is
    identified by its position and alleles, so a note attaches to exactly one variant. Folding
    the note key onto the display predicate would withdraw notes from nearly half of real
    files, silently, which is the failure direction issue #14 names.
    """
    _keyed_session(monkeypatch)

    key = variant_table._variant_key(
        pd.Series(
            {
                "Hugo_Symbol": "Unknown",
                "Chromosome": "17",
                "Start_Position": 7577120,
                "Reference_Allele": "C",
                "Tumor_Seq_Allele2": "T",
            }
        )
    )

    assert key == "Unknown:17:7577120:C>T"


def test_a_barcode_funcotator_could_not_determine_falls_back_to_the_normal(monkeypatch):
    """`Tumor_Sample_Barcode` is entirely `__UNKNOWN__` in 70 of the 303 MAFs measured.

    The Sample_Name fallback was the only site in the app that knew that spelling, and it is
    the one member of its old set that is load-bearing — so the sweep onto the shared predicate
    had to keep it rather than tidy it away.
    """
    _keyed_session(monkeypatch)
    frame = pd.DataFrame(
        {
            "Tumor_Sample_Barcode": ["__UNKNOWN__", "TUMOUR-1"],
            "Matched_Norm_Sample_Barcode": ["NORMAL-1", "NORMAL-2"],
            "Hugo_Symbol": ["TP53", "BRCA1"],
            "Chromosome": ["17", "17"],
            "Start_Position": [7577120, 41244000],
            "Reference_Allele": ["C", "A"],
            "Tumor_Seq_Allele2": ["T", "G"],
        }
    )

    named = variant_table._add_derived_columns(frame)

    assert list(named["Sample_Name"]) == ["NORMAL-1", "TUMOUR-1"]


# ---------------------------------------------------------------------------
# Each classification, in one place, gated on the file (issue #187)
# ---------------------------------------------------------------------------
#
# The panel drew `CancerVar` and `InterVar` in the clinical-badge row *and* again three rows
# down as the guideline verdict, so a germline file showed `InterVar: Pathogenic` twice, in two
# colours, under two headings. The verdict values below are deliberately unlike anything else
# the panel can draw, so counting them counts renders rather than coincidences.

_GERMLINE_EVIDENCE_COLUMN = " InterVar: InterVar and Evidence "
_SOMATIC_EVIDENCE_COLUMN = " CancerVar: CancerVar and Evidence "


def test_a_germline_classification_is_drawn_exactly_once(monkeypatch):
    """One value, one row, one badge — and the badge row no longer claims a verdict.

    Both halves matter. Keeping the two out of the clinical row's membership is what stops the
    duplicate, and asserting the count is what stops it coming back through the other row. That
    membership is `components.clinical_badges.CLINICAL_ROW` since issue #204 moved it there;
    `tests/test_clinical_badges.py` makes the same drawn-exactly-once claim about the five
    annotations that map #199 re-homed.
    """
    from components.clinical_badges import CLINICAL_ROW as _CLINICAL_BADGE_CONFIG

    drawn = _panel_text(
        monkeypatch,
        _variant(
            InterVar="Pathogenic_ZZQ",
            **{_GERMLINE_EVIDENCE_COLUMN: "PVS1=1"},
        ),
    )

    assert "\n".join(drawn).count("Pathogenic_ZZQ") == 1
    assert {entry[0] for entry in _CLINICAL_BADGE_CONFIG}.isdisjoint({"InterVar", "CancerVar"})


def test_a_somatic_tier_is_drawn_exactly_once(monkeypatch):
    drawn = _panel_text(
        monkeypatch,
        _variant(
            CancerVar="Tier_II_potential_ZZQ",
            **{
                _SOMATIC_EVIDENCE_COLUMN: (
                    "10#Tier_II_potential EVS=[2, 2, 0, 1, 0, 0, 0, 1, 1, 0, 1, 2]"
                )
            },
        ),
    )

    assert "\n".join(drawn).count("Tier_II_potential_ZZQ") == 1


def test_the_off_arm_classifier_is_not_named_at_all(monkeypatch):
    """A somatic file says nothing about InterVar — not even that it is missing.

    `InterVar` is on 0 of 95 real somatic MAFs, so the column's own absence is the gate, and a
    somatic panel that explained the absence of a germline classifier would be answering a
    question the file never raised.
    """
    drawn = _panel_text(
        monkeypatch,
        _variant(
            CancerVar="Tier_II_potential",
            **{_SOMATIC_EVIDENCE_COLUMN: "10#Tier_II_potential"},
        ),
    )

    assert "InterVar" not in "\n".join(drawn)


def test_a_tier_whose_evidence_column_is_absent_says_which_column_is_missing(monkeypatch):
    """4 of the 98 real files carrying `CancerVar` have no evidence vector — 4,044 rows.

    Issue #210 measured what the sentence claims rather than only that a column name appears in
    it: on those four files the only *other* columns with `Evidence` in the name are CIViC's,
    which is a different tool's opinion and not the reasoning behind this tier, so "the evidence
    CancerVar recorded behind this tier is not in it" holds. Asserted as the whole sentence,
    because the column name alone was in the drawn text on every somatic row in the corpus and
    so could not tell this state from any other.
    """
    drawn = _panel_text(monkeypatch, _variant(CancerVar="Tier_II_potential"))

    said = "\n".join(drawn)
    assert (
        "This file carries no `CancerVar: CancerVar and Evidence` column, so the evidence "
        "CancerVar recorded behind this tier is not in it." in said
    )


def test_evidence_with_no_tier_column_points_at_the_tier_drawn_below_it(monkeypatch):
    """15 of the 109 real files carrying a CancerVar evidence vector have no `CancerVar` column.

    The tier is a sum CancerVar itself computes, and map #184's standing decision is to read the
    tool's answer rather than re-derive it — so the panel points at the tier the section below
    badges from the string instead of adding the vector up.

    **A full twelve-entry vector, and that is the point.** This test used to feed
    `EVS=[2, 2, 0]` — three entries, which `parse_cancervar` refuses — and assert only that
    "`CancerVar`" and "AMP/ASCO/CAP" appeared somewhere in the panel. Both did, so it passed
    while the caption pointed at a tier the section had declined to draw: issue #210's defect,
    sitting under a guard that could not see it. The unparseable case is now its own test, and
    this one holds the state where the promise is true.
    """
    drawn = _panel_text(
        monkeypatch,
        _variant(
            **{
                _SOMATIC_EVIDENCE_COLUMN: (
                    " CancerVar: 10#Tier_II_potential EVS=[2, 2, 0, 1, 0, 0, 0, 1, 1, 0, 1, 2] "
                )
            }
        ),
    )

    said = "\n".join(drawn)
    assert (
        "This file carries no `CancerVar` column, so the AMP/ASCO/CAP tier below is the one "
        "CancerVar printed inside its own evidence string." in said
    )
    assert "shows none" not in said, "a drawable tier was described as absent"


def test_an_empty_verdict_cell_is_named_as_a_variant_level_absence(monkeypatch):
    """1 row of the 107,296 real rows carrying `CancerVar`, in 1 of 98 files.

    A *file*-level sentence would be wrong here: the column is there and the tool ran. What it
    reached no verdict on is this one variant.

    Issue #187 recorded this as *53 of 24,330 somatic and 55 of 25,957 germline rows, present in
    nearly every file*. Issue #210 re-walked all 336,170 rows of all 184 files and found **one**
    somatic row and **no** germline row — so the state is real and orders of magnitude rarer than
    the docstring said, and "present in nearly every file" was the reverse of true. The test
    stays: a state with one row is still a state whose sentence has to be right.
    """
    drawn = _panel_text(
        monkeypatch,
        _variant(CancerVar="", **{_SOMATIC_EVIDENCE_COLUMN: "10#Tier_II"}),
    )

    assert "CancerVar recorded no tier for this variant." in "\n".join(drawn)


def test_a_germline_file_with_no_intervar_column_names_the_marker_not_the_arm(monkeypatch):
    """2 of 61 real germline files carry `RENOVO_Class` and no `InterVar` at all — 6,597 rows.

    The sentence used to open *"This file was annotated for germline analysis"*, and issue #210
    measured what sits under it: on **6,524 of those 6,597 rows** the same file carries CancerVar's
    evidence vector, so a somatic AMP/ASCO/CAP tier badge draws one row below and the page
    contradicts its own caption. `read_arm_evidence` calls the file germline because it has
    `RENOVO_Class` and no `CancerVar` *column* — not because CancerVar did not run on it.

    So the sentence quotes the marker, which is a fact about the row, instead of the arm, which is
    an inference from it. `classifier.scale` and `classifier.column` still name the classification
    that is missing, which was the arm word's stated job.
    """
    drawn = _panel_text(monkeypatch, _variant(RENOVO_Class="Pathogenic"))

    said = "\n".join(drawn)
    assert (
        "This file carries `RENOVO_Class` but no `InterVar` column, so MAFigate has no ACMG/AMP "
        "classification to show for this variant." in said
    )
    assert "annotated for" not in said, "the arm word came back"


def test_the_marker_sentence_cannot_be_reached_on_the_somatic_arm():
    """Why the somatic form of that sentence needs no wording: it cannot fire.

    The branch runs only when the classifier's verdict column is **absent**, and
    `read_arm_evidence` answers `"somatic"` only when a somatic marker is **present**. Those
    markers are exactly `("CancerVar",)` — the column itself — so on the somatic arm the two
    conditions are contradictory and there is no such sentence to be right or wrong about. Over
    336,170 real rows the marker sentence fired 6,597 times and named `InterVar` every time.

    This is the invariant rather than a sample: add a somatic marker that is *not* the verdict
    column and the branch becomes reachable, at which point someone has to decide what it says
    about a file annotated by CancerVar that carries no tier. This test is where that decision
    gets asked for.
    """
    from components.variant_detail import _CLASSIFIERS
    from filters.arm_detection import ANNOTATOR_MARKERS

    somatic = next(c for c in _CLASSIFIERS if c.arm == "somatic")

    assert ANNOTATOR_MARKERS["somatic"] == (somatic.column,), (
        "the somatic arm can now be detected without the CancerVar column, so the "
        "marker sentence is reachable on it and needs wording measured against real files"
    )


def test_the_absence_is_named_even_when_other_badges_draw(monkeypatch):
    """The gap this ticket closes, exactly.

    ClinVar is on 168 of 177 real MAFs, so the guideline row nearly always has *something* to
    draw — and the old caption fired only when it had nothing. A missing tier beside a drawn
    ClinVar badge was therefore silent, which reads as "this variant has no classification".
    """
    drawn = _panel_text(
        monkeypatch,
        _variant(ClinVar_VCF_CLNSIG="Pathogenic", CancerVar=""),
    )

    said = "\n".join(drawn)
    assert "ClinVar" in said
    assert "CancerVar recorded no tier for this variant." in said


@pytest.mark.parametrize(
    "tier, expected, why",
    [
        ("Tier_IV_benign", "#2ca02c", "19,740 real cells — the commonest tier by far"),
        ("Tier_III_Uncertain", "#f0c420", "4,140 real cells"),
        ("Tier_II_potential", "#ff7f0e", "379 real cells"),
        ("Tier_I_strong", "#d62728", "11 real cells — the only tier that earns the red"),
        ("1#Tier_I_strong Evidence", "#d62728", "7 real cells carry the whole evidence string"),
        ("Tier_V_invented", "#7f7f7f", "outside CancerVar's vocabulary: unknown, not benign"),
    ],
)
def test_each_tier_draws_in_its_own_colour(monkeypatch, tier, expected, why):
    """`"Tier_I" in val` matched all four tier names, so every badge drew Tier I's red.

    The three `elif` branches beneath it were unreachable, which is why this went unseen: the
    colour was wrong on 99.9% of somatic variants and no value ever reached the code that would
    have coloured it correctly.
    """
    drawn = "\n".join(_panel_text(monkeypatch, _variant(CancerVar=tier)))

    assert f"background-color:{expected}" in drawn
    assert tier in drawn


def test_the_tier_colours_are_ordered_so_no_branch_is_unreachable():
    """The table itself, checked without rendering: a longer name must never sit after a
    shorter one it contains, or the shorter one shadows it exactly as before.

    The table moved to ``components.cbp_evidence`` in issue #189 — CancerVar's tier names are
    CancerVar's vocabulary, and one definition is what keeps the guideline-row badge and the CBP
    section's badge from drawing one tier in two colours. The rule it guards did not move.
    """
    from components.cbp_evidence import _TIER_COLORS

    names = [name for name, _, _ in _TIER_COLORS]
    for position, name in enumerate(names):
        shadowed = [later for later in names[position + 1:] if name in later]
        assert not shadowed, f"{name} at {position} shadows {shadowed}"


def test_the_panel_never_reads_the_arm_the_user_selected():
    """The gate is the file's columns, never `filter_params["sample_type"]` (issue #187).

    Gating on the selected arm would hide a genuinely present InterVar verdict from anyone
    sitting on the somatic arm — the wrong-arm session issue #135 exists for, whose settled
    rule is *detect and warn, never override*.

    Asserted over the parsed module rather than its text, because this module's own comments
    name the parameter in order to say it is not read; a text guard would be satisfied by the
    sentence explaining the rule it is meant to enforce.
    """
    import ast
    import inspect

    from components import variant_detail

    tree = ast.parse(inspect.getsource(variant_detail))
    constants = {
        node.value
        for node in ast.walk(tree)
        if isinstance(node, ast.Constant) and isinstance(node.value, str)
    }
    named = {
        node.attr for node in ast.walk(tree) if isinstance(node, ast.Attribute)
    } | {node.id for node in ast.walk(tree) if isinstance(node, ast.Name)}

    assert "sample_type" not in constants
    assert "filter_params" not in constants
    assert "current_arm" not in named


# --- What the crash handler can still hand back (issue #144) -------------------------
#
# `MAFigate.py`'s outermost handler used to end with *Please refresh the page or contact
# support*. Refreshing starts the session over, and a note lives in session state and
# nowhere else (issue #67) — so the app's advice of last resort was *discard your writing*,
# offered at the moment a user is most likely to take it. `written_notes` is what it offers
# instead, and every property below is a way that offer could quietly become worthless.


@pytest.fixture
def written(monkeypatch):
    """Two variants written about, and a report that holds only one of them.

    The report is seeded deliberately even though `written_notes` must never read it: it is
    what makes `test_the_rescue_reaches_writing_the_report_does_not_hold` a real claim rather
    than a restatement of the implementation.
    """
    session = FakeSessionState()
    monkeypatch.setattr(variant_table.st, "session_state", session, raising=False)

    passed = pd.DataFrame(
        {
            "Hugo_Symbol": ["TP53"],
            "Chromosome": ["17"],
            "Start_Position": [7577120],
            "Reference_Allele": ["C"],
            "Tumor_Seq_Allele2": ["T"],
        }
    )
    rejected = pd.DataFrame(
        {
            "Hugo_Symbol": ["KRAS"],
            "Chromosome": ["12"],
            "Start_Position": [25398284],
            "Reference_Allele": ["C"],
            "Tumor_Seq_Allele2": ["A"],
        }
    )
    kept_key = variant_table._variant_key(passed.iloc[0])
    dropped_key = variant_table._variant_key(rejected.iloc[0])

    session["filtered_data"] = passed
    session["variant_notes"] = {
        kept_key: "hotspot, discussed at board",
        dropped_key: "expected this one to pass",
    }
    session["custom_annotations"] = {"Actionable": {kept_key: "yes"}}
    return session


def test_the_rescue_names_the_columns_the_users_own_maf_names(written):
    """A rescue file the user cannot re-attach to their MAF is barely a rescue.

    The store holds `gene:chrom:pos:ref>alt` and nothing else, so these five columns are
    reversed back out of the key. They are the names in the user's own header row, which is
    the map's Style rule as issue #79 bounded it — the file's names, not MAFigate's.
    """
    rescued = variant_table.written_notes()

    assert list(rescued.columns) == [
        "Hugo_Symbol",
        "Chromosome",
        "Start_Position",
        "Reference_Allele",
        "Tumor_Seq_Allele2",
        "Notes",
        "Actionable",
    ]
    kept = rescued[rescued["Hugo_Symbol"] == "TP53"].iloc[0]
    assert list(kept[list(variant_table._KEY_COLUMNS)]) == ["TP53", "17", "7577120", "C", "T"]
    assert kept["Notes"] == "hotspot, discussed at board"
    assert kept["Actionable"] == "yes"


def test_the_rescue_reaches_writing_the_report_does_not_hold(written):
    """The reason the payload is the writing rather than the report.

    `filtered_data` holds only what **passed**, and both tabs draw the note dialog — so a
    rescue built on it would drop every note written on a rejected variant while reporting
    success, which is the one silent failure `with_user_columns` exists to prevent. Fails the
    moment anyone rewrites this to serialise a frame.
    """
    rescued = variant_table.written_notes()

    assert sorted(rescued["Hugo_Symbol"]) == ["KRAS", "TP53"]
    dropped = rescued[rescued["Hugo_Symbol"] == "KRAS"].iloc[0]
    assert dropped["Notes"] == "expected this one to pass"
    # And the annotation column the other variant carries is present but blank here, rather
    # than absent — the file has one shape whatever any given variant was written about.
    assert dropped["Actionable"] == ""


def test_the_rescue_survives_having_no_report_at_all(written):
    """The state that dissolved the handler's third branch.

    A file opened and not yet filtered leaves `filtered_data` at `None` while the notes live
    on — issue #67 reversed #64's clearing rule — so a report-shaped rescue would have had
    nothing to offer at exactly the moment there was something to lose. There is therefore no
    "written something, nothing to save" state, and the handler asks one question, not two.
    """
    written["filtered_data"] = None

    rescued = variant_table.written_notes()

    assert not rescued.empty
    assert sorted(rescued["Hugo_Symbol"]) == ["KRAS", "TP53"]


def test_the_rescue_is_empty_but_still_shaped_when_nothing_is_written(monkeypatch):
    """Empty is the gate: it is what tells the handler not to alarm a user with nothing to lose.

    Keeps its columns while empty so a caller reads the shape off it rather than restating it.
    """
    session = FakeSessionState()
    monkeypatch.setattr(variant_table.st, "session_state", session, raising=False)

    rescued = variant_table.written_notes()

    assert rescued.empty
    assert list(rescued.columns)[:5] == list(variant_table._KEY_COLUMNS)
    assert list(rescued.columns)[5:] == ["Notes"]


def test_a_variant_written_about_and_then_emptied_is_not_listed(written):
    """A rescue file naming variants the user wrote nothing about misreports what it holds.

    The save path already pops an emptied note rather than storing `""`, so this is defensive
    — but the stores are plain dicts that anything could put a blank into.
    """
    written["variant_notes"] = {key: "  " for key in written["variant_notes"]}
    written["custom_annotations"] = {"Actionable": {}}

    assert variant_table.written_notes().empty


@pytest.mark.parametrize("allele", ["<DEL>", "<DUP>"])
def test_a_symbolic_allele_survives_the_round_trip(monkeypatch, allele):
    """The measured case, and the reason the split runs forwards.

    Across 340,972 rows of the 198 real MAFs on this machine carrying all five key columns,
    only `Tumor_Seq_Allele2` ever contains `>` — on two rows, where the allele is symbolic.
    `rsplit(">", 1)` reads `C><DEL>` as `C><DEL` and `""`, losing both alleles at once; the
    separator that can repeat is on the side that never carries it. This fails against that
    reading and passes against `partition`.
    """
    session = FakeSessionState()
    monkeypatch.setattr(variant_table.st, "session_state", session, raising=False)

    row = pd.Series(
        {
            "Hugo_Symbol": "TP53",
            "Chromosome": "17",
            "Start_Position": "7577120",
            "Reference_Allele": "C",
            "Tumor_Seq_Allele2": allele,
        }
    )
    session["variant_notes"] = {variant_table._variant_key(row): "structural"}

    rescued = variant_table.written_notes()

    assert list(rescued.iloc[0][list(variant_table._KEY_COLUMNS)]) == [
        "TP53",
        "17",
        "7577120",
        "C",
        allele,
    ]


def test_taking_a_key_apart_cannot_raise_inside_the_handler():
    """It runs where a second exception would replace the app's account of the first.

    Unreachable through `_variant_key`, which writes exactly three colons and refuses a key
    with a blank component — so this asserts the padding rather than a case anyone has met.
    """
    assert variant_table._components_of_key("TP53") == ["TP53", "", "", "", ""]
    assert variant_table._components_of_key("TP53:17") == ["TP53", "17", "", "", ""]

# ---------------------------------------------------------------------------
# Why this variant is not in the report (issue #147)
# ---------------------------------------------------------------------------


def _standing_text(monkeypatch, row, session):
    """Every string the report-standing block draws, for one variant.

    Reads `markdown`, `caption` and `warning`, because the block writes its heading through
    one, its actionable line through the second and the filled-column note through the third
    — and a guard watching only `markdown` would pass over either of the other two silently
    disappearing.
    """
    from components import variant_detail

    drawn = []

    fake = MagicMock()
    fake.session_state = session
    fake.markdown.side_effect = lambda text, *a, **k: drawn.append(str(text))
    fake.caption.side_effect = lambda text, *a, **k: drawn.append(str(text))
    fake.warning.side_effect = lambda text, *a, **k: drawn.append(str(text))
    monkeypatch.setattr(variant_detail, "st", fake)

    variant_detail._render_report_standing(row)
    return drawn


class _Run:
    """The snapshot `filter_run` carries, as much of it as this block reads."""

    def __init__(self, params):
        self.params = params


def _rejected_setup():
    """A one-variant MAF, the parameters that reject it, and the displayed row.

    The displayed row is built the way `create_data_table` builds it — `fillna("")` then
    `astype(str)` on object columns — so these guards run against the frame the dialog is
    really handed rather than a tidy copy of the MAF.
    """
    from tests.test_attribution import a_frame, params_for

    maf = a_frame("somatic", t_alt_count=1, t_ref_count=1, tumor_f=0.001)
    params = params_for("somatic", skip_pathogenic=True)

    displayed = maf.fillna("")
    for column in displayed.columns:
        if displayed[column].dtype == "object":
            displayed[column] = displayed[column].astype(str)
    return maf, params, displayed.iloc[0]


def test_the_panel_names_every_setting_a_variant_falls_outside(monkeypatch):
    """The 70-95% case reaches the screen with all of its settings, not one of them.

    Mutation: draw only `explanation.failing[0]` -> the variant is told about one setting
    while a second is equally responsible, which is the misleading sentence this ticket
    exists to avoid.
    """
    from filters.attribution import label_for

    maf, params, row = _rejected_setup()
    drawn = _standing_text(
        monkeypatch, row, {"maf_data": maf, "filter_run": _Run(params)}
    )
    text = "\n".join(drawn)

    assert "Not in this report" in text
    assert "2 of your settings" in text
    assert label_for("depth", "somatic") in text
    assert label_for("vaf", "somatic") in text
    assert "No single change brings it back" in text


def test_the_panel_shows_what_the_file_says_beside_each_setting(monkeypatch):
    """The row's own values, so the user can see why it fell outside.

    Mutation: stop drawing `outcome.values` -> the block names settings and gives the reader
    nothing to check them against.
    """
    maf, params, row = _rejected_setup()
    text = "\n".join(
        _standing_text(monkeypatch, row, {"maf_data": maf, "filter_run": _Run(params)})
    )

    assert "t_alt_count" in text and "t_ref_count" in text
    assert "tumor_f" in text


def test_the_panel_says_nothing_when_the_row_is_not_the_one_on_screen(monkeypatch):
    """The wrong-row guard, and the reason it is not paranoia.

    The AgGrid path recovers the index with
    `pd.DataFrame(grid_response.selected_rows).iloc[0].name`, and that frame is built from a
    list of dicts returned by the browser — so its index is a RangeIndex position, which can
    exist in the report's index and resolve to a different variant. An explanation of a
    variant the user did not select is worse than none, so a mismatch draws nothing.

    Mutation: drop the identity comparison from `_row_as_the_filter_saw_it` -> this variant
    is confidently told why a different one is missing.
    """
    maf, params, row = _rejected_setup()
    impostor = row.copy()
    impostor["Hugo_Symbol"] = "EGFR"

    drawn = _standing_text(
        monkeypatch, impostor, {"maf_data": maf, "filter_run": _Run(params)}
    )
    assert drawn == [], drawn


def test_the_panel_says_nothing_with_no_run_to_describe(monkeypatch):
    """No filter run, no settings to be outside of.

    Mutation: fall back to `st.session_state.filter_params` -> the block describes the
    controls as they stand for a report that does not exist, which is the defect #137 spent
    a ticket removing from the downloaded report.

    The third case is what catches that mutation, and the first draft of this test lacked it:
    `filter_params` is seeded by `MAFigate.py` on every session, so a session without it does
    not occur and a fallback reading it survived a guard built from the first two alone.
    """
    maf, params, row = _rejected_setup()
    assert _standing_text(monkeypatch, row, {"maf_data": maf}) == []
    assert _standing_text(monkeypatch, row, {"maf_data": maf, "filter_run": None}) == []
    assert _standing_text(
        monkeypatch, row, {"maf_data": maf, "filter_params": params}
    ) == []


def test_the_panel_tells_a_passing_variant_why_it_is_in(monkeypatch):
    """The dialog is shared by both tabs, so the block has to be true of a passing row.

    Mutation: return early for anything in the report -> a user opening a variant from the
    Passed tab meets a block that only ever speaks about absences.
    """
    from tests.test_attribution import a_frame, params_for

    maf = a_frame("somatic")
    params = params_for("somatic")
    row = maf.iloc[0]

    text = "\n".join(
        _standing_text(monkeypatch, row, {"maf_data": maf, "filter_run": _Run(params)})
    )
    assert "In this report" in text
    assert "met your filter criteria" in text



if __name__ == "__main__":
    unittest.main()
