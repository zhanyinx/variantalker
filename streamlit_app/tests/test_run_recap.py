"""What the filter recap claims, and what stops it claiming something false (issue #137).

The results view said *"passed all applied filters"* three times and named none of them.
The recap names them — which puts it in the company of the eight on-screen parameter echoes
issue #28 deleted, every one of which was **wrong**: they read a catch-all sentinel that no
longer meant anything, the superseded germline VAF key, and the classification list as
though it were an include list.

So the tests here are not "does it render". They are the three ways that class of defect
gets in:

* it names a setting in words the interface does not use anywhere else;
* it reads a parameter by a key the filter does not read it by;
* it describes the settings as they are *now* rather than the run on screen.

Every assertion is mutation-checked and the mutation is named, because a green test over a
wired-up page proves nothing on its own (issues #67, #83, #90).
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

STREAMLIT_APP = Path(__file__).resolve().parent.parent
if str(STREAMLIT_APP) not in sys.path:
    sys.path.insert(0, str(STREAMLIT_APP))

from config.param_labels import PARAM_LABELS, label_of, labels_for  # noqa: E402
from config.pipeline_params import PIPELINE_PARAMS, pipeline_params  # noqa: E402
from config.presets import (  # noqa: E402
    CLINICAL_GERMLINE_PARAMS,
    CLINICAL_SOMATIC_PARAMS,
    SOFT_GERMLINE_PARAMS,
    SOFT_SOMATIC_PARAMS,
)
from config.vocabularies import GUIDELINE_SOURCES  # noqa: E402
from filters.absent_columns import CIVIC_GUARD_COLUMN  # noqa: E402
from filters.run_recap import (  # noqa: E402
    NO_CIVIC,
    NO_FREQUENCY_COLUMNS,
    NO_RESTRICTION,
    describe_run,
)
from filters.variant_filters import (  # noqa: E402
    FREQUENCY_COLUMNS,
    Diagnostics,
    _Settings,
)

ARMS = ("somatic", "germline")

#: A file carrying every column any clause can be dropped for. The default for tests about
#: *settings*, so that a test about a label cannot pass or fail for reasons about columns.
EVERY_COLUMN = (CIVIC_GUARD_COLUMN,) + tuple(FREQUENCY_COLUMNS)

#: Every parameter set the app can actually be filtering with: the two contracts it opens
#: on and the four presets. Parametrising over these is what makes the suite fail when a
#: preset carries a spelling the recap cannot read — which is not hypothetical, the germline
#: presets carry ``vaf_threshold_germline`` with no ``vaf_threshold`` at all.
LIVE_PARAM_SETS = {
    "somatic contract": PIPELINE_PARAMS["somatic"],
    "germline contract": PIPELINE_PARAMS["germline"],
    "Broad Somatic": SOFT_SOMATIC_PARAMS,
    "Stringent Somatic": CLINICAL_SOMATIC_PARAMS,
    "Broad Germline": SOFT_GERMLINE_PARAMS,
    "Stringent Germline": CLINICAL_GERMLINE_PARAMS,
}


def _values(recap) -> dict:
    return {line.label: line.value for line in recap.lines}


class TestItNamesSettingsAsTheInterfaceDoes:
    """The Style rule: no internal parameter names reach the interface."""

    @pytest.mark.parametrize("name", sorted(LIVE_PARAM_SETS))
    def test_no_line_renders_a_parameter_key(self, name):
        """``filter_cancervar: ["Tier_I_strong"]`` is precisely what must not appear.

        Checked against the *keys themselves*, so it cannot be satisfied by a label that
        merely looks tidy. The keys are drawn from the parameter sets under test rather
        than listed here, so a key added to the contract is covered without this test being
        edited.

        Mutation: label a row ``"filter_cancervar"`` → fails on that row.
        """
        params = LIVE_PARAM_SETS[name]
        recap = describe_run(params, EVERY_COLUMN)
        rendered = " ".join(f"{line.label} {line.value}" for line in recap.lines)

        for key in params:
            assert key not in rendered, f"{key!r} reached the interface: {rendered}"

    def test_every_row_in_the_table_can_be_spelled(self):
        """``spell`` raises on a ``kind`` it does not know, and nothing may reach that.

        A row added to ``PARAM_LABELS`` with a new ``kind`` and no spelling for it would
        raise inside the results view, where the page's own ``except`` turns it into
        *Error displaying results* — the report replaced by a stack trace's worth of
        apology because a label was added.

        Mutation: give any row ``kind="rate"`` → fails here, naming the row.
        """
        from filters.run_recap import spell

        for arm in ARMS:
            settings = _Settings.from_params(pipeline_params(arm))
            for row in labels_for(arm):
                assert spell(row, settings, EVERY_COLUMN), row.id

    @pytest.mark.parametrize("arm", ARMS)
    def test_every_guideline_control_is_named_by_the_shared_table(self, arm):
        """The parameter page's controls and the recap read one set of labels.

        Not "the labels happen to match" — the page's ``_GUIDELINE_CONTROLS`` is *built
        from* :data:`config.param_labels.PARAM_LABELS`, and this asserts that it still is.
        Two hand-written sets of names is how the copy beside the control stays right and
        the copy elsewhere drifts, which issue #79 measured nineteen times over.

        Mutation: write a literal label back into ``_GUIDELINE_CONTROLS`` → fails, naming
        the control whose label stopped coming from the table.
        """
        from page_modules.parameter_config import _GUIDELINE_CONTROLS

        known = {row.label for row in PARAM_LABELS}
        for key, label, _options, _help in _GUIDELINE_CONTROLS[arm]:
            assert label in known, f"{key}'s label {label!r} is not in PARAM_LABELS"

    @pytest.mark.parametrize("arm", ARMS)
    def test_every_guideline_source_this_arm_filters_on_has_a_line(self, arm):
        """A source that cuts the report and is not named is the gap this surface exists for.

        Derived from ``GUIDELINE_SOURCES``, the same constant the parameter page's
        empty-selection warning reads, rather than from a list written here — a hand-kept
        list is the failure issue #127 found in the Required Columns tab.

        Mutation: drop the ESCAT row from ``PARAM_LABELS`` → fails on the somatic arm.
        """
        recap = describe_run(pipeline_params(arm), EVERY_COLUMN)
        labels = set(_values(recap))

        for key in GUIDELINE_SOURCES[arm]:
            # `filter_clinvar` is drawn by two controls, so it is two rows; the rest are one.
            rows = [row for row in labels_for(arm) if key.replace("filter_", "") in row.id]
            assert rows, f"{key} has no row on the {arm} arm"
            for row in rows:
                assert row.label in labels

    @pytest.mark.parametrize("arm", ARMS)
    def test_it_names_no_source_this_arm_does_not_filter_on(self, arm):
        """ESCAT on the germline arm is the app's single largest measured divergence.

        ``germline_filters()`` takes no ESCAT argument, so a germline recap naming ESCAT
        would describe a clause that does not exist — and this app has shipped exactly that
        mistake before, in the widget the germline controls used to draw.

        Mutation: give the ESCAT row both arms → fails on germline.
        """
        other = "germline" if arm == "somatic" else "somatic"
        recap = describe_run(pipeline_params(arm), EVERY_COLUMN)
        labels = set(_values(recap))

        for key in GUIDELINE_SOURCES[other]:
            if key in GUIDELINE_SOURCES[arm]:
                continue
            source = key.replace("filter_", "")
            assert not any(
                source.lower() in label.lower() for label in labels
            ), f"the {arm} recap names {source}, which that arm does not filter on"


class TestItReadsWhatTheFilterRead:
    """The recap's values come off the object the filter ran on, not off the dict."""

    @pytest.mark.parametrize("name", sorted(LIVE_PARAM_SETS))
    def test_every_value_is_the_resolved_one(self, name):
        """No line may say ``None``, and the threshold lines must carry the real number.

        The germline VAF is the case that makes this more than a smoke test: the widget
        writes ``vaf_threshold_germline``, the contract writes ``vaf_threshold``, and both
        germline presets carry only the former. A recap reading the dict by key renders
        ``None`` for it on three of the six sets under test.

        Mutation: read ``params.get("vaf_threshold")`` instead of the resolved settings →
        fails on both germline presets.
        """
        params = LIVE_PARAM_SETS[name]
        recap = describe_run(params, EVERY_COLUMN)
        settings = _Settings.from_params(params)
        values = _values(recap)

        for line in recap.lines:
            assert "None" not in line.value, f"{line.label} rendered None"

        vaf = [v for label, v in values.items() if label.startswith("VAF Threshold")]
        assert len(vaf) == 1, values
        assert f"{settings.vaf_threshold:.4g}" in vaf[0]

    @pytest.mark.parametrize("arm", ARMS)
    def test_pathogenic_retention_follows_the_filter_not_the_key(self, arm):
        """The control says *retain*, the parameter says *skip*, and the presets say neither.

        ``keep_pathogenic`` is the older spelling and the SOFT/CLINICAL presets still carry
        it; the filter prefers ``skip_pathogenic`` wherever a dict has both. A recap reading
        the key it likes best is how the app once shipped a checkbox with no effect.

        Mutation: render ``params.get("skip_pathogenic")`` directly → the ``keep_pathogenic``
        case reads as off while the filter has it on.
        """
        params = dict(pipeline_params(arm))
        params.pop("skip_pathogenic", None)
        params["keep_pathogenic"] = True

        recap = describe_run(params, EVERY_COLUMN)
        assert _values(recap)[label_of("skip_pathogenic")] == "on"

        params["keep_pathogenic"] = False
        assert _values(describe_run(params, EVERY_COLUMN))[label_of("skip_pathogenic")] == "off"

    def test_an_empty_source_says_it_places_no_restriction(self):
        """The state the surface exists to catch, in the words the control's caption uses.

        An empty selection is not a filter set to nothing; it drops out of the OR. Saying
        "none" would read as *this source rejected everything*, which is the opposite.

        Mutation: render an empty list as ``""`` → the line reads as a blank setting.
        """
        params = dict(PIPELINE_PARAMS["somatic"])
        params["filter_civic"] = []

        values = _values(describe_run(params, EVERY_COLUMN))
        assert values[label_of("civic_keep")] == NO_RESTRICTION

    def test_an_unset_gene_panel_says_so_rather_than_going_unmentioned(self):
        """The dev's own example: a gene panel you believe you set, and did not.

        Listing only what fired would drop this line entirely, which is the one omission
        this surface cannot afford.

        Mutation: skip rows whose setting cuts nothing → this line disappears.
        """
        params = dict(PIPELINE_PARAMS["somatic"])
        params["filter_genes"] = []

        values = _values(describe_run(params, EVERY_COLUMN))
        assert "every gene" in values[label_of("gene_selection")]

        params["filter_genes"] = ["TP53", "BRCA1", "EGFR"]
        assert "3 genes" in _values(describe_run(params, EVERY_COLUMN))[label_of("gene_selection")]

    def test_the_frequency_filter_says_when_it_is_off(self):
        """1.0 is not "very permissive"; the filter is short-circuited and does not run.

        Mutation: render ``1.0`` as a number → the recap reports a threshold for a filter
        that made no comparison.
        """
        params = dict(PIPELINE_PARAMS["somatic"])
        params["max_freq_population"] = 1.0
        assert "off" in _values(describe_run(params, EVERY_COLUMN))[label_of("max_freq_population")]

        params["max_freq_population"] = 0.01
        assert "0.01" in _values(describe_run(params, EVERY_COLUMN))[label_of("max_freq_population")]

    def test_the_thresholds_state_their_own_boundary(self):
        """Depth is ``>=`` in the vendored clause and VAF is ``>``; the recap says which.

        A variant sitting exactly on the VAF threshold is dropped, so "or above" would be
        wrong by one boundary — the kind of quiet inaccuracy that makes a user trust a
        recap and mis-read their report.

        Mutation: spell VAF as "or above" → fails.
        """
        values = _values(describe_run(PIPELINE_PARAMS["somatic"], EVERY_COLUMN))
        assert values[label_of("min_depth")] == "50 reads or more"
        assert values[label_of("vaf_threshold")] == "above 0.01"


class TestItWillNotClaimAClauseThatDidNotRun:
    """A setting the user made is not the same fact as a restriction that applied.

    Two clauses can be dropped whole by a column the MAF does not carry, and neither leaves
    a trace anywhere else: nothing is filled, so ``filled_columns`` is empty and
    ``filled_input_note`` says nothing. They are the one way a recap built from a correct
    snapshot of correct parameters can still be false about the report beside it.
    """

    def test_it_says_when_the_civic_clause_did_not_run(self):
        """``CIViC_Evidence_Level`` is absent in 100 of 100 reference files.

        So this is not an edge case — it is what the recap says about nearly every MAF this
        app has been pointed at. ``somatic_filters`` opens ``if skip_civic or not
        civic_column_exists`` and drops CIViC from the criteria OR *and* the
        pathogenic-retention OR, so the user's selection cut nothing.

        Mutation: render the CIViC row like any other keep-list → the recap states
        ``A, B, C`` as an applied restriction on a file that was never cut on CIViC.
        """
        params = dict(PIPELINE_PARAMS["somatic"])
        without_civic = tuple(FREQUENCY_COLUMNS)

        value = _values(describe_run(params, without_civic))[label_of("civic_keep")]
        assert NO_CIVIC in value
        # The selection is still printed: the recap answers what was asked for as well as
        # what happened, and a user checking their memory needs to see both.
        assert "A, B, C" in value

        applied = _values(describe_run(params, EVERY_COLUMN))[label_of("civic_keep")]
        assert NO_CIVIC not in applied

    def test_skip_civic_reads_as_not_applied_even_where_the_column_is_present(self):
        """No widget sets ``skip_civic``, but an uploaded parameter file can.

        ``adopt_parameters`` replaces ``filter_params`` wholesale, so a user can arrive on
        this arm with ``skip_civic: True`` and no control anywhere showing it. That is
        precisely why the recap must derive this from the run rather than from the controls.
        """
        params = dict(PIPELINE_PARAMS["somatic"])
        params["skip_civic"] = True

        value = _values(describe_run(params, EVERY_COLUMN))[label_of("civic_keep")]
        assert "not applied" in value

    def test_it_says_when_the_frequency_filter_had_nothing_to_judge_by(self):
        """The app's own layer skips where no frequency column is present.

        ``frequency_mask`` returns ``None`` there and the filter appends a note at run
        time — which is drawn once and gone. The recap is what still says it later.

        Mutation: drop the column check → a threshold is reported as applied over a file no
        frequency column speaks for.
        """
        params = dict(PIPELINE_PARAMS["somatic"])
        params["max_freq_population"] = 0.01

        value = _values(describe_run(params, (CIVIC_GUARD_COLUMN,)))[
            label_of("max_freq_population")
        ]
        assert NO_FREQUENCY_COLUMNS in value

        applied = _values(describe_run(params, EVERY_COLUMN))[
            label_of("max_freq_population")
        ]
        assert NO_FREQUENCY_COLUMNS not in applied

    def test_a_frequency_filter_that_is_off_is_not_called_unapplied(self):
        """At 1.0 the filter is off because the user turned it off, not because of the file.

        Two different sentences, and saying the column one over a setting the user chose
        would send them looking for a problem with their MAF.
        """
        params = dict(PIPELINE_PARAMS["somatic"])
        params["max_freq_population"] = 1.0

        value = _values(describe_run(params, (CIVIC_GUARD_COLUMN,)))[
            label_of("max_freq_population")
        ]
        assert NO_FREQUENCY_COLUMNS not in value
        assert "off" in value

    def test_the_guarded_column_is_derived_from_the_vendored_source(self):
        """Not written down here, and not written down in ``absent_columns`` either.

        Mutation: rename the column in ``vendor/pipeline_filters.py`` → this follows it,
        where a hand-written constant would go quietly stale. (The drift guard would also
        fire, which is the point: the two agree by construction rather than by upkeep.)
        """
        assert CIVIC_GUARD_COLUMN == "CIViC_Evidence_Level"

        from filters.absent_columns import _derive_civic_guard_column

        with pytest.raises(RuntimeError):
            _derive_civic_guard_column("def check_civic_column_exists(maf):\n    return True\n")


class TestItDescribesTheRunAndNotTheControls:
    """A snapshot, because the parameters can move away from the report on screen."""

    def test_the_snapshot_does_not_follow_the_dict_it_was_taken_from(self):
        """Every widget on the parameter page writes into ``filter_params`` in place.

        So a shallow copy would leave the lists shared and the recap would silently follow
        the controls it exists to be independent of — the exact state ``data_params_hash``
        was added to detect.

        Mutation: ``dict(params)`` instead of ``deepcopy`` → the mutated list reaches the
        stored snapshot and this fails.
        """
        params = {
            "sample_type": "somatic",
            "min_depth": 50,
            "vaf_threshold": 0.01,
            "filter_cancervar": ["Tier_I_strong"],
        }
        recap = describe_run(params, EVERY_COLUMN)

        params["filter_cancervar"].append("Tier_III_Uncertain")
        params["min_depth"] = 5

        assert recap.params["filter_cancervar"] == ["Tier_I_strong"]
        assert recap.params["min_depth"] == 50
        assert _values(recap)[label_of("min_depth")] == "50 reads or more"

    def test_it_carries_the_filled_columns_the_diagnostics_reported(self):
        """The only account of a filled column that outlives the render it was made in.

        Mutation: drop ``diagnostics`` from the call → the note is ``None`` and a report
        built on stand-in values says so nowhere once its banner has been drawn.
        """
        diagnostics = Diagnostics(
            rows=10,
            criteria_only=4,
            rescue_only=1,
            both=2,
            rejected=3,
            filled_columns=("CancerVar",),
            degraded_columns=("CancerVar",),
        )
        recap = describe_run(PIPELINE_PARAMS["somatic"], EVERY_COLUMN, diagnostics)

        assert recap.note is not None
        assert "CancerVar" in recap.note
        assert "not a complete result" in recap.note

    def test_a_complete_file_gets_no_note(self):
        """The ordinary case, and the one the escalated warning almost always is not.

        Issue #136 measured the escalated fill warning firing on the file's own arm for 2
        of 173 placeable MAFs, so a note drawn unconditionally would be wrong far more
        often than right.
        """
        diagnostics = Diagnostics(
            rows=10, criteria_only=4, rescue_only=1, both=2, rejected=3
        )
        assert describe_run(PIPELINE_PARAMS["somatic"], EVERY_COLUMN, diagnostics).note is None


REFERENCE_MAF = STREAMLIT_APP / "tests" / "fixtures" / "parity" / "somatic_reference.maf"

#: The results section over a real MAF, with the parameters left where the filter run put
#: them or moved on afterwards. ``{after}`` is what happens *between* the run and the
#: render, which is the seam this whole surface is about.
_PAGE = """
import os, sys
sys.path.insert(0, {app!r})

import streamlit as st
from config.pipeline_params import pipeline_params
from page_modules import data_loading
from utils import read_maf

st.session_state.filter_params = pipeline_params("somatic")
st.session_state.maf_source_name = os.path.basename({maf!r})
st.session_state.maf_data = read_maf({maf!r})
data_loading.apply_filters_to_data(show_messages=False)
st.session_state["data_page_section"] = "results"

{after}

data_loading.show_data_loading_page()
"""


def _run_page(after: str = ""):
    pytest.importorskip("streamlit", reason="streamlit not installed")
    from streamlit.testing.v1 import AppTest

    app = AppTest.from_string(
        _PAGE.format(app=str(STREAMLIT_APP), maf=str(REFERENCE_MAF), after=after)
    )
    app.run(timeout=300)
    assert not app.exception, [str(e.value) for e in app.exception]
    # The page turns its own render errors into `st.error`, so an exception-free run is
    # not by itself a rendered report.
    assert not app.error, [element.value for element in app.error]
    return app


def _markdown(app) -> str:
    return "\n".join(block.value for block in app.markdown)


class TestItIsOnScreenAboveTheReport:
    """The wiring, read back off a running app rather than off the call site."""

    def test_the_recap_names_the_settings_beside_the_report(self):
        """Every line the recap builds reaches the page.

        Mutation: drop the ``render_filter_recap`` call from
        ``create_enhanced_data_table`` → nothing on the page names a filter, which is the
        state this ticket started from.
        """
        from utils import read_maf

        app = _run_page()
        rendered = _markdown(app)

        # Built from *this file's* columns, not from a complete set, so the expectation is
        # what the run actually did. The reference MAF carries no CIViC column — 100 of 100
        # do not — so the expected CIViC line reads *not applied*, and a recap that claimed
        # otherwise would fail here rather than passing on a substring.
        columns = list(read_maf(str(REFERENCE_MAF)).columns)
        expected = describe_run(pipeline_params("somatic"), columns)
        assert "Somatic analysis" in rendered
        for line in expected.lines:
            assert line.label in rendered, f"{line.label} is not on screen"
            assert line.value in rendered, f"{line.label}'s value is not on screen"

        assert NO_CIVIC in rendered, (
            "this MAF carries no CIViC column, so the recap must not report the CIViC "
            "selection as a restriction that applied"
        )

    def test_the_way_back_is_offered(self):
        """A recap with no route is half the request.

        Mutation: delete the button → the user is left reading what ran with no way to
        change it, which is the half the dev asked for explicitly.
        """
        labels = [button.label for button in _run_page().button]
        assert any("Change these filters" in label for label in labels), labels

    def test_a_stale_report_is_described_by_the_run_that_made_it(self):
        """The whole reason the recap is a snapshot.

        The user changes a threshold and does not re-filter. The table on screen is still
        the old cut, so the recap must still describe the old cut — and say that the
        settings have moved rather than let the two be read as one.

        Mutation: read ``st.session_state.filter_params`` at render time instead of the
        snapshot → the recap reports ``500 reads or more`` above a table cut at 50, with
        nothing on screen saying which is which.
        """
        app = _run_page(after='st.session_state.filter_params["min_depth"] = 500\n')
        rendered = _markdown(app)

        assert "50 reads or more" in rendered
        assert "500 reads or more" not in rendered
        assert any(
            "settings have changed" in warning.value.lower() for warning in app.warning
        ), [w.value for w in app.warning]
        assert any("Re-apply" in button.label for button in app.button)

    def test_an_unchanged_report_is_not_called_stale(self):
        """The ordinary case has no warning and no re-run button in the recap.

        Mutation: drop the stamp comparison and warn unconditionally → every report is
        announced as stale, which teaches the user to ignore the one that is.
        """
        app = _run_page()

        assert not any(
            "settings have changed" in warning.value.lower() for warning in app.warning
        ), [w.value for w in app.warning]
        assert not any("Re-apply" in button.label for button in app.button)


class TestWhatTheCutDid:
    """Issue #147's half of the recap: how many variants each setting kept out.

    Read off a running app, because the numbers are computed by the *filter run* and drawn
    by a different module several calls away — the seam #137 established and the one this
    reuses. A guard on ``attribute_report`` alone would pass over the block never being
    drawn, which is issue #140's finding about a message nobody reads.
    """

    def test_the_block_says_how_many_each_setting_kept_out(self):
        """The counts on screen are the ones the run computed, setting by setting.

        Mutation: drop the ``_render_attribution`` call from ``render_filter_recap`` →
        nothing on the page says what the cut did, which is the state this ticket started
        from. Mutation: have the page pass ``st.session_state.filter_params`` instead of the
        snapshot → the numbers follow the controls rather than the report.
        """
        from filters.attribution import attribute_report
        from utils import read_maf

        app = _run_page()
        rendered = _markdown(app)

        expected = attribute_report(
            read_maf(str(REFERENCE_MAF)), pipeline_params("somatic")
        )
        assert expected.excluded_by, "the reference MAF should leave variants out"

        assert f"{expected.rows:,} variants read" in rendered
        assert f"{expected.in_report:,} in the report" in rendered
        for label, count in expected.excluded_by:
            assert f"**{label}** — {count:,}" in rendered, (
                f"{label} is not on screen with its count"
            )

    def test_the_block_says_the_counts_overlap(self):
        """A reader who adds them and overshoots has been told something that looks wrong.

        The counts sum to more than the number left out, by design — so the sentence saying
        so is load-bearing rather than decoration, for the same reason
        ``_frequency_note`` states its exemption count even at zero.

        Mutation: delete the caption → the block reads as a partition and its arithmetic
        does not work.
        """
        from filters.attribution import attribute_report
        from utils import read_maf

        app = _run_page()
        captions = "\n".join(block.value for block in app.caption)
        expected = attribute_report(
            read_maf(str(REFERENCE_MAF)), pipeline_params("somatic")
        )

        assert sum(count for _, count in expected.excluded_by) > expected.left_out, (
            "the fixture must actually overlap, or this test asserts nothing"
        )
        assert f"more than {expected.left_out:,}" in captions, captions
