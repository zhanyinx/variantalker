"""
Test configuration management and constants.
"""

import types
import unittest
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import config.constants as constants
from config.constants import APP_NAME, APP_VERSION, PAGES
from config.columns import REQUIRED_COLUMNS, OPTIONAL_COLUMNS


def _flatten_columns(grouped):
    """Flatten a {category: [columns]} mapping into a single list of columns."""
    return [col for cols in grouped.values() for col in cols]


#: Everything ``config/constants.py`` is allowed to hold.
#:
#: The file was 4,021 lines and five kinds of thing, and it got that way one harmless
#: addition at a time — a constant with nowhere obvious to live went here, because "constants"
#: is a name that accepts anything. Its docstring now states a rule narrow enough to check, so
#: this list is the check: a constant about *variants* belongs to ``config/vocabularies.py``,
#: ``config/presets.py``, ``config/gene_panels.py`` or ``config/columns.py``, and only what is
#: about the application as an application belongs here.
#:
#: Adding a name to this tuple is therefore a deliberate act. If it is not the app's identity,
#: its navigation, or a presentation default, the honest move is a different module rather than
#: a longer tuple.
#:
#: ``APP_TAGLINE`` was added deliberately (issue #71). It is the app's identity in the most
#: literal sense — the one sentence it uses to say what it is — and it is here rather than in
#: a page module because two surfaces read it: the header drawn on every page, and the About
#: dialog. It replaces three hand-written variants of the same sentence.
ALLOWED_IN_CONSTANTS = ("APP_NAME", "APP_TAGLINE", "APP_VERSION", "PAGES", "TABLE_HEIGHT")


def test_constants_holds_only_the_apps_own_chrome():
    """The junk-drawer guard: no unlisted public name in ``config/constants.py``.

    Imported modules are skipped rather than counted. The smell being guarded against is a
    *constant* filed here for want of a better home; an ``import json`` is not that, and
    counting it would fail the check with a message about variant vocabularies that had
    nothing to do with the edit.
    """
    public = tuple(
        sorted(
            name
            for name, value in vars(constants).items()
            if not name.startswith("_") and not isinstance(value, types.ModuleType)
        )
    )
    assert public == tuple(sorted(ALLOWED_IN_CONSTANTS)), (
        f"config/constants.py now holds {public}, but only {tuple(sorted(ALLOWED_IN_CONSTANTS))} "
        "are about the application as an application. Anything about variants — a vocabulary, a "
        "preset, a gene panel, a column category — has a module of its own; see the docstring on "
        "config/constants.py. If the new name really is app chrome, add it here in the same "
        "commit and say why."
    )


class TestConstants(unittest.TestCase):
    """Test application constants."""

    def test_app_info_constants(self):
        """Test that basic app constants are defined."""
        self.assertIsInstance(APP_NAME, str)
        self.assertIsInstance(APP_VERSION, str)
        self.assertTrue(len(APP_NAME) > 0)
        self.assertTrue(len(APP_VERSION) > 0)

    def test_pages_configuration(self):
        """Test that page configuration is properly defined."""
        self.assertIsInstance(PAGES, dict)
        self.assertTrue(len(PAGES) > 0)

        # Check that all pages have required keys
        required_keys = {"title", "icon"}
        for page_key, page_config in PAGES.items():
            self.assertIsInstance(page_config, dict)
            for key in required_keys:
                self.assertIn(key, page_config)
                self.assertIsInstance(page_config[key], str)

    # There were two more tests here, on the parameter defaults, and issue #77 deleted
    # them with the constant they read. Both went through `DEFAULT_PARAMS`, which no app
    # module ever read, and both are asserted better elsewhere:
    #
    # * the contract's keys and values, against `nextflow.config` and for *both* arms, by
    #   `test_param_contract.py` — which skips where `nextflow.config` is absent, so in a
    #   packaged tree what holds the contract to its values is `test_filter_entry_point.py`,
    #   ungated by design, which filters at `pipeline_params(arm)` against pipeline verdicts
    #   frozen in git. That is a stronger claim than the `assertIsInstance` deleted here.
    # * the four presets, by
    #   `test_param_contract.test_soft_and_clinical_presets_still_exist_as_deviations`,
    #   which carries no gate and now covers all three keys the germline half looked for.


class TestColumns(unittest.TestCase):
    """Test column configuration."""

    def test_required_columns_exist(self):
        """Test that required columns are defined.

        REQUIRED_COLUMNS is a {category: [columns]} mapping.
        """
        self.assertIsInstance(REQUIRED_COLUMNS, dict)
        all_required = _flatten_columns(REQUIRED_COLUMNS)
        self.assertTrue(len(all_required) > 0)

        # Check for essential MAF columns
        essential_columns = ["Hugo_Symbol", "Variant_Classification"]
        for col in essential_columns:
            self.assertIn(col, all_required)

    def test_optional_columns_exist(self):
        """Test that optional columns are defined."""
        self.assertIsInstance(OPTIONAL_COLUMNS, dict)
        self.assertTrue(len(_flatten_columns(OPTIONAL_COLUMNS)) > 0)

    def test_no_column_overlap(self):
        """Test that required and optional columns don't overlap."""
        required_set = set(_flatten_columns(REQUIRED_COLUMNS))
        optional_set = set(_flatten_columns(OPTIONAL_COLUMNS))

        overlap = required_set.intersection(optional_set)
        self.assertEqual(len(overlap), 0, f"Found overlapping columns: {overlap}")

    def test_columns_are_strings(self):
        """Test that all column names are strings."""
        all_columns = _flatten_columns(REQUIRED_COLUMNS) + _flatten_columns(
            OPTIONAL_COLUMNS
        )
        for col in all_columns:
            self.assertIsInstance(col, str)
            self.assertTrue(len(col) > 0)


# A `TestConfigurationIntegration` class stood here with four more tests, and issue #77
# deleted it whole. Every one of them read `DEFAULT_PARAMS` behind an `if param in
# DEFAULT_PARAMS` guard, so each asserted the type of whatever happened to be there and
# skipped in silence otherwise — over a constant that, being a deep copy of the contract,
# could only ever hold exactly what the contract holds.
#
# What covers those keys now is set out above. The honest residue: nothing asserts the
# *types* of the contract's values as types. `test_filter_entry_point.py` and
# `test_param_migration.py` both read them ungated, but they read them for what they mean,
# so a wrong type surfaces as a filter or migration failure rather than as an
# `assertIsInstance`. That is the trade this deletion makes, and it is the right way round:
# the deleted checks skipped whenever the key they named was absent, which is the one case
# a shape guard exists to catch.


if __name__ == "__main__":
    unittest.main()


class TestSpelledIn(unittest.TestCase):
    """``config.columns.spelled_in`` — a column's name as *this file* writes it (issue #190).

    Somewhere between the pipeline and the MAF a step passes the header through R's ``make.names``,
    which turns every character illegal in an identifier into a dot. Issue #189 met it on the
    evidence column and #190 on ``GERP++_RS``, spelled ``GERP.._RS`` on 2 of 167 real MAFs.
    """

    def test_the_canonical_spelling_is_returned_when_the_file_uses_it(self):
        from config.columns import spelled_in

        self.assertEqual(spelled_in(["Chromosome", "GERP++_RS"], "GERP++_RS"), "GERP++_RS")

    def test_a_mangled_spelling_is_found(self):
        from config.columns import spelled_in

        self.assertEqual(spelled_in(["Chromosome", "GERP.._RS"], "GERP++_RS"), "GERP.._RS")

    def test_a_column_the_file_does_not_have_is_none(self):
        from config.columns import spelled_in

        self.assertIsNone(spelled_in(["Chromosome"], "GERP++_RS"))

    def test_the_canonical_spelling_wins_over_a_mangled_sibling(self):
        from config.columns import spelled_in

        columns = ["GERP.._RS", "GERP++_RS"]
        self.assertEqual(spelled_in(columns, "GERP++_RS"), "GERP++_RS")

    def test_an_ambiguous_mangled_spelling_is_refused_rather_than_guessed(self):
        """Two columns mangling to one shape, and neither is the documented name.

        An already-mangled ``GERP.._RS`` and a ``GERP+-_RS`` both reduce to ``GERP.._RS``, so a
        file carrying both gives no way to tell which the caller meant. Returning either would put
        an arbitrary number under a named predictor on a clinical screen; ``None`` sends the caller
        down its already-required *absent* path instead.

        The real pair does **not** collide, and that is why the comparison is on the whole
        normalised name rather than a substring: ANNOVAR ships ``GERP++_RS`` beside ``GERP++_NR``,
        which reduce to ``GERP.._RS`` and ``GERP.._NR`` — distinct. A ``GERP``-prefix search would
        have matched both.
        """
        from config.columns import spelled_in

        self.assertIsNone(spelled_in(["GERP.._RS", "GERP+-_RS"], "GERP++_RS"))
        self.assertEqual(
            spelled_in(["GERP.._RS", "GERP.._NR"], "GERP++_RS"), "GERP.._RS"
        )

    def test_a_name_needing_no_mangling_still_resolves(self):
        from config.columns import spelled_in

        self.assertEqual(spelled_in(["SIFT_score"], "SIFT_score"), "SIFT_score")
        self.assertIsNone(spelled_in(["SIFT_pred"], "SIFT_score"))
