"""Tests for the app's MAF loader.

Since issue #32 the app has no loader of its own. ``utils.read_maf`` is a routing
layer: it accepts the two shapes the UI can hand it — a path on disk, or the buffer
of an uploaded file — and delegates every parse to ``vendor.pipeline_utils.read_maf``,
which is ``bin/utils.py:read_maf`` held byte-identical by the drift guard.

So most of what is asserted here is *equivalence*: what the app loads must be what the
pipeline loads, column for column and dtype for dtype. That is the point of the swap —
the app's own loader built its frame from ``line.split("\\t")``, so every column arrived
as a Python string and the pipeline's filters raised ``TypeError`` on the first numeric
comparison. ``TestVendoredFiltersAcceptAppLoadedFrames`` is the record that they no
longer do.

The deleted loader's quirks — padding short rows, truncating long ones, stripping
header whitespace — are not reproduced. Each was a silent divergence from the pipeline,
and ``TestTheDeletedLoaderStaysDeleted`` keeps them from creeping back.
"""

import io
import os
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from utils import read_maf, read_maf_without_comment_lines  # noqa: E402
from vendor.pipeline_utils import read_maf as pipeline_read_maf  # noqa: E402

FIXTURE_DIR = Path(__file__).resolve().parent / "fixtures" / "parity"
SOMATIC_REFERENCE = FIXTURE_DIR / "somatic_reference.maf"
GERMLINE_REFERENCE = FIXTURE_DIR / "germline_reference.maf"


class TestReadMafRoutesToTheVendoredReader(unittest.TestCase):
    """``utils.read_maf`` must be the vendored reader plus input plumbing, nothing more."""

    def setUp(self):
        self.sample_maf_content = """# MAF version 2.4
# This is a test MAF file
Hugo_Symbol	Variant_Classification	DP	tumor_f	gnomAD_exome_AF
BRCA1	Missense_Mutation	45	0.25	0.001
TP53	Nonsense_Mutation	52	0.38	0.0001
EGFR	Missense_Mutation	30	0.15	0.05
KRAS	Silent	25	0.08	0.15
PIK3CA	Frame_Shift_Del	5	0.02	.
"""

    def test_a_path_is_read_by_the_vendored_reader_verbatim(self):
        """Given a path, the app's loader must produce the pipeline's frame exactly.

        Asserted on a real reference MAF rather than a synthetic one, because dtype
        inference is what this ticket is about and the reference is where the awkward
        columns live.
        """
        pd.testing.assert_frame_equal(
            read_maf(str(SOMATIC_REFERENCE)), pipeline_read_maf(str(SOMATIC_REFERENCE))
        )

    def test_a_path_object_is_accepted(self):
        """``pathlib.Path`` reaches the loader from ``_auto_load_file_from_path``."""
        pd.testing.assert_frame_equal(
            read_maf(SOMATIC_REFERENCE), pipeline_read_maf(str(SOMATIC_REFERENCE))
        )

    def test_an_uploaded_text_buffer_gives_the_same_frame_as_the_file(self):
        """The upload path must not be a second loader with its own dtype rules."""
        buffer = io.StringIO(SOMATIC_REFERENCE.read_text(encoding="utf-8"))
        pd.testing.assert_frame_equal(
            read_maf(buffer), pipeline_read_maf(str(SOMATIC_REFERENCE))
        )

    def test_an_uploaded_bytes_buffer_gives_the_same_frame_as_the_file(self):
        """Streamlit's ``UploadedFile.getvalue()`` returns bytes, not str."""
        buffer = io.BytesIO(SOMATIC_REFERENCE.read_bytes())
        pd.testing.assert_frame_equal(
            read_maf(buffer), pipeline_read_maf(str(SOMATIC_REFERENCE))
        )

    def test_numeric_columns_arrive_numeric(self):
        """The property the whole ticket exists for: no more all-strings frames."""
        df = read_maf(io.StringIO(self.sample_maf_content))

        self.assertEqual(df["DP"].dtype, "int64")
        self.assertEqual(df["tumor_f"].dtype, "float64")

    def test_comment_lines_are_skipped(self):
        df = read_maf(io.StringIO(self.sample_maf_content))

        self.assertEqual(len(df), 5)
        self.assertEqual(list(df.columns)[0], "Hugo_Symbol")
        self.assertFalse(any(df["Hugo_Symbol"].str.startswith("#")))

    def test_dot_is_preserved_in_an_otherwise_numeric_column(self):
        """``.`` is MAF's missing marker and pandas does not read it as NaN.

        Worth pinning because it is the mechanism behind the ``somatic_dot_numeric``
        parity fixture: one ``.`` leaves the column as ``object``, and the pipeline's
        filters then raise on it. The app inherits that, deliberately — the alternative
        is the app quietly disagreeing with the pipeline about what the file says.
        """
        df = read_maf(io.StringIO(self.sample_maf_content))

        # Asserted as "not numeric" rather than "== object" on purpose: pandas 3.0
        # defaults strings to the new `str` dtype, and the assertion this replaced
        # (`test_app_loader_cannot_reach_parity`) went red on a fresh install for
        # spelling `object` while its intent still held. requirements.txt pins pandas
        # exactly as of issue #256, so a fresh install no longer walks this test into
        # 3.0 on its own — but a pin is a decision someone will revise, and the whole
        # point of writing the assertion to its intent is that it survives the bump.
        self.assertFalse(pd.api.types.is_numeric_dtype(df["gnomAD_exome_AF"]))
        self.assertTrue(any(df["gnomAD_exome_AF"] == "."))

    def test_an_unreadable_file_raises_so_the_ui_fallback_engages(self):
        """``data_loading`` distinguishes primary from fallback by the exception."""
        with self.assertRaises(pd.errors.EmptyDataError):
            read_maf(io.StringIO(""))


class TestVendoredFiltersAcceptAppLoadedFrames(unittest.TestCase):
    """The acceptance criterion of issue #32, stated as three assertions.

    Before the swap, every one of these raised ``TypeError: '>=' not supported between
    instances of 'str' and 'int'`` — the app's loader stringified ``t_alt_count`` and
    ``tumor_f``, so the pipeline's filters could not run on app-loaded data at all.
    ``tests/parity/harness.py`` had to load both sides itself to work around it.

    The masks are compared against the same functions run on a vendored-loader frame,
    not against hardcoded counts alone: equality of *verdicts* is the claim, and it
    stays true if the fixtures are ever regenerated.
    """

    def test_common_filters_runs_and_agrees(self):
        from vendor.pipeline_filters import common_filters

        app_frame = read_maf(str(SOMATIC_REFERENCE))
        vendored_frame = pipeline_read_maf(str(SOMATIC_REFERENCE))
        args = (50, ["Silent", "IGR", "RNA"])

        mask = common_filters(app_frame, *args)

        pd.testing.assert_series_equal(mask, common_filters(vendored_frame, *args))
        self.assertEqual(mask.sum(), 74)

    def test_somatic_filters_runs_and_agrees(self):
        from vendor.pipeline_filters import somatic_filters

        app_frame = read_maf(str(SOMATIC_REFERENCE))
        vendored_frame = pipeline_read_maf(str(SOMATIC_REFERENCE))
        args = (0.05, "null", ["Tier_I_strong"], ["A", "B"], ["IA"], ["Pathogenic"], True)

        specific, patho = somatic_filters(app_frame, *args)

        expected_specific, expected_patho = somatic_filters(vendored_frame, *args)
        pd.testing.assert_series_equal(specific, expected_specific)
        pd.testing.assert_series_equal(patho, expected_patho)
        self.assertEqual((int(specific.sum()), int(patho.sum())), (25, 19))

    def test_germline_filters_runs_and_agrees(self):
        from vendor.pipeline_filters import germline_filters

        app_frame = read_maf(str(GERMLINE_REFERENCE))
        vendored_frame = pipeline_read_maf(str(GERMLINE_REFERENCE))
        args = (0.05, "null", ["Pathogenic"], ["HP Pathogenic"], ["Pathogenic"])

        specific, patho = germline_filters(app_frame, *args)

        expected_specific, expected_patho = germline_filters(vendored_frame, *args)
        pd.testing.assert_series_equal(specific, expected_specific)
        pd.testing.assert_series_equal(patho, expected_patho)
        self.assertEqual((int(specific.sum()), int(patho.sum())), (19, 18))


class TestTheFallbackPathSharesTheParser(unittest.TestCase):
    """The recovery branch behind both ``except`` clauses in ``data_loading.py``.

    Two properties, and they pull in opposite directions, which is why both are pinned:

    * its frame must be *type-identical* to the primary path's — that is the half of
      issue #32's second criterion that sits after the comma ("so no path produces a
      differently-typed frame"), and the reason the branch no longer calls
      ``pd.read_csv`` itself;
    * it must still *disagree* with the primary path about where the header is, or it is
      a retry of a failure rather than a recovery from one.
    """

    CLEAN = "# MAF version 2.4\n# note\nHugo_Symbol\tDP\nBRCA1\t45\nTP53\t52\n"
    #: A comment line among the data rows — the shape the fallback exists for. The
    #: vendored reader counts all three ``#`` lines and reads the header from line 3,
    #: which is the first *data* row.
    STRAY_COMMENT = "# MAF version 2.4\nHugo_Symbol\tDP\nBRCA1\t45\n# stray\nTP53\t52\n"

    def test_it_agrees_with_the_primary_path_on_a_well_formed_file(self):
        pd.testing.assert_frame_equal(
            read_maf_without_comment_lines(self.CLEAN), read_maf(io.StringIO(self.CLEAN))
        )

    def test_it_accepts_bytes_as_well_as_str(self):
        """``UploadedFile.getvalue()`` hands the upload fallback bytes."""
        pd.testing.assert_frame_equal(
            read_maf_without_comment_lines(self.CLEAN.encode("utf-8")),
            read_maf_without_comment_lines(self.CLEAN),
        )

    def test_it_recovers_the_header_the_primary_path_loses(self):
        """The primary path mis-locates the header here; the fallback must not.

        This is the whole justification for the branch. Stripping only the *leading*
        comment block would not have achieved it: the vendored reader re-counts whatever
        it is handed, so line ``k + m`` of the original and line ``m`` of a
        ``k``-stripped remainder are the same physical line, and the retry would repeat
        the primary path's mistake exactly. Removing every comment line leaves the count
        at zero.
        """
        primary = read_maf(io.StringIO(self.STRAY_COMMENT))
        self.assertEqual(
            list(primary.columns),
            ["BRCA1", "45"],
            "the pipeline's comment count no longer overshoots on this shape — "
            "re-examine whether this fallback still has a job",
        )

        recovered = read_maf_without_comment_lines(self.STRAY_COMMENT)
        self.assertEqual(list(recovered.columns), ["Hugo_Symbol", "DP"])
        self.assertEqual(list(recovered["Hugo_Symbol"]), ["BRCA1", "TP53"])

    def test_recovered_frames_are_typed_like_the_primary_path(self):
        """Recovery must not cost the dtypes — the entire point of the swap."""
        recovered = read_maf_without_comment_lines(self.STRAY_COMMENT)

        self.assertEqual(recovered["DP"].dtype, "int64")


class TestTheDeletedLoaderStaysDeleted(unittest.TestCase):
    """One loader, and none of the old one's silent reshaping."""

    def test_no_second_parser_survives_in_the_utils_package(self):
        import utils
        import utils.main_utils as main_utils

        self.assertNotIn("read_maf_simple", utils.__all__)
        self.assertFalse(
            hasattr(main_utils, "read_maf"),
            "main_utils.read_maf is back — a second loader with its own parse rules",
        )

    def test_the_drifted_column_list_is_gone(self):
        """``main_utils.KEEP`` was a stale copy of the pipeline's ``KEEP``.

        39 entries against the pipeline's 45, and nothing read it. The live copy is
        ``vendor.pipeline_utils.KEEP``, which the drift guard holds to ``bin/``.
        """
        import utils
        import utils.main_utils as main_utils

        self.assertFalse(hasattr(main_utils, "KEEP"))
        self.assertFalse(hasattr(utils, "KEEP"))

        from vendor.pipeline_utils import KEEP

        self.assertEqual(len(KEEP), 45)

    def test_a_malformed_row_is_read_the_way_the_pipeline_reads_it(self):
        """The old loader padded short rows and truncated long ones, silently.

        Truncation looked harmless and was not: it produced a plausible frame that the
        file did not state, and one the pipeline would never have produced. pandas does
        something different and also lossy — an extra field makes it promote the first
        column into the index, shifting the row — but it does it *identically on both
        sides*, which is the only property the app can offer here. Pinned by
        equivalence rather than by a preferred answer, because the answer is the
        pipeline's to choose.
        """
        malformed = (
            "# MAF version 2.4\n"
            "Hugo_Symbol\tVariant_Classification\n"
            "BRCA1\tMissense_Mutation\tExtraColumn\n"
        )
        with tempfile.NamedTemporaryFile(
            "w", suffix=".maf", delete=False, encoding="utf-8"
        ) as handle:
            handle.write(malformed)
            path = handle.name
        try:
            frame = read_maf(io.StringIO(malformed))
            pd.testing.assert_frame_equal(frame, pipeline_read_maf(path))
        finally:
            os.unlink(path)

        self.assertEqual(list(frame.index), ["BRCA1"])


class TestUtilsIntegration(unittest.TestCase):
    """Integration tests for utility functions."""

    def test_read_maf_with_realistic_data(self):
        """Test with more realistic MAF data structure."""
        realistic_maf = """# MAF version 2.4
# Source: Test Suite
Hugo_Symbol	Entrez_Gene_Id	Variant_Classification	Variant_Type	DP	tumor_f	gnomAD_exome_AF
BRCA1	672	Missense_Mutation	SNP	45	0.25	0.001
TP53	7157	Nonsense_Mutation	SNP	52	0.38	0.0001
EGFR	1956	Missense_Mutation	SNP	30	0.15	0.05
"""

        maf_file = io.StringIO(realistic_maf)
        df = read_maf(maf_file)

        self.assertEqual(len(df), 3)
        self.assertEqual(len(df.columns), 7)
        self.assertTrue(
            all(col in df.columns for col in ["Hugo_Symbol", "DP", "tumor_f"])
        )


if __name__ == "__main__":
    unittest.main()
