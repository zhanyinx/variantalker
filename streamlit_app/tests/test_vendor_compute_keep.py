"""Behaviour of the vendored output-column computation, `vendor.compute_keep`.

`tests/test_vendor_drift.py` proves `compute_keep` is still a faithful copy of the
statements in `bin/filter_variants.py:main()`. This module proves the one thing the
copy is deliberately *not* faithful about: the pipeline aliases its module-level
column list and then mutates it, and `compute_keep` takes a local copy first.

That difference is the whole reason this symbol exists as a separate ticket. In a
pipeline process the aliasing is harmless — `main()` runs once and exits. In the app
the module stays loaded for the life of the session, so the mutations accumulate:

* a germline run removes three entries from the shared list and substitutes a fourth,
  so a **somatic** run that follows silently returns germline columns;
* a second **germline** run raises `ValueError: list.remove(x): x not in list`,
  because the entries it removes are already gone.

Both are reproduced below against the fixed function, which must not exhibit either.
"""

import sys
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

STREAMLIT_APP = Path(__file__).resolve().parents[1]
if str(STREAMLIT_APP) not in sys.path:
    sys.path.insert(0, str(STREAMLIT_APP))

from vendor.pipeline_filters import compute_keep  # noqa: E402
from vendor.pipeline_utils import KEEP  # noqa: E402

#: Columns a MAF carries beyond `KEEP` that the computation reacts to. Both gnomAD
#: pairs use the `_raw` spelling, which is the branch the pipeline prefers.
EXTRA_MAF_COLUMNS = ["gnomAD_exome_AF_raw", "gnomAD_genome_AF_raw"]


def args_shim(sample_type: str, skip_civic: bool = False) -> SimpleNamespace:
    """The pipeline's `argparse` namespace, reduced to what the statements read.

    `compute_keep` is vendored verbatim, so it reads `args.sample_type` and
    `args.skip_civic` exactly as `main()` does. Supplying that namespace is the
    caller's job — this is the call-site shim the app will build too.
    """
    return SimpleNamespace(sample_type=sample_type, skip_civic=skip_civic)


def maf_frame(columns=None) -> pd.DataFrame:
    """An empty frame carrying `columns`, defaulting to a full CIViC + gnomAD MAF.

    Only `out.columns` is ever read, so the rows are irrelevant.
    """
    if columns is None:
        columns = list(KEEP) + EXTRA_MAF_COLUMNS
    return pd.DataFrame({name: pd.Series(dtype="object") for name in columns})


@pytest.fixture
def frame() -> pd.DataFrame:
    return maf_frame()


@pytest.fixture(autouse=True)
def keep_is_pristine():
    """Fail loudly if any test in this module leaves the shared `KEEP` mutated.

    Without this, a regression in `compute_keep` would leak across test functions and
    show up as an unrelated test failing — which is precisely how the bug behaves in
    the app, and precisely what makes it hard to attribute.
    """
    before = list(KEEP)
    yield
    assert KEEP == before, (
        "vendor.pipeline_utils.KEEP was mutated by this test. compute_keep must copy "
        "the list before running the pipeline's mutating statements."
    )


def test_germline_call_returns_germline_columns(frame):
    """The germline substitutions happen at all — the control for the tests below."""
    keep = compute_keep(args_shim("germline"), frame)

    assert "InterVar" in keep
    assert "CancerVar" not in keep
    assert "RENOVO_Class" in keep and "RENOVO_pls" in keep
    assert "Tumor_Sample_Barcode" not in keep
    assert not [c for c in keep if c.startswith("CIViC")]


def test_somatic_call_returns_somatic_columns(frame):
    """The somatic arm keeps what germline strips — the other control."""
    keep = compute_keep(args_shim("somatic"), frame)

    assert "CancerVar" in keep
    assert "InterVar" not in keep
    assert "Tumor_Sample_Barcode" in keep
    assert "cosmic" in keep and "Freq_ExAC_ALL" in keep
    assert [c for c in keep if c.startswith("CIViC")]


def test_germline_then_somatic_returns_somatic_columns(frame):
    """A germline call must not contaminate the somatic call that follows it.

    The pipeline's `keep = KEEP` makes this fail: the germline arm removes
    `Tumor_Sample_Barcode`, `cosmic` and `Freq_ExAC_ALL` from the shared list and
    substitutes `InterVar` for `CancerVar` in place, so the somatic caller receives a
    germline column set and nothing says so.
    """
    compute_keep(args_shim("germline"), frame)
    somatic = compute_keep(args_shim("somatic"), frame)

    assert somatic == compute_keep(args_shim("somatic"), maf_frame()), (
        "the somatic column set depends on whether a germline call came first"
    )
    assert "CancerVar" in somatic, "somatic run lost CancerVar to a previous germline run"
    assert "InterVar" not in somatic, "somatic run emitted the germline InterVar column"
    for column in ("Tumor_Sample_Barcode", "cosmic", "Freq_ExAC_ALL"):
        assert column in somatic, f"somatic run lost {column} to a previous germline run"
    for column in ("RENOVO_Class", "RENOVO_pls"):
        assert column not in somatic, f"somatic run emitted the germline column {column}"


def test_two_germline_calls_do_not_raise(frame):
    """Two germline calls in one process must both succeed, and agree.

    With the pipeline's aliasing this raises `ValueError: list.remove(x): x not in
    list` on the second call, because `Tumor_Sample_Barcode` was removed by the first.
    """
    first = compute_keep(args_shim("germline"), frame)
    second = compute_keep(args_shim("germline"), frame)

    assert first == second


def test_repeated_calls_across_arms_are_all_independent(frame):
    """Arm order must not matter, in either direction or over many calls."""
    somatic_first = compute_keep(args_shim("somatic"), frame)
    germline_first = compute_keep(args_shim("germline"), frame)

    for _ in range(3):
        assert compute_keep(args_shim("germline"), frame) == germline_first
        assert compute_keep(args_shim("somatic"), frame) == somatic_first


def test_returned_list_is_not_the_shared_keep(frame):
    """The caller gets a list it may mutate without touching the module-level one."""
    keep = compute_keep(args_shim("somatic"), frame)

    assert keep is not KEEP
    keep.append("a_column_the_caller_added")
    assert "a_column_the_caller_added" not in KEEP


def test_mutating_a_returned_list_does_not_affect_the_next_call(frame):
    """The copy is per call, not one copy shared by every call."""
    first = compute_keep(args_shim("somatic"), frame)
    first.clear()

    assert compute_keep(args_shim("somatic"), frame), "second call returned an empty list"


def test_skip_civic_drops_the_civic_columns(frame):
    """`--skip_civic` removes the CIViC block on the somatic arm."""
    keep = compute_keep(args_shim("somatic", skip_civic=True), frame)

    assert not [c for c in keep if c.startswith("CIViC")]
    assert "CancerVar" in keep


def test_absent_civic_columns_drop_the_civic_block_without_the_flag(frame):
    """A MAF with no CIViC columns drops them even when `--skip_civic` was not given."""
    columns = [c for c in KEEP if not c.startswith("CIViC")] + EXTRA_MAF_COLUMNS
    keep = compute_keep(args_shim("somatic"), maf_frame(columns))

    assert not [c for c in keep if c.startswith("CIViC")]


@pytest.mark.parametrize(
    ("present", "expected", "unexpected"),
    [
        (["gnomAD_exome_AF_raw"], "gnomAD_exome_AF_raw", "gnomAD_exome_AF"),
        (["gnomAD_exome_AF"], "gnomAD_exome_AF", "gnomAD_exome_AF_raw"),
        (["gnomAD_genome_AF_raw"], "gnomAD_genome_AF_raw", "gnomAD_genome_AF"),
        (["gnomAD_genome_AF"], "gnomAD_genome_AF", "gnomAD_genome_AF_raw"),
    ],
)
def test_gnomad_columns_are_appended_when_present(present, expected, unexpected):
    """Each gnomAD pair appends whichever spelling the MAF carries, `_raw` first."""
    keep = compute_keep(args_shim("somatic"), maf_frame(list(KEEP) + present))

    assert expected in keep
    assert unexpected not in keep


def test_absent_gnomad_columns_are_not_appended():
    """A MAF with no gnomAD columns gets neither name."""
    keep = compute_keep(args_shim("somatic"), maf_frame(list(KEEP)))

    assert not [c for c in keep if c.startswith("gnomAD")]


def test_gnomad_columns_do_not_accumulate_across_calls(frame):
    """The append statements must not stack a second gnomAD pair onto a shared list.

    This is the aliasing bug's quietest symptom: `keep.append` mutates in place, so an
    aliased list grows by two entries per call and the column set gains duplicates
    rather than losing entries.
    """
    first = compute_keep(args_shim("somatic"), frame)
    compute_keep(args_shim("somatic"), frame)
    third = compute_keep(args_shim("somatic"), frame)

    assert third == first
    assert len(third) == len(set(third)), f"duplicate columns: {sorted(third)}"
