"""The harness's loading premise: its frame is the frame the app loads.

``harness.py`` hands **one** DataFrame to both sides of every comparison, loaded by the
vendored ``pipeline_utils.read_maf``. That is what makes each divergence it reports a
divergence in *filter logic* rather than in I/O — and it is only honest while the app
loads the same way.

It began as a workaround rather than a premise. Until issue #32 the app had a loader of
its own, building its frame from ``line.split("\\t")``: every column a Python string,
nothing dtype-inferred, so the pipeline's ``common_filters`` raised
``TypeError: '>=' not supported between instances of 'str' and 'int'`` on the first
comparison (issue #10). The harness could not have used it, and ``test_parity.py``
recorded that impossibility as ``test_app_loader_cannot_reach_parity``. Issue #32
deleted that loader — ``utils.read_maf`` now delegates to the same vendored reader — so
the workaround became an identity, and the test became this one.

Asserted per fixture, because it is a claim about every file the parity suite measures,
not about one of them.

Deliberately **not** under ``test_parity.py``'s ``skipif(not PIPELINE_AVAILABLE)``, for
the same reason as ``test_baseline_integrity.py``: this compares two in-process
readers over committed fixtures, with no ``bin/`` and no pipeline subprocess anywhere in
it. Left in that module it would have been skipped in exactly the pipeline-less CI job
where issue #24 found the parity suite silently disappearing — and this is the assertion
that the rest of the suite's frames can be trusted.

That the vendored *filters* then run on such a frame is a different claim, pinned by
``tests/test_utils.py::TestVendoredFiltersAcceptAppLoadedFrames``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

STREAMLIT_APP = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(STREAMLIT_APP))

from tests.parity import harness as H  # noqa: E402
from tests.parity.contract import CASES  # noqa: E402


@pytest.mark.parametrize("fixture", sorted({c.fixture for c in CASES}))
def test_harness_loads_the_frame_the_app_would_load(fixture):
    """The app's loader and the harness's must return the same frame, exactly."""
    from utils import read_maf as app_read_maf

    path = H.FIXTURE_DIR / fixture
    pd.testing.assert_frame_equal(app_read_maf(str(path)), H.read_maf(path))
