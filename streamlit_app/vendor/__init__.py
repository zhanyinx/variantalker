"""Vendored copies of the VarianTalker pipeline's variant-filter code.

Everything in this package except ``_sync.py`` and this docstring is a copy of
code under ``bin/``, held byte-for-byte (``pipeline_utils.py``) or symbol-for-symbol
(``pipeline_filters.py``) so that the app's PASS/NOPASS decision is the pipeline's
decision rather than a re-derivation of it.

**Do not edit these modules.** Change ``bin/`` (which this effort treats as frozen)
and re-run ``python streamlit_app/vendor/_sync.py``. Any hand-edit is caught by
``streamlit_app/tests/test_vendor_drift.py``, which also runs in CI, as a
pre-commit hook, and as a gate on the installer builds.
"""
