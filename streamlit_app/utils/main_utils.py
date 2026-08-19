"""
Main utility functions for MAFigate
Moved from standalone utils.py

The MAF loader used to live here too, alongside a ``KEEP`` column list. Both are gone
(issue #32): the loader was a second parser competing with the pipeline's, and ``KEEP``
was a copy of ``bin/utils.py:KEEP`` that had drifted to 39 entries against the
pipeline's 45 while nothing read it. The live copies are ``utils.read_maf`` — which
delegates to the vendored reader — and ``vendor.pipeline_utils.KEEP``.

Two serialisers lived here too, and went with the results page's third door (issue #83).
``convert_df`` was read only by ``create_export_buttons``, the CSV/TSV/Excel trio nothing
called; ``download_csv`` built a base64 ``<a href>``, the idiom ``st.download_button``
replaced, and was read by nothing at all. Deleting them took ``base64``, ``pandas`` and
``streamlit`` out of this module's imports with them.
"""

import hashlib


def params_hash(params: dict) -> str:
    """Stamp a parameter set, so a later render can tell whether it has changed.

    One recipe in one place. The data page stamps ``data_params_hash`` with it when the
    filters run and compares against it to warn that the settings have moved on; the
    sidebar's load status compares against the same stamp so it cannot call a stale
    result "current" while the data page calls it stale (issue #58).
    """
    return hashlib.md5(str(sorted(params.items())).encode()).hexdigest()


def settings_have_moved(params: dict, stamp: str) -> bool:
    """Whether ``params`` differ from the set ``stamp`` was taken of.

    The comparison, not just the recipe. Three surfaces ask this question — the sidebar's
    load status, the load section's warning, and the filter recap above the results (issue
    #137) — and each had written out the same three-part conjunction, including the part
    that is easy to get wrong: **no stamp means no run to be stale against**, not "stale".
    A fourth spelling of it is how one of them would come to call a report current while
    another calls it stale, which is the disagreement issue #58 fixed between the first
    two.

    Streamlit-free, like the rest of this module: the caller passes the two session values
    in rather than this reaching for them.
    """
    return bool(stamp and params and params_hash(params) != stamp)
