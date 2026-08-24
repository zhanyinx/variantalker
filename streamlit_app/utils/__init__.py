"""MAFigate Utils Module

The app's MAF loader lives here, and it is deliberately thin. The *loader* is the
pipeline's — ``vendor.pipeline_utils.read_maf``, held byte-identical to
``bin/utils.py`` by ``tests/test_vendor_drift.py``. This module only adapts the input
shapes the UI can produce (a path from "Open With", a buffer from ``st.file_uploader``)
onto the path-taking reader.

Why it matters that there is exactly one loader: the app used to build its frame from
``line.split("\\t")``, which made every column a Python string. The pipeline's filters
compare ``t_alt_count + t_ref_count`` against a depth and ``tumor_f`` against a VAF, so
on such a frame they raised ``TypeError`` before deciding anything — the app could not
have run the pipeline's decision even if it wanted to (issue #32).
"""

import io
import os
import tempfile
from contextlib import contextmanager

import pandas as pd

from vendor.pipeline_utils import read_maf as _pipeline_read_maf


def read_maf(file) -> pd.DataFrame:
    """Read a MAF exactly as the pipeline reads it.

    Args:
        file: A path (``str`` or ``os.PathLike``), or a file-like buffer whose
            ``getvalue()``/``read()`` returns ``str`` or ``bytes`` — Streamlit's
            ``UploadedFile`` returns bytes.

    Returns:
        pandas.DataFrame: the frame ``vendor.pipeline_utils.read_maf`` returns for the
        same bytes. Column names are *not* normalised and ``.`` is *not* read as NaN,
        because the pipeline does neither; a frame that differed in either would make
        the app disagree with the pipeline about what the file says.

    Raises:
        pandas.errors.EmptyDataError, pandas.errors.ParserError: as the vendored reader
            raises them. ``page_modules/data_loading.py`` catches these and retries via
            :func:`read_maf_without_comment_lines`.
    """
    if hasattr(file, "getvalue") or hasattr(file, "read"):
        with _spilled_to_disk(_content_of(file)) as path:
            return _pipeline_read_maf(path)
    return _pipeline_read_maf(file)


def read_maf_without_comment_lines(content) -> pd.DataFrame:
    """Retry a MAF that :func:`read_maf` raised on, with every comment line removed.

    Args:
        content: the file's full text, as ``str`` or ``bytes``.

    This is the recovery branch behind both ``except`` clauses in
    ``page_modules/data_loading.py``. It parses with the same vendored reader — so its
    dtypes, missing-value handling and column names are the primary path's, which is the
    point: those branches used to call ``pd.read_csv`` themselves and could hand the
    filters an all-strings frame the primary path would never produce.

    What it does differently is *where the header is*, and the difference is the whole
    reason to have it. The vendored reader inherits the pipeline's rule: count every
    ``#`` line in the file, then treat line ``n`` as the header. That is correct only
    while the comments form one leading block. Drop ``#`` lines anywhere else and the
    count overshoots, the real header is consumed as a comment, and the first data row
    becomes the header — measured on ``# lead / A\\tB / 1\\t2 / # stray / 3\\t4``, the
    primary path returns columns ``["1", "2"]``.

    Removing *all* comment lines rather than just the leading block leaves the reader
    counting zero, so it reads the header where the header actually is. Stripping only
    the leading block would not have been enough: the reader re-counts what it is given,
    so line ``k + m`` of the original and line ``m`` of a ``k``-stripped remainder are
    the same physical line, and the retry would be a no-op that merely repeats the
    primary path's mistake.

    The one behavioural cost is that a data row whose first field begins with ``#`` is
    dropped rather than kept as data. That is the right trade for a recovery path: such a
    row is a comment by every convention MAF has, and the alternative is a frame with a
    ``#``-prefixed ``Hugo_Symbol``.

    Pinned by ``tests/test_utils.py::TestTheFallbackPathSharesTheParser``; measured under
    issue #32 (a private-tree note — the runnable half of it now lives beside its fixtures
    as ``tests/run_app_check.py``, per issue #345).
    """
    text = content.decode("utf-8") if isinstance(content, bytes) else content
    lines = [line for line in text.strip().split("\n") if not line.startswith("#")]

    return read_maf(io.StringIO("\n".join(lines)))


def _content_of(buffer) -> bytes:
    """The whole buffer, as bytes, regardless of how it was opened.

    ``getvalue()`` is preferred over ``read()`` because it ignores the cursor: an
    uploaded file may already have been read once by validation code upstream.
    """
    raw = buffer.getvalue() if hasattr(buffer, "getvalue") else buffer.read()
    return raw.encode("utf-8") if isinstance(raw, str) else raw


@contextmanager
def _spilled_to_disk(content: bytes):
    """Hand the vendored reader a path, because a path is all it takes.

    It opens the file twice — once to count the comment lines, once to parse — so a
    buffer cannot be passed through. Spilling to disk keeps this module a routing layer:
    re-implementing the comment count against an in-memory buffer would be a second
    loader again, free to drift from ``bin/`` exactly as the deleted one did.
    """
    handle = tempfile.NamedTemporaryFile(suffix=".maf", delete=False)
    try:
        handle.write(content)
        handle.close()
        yield handle.name
    finally:
        os.unlink(handle.name)


# Export all functions
__all__ = [
    "read_maf",
    "read_maf_without_comment_lines",
]
