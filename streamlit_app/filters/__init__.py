"""The app's filter package: one entry point, and the vocabulary of its output.

Everything the filter decides is decided by ``vendor/``, the pipeline's own code. What is
re-exported here is the seam (:func:`apply_filters`), the two columns it labels rows
with, the four reasons it can give, and the one refusal it can raise instead — nothing
that computes a verdict.
"""

from .numeric_columns import UnreadableNumericColumns
from .variant_filters import (
    MAFIGATE_FILTER,
    MAFIGATE_REASON,
    NOPASS,
    PASS,
    REASON_BOTH,
    REASON_CRITERIA,
    REASON_REJECTED,
    REASON_RESCUE,
    Diagnostics,
    apply_filters,
    decomposition,
)

__all__ = [
    "apply_filters",
    "decomposition",
    "Diagnostics",
    "MAFIGATE_FILTER",
    "MAFIGATE_REASON",
    "PASS",
    "NOPASS",
    "REASON_CRITERIA",
    "REASON_RESCUE",
    "REASON_BOTH",
    "REASON_REJECTED",
    "UnreadableNumericColumns",
]
