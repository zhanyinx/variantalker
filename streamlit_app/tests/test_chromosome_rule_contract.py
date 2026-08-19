"""Drift guard for the chromosome-spelling rule (issue #149).

``config/chromosome_spelling.py`` does not invent a rule. It holds a **copy** of one the
pipeline already applies, in ``bin/add_guidelines_escat_to_funcotator.py``::

    if not any(["chr" in str(x) for x in maf["Chromosome"].values]):
        maf["Chromosome"] = "chr" + maf["Chromosome"].astype(str)

That is the same relationship ``config/pipeline_params.py`` has to ``nextflow.config``
and ``vendor/`` has to ``bin/``, and it fails the same way: silently. If the pipeline's
rule moves and the app's does not, MAFigate spells chromosomes differently from the
report the file it is reading came out of, and nothing raises — the app would simply be
quietly wrong about a value on every screen that draws one.

**Why this is a fifth module allowed to need the pipeline**, when
``test_suite_organisation.MODULES_ALLOWED_TO_NEED_THE_PIPELINE`` says the list is meant
to shrink: because the claim being checked is *exactly* "the app agrees with the
pipeline", and there is no way to check that without the pipeline. It is kept in its own
module rather than folded into ``tests/test_chromosome_spelling.py`` so that the twenty
assertions about the app's own behaviour keep holding in a checkout with no pipeline in
it, which is the whole point of the ``bin/``-free net. It is not folded into
``test_param_contract.py`` either: that module states, and keeps, the property of
importing neither pandas nor the app's filter code, and this comparison is behavioural
and needs both.

The comparison is behavioural rather than textual on purpose. The two are written in
different idioms — a list comprehension over ``.values`` against ``Series.str.contains``
— so no AST or hash equality could hold between them, and pinning the *text* would fail
on a reformat while passing on a reversal. What is pinned is the **decision**: for the
same column, both say prefix, or both say leave it.

**With one declared exception, asserted here rather than only described elsewhere.** On a
column that says nothing at all, the two rules *disagree*, and deliberately: the pipeline
would write ``chrnan`` into every cell, and ``config/chromosome_spelling.py`` refuses to,
because ``config.missing_values.says_nothing`` does not recognise that string and the
detail panel would print it at a clinician where it prints an em dash. That is the same
deviation the module documents for a blank cell in an otherwise-populated column, and it
is the only one. It is asserted by
:func:`test_the_one_place_the_app_declines_to_follow_the_pipeline` — declared exceptions
live *at the guard*, the way ``vendor/README.md`` keeps its one unfaithfulness beside the
thing it is unfaithful about, so that a reader of this file cannot come away believing the
agreement is total.
"""

from __future__ import annotations

import ast
import sys
from pathlib import Path

import pandas as pd
import pytest

STREAMLIT_APP = Path(__file__).resolve().parents[1]
if str(STREAMLIT_APP) not in sys.path:
    sys.path.insert(0, str(STREAMLIT_APP))

REPO_ROOT = STREAMLIT_APP.parent
ANNOTATION_SCRIPT = REPO_ROOT / "bin" / "add_guidelines_escat_to_funcotator.py"

needs_annotation_script = pytest.mark.skipif(
    not ANNOTATION_SCRIPT.is_file(),
    reason="pipeline bin/ not present (packaged app or partial checkout)",
)

#: Columns to put to both rules, on which they must **agree**. Every shape the corpus
#: holds, plus the one it does not: 174 of 187 real MAFs are all-prefixed, 9 all-bare, and
#: **none** mixes the two — so the mixed case is where an unwritten disagreement would
#: hide, which is why it is here.
#:
#: A column that says nothing is deliberately not in this list; it is the declared
#: exception, and has its own test below.
COLUMNS = [
    ["1", "2", "X"],
    ["chr1", "chr2", "chrX"],
    ["chr1", "2"],
    [1, 2, 3],
    ["MT"],
    ["chrY"],
]

#: The one column shape on which the two rules part company, kept as data so the exception
#: is a list of one rather than a sentence in a docstring.
SAYS_NOTHING = [None, float("nan"), ""]


def _pipeline_rule() -> "callable":
    """The pipeline's own test, located in its source and transcribed.

    Located rather than assumed: if the ``if`` this guard is about is gone from
    ``bin/``, the rule has moved and the app's copy needs re-deriving by a human — so
    this fails loudly instead of silently comparing against a rule nobody applies.
    """
    tree = ast.parse(ANNOTATION_SCRIPT.read_text(encoding="utf-8"))
    conditions = [
        ast.unparse(node.test)
        for node in ast.walk(tree)
        if isinstance(node, ast.If)
        and "Chromosome" in ast.unparse(node.test)
        and "chr" in ast.unparse(node.test)
    ]
    assert conditions, (
        f"no chromosome-prefix rule found in {ANNOTATION_SCRIPT.relative_to(REPO_ROOT)}. "
        "The pipeline's rule has moved or gone; config/chromosome_spelling.py holds a "
        "copy of it and needs re-deriving against wherever it now lives."
    )
    assert len(conditions) == 1, (
        f"the pipeline now tests the chromosome column in {len(conditions)} places "
        f"({conditions}) — the app copies one rule and there is no longer one rule"
    )
    assert conditions[0] == "not any(['chr' in str(x) for x in maf['Chromosome'].values])", (
        f"the pipeline's chromosome rule now reads {conditions[0]!r}. The app's copy in "
        "config/chromosome_spelling.py was derived from the previous wording; re-derive "
        "it, and update the quotation in that module's docstring, before relaxing this."
    )

    def decides_to_prefix(values):
        return not any("chr" in str(x) for x in values)

    return decides_to_prefix


@needs_annotation_script
@pytest.mark.parametrize("values", COLUMNS, ids=lambda v: "-".join(str(x) for x in v))
def test_the_app_and_the_pipeline_agree_on_every_column_shape(values):
    """The app's answer is the pipeline's answer, for the same column.

    Mutation: make the app's rule per-row rather than whole-column and ``chr1, 2``
    disagrees; widen its substring test to ``startswith`` and nothing here moves, which
    is correct — the two agree on every value either has ever seen.
    """
    from config.chromosome_spelling import CHROMOSOME, normalise_chromosome_spelling

    pipeline_says = _pipeline_rule()(values)

    frame = pd.DataFrame({CHROMOSOME: list(values)})
    app_says = normalise_chromosome_spelling(frame)

    assert app_says is pipeline_says, (
        f"for {values!r} the app says {'prefix' if app_says else 'leave it'} and the "
        f"pipeline says {'prefix' if pipeline_says else 'leave it'} — a MAF would be "
        "spelled one way in the report it came from and another way on screen"
    )


@needs_annotation_script
def test_the_one_place_the_app_declines_to_follow_the_pipeline():
    """The declared exception, asserted rather than only described.

    On a column that says nothing, the pipeline's expression writes ``chrnan`` into every
    cell and the app writes nothing at all. Both halves are pinned: that the two really do
    disagree here — so the exception is not quietly obsolete, the failure mode a declared
    exception has — and that the app's side of the disagreement is the one the display can
    read, which is the whole reason for it.

    Mutation: make the app follow the pipeline exactly and the first assertion fails; drop
    the app's blank guard for some *other* reason and the second one names what broke.
    """
    from config.chromosome_spelling import CHROMOSOME, normalise_chromosome_spelling
    from config.missing_values import says_nothing

    assert _pipeline_rule()(SAYS_NOTHING) is True, (
        "the pipeline would no longer prefix a column that says nothing, so the app has "
        "nothing left to decline — delete this exception rather than relaxing it"
    )

    frame = pd.DataFrame({CHROMOSOME: list(SAYS_NOTHING)})
    assert normalise_chromosome_spelling(frame) is False
    assert all(says_nothing(value) for value in frame[CHROMOSOME]), (
        "a cell that said nothing now says something — the detail panel draws this field "
        f"unconditionally, and would print it at a clinician: {list(frame[CHROMOSOME])!r}"
    )


@needs_annotation_script
def test_the_app_writes_what_the_pipeline_would_have_written():
    """Agreeing on *whether* is not agreeing on *what*.

    The pipeline's assignment is ``"chr" + column.astype(str)``; this asserts the values
    the app leaves behind are the ones that expression produces. Mutation: prefix with
    ``Chr`` or append instead of prepend, and the decision test above still passes while
    this one fails.
    """
    from config.chromosome_spelling import CHROMOSOME, normalise_chromosome_spelling

    values = ["1", "2", "X"]
    frame = pd.DataFrame({CHROMOSOME: list(values)})
    normalise_chromosome_spelling(frame)

    as_the_pipeline_writes_it = "chr" + pd.Series(values).astype(str)
    assert list(frame[CHROMOSOME]) == list(as_the_pipeline_writes_it)
