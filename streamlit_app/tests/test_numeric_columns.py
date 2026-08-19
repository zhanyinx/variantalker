"""Refusing MAF shapes the pipeline crashes on (issue #38).

Three things are asserted here, and they are asserted in this order because each one is
worthless without the one before it:

1. **The derivation.** Which columns the pipeline reads as *numbers* is not written down
   anywhere in the app — it is read out of the vendored source by ``ast`` at import. So
   the first thing to check is that the reader reads correctly, and the only way to check
   a derivation rather than its current answer is to mutate the input and watch the answer
   follow. A test that only pinned ``("t_alt_count", "t_ref_count", "tumor_f")`` would pass
   just as well against a hard-coded tuple, which is the thing this ticket refuses to ship.

2. **The biconditional.** The app must refuse **if and only if** the pipeline's own code
   raises. One direction keeps the app from being fussier than the pipeline — refusing a
   file the pipeline reports on happily; the other keeps it from being laxer — returning a
   confident verdict on a file the pipeline crashes on, which is the worse of the two and
   is what ``somatic_dot_numeric.maf`` used to catch. Both directions are checked against
   the vendored functions themselves, so the oracle is the pipeline's code rather than a
   restatement of what it is believed to do.

3. **The carve-out.** Population-frequency columns are *coerced*, never refused. Every
   fixture in the parity set carries ``.`` in at least one of them — this is the pipeline's
   own blank convention for a column it never reads as a number — and a validator that
   annexed them would refuse all seven, including the six the pipeline handles perfectly.

The oracle for (2) is the **vendored code, called in process**. Not ``bin/``: the app has
to be assertable on a checkout that has no pipeline, and since issue #33 the two are the
same code anyway. ``tests/parity/`` still runs the real subprocess where it can.

Nothing here is skipped without ``bin/`` — that is the point of the file. Two groups are
marked ``slow`` so ``make test-fast`` can drop them: the injected biconditional cases,
which re-read a real reference MAF per case, and the packaged-build probes, which copy
three packages and spawn an interpreter. Everything else runs in milliseconds.
"""

from __future__ import annotations

import ast
import json
import os
import re
import shutil
import subprocess
import sys
import textwrap
from pathlib import Path

import pandas as pd
import pytest

STREAMLIT_APP = Path(__file__).resolve().parent.parent
if str(STREAMLIT_APP) not in sys.path:
    sys.path.insert(0, str(STREAMLIT_APP))

from config.pipeline_params import pipeline_params  # noqa: E402
from filters.numeric_columns import (  # noqa: E402
    NUMERIC_COLUMNS,
    UnreadableNumericColumns,
    derive_numeric_columns,
    unreadable_values,
    vendored_source,
)
from filters.variant_filters import FREQUENCY_COLUMNS, apply_filters  # noqa: E402
from vendor.pipeline_filters import common_filters as pipeline_common_filters  # noqa: E402
from vendor.pipeline_filters import germline_filters as pipeline_germline_filters  # noqa: E402
from vendor.pipeline_filters import somatic_filters as pipeline_somatic_filters  # noqa: E402
from vendor.pipeline_utils import read_maf  # noqa: E402

FIXTURE_DIR = STREAMLIT_APP / "tests" / "fixtures" / "parity"
MANIFEST = json.loads((FIXTURE_DIR / "MANIFEST.json").read_text())

#: Every committed MAF fixture with the arm it must be filtered under. Read from the
#: manifest rather than listed, so a fixture added to the set is covered by the
#: biconditional without anyone remembering to add it here.
FIXTURES = sorted(
    (name, record["arm"])
    for name, record in MANIFEST["fixtures"].items()
    if name.endswith(".maf")
)

#: Tokens that make a column unreadable. Each is a real shape: ``.`` is MAF's own blank
#: marker for text columns and the one ``somatic_dot_numeric.maf`` carries; a lone space
#: is what a hand-edited TSV grows; ``1,5`` is a European decimal comma; ``>100`` is a
#: capped depth some callers emit as free text.
POISON = [".", " ", "1,5", ">100", "n/a "]

#: Tokens that read cleanly as missing and leave the column numeric. The distinction is
#: the whole subtlety of the ticket: an empty cell is *not* the trigger, and refusing on
#: one would refuse the synthetic fixtures, whose whole purpose is the NaN ruling.
CLEAN = ["", "NA", "NaN", "N/A", "null"]


def _params(arm: str, **overrides) -> dict:
    params = pipeline_params(arm)
    params.update(overrides)
    return params


def _frame(fixture: str) -> pd.DataFrame:
    return read_maf(str(FIXTURE_DIR / fixture))


def _read_tsv(text: str, comments: int = 0) -> pd.DataFrame:
    """Parse ``text`` exactly as ``pipeline_utils.read_maf`` parses a file.

    One helper for both of this module's constructors — the hand-built frames of
    :func:`read_maf_like` and the one-cell edits of :func:`_inject`. They differ only in
    where the text comes from, and the part that has to be right is the same for both: the
    reader's own ``header``/``sep``/``low_memory`` arguments, because a frame built any other
    way can carry a dtype ``read_maf`` would never have produced, and dtype is the subject.
    """
    import io

    return pd.read_csv(io.StringIO(text), header=comments, sep="\t", low_memory=False)


# ---------------------------------------------------------------------------
# The oracle: does the pipeline's own code raise on this frame?
# ---------------------------------------------------------------------------


def pipeline_raises(frame: pd.DataFrame, arm: str) -> bool:
    """Whether the vendored filters raise ``TypeError`` on ``frame``.

    The three calls ``apply_filters`` makes, in the order it makes them, with the gene
    clause disabled — a gene list changes which rows survive and never whether a
    comparison is possible. ``KeyError`` is deliberately **not** caught: an absent column
    is issue #39's decision, not this one's, and swallowing it here would let a test pass
    for the wrong reason.
    """
    contract = pipeline_params(arm)
    try:
        pipeline_common_filters(
            frame, contract["min_depth"], contract["filter_variant_classification"]
        )
        if arm == "somatic":
            pipeline_somatic_filters(
                frame,
                vaf=contract["vaf_threshold"],
                somatic_genes="null",
                cancervar_keep=contract["filter_cancervar"],
                civic_keep=contract["filter_civic"],
                escat_keep=contract["filter_escat"],
                clinvar_keep=contract["filter_clinvar"],
                skip_civic=contract.get("skip_civic", False),
            )
        else:
            pipeline_germline_filters(
                frame,
                vaf=contract["vaf_threshold"],
                germline_genes="null",
                intervar_keep=contract["filter_intervar"],
                renovo_keep=contract["filter_renovo"],
                clinvar_keep=contract["filter_clinvar"],
            )
    except TypeError:
        return True
    return False


def app_refuses(frame: pd.DataFrame, arm: str) -> bool:
    """Whether the entry point refuses ``frame`` with the typed exception.

    A ``TypeError`` escaping here is a biconditional failure in the dangerous direction
    inverted — the validator let a frame through that the vendored code then choked on —
    so it is turned into an explicit failure rather than an error with a bare traceback.
    """
    try:
        apply_filters(frame, _params(arm))
    except UnreadableNumericColumns:
        return True
    except TypeError as exc:  # pragma: no cover - only on a broken validator
        pytest.fail(
            f"the validator passed a frame the vendored code then raised on: {exc}"
        )
    return False


# ---------------------------------------------------------------------------
# 1. The derivation
# ---------------------------------------------------------------------------


def test_the_derivation_names_the_columns_the_pipeline_compares_as_numbers():
    """``t_alt_count + t_ref_count >= coverage`` and ``tumor_f > vaf``. Nothing else.

    Pinned so the *answer* is visible in the suite, but pinning it is the weakest thing
    here — the mutation tests below are what establish that this is read rather than
    remembered.
    """
    assert NUMERIC_COLUMNS == ("t_alt_count", "t_ref_count", "tumor_f")


@pytest.mark.parametrize(
    "column",
    [
        "Variant_Classification",
        "CancerVar",
        "ClinVar_VCF_CLNSIG",
        "ESCAT",
        "InterVar",
        "RENOVO_Class",
        "CIViC_Evidence_Level",
        "Hugo_Symbol",
    ],
)
def test_columns_the_pipeline_reads_as_text_are_not_derived(column):
    """``.isin``, ``.apply`` and ``.str.upper`` are not numeric contexts.

    Every one of these columns is read by the vendored filters, and every one of them is
    routinely non-numeric in a real MAF. Deriving any of them would refuse every file in
    the reference set.
    """
    assert column not in NUMERIC_COLUMNS


@pytest.mark.parametrize("column", FREQUENCY_COLUMNS)
def test_frequency_columns_are_not_derived(column):
    """The pipeline has no frequency filter, so no frequency column is read numerically.

    ``compute_keep`` names ``gnomAD_exome_AF`` and its three siblings, but only to test
    membership in ``out.columns`` — a column *set* decision, not a value comparison. A
    derivation that walked every mention rather than every numeric context would pick them
    up and refuse all seven fixtures.
    """
    assert column not in NUMERIC_COLUMNS


#: Mutations of a minimal source, each stating one rule of the derivation. Written as
#: whole functions rather than as edits to the vendored text so that what is being claimed
#: is legible: the frame parameter is found by its ``pd.DataFrame`` annotation, and a
#: column counts only when its own values reach an arithmetic or ordering operator.
_MUTATIONS = [
    pytest.param(
        """
        def f(maf: pd.DataFrame, x):
            return maf["A"] > x
        """,
        ("A",),
        id="ordering-compare",
    ),
    pytest.param(
        """
        def f(maf: pd.DataFrame, x):
            return (maf["A"] + maf["B"]) >= x
        """,
        ("A", "B"),
        id="sum-inside-a-compare",
    ),
    pytest.param(
        """
        def f(maf: pd.DataFrame, x):
            return -maf["A"] / 2 < x
        """,
        ("A",),
        id="unary-and-division",
    ),
    pytest.param(
        """
        def f(maf: pd.DataFrame, x):
            return maf["A"].isin(x)
        """,
        (),
        id="isin-is-not-numeric",
    ),
    pytest.param(
        """
        def f(maf: pd.DataFrame, x):
            return maf["A"].apply(lambda v: v > 1)
        """,
        (),
        id="apply-hides-the-column-behind-a-call",
    ),
    pytest.param(
        """
        def f(maf: pd.DataFrame, x):
            return maf["A"].str.len() > x
        """,
        (),
        id="a-derived-length-is-not-the-columns-own-values",
    ),
    pytest.param(
        """
        def f(maf: pd.DataFrame, x):
            return maf["A"] == x
        """,
        (),
        id="equality-is-not-an-ordering",
    ),
    pytest.param(
        """
        def f(frame, x):
            return frame["A"] > x
        """,
        (),
        id="an-unannotated-parameter-is-not-known-to-be-a-frame",
    ),
    pytest.param(
        """
        def f(maf: pd.DataFrame, other: pd.DataFrame, x):
            return (maf["A"] > x) & (other["B"] > x)
        """,
        ("A", "B"),
        id="every-annotated-frame-counts",
    ),
    pytest.param(
        """
        def f(maf: pd.DataFrame, x):
            column = "A"
            return maf[column] > x
        """,
        (),
        id="only-literal-column-names-are-readable",
    ),
]


@pytest.mark.parametrize("source,expected", _MUTATIONS)
def test_the_derivation_follows_the_source_under_mutation(source, expected):
    assert derive_numeric_columns(textwrap.dedent(source)) == expected


def test_the_derivation_follows_a_renamed_column_in_the_vendored_source():
    """Rename the depth column throughout the real vendored source; the answer follows.

    The strongest form of the claim, because the input is the actual file rather than a
    miniature: if this ever came back ``t_alt_count`` it would mean the tuple is coming
    from somewhere other than the source.
    """
    renamed, substitutions = re.subn(r'"t_alt_count"', '"depth_alt"', vendored_source())
    assert substitutions, "t_alt_count is not in the vendored source; nothing was renamed"

    assert derive_numeric_columns(renamed) == ("depth_alt", "t_ref_count", "tumor_f")


def test_the_derivation_picks_up_a_numeric_gate_added_to_the_vendored_source():
    """A gate the pipeline does not have today, spliced into the source it does have.

    The splice is matched by regex rather than by an exact literal including its
    indentation. ``bin/`` is frozen for the parity effort but not embalmed, and
    ``test_vendor_drift.py``'s standing rule is that nothing may go red over formatting
    churn — an exact-text ``.replace`` would silently become a no-op if the line were
    re-wrapped, and this test would then fail for a reason that is not about the derivation
    at all. The assertion below that the substitution actually happened is what keeps the
    looser match honest.
    """
    source, substitutions = re.subn(
        r'filter_vaf\s*=\s*maf\["tumor_f"\]\s*>\s*vaf',
        'filter_vaf = (maf["tumor_f"] > vaf) & (maf["DP"] >= 8)',
        vendored_source(),
        count=1,
    )
    assert substitutions == 1, (
        "the VAF gate was not found in the vendored source, so nothing was mutated and "
        "this test would pass without testing anything"
    )

    assert "DP" in derive_numeric_columns(source)
    assert "DP" not in NUMERIC_COLUMNS


def test_the_derivation_drops_a_column_the_pipeline_stops_comparing():
    """Turn every ``tumor_f`` comparison into a membership test; the column drops out.

    Both comparisons have to go: the VAF gate and the ``> -1`` all-true gene default. That
    is why the substitution count is asserted — mutating only one of them would leave
    ``tumor_f`` derived for the *right* reason and the test would look like it had shown
    something it had not.
    """
    source, substitutions = re.subn(
        r'maf\["tumor_f"\]\s*>\s*(vaf|-1)',
        r'maf["tumor_f"].isin([\1])',
        vendored_source(),
    )
    assert substitutions >= 2, f"expected both tumor_f gates, mutated {substitutions}"

    assert "tumor_f" not in derive_numeric_columns(source)


def test_the_derived_columns_are_all_read_by_the_vendored_source():
    """No derived column is an artefact: each one is subscripted in the source it came
    from. A cheap sanity net under the parser, independent of the rules above."""
    source = vendored_source()
    for column in NUMERIC_COLUMNS:
        assert f'"{column}"' in source


def test_the_vendored_source_is_the_module_the_app_actually_calls():
    """The source read for the derivation must parse to the module being imported.

    Guards the one way this could be quietly wrong: reading some *other* copy of the
    filter code — ``bin/``, a stale file, an installed egg — and deriving a contract the
    running app does not obey.
    """
    from vendor import pipeline_filters

    parsed = {
        node.name
        for node in ast.parse(vendored_source()).body
        if isinstance(node, ast.FunctionDef)
    }
    assert "common_filters" in parsed
    assert parsed <= set(dir(pipeline_filters))


# ---------------------------------------------------------------------------
# 2. The validator, unit by unit
# ---------------------------------------------------------------------------


def _minimal(**columns) -> pd.DataFrame:
    """A frame carrying only the derived columns, at readable values unless overridden."""
    base = {"t_alt_count": [30, 40], "t_ref_count": [30, 40], "tumor_f": [0.5, 0.6]}
    base.update(columns)
    return pd.DataFrame(base)


def test_a_clean_frame_is_not_refused():
    assert unreadable_values(_minimal()) == {}


def test_an_absent_column_is_skipped():
    """The validator does not annex the missing-column decision, which is issue #39's.

    Dropping a derived column entirely must leave the validator silent — it has nothing to
    say about a column that is not there, and saying something would pre-empt a ruling this
    ticket has no measurement for.
    """
    frame = _minimal().drop(columns=["tumor_f"])
    assert unreadable_values(frame) == {}


@pytest.mark.parametrize("token", CLEAN)
def test_tokens_that_read_as_missing_are_not_refused(token):
    """``''`` and ``NA`` leave the column numeric with a NaN in it. The pipeline then
    drops the row on the comparison, which is a verdict — not a crash — and the app must
    reach the same one rather than refusing the file."""
    frame = read_maf_like({"t_alt_count": ["30", token], "tumor_f": ["0.5", "0.6"]})
    assert unreadable_values(frame) == {}


@pytest.mark.parametrize("token", POISON)
def test_a_single_unreadable_token_refuses_the_column(token):
    """One bad cell refuses its column, and the message names that exact cell.

    The offender is compared **exactly**, not by substring: the earlier spelling was
    ``token.strip() in " ".join(...)``, which for the lone-space token asks whether ``""``
    appears in a string — true of every string, so the one token whose handling is least
    obvious was the one token going unchecked.
    """
    frame = read_maf_like({"t_alt_count": ["30", token], "tumor_f": ["0.5", "0.6"]})
    assert unreadable_values(frame) == {"t_alt_count": {token: 1}}


def read_maf_like(columns: dict) -> pd.DataFrame:
    """A frame with the dtypes ``read_maf`` would have inferred from these cells.

    Round-tripping through the reader rather than constructing the frame directly, because
    the dtype is the whole subject: ``pd.DataFrame({"a": ["30", "."]})`` and a TSV holding
    ``30`` and ``.`` produce the same object column, but ``["30", ""]`` and a TSV holding
    ``30`` and a blank do not — the reader turns the blank into ``NaN`` and keeps the column
    numeric, and that difference is exactly what decides refusal.
    """
    frame = pd.DataFrame(columns)
    for missing in NUMERIC_COLUMNS:
        if missing not in frame:
            frame[missing] = ["30"] * len(frame)
    return _read_tsv(frame.to_csv(index=False, sep="\t"))


def test_every_offending_column_and_value_is_named_in_one_message():
    """One message, not one per column and not a sample of them.

    A user handed "column X is unreadable" fixes X, reloads, and is told about Y. The
    columns are unreadable *together* — the file is one artefact — so the refusal names
    everything at once.
    """
    frame = read_maf_like(
        {
            "t_alt_count": ["30", "."],
            "t_ref_count": [">100", "40"],
            "tumor_f": ["0.5", "1,5"],
        }
    )
    with pytest.raises(UnreadableNumericColumns) as caught:
        apply_filters(frame, _params("somatic"))

    message = str(caught.value)
    for column in ("t_alt_count", "t_ref_count", "tumor_f"):
        assert column in message
    for value in (".", ">100", "1,5"):
        assert repr(value) in message


def test_the_message_does_not_name_the_numbers_the_bad_cell_dragged_down_with_it():
    """Only what the user can act on.

    One ``.`` leaves every other cell in the column a *string*, so "every value that is not
    a number" is the whole column — 150,000 distinct numeric strings on a real MAF. The
    message names the tokens that are not numbers in any reading and leaves the rest out.
    """
    frame = read_maf_like({"t_alt_count": ["30", "."], "tumor_f": ["0.5", "0.6"]})
    with pytest.raises(UnreadableNumericColumns) as caught:
        apply_filters(frame, _params("somatic"))

    message = str(caught.value)
    assert repr(".") in message
    assert repr("30") not in message


def test_a_readable_token_in_a_poisoned_column_is_not_reported_as_an_offender():
    frame = read_maf_like({"tumor_f": ["0.5", "."]})
    assert unreadable_values(frame) == {"tumor_f": {".": 1}}


def test_the_refusal_carries_the_offenders_for_a_caller_that_wants_them():
    frame = read_maf_like({"t_alt_count": ["30", "."], "tumor_f": ["0.5", "."]})
    with pytest.raises(UnreadableNumericColumns) as caught:
        apply_filters(frame, _params("somatic"))

    offenders = caught.value.offenders
    assert set(offenders) == {"t_alt_count", "tumor_f"}
    assert offenders["t_alt_count"] == {".": 1}


def test_the_refusal_is_a_typed_exception_not_a_diagnostic():
    """There is no verdict to return, so it cannot travel on the diagnostics channel.

    ``Diagnostics.warnings`` describes a report that exists. A refused MAF produces no
    report, so a warning would have to accompany a frame the app cannot have computed.
    """
    frame = read_maf_like({"tumor_f": ["0.5", "."]})
    with pytest.raises(UnreadableNumericColumns):
        apply_filters(frame, _params("somatic"))
    assert issubclass(UnreadableNumericColumns, ValueError)


def test_a_column_of_python_numbers_is_not_refused_for_its_dtype():
    """Object dtype alone is not the trigger — the pipeline compares these fine.

    A dtype check would be the obvious implementation and would be *fussier than the
    pipeline*: ``pd.Series([1.0, 2.0], dtype=object) > 0.5`` is a perfectly good mask. The
    trigger is a value that is not a number, not a container that has forgotten it holds
    numbers.
    """
    frame = _minimal(tumor_f=pd.Series([0.5, 0.6], dtype=object))
    assert not pd.api.types.is_numeric_dtype(frame["tumor_f"])
    assert unreadable_values(frame) == {}


#: Tokens that a naive "name what ``to_numeric`` cannot read" rule might suppress, because
#: ``pd.to_numeric`` *can* read them. Included to show the suppression is not reachable:
#: ``read_csv`` reads them too, so the column stays numeric and there is nothing to refuse.
COERCIBLE_ODDITIES = ["inf", "-inf", "1e5", "+30", "0030"]


@pytest.mark.parametrize("token", COERCIBLE_ODDITIES + POISON)
def test_no_token_that_poisons_a_column_goes_unnamed(token):
    """The gap between "refuses" and "names a value" is empty for anything a reader produces.

    The message names the values ``pd.to_numeric`` cannot read, and drops the rest — which
    on its face could suppress a real culprit that ``to_numeric`` happens to accept, like
    ``inf`` or ``1e5``. It cannot, and this is where that is established rather than argued:
    a token ``to_numeric`` accepts is a token ``read_csv`` accepted first, so the column
    never became object dtype and there was nothing to refuse. Either the file loads
    normally, or the offending token is named.

    Stated as a property over both lists so it keeps holding if either grows. If some future
    token does land in the gap — refused with an empty offender map out of a real file — this
    goes red and the reporting rule has to be widened.
    """
    frame = read_maf_like({"t_alt_count": ["30", token], "tumor_f": ["0.5", "0.6"]})
    offenders = unreadable_values(frame)

    if pd.api.types.is_numeric_dtype(frame["t_alt_count"]):
        assert offenders == {}, f"{token!r} kept the column numeric but was still refused"
    else:
        assert offenders == {"t_alt_count": {token: 1}}, (
            f"{token!r} poisoned the column but the refusal does not name it: {offenders}"
        )


def test_a_column_of_numeric_strings_is_refused_with_no_value_at_fault():
    """``"5" + 3`` raises just as loudly as ``"." + 3``, so a coercible string is still
    a refusal. Named because the tempting implementation — ``pd.to_numeric(coerce)`` and
    look for new NaNs — gets this one wrong in the dangerous direction.

    No individual cell is wrong here, so the message says that rather than naming the whole
    column back at the user.
    """
    frame = _minimal(tumor_f=pd.Series(["0.5", "0.6"], dtype=object))
    assert unreadable_values(frame) == {"tumor_f": {}}
    assert pipeline_raises(_somatic_shaped(frame), "somatic")

    with pytest.raises(UnreadableNumericColumns) as caught:
        apply_filters(_somatic_shaped(frame), _params("somatic"))
    assert "no single value is at fault" in str(caught.value)


def _somatic_shaped(frame: pd.DataFrame) -> pd.DataFrame:
    """``frame`` plus the text columns the somatic filters read, at neutral values.

    So that :func:`pipeline_raises` can be used as the oracle on a hand-built frame: without
    them the vendored code raises ``KeyError`` on a missing column, which is issue #39's
    subject and not an answer to the question being asked here.
    """
    filled = frame.copy()
    for column, value in [
        ("Variant_Classification", "Missense_Mutation"),
        ("CancerVar", "Tier_III_Uncertain"),
        ("ClinVar_VCF_CLNSIG", "Uncertain_significance"),
        ("ESCAT", "IIIA"),
        ("Hugo_Symbol", "TP53"),
    ]:
        if column not in filled:
            filled[column] = value
    return filled


def test_the_validator_runs_before_anything_else_in_the_entry_point():
    """Refusal must not depend on the parameters, because the pipeline's raise does not.

    Asserted with a parameter dict that would fail its own validation later on: if the
    refusal came after parameter resolution the caller would be told about ``sample_type``
    and never see the values that actually stopped the file.
    """
    frame = read_maf_like({"tumor_f": ["0.5", "."]})
    with pytest.raises(UnreadableNumericColumns):
        apply_filters(frame, {"sample_type": "not-an-arm"})


# ---------------------------------------------------------------------------
# 3. The biconditional
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("fixture,arm", FIXTURES, ids=[f for f, _ in FIXTURES])
def test_refusal_matches_the_pipeline_on_every_fixture(fixture, arm):
    """Refuse if and only if the vendored code raises — over the committed set.

    Seven files, six of which the pipeline filters happily and one (``somatic_dot_numeric``)
    it raises on. Both arms are represented.
    """
    frame = _frame(fixture)
    assert app_refuses(frame, arm) == pipeline_raises(frame, arm)


def _injected_cases():
    """Every derived column poisoned in turn, on both arms, against clean controls."""
    cases = []
    for fixture, arm in [
        ("somatic_reference.maf", "somatic"),
        ("germline_reference.maf", "germline"),
    ]:
        for column in NUMERIC_COLUMNS:
            for token in POISON + CLEAN:
                cases.append(
                    pytest.param(
                        fixture,
                        arm,
                        column,
                        token,
                        id=f"{arm}-{column}-{token.strip() or 'blank'!r}",
                    )
                )
    return cases


@pytest.mark.slow
@pytest.mark.parametrize("fixture,arm,column,token", _injected_cases())
def test_refusal_matches_the_pipeline_on_injected_values(fixture, arm, column, token):
    """The biconditional over shapes the fixture set does not contain.

    One cell of a real reference MAF is overwritten and the file re-read, so the dtype is
    whatever ``read_maf`` would really have inferred — the injection cannot accidentally
    prove the point by constructing a dtype the reader never produces.
    """
    frame = _inject(FIXTURE_DIR / fixture, column, token)
    assert app_refuses(frame, arm) == pipeline_raises(frame, arm)


def _inject(path: Path, column: str, token: str, row: int = 0) -> pd.DataFrame:
    """``path`` with one cell replaced, re-read through the pipeline's own reader."""
    lines = path.read_text().splitlines()
    header_index = next(i for i, line in enumerate(lines) if not line.startswith("#"))
    position = lines[header_index].split("\t").index(column)

    target = header_index + 1 + row
    cells = lines[target].split("\t")
    cells[position] = token
    lines[target] = "\t".join(cells)

    return _read_tsv("\n".join(lines) + "\n", comments=header_index)


def test_one_bad_value_poisons_the_whole_column_so_there_is_no_per_row_rescue():
    """Dropping the offending rows does not make the column readable again.

    The reason the refusal is whole-file rather than a row filter. ``read_maf`` decides the
    dtype once, when it reads the file; deleting the row afterwards leaves an object column
    of Python strings behind, which the pipeline still cannot compare.
    """
    frame = _inject(FIXTURE_DIR / "somatic_reference.maf", "tumor_f", ".")
    survivors = frame[frame["tumor_f"] != "."]

    assert len(survivors) == len(frame) - 1
    assert not pd.api.types.is_numeric_dtype(survivors["tumor_f"])
    assert pipeline_raises(survivors, "somatic")
    assert app_refuses(survivors, "somatic")


# ---------------------------------------------------------------------------
# 4. The frequency carve-out
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("fixture", [f for f, _ in FIXTURES])
def test_frequency_columns_are_never_refused(fixture):
    """No fixture is refused *for* a frequency column — and every one has cause to be.

    Each committed MAF carries ``.`` in at least one frequency column. That is the shape
    that refuses a depth column, so this is not a vacuous check: it is the same token in a
    column the pipeline does not read as a number, and the answer has to be different.
    """
    frame = _frame(fixture)
    present = [c for c in FREQUENCY_COLUMNS if c in frame.columns]
    dotted = [c for c in present if (frame[c].astype(str) == ".").any()]
    assert dotted, f"{fixture} carries no '.' in a frequency column — check is vacuous"

    assert set(unreadable_values(frame)) & set(present) == set()


@pytest.mark.parametrize("fixture", [f for f, _ in FIXTURES])
def test_the_frequency_coercion_preserves_every_number_on_every_fixture(fixture):
    """The coercion is a no-op on the values that are really there.

    ``pd.to_numeric(errors="coerce")`` exists in the frequency mask to turn ``.`` into
    *missing* — "not in the panel" — and the risk of a silent coercion is that it also eats
    a value that meant something. Two claims per fixture, per column, and they are the two
    halves of "no-op": every value the coercion keeps is exactly the number its text said,
    and every value it drops was a missing-marker and not data.

    Both halves are guarded against vacuity at the end, because the obvious spelling of the
    first one is vacuous here and silently so. A ``.`` makes the whole column object dtype
    holding *strings*, so a check phrased as "values that are already ``float`` survive
    unchanged" matches nothing at all on any of these seven files and passes by having
    nothing to say.

    The guard is per **fixture**, not per column: some columns really are empty end to end
    on some fixtures — ``gnomAD_exome_AF`` is blank in every row of ``somatic_civic.maf`` —
    and an all-missing column has no value to preserve and nothing to prove.
    """
    frame = _frame(fixture)
    kept_total = lost_total = 0

    for column in FREQUENCY_COLUMNS:
        if column not in frame.columns:
            continue
        raw = frame[column]
        coerced = pd.to_numeric(raw, errors="coerce")

        kept = coerced.notna()
        assert [float(text) for text in raw[kept].astype(str)] == coerced[kept].tolist(), (
            f"{fixture}:{column}: the coercion changed a value it kept"
        )

        lost = coerced.isna() & raw.notna()
        assert set(raw[lost].astype(str)) <= {".", "", "NA"}, (
            f"{fixture}:{column}: the coercion discarded "
            f"{sorted(set(raw[lost].astype(str)))[:5]}, which are not missing-markers"
        )

        kept_total += int(kept.sum())
        lost_total += int(lost.sum())

    assert kept_total, f"{fixture}: no frequency value survived the coercion anywhere"
    assert lost_total, (
        f"{fixture}: the coercion dropped nothing anywhere, so it was never exercised — "
        "this fixture is supposed to carry '.' in a frequency column"
    )


@pytest.mark.parametrize("fixture,arm", FIXTURES, ids=[f for f, _ in FIXTURES])
def test_moving_the_frequency_slider_does_not_crash_on_any_fixture(fixture, arm):
    """The app-only extra, exercised off its default.

    The frequency layer is neutral at 1.0 and every other test runs it there, so a crash in
    it would be latent until the first user moved the slider. Any fixture the pipeline
    itself raises on is refused first, which is the correct answer at any threshold.
    """
    frame = _frame(fixture)
    if pipeline_raises(frame, arm):
        with pytest.raises(UnreadableNumericColumns):
            apply_filters(frame, _params(arm, max_freq_population=0.01))
        return

    labelled, diagnostics = apply_filters(frame, _params(arm, max_freq_population=0.01))
    assert len(labelled) == len(frame)
    assert sum(diagnostics.cells().values()) == len(frame)


# ---------------------------------------------------------------------------
# 5. What the parity harness counts as agreement
# ---------------------------------------------------------------------------

#: ``(pipeline_error, refused_columns, in_parity)``. The rule
#: ``harness.ParityResult.errors_in_parity`` applies, stated as a table.
#:
#: Tested here rather than in ``tests/parity/test_parity.py`` on purpose: that module skips
#: entirely without ``bin/``, and this is pure logic about the app's refusal that needs no
#: pipeline at all. The repo already draws that line — ``test_loader_premise.py`` and
#: ``test_baseline_integrity.py`` carry no ``skipif`` for the same reason.
_ERROR_PARITY = [
    pytest.param(None, None, True, id="neither-side-failed"),
    pytest.param("TypeError: '>=' not supported", None, False, id="pipeline-only"),
    pytest.param(
        "TypeError: '>=' not supported",
        ["t_alt_count"],
        True,
        id="refusal-against-the-matching-TypeError",
    ),
    pytest.param(
        "KeyError: 'InterVar'",
        ["t_alt_count"],
        False,
        id="refusal-against-an-absent-column-KeyError",
    ),
    pytest.param(None, ["t_alt_count"], False, id="app-refused-where-pipeline-did-not"),
]


@pytest.mark.parametrize("pipeline_error,refused,expected", _ERROR_PARITY)
def test_only_a_refusal_against_a_type_error_counts_as_agreement(
    pipeline_error, refused, expected
):
    """The refusal stands in for the pipeline's ``TypeError`` and for nothing else.

    The ``KeyError`` row is the one worth having. An absent filter-input column makes the
    vendored code raise ``KeyError``, which is issue #39's ruling to make — and this
    validator deliberately *skips* absent columns rather than annexing it. So a refusal
    sitting opposite a ``KeyError`` means the two sides stopped the file for unrelated
    reasons, and counting that as parity would let a real divergence read as agreement.
    """
    from tests.parity.harness import ParityResult

    result = ParityResult(
        case="synthetic", arm="somatic", fixture="none.maf", rows=1,
        pipeline_error=pipeline_error, app_refused_columns=refused,
        app_error=None if refused is None else "UnreadableNumericColumns",
    )
    assert result.errors_in_parity is expected


# ---------------------------------------------------------------------------
# 6. The packaged build
# ---------------------------------------------------------------------------


#: The packages the app needs to reach the derivation. Notably **not** ``bin/``, which no
#: installer ships — see ``test_parity.py::test_harness_is_excluded_from_packaged_builds``.
SHIPPED_PACKAGES = ("vendor", "filters", "config")

PACKAGED_PROBE = """
import sys
sys.path.insert(0, sys.argv[1])
from filters.numeric_columns import NUMERIC_COLUMNS
print(",".join(NUMERIC_COLUMNS))
"""


def _shipped_tree(root: Path) -> Path:
    """The app's packages, copied without ``bin/``, ``tests/`` or any bytecode."""
    app = root / "app"
    app.mkdir()
    for package in SHIPPED_PACKAGES:
        shutil.copytree(
            STREAMLIT_APP / package,
            app / package,
            ignore=shutil.ignore_patterns("__pycache__", "*.pyc"),
        )
    return app


def _derive_in(tmp_path: Path, import_root: str, cwd: Path) -> str:
    """Run the derivation in a fresh interpreter with ``import_root`` on ``sys.path``."""
    probe = tmp_path / "probe.py"
    probe.write_text(PACKAGED_PROBE)
    proc = subprocess.run(
        [sys.executable, "-B", str(probe), import_root],
        cwd=cwd,
        capture_output=True,
        text=True,
        # Empty PYTHONPATH so the probe cannot fall back to this checkout and answer from
        # the very tree the test is trying to exclude.
        env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1", "PYTHONPATH": ""},
    )
    assert proc.returncode == 0, proc.stderr
    return proc.stdout.strip()


@pytest.mark.slow
@pytest.mark.integration
def test_the_derivation_works_from_a_plain_source_tree(tmp_path):
    """No ``bin/``, no ``tests/``, no ``__pycache__`` — what the .dmg and .exe actually ship.

    The failure this rules out is the derivation reading ``bin/filter_variants.py``. That
    file is the source of truth for the *drift guard* and is the obvious thing to parse, but
    no installer ships it, so an app built that way would fail to start — or, worse, refuse
    nothing at all. ``-B`` and the pruned copy also mean nothing here can be answered from
    cached bytecode.
    """
    app = _shipped_tree(tmp_path)
    assert _derive_in(tmp_path, str(app), app) == ",".join(NUMERIC_COLUMNS)


@pytest.mark.slow
@pytest.mark.integration
def test_the_derivation_works_when_the_package_is_imported_from_a_zip(tmp_path):
    """The other half of "works in a packaged build": source behind a zip importer.

    A packaging shape with no ``.py`` on disk to open, which is the one the plain-tree case
    cannot speak to. Measured, so the claim is not stronger than the evidence: **both**
    routes in :func:`~filters.numeric_columns.vendored_source` work here —
    ``zipimporter.get_source`` and ``inspect.getsource`` each return the same 8,525
    characters — so this does not discriminate between them, and the loader is asked first
    for the reason given there rather than because ``inspect`` fails under a zip.

    What it does establish is the thing the acceptance criterion asks for: the derivation
    survives being imported from an archive, which is what would have gone wrong had it
    reached for ``__file__`` and joined a path onto it.
    """
    app = _shipped_tree(tmp_path)
    archive = shutil.make_archive(str(tmp_path / "bundle"), "zip", root_dir=str(app))

    # cwd deliberately outside the tree, so the only importable copy is the archive.
    assert _derive_in(tmp_path, archive, tmp_path) == ",".join(NUMERIC_COLUMNS)
