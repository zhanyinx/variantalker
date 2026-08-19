"""Drift guard for the predictor cutoffs the panel prints (issue #190).

``components/predictor_context.py`` invents no thresholds. It holds a **copy** of the ones two
vendored classifiers already apply::

    resources/CancerVar/CancerVar.py:check_PreP   sift_cutoff = 0.05, cutoff_conserv = 2
    resources/InterVar/Intervar.py:check_PP3      metasvm_cutoff = 0.0, cutoff_conserv = 2,
                                                 dbscSNV_cutoff = 0.6

That is the same relationship ``config/chromosome_spelling.py`` has to
``bin/add_guidelines_escat_to_funcotator.py``, and it fails the same way: silently. The section
prints *"damaging — below CancerVar's 0.05 cutoff"* beside a score. If the vendored script's cutoff
moves and the app's copy does not, the panel keeps stating a threshold nobody applies, on a screen
a clinician reads to check a verdict. Nothing raises; the sentence is simply false.

Map #184 rule 2 makes ``resources/`` read-only for this effort, so the app can only ever hold a copy
— which is exactly why printing the number needs this guard. Issue #190's decision was that the
threshold *is* shown, on the condition that a test reads it back out of the authority.

**Why this compares the assignments textually**, where
``tests/test_chromosome_rule_contract.py`` deliberately compares behaviour: there, the two sides are
written in different idioms and no textual equality could hold. Here the authority *is* a literal —
five ``name = number`` assignments inside two named functions — so the faithful comparison is of
the numbers themselves. Parsed from the AST rather than grepped, so a cutoff mentioned in a comment
or in a neighbouring function cannot satisfy the guard: this repo has already had a text guard
satisfied by its own file header.

**Why the function is named as well as the variable.** ``cutoff_conserv = 2`` appears in four
functions across the two scripts, and ``check_BP4``'s copy is the same number by coincidence rather
than by construction. The guard asks each function for its own.
"""

import ast
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from components.predictor_context import (  # noqa: E402
    CANCERVAR_GERP_CONSERVED_FROM,
    CANCERVAR_SIFT_DAMAGING_BELOW,
    INTERVAR_DBSCSNV_SPLICE_ALTERING_ABOVE,
    INTERVAR_GERP_CONSERVED_ABOVE,
    INTERVAR_METASVM_CUTOFF,
)

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CANCERVAR = os.path.join(REPO, "resources", "CancerVar", "CancerVar.py")
INTERVAR = os.path.join(REPO, "resources", "InterVar", "Intervar.py")

#: ``(script, function, variable, the app's copy)``. One row per number the panel prints.
CONTRACT = (
    (CANCERVAR, "check_PreP", "sift_cutoff", CANCERVAR_SIFT_DAMAGING_BELOW),
    (CANCERVAR, "check_PreP", "cutoff_conserv", CANCERVAR_GERP_CONSERVED_FROM),
    (INTERVAR, "check_PP3", "metasvm_cutoff", INTERVAR_METASVM_CUTOFF),
    (INTERVAR, "check_PP3", "cutoff_conserv", INTERVAR_GERP_CONSERVED_ABOVE),
    (INTERVAR, "check_PP3", "dbscSNV_cutoff", INTERVAR_DBSCSNV_SPLICE_ALTERING_ABOVE),
)

#: The same three numbers as ``check_PP3`` reads them, checked in ``check_BP4`` too. BP4 is the
#: reversed comparison over the same cutoffs, and the panel's sentences name it, so a script that
#: moved one and not the other would leave half the section true.
BP4_CONTRACT = (
    (INTERVAR, "check_BP4", "metasvm_cutoff", INTERVAR_METASVM_CUTOFF),
    (INTERVAR, "check_BP4", "cutoff_conserv", INTERVAR_GERP_CONSERVED_ABOVE),
    (INTERVAR, "check_BP4", "dbscSNV_cutoff", INTERVAR_DBSCSNV_SPLICE_ALTERING_ABOVE),
)

pytestmark = pytest.mark.skipif(
    not (os.path.exists(CANCERVAR) and os.path.exists(INTERVAR)),
    reason="the pipeline's vendored classifiers under resources/ are not present",
)


def _assignments(path: str, function: str) -> dict:
    """Every ``name = <number>`` assigned directly in one function, as a dict.

    Read from the AST so that a number in a comment, a docstring, or another function cannot
    answer for this one. Unary minus is handled because a cutoff may legitimately be negative.
    """
    with open(path, encoding="utf-8") as handle:
        tree = ast.parse(handle.read())

    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name == function:
            found = {}
            for statement in ast.walk(node):
                if not isinstance(statement, ast.Assign):
                    continue
                value = statement.value
                if isinstance(value, ast.UnaryOp) and isinstance(value.op, ast.USub):
                    if not isinstance(value.operand, ast.Constant):
                        continue
                    number = -value.operand.value
                elif isinstance(value, ast.Constant):
                    number = value.value
                else:
                    continue
                if not isinstance(number, (int, float)) or isinstance(number, bool):
                    continue
                for target in statement.targets:
                    if isinstance(target, ast.Name):
                        found[target.id] = number
            return found
    return {}


@pytest.mark.parametrize(
    ("path", "function", "variable", "ours"), CONTRACT + BP4_CONTRACT,
    ids=[f"{f}:{v}" for _, f, v, _ in CONTRACT + BP4_CONTRACT],
)
def test_the_panel_prints_the_cutoff_the_script_applies(path, function, variable, ours):
    """Each printed threshold equals the one its own classifier assigns."""
    assignments = _assignments(path, function)
    assert variable in assignments, (
        f"{os.path.basename(path)}:{function} no longer assigns `{variable}`. "
        f"components/predictor_context.py prints {ours} as that cutoff on the variant panel; "
        "find where the script's threshold went and move the app's copy with it."
    )
    assert float(assignments[variable]) == float(ours), (
        f"{os.path.basename(path)}:{function} sets `{variable}` to {assignments[variable]}, "
        f"and components/predictor_context.py prints {ours}. The panel states this number to a "
        "clinician as the tool's own cutoff, so one of the two is now lying."
    )


def test_the_functions_the_section_speaks_for_still_exist():
    """The two criteria have not been renamed out from under the section.

    A ``check_PreP`` or ``check_PP3`` that vanished would make every parametrised case above skip
    its assertion by finding an empty dict — so the existence of the functions is asserted once,
    separately, rather than inferred from the cutoffs being found.
    """
    with open(CANCERVAR, encoding="utf-8") as handle:
        cancervar = ast.parse(handle.read())
    with open(INTERVAR, encoding="utf-8") as handle:
        intervar = ast.parse(handle.read())

    def names(tree):
        return {n.name for n in ast.walk(tree) if isinstance(n, ast.FunctionDef)}

    assert "check_PreP" in names(cancervar), "CancerVar's CBP10 function is gone"
    for function in ("check_PP3", "check_BP4"):
        assert function in names(intervar), f"InterVar's {function} is gone"


def test_the_guard_would_notice_a_moved_cutoff(tmp_path):
    """The parser reads a real assignment, and would report a changed one.

    Asserted rather than assumed, because a guard whose parser silently returns ``{}`` passes every
    comparison it is asked to make. Written against a stand-in script so it does not depend on the
    vendored one staying as it is.
    """
    script = tmp_path / "fake.py"
    script.write_text(
        "sift_cutoff = 999  # module level, must not answer for the function\n"
        "def check_PreP(line):\n"
        "    # sift_cutoff = 0.05 in a comment must not answer either\n"
        "    sift_cutoff = 0.07\n"
        "    metasvm_cutoff = -1.5\n"
        "    return sift_cutoff\n"
        "def other(line):\n"
        "    sift_cutoff = 0.05\n"
        "    return sift_cutoff\n",
        encoding="utf-8",
    )
    found = _assignments(str(script), "check_PreP")
    assert found["sift_cutoff"] == 0.07, "the parser read the wrong scope"
    assert found["metasvm_cutoff"] == -1.5, "a negative cutoff was not read"
    assert _assignments(str(script), "missing") == {}


def test_the_two_tools_read_gerp_differently_and_the_app_keeps_both():
    """One column, two cutoffs — and the app must not collapse them into one constant.

    CancerVar's CBP10 calls ``GERP++_RS`` conserved at ``>= 2`` and InterVar's PP3 at ``> 2``. The
    numbers are equal and the comparisons are not, so the two constants are separately named and
    separately attributed on screen. This asserts the comparisons are still what the section
    claims, since equal numbers make the difference invisible to the cases above.
    """
    with open(CANCERVAR, encoding="utf-8") as handle:
        cancervar = handle.read()
    with open(INTERVAR, encoding="utf-8") as handle:
        intervar = handle.read()

    assert 'if float(cls[Funcanno_flgs["GERP++_RS"]]) >= cutoff_conserv:' in cancervar, (
        "CancerVar's GERP++ comparison is no longer `>=`. The panel says "
        "'at or above CancerVar's 2 cutoff' for the somatic arm."
    )
    assert 'if float(cls[Funcanno_flgs["GERP++_RS"]]) > cutoff_conserv:' in intervar, (
        "InterVar's PP3 GERP++ comparison is no longer `>`. The panel says "
        "'above InterVar's 2 cutoff' for the germline arm."
    )
