"""Nothing the app says is drawn into a frame ``st.rerun`` throws away.

A guard rather than a unit module: nothing here calls the app. It reads the source of every
module the app renders from and asserts a *shape* — that no message-drawing call stands
before an ``st.rerun`` in the same block — so it holds in a checkout with no Streamlit
installed at all.

``st.rerun`` raises immediately and the run's output is discarded, so a sentence written on
the way to one reaches nobody. That is not a hazard this app avoids; it is a defect found
across it. Issue #133 found five instances on the parameter page and issue #140 the last
three, in two places: the variant dialog's *Saved for this session.* — the app's only
statement of how long a note lasts — and the page router's two failure branches, which
bounced a user to Home with no explanation at all.

Every one of them read as working, and that is the reason this file exists rather than a
larger test. ``AppTest`` reports an element drawn before a rerun **when it sits inside a
container**, measured with a one-element probe in issue #133, so a behavioural test can
assert the message and pass while no user has ever seen it. The instrument agrees with the
bug. Reading the source does not.

There is no table of sanctioned exceptions here, unlike the sweep in
``test_parameter_adoption.py`` next door, and the difference is not an oversight. That rule
governs writes that are *legitimate in the right place*, so the right shape is a list of
places with reasons. This one governs an act with no legitimate case: a message the user
cannot read is never what the author meant. The two ways to keep one are to park it for the
run after the rerun, or not to rerun.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

STREAMLIT_APP = Path(__file__).resolve().parent.parent

#: The calls that put something on the screen.
#:
#: ``st.toast`` is in the list, and it is the member most likely to be argued about. It is
#: not exempt: in Streamlit 1.47 a toast is an ordinary element in the delta tree — it is
#: *not* routed to the event container the way ``st.dialog`` and ``st.html`` are — and the
#: frontend clears it when its node is pruned, which happens as soon as the following run
#: finishes. So a toast drawn before a rerun is visible for exactly as long as the rerun
#: takes. It is a fine thing to draw *after* one, which is what issue #140 does with it.
#:
#: Wider than the five alert boxes issue #140 set out to find, and review is why: this
#: change's own new sentence — *You can pick another page in the sidebar* — is an
#: ``st.write``, so a guard listing only the alert boxes could be walked past by choosing a
#: different call, including by the commit that added it. The rule is not about alert boxes;
#: a discarded frame discards everything drawn into it. Nothing in the app is caught by the
#: three extra names today, so they cost nothing and close the hole.
DRAWERS = frozenset(
    {
        "success",
        "info",
        "warning",
        "error",
        "toast",
        "exception",
        "write",
        "markdown",
        "caption",
    }
)


def _app_modules() -> list[Path]:
    """Every module of the app a discarded message could hide in.

    App-wide rather than one file, for the reason the sweep in
    ``test_parameter_adoption.py`` gives: the rule is about the app. ``MAFigate.py`` is in
    the list because that is where two of the three known instances were, and a guard over
    ``page_modules/`` and ``components/`` alone would have found neither.
    """
    modules = [STREAMLIT_APP / "MAFigate.py"]
    for package in ("page_modules", "components"):
        modules.extend(sorted((STREAMLIT_APP / package).glob("*.py")))
    return modules


def _st_call(node: ast.stmt) -> str | None:
    """The ``st.<name>`` a bare expression statement calls, if it is one."""
    if not isinstance(node, ast.Expr) or not isinstance(node.value, ast.Call):
        return None
    func = node.value.func
    if (
        isinstance(func, ast.Attribute)
        and isinstance(func.value, ast.Name)
        and func.value.id == "st"
    ):
        return func.attr
    return None


def _on_one_path(block: list) -> list:
    """The statements in ``block`` that all run, in order, flattening ``with`` bodies.

    A ``with`` body runs unconditionally when the ``with`` is reached, so a message inside
    one and a rerun after it are on the same path. Flattening matters more here than
    anywhere else in this file, because a container is **exactly** the shape that hides this
    defect from the other instrument: ``AppTest`` reports an element drawn before a rerun
    when it sits inside a container and not at page level (issue #133), so

    .. code-block:: python

        with st.container():
            st.success("Saved.")
        st.rerun()

    is the one arrangement that reads as working from both a test and the source. Found by
    review; no site in the app has this shape today.

    ``if``/``else`` branches are deliberately *not* flattened — a message in one branch and a
    rerun in another are not on the same path, and reporting them would be the false positive
    that teaches a reader to ignore this file. ``try`` bodies are not flattened either, for a
    weaker reason: a statement there may or may not be reached, so a claim about it would be
    a guess. Both are still swept in their own right by :func:`_statement_blocks`.
    """
    flattened = []
    for statement in block:
        if isinstance(statement, (ast.With, ast.AsyncWith)):
            flattened.extend(_on_one_path(statement.body))
        else:
            flattened.append(statement)
    return flattened


def _statement_blocks(tree: ast.AST):
    """Every list of sibling statements in the module.

    Siblings rather than *anywhere lexically before*, and the difference is the whole
    accuracy of this guard. A looser sweep over one function at a time reports thirteen
    sites in this app and ten are false: a message in one branch of an ``if`` and a rerun in
    another are not on the same path, and calling them a defect would teach a reader to
    ignore this file. Siblings in one block share a path by construction.
    """
    for node in ast.walk(tree):
        for field in ("body", "orelse", "finalbody"):
            block = getattr(node, field, None)
            if isinstance(block, list) and block and isinstance(block[0], ast.stmt):
                yield block


def _messages_drawn_into_a_discarded_frame(tree: ast.AST) -> list[tuple[int, str]]:
    """``(line, drawer)`` for everything this tree draws before rerunning.

    Takes a parsed tree rather than a path so that the tests below can drive *this* function
    over fixture source instead of re-implementing it. That is not tidiness: an earlier
    version of this file had them re-implement the loop, which left the slicing here — the
    one line that decides what counts as "before" — exercised by nothing, so a bug in it
    would have made all eighteen parametrised cases and their own anti-vacuity test green
    together. Caught by review.
    """
    found = set()
    for block in _statement_blocks(tree):
        calls = [(statement, _st_call(statement)) for statement in _on_one_path(block)]
        reruns = [i for i, (_, name) in enumerate(calls) if name == "rerun"]
        if not reruns:
            continue
        for statement, name in calls[: reruns[0]]:
            if name in DRAWERS:
                # A set, because flattening makes one `with` body reachable from two blocks
                # — its own, and its parent's — and a sentence reported twice reads as two
                # defects.
                found.add((statement.lineno, name))
    return sorted(found)


@pytest.mark.parametrize("module", _app_modules(), ids=lambda path: path.name)
def test_no_message_is_drawn_into_a_frame_that_is_discarded(module):
    """A message before a rerun, in the same block, is a message nobody reads."""
    drawn = _messages_drawn_into_a_discarded_frame(ast.parse(module.read_text()))
    assert drawn == [], (
        f"{module.name} draws a message immediately before `st.rerun`, which discards the "
        "frame it was drawn into, so nobody reads it: "
        + ", ".join(f"st.{name} at line {line}" for line, name in drawn)
        + ". Park the sentence for the run after the rerun, or do not rerun."
    )


def test_the_sweep_can_see_a_message_that_is_discarded():
    """The guard above, made to fail on purpose.

    Without this, a bug anywhere in the sweep would make every parametrised case above
    vacuously green, and a green run would mean the parser found nothing rather than that
    the app does nothing wrong — the failure mode issue #133's ``AppTest`` probe and issue
    #77's ``deepcopy(X) == X`` were both instances of.

    It drives :func:`_messages_drawn_into_a_discarded_frame` itself, the same function the
    parametrised cases call, rather than re-running its loop over a fixture. The first
    version of this file did the latter and review caught it: the two would have gone green
    together over a broken wrapper.

    The source is the shape the app used to have at ``components/variant_table.py``, kept as
    a string so that reinstating the defect is not required to know this works.
    """
    module = ast.parse(
        "import streamlit as st\n"
        "\n"
        "def save():\n"
        "    if st.button('Save'):\n"
        "        st.success('Saved for this session.')\n"
        "        st.rerun()\n"
    )
    assert _messages_drawn_into_a_discarded_frame(module) == [(5, "success")]


def test_the_sweep_sees_through_a_container():
    """The shape that fools the *other* instrument does not fool this one.

    ``AppTest`` reports an element drawn before a rerun when it sits inside a container, and
    not at page level (issue #133, measured with a one-element probe). So a message wrapped
    in ``st.container`` and then discarded is the arrangement that reads as working
    everywhere — which is exactly why the sweep flattens ``with`` bodies into the block
    around them. Review found this gap; nothing in the app has the shape.
    """
    module = ast.parse(
        "import streamlit as st\n"
        "\n"
        "def save():\n"
        "    if st.button('Save'):\n"
        "        with st.container():\n"
        "            st.success('Saved for this session.')\n"
        "        st.rerun()\n"
    )
    assert _messages_drawn_into_a_discarded_frame(module) == [(6, "success")]


def test_a_message_in_the_other_branch_is_not_reported():
    """The false positive this guard refuses to make.

    A message in one branch of an ``if`` and a rerun in another are not on the same path, so
    the first is drawn on runs where the second never happens. Reporting it would be the
    kind of noise that teaches a reader to stop reading the failure — and the app has ten
    such pairs, every one of them fine.
    """
    module = ast.parse(
        "import streamlit as st\n"
        "\n"
        "def save():\n"
        "    if st.button('Save'):\n"
        "        st.rerun()\n"
        "    else:\n"
        "        st.info('Nothing to save yet.')\n"
    )
    assert _messages_drawn_into_a_discarded_frame(module) == []


def test_a_parked_message_is_not_reported():
    """The sanctioned shape passes, so the guard is a rule and not a ban on reruns.

    Parking is what issues #133 and #140 both landed on: the sentence goes into session
    state and the run *after* the rerun draws it. A guard that could not tell the two apart
    would push authors back towards saying nothing.
    """
    module = ast.parse(
        "import streamlit as st\n"
        "\n"
        "def save():\n"
        "    if st.button('Save'):\n"
        "        park_note_confirmation('Saved for this session.')\n"
        "        st.rerun()\n"
        "\n"
        "def page():\n"
        "    st.toast('Saved for this session.')\n"
    )
    assert _messages_drawn_into_a_discarded_frame(module) == []
