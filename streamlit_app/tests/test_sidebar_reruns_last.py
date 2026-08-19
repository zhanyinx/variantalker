"""The sidebar takes its reruns after the file chooser has been drawn, or not at all.

A guard rather than a unit module: nothing here calls the app. It reads the source of
``components/sidebar.py`` and asserts a *shape*, so it holds in a checkout with no Streamlit
installed at all — the same reason ``tests/test_sidebar_doors.py`` is written this way.

The rule comes from issue #283, and it is about a Streamlit mechanic rather than about this
app's taste. ``st.rerun`` is **not** a premature stop — ``scriptrunner/exec_code.py`` sets
``premature_stop = False`` on ``RerunException`` — so Streamlit still runs
``on_script_finished(ctx.widget_ids_this_run)`` when a run ends that way, and that prunes the
state of every widget the aborted run did not register. The sidebar's file chooser is the last
thing the app draws: it lives in the status slot that ``render_into_status_slot`` fills after
the page has run. So *any* rerun sited earlier in this module drops the chooser's state, the
rerun that follows re-registers it at ``None`` while the browser is still holding the file, and
the user's next interaction anywhere in the app reads as a fresh choice and fires the chooser's
``on_change`` — which navigates to the data page, before the page whose control was touched has
run. That was issue #277: every selection on the parameters page was thrown away.

Worth a guard rather than a comment because the reruns this rule forbids were each put there
deliberately and read as obviously correct. Two of them turned out to have no work to do at all
(``create_sidebar_navigation`` sets ``current_page`` before ``MAFigate.main`` routes, and
deleting a widget key takes effect in the same run), and the third — the route back to your
results, drawn *above* the chooser — genuinely needs one, which is why the rule is *where*
rather than *whether*. The next author who needs a rerun in this column should have to answer
*is the chooser drawn yet?* first, and the honest answer is available only at the bottom of
``render_into_status_slot``.

A rule over ``st.rerun`` alone is enough here, and deliberately so: the two other ways a run
can end early are not this hazard. ``st.stop`` raises ``StopException``, which *is* a premature
stop and prunes nothing, and an uncaught exception is a broken render rather than a navigation.

**What the rule protects is every widget the run draws, not only the chooser.** The chooser is
simply the last of them, so a rerun below it has necessarily let the page body run too — which
is why ``_keep_a_section_selected`` on the data page, the app's only other ``on_change`` and
exposed to exactly the same prune-and-re-fire shape, is covered by this file without being
named in it. That is the ``also check`` issue #283 asked for, discharged by the shape of the
rule rather than by a second test.

Two limits, named because a guard that looks wider than it is reads as a promise: this file
parses ``components/sidebar.py`` and nothing else, so a rerun moved into a helper module
escapes it (contrast ``tests/test_discarded_frames.py``, which sweeps ``page_modules/`` and
``components/``) — and that sibling guard's own rule, *no message drawn before an ``st.rerun``
in the same block*, cannot see the one rerun this file permits, because ``render_load_status``
draws through a call rather than as a bare ``st.<name>`` expression. Nothing is lost by that
today: no sentence is parked and drained down there at all since the file history went, and
the one that used to be — its *cleared* confirmation — was set by an ``on_click``, and a run
that fires a callback is not a run in which ``_nav_button`` returns True. A control parked
below the status block would need checking by hand.
"""

from __future__ import annotations

import ast
from pathlib import Path

SIDEBAR = Path(__file__).resolve().parent.parent / "components" / "sidebar.py"

#: The one function allowed to call ``st.rerun``, and the call that has to come first.
#:
#: Named as a pair, because the permission is not "this function" — it is "after the column is
#: complete". ``render_into_status_slot`` is where the chooser reaches the screen, by way of
#: ``render_load_status``, so a rerun *above* that call inside this very function would be the
#: same defect with a shorter reach.
DRAINING_FUNCTION = "render_into_status_slot"
DRAWS_THE_CHOOSER = "render_load_status"


def _tree() -> ast.Module:
    return ast.parse(SIDEBAR.read_text(encoding="utf-8"))


def _calls_named(node: ast.AST, name: str) -> list[ast.Call]:
    """Every call to ``name`` under ``node``, however it is reached.

    Matched on the final name rather than on the full dotted path, so ``st.rerun()`` and a
    bare ``rerun()`` from ``from streamlit import rerun`` both count, and so does a rerun
    called on a module passed in under another name. Both rules in this file quantify over
    call sites — one over the reruns, one over the call that draws the column — and a matcher
    that recognised only one spelling would let the other be renamed past it.
    """
    return [
        child
        for child in ast.walk(node)
        if isinstance(child, ast.Call)
        and (
            (isinstance(child.func, ast.Attribute) and child.func.attr == name)
            or (isinstance(child.func, ast.Name) and child.func.id == name)
        )
    ]


def _rerun_calls(node: ast.AST) -> list[ast.Call]:
    """Every ``st.rerun(...)`` under ``node``."""
    return _calls_named(node, "rerun")


def _named_functions(tree: ast.Module) -> dict[str, ast.FunctionDef]:
    return {
        node.name: node
        for node in ast.walk(tree)
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
    }


def test_the_sidebars_only_rerun_is_the_one_the_status_block_drains():
    """No rerun anywhere in the module except after the chooser has been drawn."""
    tree = _tree()
    calls = _rerun_calls(tree)

    # Vacuity, asserted rather than assumed: every claim below is quantified over these calls,
    # so a module with none of them would pass this file while saying nothing. The sidebar has
    # exactly one and needs it — the route back to your results is drawn after the page has
    # already rendered, so nothing else would put the new page on screen.
    assert calls, (
        f"No st.rerun call sites found in {SIDEBAR.name}. Either the sidebar's navigation has "
        "been rewritten — in which case this guard needs rewriting with it — or the parse is "
        "broken, and a guard that catches nothing reads exactly like a guard with nothing to "
        "catch."
    )

    functions = _named_functions(tree)
    drainer = functions.get(DRAINING_FUNCTION)
    assert drainer is not None, (
        f"{SIDEBAR.name} no longer defines {DRAINING_FUNCTION}. That function is where the "
        "file chooser reaches the screen, so it is the only place in this module where a rerun "
        "is safe; without it this rule cannot be checked."
    )

    draws = _rerun_calls(drainer)
    allowed = {id(call) for call in draws}
    strays = [call for call in calls if id(call) not in allowed]
    assert not strays, (
        "These call st.rerun before the sidebar's file chooser has been drawn: "
        + ", ".join(f"{SIDEBAR.name}:{call.lineno}" for call in strays)
        + ". Streamlit prunes the state of every widget a rerunning run did not register, so "
        "the chooser comes back holding None while the browser still holds the file, and the "
        "user's next interaction anywhere in the app fires its on_change and navigates away "
        "(issues #277, #283). Park the request in a session key instead and let "
        f"{DRAINING_FUNCTION} take it, which is what NAV_RERUN_PENDING is for."
    )


def test_that_rerun_comes_after_the_chooser_and_not_before_it():
    """Inside the draining function, the rerun follows the call that draws the column.

    Quantified over the reruns in that one function, so it says nothing about a module with
    none — and that is deliberate rather than an oversight: the test above refuses a module
    with no rerun at all *and* refuses one sited anywhere else, so between them the set this
    loop walks is never empty. Asserting non-emptiness here as well would fail a future
    sidebar that legitimately needs no deferred rerun because it draws the chooser above the
    button, which is a fix rather than a regression.
    """
    tree = _tree()
    drainer = _named_functions(tree)[DRAINING_FUNCTION]

    renders = _calls_named(drainer, DRAWS_THE_CHOOSER)
    assert renders, (
        f"{DRAINING_FUNCTION} no longer calls {DRAWS_THE_CHOOSER}, which is what draws the "
        "status block and the file chooser inside it. The rule in this file is about a rerun "
        "sited after that call, so it cannot be checked without it."
    )
    drawn_by = max(call.lineno for call in renders)

    for call in _rerun_calls(drainer):
        assert call.lineno > drawn_by, (
            f"{SIDEBAR.name}:{call.lineno} reruns before {DRAWS_THE_CHOOSER} has drawn the "
            f"file chooser on line {drawn_by}. Being inside {DRAINING_FUNCTION} is not the "
            "permission — having drawn the column is (issues #277, #283)."
        )
