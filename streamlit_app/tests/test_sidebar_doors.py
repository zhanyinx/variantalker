"""The sidebar offers no second door onto a page its nav radio already lists.

A guard rather than a unit module: nothing here calls the app. It reads the source of
``components/sidebar.py`` and asserts a *shape*, so it holds in a checkout with no Streamlit
installed at all.

Issue #161 deleted a ``❓ Need Help? Click here`` button that set ``current_page = "help"`` —
the same page the nav radio's fourth entry opens, two elements above it and behind the same
❓ glyph. What is left in the column is one button, ``📊 Back to your results``, and the
difference between the two is the rule this file holds: that button asks the data page for a
*section* as well as for the page (``to_results``), which the radio cannot express. Naming a
page the radio already lists, and nothing else, is not a second invitation — it is the same
one, and on the page itself it is worse than the same one, because it offers a trip to where
you already are.

The rule is worth a guard rather than a comment because the button was **not** an accident.
Issue #58 absorbed the ``📊 Session Info`` expander into the load-status block and
deliberately kept it, on the argument that a radio entry reads as a place while a button
reads as an offer of help. That argument will be made again, and the next author making it
should have to answer *what can this reach that the radio cannot?* first.

Three syntaxes, not one, and that is deliberate: a rule over ``_nav_button`` calls alone is
walked past by writing the button out longhand, or by hanging the navigation off an
``on_click`` callback. All three are checked below, and each one has a caller in this module
today — the sanctioned route, the clear-list button, and the file chooser's callback — so
none of the three sweeps is looking at an empty set.
"""

from __future__ import annotations

import ast
from pathlib import Path

SIDEBAR = Path(__file__).resolve().parent.parent / "components" / "sidebar.py"

#: The one function allowed to turn a sidebar button into a navigation.
#:
#: Named rather than derived, because the shape it is excluded for — ``if ui.button(...):``
#: followed by a write to ``current_page`` — is exactly the shape the rule is about. It is
#: the sanctioned route precisely because it is the one place that takes ``to_results``, so
#: every button routed through it can be asked what it reaches.
SANCTIONED_ROUTE = "_nav_button"

#: The session key that moves the user between pages.
CURRENT_PAGE = "current_page"


def _tree() -> ast.Module:
    return ast.parse(SIDEBAR.read_text(encoding="utf-8"))


def _page_labels(tree: ast.Module) -> list[str]:
    """The pages the nav radio lists, read off ``PAGE_LABELS`` in the source.

    Read rather than imported: importing this module needs Streamlit, and a guard that
    stops running when a dependency is absent is issue #24's failure returning.
    """
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign) and any(
            isinstance(t, ast.Name) and t.id == "PAGE_LABELS" for t in node.targets
        ):
            assert isinstance(node.value, ast.Dict), "PAGE_LABELS is no longer a dict literal."
            return [k.value for k in node.value.keys if isinstance(k, ast.Constant)]
    raise AssertionError(
        "components/sidebar.py no longer defines PAGE_LABELS. Every rule in this file is "
        "about the pages the nav radio lists, so it cannot be checked without that table."
    )


def _functions(tree: ast.Module) -> dict[str, ast.FunctionDef]:
    return {
        node.name: node
        for node in ast.walk(tree)
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
    }


def _writes_current_page(node: ast.AST) -> bool:
    """Whether anything under ``node`` assigns the page the user is on.

    Both spellings, because they are equivalent and the app uses each of them somewhere:
    ``st.session_state.current_page = ...`` and ``st.session_state["current_page"] = ...``.
    """
    for child in ast.walk(node):
        if not isinstance(child, (ast.Assign, ast.AugAssign, ast.AnnAssign)):
            continue
        targets = child.targets if isinstance(child, ast.Assign) else [child.target]
        for target in targets:
            if isinstance(target, ast.Attribute) and target.attr == CURRENT_PAGE:
                return True
            if (
                isinstance(target, ast.Subscript)
                and isinstance(target.slice, ast.Constant)
                and target.slice.value == CURRENT_PAGE
            ):
                return True
    return False


def _button_calls(node: ast.AST) -> list[ast.Call]:
    """Every ``<anything>.button(...)`` under ``node``.

    The receiver is not pinned to ``st``: this module draws buttons on ``st.sidebar``, on a
    passed-in ``ui``, and on an expander's ``box``, and a rule that only recognised one of
    those would be blind to the other two.
    """
    return [
        child
        for child in ast.walk(node)
        if isinstance(child, ast.Call)
        and isinstance(child.func, ast.Attribute)
        and child.func.attr == "button"
    ]


def _asks_for_more_than_the_page(call: ast.Call) -> bool:
    """Whether a ``_nav_button`` call asks for something the nav radio cannot express.

    ``to_results`` today, positionally or by keyword. An explicit ``False`` does not count;
    anything else does, including an expression — the live call site passes
    ``to_results=filtered_data is not None``, because a file with no results behind it has no
    Results section to land on, and pinning this to a literal ``True`` would fail the one
    button the rule exists to permit.
    """
    for keyword in call.keywords:
        if keyword.arg == "to_results":
            return not (isinstance(keyword.value, ast.Constant) and keyword.value.value is False)
    if len(call.args) >= 4:
        return not (isinstance(call.args[3], ast.Constant) and call.args[3].value is False)
    return False


def test_no_sidebar_button_merely_repeats_a_nav_radio_entry():
    """A button naming a page the radio lists must reach more than the page."""
    tree = _tree()
    pages = _page_labels(tree)
    calls = [
        node
        for node in ast.walk(tree)
        if isinstance(node, ast.Call)
        and isinstance(node.func, ast.Name)
        and node.func.id == SANCTIONED_ROUTE
    ]

    # Vacuity, asserted rather than assumed: every claim below is quantified over these
    # calls, so a rename that parsed to none of them would leave this test green while
    # saying nothing. The app has one such call and needs at least one — it is the way back
    # to your results.
    assert calls, (
        f"No {SANCTIONED_ROUTE} call sites found in {SIDEBAR.name}. Either the sidebar's "
        "navigation has been rewritten, in which case this guard needs rewriting with it, "
        "or the parse is broken — and a guard that catches nothing reads exactly like a "
        "guard with nothing to catch."
    )

    for call in calls:
        page = call.args[2] if len(call.args) >= 3 else None
        assert isinstance(page, ast.Constant), (
            f"{SIDEBAR.name}:{call.lineno} calls {SANCTIONED_ROUTE} with a destination this "
            "guard cannot read. Pass the page as a literal so the rule below can be checked "
            "— a computed destination is how a second door onto Help gets past this file."
        )
        if page.value in pages and not _asks_for_more_than_the_page(call):
            raise AssertionError(
                f"{SIDEBAR.name}:{call.lineno} puts a button in the sidebar that goes to "
                f"{page.value!r} — a page the nav radio already lists, in the same column, "
                "with nothing asked for beyond the page itself. That is the button issue "
                "#161 deleted: on every page it is the radio entry above it, and on "
                f"{page.value!r} itself it offers a trip to where the user already is. A "
                "button here earns its place by reaching what the radio cannot — a section, "
                "as `to_results` asks for."
            )


def test_every_navigating_sidebar_button_goes_through_the_sanctioned_route():
    """So the rule above cannot be walked past by writing the button out longhand."""
    tree = _tree()
    offenders = []
    for name, function in _functions(tree).items():
        if name == SANCTIONED_ROUTE:
            continue
        for node in ast.walk(function):
            if not isinstance(node, ast.If) or not _button_calls(node.test):
                continue
            if _writes_current_page(ast.Module(body=node.body, type_ignores=[])):
                offenders.append(f"{name} (line {node.lineno})")

    assert not offenders, (
        f"These branch on a button and then set {CURRENT_PAGE} themselves: "
        f"{', '.join(offenders)}. Route them through {SANCTIONED_ROUTE}, which is where the "
        "rerun the nav radio needs lives, and where a destination can be checked against the "
        "pages the radio already lists."
    )


def test_no_sidebar_button_navigates_from_a_callback_instead():
    """The third syntax: navigation hung off ``on_click`` rather than a return value.

    The file chooser really does navigate from a callback — choosing a MAF carries you to the
    page that loads it (issue #64) — and that is not what this rule is about: it is an
    uploader's ``on_change``, not a button, and the navigation is a consequence of opening a
    file rather than the whole of what the control does. A *button* whose callback moves the
    user is the deleted Help button with one more level of indirection.
    """
    tree = _tree()
    functions = _functions(tree)
    offenders = []
    for call in _button_calls(tree):
        for keyword in call.keywords:
            if keyword.arg != "on_click" or not isinstance(keyword.value, ast.Name):
                continue
            callback = functions.get(keyword.value.id)
            if callback is not None and _writes_current_page(callback):
                offenders.append(f"line {call.lineno} -> {keyword.value.id}")

    assert not offenders, (
        f"These sidebar buttons navigate from their callback: {', '.join(offenders)}. The "
        f"rule in this file is about where a sidebar button can take you, so a button that "
        f"sets {CURRENT_PAGE} in an `on_click` is inside it, not outside it."
    )
