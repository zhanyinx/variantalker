"""The arrangement of ``components/``, held against the account it gives of itself.

A guard rather than a unit module: nothing here calls the app. Every check reads the
package's source with ``ast`` and asserts a *shape* — which module may import which, what
the package re-exports, whether its own module table still lists the modules on disk. So
these hold in a checkout with no Streamlit and no ``st_aggrid``, which is the point: the
arrangement is what issue #76 decided, and it should be checkable without booting anything.

Issue #76 split ``components/ui_components.py`` (1,536 lines) and ``visualizations.py``
(1,142) into seven modules. Two of the three things that split protected are invisible to
any behavioural test — a cycle and a re-export both leave a working app — and the third,
the hand-written module table, is documentation that can quietly stop being true.
"""

from __future__ import annotations

import ast
from pathlib import Path

COMPONENTS_DIR = Path(__file__).resolve().parent.parent / "components"


def _component_modules() -> set[str]:
    return {path.stem for path in COMPONENTS_DIR.glob("*.py") if path.stem != "__init__"}


def _component_imports() -> dict[str, set[str]]:
    """``module -> the sibling modules it imports``, read statically."""
    modules = _component_modules()
    graph: dict[str, set[str]] = {}
    for path in COMPONENTS_DIR.glob("*.py"):
        if path.stem == "__init__":
            continue
        siblings = set()
        for node in ast.walk(ast.parse(path.read_text())):
            if isinstance(node, ast.ImportFrom):
                # `from .charts import x` (level 1) or `from components.charts import x`
                name = node.module or ""
                if node.level:
                    head = name.split(".")[0]
                elif name.startswith("components."):
                    head = name.split(".")[1]
                else:
                    continue
                if head in modules:
                    siblings.add(head)
        graph[path.stem] = siblings
    return graph


def test_the_components_depend_on_each_other_in_one_direction_only():
    """No cycles, which is what let the deferred imports go.

    Two imports in the old ``ui_components.py`` were written inside the functions that
    needed them, to dodge a cycle through this package's own ``__init__``. Both are
    module-level now because the split left an acyclic graph, and a cycle reintroduced
    here would send the next reader back to hiding imports in function bodies rather than
    to the module that should not have been imported.
    """
    graph = _component_imports()

    visiting, done = set(), set()

    def walk(module, trail):
        if module in done:
            return
        assert module not in visiting, (
            f"import cycle in components/: {' -> '.join(trail + [module])}"
        )
        visiting.add(module)
        for sibling in sorted(graph.get(module, ())):
            walk(sibling, trail + [module])
        visiting.discard(module)
        done.add(module)

    for module in sorted(graph):
        walk(module, [])


def test_the_package_re_exports_nothing():
    """The mechanism that hid seven dead functions, kept shut.

    ``components/__init__.py`` opened with ``from .ui_components import *`` and listed
    twelve names in ``__all__``. Seven named functions with no caller anywhere, leaving
    five that were really read — and the ``__all__`` entry was the reason nothing reported
    them, since a name re-exported from a package reads as API kept for its callers rather
    than as code nobody calls. ``config/__init__.py`` lost the same wildcard in issue #54.
    """
    tree = ast.parse((COMPONENTS_DIR / "__init__.py").read_text())

    imports = [n for n in ast.walk(tree) if isinstance(n, (ast.Import, ast.ImportFrom))]
    assert not imports, (
        "components/__init__.py imports again. Import from the module that defines a "
        "thing: a re-export here makes a dead name look like a kept one (issue #76)."
    )

    assigned = [
        target.id
        for node in tree.body
        if isinstance(node, ast.Assign)
        for target in node.targets
        if isinstance(target, ast.Name)
    ]
    assert "__all__" not in assigned, (
        "components/__init__.py declares __all__ again — see issue #76"
    )


def test_the_packages_own_module_table_names_the_modules_that_exist():
    """``__init__.py`` is now documentation, so it can rot without anything failing.

    Having deleted the re-exports, the only thing left in that file is a table saying what
    each module holds and a diagram of which imports which. That is the map a reader
    reaches for first, and nothing executes it. A module added or renamed without touching
    the table leaves the file describing a package that no longer exists — the same defect
    as a comment that has stopped being true, which this repo treats as a defect.
    """
    source = (COMPONENTS_DIR / "__init__.py").read_text()
    listed = {
        line.strip().strip("|").split("|")[0].strip().strip("`")
        for line in source.splitlines()
        if line.startswith("| `")
    }
    on_disk = _component_modules()

    assert listed == on_disk, (
        "components/__init__.py's module table and the modules on disk disagree.\n"
        f"  in the table, not on disk: {sorted(listed - on_disk)}\n"
        f"  on disk, not in the table: {sorted(on_disk - listed)}"
    )
