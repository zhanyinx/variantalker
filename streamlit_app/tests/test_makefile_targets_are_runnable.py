"""Every script the Makefile runs is in the tree that ships with the Makefile.

Issue #345, and the reason it is a guard rather than a one-line fix. `make app-load-check`
ran `../docs/wayfinder/issue-32/run_app_check.py` for six months. It worked on every
developer's machine and failed in **every clone of the public repository**, because
`.publicignore` strips `/docs/wayfinder/` entirely — so the one target
`streamlit_app/README.md` calls the quickest evidence a change has not broken the app
stopped with a missing-file traceback for exactly the people who had just cloned the repo
and wanted to check their setup worked. Nobody noticed until somebody cloned it and looked.

The two tests here are deliberately the same claim read from opposite ends, because neither
tree can see the bug on its own:

* :func:`test_every_script_a_recipe_runs_is_present` fails in an **exported** tree, where
  the file is genuinely gone. That is where the bug was visible, and where nobody runs the
  suite before publishing.
* :func:`test_no_recipe_reaches_into_a_directory_the_export_strips` fails in the **private**
  tree, the day someone points a recipe at a deny-listed path. That is where the bug can
  still be cheaply fixed, and it is the check that would have refused #345 at birth.

A recipe reaching outside `streamlit_app/` is not itself the defect — `..` is ordinary in a
monorepo, and `make test-cov` legitimately builds a throwaway tree above this one. Reaching
into a directory the export strips is the defect, so that is what the second test names.

WHY NOT ASK THE EXPORT GATE. The export script reads the same deny-list and would be
the tidier authority, but it is itself deny-listed: importing it here would make this module
fail in the tree it most needs to pass in. So the rules are re-read from `.publicignore`
directly, and where that file is absent — which is exactly the exported tree — the second
test skips rather than passing vacuously. That shape is `test_suite_organisation.py`'s, from
issue #348: a check that cannot run says so instead of going quiet.
"""

import re
import unittest
from pathlib import Path

APP_DIR = Path(__file__).resolve().parents[1]
REPO_ROOT = APP_DIR.parent
MAKEFILE = APP_DIR / "Makefile"
DENY_LIST = REPO_ROOT / ".publicignore"

#: Two scripts this reading must find, or it is not reading the recipes. `run_app_check.py`
#: is the one this module exists for; `check_dependencies.py` is the oldest and least likely
#: to move, so a failure naming only the first is a moved file rather than a broken parse.
ANCHOR_SCRIPTS = {"tests/run_app_check.py", "check_dependencies.py"}

#: Below the eight distinct scripts the recipes name today, above anything a broken regex
#: would return.
MINIMUM_SCRIPTS = 6


def _scripts_the_recipes_run() -> set[str]:
    """Every ``.py`` path a recipe hands to ``$(PYTHON)``, as written in the Makefile.

    Only literal paths. ``-m pytest``, ``-m pip``, ``-m streamlit`` and ``-c "…"`` name no
    file, and anything still holding a ``$`` after the match is a make or shell expansion
    this module cannot resolve and must not guess at.
    """
    text = MAKEFILE.read_text(encoding="utf-8")
    found = set()
    for match in re.finditer(r"\$\(PYTHON\)\s+(-m\s+\S+|\S+)", text):
        argument = match.group(1)
        if argument.startswith("-") or not argument.endswith(".py"):
            continue
        if "$" in argument:
            continue
        found.add(argument)
    return found


def _stripped_directories() -> list[str]:
    """The directory rules in ``.publicignore``, as repository-relative prefixes.

    A rule ending in ``/`` strips that directory and everything under it; anything else
    names one exact file, which a recipe path can also hit, so both shapes are returned —
    directories with their trailing slash intact so a prefix test cannot match a sibling
    whose name merely starts the same way.
    """
    rules = []
    for line in DENY_LIST.read_text(encoding="utf-8").splitlines():
        rule = line.strip()
        if not rule or rule.startswith("#"):
            continue
        rules.append(rule.lstrip("/"))
    return rules


class TestTheRecipesCanRun(unittest.TestCase):
    def setUp(self):
        self.scripts = _scripts_the_recipes_run()

    def test_the_reading_found_the_recipes_at_all(self):
        """Anchored first, so the two assertions below cannot pass over an empty set."""
        self.assertGreaterEqual(
            len(self.scripts),
            MINIMUM_SCRIPTS,
            f"read only {len(self.scripts)} script path(s) out of the Makefile "
            f"({sorted(self.scripts)}); its shape changed and this module is no longer "
            "reading the recipes",
        )
        missing_anchors = ANCHOR_SCRIPTS - self.scripts
        self.assertFalse(
            missing_anchors,
            f"the Makefile no longer hands {sorted(missing_anchors)} to $(PYTHON). If a "
            "target was renamed or moved, move these anchors with it",
        )

    def test_every_script_a_recipe_runs_is_present(self):
        """The claim from the exported tree's side: nothing a target runs is missing.

        This is the assertion that was false in every public clone from the day
        `app-load-check` was written until issue #345.
        """
        absent = sorted(script for script in self.scripts if not (APP_DIR / script).is_file())
        self.assertFalse(
            absent,
            f"the Makefile runs {absent}, which is not in this tree. A target that cannot "
            "find its own script stops with an interpreter traceback, and the README goes "
            "on advertising it (issue #345)",
        )

    def test_no_recipe_reaches_into_a_directory_the_export_strips(self):
        """The claim from the private tree's side: no recipe points at content that cannot travel.

        Skipped rather than passed where the deny-list is absent, which is precisely the
        exported tree — there the test above is the one that has teeth.
        """
        if not DENY_LIST.is_file():
            self.skipTest(
                "no deny-list in this tree, so it is an exported one; "
                "test_every_script_a_recipe_runs_is_present carries the claim here"
            )

        rules = _stripped_directories()
        self.assertTrue(rules, f"{DENY_LIST} parsed to no rules, so this check reads nothing")

        offenders = {}
        for script in self.scripts:
            relative = (APP_DIR / script).resolve().relative_to(REPO_ROOT).as_posix()
            for rule in rules:
                if relative == rule or (rule.endswith("/") and relative.startswith(rule)):
                    offenders[script] = rule
        self.assertFalse(
            offenders,
            "the Makefile runs script(s) the public export strips, so these targets work "
            f"here and fail in every public clone: {offenders}. Move the script into "
            "streamlit_app/, or take the target out of the README that advertises it",
        )


if __name__ == "__main__":
    unittest.main()
