"""Installing and launching MAFigate must land in the same interpreter.

Issue #162: the app's Summary tab drew six wrong charts for as long as ``plotly`` was
declared in ``requirements.txt`` and installed in no interpreter on the machine. The
declaration was never the gap. The gap was that ``setup.sh`` installed with a bare
``pip`` and ``run_mafigate.sh`` launched a bare ``streamlit``, and those are two separate
entries on ``PATH`` that need not belong to the same python. On the machine where this
was measured they did not: ``pip`` was miniconda's and ``streamlit`` was Homebrew's, so
every declared dependency installed successfully into an interpreter the app never ran
in, and nothing anywhere failed.

So the guards here are about *agreement*, not about any one package:

* both scripts, and the Makefile, go through one interpreter variable, and reach pip,
  streamlit and pytest as ``-m`` modules of it rather than as PATH entries of their own;
* the launcher checks before it boots, because the app degrades quietly;
* ``check_dependencies.py`` actually parses ``requirements.txt`` — a checker that silently
  reads nothing would report a clean environment on a bare machine, which is the failure
  mode of every guard that cannot fire.

**Where the presence check itself is asserted.** Mostly not here. Whether *this* machine
has plotly is an environment fact, and it is made loud in the two places a human meets it:
``make check-deps`` and a launcher that refuses to start. A test that simply demanded a
full install would fail for anyone running the suite in the narrow environment the vendor
drift CI job uses, which is deliberately pytest and pandas only. The one exception is
:func:`test_an_interpreter_that_can_run_the_app_has_plotly_too`, which closes the exact
hole #162 fell through — see its own docstring.
"""

from __future__ import annotations

import re
import unittest
from pathlib import Path

import check_dependencies
from check_dependencies import declared_requirements, missing_requirements

APP_DIR = Path(check_dependencies.__file__).resolve().parent

SETUP = (APP_DIR / "setup.sh").read_text()
LAUNCHER = (APP_DIR / "run_mafigate.sh").read_text()
MAKEFILE = (APP_DIR / "Makefile").read_text()

#: The commands that must never be invoked as a bare PATH entry, because each one exists
#: in more than one interpreter on a normal machine.
AMBIGUOUS_COMMANDS = ("pip", "pip3", "streamlit", "pytest")


def _command_lines(text):
    """The lines of a script that run something, comments dropped.

    Both of this module's "does the script do X" guards read these rather than the raw
    file, because both scripts *describe* what they do in a header comment as well as
    doing it. Two early drafts of these tests read the description by accident: one
    ordered the install and the verify by their position in the whole text and failed on a
    correct script, and one passed while a mutation deleted the launcher's dependency
    check outright, because the filename still appeared in a comment two lines up.
    """
    lines = []
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if line and not line.startswith("#"):
            lines.append(line)
    return lines


def _invoked_commands(text, recipe_lines_only=False):
    """The first word of every line that runs something, echoes excluded.

    Deliberately crude, and crude in the safe direction: it reads the *first* word, so a
    line that mentions ``pip install`` inside an ``echo`` — which both scripts and the
    Makefile's ``deps-check`` help text do — is not mistaken for an invocation of it.
    """
    commands = []
    for raw_line in text.splitlines():
        if recipe_lines_only and not raw_line.startswith("\t"):
            continue
        line = raw_line.strip()
        # Makefile recipe prefixes, then shell noise.
        line = line.lstrip("@-+")
        if not line or line.startswith("#"):
            continue
        first = line.split()[0]
        if first in ("echo", "if", "elif", "else", "fi", "then", "exit", "cd", "set"):
            continue
        if first == "exec":
            parts = line.split()
            first = parts[1] if len(parts) > 1 else first
        commands.append(first.strip('"'))
    return commands


class TestOneInterpreter(unittest.TestCase):
    """Setup, launch and the suite all reach the same python."""

    def test_setup_installs_as_a_module_of_the_chosen_interpreter(self):
        self.assertIn('"$PYTHON" -m pip install -r requirements.txt', SETUP)

    def test_the_launcher_starts_streamlit_as_a_module_of_the_chosen_interpreter(self):
        """``-m streamlit`` is what makes the checked interpreter the running one.

        A bare ``streamlit`` is resolved by PATH, so it can be a different python from the
        one that was just checked — which is precisely how an app with no plotly passed a
        launcher's dependency check in #162.
        """
        self.assertIn('exec "$PYTHON" -m streamlit run MAFigate.py', LAUNCHER)

    def test_both_scripts_read_the_same_variable_with_the_same_default(self):
        """One variable, or the two scripts can disagree about where the app lives."""
        for name, text in (("setup.sh", SETUP), ("run_mafigate.sh", LAUNCHER)):
            with self.subTest(script=name):
                self.assertIn('PYTHON="${MAFIGATE_PYTHON:-python3}"', text)

    def test_neither_script_invokes_an_ambiguous_command_by_itself(self):
        for name, text in (("setup.sh", SETUP), ("run_mafigate.sh", LAUNCHER)):
            for command in _invoked_commands(text):
                with self.subTest(script=name, command=command):
                    self.assertNotIn(
                        command,
                        AMBIGUOUS_COMMANDS,
                        f"{name} runs a bare `{command}`, which PATH resolves — possibly "
                        "to a different python than the one running the app (issue #162)",
                    )

    def test_no_makefile_recipe_invokes_an_ambiguous_command_by_itself(self):
        """``make install`` filled miniconda while ``make run`` launched Homebrew."""
        for command in _invoked_commands(MAKEFILE, recipe_lines_only=True):
            with self.subTest(command=command):
                self.assertNotIn(
                    command,
                    AMBIGUOUS_COMMANDS,
                    f"a Makefile recipe runs a bare `{command}`; use $(PYTHON) -m "
                    f"{command} so installing, launching and testing agree (issue #162)",
                )

    def test_the_makefile_installs_and_launches_through_the_python_variable(self):
        self.assertIn("$(PYTHON) -m pip install -r requirements.txt", MAKEFILE)
        self.assertIn("$(PYTHON) -m streamlit run MAFigate.py", MAKEFILE)


class TestTheLauncherChecksBeforeItBoots(unittest.TestCase):
    """The app degrades quietly, so the check cannot wait for the chart."""

    def test_the_launcher_runs_the_dependency_check(self):
        """Asserted over the launcher's *commands*, not its text — see _command_lines."""
        invocations = [
            line for line in _command_lines(LAUNCHER) if "check_dependencies.py" in line
        ]
        self.assertTrue(
            invocations,
            "run_mafigate.sh does not run check_dependencies.py, so it can boot an app "
            "whose charts are silently broken (issue #162)",
        )

    def test_setup_verifies_that_what_it_asked_for_arrived(self):
        """A pip install that reports success is not evidence the app can run.

        Ordered over the *commands* rather than over the file's text: both filenames are
        named in the header comment too, in the other order, and an earlier draft of this
        test read those instead and failed on a correct script.
        """
        commands = _command_lines(SETUP)
        installs = [i for i, line in enumerate(commands) if "-m pip install" in line]
        verifies = [i for i, line in enumerate(commands) if "check_dependencies.py" in line]
        self.assertTrue(installs, "setup.sh installs nothing")
        self.assertTrue(verifies, "setup.sh does not verify the install")
        self.assertLess(installs[0], verifies[0], "setup.sh must verify *after* installing")

    def test_the_launcher_refusal_can_be_overridden(self):
        """A launcher nobody can override is its own kind of trap."""
        self.assertIn("MAFIGATE_SKIP_DEPS_CHECK", LAUNCHER)

    def test_the_launcher_serves_on_loopback_only(self):
        """Who can reach MAFigate, asserted where the answer is written.

        ``0.0.0.0`` — every interface — stood here until issue #182 and was never a decision:
        it left the developer script publishing an app that opens patient data to anyone who
        could route to the machine, while both packaged launchers bound loopback. The flag is
        asserted rather than trusted because nothing else in the suite reads this file, and a
        launcher's reach is not visible from inside the app it starts.

        Read through :func:`_command_lines`, not the raw file, for the reason that function's
        own docstring gives — and this test proved it the hard way: the first draft asserted
        ``0.0.0.0`` was absent from ``LAUNCHER`` and failed on a correct script, because the
        comment explaining the change quotes the address it removed.

        The upload half of this test is unchanged and unrelated: ``tests/test_help_claims.py``
        reasons that 200 MB is the real cap because nothing here raises
        ``server.maxUploadSize``, and that stays a Help-page claim this line can falsify.
        """
        commands = "\n".join(_command_lines(LAUNCHER))
        self.assertIn("--server.address 127.0.0.1", commands)
        self.assertNotIn("0.0.0.0", commands)
        self.assertNotIn("maxUploadSize", LAUNCHER)
        self.assertNotIn("maxUploadSize", MAKEFILE)


class TestTheCheckerReadsTheRequirements(unittest.TestCase):
    """A checker that reads nothing reports every environment clean."""

    def test_it_finds_every_runtime_requirement(self):
        """The exact set, so a parser that drops a line fails rather than passing quietly."""
        self.assertEqual(
            sorted(declared_requirements()),
            ["PyYAML", "numpy", "pandas", "plotly", "streamlit", "streamlit-aggrid"],
        )

    def test_it_ignores_the_commented_development_extras(self):
        """``pytest``, ``black`` and friends are ``pip install`` hints, not requirements.

        The app runs without them, so a checker that demanded them would refuse to launch
        a perfectly good install.
        """
        declared = declared_requirements()
        for extra in ("pytest", "pytest-cov", "black", "isort", "pre-commit"):
            with self.subTest(extra=extra):
                self.assertNotIn(extra, declared)

    def test_it_takes_the_name_without_the_version_floor(self):
        """``plotly>=5.15.0`` names the distribution ``plotly``."""
        self.assertIn("plotly>=5.15.0", (APP_DIR / "requirements.txt").read_text())
        self.assertIn("plotly", declared_requirements())

    def test_it_reports_a_distribution_that_is_not_installed(self):
        absent = "mafigate-nothing-is-called-this"
        self.assertEqual(missing_requirements([absent]), [absent])

    def test_it_does_not_report_one_that_is(self):
        """Checked against pandas, which the suite cannot run without in any case."""
        self.assertEqual(missing_requirements(["pandas"]), [])

    def test_the_requirements_file_is_read_from_beside_the_checker(self):
        """Not from the working directory, so ``make -C streamlit_app`` sees the same file."""
        self.assertEqual(check_dependencies.REQUIREMENTS.parent, APP_DIR)


class TestTheSilentSkipIsClosed(unittest.TestCase):
    """The one environment assertion, scoped to where it cannot be noise."""

    def test_an_interpreter_that_can_run_the_app_has_plotly_too(self):
        """plotly was the only declared dependency that could vanish without a failure.

        Every other one is imported at module scope somewhere the suite touches, so its
        absence is a collection error nobody can miss. ``plotly`` is imported behind a
        ``try`` in ``components/charts.py`` and behind a ``skipTest`` in
        ``tests/test_components.py``, so its absence read as six passing skips — and did,
        on the developer's own machine, for every run before #162.

        Asserted only where the app could actually run, which keeps it out of the vendor
        drift CI job's deliberately narrow environment (pytest and pandas only). In a full
        app environment, a missing plotly is now one loud failure instead of six quiet
        skips.
        """
        others = [name for name in declared_requirements() if name != "plotly"]
        if missing_requirements(others):
            self.skipTest(
                "this interpreter cannot run MAFigate at all "
                f"(missing {', '.join(missing_requirements(others))}), "
                "so plotly's absence is not the interesting fact here"
            )

        self.assertEqual(
            missing_requirements(["plotly"]),
            [],
            "this interpreter can run MAFigate but has no plotly, so all six Summary "
            "charts will fall back to st.bar_chart and the VAF chart's x-axis will read "
            "{'left': ..., 'right': ...} (issue #162). Run ./setup.sh.",
        )


class TestTheDocumentedSkipStillMatchesTheTests(unittest.TestCase):
    """The chart class's note says how many tests it hands on; it drifted by one."""

    def test_the_chart_class_counts_its_own_tests_correctly(self):
        """Written down because the note was stale: it said five, and there are six.

        The count is the whole point of that paragraph — "whoever widens CI to the full
        suite inherits six live tests" — so a sixth test arriving without the sentence
        moving makes the note misreport what it exists to report.
        """
        from tests import test_components

        chart_tests = [
            name for name in dir(test_components.TestChartLayout) if name.startswith("test")
        ]
        self.assertEqual(len(chart_tests), 6)

        note = test_components.TestChartLayout.__doc__
        self.assertIsNotNone(note)
        counted = re.search(r"inherits (\w+) live tests", note)
        self.assertIsNotNone(counted, "the chart class no longer states how many tests it has")
        self.assertEqual(counted.group(1), "six")


if __name__ == "__main__":
    unittest.main()
