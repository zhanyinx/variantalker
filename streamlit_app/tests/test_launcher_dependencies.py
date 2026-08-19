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
* that one interpreter is a virtual environment inside the checkout — ``setup.sh`` builds
  it and the launcher resolves it through the same sourced file, so the two cannot name
  different directories, and installing MAFigate cannot disturb an interpreter the
  collaborator already depends on (issue #258);
* the launcher checks before it boots, because the app degrades quietly;
* ``check_dependencies.py`` actually parses ``requirements.txt`` — a checker that silently
  reads nothing would report a clean environment on a bare machine, which is the failure
  mode of every guard that cannot fire;
* and every line of that file pins one exact version, so the interpreter the two scripts
  agree on is filled with one known set of software rather than with whatever PyPI served
  that day (issue #256 — see :class:`TestEveryRequirementIsPinned`).

**Where the presence check itself is asserted.** Mostly not here. Whether *this* machine
has plotly is an environment fact, and it is made loud in the two places a human meets it:
``make check-deps`` and a launcher that refuses to start. A test that simply demanded a
full install would fail for anyone running the suite in the narrow environment the vendor
drift CI job uses, which is deliberately pytest and pandas only. The one exception is
:func:`test_an_interpreter_that_can_run_the_app_has_plotly_too`, which closes the exact
hole #162 fell through — see its own docstring.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

import pytest

import check_dependencies
from check_dependencies import declared_requirements, missing_requirements

APP_DIR = Path(check_dependencies.__file__).resolve().parent

#: The file that decides which interpreter both scripts use, sourced by both of them.
#: Named here as a path as well as read, because two guards below run it.
ENV_SCRIPT = APP_DIR / "mafigate_python.sh"

SETUP = (APP_DIR / "setup.sh").read_text()
LAUNCHER = (APP_DIR / "run_mafigate.sh").read_text()
MAKEFILE = (APP_DIR / "Makefile").read_text()
ENV = ENV_SCRIPT.read_text()

REQUIREMENTS = APP_DIR / "requirements.txt"

#: The commands that must never be invoked as a bare PATH entry, because each one exists
#: in more than one interpreter on a normal machine.
AMBIGUOUS_COMMANDS = ("pip", "pip3", "streamlit", "pytest")

#: Distributions ``requirements.txt`` is known to declare. Asserted by the pin guard as a
#: *parse anchor*, following ``test_vendor_drift.DMG_ANCHOR_SOURCES``: if the file's shape
#: changes so that :func:`_declaring_lines` stops recognising it, these fail loudly instead
#: of every pin assertion passing vacuously over an empty list. Three names rather than one
#: so that the anchor cannot be satisfied by a file that has lost all but its first line.
PIN_ANCHOR_DISTRIBUTIONS = {"streamlit", "pandas", "plotly"}

#: A line that pins one exact version: a distribution name, optionally an extras bracket,
#: ``==``, one version, nothing else. Matched over the *whole* line, which is what rejects
#: every other shape a specifier can take — a floor (``>=1.5.0``), a compatible release
#: (``~=1.5``), a wildcard (``==1.5.*``), a comma-separated range, a ``-r`` include, an
#: option line. Under any of those, two people installing the same MAFigate release a month
#: apart can end up running different software, which is the whole of issue #256.
#:
#: Deliberately tolerant of two shapes that *are* exact pins even though this file writes
#: neither: whitespace around the operator (``pandas == 2.3.1``) and an extras bracket
#: (``streamlit[snowflake]==1.47.0``), both PEP 508-legal. A guard that called those unpinned
#: would be making a false accusation, and its own message would tell the reader to write
#: what they had already written.
#:
#: An environment marker (``; python_version < "3.12"``) is refused, and that is a choice
#: rather than an oversight: the version is exact but *which* requirements apply stops being
#: a property of the file alone, so the reproducibility this guard exists for would hold only
#: per-platform. If MAFigate ever needs one, this is the line to argue with.
EXACT_PIN = re.compile(r"^[A-Za-z0-9._-]+(\[[A-Za-z0-9._,-]+\])?\s*==\s*[A-Za-z0-9._+!-]+$")

#: The development extras ``requirements.txt`` names as commented ``pip install`` hints
#: rather than as requirements. One tuple, read by the checker's test and by the pin guard's,
#: so the two cannot fall out of step about which names are hints.
DEVELOPMENT_EXTRAS = ("pytest", "pytest-cov", "black", "isort", "pre-commit")


def _is_declaring(raw_line):
    """Whether one line of ``requirements.txt`` declares a requirement rather than prose."""
    return bool(raw_line.split("#", 1)[0].strip())


def _declaring_lines(text):
    """Every line of ``requirements.txt`` that declares a requirement, comment dropped.

    The same reduction ``check_dependencies.declared_requirements`` performs — ``#`` and
    everything after it goes, then whitespace — but stopping *before* the distribution name
    is split off, because what the pin guard has to read is the specifier and not the name.
    That is also why this is a second parse rather than a call into the checker: the guard
    would be unfailable if it asked the module it guards what the file says. The two are held
    together by
    :meth:`TestEveryRequirementIsPinned.test_the_pin_guard_reads_what_the_checker_reads`,
    so a line shape that one of them stops recognising cannot leave the other quietly
    reading less than the whole file.
    """
    return [
        raw_line.split("#", 1)[0].strip()
        for raw_line in text.splitlines()
        if _is_declaring(raw_line)
    ]


def _distribution_name(specifier):
    """The distribution a specifier names — ``plotly==6.9.0`` is ``plotly``.

    A second, independent spelling of ``check_dependencies._NAME_ENDS``, and deliberately
    so, for the reason :func:`_declaring_lines` gives. This is the finer of the two copies —
    it is a character class rather than a scan — and what holds them together is the
    cross-check, which compares the names the two parses come away with rather than trusting
    either spelling.
    """
    return re.split(r"[<>=!~;\[,\s]", specifier, maxsplit=1)[0]


def _reason_for(prefix, text):
    """The comment that explains one requirement: its own, plus its continuation lines.

    A continuation is a comment-only line *indented* under the requirement, which is how
    this file writes a reason too long for one line. The indentation is what distinguishes
    one from a comment of the file's own — the ``openpyxl`` paragraph and the development
    extras both start at column zero — so a reason cannot be answered for by prose that
    belongs to the file rather than to the pin.
    """
    lines = text.splitlines()
    for index, line in enumerate(lines):
        if not line.startswith(prefix):
            continue
        reason = [line.partition("#")[2]]
        for continuation in lines[index + 1 :]:
            if not continuation[:1].isspace() or not continuation.strip().startswith("#"):
                break
            reason.append(continuation.strip().lstrip("#"))
        return "\n".join(reason)
    return ""


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


def _fake_app_dir(stubs=()):
    """A throwaway app directory holding the real env script and stub interpreters.

    The stubs are files, not real virtual environments, because what is under test here is
    the *resolution* — which path the two scripts end up calling — and a real venv would
    make these tests slow without making them stricter. The one test that does want a real
    venv builds one and says so.
    """
    directory = Path(tempfile.mkdtemp(prefix="mafigate-venv-"))
    shutil.copy(ENV_SCRIPT, directory / ENV_SCRIPT.name)
    for relative in stubs:
        stub = directory / relative
        stub.parent.mkdir(parents=True, exist_ok=True)
        stub.write_text("#!/bin/sh\nexit 0\n")
        stub.chmod(0o755)
    return directory


#: The two lines both scripts run before they do anything else, and nothing more. Run as a
#: script rather than asserted as text so that what is checked is the resolution's
#: *behaviour*: a comment claiming the venv is preferred cannot satisfy it.
_RESOLVE = (
    "set -eu\n"
    'PYTHON="${MAFIGATE_PYTHON:-python3}"\n'
    ". ./mafigate_python.sh\n"
    "if mafigate_use_venv; then promoted=yes; else promoted=no; fi\n"
    'printf "%s %s" "$promoted" "$PYTHON"\n'
)


def _resolve_python(app_dir, **environment):
    """``(did it use the venv, which interpreter)``, as the scripts themselves would get it."""
    env = {key: value for key, value in os.environ.items() if key != "MAFIGATE_PYTHON"}
    env.update(environment)
    completed = subprocess.run(
        ["bash", "-c", _RESOLVE],
        cwd=app_dir,
        env=env,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        raise AssertionError(
            f"resolving the interpreter exited {completed.returncode}: {completed.stderr!r}"
        )
    promoted, python = completed.stdout.split(" ", 1)
    return promoted == "yes", python


class TestTheOneInterpreterIsAVenvInTheCheckout(unittest.TestCase):
    """Issue #258: ``setup.sh`` builds the interpreter rather than borrowing one.

    The interpreter-agreement claim above is unchanged by this — one variable still decides
    where MAFigate installs and where it runs. What changes is what that variable defaults
    to: a virtual environment in the checkout, so a collaborator's first run of MAFigate
    cannot leave a dependency conflict in the python their other work needs.

    Two properties carry the weight, and both can fail:

    * the setup script and the launcher must resolve **the same** venv, which is asserted by
      requiring that neither of them names the directory at all — the path is written once,
      in the file both of them source;
    * ``MAFIGATE_PYTHON`` must still short-circuit the whole thing, because the venv is a
      default and not a cage. That one is asserted by *running* the resolution.
    """

    def test_setup_builds_the_venv_through_the_chosen_interpreter(self):
        """``-m venv``, for the same reason as ``-m pip``: PATH must not get a vote."""
        self.assertIn('"$PYTHON" -m venv "$MAFIGATE_VENV"', ENV)
        self.assertIn("mafigate_create_venv", "\n".join(_command_lines(SETUP)))

    def test_both_scripts_resolve_the_venv_through_the_same_file(self):
        for name, text in (("setup.sh", SETUP), ("run_mafigate.sh", LAUNCHER)):
            with self.subTest(script=name):
                self.assertIn(f". ./{ENV_SCRIPT.name}", _command_lines(text))

    def test_nothing_but_the_env_script_spells_the_venv_path(self):
        """Two copies of a path are two things that can drift — the shape of #162 again.

        The Makefile is held to this too, and that is not pedantry: it was the first thing
        to break the rule. It needs the same two facts the scripts need, and it first got
        them from its own ``$(wildcard $(CURDIR)/.venv/bin/python …)`` — a copy that would go
        on resolving happily after the environment moved, and point ``make test`` at a python
        ``./setup.sh`` had not filled. It reads them out of the env script now.
        """
        sources = (
            ("setup.sh", SETUP),
            ("run_mafigate.sh", LAUNCHER),
            ("Makefile", MAKEFILE),
        )
        for name, text in sources:
            for line in _command_lines(text):
                with self.subTest(source=name, line=line):
                    self.assertNotIn(
                        ".venv",
                        line,
                        f"{name} names the venv directory itself; it must read it from "
                        f"{ENV_SCRIPT.name}, which is the one place that decides where the "
                        "environment is",
                    )

    def test_the_venv_lives_in_the_checkout(self):
        """Not in the ``~/.mafigate`` the installers use: two checkouts would collide there."""
        self.assertIn('MAFIGATE_VENV=".venv"', ENV)
        commands = "\n".join(_command_lines(ENV))
        for outside in ("$HOME", "~/", ".mafigate"):
            with self.subTest(path=outside):
                self.assertNotIn(outside, commands)

    def test_the_env_script_invokes_no_ambiguous_command_by_itself(self):
        """The third script now on this path, held to the same rule as the other two."""
        for command in _invoked_commands(ENV):
            with self.subTest(command=command):
                self.assertNotIn(command, AMBIGUOUS_COMMANDS)

    def test_the_venv_is_git_ignored(self):
        """Asked of git, not of ``.gitignore``'s text, so a shadowed pattern still fails."""
        completed = subprocess.run(
            ["git", "check-ignore", "-q", "streamlit_app/.venv/bin/python"],
            cwd=APP_DIR.parent,
            capture_output=True,
        )
        if completed.returncode == 128:
            self.skipTest("not a git checkout, so there is no ignore rule to read")
        self.assertEqual(
            completed.returncode,
            0,
            "streamlit_app/.venv is not git-ignored, so setup.sh's own environment would "
            "show up as thousands of untracked files (issue #258)",
        )

    def test_the_resolution_prefers_the_venv_it_finds(self):
        app_dir = _fake_app_dir([".venv/bin/python"])
        try:
            promoted, python = _resolve_python(app_dir)
        finally:
            shutil.rmtree(app_dir)
        self.assertTrue(promoted)
        self.assertEqual(python, ".venv/bin/python")

    def test_the_resolution_finds_the_windows_layout_too(self):
        """``python3 -m venv`` writes ``Scripts/python.exe`` when the bash is Git Bash."""
        app_dir = _fake_app_dir([".venv/Scripts/python.exe"])
        try:
            promoted, python = _resolve_python(app_dir)
        finally:
            shutil.rmtree(app_dir)
        self.assertTrue(promoted)
        self.assertEqual(python, ".venv/Scripts/python.exe")

    def test_no_venv_leaves_the_interpreter_alone(self):
        """So the launcher can say so and let the dependency check refuse, rather than
        inventing a path that does not exist."""
        app_dir = _fake_app_dir()
        try:
            promoted, python = _resolve_python(app_dir)
        finally:
            shutil.rmtree(app_dir)
        self.assertFalse(promoted)
        self.assertEqual(python, "python3")

    def test_an_explicit_interpreter_short_circuits_the_venv(self):
        """Measured with a venv present, which is the only case that can fail."""
        app_dir = _fake_app_dir([".venv/bin/python"])
        try:
            promoted, python = _resolve_python(
                app_dir, MAFIGATE_PYTHON="/somewhere/else/bin/python"
            )
        finally:
            shutil.rmtree(app_dir)
        self.assertFalse(promoted)
        self.assertEqual(python, "/somewhere/else/bin/python")

    def test_a_dangling_venv_interpreter_is_not_used(self):
        """A venv outlives the python that built it — a Homebrew upgrade is enough.

        The symlink is then still there and points at nothing, so anything that handed the
        launcher this path would hand it an interpreter that cannot start. setup.sh rebuilds
        in place from here, which is what makes re-running it the repair rather than a
        ``rm -rf``.
        """
        app_dir = _fake_app_dir()
        (app_dir / ".venv" / "bin").mkdir(parents=True)
        (app_dir / ".venv" / "bin" / "python").symlink_to("/nonexistent/python")
        try:
            promoted, python = _resolve_python(app_dir)
        finally:
            shutil.rmtree(app_dir)
        self.assertFalse(promoted)
        self.assertEqual(python, "python3")

    def test_a_venv_interpreter_that_cannot_be_run_is_not_used(self):
        """The case that separates *executable* from *present*, which is why -x is the test.

        The dangling symlink above does not: ``[ -e ]`` follows links, so it rejects that
        one too, and mutating -x to -e was measured surviving on it alone. A file that is
        there but not runnable — a checkout that lost its permission bits, a directory that
        happens to be called ``python`` — is the case only -x rejects, and a resolution that
        accepted it would report a working environment and then fail at the first launch.
        """
        app_dir = _fake_app_dir()
        stub = app_dir / ".venv" / "bin" / "python"
        stub.parent.mkdir(parents=True)
        stub.write_text("#!/bin/sh\nexit 0\n")
        stub.chmod(0o644)
        try:
            promoted, python = _resolve_python(app_dir)
        finally:
            shutil.rmtree(app_dir)
        self.assertFalse(promoted)
        self.assertEqual(python, "python3")

    def test_the_makefile_reaches_the_same_venv(self):
        """``make test`` must not test an interpreter ``./setup.sh`` never filled.

        Read by asking make what it would run — ``-n`` prints the recipe with ``$(PYTHON)``
        expanded — rather than by matching the variable's definition, because what matters
        is the answer and there are three ways to spell the question wrongly. All three
        cases are checked from one copied Makefile, with a stub environment beside it.
        """
        app_dir = _fake_app_dir([".venv/bin/python"])
        shutil.copy(APP_DIR / "Makefile", app_dir / "Makefile")
        try:
            resolved = {}
            for case, overrides in (
                ("venv", {}),
                ("override", {"MAFIGATE_PYTHON": "/somewhere/else/bin/python"}),
            ):
                env = {k: v for k, v in os.environ.items() if k != "MAFIGATE_PYTHON"}
                env.update(overrides)
                completed = subprocess.run(
                    ["make", "-n", "check-deps"],
                    cwd=app_dir,
                    env=env,
                    capture_output=True,
                    text=True,
                )
                self.assertEqual(completed.returncode, 0, completed.stderr)
                resolved[case] = completed.stdout
        finally:
            shutil.rmtree(app_dir)

        self.assertIn(
            f"{app_dir}/.venv/bin/python check_dependencies.py",
            resolved["venv"],
            "the Makefile ignores the virtual environment ./setup.sh builds, so `make "
            "install` and `make test` would work on a different interpreter from the two "
            f"scripts (issue #258). It would run: {resolved['venv']!r}",
        )
        self.assertIn(
            "/somewhere/else/bin/python check_dependencies.py",
            resolved["override"],
            f"MAFIGATE_PYTHON does not reach the Makefile: {resolved['override']!r}",
        )

    @pytest.mark.slow
    def test_the_venv_it_builds_is_the_one_it_then_installs_into(self):
        """The end-to-end property, minus the network: build, resolve, ask the interpreter.

        ``sys.prefix`` is the interpreter's own answer to "which environment am I", so this
        fails if the venv is built in one place and resolved from another — the coincidence
        the two guards above forbid by construction, checked here by consequence. It also
        pins the half of "writes nothing outside the checkout" that needs no network: the
        environment lands under the app directory.
        """
        app_dir = _fake_app_dir()
        script = (
            "set -eu\n"
            'PYTHON="${MAFIGATE_PYTHON:-python3}"\n'
            ". ./mafigate_python.sh\n"
            "mafigate_create_venv\n"
            "mafigate_use_venv\n"
            '"$PYTHON" -c "import sys; print(sys.prefix)"\n'
        )
        env = {key: value for key, value in os.environ.items() if key != "MAFIGATE_PYTHON"}
        try:
            completed = subprocess.run(
                ["bash", "-c", script],
                cwd=app_dir,
                env=env,
                capture_output=True,
                text=True,
            )
            if completed.returncode != 0:
                if "ensurepip" in completed.stderr or "venv" in completed.stderr:
                    self.skipTest(
                        "this python cannot create a venv "
                        f"(the failure setup.sh reports): {completed.stderr.strip()!r}"
                    )
                raise AssertionError(f"building the venv failed: {completed.stderr!r}")
            prefix = Path(completed.stdout.strip()).resolve()
            expected = (app_dir / ".venv").resolve()
        finally:
            shutil.rmtree(app_dir)
        self.assertEqual(prefix, expected)


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
        for extra in DEVELOPMENT_EXTRAS:
            with self.subTest(extra=extra):
                self.assertNotIn(extra, declared)

    def test_it_takes_the_name_without_the_version_pin(self):
        """``plotly==6.9.0`` names the distribution ``plotly``.

        Written against whatever version stands there rather than against the literal
        string, so that bumping a pin does not have to come here: the assertion is that the
        line carries a specifier *and* that the checker hands back the bare name, which is
        what ``importlib.metadata`` can be asked about. The first half is what keeps the
        second from passing over a file whose specifiers have all gone.
        """
        declaring = _declaring_lines(REQUIREMENTS.read_text())
        pinned = [line for line in declaring if line.startswith("plotly==")]
        self.assertEqual(
            len(pinned),
            1,
            f"requirements.txt does not pin plotly with `==`: {declaring}",
        )
        self.assertIn("plotly", declared_requirements())
        self.assertNotIn(pinned[0], declared_requirements())

    def test_it_reports_a_distribution_that_is_not_installed(self):
        absent = "mafigate-nothing-is-called-this"
        self.assertEqual(missing_requirements([absent]), [absent])

    def test_it_does_not_report_one_that_is(self):
        """Checked against pandas, which the suite cannot run without in any case."""
        self.assertEqual(missing_requirements(["pandas"]), [])

    def test_the_requirements_file_is_read_from_beside_the_checker(self):
        """Not from the working directory, so ``make -C streamlit_app`` sees the same file."""
        self.assertEqual(check_dependencies.REQUIREMENTS.parent, APP_DIR)


class TestEveryRequirementIsPinned(unittest.TestCase):
    """One release, one set of software — asserted line by line.

    Issue #256. ``requirements.txt`` declared floors (``streamlit>=1.47.0``,
    ``pandas>=1.5.0``) and every channel — the installers and the clone route alike — builds
    its environment from this one file at install time. So the venv was whatever PyPI served
    that day: two people installing the same MAFigate release a month apart got different
    software behind the same version number, in a tool that filters clinical variants. A
    version number that does not identify what is running is not a version number.

    The parse is proved before anything is asserted, because this is a whole-file sweep and
    a sweep that reads nothing reports every file clean. Two independent proofs, and neither
    is a count against a threshold: the anchor names three distributions the file is known
    to carry, and the cross-check demands the pin guard and ``check_dependencies.py`` come
    away from the file with the same set of names. Prior art is the vendor-drift suite, which
    fails when it cannot parse the DMG copy list rather than passing over an empty one.
    """

    def setUp(self):
        self.text = REQUIREMENTS.read_text()
        self.declaring = _declaring_lines(self.text)

    def _assert_the_parse_worked(self):
        """The anchor, asserted by each test *before* it sweeps anything.

        Called from every test in this class rather than left to the sibling below, and the
        distinction matters: a sweep whose parse returned nothing runs zero iterations and
        goes green on its own, so a class-level proof living in another test leaves each
        sweep individually satisfiable by an empty file. Only the whole class failing tells
        you that, and a single test is what a person reads. This is the vendor-drift
        suite's shape — ``test_both_installers_ship_the_vendor_package`` asserts
        ``MAFigate.py`` is in the parsed copy list inside the same test that then checks
        what the list contains, rather than trusting a neighbour to have checked.
        """
        parsed = {_distribution_name(line) for line in self.declaring}
        self.assertLessEqual(
            PIN_ANCHOR_DISTRIBUTIONS,
            parsed,
            "could not parse the requirement list out of requirements.txt — the file's "
            f"shape changed and this guard is no longer reading it. Parsed: {sorted(parsed)}",
        )
        return parsed

    def test_the_pin_guard_reads_what_the_checker_reads(self):
        """The parse proof, plus the subtler half of it.

        The anchor catches a file this function has stopped understanding. The cross-check
        against ``declared_requirements`` catches the case where both parses still work but
        see different lines — which is how a requirement could be pinned as far as one of
        them can tell and absent as far as the other can.

        Said plainly, because it is easy to overrate: the two parses apply the same
        comment-and-whitespace reduction, so today the cross-check is close to a tautology
        and the anchor is doing the work. It earns its place only if one of the two parses is
        ever taught something the other is not — which is exactly when a copy diverges, and
        exactly the day nobody would think to look.
        """
        parsed = self._assert_the_parse_worked()
        self.assertEqual(
            parsed,
            set(declared_requirements()),
            "the pin guard and check_dependencies.py disagree about which lines of "
            "requirements.txt declare a requirement, so one of them is reading less than "
            "the whole file",
        )

    def test_every_runtime_requirement_pins_an_exact_version(self):
        """The claim itself: no floors, no ranges, no wildcards, no markers."""
        self._assert_the_parse_worked()
        for specifier in self.declaring:
            with self.subTest(requirement=specifier):
                self.assertRegex(
                    specifier,
                    EXACT_PIN,
                    f"requirements.txt declares `{specifier}`, which this guard cannot read "
                    "as one exact version — either it is looser than `==` (a floor, a range, "
                    "a wildcard) or it carries an environment marker, which makes the "
                    "resolved set depend on the machine rather than on this file. Every "
                    "channel builds its environment from here at install time, so either way "
                    "two people installing the same MAFigate release can get different "
                    "software (issue #256). Write it as `name==version` — extras and spaces "
                    "around the operator are fine — and keep the comment saying why the pin "
                    "sits where it does.",
                )

    def test_the_development_extras_are_still_only_comments(self):
        """This ticket pins what ships, not what a contributor chooses to install.

        ``pytest`` and friends are written here as commented ``pip install`` hints. Pinning
        them would make them requirements, and the app runs without them — a checker that
        demanded them would refuse to launch a perfectly good install.

        Both halves are the claim, and neither stands alone: that a name is not a requirement
        is satisfied by a file that never mentions it, and that a name is in the file is
        satisfied by a file that declares it outright. The sibling
        :meth:`TestTheCheckerReadsTheRequirements.test_it_ignores_the_commented_development_extras`
        makes the same not-a-requirement point about the *checker's* reading; this one is
        about the file, and the two share :data:`DEVELOPMENT_EXTRAS` so they cannot disagree
        about which names are hints.
        """
        declared = self._assert_the_parse_worked()
        for extra in DEVELOPMENT_EXTRAS:
            with self.subTest(extra=extra):
                self.assertIn(
                    extra,
                    self.text,
                    f"the `pip install {extra}` hint has gone from requirements.txt, so "
                    "this test can no longer tell a comment from a requirement",
                )
                self.assertNotIn(extra, declared)

    def test_every_pin_still_says_why_it_is_there(self):
        """A pin with no reason beside it is a number nobody can safely change.

        The comments are the only record of *why* each version was chosen, and a pin makes
        them load-bearing in a way a floor did not: whoever bumps one has to know what the
        old number was protecting. Asserted as "there is a comment", not as particular
        prose, so the wording stays the author's — with one exception below, for the one
        rationale that is a list of specific features.
        """
        self._assert_the_parse_worked()
        for raw_line in self.text.splitlines():
            if not _is_declaring(raw_line):
                continue
            with self.subTest(line=raw_line.strip()):
                _, _, comment = raw_line.partition("#")
                self.assertTrue(
                    comment.strip(),
                    f"`{raw_line.strip()}` pins a version with no comment saying why, so "
                    "the next person to bump it cannot know what the number was for",
                )

    def test_the_streamlit_pin_still_names_the_features_that_set_its_floor(self):
        """1.40 and 1.46 are why the Streamlit pin cannot move downwards.

        The data page's section switch is an ``st.segmented_control`` (1.40) sized with
        ``width=`` (1.46). Those two numbers were the floor's justification and they are now
        the pin's: a bump is safe above them and breaks the page below them. They are the
        one piece of this file's prose worth asserting, because they are the only part that
        a future operator has to *check* rather than merely read.

        Read off the Streamlit pin's *own* comment block rather than off the whole file, so
        that the two numbers cannot be answered for by a mention somewhere else in it — the
        header, say, which lists the floors this file used to carry.
        """
        reason = _reason_for("streamlit==", self.text)
        self.assertIn(
            "segmented_control",
            reason,
            "could not find the Streamlit pin's comment block — the file's shape changed "
            f"and this test is no longer reading it. Read: {reason!r}",
        )
        for feature_floor in ("1.40", "1.46"):
            with self.subTest(floor=feature_floor):
                self.assertIn(
                    feature_floor,
                    reason,
                    f"requirements.txt no longer records that {feature_floor} is a Streamlit "
                    "feature floor, so the pin above it looks arbitrary and a bump cannot "
                    "be checked against anything",
                )


class TestThePinsAndThePythonFloorAgree(unittest.TestCase):
    """The pins decide which pythons can install MAFigate, so they decide what to promise.

    A consequence of issue #256 that is easy to miss, because nothing fails on the machine
    where the pinning is done: ``numpy==2.3.1`` declares Requires-Python ``>=3.11`` and
    ``streamlit-aggrid==1.1.6`` declares ``>=3.10``, so pip refuses the file outright on 3.9
    and 3.10. Under the old floors it would instead have resolved *something* older and the
    app would have started. The README promised 3.9 or newer at the time, and that sentence
    is what a person following the clone route reads before they find out.

    So the floor is a claim in two places at once — the pins that enforce it and the prose
    that states it — and this holds them together. Kept in this module rather than in
    ``test_delivery_channels_copy.py``, which owns the README's channel promises: the fact
    being asserted is a property of ``requirements.txt``, and it belongs beside the guard that
    made it true.
    """

    README = (APP_DIR / "README.md").read_text()

    #: The minor the pinned set requires, and the one both installers bundle. A literal,
    #: because the alternative — reading Requires-Python off the installed distributions —
    #: can only run in an environment that already has them, which is exactly the narrow
    #: CI job where it would skip in silence.
    PYTHON_FLOOR = "3.11"

    #: What the floor used to be, swept as a literal so that a half-finished bump fails.
    SUPERSEDED_FLOOR = "3.9"

    def test_the_requirements_file_states_the_floor_its_pins_impose(self):
        text = REQUIREMENTS.read_text()
        self.assertIn(
            self.PYTHON_FLOOR,
            text,
            "requirements.txt does not say which Python its pins require, so nobody "
            "changing a pin can tell whether they have moved the floor (issue #256)",
        )

    def test_the_readme_promises_the_floor_the_pins_actually_impose(self):
        """Read off the prerequisites block, which is the sentence the clone route shows."""
        blocks = [
            line for line in self.README.splitlines() if line.startswith("> **Prerequisites.")
        ]
        self.assertTrue(
            blocks,
            "the README has no prerequisites blockquote, so this test is reading nothing",
        )
        prerequisites = "\n".join(blocks)
        self.assertIn(
            self.PYTHON_FLOOR,
            prerequisites,
            f"the README's prerequisites must name Python {self.PYTHON_FLOOR}, which is what "
            f"the pins require — it reads: {prerequisites!r}",
        )
        self.assertNotIn(
            self.SUPERSEDED_FLOOR,
            prerequisites,
            f"the README still offers Python {self.SUPERSEDED_FLOOR}, on which "
            "`pip install -r requirements.txt` now fails outright rather than installing "
            f"anything — it reads: {prerequisites!r}",
        )


class TestThePinPatternItself(unittest.TestCase):
    """:data:`EXACT_PIN` against the shapes a requirements line can take.

    The guard above is only as good as this pattern, and a pattern is the kind of thing that
    is wrong in one direction while looking right in the other. Both directions are tabled
    here, because each catches a different mistake: a pattern too loose lets a floor through
    and the whole of issue #256 comes back, and a pattern too strict calls a legal exact pin
    unpinned — telling the reader, in the failure message, to write what they had written.

    Written against the pattern directly rather than by mutating the real file, so the cases
    that this file does not use (extras, whitespace, epochs, local versions, a ``-r``
    include) are covered without asking anyone to put them there.
    """

    #: Exact pins, including three shapes ``requirements.txt`` does not currently write:
    #: whitespace around the operator, an extras bracket, and PEP 440's epoch and local
    #: version forms. All are one exact version and must not be called otherwise.
    PINNED = (
        "pandas==2.3.1",
        "pandas == 2.3.1",
        "pandas ==2.3.1",
        "streamlit[snowflake]==1.47.0",
        "streamlit-aggrid==1.1.6",
        "PyYAML==6.0.2",
        "foo==1.0.0rc1",
        "foo==1!2.0",
        "foo==1.0.0+local",
    )

    #: Everything that leaves the installed version up to the day, the machine, or another
    #: file — plus ``===``, which is arbitrary equality rather than version equality, and a
    #: bare name, which pins nothing at all.
    NOT_PINNED = (
        "pandas>=1.5.0",
        "pandas>1.5",
        "pandas<=3",
        "pandas~=2.3",
        "pandas!=2.3.0",
        "pandas==2.3.*",
        "pandas===2.3.1",
        "pandas>=1.5,<3",
        "pandas==2.3.1,!=2.3.0",
        'pandas==2.3.1; python_version < "3.12"',
        "-r base.txt",
        "--index-url https://example.invalid",
        "pandas",
        "git+https://example.invalid/pandas.git",
    )

    def test_it_accepts_every_exact_pin(self):
        for specifier in self.PINNED:
            with self.subTest(specifier=specifier):
                self.assertRegex(
                    specifier,
                    EXACT_PIN,
                    f"`{specifier}` pins one exact version and the guard would refuse it, "
                    "so it would report a correctly pinned file as unpinned",
                )

    def test_it_refuses_everything_that_is_not_one(self):
        for specifier in self.NOT_PINNED:
            with self.subTest(specifier=specifier):
                self.assertNotRegex(
                    specifier,
                    EXACT_PIN,
                    f"`{specifier}` does not fix one version, and the guard would accept it "
                    "— so requirements.txt could go back to resolving at install time with "
                    "the suite still green (issue #256)",
                )


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
