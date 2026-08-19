"""What the Windows installer inputs must say (issue #261).

Two claims the .exe makes to a clinician, neither of which any other guard can see,
because both live in files no Python ever imports:

* **It brings its own interpreter.** Windows used to demand a system Python and check for
  it *softly*: setup popped *"Python 3 was not detected… continue anyway?"* and installing
  proceeded on Yes. So the install could complete, report success, and leave an app that
  fails on first launch — the worst shape a prerequisite can have. macOS already bundles a
  portable CPython; Windows now bundles the same pinned release, and the launcher stops
  looking anywhere else.
* **Uninstalling does not destroy the user's work.** ``[UninstallDelete]`` used to name
  ``%USERPROFILE%\\.mafigate`` — the whole directory, which is not just the venv. It also
  holds the parameter cache, so uninstalling silently deleted a clinician's saved filter
  setups, with no dialog saying so.

Why these are static text assertions
------------------------------------
Nothing here can be exercised: an Inno Setup script is compiled by a Windows-only
compiler, and ``launch.bat`` runs under ``cmd.exe``. So the guards read the inputs. That
makes them worth exactly as much as their parsing, which is why every one of them first
proves it found the structure it means to check — the vendor-drift suite's rule, and the
same reason :func:`_iss_sections` fails loudly rather than returning an empty section that
would satisfy every "does not contain" assertion below.

Each assertion here was shown failing against the files as they stood before this ticket:
the uninstall rule named the parent, the launcher probed ``PATH`` and sent the user to
python.org, the setup hook popped the soft prompt, and no ``Source:`` line shipped an
interpreter. A guard that is green on arrival has parsed the wrong thing.

**What is not asserted here.** That the .exe installs and launches on a machine with no
system Python. That is not unit-testable at all — it needs a clean Windows box — and the
PR that closes #261 records the hand-verification instead.
"""

from __future__ import annotations

import re
import subprocess
import unittest
from pathlib import Path

import pytest

STREAMLIT_APP = Path(__file__).resolve().parents[1]
BUILD_DIR = STREAMLIT_APP / "build"
WINDOWS_DIR = BUILD_DIR / "windows"

ISS_SCRIPT = WINDOWS_DIR / "installer.iss"
LAUNCH_BAT = WINDOWS_DIR / "launch.bat"
BUILD_BAT = WINDOWS_DIR / "build_installer.bat"
DMG_SCRIPT = BUILD_DIR / "mac" / "build_dmg.sh"

#: The one place the bundled interpreter's release is written down. Both builds read it, so
#: they cannot drift to different interpreters — which is the failure this file exists to
#: make impossible, since a Windows build pinning a different CPython from the DMG would
#: look perfectly healthy on both machines while shipping two different applications.
PIN_FILE = BUILD_DIR / "python_release.env"

#: Gated on the *pre-existing* inputs only. Deliberately not on :data:`PIN_FILE`: a gate
#: that included it would have turned every assertion below into a skip on the tree this
#: ticket started from, which is the vacuous-green failure the whole module is written
#: against. ``build/`` is stripped from the packaged app, so the skip is real there.
needs_build_scripts = pytest.mark.skipif(
    not (ISS_SCRIPT.is_file() and LAUNCH_BAT.is_file() and BUILD_BAT.is_file()),
    reason="installer scripts not present (build/ is stripped from the packaged app)",
)


# ---------------------------------------------------------------------------
# Reading the inputs, and failing when they cannot be read
# ---------------------------------------------------------------------------


def _iss_sections(iss: str) -> dict[str, list[str]]:
    """The Inno Setup script as ``{section name: [lines]}``, comments dropped.

    Blank lines and ``;`` comments go, so a section that contains nothing but commentary
    reads as empty — which is what the callers below want: a directive commented out is a
    directive that is not there.
    """
    sections: dict[str, list[str]] = {}
    current: str | None = None
    for raw_line in iss.splitlines():
        line = raw_line.strip()
        if not line or line.startswith(";"):
            continue
        if line.startswith("[") and line.endswith("]"):
            current = line[1:-1]
            sections.setdefault(current, [])
            continue
        if current is not None:
            sections[current].append(line)
    return sections


def _directive(line: str) -> dict[str, str]:
    """One Inno Setup directive line as ``{parameter: value}``, quotes stripped.

    Inno separates a line's parameters with ``;`` and a parameter's name from its value
    with ``:``, so this is the whole grammar. It would misread a value containing either
    character; none in this script does, and :func:`test_the_directives_parse_at_all`
    checks that the shapes this module relies on really came out.
    """
    parts: dict[str, str] = {}
    for chunk in line.split(";"):
        if ":" not in chunk:
            continue
        key, _, value = chunk.partition(":")
        parts[key.strip()] = value.strip().strip('"')
    return parts


def _batch_command_lines(text: str) -> list[str]:
    """The lines of a build script that *do* something — comments and ``echo`` dropped.

    Reads ``.bat`` (``REM``, ``::``) and ``sh`` (``#``) alike, because the pin is asserted
    over both build scripts and the two comment syntaxes are the only difference that
    matters at this resolution.

    The same reasoning as ``tests/test_launcher_dependencies.py``'s ``_command_lines``:
    these scripts describe themselves in comments as well as acting, and a guard that
    reads the description passes while the thing it guards is deleted two lines below.
    ``echo`` goes for the mirror-image reason — a script may legitimately *print* the word
    ``python``.

    The one assertion that reads the raw file instead is
    :meth:`TestTheLaunchScriptUsesIt.test_it_does_not_send_the_user_to_python_org`, and its
    own docstring says why.
    """
    lines = []
    for raw_line in text.splitlines():
        line = raw_line.strip().lstrip("@")
        if not line:
            continue
        lowered = line.lower()
        if lowered.startswith("rem") or line.startswith("::") or line.startswith("#"):
            continue
        if lowered.startswith("echo"):
            continue
        lines.append(line)
    return lines


def _tar_exclusions(text: str) -> set[str]:
    """Every pattern a build script passes to ``tar --exclude``, both syntaxes normalised.

    ``build_dmg.sh`` writes ``--exclude '*/tkinter'`` and ``build_installer.bat`` writes
    ``--exclude=*/tkinter``; this reads either. Taken from the command lines and not the raw
    text, which matters here: the Windows script's *comment* says "the --exclude list is
    build_dmg.sh's", and over raw text that sentence contributes a pattern called ``list``.
    """
    patterns = set()
    for line in _batch_command_lines(text):
        patterns.update(re.findall(r"--exclude[= ]'?([^'\s]+)'?", line))
    return patterns


#: An assignment to one of the pin's two names, in either script's syntax.
#:
#: Searched anywhere in the line rather than anchored to its start, which mutation testing
#: insisted on: ``cmd.exe`` writes its assignments after a condition — ``if /i "%%a"=="…"
#: set PYTHON_VERSION=%%b`` — so an anchored pattern read the batch build's only two
#: assignments as no assignments at all, and passed a script that hardcoded the version.
PIN_ASSIGNMENT = re.compile(r"(?:^|\bset\s+)PYTHON_(VERSION|TAG)=(\S*)")


def _pin_assignments(text: str) -> list[tuple[str, str]]:
    """Every value a build script assigns to ``PYTHON_VERSION``/``PYTHON_TAG`` itself.

    Not every assignment is a second pin: ``build_installer.bat`` legitimately clears both
    names and then fills them from the file it just read. What matters is whether the value
    is a *literal* — see :meth:`TestThePinIsTypedOnce.test_neither_build_pins_an_interpreter_
    of_its_own`, which is the assertion this exists for.
    """
    found = []
    for line in _batch_command_lines(text):
        for match in PIN_ASSIGNMENT.finditer(line):
            found.append((match.group(1), match.group(2).strip()))
    return found


def _pinned_release() -> dict[str, str]:
    """The pinned CPython release, read the way both builds read it.

    ``KEY=value`` and nothing else, because the file has two consumers with very little
    grammar between them: ``bash`` sources it, and ``cmd.exe`` reads it with
    ``for /f "tokens=1,2 delims=="``. :meth:`TestThePinIsTypedOnce.test_the_pin_file_is_
    readable_by_a_shell_and_by_cmd` is what keeps it inside that intersection.
    """
    values: dict[str, str] = {}
    for raw_line in PIN_FILE.read_text().splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        key, _, value = line.partition("=")
        values[key.strip()] = value.strip()
    return values


@needs_build_scripts
class TestTheDirectivesParseAtAll(unittest.TestCase):
    """Before anything else: this module can find the structures it asserts over.

    Every other test here is of the form "X is not in the script". Each one of them is
    satisfied by a parser that read nothing, so the parse gets its own red canary — the
    prior art being the vendor-drift suite, which fails when it cannot parse the DMG
    script's copy list rather than passing over an empty list.
    """

    def test_the_sections_are_found(self):
        sections = _iss_sections(ISS_SCRIPT.read_text())

        self.assertTrue(
            any(line.startswith("AppName=") for line in sections.get("Setup", [])),
            "could not parse installer.iss into sections — [Setup] has no AppName line, so "
            f"the section headers are not where this module looks. Parsed: {sorted(sections)}",
        )
        self.assertTrue(
            sections.get("Files"),
            "installer.iss's [Files] section parsed as empty, so every check below about "
            "what the installer ships would pass over nothing",
        )

    def test_the_directives_parse_at_all(self):
        """A ``Source:``/``DestDir:`` pair really comes out of a [Files] line."""
        sources = [
            _directive(line) for line in _iss_sections(ISS_SCRIPT.read_text())["Files"]
        ]
        app_source = [d for d in sources if d.get("Source", "").endswith("MAFigate.py")]

        self.assertEqual(
            len(app_source),
            1,
            "installer.iss has no single [Files] line whose Source is MAFigate.py, so "
            f"_directive is not reading this script's grammar. Parsed: {sources}",
        )
        self.assertEqual(app_source[0].get("DestDir"), "{app}\\streamlit_app")

    def test_the_launcher_parses_into_commands(self):
        commands = _batch_command_lines(LAUNCH_BAT.read_text())

        self.assertTrue(
            any("MAFigate.py" in line for line in commands),
            "no command line in launch.bat runs MAFigate.py — the script's shape changed "
            f"and _batch_command_lines is no longer reading it. Parsed: {commands}",
        )

    def test_the_build_script_parses_into_commands(self):
        commands = _batch_command_lines(BUILD_BAT.read_text())

        self.assertTrue(
            [line for line in commands if line.startswith("%ISCC%") and "installer.iss" in line],
            "build_installer.bat has no line that compiles installer.iss, so the ordering "
            f"check below has nothing to order against. Parsed: {commands}",
        )


@needs_build_scripts
class TestUninstallKeepsTheUsersWork(unittest.TestCase):
    """``~/.mafigate`` is not the venv. It is also everything the user saved.

    Asserted structurally rather than by naming the two documents: every entry must sit
    *under* the directory and its first component must be one of the machine-generated
    ones. That way a third thing the app decides to keep between sessions survives an
    uninstall by construction, without anyone remembering to add it here.
    """

    #: The only things in ``~/.mafigate`` an uninstall may reclaim: bytes the app generated
    #: and can generate again. Everything else in there was put there by the user.
    MACHINE_GENERATED = {"venv", "logs"}

    def _names(self) -> list[str]:
        entries = _iss_sections(ISS_SCRIPT.read_text()).get("UninstallDelete", [])
        self.assertTrue(
            entries,
            "installer.iss has no [UninstallDelete] section at all. That is not this "
            "ticket's fix — the venv and the logs should still be reclaimed — so it more "
            "likely means this module is no longer reading the file it thinks it is.",
        )
        names = [_directive(line).get("Name", "") for line in entries]
        self.assertTrue(all(names), f"an [UninstallDelete] line has no Name: {entries}")
        return names

    def test_the_users_directory_is_where_the_app_really_keeps_their_work(self):
        """The premise, read off the app rather than transcribed.

        The whole argument below is "this directory holds things the user cannot
        regenerate", so if the app ever moved its cache elsewhere this file's rule would be
        protecting an empty directory while quietly permitting the deletion of the real one.

        Read from the *source* and not imported, which is not squeamishness: a test rig is
        free to redirect where the app keeps things — ``config/file_history.py`` was covered
        here too, and a suite-wide fixture pointed its store at a throwaway home so a test
        run could not touch the developer's own — so an imported constant describes the rig
        rather than the app, and an equality against it would be asserting where *pytest*
        keeps things. That module and its fixture are gone with the recents list; the
        parameter cache is what is left in ``~/.mafigate``, and it is read the same way.
        """
        for module in (STREAMLIT_APP / "page_modules" / "param_store.py",):
            with self.subTest(module=module.name):
                self.assertIn(
                    'Path.home() / ".mafigate"',
                    module.read_text(),
                    f"{module.name} no longer keeps the user's work in ~/.mafigate, so the "
                    "uninstall rule below is guarding the wrong directory",
                )

    def test_no_rule_names_the_parent_directory(self):
        """The bug: ``Type: filesandordirs; Name: "{userprofile}\\.mafigate"``."""
        for name in self._names():
            with self.subTest(name=name):
                self.assertFalse(
                    name.rstrip("\\").endswith(".mafigate"),
                    f"installer.iss deletes {name} on uninstall — that is the whole "
                    "directory, so it takes the parameter cache with it. Name the venv and "
                    "the logs instead.",
                )

    def test_every_rule_names_something_the_app_can_regenerate(self):
        for name in self._names():
            with self.subTest(name=name):
                match = re.match(r"\{userprofile\}\\\.mafigate\\([^\\]+)", name)
                self.assertIsNotNone(
                    match,
                    f"{name} is not a path under the app's own directory; an uninstall "
                    "rule reaching anywhere else needs its own argument",
                )
                self.assertIn(
                    match.group(1),
                    self.MACHINE_GENERATED,
                    f"installer.iss deletes {name} on uninstall, which is not one of the "
                    f"machine-generated {sorted(self.MACHINE_GENERATED)}. Anything else in "
                    "there is the user's own work and survives a reinstall.",
                )

    def test_the_venv_and_the_logs_are_still_reclaimed(self):
        """The other half: keeping the user's work is not an excuse to leave 400 MB."""
        names = " ".join(self._names())

        for reclaimed in sorted(self.MACHINE_GENERATED):
            with self.subTest(reclaimed=reclaimed):
                self.assertIn(f".mafigate\\{reclaimed}", names)


@needs_build_scripts
class TestTheInstallerShipsAnInterpreter(unittest.TestCase):
    """The .exe carries its own CPython, staged by the build script beforehand."""

    def _files(self) -> list[dict[str, str]]:
        return [_directive(line) for line in _iss_sections(ISS_SCRIPT.read_text())["Files"]]

    def test_the_bundled_interpreter_is_installed_into_the_application_directory(self):
        shipped = [
            directive
            for directive in self._files()
            if directive.get("DestDir") == "{app}\\python"
        ]

        self.assertTrue(
            shipped,
            "installer.iss has no [Files] line installing anything into {app}\\python, so "
            "the .exe ships no interpreter and the launcher has nothing to point at",
        )
        for directive in shipped:
            with self.subTest(source=directive.get("Source")):
                self.assertTrue(
                    directive["Source"].startswith("python\\"),
                    f"{directive['Source']} is not the staged interpreter tree that "
                    "build_installer.bat extracts beside installer.iss",
                )
                self.assertIn("recursesubdirs", directive.get("Flags", ""))
                self.assertIn("createallsubdirs", directive.get("Flags", ""))

    def test_the_build_downloads_and_caches_the_pinned_release(self):
        """The DMG's mechanism with a different triple, not a new one.

        The triple is ``x86_64-pc-windows-msvc-install_only_stripped`` and **not** the
        ``…-msvc-shared-install_only_stripped`` issue #261 specified: ``-shared`` is a
        component of this project's ``full`` archives only, every Windows ``install_only``
        build already being the shared one, so the ticket's name is not a published asset of
        any release and would 404 on the build machine. Pinned here so a later reader
        following the ticket back does not "fix" it into a broken download.

        Read off the commands rather than the file, and that is not incidental: the comment
        in ``build_installer.bat`` explaining why ``-shared`` is absent has to *name* the
        wrong triple to warn about it, which fails a raw-text assertion over a correct
        script. Caught here rather than reasoned about, exactly as
        ``test_launcher_dependencies.py`` was.
        """
        build = "\n".join(_batch_command_lines(BUILD_BAT.read_text()))

        self.assertIn("python-build-standalone", build)
        self.assertIn("x86_64-pc-windows-msvc-install_only_stripped", build)
        self.assertNotIn("msvc-shared-install_only", build)
        self.assertIn(".python_cache", build)
        self.assertIn("curl", build.lower())

    def test_the_interpreter_is_staged_before_the_installer_is_compiled(self):
        """Order, because the wrong one ships an .exe with an empty python directory.

        Inno Setup fails a missing ``Source:`` outright, so this is belt to that braces —
        but the failure it prevents is a build script that *looks* right, with the
        download somewhere after the compile step, and the check costs one line.
        """
        commands = _batch_command_lines(BUILD_BAT.read_text())
        extracts = [i for i, line in enumerate(commands) if line.startswith("tar ")]
        compiles = [
            i
            for i, line in enumerate(commands)
            if line.startswith("%ISCC%") and "installer.iss" in line
        ]

        self.assertTrue(extracts, f"build_installer.bat extracts nothing: {commands}")
        self.assertTrue(compiles, f"build_installer.bat compiles nothing: {commands}")
        self.assertLess(extracts[-1], compiles[0])

    def test_it_trims_the_same_components_the_dmg_does(self):
        """"The same way the DMG build does" includes what the DMG throws away.

        Compared against ``build_dmg.sh``'s own exclusion list rather than restated, and as
        an *equality* so it holds in both directions — a one-sided version passed a DMG that
        had stopped trimming ``idlelib``, since dropping a component from the reference set
        only relaxes what the other build must do. Without this the
        Windows tree shipped 73 MB unpacked where the DMG ships 56 — CPython's test suite,
        ``idlelib`` and ``tkinter``, none of it reachable from this app — and nothing would
        have said so: ``test_vendor_drift.py``'s no-tests-in-the-bundle rule reads
        ``Source:`` line *text*, and ``Source: "python\\*"`` names no test path.
        """
        dmg = _tar_exclusions(DMG_SCRIPT.read_text())
        windows = _tar_exclusions(BUILD_BAT.read_text())

        self.assertTrue(
            dmg,
            "could not parse an --exclude list out of build_dmg.sh, so this test would "
            "demand nothing of the Windows build",
        )
        self.assertEqual(
            windows,
            dmg,
            "the two builds no longer trim the same components. Only in "
            f"build_dmg.sh: {sorted(dmg - windows)}; only in build_installer.bat: "
            f"{sorted(windows - dmg)}",
        )

    def test_x64_stays_the_only_windows_build(self):
        """ARM Windows runs the x64 build under emulation; no native variant is built.

        Written down because it is a decision rather than an oversight, and because the
        cheapest way to "support" ARM would be to add a second triple here and ship two
        interpreters that no one has ever launched.
        """
        build = BUILD_BAT.read_text()

        self.assertNotIn("aarch64-pc-windows", build)
        self.assertNotIn("arm64-pc-windows", build)
        self.assertIn(
            "ArchitecturesAllowed=x64compatible",
            ISS_SCRIPT.read_text(),
            "x64compatible is what lets the one build run on ARM Windows under emulation",
        )


@needs_build_scripts
class TestTheLaunchScriptUsesIt(unittest.TestCase):
    """One interpreter, in the application directory, and no search anywhere else."""

    def test_it_points_at_the_bundled_interpreter(self):
        commands = _batch_command_lines(LAUNCH_BAT.read_text())

        bundled = re.compile(r"set PYTHON=%~dp0python\\python\.exe")

        self.assertTrue(
            [line for line in commands if bundled.fullmatch(line)],
            "launch.bat does not set PYTHON to the interpreter installed beside it "
            f"(%~dp0python\\python.exe). Parsed: {commands}",
        )

    def test_it_contains_no_path_probe(self):
        """``where python`` was the first half of the bug: PATH is not ours to trust.

        Read off the commands rather than the file, so the comment explaining that the
        probe is gone does not fail the test that checks it is gone.
        """
        commands = _batch_command_lines(LAUNCH_BAT.read_text())

        for line in commands:
            lowered = line.lower()
            with self.subTest(line=line):
                self.assertNotIn("where python", lowered)
                self.assertNotIn("localappdata%\\programs\\python", lowered)
                self.assertNotRegex(lowered, r"c:\\python\d")

    def test_it_does_not_send_the_user_to_python_org(self):
        """Asserted over the *raw* file, uniquely in this module.

        The instruction was an ``echo``, so reading only the command lines would leave the
        thing this ticket removes invisible to the test that removes it. The cost is that
        ``launch.bat`` may not name the URL in a comment either — deliberate, and the
        lesson ``test_launcher_dependencies.py`` learned the other way round, where a
        comment quoting a removed address failed a correct script.
        """
        self.assertNotIn("python.org", LAUNCH_BAT.read_text().lower())

    def test_it_rebuilds_a_venv_that_cannot_run(self):
        """The upgrade path, which is the one case "points at the bundled interpreter" misses.

        An installation upgraded from a version of MAFigate that used a Python of the user's
        own keeps a venv built by that interpreter, and ``launch.bat`` would reuse it — so
        the app would run on an interpreter this installer never shipped, or fail outright
        once that Python is uninstalled. The repair is gated on the venv's own interpreter
        failing to run, which is why it cannot loop on a healthy one.
        """
        commands = "\n".join(_batch_command_lines(LAUNCH_BAT.read_text()))

        self.assertIn('"%VENV_DIR%\\Scripts\\python.exe" -c "import sys"', commands)
        self.assertIn('rmdir /s /q "%VENV_DIR%"', commands)
        self.assertIn("if errorlevel 1", commands)

    def test_it_still_says_so_when_the_interpreter_is_missing(self):
        """Dropping the error is not the same as dropping the prerequisite.

        A bundled interpreter can still be absent — a half-copied install, an antivirus
        quarantine — and the launcher's answer must be "reinstall MAFigate", not a batch
        script dying on a path that does not exist.
        """
        commands = "\n".join(_batch_command_lines(LAUNCH_BAT.read_text()))

        self.assertIn('if not exist "%PYTHON%"', commands)


@needs_build_scripts
class TestSetupAsksForNothing(unittest.TestCase):
    """No probe, no soft prompt, and a welcome screen that states no prerequisite."""

    def test_there_is_no_python_probe_in_the_setup_hook(self):
        """The soft check is worse than none: it lets a doomed install report success.

        Read off the raw file for the three Pascal tokens, and *not* off the parsed
        ``[Code]`` section, because that section is now absent entirely — so a
        section-scoped assertion would be three ``assertNotIn`` calls over an empty list,
        satisfied by a broken parser as readily as by a fixed installer. A test whose
        subject is an absence cannot depend on locating it.

        The section lookup stays as the second half: if a ``[Code]`` section ever returns
        for some unrelated reason, it may not go looking for a Python.
        """
        iss = ISS_SCRIPT.read_text()
        sections = _iss_sections(iss)
        self.assertIn(
            "Setup",
            sections,
            "installer.iss did not parse into sections at all, so the [Code] check below "
            "is reading nothing",
        )

        lowered = iss.lower()
        self.assertNotIn("initializesetup", lowered)
        self.assertNotIn("mb_yesno", lowered)
        self.assertNotIn("msgbox", lowered)
        self.assertNotIn("python", " ".join(sections.get("Code", [])).lower())

    def test_the_welcome_screen_states_no_python_requirement(self):
        messages = _iss_sections(ISS_SCRIPT.read_text()).get("Messages", [])
        welcome = [line for line in messages if line.startswith("WelcomeLabel2=")]

        self.assertEqual(
            len(welcome),
            1,
            "installer.iss's [Messages] section has no WelcomeLabel2 line, so this test "
            f"is asserting over nothing. Parsed: {messages}",
        )
        self.assertNotIn("python", welcome[0].lower())


@needs_build_scripts
class TestNeitherBatchFileFallsIntoCmdsParsingTraps(unittest.TestCase):
    """Two ways a ``.bat`` reads correctly and behaves otherwise, both met while writing this.

    ``cmd.exe`` expands ``%VAR%`` when it parses a *whole* parenthesised block, before any of
    it runs. Two consequences, and neither shows up until it does:

    * ``if %errorlevel% neq 0`` inside a block tests the value from **before** the block
      started, so the failure branch above it is unreachable. Both of ``launch.bat``'s error
      checks were written that way, which made a failed venv creation and a failed dependency
      install both pass silently.
    * an *unquoted* ``%VAR%`` holding a path with a parenthesis closes the block early. Both
      files interpolate paths the machine chooses — ``%~dp0`` is the install directory in one
      and the checkout in the other — and Windows is a platform where ``C:\\Program Files
      (x86)`` is a perfectly ordinary answer.

    Asserted over the two files this repo has rather than left as knowledge, because the
    idiom that breaks is the idiom everyone writes first.

    The second rule is deliberately blanket: it flags *any* unquoted ``%VAR%`` echoed inside
    a block, without reasoning about what the variable holds. It therefore catches
    ``%APP_NAME%``, which is the fixed string ``MAFigate`` and perfectly safe — and that is
    the right trade, because ``!APP_NAME!`` costs nothing and a rule with a
    which-variables-are-paths exemption list is a rule that will be wrong about the next
    variable someone adds.
    """

    BATCH_FILES = (LAUNCH_BAT, BUILD_BAT)

    def _blocks(self, text: str) -> list[str]:
        """Every line that sits inside an unclosed ``(`` at the time it is parsed.

        Crude but in the safe direction: it tracks nothing but paren depth, with ``^(``/``^)``
        escapes removed first, so a line it reports really is inside a block.
        """
        inside, depth = [], 0
        for raw_line in text.splitlines():
            line = raw_line.strip()
            if line.lower().startswith("rem") or line.startswith("::"):
                continue
            bare = line.replace("^(", "").replace("^)", "")
            if depth > 0:
                inside.append(line)
            depth = max(0, depth + bare.count("(") - bare.count(")"))
        return inside

    def test_the_blocks_are_found_at_all(self):
        """The parse canary: ``launch.bat`` does have blocks, so an empty read is a bug."""
        self.assertTrue(
            self._blocks(LAUNCH_BAT.read_text()),
            "no parenthesised block was found in launch.bat, so both checks below would "
            "pass over nothing",
        )

    def test_no_block_reads_errorlevel_by_expansion(self):
        for path in self.BATCH_FILES:
            for line in self._blocks(path.read_text()):
                with self.subTest(script=path.name, line=line):
                    self.assertNotIn(
                        "%errorlevel%",
                        line.lower(),
                        f"{path.name} tests %errorlevel% inside a block, where cmd has "
                        "already substituted the value from before the block ran. Use "
                        "`if errorlevel 1`.",
                    )

    def test_no_block_echoes_an_unquoted_path_variable(self):
        for path in self.BATCH_FILES:
            for line in self._blocks(path.read_text()):
                if not line.lower().startswith("echo"):
                    continue
                with self.subTest(script=path.name, line=line):
                    self.assertNotRegex(
                        line,
                        r"(?<!\")%[A-Za-z_~][^%\s]*%(?!\")",
                        f"{path.name} echoes an unquoted %VAR% inside a block; a value "
                        "holding a parenthesis closes the block at parse time. Use "
                        "!VAR! where delayed expansion is on, or move the echo out.",
                    )


@needs_build_scripts
class TestThePinIsTypedOnce(unittest.TestCase):
    """One release, one place. Two would be two applications with one version number."""

    def test_the_pin_file_names_a_release(self):
        self.assertTrue(
            PIN_FILE.is_file(),
            f"{PIN_FILE.name} does not exist, so each build names its own CPython and "
            "the two can drift to different interpreters without anything noticing",
        )
        pinned = _pinned_release()

        self.assertRegex(pinned.get("PYTHON_VERSION", ""), r"^\d+\.\d+\.\d+$")
        self.assertRegex(pinned.get("PYTHON_TAG", ""), r"^\d{8}$")

    def test_the_pin_file_is_readable_by_a_shell_and_by_cmd(self):
        """Its grammar is the intersection of two consumers, so it is asserted.

        ``bash`` sources this file and ``cmd.exe`` reads it with ``delims==``. A quoted
        value, a space around the ``=``, or a continuation would break one of the two —
        and it would break it on the build machine, which is the one place nobody is
        running this suite.
        """
        for raw_line in PIN_FILE.read_text().splitlines():
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            with self.subTest(line=line):
                self.assertRegex(line, r"^[A-Z][A-Z0-9_]*=[A-Za-z0-9.+_-]+$")

    def test_neither_build_pins_an_interpreter_of_its_own(self):
        """The assertion with teeth: no *literal* version reaches either build script.

        :meth:`test_neither_build_repeats_the_version` alone was not enough, and mutation
        testing is what said so — replacing the DMG's ``source`` line with
        ``PYTHON_VERSION="3.12.0"`` passed it, because the file was still *mentioned* and
        the string it looked for was the pin the mutant had stopped using. A build that
        names a different CPython is the whole failure being guarded against, so the rule is
        about the shape of the assignment rather than about any one value: empty, or read
        from something (``$``, ``%``), never typed.
        """
        for script in (DMG_SCRIPT, BUILD_BAT):
            for name, value in _pin_assignments(script.read_text()):
                with self.subTest(script=script.name, name=name):
                    self.assertTrue(
                        value == "" or "$" in value or "%" in value,
                        f"{script.name} sets PYTHON_{name} to the literal {value!r} instead "
                        f"of reading {PIN_FILE.name}, so the two builds can name different "
                        "interpreters while both look healthy",
                    )

    #: How each build reads the pin. Two pinned literals, because "reads the file" is the
    #: claim and a script that merely *names* it satisfies nothing: mutation testing passed
    #: a batch build whose loop had been repointed at ``nothing.env``, leaving the pin
    #: mentioned, unread, and both values empty.
    READS_THE_PIN = {
        "build_dmg.sh": 'source "${PYTHON_PIN}"',
        "build_installer.bat": 'in ("%PYTHON_PIN%")',
    }

    def test_both_builds_read_the_pin_file(self):
        for script in (DMG_SCRIPT, BUILD_BAT):
            commands = "\n".join(_batch_command_lines(script.read_text()))
            with self.subTest(script=script.name):
                self.assertIn(
                    PIN_FILE.name,
                    commands,
                    f"{script.name} does not name {PIN_FILE.name} in anything it runs — a "
                    "comment about the pin is not a read of it",
                )
                self.assertIn(self.READS_THE_PIN[script.name], commands)

    def test_neither_build_repeats_the_version(self):
        pinned = _pinned_release()
        for script in (DMG_SCRIPT, BUILD_BAT):
            text = script.read_text()
            with self.subTest(script=script.name):
                self.assertNotIn(
                    pinned["PYTHON_VERSION"],
                    text,
                    f"{script.name} types the pinned version itself; the two builds can "
                    "then be edited apart",
                )
                self.assertNotIn(pinned["PYTHON_TAG"], text)


#: The ``[Setup]`` directives whose value is a path the **compiler** reads. Inno resolves
#: each relative to the directory holding the script, and aborts the compile when one is
#: missing — there is no ``skipifsourcedoesntexist`` for these, which is the whole trouble:
#: the ``[Files]`` entries below express tolerance and these cannot.
#:
#: ``UninstallDisplayIcon`` is deliberately absent. Its value is a ``{app}\…`` path resolved
#: on the *user's* machine at install time, so it names nothing this repository has.
COMPILE_TIME_PATH_DIRECTIVES = (
    "SetupIconFile",
    "LicenseFile",
    "InfoBeforeFile",
    "InfoAfterFile",
    "WizardImageFile",
    "WizardSmallImageFile",
)


def _iss_source_paths(iss: str):
    """Every path ``installer.iss`` asks the compiler to read, as (directive, value) pairs.

    ``[Setup]`` directives and ``[Files]`` ``Source:`` entries alike, since both are read at
    compile time and both abort the build when absent. Lines inside an ``#ifexist`` block are
    skipped: that is exactly how a script says *use this if it is there*, and holding a guarded
    directive to existing would forbid the one construction that makes absence safe.
    """
    found = []
    guarded = False
    section = None
    for raw_line in iss.splitlines():
        line = raw_line.strip()
        if line.startswith("#ifexist"):
            guarded = True
            continue
        if line.startswith("#endif"):
            guarded = False
            continue
        if guarded or not line or line.startswith(";"):
            continue
        if line.startswith("[") and line.endswith("]"):
            section = line[1:-1]
            continue
        if section == "Setup" and "=" in line:
            name, _, value = line.partition("=")
            if name.strip() in COMPILE_TIME_PATH_DIRECTIVES:
                found.append((name.strip(), value.strip().strip('"')))
        elif section == "Files":
            parts = _directive(line)
            source = parts.get("Source")
            if source and "skipifsourcedoesntexist" not in parts.get("Flags", ""):
                found.append(("Source", source))
    return found


def _is_a_build_product(path: Path) -> bool:
    """Does this repository ignore *path* — i.e. is it something a build makes?

    ``Source: "python\\*"`` names the interpreter tree ``build_installer.bat`` extracts before
    it calls the compiler, so it is legitimately absent from a checkout. Asked of git rather
    than listed here, because ``build/.gitignore`` is already the one place every build product
    of these two scripts is declared (issue #260 consolidated it there) — and a hand-written
    exemption list is the thing that goes stale the next time one is added.

    That a staged path is staged *in time* is a separate claim, and
    :meth:`test_the_interpreter_is_staged_before_the_installer_is_compiled` above is the one
    that makes it.

    Asked in both spellings, with and without a trailing slash, and that is not belt and
    braces: ``build/.gitignore`` writes ``windows/python/`` with a slash, which matches
    directories only — and ``git check-ignore`` cannot see that a path it is asked about *is*
    a directory when nothing is there. Without the slash the staged interpreter reads as not
    ignored, which is precisely the case this function exists to recognise.
    """
    for candidate in (f"{path}/", str(path)):
        completed = subprocess.run(
            ["git", "check-ignore", "-q", candidate],
            cwd=STREAMLIT_APP,
            capture_output=True,
        )
        if completed.returncode not in (0, 1):
            raise AssertionError(
                f"git check-ignore could not answer for {candidate}: {completed.stderr!r}. "
                "Failing rather than assuming, since assuming *not ignored* turns a build "
                "product into a false red and assuming *ignored* turns a missing file into a "
                "silent pass."
            )
        if completed.returncode == 0:
            return True
    return False


@needs_build_scripts
class TestEveryPathTheCompilerReadsExists(unittest.TestCase):
    """The class of bug that only a compiler had ever been able to find.

    ``mafigate-v1.0.0-rc2`` — the second rehearsal of the release workflow, and the first
    Windows compile in this project's history to get past parsing — died on
    ``LicenseFile=..\\..\\LICENSE``. Two levels up from ``build/windows/`` is
    ``streamlit_app/``; the licence is at the repository root, one further. The path had been
    wrong since it was written and no test could see it, because every test here read the
    script as *text*.

    ``SetupIconFile=icon.ico`` was the same shape and would have failed on the next line: no
    ``icon.ico`` is tracked, and while the ``[Files]`` and ``[Icons]`` entries for it were
    written with ``skipifsourcedoesntexist``, a ``[Setup]`` directive has no such flag. It is
    now wrapped in ``#ifexist``, which this parser understands as the tolerance it is.

    Resolving a path against the filesystem is cheap, needs no Windows, and would have caught
    both. It does not replace compiling — ``[Icons]``, ``[Run]`` and the wizard are still only
    exercised by a real build — but it retires the two ways this file can be wrong about the
    tree it sits in.
    """

    def setUp(self):
        self.paths = _iss_source_paths(ISS_SCRIPT.read_text(encoding="utf-8"))

    def test_the_parse_found_the_directives_it_claims_to(self):
        """A parser that returned nothing would make every assertion below vacuously true."""
        directives = {name for name, _ in self.paths}
        self.assertIn(
            "LicenseFile",
            directives,
            f"no LicenseFile directive was parsed out of installer.iss; found {directives}",
        )
        self.assertIn(
            "Source",
            directives,
            "no [Files] Source: entries were parsed, so this guard is reading nothing",
        )

    def test_every_unguarded_path_resolves_to_something_in_the_tree(self):
        missing = []
        for name, value in self.paths:
            if "{#" in value or value.startswith("{"):
                # Preprocessor output or an install-time constant; not a path in this tree.
                continue
            # A wildcard names a directory that must exist; the compiler is content with an
            # empty one, so the directory is what is checked.
            target = WINDOWS_DIR / value
            if "*" in value:
                target = WINDOWS_DIR / value.rsplit("\\", 1)[0]
            resolved = Path(str(target).replace("\\", "/"))
            if resolved.exists() or _is_a_build_product(resolved):
                continue
            missing.append(f"{name}={value} → {resolved}")
        self.assertFalse(
            missing,
            "installer.iss asks the compiler to read paths that are not in this tree, and "
            "Inno Setup aborts on each rather than skipping it:\n  " + "\n  ".join(missing),
        )

    def test_a_guarded_directive_is_not_held_to_existing(self):
        """The tolerance has to be expressible, or the fix for a missing icon is a red test."""
        guarded = "#ifexist \"nonesuch.ico\"\nSetupIconFile=nonesuch.ico\n#endif\n"
        self.assertEqual(_iss_source_paths("[Setup]\n" + guarded), [])
        self.assertEqual(
            _iss_source_paths("[Setup]\nSetupIconFile=nonesuch.ico\n"),
            [("SetupIconFile", "nonesuch.ico")],
        )


if __name__ == "__main__":
    unittest.main()
