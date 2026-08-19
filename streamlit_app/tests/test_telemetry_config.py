"""Streamlit's usage reporting is off wherever MAFigate is launched from.

Issue #259. ``browser.gatherUsageStats`` defaults to **on**, so until this repo tracks a
config that turns it off, the strongest promise the app makes — the one the README's
channel table carries, and the Windows installer's welcome screen alongside it — is false.
Shipping the config is one line of TOML. Proving that every launch route *finds* it is the
work, and it is not a formality: Streamlit resolves config from three places, two of which
move with where you launch from, and the four routes here do not agree about that.

**What Streamlit actually reads**, in ascending precedence (``streamlit.config``'s
``get_config_files``):

1. ``~/.streamlit/config.toml`` — the user's, not this repo's to write;
2. ``$PWD/.streamlit/config.toml`` — the *project* config, which follows the working
   directory;
3. ``<directory of the main script>/.streamlit/config.toml`` — the *script-level* config,
   which follows ``MAFigate.py`` instead, and overrides both of the above.

The third is what lets one tracked file cover every route, and it is why this repo holds a
single config rather than one per launcher. **Half** the routes never run from the app
directory at all: ``build/windows/launch.bat`` runs from the installed application root,
one level above ``streamlit_app\\``, and the macOS ``.app`` inherits whatever working
directory Finder hands it — so a project config alone would miss both, and the config would
be "shipped" while being read on neither packaged channel. ``requirements.txt`` pins
``streamlit==1.47.0`` (issue #256), which resolves script-level config, so this is a
guarantee about the software the app runs on rather than a hope — and a pin means the
guarantee changes only when someone deliberately changes it.

What the file's placement therefore has to survive is not the working directory but the
**copy lists**: ``build/mac/build_dmg.sh`` copies a fixed list of names into the bundle and
``installer.iss`` names one ``Source:`` line per directory, so a dotted directory nobody
added by hand is simply absent from both installers. The two simulations below are what
turn that into a test failure rather than a promise that quietly holds only in a checkout.

**Resolution is measured, not modelled.** Each route test derives a working directory and a
main script from the launcher itself, then asks *Streamlit's own config machinery*, in a
subprocess with that working directory, what ``browser.gatherUsageStats`` came out as.
Nothing here re-implements the search order: a guard that modelled it would agree with its
own model rather than with the framework, and a framework behaviour is exactly what is
being relied on.

**The environment is scrubbed on purpose.** All three packaged launchers also export
``STREAMLIT_BROWSER_GATHER_USAGE_STATS=false``, and an environment variable outranks every
file — so a probe that inherited one would report *off* for every route, including a route
whose config had gone missing, which is the same as reporting nothing at all. The probe
therefore drops every ``STREAMLIT_*`` variable and points ``HOME`` at an empty directory,
leaving exactly one thing that can answer: the tracked file, found from where the route
stands. Those exports stay in the launchers as a belt-and-braces fallback; they are simply
not what this file accepts as proof.

**It can go red.** ``test_the_default_really_is_on`` pins the premise the whole module
rests on, and ``test_the_reachability_check_can_identify_an_unreachable_route`` runs the
route assertion against an app directory with no config and fails if it *passes*. Both are
here because a reachability test that cannot name an unreachable route is decoration, and
this suite has shipped that kind before — see ``tests/README.md`` on issue #24.
"""

from __future__ import annotations

import ast
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Callable, NamedTuple

import pytest

# Imported rather than re-derived, and it is the only import here from another test module.
# ``test_delivery_channels_copy.py`` owns this fact — *does anything tracked turn Streamlit's
# usage reporting off?* — and answers it over all three spellings a launch could use: the
# config file, the ``--browser.gatherUsageStats`` flag, or the environment variable. Two
# definitions of one switch is how two guards come to disagree about one tree: were this
# module to key its copy assertions on the config file alone, moving the setting to a
# launcher flag would leave the sibling green while these tests demanded the deletion of
# copy that had just become *more* true. The sibling owns the fact; this module owns the
# surfaces and the routes.
from tests.test_delivery_channels_copy import telemetry_is_turned_off_in_the_tree

STREAMLIT_APP = Path(__file__).resolve().parents[1]
REPO_ROOT = STREAMLIT_APP.parent

#: The one tracked config, beside ``MAFigate.py`` so the script-level lookup finds it.
CONFIG = STREAMLIT_APP / ".streamlit" / "config.toml"

MAIN_SCRIPT_NAME = "MAFigate.py"

#: The option under test, spelled as Streamlit spells it.
OPTION = "browser.gatherUsageStats"

BUILD = STREAMLIT_APP / "build"
DMG_SCRIPT = BUILD / "mac" / "build_dmg.sh"
ISS_SCRIPT = BUILD / "windows" / "installer.iss"
MAC_LAUNCHER = BUILD / "mac" / "MAFigate.app" / "Contents" / "MacOS" / "launch.sh"
WINDOWS_LAUNCHER = BUILD / "windows" / "launch.bat"

#: What an installer must place beside ``MAFigate.py`` for a launch to find the config.
#: The simulations copy exactly these two entries out of each build script's copy list and
#: nothing else — they are what decides config resolution, and dragging ``vendor/`` and
#: ``resources/`` along would cost seconds per test without changing a single answer.
RESOLUTION_ENTRIES = (MAIN_SCRIPT_NAME, ".streamlit")


class Launch(NamedTuple):
    """Where a launch stands: the two paths Streamlit's config search depends on.

    They travel together everywhere below — derived together, measured together, reported
    together — because neither answers the question alone. That is the whole finding of this
    module: the working directory decides the *project* config, the main script decides the
    *script-level* one, and MAFigate's routes disagree about the first.
    """

    working_directory: Path
    main_script: Path


def _text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def _plain(text: str) -> str:
    """Markdown emphasis and code ticks stripped, lowercased, for phrase matching."""
    return re.sub(r"[*`_]", "", text).lower()


def _string_literals(path: Path) -> str:
    """Every string constant in a Python module, joined.

    What a page *says*, as opposed to what its file contains. The difference is the whole
    reason this exists — see
    :func:`test_the_help_page_claims_no_external_query_only_while_that_is_true`.
    """
    literals = [
        node.value
        for node in ast.walk(ast.parse(_text(path)))
        if isinstance(node, ast.Constant) and isinstance(node.value, str)
    ]
    assert literals, f"no string literals parsed out of {path} — this read is not working"
    return "\n".join(literals)


# --------------------------------------------------------------------------------------
# Asking Streamlit
# --------------------------------------------------------------------------------------

#: Three lines of what ``streamlit run`` itself does. The middle one is the load-bearing
#: one and is copied from ``streamlit/web/cli.py``, which sets ``_main_script_path`` before
#: the options are parsed: that assignment is *how* the script-level config is found, so a
#: probe without it would measure only the working directory and pass over the mechanism
#: this module depends on.
_PROBE = "\n".join(
    (
        "import os, sys",
        "import streamlit.config as config",
        "config._main_script_path = os.path.abspath(sys.argv[1])",
        "config.get_config_options(force_reparse=True)",
        f"print(repr(config.get_option({OPTION!r})))",
    )
)


def resolved_usage_stats(launch: Launch) -> bool:
    """What Streamlit makes of :data:`OPTION` for a launch standing at ``launch``.

    A subprocess, because the answer depends on the process's working directory and on
    module state Streamlit computes as it goes. ``HOME`` points at an empty directory and
    every ``STREAMLIT_*`` variable is dropped, so the only thing left that can turn the
    option off is a config file reachable from the working directory or from the main
    script.
    """
    environment = {
        key: value for key, value in os.environ.items() if not key.startswith("STREAMLIT_")
    }
    with tempfile.TemporaryDirectory(prefix="mafigate-empty-home-") as empty_home:
        environment["HOME"] = empty_home
        environment.pop("USERPROFILE", None)
        completed = subprocess.run(
            [sys.executable, "-c", _PROBE, str(launch.main_script)],
            cwd=str(launch.working_directory),
            env=environment,
            capture_output=True,
            text=True,
        )
    if completed.returncode != 0:
        raise AssertionError(
            "could not ask Streamlit what it resolved — the probe itself failed:\n"
            f"{completed.stderr.strip()[-2000:]}"
        )
    answer = completed.stdout.strip().splitlines()[-1] if completed.stdout.strip() else ""
    if answer not in ("True", "False"):
        raise AssertionError(
            f"Streamlit answered {answer!r} rather than a boolean for {OPTION}; the probe "
            "is no longer reading what it thinks it is"
        )
    return answer == "True"


def assert_route_finds_the_setting(label: str, launch: Launch):
    """The one assertion every route makes, so the self-check below can invert it."""
    assert launch.main_script.is_file(), (
        f"{label} launches {launch.main_script}, which does not exist — the route's own "
        "shape could not be read, so nothing below was measured"
    )
    assert resolved_usage_stats(launch) is False, (
        f"launched through {label}, Streamlit resolves {OPTION} to on.\n"
        f"  working directory: {launch.working_directory}\n"
        f"  main script:       {launch.main_script}\n"
        f"Neither {launch.working_directory / '.streamlit' / 'config.toml'} nor "
        f"{launch.main_script.parent / '.streamlit' / 'config.toml'} turned it off, so this "
        "route reports usage to Streamlit and the README's data-posture claim is false "
        "on it."
    )


# --------------------------------------------------------------------------------------
# Where each route stands
# --------------------------------------------------------------------------------------


def _makefile_recipe(makefile: str, target: str) -> str:
    """The indented body of one Make target."""
    lines, body, inside = makefile.splitlines(), [], False
    for line in lines:
        if line.startswith(f"{target}:"):
            inside = True
            continue
        if inside:
            if line.startswith("\t"):
                body.append(line.strip())
            elif line.strip():
                break
    assert body, f"the Makefile has no `{target}` recipe to read"
    return "\n".join(body)


def _launched_script(command_text: str, source: Path) -> str:
    """The argument a ``streamlit run`` line hands over, as written."""
    launched = re.search(r"streamlit run [\"']?([^\s\"']+)", command_text)
    assert launched, (
        f"{source} no longer has a `streamlit run` line this test can read, so the route "
        "it describes is no longer being measured"
    )
    return launched.group(1)


def clone_launcher_route(_tmp_path: Path) -> Launch:
    """``./run_mafigate.sh`` — changes into its own directory before launching."""
    script = STREAMLIT_APP / "run_mafigate.sh"
    text = _text(script)
    assert 'cd "$(dirname "$0")"' in text, (
        "run_mafigate.sh no longer changes into its own directory; this route's working "
        "directory was derived from that line and can no longer be read off it"
    )
    working_directory = script.parent
    return Launch(
        working_directory, (working_directory / _launched_script(text, script)).resolve()
    )


def make_run_route(_tmp_path: Path) -> Launch:
    """``make run`` — Make runs a recipe in the directory holding the Makefile."""
    makefile = STREAMLIT_APP / "Makefile"
    recipe = _makefile_recipe(_text(makefile), "run")
    working_directory = makefile.parent
    return Launch(
        working_directory, (working_directory / _launched_script(recipe, makefile)).resolve()
    )


def windows_installed_route(tmp_path: Path) -> Launch:
    """The Start Menu shortcut — ``launch.bat``, from the installation root.

    The shortcut Inno Setup writes runs ``{app}\\launch.bat``, and the .bat never changes
    directory, so the working directory is the installation root — a level *above* the app.
    That is the route a project-level config would miss, and it is the reason the file sits
    beside ``MAFigate.py`` instead.
    """
    text = _text(WINDOWS_LAUNCHER)
    assert "set APP_DIR=%~dp0streamlit_app" in text, (
        "launch.bat no longer derives its app directory from its own location; this "
        "route's paths were read off that line"
    )
    assert not re.search(r"(?im)^\s*cd[ /]", text), (
        "launch.bat now changes directory — the working directory asserted here is no "
        "longer the one it launches from, so this route needs re-deriving"
    )
    # Inno Setup starts a shortcut in the directory of the file it points at unless an
    # `[Icons]` entry says otherwise, and none does — which is what makes `{app}` the
    # working directory here. An added `WorkingDir:` would move it silently.
    icons = _text(ISS_SCRIPT)
    assert "WorkingDir:" not in icons, (
        "installer.iss now sets a shortcut WorkingDir, so this route's working directory "
        "is whatever that says rather than the installation root asserted here"
    )
    root = simulate_windows_install(tmp_path)
    return Launch(root, root / "streamlit_app" / MAIN_SCRIPT_NAME)


def macos_bundle_route(tmp_path: Path) -> Launch:
    """The ``.app`` — ``launch.sh`` inside the bundle, from Finder's working directory.

    The bundle's executable never changes directory either, and what it inherits from
    Finder is ``/``. Rather than read the launching process's mind, this route stands in an
    empty directory that is *not* the app's: the assertion that matters is that a launch
    from somewhere else still finds the config, and any such directory proves it.
    """
    text = _text(MAC_LAUNCHER)
    for assignment in (
        'RESOURCES_DIR="${BUNDLE_DIR}/Contents/Resources"',
        'APP_DIR="${RESOURCES_DIR}/streamlit_app"',
    ):
        assert assignment in text, (
            f"the macOS launcher no longer spells {assignment!r}, so where it launches "
            "from can no longer be read off it"
        )
    assert not re.search(r"(?m)^\s*cd\s", text), (
        "the macOS launcher now changes directory — this route needs re-deriving"
    )
    app_dir = simulate_macos_bundle(tmp_path)
    elsewhere = tmp_path / "finder-working-directory"
    elsewhere.mkdir()
    return Launch(elsewhere, app_dir / MAIN_SCRIPT_NAME)


# --------------------------------------------------------------------------------------
# Standing the installers up
# --------------------------------------------------------------------------------------


class IssEntry(NamedTuple):
    """One ``[Files]`` directive: what Inno Setup copies, where, and what it leaves out."""

    source: str
    destination: str
    excludes: frozenset
    recursive: bool


def _iss_copy_list(iss: str) -> list[IssEntry]:
    """Every ``Source:`` directive in the Inno Setup script, as written.

    All four fields come off the one line they are written on. An earlier draft read the
    ``Excludes:`` clause back out of the file in a second pass, keyed on the source and
    destination it had already parsed — a second read that could disagree with the first.
    """
    entries = []
    for line in iss.splitlines():
        stripped = line.strip()
        if not stripped.startswith("Source:"):
            continue
        source = re.search(r'Source:\s*"([^"]+)"', stripped)
        destination = re.search(r'DestDir:\s*"([^"]+)"', stripped)
        if not (source and destination):
            continue
        excludes = re.search(r"Excludes:\s*(.+?)(?:;|$)", stripped)
        excluded = (excludes.group(1) if excludes else "").split(",")
        entries.append(
            IssEntry(
                source=source.group(1),
                destination=destination.group(1),
                excludes=frozenset(name.strip().strip('"') for name in excluded),
                recursive="recursesubdirs" in stripped,
            )
        )
    return entries


def simulate_windows_install(tmp_path: Path) -> Path:
    """The tree ``installer.iss`` installs, for the entries config resolution depends on.

    Each entry's own ``DestDir`` is honoured, because *where* the config lands is the whole
    question: a ``.streamlit`` shipped to ``{app}`` rather than ``{app}\\streamlit_app``
    installs a file that no launch would ever read, and a test that assumed the
    destination would pass over exactly that mistake.
    """
    entries = _iss_copy_list(_text(ISS_SCRIPT))
    root = tmp_path / "windows-install"
    copied = []
    for entry in entries:
        relative = entry.source.replace("\\", "/")
        contents_only = relative.endswith("/*")
        origin = ISS_SCRIPT.parent / (relative[:-2] if contents_only else relative)
        if origin.name not in RESOLUTION_ENTRIES or not origin.exists():
            continue
        assert origin.name not in entry.excludes, (
            f"installer.iss excludes {origin.name} from the files it copies"
        )
        under_app = entry.destination.replace("{app}", "").replace("\\", "/")
        target = root.joinpath(*[part for part in under_app.split("/") if part])
        if contents_only:
            # ``dir\*`` copies that directory's files; only ``recursesubdirs`` reaches
            # deeper. Copying the whole tree regardless would make this simulation more
            # generous than the installer it stands in for, and the day the config grew a
            # subdirectory the route test would pass over an .exe that had not shipped it.
            shutil.copytree(
                origin,
                target,
                dirs_exist_ok=True,
                ignore=None if entry.recursive else _ignore_subdirectories,
            )
        else:
            target.mkdir(parents=True, exist_ok=True)
            shutil.copy2(origin, target / origin.name)
        copied.append(origin.name)
    assert MAIN_SCRIPT_NAME in copied, (
        "could not read an app-source copy list out of installer.iss — its Source: lines "
        f"changed shape and this simulation is no longer standing the installer up. "
        f"Parsed: {entries}"
    )
    return root


def _ignore_subdirectories(directory, names):
    """A :func:`shutil.copytree` filter that drops directories, keeping files."""
    return [name for name in names if (Path(directory) / name).is_dir()]


def _dmg_copy_list(dmg: str) -> tuple[list[str], str]:
    """The literal names ``build_dmg.sh`` copies, and the bundle path it copies them into.

    Shell expansions — the loop body's own ``cp -R "${PROJECT_DIR}/${item}"`` — are
    dropped rather than returned as names, so a caller's "did this parse anything?" check
    stays falsifiable.

    This script is read by two other places for their own purposes: ``test_vendor_drift.py``
    wants the names alone, and ``tools/export_public.py`` scans the copy list for offenders.
    Each keeps its own reading and its own anchor assertion, and this one needs a shape
    neither has — the destination as well as the names, since where the config lands is the
    question. A shared parser returning the union of three shapes would be the abstraction
    none of the three asked for.
    """
    names = []
    for line in dmg.splitlines():
        stripped = line.strip()
        if stripped.startswith("for item in ") and "; do" in stripped:
            names.extend(stripped[len("for item in ") : stripped.index("; do")].split())
        match = re.search(r'cp -R\s+"\$\{PROJECT_DIR\}/([^"$]+)"', stripped)
        if match:
            names.append(match.group(1))
    destination = re.search(r'DEST="\$\{APP_BUNDLE\}/([^"]+)"', dmg)
    assert destination, (
        "could not read where build_dmg.sh copies the app source; the DEST assignment "
        "changed shape and this simulation is no longer standing the bundle up"
    )
    return names, destination.group(1)


def simulate_macos_bundle(tmp_path: Path) -> Path:
    """The ``streamlit_app`` directory ``build_dmg.sh`` assembles inside the ``.app``."""
    dmg = _text(DMG_SCRIPT)
    names, destination = _dmg_copy_list(dmg)
    assert MAIN_SCRIPT_NAME in names, (
        "could not read the app-source copy list out of build_dmg.sh — the script's shape "
        f"changed and this simulation is no longer reading it. Parsed: {names}"
    )
    app_dir = tmp_path / "MAFigate.app" / destination
    app_dir.mkdir(parents=True)
    for name in names:
        origin = STREAMLIT_APP / name
        if name not in RESOLUTION_ENTRIES or not origin.exists():
            continue
        _assert_the_cleanup_spares(dmg, name)
        if origin.is_dir():
            shutil.copytree(origin, app_dir / name)
        else:
            shutil.copy2(origin, app_dir / name)
    return app_dir


def _assert_the_cleanup_spares(dmg: str, name: str):
    """``build_dmg.sh`` prunes the copied tree afterwards; it must not prune ``name``.

    The script copies and then deletes — ``tests``, ``build``, ``__pycache__``,
    ``.DS_Store`` — so a copy list naming something the cleanup removes ships nothing,
    and would do it silently. Cheap to check, and the alternative is a simulation that
    models only half of what the build does.
    """
    for line in dmg.splitlines():
        stripped = line.strip()
        if "${DEST}" not in stripped:
            continue
        if not (stripped.startswith("rm ") or "-delete" in stripped or "-exec rm" in stripped):
            continue
        assert name not in stripped, (
            f"build_dmg.sh copies {name} into the bundle and then removes it again: "
            f"{stripped!r}"
        )


# --------------------------------------------------------------------------------------
# The routes
# --------------------------------------------------------------------------------------


class Route(NamedTuple):
    """One way MAFigate gets started, and where that leaves Streamlit standing."""

    label: str
    #: The tracked file that launches Streamlit, repo-relative. Read by the completeness
    #: test below, which is what stops a sixth route arriving unmeasured.
    launcher: str
    derive: Callable[[Path], Launch]


ROUTES = (
    Route("./run_mafigate.sh", "streamlit_app/run_mafigate.sh", clone_launcher_route),
    Route("make run", "streamlit_app/Makefile", make_run_route),
    Route(
        "the Windows installer's launch.bat",
        "streamlit_app/build/windows/launch.bat",
        windows_installed_route,
    ),
    Route(
        "the macOS .app launcher",
        "streamlit_app/build/mac/MAFigate.app/Contents/MacOS/launch.sh",
        macos_bundle_route,
    ),
)


def _tracked_files() -> list[str]:
    """Every file git tracks, as repo-relative paths.

    Through git rather than a filesystem walk, for the reason that has bitten this repo
    before: ``.claude/worktrees`` holds whole checkouts of the same tree, so an ``rglob``
    from the repo root answers questions about a branch this test is not running on.
    """
    completed = subprocess.run(
        ["git", "ls-files"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=True,
    )
    tracked = [line for line in completed.stdout.splitlines() if line]
    assert tracked, "git ls-files listed nothing — this guard would pass over any tree"
    return tracked


def config_is_tracked() -> bool:
    """Is the telemetry config a file the repo ships, rather than one on this machine?"""
    return CONFIG.relative_to(REPO_ROOT).as_posix() in _tracked_files()


# --------------------------------------------------------------------------------------
# The tests
# --------------------------------------------------------------------------------------


def test_the_telemetry_config_is_tracked():
    """An untracked config protects the machine it was written on and nobody else."""
    assert config_is_tracked(), (
        f"{CONFIG.relative_to(REPO_ROOT)} is not tracked. Streamlit reports usage by "
        "default, so every claim in this repo that nothing leaves the user's machine "
        "rests on this file being one the clone and both installers carry"
    )


def test_the_config_does_not_publish_the_server():
    """A config file is a fine place to put a setting; it is a terrible place to lose one.

    ``run_mafigate.sh`` read ``--server.address 0.0.0.0`` until issue #182, which published
    MAFigate on every interface and let anyone who could route to the machine open patient
    data through it. That flag is now loopback on every launcher — and this file is the
    first place since where a server setting could arrive without passing through one of
    them, so it may not be where the address comes back.

    Read raw rather than through :func:`_plain`, which strips ``*`` for markdown's sake and
    so deleted the wildcard needle before it was ever searched for — the first draft of this
    test could not fail on ``address = "*"``, and did not, when that line was appended to
    the config and the assertion still passed.
    """
    text = _text(CONFIG).lower()
    for published in ("0.0.0.0", '"*"', "'*'"):
        assert published not in text, (
            f"{CONFIG.relative_to(REPO_ROOT)} sets a server address of {published}, which "
            "publishes MAFigate beyond the machine it runs on (issue #182)"
        )


def test_the_config_does_not_move_the_upload_cap():
    """The Help page states 200 MB as a fact, and this file is a new way to falsify it.

    ``test_help_claims.py`` holds the page to Streamlit's real cap, and its reasoning is
    that nothing in the repo raises ``server.maxUploadSize``: the two places it could have
    been raised were ``run_mafigate.sh`` and the ``Makefile``, and
    ``test_launcher_dependencies.py`` watches both. Adding a config file adds a third, one
    neither of those tests reads — so raising the cap here would leave the Help page
    telling a clinician a number the server no longer enforces, with nothing red.
    """
    assert "maxuploadsize" not in _plain(_text(CONFIG)), (
        f"{CONFIG.relative_to(REPO_ROOT)} sets server.maxUploadSize. The Help page states "
        "Streamlit's 200 MB default as the cap it enforces, so moving it here makes that "
        "sentence false — raise it on the command line, or correct the page in the same "
        "change"
    )


@pytest.mark.parametrize("route", ROUTES, ids=lambda route: route.label)
def test_every_launch_route_resolves_usage_reporting_off(route: Route, tmp_path: Path):
    """The claim, per route, answered by Streamlit rather than by this file."""
    assert_route_finds_the_setting(route.label, route.derive(tmp_path))


def test_the_default_really_is_on(tmp_path: Path):
    """The premise every assertion above rests on, measured rather than cited.

    If Streamlit's default were already off, each route test would pass on a tree with no
    config in it at all and this module would be measuring nothing.
    """
    app_dir = tmp_path / "streamlit_app"
    app_dir.mkdir()
    (app_dir / MAIN_SCRIPT_NAME).write_text("")
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    assert resolved_usage_stats(Launch(elsewhere, app_dir / MAIN_SCRIPT_NAME)) is True, (
        "Streamlit no longer reports usage by default, which is the thing this whole "
        "module exists to turn off — re-read it before deleting anything"
    )


def test_the_reachability_check_can_identify_an_unreachable_route(tmp_path: Path):
    """The route assertion, run against a route that cannot find the config, must fail.

    Constructed rather than argued. The app directory here holds ``MAFigate.py`` and no
    ``.streamlit``, and the working directory is a sibling — which is precisely the shape
    of a bundle whose copy list forgot the config, the failure the two simulations exist
    to catch. If this test ever passes, every route test above is vacuous.
    """
    app_dir = tmp_path / "bundle" / "streamlit_app"
    app_dir.mkdir(parents=True)
    (app_dir / MAIN_SCRIPT_NAME).write_text("")
    elsewhere = tmp_path / "bundle"
    with pytest.raises(AssertionError, match="resolves browser.gatherUsageStats to on"):
        assert_route_finds_the_setting(
            "a bundle whose copy list forgot the config",
            Launch(elsewhere, app_dir / MAIN_SCRIPT_NAME),
        )


#: A Streamlit executable, invoked with ``run``. Loose about what stands between the two
#: because the four launchers spell it three ways — ``-m streamlit run``,
#: ``"${VENV_STREAMLIT}" run``, ``"%VENV_STREAMLIT%" run`` — and looser still about
#: newlines, because a launcher may build its command as a list with the executable and
#: ``"run"`` on separate lines. ``build/launcher.py`` did exactly that until issue #257
#: deleted it, which is why the match is made against the file rather than line by line: the
#: next such launcher should be caught by this sweep, not missed by it.
LAUNCH_INVOCATION = re.compile(r"""streamlit[\w%${}"'\\/.()-]*["']?[\s,]*["']?run\b""", re.I)

#: Lines that talk about the command instead of running it. Three tracked files do:
#: ``setup.sh`` closes by printing it as advice, ``check_dependencies.py`` mentions
#: ``make run`` in a comment, and ``pyproject.toml`` notes how the app starts. Dropping
#: them is what lets the match above be as loose as it is.
COMMENTARY_PREFIXES = ("#", "//", ";", "REM ", "rem ", "echo ", "@echo", "*", '"""', "'''")

#: Where a launcher could live: something executable, inside the app, outside the suite.
LAUNCHER_SUFFIXES = (".sh", ".bat", ".py", "")


def _commands_only(text: str) -> str:
    """``text`` with its commentary lines removed."""
    return "\n".join(
        line for line in text.splitlines() if not line.strip().startswith(COMMENTARY_PREFIXES)
    )


def test_no_launcher_in_the_tree_is_missing_from_the_routes():
    """A fifth way to start the app must not arrive unmeasured.

    The route list above is written by hand, because there is no reading a working
    directory out of an arbitrary shell script. What keeps it honest is this: every tracked
    file under ``streamlit_app/`` that invokes Streamlit with ``run`` has to be one of the
    four. A new launcher then arrives with a red test naming it rather than with a claim in
    the README quietly covering one route fewer than it says — and a *deleted* one arrives
    the same way, which is how this sweep reported ``build/launcher.py``'s removal in #257.
    """
    launching = set()
    for relative in _tracked_files():
        if not relative.startswith("streamlit_app/"):
            continue
        if relative.startswith("streamlit_app/tests/") or relative.endswith(".md"):
            continue
        if Path(relative).suffix not in LAUNCHER_SUFFIXES:
            continue
        path = REPO_ROOT / relative
        if not path.is_file():
            continue
        if LAUNCH_INVOCATION.search(_commands_only(_text(path))):
            launching.add(relative)
    assert launching == {route.launcher for route in ROUTES}, (
        "the set of tracked files that launch Streamlit no longer matches the routes this "
        "module measures.\n"
        f"  launching, unmeasured: {sorted(launching - {r.launcher for r in ROUTES})}\n"
        f"  measured, no longer launching: {sorted({r.launcher for r in ROUTES} - launching)}"
    )


# --------------------------------------------------------------------------------------
# The three surfaces next door
# --------------------------------------------------------------------------------------

HELP_PAGE = STREAMLIT_APP / "page_modules" / "help.py"
PARAMETER_CONFIG = STREAMLIT_APP / "page_modules" / "parameter_config.py"

#: The Windows installer's welcome screen, before a single file has been copied.
NO_UPLOAD_CLAIM = "nothing is uploaded to external servers"

#: The Help page's answer to "where does my data go?".
NO_EXTERNAL_QUERY_CLAIM = "queries no external service"

#: The parameter cache's *Technical Details* line.
CACHE_LOCALITY_CLAIM = "cache is stored locally on your computer only"

#: The README's canonical data-posture sentence (issue #229), pinned in
#: ``test_delivery_channels_copy.py``. Named here only to keep it *out* of these three
#: surfaces: one statement, in one table, and every other page pointing at it.
CANONICAL_CLAIM = "leaves your machine"


def test_the_installer_promises_no_upload_only_while_nothing_is_uploaded():
    """The welcome screen makes this promise before the user has agreed to anything.

    ``WelcomeLabel2`` reads *"All data stays on your computer - nothing is uploaded to
    external servers"*, and while Streamlit reported usage that was false — the process
    the installer is about to put on the machine phoned home on every launch. It is true
    now, and it is true *because* something in the tree turns that reporting off, which is
    what this asserts: the sentence and that fact stand or fall together.

    The wording is pinned, at the usual cost — a maintainer who rewords it gets a red
    test. Issue #229 pinned one sentence for the README for the reason that applies here
    too: this is the same promise on a second surface, and two surfaces wording the same
    promise separately is how one of them comes to outlive the thing that made it true.
    """
    welcome = _plain(_text(ISS_SCRIPT))
    if telemetry_is_turned_off_in_the_tree():
        assert NO_UPLOAD_CLAIM in welcome, (
            f"the installer's welcome screen must state {NO_UPLOAD_CLAIM!r} — the config "
            "makes it true, and this is the surface a Windows user reads first"
        )
    else:
        assert NO_UPLOAD_CLAIM not in welcome, (
            f"the installer's welcome screen promises {NO_UPLOAD_CLAIM!r} while nothing "
            "tracked turns Streamlit's usage reporting off"
        )


def test_the_help_page_claims_no_external_query_only_while_that_is_true():
    """*"queries no external service"* — the reader's question is about the process.

    Whose telemetry it is does not reach the person opening a patient MAF: what they are
    told is that the software they just launched makes no outbound call, and until this
    config landed, it made one on every launch. The same two-state switch as the
    installer's, on the sentence the Help page's *Technical Questions* answer opens with.

    **Read off the page's string literals, not its file.** The first draft read the whole
    source, and this change is what made that vacuous: the commit adding it also added a
    code comment quoting the sentence, so the guard began certifying its own commentary. It
    was caught by rewriting the *rendered* answer to say the app "queries three external
    services" and watching the test stay green. Literals through :mod:`ast` cannot be
    satisfied by a comment.

    Not a render, unlike the rest of this page's claims (``test_help_claims.py``, which owns
    them and holds the rendering fixture). A literal is the weaker instrument for *"the page
    must say X"* — a string in a branch nothing draws would count — and the right one for
    *"the page must not say X"*, which is the direction that matters while the claim is
    false.
    """
    page = _plain(_string_literals(HELP_PAGE))
    if telemetry_is_turned_off_in_the_tree():
        assert NO_EXTERNAL_QUERY_CLAIM in page, (
            f"the Help page must keep {NO_EXTERNAL_QUERY_CLAIM!r} — with usage reporting "
            "off it is now true without qualification, and it is the answer to the one "
            "question a clinician asks about patient data"
        )
    else:
        assert NO_EXTERNAL_QUERY_CLAIM not in page, (
            f"the Help page tells the user MAFigate {NO_EXTERNAL_QUERY_CLAIM!r} while "
            "nothing tracked stops Streamlit reporting usage from the same process"
        )


def test_the_parameter_cache_line_stays_about_the_cache():
    """The third surface #259 asked about, and the one whose truth never depended on this.

    *"Cache is stored locally on your computer only"* was the app's last false locality
    claim while a hosted channel existed — on a shared server it named someone else's
    disk (issue #224). The channel is gone (#229), so the sentence is true and stays as
    written; no two-state switch here, because Streamlit's usage reporting never carried
    a parameter cache anywhere.

    What is worth holding is the boundary. This line is about a JSON file in the user's
    home directory, not about where the MAF goes — that promise has exactly one home, the
    README's channel table, and a copy of it here would be the third wording of one claim.
    """
    text = _plain(_text(PARAMETER_CONFIG))
    assert CACHE_LOCALITY_CLAIM in text, (
        f"the parameter cache's privacy line must state {CACHE_LOCALITY_CLAIM!r} — it is "
        "the only place the app says where the cache it writes without asking ends up"
    )
    assert CANONICAL_CLAIM not in text, (
        "the parameter cache's technical details restate the data-posture promise, which "
        "belongs in the README's channel table and nowhere else (issue #229)"
    )
