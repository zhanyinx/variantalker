"""The mac launcher's venv survives the app having moved (issues #411, #412, #414).

The field failure this guards against: v1.0.0's first launch worked, and the second
died with ``Error: [Errno 2] No such file or directory: '~/.mafigate/venv/bin/python3'``.
The venv was built with default symlinks pinning the bundled interpreter at its
ephemeral first-launch path (Gatekeeper translocation gives the app a new path per
launch); on relaunch the ``[ ! -x venv/bin/python ]`` branch set ``NEED_NEW_VENV=true``
without clearing the stale venv — its two sibling branches do clear — and
``python3 -m venv`` over the stale directory skips existing symlinks, even dangling
ones, so ensurepip then executed the dead ``bin/python3``. Reproduced headlessly on
issue #412 with exactly the harness shape used here.

The functional tests load ``launch.sh`` through its ``MAFIGATE_LAUNCH_LIB`` seam and
drive ``ensure_venv()`` the way the real launcher does — dialogs stubbed, ``HOME``
sandboxed — against a fake bundle whose interpreter is a movable copy of this test
run's own base (``sys.base_prefix``). Moving that copy between two runs dangles the
venv's symlinks exactly like the field failure; no network, no fetch (the house rule:
tests make no outbound calls). Fidelity to the exact shipped python-build-standalone
interpreter was proven once, by hand, on issue #412 — the mechanism is python-agnostic.

If this interpreter's base cannot be copied into a working python, every case would
fail for a reason that is not the launcher's, so the staging step skips loudly with
the reason instead (the pattern of ``test_launcher_dependencies``'s one real-venv test).
"""

import atexit
import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pytest

if sys.platform == "win32":
    pytest.skip(
        "launch.sh is the mac launcher; the Windows launcher is immune by construction "
        "(issue #413) and has its own guards in test_windows_installer.py",
        allow_module_level=True,
    )

TESTS_DIR = Path(__file__).resolve().parent
APP_DIR = TESTS_DIR.parent
LAUNCH_SH = APP_DIR / "build" / "mac" / "MAFigate.app" / "Contents" / "MacOS" / "launch.sh"

#: The stamp exactly as the v1.0.0 launcher wrote it: version and architecture only,
#: no base-interpreter path — which is why a same-build move could never trip the
#: stamp branch (issue #412's finding, issue #411's fix).
V100_STAMP_EXPR = (
    "import sys, platform; "
    'print(f"{sys.version_info.major}.{sys.version_info.minor}.'
    '{sys.version_info.micro}-{platform.machine()}")'
)

SENTINEL = "mafigate-test-sentinel"

#: The driver mirrors the real launcher's own flow — resolve the interpreter with
#: find_python, then ensure_venv — with the osascript dialog helpers stubbed to
#: prints, which is the same stubbing the issue #412 reproduction used.
DRIVER = """#!/bin/bash
set -e
export MAFIGATE_LAUNCH_LIB=1
source "$(cd "$(dirname "$0")" && pwd)/launch.sh"
show_error() { echo "STUB_ERROR: $1"; exit 1; }
show_notification() { echo "STUB_NOTIFICATION: $1"; }
show_progress() { true; }
PYTHON=$(find_python) || show_error "find_python found no interpreter in the bundle"
ensure_venv
echo "ENSURE_VENV_DONE"
"""

_STAGED = {}


def _pristine_base() -> Path:
    """A movable copy of this interpreter's base, staged once per test run.

    Copied with ``symlinks=True`` so venvs built from the copy symlink into it —
    moving the copy then dangles them, which is the field mechanism under test.
    """
    if "outcome" in _STAGED:
        kind, value = _STAGED["outcome"]
        if kind == "skip":
            raise unittest.SkipTest(value)
        return value
    root = Path(tempfile.mkdtemp(prefix="mafigate-launcher-base-"))
    atexit.register(shutil.rmtree, root, ignore_errors=True)
    base = root / "python"
    try:
        shutil.copytree(sys.base_prefix, base, symlinks=True)
        python3 = base / "bin" / "python3"
        if not python3.exists():
            versioned = sorted(
                candidate
                for candidate in (base / "bin").glob("python3.*")
                if re.fullmatch(r"python3\.\d+", candidate.name)
            )
            if not versioned:
                raise RuntimeError(f"no python3 binary under {base / 'bin'}")
            python3.symlink_to(versioned[0].name)
        probe = root / "probe-venv"
        completed = subprocess.run(
            [str(python3), "-m", "venv", str(probe)],
            capture_output=True,
            text=True,
            timeout=600,
        )
        if completed.returncode != 0:
            raise RuntimeError(
                f"the copied base cannot build a venv: {completed.stderr.strip()!r}"
            )
        # A framework build's copy resolves back to the original installation
        # (its bin/python3 is a stub that execs a compiled-in absolute path), so
        # venvs built from such a copy anchor to the original, not the copy.
        # Record it: one anchor assertion below is unaskable on such a base.
        probe_home = ""
        for line in (probe / "pyvenv.cfg").read_text().splitlines():
            key, _, value = line.partition("=")
            if key.strip() == "home":
                probe_home = os.path.realpath(value.strip())
        _STAGED["copy_anchors_to_itself"] = probe_home.startswith(
            os.path.realpath(str(base))
        )
    except (OSError, RuntimeError, subprocess.SubprocessError) as exc:
        reason = (
            f"cannot stage a movable copy of this interpreter ({sys.base_prefix}): {exc}"
        )
        _STAGED["outcome"] = ("skip", reason)
        raise unittest.SkipTest(reason)
    _STAGED["outcome"] = ("base", base)
    return base


def _clean_env(home: Path) -> dict:
    """The launcher's environment: a sandboxed HOME, no interpreter overrides."""
    env = {
        key: value
        for key, value in os.environ.items()
        if not key.startswith(("PYTHON", "VIRTUAL_ENV", "MAFIGATE", "__PYVENV"))
    }
    env["HOME"] = str(home)
    return env


def _build_bundle(install_dir: Path) -> Path:
    """A fake MAFigate.app: the real launch.sh beside a copied interpreter."""
    bundle = install_dir / "MAFigate.app"
    macos_dir = bundle / "Contents" / "MacOS"
    resources_dir = bundle / "Contents" / "Resources"
    macos_dir.mkdir(parents=True)
    resources_dir.mkdir(parents=True)
    shutil.copytree(_pristine_base(), resources_dir / "python", symlinks=True)
    shutil.copy2(LAUNCH_SH, macos_dir / "launch.sh")
    driver = macos_dir / "run_ensure_venv.sh"
    driver.write_text(DRIVER)
    driver.chmod(0o755)
    return bundle


def _run_launcher(bundle: Path, home: Path) -> subprocess.CompletedProcess:
    return subprocess.run(
        ["bash", str(bundle / "Contents" / "MacOS" / "run_ensure_venv.sh")],
        capture_output=True,
        text=True,
        env=_clean_env(home),
        timeout=600,
    )


def _transcript(completed: subprocess.CompletedProcess) -> str:
    return f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"


class _LauncherCase(unittest.TestCase):
    def setUp(self):
        root = Path(tempfile.mkdtemp(prefix="mafigate-launcher-case-"))
        self.addCleanup(shutil.rmtree, root, True)
        self.root = root
        self.home = root / "home"
        self.home.mkdir()

    @property
    def venv_dir(self) -> Path:
        return self.home / ".mafigate" / "venv"

    def venv_python_runs(self) -> subprocess.CompletedProcess:
        return subprocess.run(
            [str(self.venv_dir / "bin" / "python3"), "-c", "import sys"],
            capture_output=True,
            text=True,
            env=_clean_env(self.home),
        )

    def cfg_home(self) -> str:
        for line in (self.venv_dir / "pyvenv.cfg").read_text().splitlines():
            key, _, value = line.partition("=")
            if key.strip() == "home":
                return value.strip()
        raise AssertionError(
            f"{self.venv_dir / 'pyvenv.cfg'} has no 'home' key:\n"
            f"{(self.venv_dir / 'pyvenv.cfg').read_text()}"
        )

    def seed_v100_venv(self, bundle: Path) -> None:
        """The venv exactly as the shipped v1.0.0 launcher left it: built straight
        from the bundle's interpreter, stamped version-arch only, deps marked done.

        The interpreter symlinks are re-pointed at the bundle's own path afterwards,
        because that is the state the field machines are in — and because a venv
        built from a *copied* framework python resolves its base back to the
        original installation (measured here with the Homebrew build), which would
        quietly turn the whole scenario move-immune and this guard vacuous.
        """
        bundle_python = bundle / "Contents" / "Resources" / "python" / "bin" / "python3"
        self.venv_dir.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run(
            [str(bundle_python), "-m", "venv", str(self.venv_dir)],
            check=True,
            capture_output=True,
            text=True,
            env=_clean_env(self.home),
            timeout=600,
        )
        for entry in (self.venv_dir / "bin").iterdir():
            if entry.is_symlink() and os.path.isabs(os.readlink(entry)):
                entry.unlink()
                entry.symlink_to(bundle_python)
        stamp = subprocess.run(
            [str(bundle_python), "-c", V100_STAMP_EXPR],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
        (self.venv_dir / ".python_stamp").write_text(stamp + "\n")
        (self.venv_dir / ".deps_installed").touch()
        (self.venv_dir / SENTINEL).touch()


class TestARelaunchSurvivesTheAppMoving(_LauncherCase):
    """The field defect: first launch at one path, second launch at another."""

    @pytest.mark.slow
    def test_the_second_launch_from_a_new_path_still_yields_a_working_venv(self):
        site_a = self.root / "site-a"
        site_a.mkdir()
        bundle = _build_bundle(site_a)

        first = _run_launcher(bundle, self.home)
        self.assertEqual(
            first.returncode,
            0,
            f"the FIRST launch failed, which is not the field defect:\n{_transcript(first)}",
        )
        (self.venv_dir / SENTINEL).touch()

        site_b = self.root / "site-b"
        shutil.move(str(site_a), str(site_b))
        moved_bundle = site_b / "MAFigate.app"

        second = _run_launcher(moved_bundle, self.home)
        self.assertEqual(
            second.returncode,
            0,
            "the second launch after the app moved died — the v1.0.0 field defect "
            f"(issue #412):\n{_transcript(second)}",
        )
        alive = self.venv_python_runs()
        self.assertEqual(
            alive.returncode,
            0,
            f"the venv the second launch left behind cannot run:\n{alive.stderr}",
        )
        dangling = [
            str(path)
            for path in (self.venv_dir / "bin").iterdir()
            if path.is_symlink() and not path.exists()
        ]
        self.assertEqual(
            dangling,
            [],
            "the second launch left dangling symlinks in the venv, so it wrote over "
            "the stale environment instead of clearing it",
        )

    @pytest.mark.slow
    def test_a_broken_v100_venv_heals_with_no_manual_cleanup(self):
        """The upgrade path: a user whose v1.0.0 venv already dangles relaunches."""
        site_a = self.root / "site-a"
        site_a.mkdir()
        bundle = _build_bundle(site_a)
        self.seed_v100_venv(bundle)

        site_b = self.root / "site-b"
        shutil.move(str(site_a), str(site_b))
        moved_bundle = site_b / "MAFigate.app"

        healed = _run_launcher(moved_bundle, self.home)
        self.assertEqual(
            healed.returncode,
            0,
            "a launch over a dangling v1.0.0 venv still dies instead of clearing and "
            f"rebuilding it:\n{_transcript(healed)}",
        )
        self.assertFalse(
            (self.venv_dir / SENTINEL).exists(),
            "the stale venv's contents survived the rebuild — cleared, not overwritten, "
            "is the contract (issue #414)",
        )
        alive = self.venv_python_runs()
        self.assertEqual(alive.returncode, 0, f"the healed venv cannot run:\n{alive.stderr}")


class TestAWorkingV100VenvMigrates(_LauncherCase):
    """A still-working v1.0.0 venv is still move-fragile, so it must rebuild too."""

    @pytest.mark.slow
    def test_an_old_format_stamp_forces_a_whole_venv_rebuild(self):
        site = self.root / "site"
        site.mkdir()
        bundle = _build_bundle(site)
        self.seed_v100_venv(bundle)

        migrated = _run_launcher(bundle, self.home)
        self.assertEqual(migrated.returncode, 0, _transcript(migrated))
        self.assertFalse(
            (self.venv_dir / SENTINEL).exists(),
            "a venv stamped in the v1.0.0 format (version-arch, no base path) was kept, "
            "so existing installs would stay pinned to a path the next app move kills "
            "(issue #411's migration promise)",
        )
        self.assertFalse(
            (self.venv_dir / ".deps_installed").exists(),
            "the rebuild must delete ~/.mafigate/venv whole — .deps_installed lives "
            "inside it and a survivor would skip the dependency install into the "
            "fresh venv",
        )
        alive = self.venv_python_runs()
        self.assertEqual(alive.returncode, 0, f"the rebuilt venv cannot run:\n{alive.stderr}")


class TestTheVenvIsAnchoredUnderMafigateHome(_LauncherCase):
    """The root fact of move-immunity: the venv's base lives under ~/.mafigate."""

    def _first_launch_cfg_home(self) -> tuple:
        site = self.root / "site"
        site.mkdir()
        bundle = _build_bundle(site)
        first = _run_launcher(bundle, self.home)
        self.assertEqual(first.returncode, 0, _transcript(first))
        return bundle, Path(os.path.realpath(self.cfg_home()))

    @pytest.mark.slow
    def test_pyvenv_cfg_home_never_points_into_the_app_bundle(self):
        bundle, home_value = self._first_launch_cfg_home()
        self.assertFalse(
            home_value.is_relative_to(Path(os.path.realpath(bundle))),
            f"pyvenv.cfg pins home = {home_value}, inside the app bundle — the path "
            "Gatekeeper translocation moves between launches, which is the v1.0.0 "
            "field defect (issue #412)",
        )

    @pytest.mark.slow
    def test_pyvenv_cfg_home_points_under_the_mafigate_directory(self):
        _pristine_base()
        if not _STAGED.get("copy_anchors_to_itself"):
            self.skipTest(
                "a venv built from a copy of this interpreter anchors back to the "
                "original installation (a framework-build property), so where the "
                "launcher anchors the venv cannot be asked here; "
                "test_pyvenv_cfg_home_never_points_into_the_app_bundle and the "
                "re-pointed-symlink cases still hold the contract on this platform, "
                "and this assertion runs where copies self-anchor (the ubuntu leg)"
            )
        _, home_value = self._first_launch_cfg_home()
        anchor = Path(os.path.realpath(self.home / ".mafigate"))
        self.assertTrue(
            home_value.is_relative_to(anchor),
            f"pyvenv.cfg pins home = {home_value}, a path that moves with the app "
            "bundle; it must point at the stable interpreter copy under ~/.mafigate "
            "(issue #411) or every app move re-breaks the venv",
        )


class TestEveryRecreateBranchClearsTheStaleVenv(unittest.TestCase):
    """The static rule behind the field defect, parse-anchored on launch.sh itself.

    v1.0.0's ``[ ! -x venv/bin/python ]`` branch set ``NEED_NEW_VENV=true`` without
    the ``rm -rf`` its two siblings had, and ``python3 -m venv`` over the leftovers
    skips existing symlinks. So the rule is structural: any branch that decides to
    recreate must clear first. Comments are stripped before matching — a well-worded
    comment must never satisfy (or fail) this guard in place of the code.
    """

    @staticmethod
    def _executable(line: str) -> str:
        stripped = line.strip()
        return "" if stripped.startswith("#") else stripped

    def _branches(self, lines):
        """Each (condition line, index of its NEED_NEW_VENV=true) pair in the file."""
        branches = []
        condition = None
        condition_index = None
        for index, raw in enumerate(lines):
            line = self._executable(raw)
            if not line:
                continue
            if line.startswith(("if ", "elif ")) or line == "else":
                condition = line
                condition_index = index
            if "NEED_NEW_VENV=true" in line:
                branches.append((condition, condition_index, index))
        return branches

    def test_every_branch_that_sets_need_new_venv_clears_the_venv_first(self):
        self.assertTrue(LAUNCH_SH.is_file(), f"the mac launcher moved: {LAUNCH_SH}")
        lines = LAUNCH_SH.read_text().splitlines()
        branches = self._branches(lines)

        # Anti-vacuity: the recreate chain has three arms today (missing python,
        # stamp mismatch, broken interpreter); a parse that stops seeing them must
        # fail loudly here rather than sweep nothing.
        self.assertGreaterEqual(
            len(branches),
            3,
            "the parse found fewer NEED_NEW_VENV=true branches than launch.sh is "
            f"known to have — it is checking nothing: {branches!r}",
        )
        self.assertTrue(
            any(condition and "-x" in condition for condition, _, _ in branches),
            "the parse did not find the missing-interpreter branch "
            "([ ! -x venv/bin/python ]) — the one the field defect lived in",
        )

        uncleared = []
        for condition, start, set_index in branches:
            self.assertIsNotNone(
                condition,
                f"NEED_NEW_VENV=true on line {set_index + 1} has no enclosing "
                "if/elif/else the parse can see",
            )
            body = [self._executable(raw) for raw in lines[start : set_index + 1]]
            if not any('rm -rf "${VENV_DIR}"' in line for line in body):
                uncleared.append(f"line {set_index + 1}: {condition}")
        self.assertEqual(
            uncleared,
            [],
            "these branches decide to recreate the venv without clearing the stale "
            "one first, which is exactly how v1.0.0's second launch died (issues "
            f"#412, #414): {uncleared}",
        )


class TestTheLibrarySeamHolds(unittest.TestCase):
    """The functional tests above stand on this seam; name its loss precisely."""

    def test_launch_sh_defines_ensure_venv_behind_the_lib_guard(self):
        text = LAUNCH_SH.read_text()
        self.assertIn(
            "ensure_venv()",
            text,
            "launch.sh no longer defines ensure_venv(), so the launcher's venv logic "
            "cannot be exercised without booting the app",
        )
        self.assertIn(
            "MAFIGATE_LAUNCH_LIB",
            text,
            "launch.sh lost its MAFIGATE_LAUNCH_LIB seam, so sourcing it for the "
            "functional tests would boot the app for real",
        )
