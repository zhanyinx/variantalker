#!/usr/bin/env python3
"""Writes the build stamp: which channel produced this artifact, and from which commit.

Issue #263. ``build/version.py``'s neighbour, and the same shape for the same reasons — a
build-time writer of a generated, gitignored file that the app then reads. What it stamps is
the pair a version number cannot supply: the *route* the artifact took to a user, and the
commit it was cut from, so a bug report identifies a build rather than a release.

.. code-block:: sh

    python3 build/build_stamp.py --channel macos-dmg        # writes config/build_stamp.py
    python3 build/build_stamp.py --channel windows-installer

Each build script calls it once, before it copies ``config`` (the .dmg) or compiles (the
.exe); ``config/build_identity.py`` reads what it wrote, and reports a source checkout when
it finds nothing. See that module for why the absence is the answer rather than an error.

**Stdlib only, and it never imports the app.** ``build/version.py``'s rule, unchanged: this
runs on a bare build machine before any dependency is installed, and on Windows before the
venv exists. The channel list below is therefore a second copy of
``config.build_identity.INSTALLER_CHANNELS`` rather than an import of it —
``tests/test_build_identity.py`` holds the two equal in both directions, which is the
arrangement the two ``tar --exclude`` lists already live under.

**Where it differs from version.py: a missing ``git`` is a warning, not a failure.** That tool
refuses to guess a version, because the version names the release, the artifact's filename and
the tag — an installer stamped with a guess is worse than no installer. A commit is diagnostic
metadata: someone building from a downloaded archive has no ``.git`` to read, and refusing to
build for want of a bug-report field would trade a working installer for a tidier one. So the
commit is written empty, the app fills in its own word for *not stated*, and this prints a
warning that says so.

A dirty working tree is stamped ``-dirty``, counting tracked modifications only. Untracked
files are excluded deliberately: this repository is routinely checked out beside untracked
reference data, and marking every such build dirty would make the flag mean nothing. What the
flag has to mean is *this artifact does not correspond to a commit anyone can check out*.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

#: The channels a build may stamp. There is no source channel here on purpose: a checkout is
#: what an *absent* stamp means, so stamping one would give that channel two representations
#: and leave the default path — the one every clone and every CI run takes — unexercised.
#: There is no hosted channel either; that deployment was removed and ruled out of scope
#: permanently (#229), so it is not a route an artifact can arrive by.
INSTALLER_CHANNELS = ("macos-dmg", "windows-installer")

#: The default app root: the ``streamlit_app`` directory this file sits under, computed from
#: ``__file__`` rather than the working directory because ``build_installer.bat`` runs from
#: ``build\\windows``. ``build/version.py``'s ``DEFAULT_APP_ROOT``, and the same reason.
DEFAULT_APP_ROOT = Path(__file__).resolve().parents[1]

#: Where the stamp goes, relative to the app root. Inside ``config`` because that is a package
#: both installers already copy whole — the .dmg's copy list names it, the Inno script ships
#: ``config\\*`` with ``recursesubdirs`` — so the stamp reaches the user without either build
#: growing a line about one file, which is a line that gets forgotten.
STAMP_RELATIVE_PATH = Path("config") / "build_stamp.py"

#: How much of the commit hash to write. Short enough to be read aloud off a screen and
#: retyped into an issue, long enough to be unambiguous in a repository this size.
COMMIT_LENGTH = 7


def _git(repo: Path, *arguments: str) -> str | None:
    """``git`` in *repo*, or ``None`` if it could not answer.

    Every way it can fail to answer collapses to one value, because the caller does the same
    thing with all of them: no ``git`` on PATH (a downloaded archive on a fresh machine), not
    a repository (likewise), or a non-zero exit. Nothing is raised — see the module docstring
    on why a missing commit does not stop a build.
    """
    try:
        completed = subprocess.run(
            ["git", "-C", str(repo), *arguments],
            capture_output=True,
            text=True,
            check=False,
        )
    except (OSError, ValueError):
        return None
    if completed.returncode != 0:
        return None
    return completed.stdout.strip()


def short_commit(repo: Path) -> str:
    """The short commit *repo* is checked out at, ``-dirty`` if tracked files are modified.

    ``""`` when git cannot answer. The empty string rather than a placeholder word: the app
    owns the word for *not stated* (``config.build_identity.UNRECORDED_BUILD``) and a second
    spelling of it here is a second thing to keep true.
    """
    commit = _git(repo, "rev-parse", f"--short={COMMIT_LENGTH}", "HEAD")
    if not commit:
        return ""

    # `--untracked-files=no` for the reason the module docstring gives. A failure here reads
    # as clean rather than dirty: `git` answered the harder question already, so treating an
    # unanswered follow-up as dirty would put the flag on builds it is not true of.
    modified = _git(repo, "status", "--porcelain", "--untracked-files=no")
    return f"{commit}-dirty" if modified else commit


def stamp_source(channel: str, commit: str) -> str:
    """The generated module's text: two literals, and a note to whoever opens it.

    Literals rather than anything computed, so the app reads its own identity without running
    a subprocess at import — the app makes no outbound calls and shells out to nothing, and
    a stamp that asked ``git`` on every launch would be a new dependency on the *user's*
    machine having git and this being a checkout.
    """
    return (
        '"""GENERATED FILE — DO NOT EDIT, AND DO NOT COMMIT.\n'
        "\n"
        "Written by build/build_stamp.py when an installer is built (issue #263), and read by\n"
        "config/build_identity.py. Absent from every source checkout, which is how the app\n"
        "knows it is running from one — so a committed copy would make every clone claim to\n"
        "be an installer build, and editing this file changes nothing that survives the next\n"
        'build.\n'
        '"""\n'
        "\n"
        f'BUILD_CHANNEL = "{channel}"\n'
        f'BUILD_COMMIT = "{commit}"\n'
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Write MAFigate's build stamp — the channel an artifact was built for "
        "and the commit it was cut from — for the About dialog to report.",
    )
    parser.add_argument(
        "--channel",
        required=True,
        choices=INSTALLER_CHANNELS,
        help="the route this artifact reaches a user by. Required and closed: a build that "
        "guessed its own channel would tell a bug reporter something unverifiable.",
    )
    parser.add_argument(
        "--app-root",
        type=Path,
        default=DEFAULT_APP_ROOT,
        help="the streamlit_app directory to write the stamp inside, and the tree whose "
        "commit is read (default: the one this script lives in)",
    )
    args = parser.parse_args(argv)

    app_root = Path(args.app_root)
    destination = app_root / STAMP_RELATIVE_PATH

    commit = short_commit(app_root)
    if not commit:
        print(
            "warning: git could not name this checkout's commit, so the build stamp carries "
            "no build identifier. The installer is still complete; a bug report from it will "
            "name its channel and its version but not the commit it was cut from.",
            file=sys.stderr,
        )

    try:
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(stamp_source(args.channel, commit), encoding="utf-8")
    except OSError as error:
        print(f"error: could not write {destination}: {error}", file=sys.stderr)
        return 2

    print(f"{destination}: channel {args.channel}, build {commit or '(none)'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
