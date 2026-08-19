#!/usr/bin/env python3
"""The one reader of ``APP_VERSION``, for everything that has to stamp a version on a file.

Issue #260. A version number was written in five places — the app's constant, the Inno
Setup script twice, the DMG script's argument default, the macOS bundle's ``Info.plist``
twice, and the Makefile — and they had already drifted by a whole major version: the app
reported one number and every installer stamped a lower one. Nothing reconciled them, and
moment two puts that number on an artifact a clinician downloads and on a release tag.

So the number lives in ``config/constants.py`` and reaches a build only through this file:

.. code-block:: sh

    python3 build/version.py                    # prints APP_VERSION
    python3 build/version.py --tag              # prints mafigate-v<version>
    python3 build/version.py --write-iss        # writes build/windows/version.iss

No worked example prints a number here, deliberately. A literal in this file's own docstring
would be a copy the guard's prose sweep cannot see — it looks for a version next to the word
*MAFigate* — which is the exact shape of hole this file exists to close.

**Stdlib only, and by parsing rather than importing.** The same rule ``vendor/_sync.py``
follows and for the same reason: this runs on a bare build machine before any dependency is
installed, and on Windows before the venv exists. Importing ``config.constants`` would work
today — the module is nothing but literals — but it would make the build's version depend
on the app's import graph staying that way, and one ``import streamlit`` at the top of a
future ``constants.py`` would break the installer rather than the app. :func:`ast.parse`
reads the assignment and nothing else runs.

**Why the tag has its own namespace.** The public repository's ``v1.x`` tags belong to the
Nextflow pipeline, released on its own cadence (latest 2025-01-31). A MAFigate release
sharing that namespace would put a second product in the release feed pipeline users watch,
and a pipeline bugfix would read as a new MAFigate. Hence ``mafigate-v`` (#229, #260).

``tests/test_installer_version.py`` holds every installer input against this file, and
proves this file *derives* rather than repeats by running it against a throwaway tree whose
``APP_VERSION`` says something else.
"""

from __future__ import annotations

import argparse
import ast
import sys
from pathlib import Path

#: The release tag namespace. See the module docstring for why it is not ``v``.
TAG_PREFIX = "mafigate-v"

#: The default app root: the ``streamlit_app`` directory this file sits under. Computed
#: from ``__file__`` rather than from the working directory, so a build script can call it
#: from anywhere — ``build_installer.bat`` runs from ``build\windows``.
DEFAULT_APP_ROOT = Path(__file__).resolve().parents[1]

#: Where the generated Inno Setup include goes, beside the script that includes it.
DEFAULT_ISS_PATH = Path(__file__).resolve().parent / "windows" / "version.iss"

CONSTANT_NAME = "APP_VERSION"


class VersionError(Exception):
    """The constant could not be read. Always fatal: a build must not guess a version."""


def app_version(app_root: Path = DEFAULT_APP_ROOT) -> str:
    """The value of ``APP_VERSION`` in *app_root*'s ``config/constants.py``.

    Only a module-level assignment of a plain string counts. A version built by calling
    something — ``APP_VERSION = read_version()`` — is not a literal this can read, and
    raising beats returning the expression's source text as if it were a version.
    """
    constants = Path(app_root) / "config" / "constants.py"
    if not constants.is_file():
        raise VersionError(f"no config/constants.py under {app_root}")

    tree = ast.parse(constants.read_text(encoding="utf-8"), filename=str(constants))
    for node in tree.body:
        if not isinstance(node, ast.Assign):
            continue
        names = [target.id for target in node.targets if isinstance(target, ast.Name)]
        if CONSTANT_NAME not in names:
            continue
        if not (isinstance(node.value, ast.Constant) and isinstance(node.value.value, str)):
            raise VersionError(
                f"{CONSTANT_NAME} in {constants} is not a string literal, so no build can "
                "read it without running the app"
            )
        version = node.value.value.strip()
        # Deliberately permissive about the shape — a release candidate or a build suffix
        # is a legitimate version — and strict about the two states that would silently
        # produce an unversioned artifact.
        if not version or not version[0].isdigit():
            raise VersionError(
                f"{CONSTANT_NAME} is {node.value.value!r}, which is not a version. An "
                "installer stamped with it would name no release at all."
            )
        return version

    raise VersionError(f"{constants} has no module-level {CONSTANT_NAME} assignment")


def release_tag(version: str) -> str:
    """``mafigate-v<version>`` — the tag a release is cut from."""
    return f"{TAG_PREFIX}{version}"


def inno_include(version: str) -> str:
    """The Inno Setup preprocessor include that hands ISCC the version.

    Generated rather than committed: a tracked copy of this file would be the fifth place
    the number is written, and — being generated — the one nobody would think to update.
    """
    return (
        "; GENERATED FILE — DO NOT EDIT, AND DO NOT COMMIT.\n"
        "; Written by build/version.py from config/constants.py's APP_VERSION.\n"
        "; installer.iss includes it; build_installer.bat regenerates it before every\n"
        "; compile, so editing this file changes nothing that survives the next build.\n"
        f'#define AppVersion "{version}"\n'
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Print MAFigate's version, its release tag, or generate the "
        "installer's version include — all from config/constants.py.",
    )
    parser.add_argument(
        "--app-root",
        type=Path,
        default=DEFAULT_APP_ROOT,
        help="the streamlit_app directory to read config/constants.py from "
        "(default: the one this script lives in)",
    )
    output = parser.add_mutually_exclusive_group()
    output.add_argument(
        "--tag",
        action="store_true",
        help=f"print the release tag ({TAG_PREFIX}<version>) instead of the version",
    )
    output.add_argument(
        "--write-iss",
        nargs="?",
        const=DEFAULT_ISS_PATH,
        type=Path,
        metavar="PATH",
        help="write the Inno Setup version include (default: build/windows/version.iss)",
    )
    args = parser.parse_args(argv)

    try:
        version = app_version(args.app_root)
    except (VersionError, SyntaxError, OSError) as error:
        print(f"error: could not read {CONSTANT_NAME}: {error}", file=sys.stderr)
        return 2

    if args.write_iss is not None:
        destination = Path(args.write_iss)
        try:
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_text(inno_include(version), encoding="utf-8")
        except OSError as error:
            print(f"error: could not write {destination}: {error}", file=sys.stderr)
            return 2
        print(f"{destination}: AppVersion {version}")
        return 0

    print(release_tag(version) if args.tag else version)
    return 0


if __name__ == "__main__":
    sys.exit(main())
