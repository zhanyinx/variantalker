#!/usr/bin/env python3
"""The one reader of the pipeline's ``manifest.version``, and of the release tag it derives.

Issue #375, on map #371. The pipeline's version number sat at one value through six
published releases and nothing anywhere noticed — the defect this file exists to make
impossible is not a wrong number but an *unwatched* one. So the number keeps living where
Nextflow reads it, in ``nextflow.config``'s ``manifest {}`` block, and everything that has
an opinion about it reads it through here:

.. code-block:: sh

    python3 tools/manifest_version.py            # prints manifest.version
    python3 tools/manifest_version.py --tag      # prints the release tag it derives

No worked example prints a number in this docstring, deliberately — the same rule
``streamlit_app/build/version.py`` follows. A literal here would be a copy of the version
that ``tools/tests/test_version_contract.py``'s prose sweep is built to forbid, written
inside the tool the sweep exists to serve.

**Stdlib only, and by parsing rather than by asking Nextflow.** ``nextflow config`` would
answer this question authoritatively, and cannot be used: it needs a JVM and a Nextflow
install, and this runs in a CI job that has neither, on a tree whose config *does not parse
under Nextflow's default v2 parser* (measured on map #376 — ``NXF_SYNTAX_PARSER=v1`` is
required to run this pipeline at all). A guard that can only run where the pipeline runs is
a guard that does not run.

**What this does not do.** It reads; it decides nothing. The three policies built on it —
that a pushed tag names this tree's version, that the tree does not sit on an
already-released version, and that no second place states a pipeline version — live in
``tools/tests/test_version_contract.py`` and in ``.github/workflows/version-contract.yml``,
which is where they can fail something.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

#: The release tag namespace, ruled on issue #374: ``"variantalker-v"`` + ``manifest.version``,
#: compared as an exact string.
#:
#: **Fixed text, and it must stay fixed text.** It reads like ``manifest.name`` and is not
#: derived from it on purpose: map #376 owns that field as the pipeline's *display label*, and
#: if the label is ever reworded the tag namespace must not follow it — tags already pushed
#: cannot be renamed, and a namespace that moves with a label would split one release line in
#: two. The prefix is also what lets two products share one repository's release feed, the
#: same way ``mafigate-v`` does for the app (#229, #260).
TAG_PREFIX = "variantalker-v"

#: The config file the manifest lives in, relative to the repository root.
MANIFEST_FILE = Path("nextflow.config")

#: The default repository root: the directory this file's ``tools/`` sits in. Computed from
#: ``__file__`` rather than the working directory, so a workflow step can call this from
#: anywhere.
DEFAULT_REPO_ROOT = Path(__file__).resolve().parents[1]

#: ``manifest`` opening its block, at the start of a line. Anchored so a ``manifest`` written
#: inside a string or a comment elsewhere in the config cannot be mistaken for the block.
MANIFEST_BLOCK = re.compile(r"^\s*manifest\s*\{", re.MULTILINE)

#: ``version = "..."`` inside that block. Single or double quoted, since Nextflow's config
#: takes either and this file must not have an opinion the language does not.
VERSION_ASSIGNMENT = re.compile(
    r"""^\s*version\s*=\s*(?P<quote>["'])(?P<version>.*?)(?P=quote)\s*$""",
    re.MULTILINE,
)


class VersionError(Exception):
    """The version could not be read. Always fatal: nothing here may guess a version."""


def _manifest_block(text: str, path: Path) -> str:
    """The body of ``nextflow.config``'s ``manifest {}`` block.

    Brace-counted rather than read to the first ``}``, so a value that ever contains a brace
    cannot silently truncate the block and hide the assignment below it. Exactly one block
    must exist: two would mean two manifests, and Nextflow's rule for which one wins is not
    a rule this file should be quietly encoding.
    """
    openings = list(MANIFEST_BLOCK.finditer(text))
    if not openings:
        raise VersionError(f"{path} has no manifest {{}} block")
    if len(openings) > 1:
        raise VersionError(
            f"{path} has {len(openings)} manifest {{}} blocks. Which one Nextflow reads is "
            "not something this reader should be deciding on your behalf."
        )

    start = openings[0].end()  # just past the opening brace
    depth = 1
    for index in range(start, len(text)):
        character = text[index]
        if character == "{":
            depth += 1
        elif character == "}":
            depth -= 1
            if depth == 0:
                return text[start:index]

    raise VersionError(f"{path}'s manifest {{}} block is never closed")


def manifest_version(repo_root: Path = DEFAULT_REPO_ROOT) -> str:
    """The value of ``version`` in *repo_root*'s ``nextflow.config`` manifest block.

    Only a quoted literal counts. A version assembled at config time — ``version =
    "1." + minor`` — is not something a guard can compare to a tag without running Groovy,
    and raising beats returning the expression's source text as if it were a version.
    """
    config = Path(repo_root) / MANIFEST_FILE
    if not config.is_file():
        raise VersionError(f"no {MANIFEST_FILE} under {repo_root}")

    block = _manifest_block(config.read_text(encoding="utf-8"), config)

    matches = list(VERSION_ASSIGNMENT.finditer(block))
    if not matches:
        raise VersionError(
            f"{config}'s manifest block has no quoted `version = \"...\"` assignment. Either "
            "the pipeline no longer states a version, or it states one this reader cannot "
            "compare to a tag."
        )
    if len(matches) > 1:
        raise VersionError(
            f"{config}'s manifest block assigns version {len(matches)} times, so the "
            "pipeline's own version depends on which line Nextflow reads last"
        )

    version = matches[0].group("version").strip()
    # Deliberately permissive about the shape — the contract fixes three components today
    # (#372), and a release candidate or build suffix is still a version — and strict about
    # the two states that would make the banner announce nothing at all.
    if not version or not version[0].isdigit():
        raise VersionError(
            f"manifest.version is {matches[0].group('version')!r}, which is not a version. "
            "The pipeline would announce itself with it on every run."
        )
    return version


def release_tag(version: str) -> str:
    """The tag a release of *version* is cut from (#374)."""
    return f"{TAG_PREFIX}{version}"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Print the pipeline's manifest version, or the release tag it derives, "
        "read from nextflow.config's manifest block.",
    )
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=DEFAULT_REPO_ROOT,
        help="the repository root to read nextflow.config from "
        "(default: the one this script lives in)",
    )
    parser.add_argument(
        "--tag",
        action="store_true",
        help="print the release tag instead of the bare version",
    )
    args = parser.parse_args(argv)

    try:
        version = manifest_version(args.repo_root)
    except (VersionError, OSError) as error:
        print(f"error: could not read manifest.version: {error}", file=sys.stderr)
        return 2

    print(release_tag(version) if args.tag else version)
    return 0


if __name__ == "__main__":
    sys.exit(main())
