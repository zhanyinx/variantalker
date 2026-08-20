#!/usr/bin/env python3
"""Refuse to publish a draft release that no longer matches the tree it claims.

Issue #328. `mafigate-v1.0.0` was drafted from a commit that was the tip of public ``main``
when the tag was pushed. Twenty-two commits later it was not, and four of those changed code
that ships — including the fix for a dialog that opened the *wrong variant* (#310). The draft
was still sitting there, correct about its own commit, correctly named, correctly hashed, with
every guard in the suite green. It read exactly as it would have read an hour after being
built, and the handover comment on #265 said to publish it.

That is the gap this closes. The release workflow already refuses the mirror-image mistake — a
rebuild of an already-**published** release, because ``gh release edit --draft`` would retract
it silently — but nothing noticed a draft going stale underneath, and the window is as long as
the draft sits there.

**Why a script and not a paragraph in the build instructions.** The instruction to publish was
written by someone who had just run every check that existed and had them all pass. A step a
human is asked to remember is exactly what failed here, so this refuses instead of reminding.

Run it immediately before pressing Publish::

    make release-preflight                  # the tag build/version.py derives
    make release-preflight TAG=mafigate-v1.0.0-rc1

Exit codes: ``0`` safe to publish, ``1`` refused, ``2`` bad usage, ``3`` could not find out
(no ``gh``, no network, no such tag) — which is deliberately *not* ``0``: "I could not check"
must never read like "I checked".

Stdlib only, and it shells out to ``gh`` rather than importing anything, for the reason
``build/version.py`` gives about itself: a release machine should need nothing installed. The
decision is a pure function — :func:`verdict` — so ``tests/test_release_preflight.py`` can
drive it over states that are otherwise unreachable, this incident among them.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

#: Where the release lives. Written out rather than derived from a remote, because the
#: question this script answers is about *that* repository whichever clone it is run from —
#: and a preflight that read its target from the local git config would silently check the
#: wrong repository on a fork.
DEFAULT_REPO = "zhanyinx/variantalker"

#: Anything but this and the script has not answered the question.
OK, REFUSED, USAGE, UNKNOWN = 0, 1, 2, 3


def verdict(*, tag, tag_commit, branch, branch_commit, behind, is_draft, shipped_changes):
    """Is it safe to publish? Returns ``(exit code, [lines])``.

    Pure, and every input is something the caller had to go and find out — which is what makes
    the states below testable. *behind* is how many commits the branch tip is ahead of the
    tag's commit; *shipped_changes* the paths among them that reach an installer, which decides
    only how the refusal is *worded*, never whether it refuses. Encoding "this drift was
    harmless" is a judgement, and a preflight that makes it can be wrong in the one direction
    that matters.
    """
    lines = [f"tag      {tag} → {tag_commit[:7]}", f"{branch:<9}{branch_commit[:7]}"]

    if is_draft is None:
        return UNKNOWN, lines + [
            f"REFUSED: no release exists for {tag}. Nothing to publish, so nothing was checked."
        ]

    if not is_draft:
        return OK, lines + [
            f"{tag} is already published — this check has nothing left to protect. If you meant"
            " to replace it, delete the release first; re-running the workflow against a"
            " published tag is refused by design."
        ]

    if tag_commit == branch_commit:
        return OK, lines + [
            f"The draft was built from the tip of {branch}. Safe to publish.",
        ]

    lines.append(
        f"REFUSED: the draft is {behind} commit(s) behind {branch}, so it was built from a tree"
        " that no longer exists as the published one."
    )
    if shipped_changes:
        lines.append("")
        lines.append("Code that ships changed in between:")
        lines.extend(f"  {path}" for path in shipped_changes)
        lines.append("")
        lines.append(
            "Publishing this draft ships the older behaviour under the newer version number."
        )
    else:
        lines.append("")
        lines.append(
            "None of those commits touched a path either installer copies, so the artifacts"
            " would probably be identical — but 'probably' is not what a version number means."
        )
    lines.append("")
    lines.append(
        "Delete the draft and its tag, export, re-tag, and let the workflow build it again."
        " Deleting the draft first matters: with it still present the workflow clobbers"
        " artifacts onto it instead of building a fresh one."
    )
    return REFUSED, lines


#: The path prefixes that reach an installer. Both copy the app tree and `resources/`, and
#: neither copies `tests/` — `test_vendor_drift.py` fails an installer that carries a path
#: equal to `tests`, starting `tests/`, or containing `fixtures`. Used only to *word* a
#: refusal, so being approximate here cannot make the script permissive.
SHIPPED_PREFIXES = ("streamlit_app/",)
UNSHIPPED_MARKERS = ("streamlit_app/tests/", "streamlit_app/build/", "/fixtures/")


def shipped(paths):
    """The subset of *paths* that reaches a built installer. See :data:`SHIPPED_PREFIXES`."""
    return sorted(
        path
        for path in paths
        if path.startswith(SHIPPED_PREFIXES)
        and not any(marker in path for marker in UNSHIPPED_MARKERS)
    )


def _gh(*args):
    """``gh`` with *args*, returning stdout, or ``None`` when it could not answer."""
    try:
        completed = subprocess.run(
            ["gh", *args], capture_output=True, text=True, timeout=60
        )
    except (OSError, subprocess.SubprocessError):
        return None
    if completed.returncode != 0:
        return None
    return completed.stdout.strip()


def _derived_tag():
    """``build/version.py --tag``, the one route from ``APP_VERSION`` to a tag name (#260)."""
    version_py = Path(__file__).resolve().parent / "version.py"
    completed = subprocess.run(
        [sys.executable, str(version_py), "--tag"], capture_output=True, text=True
    )
    if completed.returncode != 0:
        return None
    return completed.stdout.strip()


def gather(tag, repo):
    """Everything :func:`verdict` needs, from the GitHub API. ``None`` when it cannot be had."""
    tag_commit = _gh("api", f"repos/{repo}/git/refs/tags/{tag}", "--jq", ".object.sha")
    if not tag_commit:
        return None
    branch = _gh("api", f"repos/{repo}", "--jq", ".default_branch")
    if not branch:
        return None
    branch_commit = _gh(
        "api", f"repos/{repo}/git/refs/heads/{branch}", "--jq", ".object.sha"
    )
    if not branch_commit:
        return None

    # A release may legitimately not exist yet; `gh` failing here is not fatal, so the
    # missing case is represented rather than treated as an outage.
    raw = _gh("api", f"repos/{repo}/releases/tags/{tag}", "--jq", ".draft")
    if raw is None:
        # Draft releases are not served by the by-tag endpoint, so ask the list too.
        listed = _gh("api", f"repos/{repo}/releases", "--jq", ".[] | select(.tag_name==\"" + tag + "\") | .draft")
        raw = listed.splitlines()[0] if listed else None
    is_draft = None if raw in (None, "") else raw.strip() == "true"

    behind, changed = 0, []
    if tag_commit != branch_commit:
        compared = _gh(
            "api", f"repos/{repo}/compare/{tag_commit}...{branch_commit}",
            "--jq", "{ahead: .ahead_by, files: [.files[].filename]}",
        )
        if compared is None:
            return None
        parsed = json.loads(compared)
        behind = parsed.get("ahead") or 0
        changed = parsed.get("files") or []

    return {
        "tag": tag,
        "tag_commit": tag_commit,
        "branch": branch,
        "branch_commit": branch_commit,
        "behind": behind,
        "is_draft": is_draft,
        "shipped_changes": shipped(changed),
    }


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Refuse to publish a draft release built from a tree that has moved on.",
        epilog="Run it immediately before pressing Publish. See build/BUILD_INSTRUCTIONS.md.",
    )
    parser.add_argument("--tag", help="the release tag (default: build/version.py --tag)")
    parser.add_argument("--repo", default=DEFAULT_REPO, help=f"default: {DEFAULT_REPO}")
    args = parser.parse_args(argv)

    tag = args.tag or _derived_tag()
    if not tag:
        print(
            "could not derive the tag from build/version.py; pass --tag", file=sys.stderr
        )
        return USAGE

    facts = gather(tag, args.repo)
    if facts is None:
        print(
            f"could not read {args.repo} through `gh` — no network, not authenticated, or no"
            f" tag {tag}. This is not a pass: nothing was checked.",
            file=sys.stderr,
        )
        return UNKNOWN

    code, lines = verdict(**facts)
    print("\n".join(lines))
    return code


if __name__ == "__main__":
    sys.exit(main())
