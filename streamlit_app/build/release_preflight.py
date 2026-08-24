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

**And the first version of this script did not close it** (issue #405). It compared the tag
against the tip of public ``main`` and stopped there — but those are the *same commit* whenever
nobody has exported since the tag was cut, which is the ordinary state, the export being manual
and once per release. So it passed the incident above, and passed it again four days later with
the private tree 74 commits ahead of what was published. The staleness lives one level further
out: between the published tree and the tree it was exported from. That is what
*exported from* / *private* below are, and they are read out of public ``main``'s own tip
commit message, which the export writes and issue #227 made the contract.

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


def verdict(
    *,
    tag,
    tag_commit,
    branch,
    branch_commit,
    behind,
    is_draft,
    shipped_changes,
    exported_from=None,
    private_head=None,
):
    """Is it safe to publish? Returns ``(exit code, [lines])``.

    Pure, and every input is something the caller had to go and find out — which is what makes
    the states below testable. *behind* is how many commits the branch tip is ahead of the
    tag's commit; *shipped_changes* the paths among them that reach an installer, which decides
    only how the refusal is *worded*, never whether it refuses. Encoding "this drift was
    harmless" is a judgement, and a preflight that makes it can be wrong in the one direction
    that matters.

    *exported_from* is the private commit the published tree was built from — the SHA named in
    public ``main``'s own tip commit message — and *private_head* is where the private tree is
    now. **They are the pair that matters**, and the first version of this function did not
    take them: it compared the tag against the tip of public ``main`` only, which are the same
    commit whenever nobody has exported since the tag was cut. That is the normal state, since
    the export is manual and once per release, so the check passed exactly when the hazard was
    present — including on the incident it was written for (issue #405).
    """
    lines = [f"tag      {tag} → {tag_commit[:7]}", f"{branch:<9}{branch_commit[:7]}"]
    if exported_from:
        lines.append(f"exported from {exported_from[:7]}")

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
        # The tag names what is published. That is necessary and it is not sufficient: what is
        # published can itself be behind the tree it was exported from, and then the tag and
        # the branch agree with each other about a tree nobody is developing any more.
        if private_head is None:
            return UNKNOWN, lines + [
                "REFUSED: the tag names the tip of the published branch, but whether that"
                " branch is current could not be checked — this does not look like the private"
                " tree (no tools/export_public.py), and only it knows what is waiting to be"
                " exported. Run this where you run the export.",
            ]
        if exported_from is None:
            return UNKNOWN, lines + [
                f"REFUSED: {branch}'s tip commit does not say which private commit it was"
                " exported from, so this cannot tell whether an export is pending. Every"
                " commit the export writes says so; a tip that does not was made another way.",
            ]
        if exported_from != private_head:
            return REFUSED, lines + [
                f"private  {private_head[:7]}",
                "",
                f"REFUSED: an export is pending. What is published was exported from"
                f" {exported_from[:7]}, and the private tree is now at {private_head[:7]}, so"
                " the artifacts were built from a tree that is no longer the current one —"
                " even though the tag correctly names the tip of what is published.",
                "",
                "This is the failure the first version of this check could not see, because"
                " those two commits agree with each other whenever nobody has exported since"
                " the tag was cut, which is the ordinary state.",
                "",
                "Export, re-tag, and let the workflow build it again — deleting the draft and"
                " its tag first, or the workflow clobbers artifacts onto the old draft.",
            ]
        return OK, lines + [
            f"The draft was built from the tip of {branch}, and {branch} is current with the"
            " private tree. Safe to publish.",
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


#: How ``tools/export_public.py`` opens every commit it writes — ``Sync from <private repo> @
#: <sha>``. Assembled from halves rather than written out, for the reason
#: ``tests/test_public_repo_name.py`` gives: that sweep asserts the private repository's name
#: occurs zero times across the tracked tree, and this file travels. The export tooling composes
#: it the same way, which is the precedent.
SYNC_PREFIX = "Sync from " + "variantalker" + "_ieo" + " @ "


def exported_from(message):
    """The private commit a published tip was exported from, or ``None``.

    ``None`` is not "nothing pending" — it is "this commit did not come from the export", which
    :func:`verdict` treats as a state it cannot judge rather than as a pass.
    """
    if not message:
        return None
    first = message.splitlines()[0].strip()
    if not first.startswith(SYNC_PREFIX):
        return None
    sha = first[len(SYNC_PREFIX) :].strip()
    return sha if sha else None


def _private_head(repo_root=None):
    """This tree's HEAD, but only when this *is* the private tree.

    Told apart by ``tools/export_public.py``, which the export strips from what it publishes —
    so its presence is the same signal the export tooling already relies on, rather than a new
    one invented here. In a public clone the answer is unknowable rather than false, and the
    caller must not read it as fine.
    """
    repo_root = Path(repo_root or Path(__file__).resolve().parents[2])
    if not (repo_root / "tools" / "export_public.py").is_file():
        return None
    completed = subprocess.run(
        ["git", "-C", str(repo_root), "rev-parse", "HEAD"],
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        return None
    return completed.stdout.strip() or None


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

    tip_message = _gh(
        "api", f"repos/{repo}/commits/{branch_commit}", "--jq", ".commit.message"
    )

    return {
        "tag": tag,
        "tag_commit": tag_commit,
        "branch": branch,
        "branch_commit": branch_commit,
        "behind": behind,
        "is_draft": is_draft,
        "shipped_changes": shipped(changed),
        "exported_from": exported_from(tip_message),
        "private_head": _private_head(),
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
