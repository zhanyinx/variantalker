"""The repository this tree names out loud, swept as a literal across every tracked file.

MAFigate is published from a *public* repository, ``zhanyinx/variantalker``, and developed in
a private one. Before issue #243 the tree named the private one 62 times in 16 files, in three
classes:

* **two live rendered strings** — the page config's *"Get Help"* and *"Report a bug"* menu
  items, which are the only two a user can click and the only two that 404 for them;
* **~51 provenance references** in docstrings and READMEs, each a Markdown link to a private
  issue. They could not be rewritten to point somewhere public: the public repository has had
  zero issues filed in its life, so no public URL exists for any of them. The URL is gone and
  the bare number stayed as prose — ``issue #199``, not a link — which is what the rest of
  these docstrings already did;
* **nine clone and directory references**, including four in the database tooling's Docker
  README that name the private repo *as a directory*, misnaming the folder a public cloner
  will actually have.

The claim here is a **zero**, over ``git ls-files`` rather than a walk of the filesystem, so
that untracked scratch files and the ignored trees cannot fail it and — more importantly — so
that a file cannot escape the sweep by being new.

**One exclusion, and only one:** ``docs/wayfinder/``. That tree is the measurement record of
the effort itself: reports written *to* the private tracker, quoting its issue numbers as the
subject of the work rather than as a citation inside shipped code. Nothing in it is packaged,
rendered, or read by a user. No other exemption exists — in particular there is no test-tree
exemption, though 30 of the 51 provenance links sat in test and fixture READMEs no user reads.
A guard that can be argued into narrowing is a guard that will be, and this repo has watched
that happen (``tests/README.md``'s note on issue #24).

**What this guard cannot see, and no guard does.**
``build/windows/installer.iss``'s ``AppPublisherURL`` is a placeholder — ``your-org`` where the
account name belongs — so it is broken today, it renders in Windows' Add/Remove Programs, and it
is permanently invisible to a sweep for the *private* name. Widening this file to every GitHub
URL was considered and declined in issue #234: that guard needs exemptions for genuine
third-party links, and an exemption list is the thing this file exists without. So the
placeholder is a copy fact owned by issue #240, and the exclusion named above is the only one
this sweep takes, not the only wrong URL in the tree.

**Why this file does not contain the string it forbids.** It composes the literal from two
halves at import. Spelled out, this module would be its own first violation and the guard
could never be green. The cost is real and worth stating: a human grepping for the private
name will not find the guard that protects them from it. That is what this paragraph is for.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from tests.fakes import page_config_kwargs

STREAMLIT_APP = Path(__file__).resolve().parents[1]
REPO_ROOT = STREAMLIT_APP.parent

#: The forbidden literal, composed rather than written — see the module docstring.
PRIVATE_NAME = "variantalker" + "_ieo"

#: The public repository, and the two places the app sends a user who clicks for help.
PUBLIC_REPO = "https://github.com/zhanyinx/variantalker"
PUBLIC_APP_README = f"{PUBLIC_REPO}/blob/main/streamlit_app/README.md"
PUBLIC_ISSUES = f"{PUBLIC_REPO}/issues"

#: The one excluded tree, as a ``git ls-files`` path prefix.
EXCLUDED_PREFIX = "docs/wayfinder/"


def _tracked_files() -> list[str]:
    """Every tracked path in the repository, relative to its root.

    A failure to run ``git`` is reported as a failure of the guard rather than a reason to
    skip it. A skip here would be indistinguishable from a pass in a CI log, and the whole
    point of this file is a claim that cannot go quiet.
    """
    try:
        completed = subprocess.run(
            ["git", "-C", str(REPO_ROOT), "ls-files", "-z"],
            capture_output=True,
            check=False,
        )
    except OSError as exc:  # pragma: no cover - depends on the machine, not the code
        pytest.fail(
            f"could not run git to list the tracked tree ({exc}). This guard's claim is "
            "about what the repository ships, so it needs the repository; it does not "
            "excuse itself when it cannot check."
        )

    assert completed.returncode == 0, (
        "`git ls-files` failed in "
        f"{REPO_ROOT} (exit {completed.returncode}): "
        f"{completed.stderr.decode('utf-8', 'replace').strip()}"
    )
    out = completed.stdout.decode("utf-8", "replace")
    return [path for path in out.split("\0") if path]


def _hits(path: Path) -> list[tuple[int, str]]:
    """The 1-based line numbers in one file that contain the forbidden literal.

    Read as bytes and decoded permissively, because the tracked tree holds MAFs, a
    spreadsheet, a jar and several PNGs. A decode error must not be able to hide a hit, and
    binary files are not skipped for the same reason: the literal is ASCII and would survive
    into one.
    """
    try:
        text = path.read_bytes().decode("utf-8", "replace")
    except OSError as exc:  # pragma: no cover - depends on the checkout, not the code
        pytest.fail(
            f"cannot read the tracked file {path} ({exc}), so this guard cannot say whether it "
            "names the private repository. An unreadable file is a hole in a claim of zero, "
            "not a file with nothing in it."
        )
    return [
        (number, line.strip())
        for number, line in enumerate(text.splitlines(), start=1)
        if PRIVATE_NAME in line
    ]


@pytest.fixture(scope="module")
def swept() -> list[str]:
    """The paths this guard sweeps: everything tracked, outside the one excluded tree."""
    return [path for path in _tracked_files() if not path.startswith(EXCLUDED_PREFIX)]


def test_the_sweep_reads_a_real_tree(swept):
    """A vacuity check, first, because every claim below is a *count of zero*.

    ``git ls-files`` returning nothing — a detached invocation, a path that is not a
    checkout, a future flag change — would satisfy this file completely while checking
    nothing at all. Two anchors: the list is large, and it contains the app's own entry
    point. Named rather than counted alone, so that a truncated list cannot pass under a
    threshold.
    """
    assert len(swept) > 100, (
        f"the tracked-file sweep found only {len(swept)} paths outside "
        f"{EXCLUDED_PREFIX!r}, which is too few for this repository. The sweep is broken, "
        "and a broken sweep reports zero violations."
    )
    assert "streamlit_app/MAFigate.py" in swept, (
        "the sweep does not include the app's entry point, so whatever it is reading is "
        f"not this repository. {len(swept)} paths were listed."
    )


def test_the_matcher_can_actually_find_the_literal(tmp_path):
    """And the other half of the vacuity check: that ``_hits`` reports a hit it is shown.

    Asserted against a fabricated file rather than against ``docs/wayfinder/`` — where the
    literal genuinely occurs hundreds of times — so that this stays true if that tree is
    ever rewritten too, and so that the machinery is proved by something this test controls.
    """
    planted = tmp_path / "planted.md"
    planted.write_text(
        "a first line with nothing in it\n"
        f"clone https://github.com/zhanyinx/{PRIVATE_NAME}.git\n"
        "a third line\n"
    )

    assert _hits(planted) == [(2, f"clone https://github.com/zhanyinx/{PRIVATE_NAME}.git")]

    clean = tmp_path / "clean.md"
    clean.write_text(f"clone {PUBLIC_REPO}.git\n")
    assert _hits(clean) == []


def test_no_tracked_file_outside_the_wayfinder_names_the_private_repo(swept):
    """The zero itself.

    The message names every hit with its line, because the fix is per-line and a bare count
    would send the next person back to grep — for a string they may not think to compose.
    """
    found = {path: hits for path in swept if (hits := _hits(REPO_ROOT / path))}

    total = sum(len(hits) for hits in found.values())
    report = "\n".join(
        f"  {path}:{number}: {line}"
        for path, hits in sorted(found.items())
        for number, line in hits
    )
    assert not found, (
        f"{total} occurrence(s) of the private repository name in {len(found)} tracked "
        f"file(s) outside {EXCLUDED_PREFIX!r} (issue #243):\n{report}\n"
        "Live strings and clone lines name the public repository; provenance references "
        "drop the URL and keep the bare issue number as prose."
    )


@pytest.fixture
def menu_items():
    """The ``menu_items`` dict the app hands ``st.set_page_config``.

    Driven through the real function rather than read off the source, and driven by the reader
    ``test_app_identity.py`` shares with this module rather than by a second copy of it.
    """
    return page_config_kwargs()["menu_items"]


def test_get_help_lands_on_the_apps_own_readme(menu_items):
    """Not the repository root.

    The root README is the *Nextflow pipeline's* — DSL2 badges, a channel diagram, MAFigate
    one line down its contents — so a straight host swap would land the least technical
    audience somewhere that does not resemble the app they just opened. The app's own README
    is one directory in.
    """
    assert menu_items["Get Help"] == PUBLIC_APP_README, (
        f"'Get Help' points at {menu_items['Get Help']!r}. It must reach the app's own "
        f"README in the public repository: {PUBLIC_APP_README}"
    )


def test_report_a_bug_lands_on_the_public_issues_page(menu_items):
    assert menu_items["Report a bug"] == PUBLIC_ISSUES, (
        f"'Report a bug' points at {menu_items['Report a bug']!r}, which is not the public "
        f"repository's issue tracker: {PUBLIC_ISSUES}"
    )
