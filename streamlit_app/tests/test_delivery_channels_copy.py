"""What the README may promise about each delivery channel, and when.

Issue #245. The public repo receives this tree at *moment one*, when the clone route
works and the desktop installers do not exist yet. The README that lands with it is the
only thing telling an external collaborator which channel is theirs, so it has two ways
to be wrong: it can describe an installer nobody can download, and it can make a privacy
promise the software does not yet keep.

Both failures are invisible to every other guard in this suite, because both are prose.
So the guards here tie the prose to a **checkable fact about the tree**, and each one is
a two-state switch rather than a fixed assertion:

* *Installers* — the release workflow that builds them either exists in ``.github/
  workflows`` or does not. While it does not, the README may not link to a release
  artifact, and the table must say the installers are not available yet.
* *The data-posture claim* — a Streamlit config turning usage statistics off either is
  tracked in this repo or is not. While it is not, Streamlit reports usage by default,
  the unqualified *"nothing leaves your machine"* is false, and the table must carry the
  caveat instead. When moment two lands that config the switch flips: the caveat must go
  and the unqualified claim must appear.

That flip is the point. A one-directional assertion would let the caveat outlive the
thing it caveats — copy that was true when written and quietly rots on the branch nobody
re-reads, which is the failure this ticket's canonical-table rule exists to prevent.

**Pinned wording, deliberately.** :data:`UNQUALIFIED_CLAIM` is the exact sentence the
delivery matrix decided on (issue #229), and these tests pin it rather than accepting any
paraphrase. That is a real cost — a maintainer who rewords it gets a red test — and it is
paid on purpose: #229 chose one wording so that the surfaces making this promise cannot
drift apart from each other, and pinning is the only way a test can hold a choice like
that. The sentence itself is *not* in the tree today, and no test here implies it is: two
neighbouring surfaces promise something adjacent in their own words — the Windows
installer's welcome screen and the parameter cache's privacy line — and both are moment
two's to reconcile, not this module's to police.

**Proving the parse.** Every assertion below runs against a table this module has to find
first, so :func:`channel_table` fails loudly when it finds no table, an unrecognised
header, or an empty cell. Prior art is the vendor-drift suite, which fails when it cannot
parse the DMG build script's copy list rather than passing over an empty one — a guard
that silently reads nothing reports success on a tree that has lost the thing it guards.
"""

from __future__ import annotations

import re
import subprocess
import unittest
from pathlib import Path

import pytest

STREAMLIT_APP = Path(__file__).resolve().parents[1]
REPO_ROOT = STREAMLIT_APP.parent

APP_README = STREAMLIT_APP / "README.md"
CONTRIBUTING = REPO_ROOT / "CONTRIBUTING.md"

#: The exact data-posture sentence the delivery matrix settled on. It becomes true — and
#: therefore sayable — only once the telemetry config lands; see the module docstring.
UNQUALIFIED_CLAIM = "nothing leaves your machine"

#: What the caveat has to name while the claim above cannot yet be made. Streamlit's
#: usage reporting is the *only* thing standing between the app and the strong claim, so
#: naming it is what makes the caveat a caveat rather than a hedge.
TELEMETRY_CAVEAT = "usage statistic"

#: The phrases that count as *making* the data-posture claim rather than pointing at it.
#: ``"leaves your machine"`` subsumes :data:`UNQUALIFIED_CLAIM`, and is listed alone so
#: that a near-miss rewording — "your data never leaves your machine" — is caught by the
#: canonical-statement rule too.
DATA_POSTURE_MARKERS = ("leaves your machine", TELEMETRY_CAVEAT)

#: Any of these, in the installers' status cell, states the channel is not ready. Several
#: spellings are accepted because the sentence is a maintainer's to write; what is not
#: negotiable is that the cell says *something* of this kind while no release exists.
NOT_AVAILABLE_YET = (
    "not built yet",
    "not yet built",
    "not available yet",
    "not yet available",
    "coming",
    "no release yet",
)

#: A link into a GitHub release. None of these may appear while no workflow builds one.
RELEASE_ARTIFACT_LINKS = (
    "/releases/download",
    "/releases/latest",
    "/releases/tag/mafigate-v",
)

#: The private repo. The export publishes this tree verbatim, so the literal may not
#: appear in copy this ticket writes. Issue #243 owns the tree-wide sweep; this is the
#: narrow guard over the two files written here, and the two agree by asserting a zero.
#:
#: Assembled from halves rather than written out, because #243's sweep asserts the literal
#: occurs zero times across ``git ls-files``: spelling it here would make this guard the
#: one hit that fails its sibling, which is a silly way for a guard to be the bug.
PRIVATE_REPO = "variantalker" + "_ieo"


def _tracked_files():
    """Every file git tracks, as repo-relative paths.

    Read through git rather than :meth:`Path.rglob` for one reason that has bitten this
    repo: ``.claude/worktrees`` holds full checkouts of the same tree, so a filesystem
    walk from the repo root finds another session's files and answers questions about a
    branch this test is not running on.
    """
    completed = subprocess.run(
        ["git", "ls-files"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=True,
    )
    tracked = [line for line in completed.stdout.splitlines() if line]
    if not tracked:
        raise AssertionError(
            "git ls-files listed nothing — this guard would pass over any tree at all"
        )
    return tracked


#: The three ways Streamlit can be told to stop reporting usage: a config file, the
#: command-line flag the launchers already pass their server address through, and the
#: environment variable. Moment two picks one — the choice is genuinely open, because the
#: config has to be reachable from working directories the four launch routes do not
#: share — so this guard accepts any of them.
#:
#: Accepting all three is what stops the switch jamming shut. If this recognised only a
#: config file and moment two used a flag, the caveat below would be demanded forever
#: against software that no longer needs it, and the flip this module exists to perform
#: would silently never happen. A guard that cannot be satisfied is a guard that will be
#: deleted by whoever meets it.
TELEMETRY_OFF_PATTERNS = (
    r"gatherUsageStats\s*=\s*false",
    r"--browser\.gatherUsageStats[= ]+false",
    r"STREAMLIT_BROWSER_GATHER_USAGE_STATS[= ]+[\"']?(?:false|0)",
)

#: Where a launch route could carry that setting: the two shell entry points, the
#: Makefile, and the bundled launcher.
LAUNCH_ROUTES = (
    "streamlit_app/setup.sh",
    "streamlit_app/run_mafigate.sh",
    "streamlit_app/Makefile",
    "streamlit_app/build/launcher.py",
)


def telemetry_is_turned_off_in_the_tree():
    """Does anything tracked here turn Streamlit's usage reporting off?

    ``browser.gatherUsageStats`` defaults to on, so until something in the tree turns it
    off the strong claim is false on every local channel. Read over the tracked files
    rather than the filesystem, and over all three spellings, for the reasons above each
    constant.
    """
    candidates = []
    for relative in _tracked_files():
        path = Path(relative)
        if path.name == "config.toml" and path.parent.name == ".streamlit":
            candidates.append(path)
        elif relative in LAUNCH_ROUTES:
            candidates.append(path)
    for path in candidates:
        text = (REPO_ROOT / path).read_text(encoding="utf-8", errors="replace")
        if any(re.search(pattern, text, re.IGNORECASE) for pattern in TELEMETRY_OFF_PATTERNS):
            return True
    return False


def a_release_workflow_exists():
    """Does a workflow build the installers on a ``mafigate-v*`` tag?

    This is the tree's own answer to "do the installers exist yet". Moment two adds the
    workflow (#229 decision 5); until then no tag has ever produced an artifact, so any
    README linking to one is linking to a 404.
    """
    workflows = REPO_ROOT / ".github" / "workflows"
    if not workflows.is_dir():
        return False
    return any(
        "mafigate-v" in path.read_text(encoding="utf-8")
        for path in workflows.glob("*.yml")
    )


def _plain(text):
    """Markdown emphasis and code ticks stripped, lowercased, for phrase matching."""
    return re.sub(r"[*`_]", "", text).lower()


def _table_runs(markdown):
    """Every run of consecutive table lines, as ``(line indices, lines)`` pairs.

    The indices come back with the text because two guards here need to know *where* the
    channel table is, not merely what it says: "this claim is made in the table and
    nowhere else" is a statement about position, and answering it by asking whether a
    line starts with a pipe would license the claim inside any of the README's other
    tables — the annotations table, for one.
    """
    runs, current = [], []
    for index, line in enumerate(markdown.splitlines()):
        if line.lstrip().startswith("|"):
            current.append((index, line.strip()))
        elif current:
            runs.append(current)
            current = []
    if current:
        runs.append(current)
    return runs


def _cells(line):
    return [cell.strip() for cell in line.strip().strip("|").split("|")]


def channel_table_span(markdown):
    """The line indices the channel table occupies."""
    for run in _table_runs(markdown):
        if _is_channel_header(run[0][1]):
            return {index for index, _ in run}
    raise AssertionError("no per-channel table found in the app README")


def _is_channel_header(line):
    header = _plain(line)
    return "installer" in header and "clone" in header


def channel_table(markdown):
    """The per-channel table, as ``{row label: {"installers": cell, "clone": cell}}``.

    Fails rather than returns empty when the table is missing or malformed: this module's
    every other assertion reads what comes back from here, so a silent empty dict would
    turn the whole file green against a README that lost its table.
    """
    for indexed_run in _table_runs(markdown):
        run = [line for _, line in indexed_run]
        if not _is_channel_header(run[0]):
            continue
        columns = _cells(run[0])
        if len(columns) != 3:
            raise AssertionError(
                f"the channel table must have a label column and two channels, got: {columns}"
            )
        installer_first = "installer" in _plain(columns[1])
        table = {}
        for line in run[1:]:
            row = _cells(line)
            if len(row) != 3 or all(set(cell) <= {"-", ":", " "} for cell in row):
                continue
            label, first, second = row
            if not label or not first.strip() or not second.strip():
                raise AssertionError(f"the channel table has an empty cell in row: {row}")
            table[_plain(label).strip()] = {
                "installers": first if installer_first else second,
                "clone": second if installer_first else first,
            }
        if len(table) < 3:
            raise AssertionError(
                f"the channel table has too few rows to be canonical: {sorted(table)}"
            )
        return table

    raise AssertionError(
        "no per-channel table found in the app README — it must have a header row naming "
        "both the installers and the clone route"
    )


def _prerequisites_block(markdown):
    """The blockquote stating what the clone route needs, or a failure saying it is gone.

    Fails rather than returning "" for the usual reason: an empty string satisfies no
    assertion in the test that reads it, but only because it satisfies none — the guard
    would go red for the wrong reason and be read as a broken test rather than as a
    README that has lost its prerequisites.
    """
    blocks, current = [], []
    for line in markdown.splitlines():
        if line.lstrip().startswith(">"):
            current.append(line.lstrip().lstrip(">").strip())
        elif current:
            blocks.append(" ".join(current))
            current = []
    if current:
        blocks.append(" ".join(current))
    for block in blocks:
        if "prerequisite" in block.lower():
            return block
    raise AssertionError(
        "the app README has no prerequisites blockquote — the clone route's floor "
        f"(git, a python3, a POSIX shell) is stated nowhere. Blockquotes found: {blocks}"
    )


def _row(table, *keywords):
    """The one row whose label carries any of ``keywords``, or a failure naming what is there."""
    for label, cells in table.items():
        if any(keyword in label for keyword in keywords):
            return cells
    raise AssertionError(
        f"the channel table has no row about {keywords!r}; its rows are {sorted(table)}"
    )


class TheInstallerChannelIsNotPromisedEarly(unittest.TestCase):
    """Nothing may offer a download that no release has produced."""

    def setUp(self):
        self.readme = APP_README.read_text(encoding="utf-8")

    def test_the_table_says_the_installers_are_not_available_yet(self):
        if a_release_workflow_exists():
            pytest.skip("a release workflow exists — the installers are moment two's to describe")
        status = _plain(_row(channel_table(self.readme), "status", "available")["installers"])
        self.assertTrue(
            any(marker in status for marker in NOT_AVAILABLE_YET),
            f"the installers' status cell must say they are not ready yet, got: {status!r}",
        )

    def test_the_readme_links_to_no_release_artifact(self):
        if a_release_workflow_exists():
            pytest.skip("a release workflow exists — release links are legitimate from here on")
        for link in RELEASE_ARTIFACT_LINKS:
            self.assertNotIn(
                link,
                self.readme,
                f"the README links to {link!r}, which no release has produced yet",
            )


class TheDataPostureClaimMatchesTheSoftware(unittest.TestCase):
    """The strongest promise the app makes, held against whether it is currently true.

    Both directions are asserted. The caveat is required while usage reporting is still
    on, and forbidden once anything in the tree turns it off — so moment two cannot fix
    the problem and leave the README still apologising for it.
    """

    def setUp(self):
        self.readme = APP_README.read_text(encoding="utf-8")
        table = channel_table(self.readme)
        self.posture = _row(table, "file goes", "data", "privacy")
        self.cells = _plain(self.posture["installers"] + " " + self.posture["clone"])
        #: Every cell of the table, because the claim is forbidden *anywhere* in it while
        #: it is false. Checking only the posture row would let it be made one row up, in
        #: an audience or status cell, and pass.
        self.whole_table = _plain(
            " ".join(
                cell for row in table.values() for cell in (row["installers"], row["clone"])
            )
        )

    def test_the_claim_is_qualified_while_streamlit_still_reports_usage(self):
        if telemetry_is_turned_off_in_the_tree():
            pytest.skip("telemetry is turned off in the tree — see the test below")
        self.assertIn(
            TELEMETRY_CAVEAT,
            self.cells,
            "nothing tracked turns Streamlit's usage reporting off, so the table's "
            "data-posture row must name it rather than promising more than is true",
        )
        self.assertNotIn(
            UNQUALIFIED_CLAIM,
            self.whole_table,
            f"the table makes the unqualified claim {UNQUALIFIED_CLAIM!r} while Streamlit "
            "still reports usage statistics by default",
        )

    def test_the_caveat_retires_when_telemetry_is_turned_off(self):
        if not telemetry_is_turned_off_in_the_tree():
            pytest.skip("telemetry is not turned off in the tree yet — see the test above")
        self.assertIn(
            UNQUALIFIED_CLAIM,
            self.cells,
            "the tree now turns Streamlit's usage reporting off, so the table may state "
            f"the claim {UNQUALIFIED_CLAIM!r} without qualification",
        )
        self.assertNotIn(
            TELEMETRY_CAVEAT,
            self.cells,
            "the tree now turns usage reporting off, so this caveat is stale "
            "copy and must go",
        )

    def test_the_claim_is_made_in_the_table_and_nowhere_else(self):
        """One canonical statement, per this ticket's rule. Restated copy rots on half
        its branches, so every other surface points at the table instead of repeating it.

        Scoped to the channel table's own lines, not to every line that starts with a
        pipe. The README has other tables — the annotations one, the gene sets — and
        licensing all of them would let the data-posture claim be restated in a table
        about something else entirely, which is the same rot in a different row.
        """
        span = channel_table_span(self.readme)
        for index, line in enumerate(self.readme.splitlines()):
            plain = _plain(line)
            for marker in DATA_POSTURE_MARKERS:
                if marker in plain:
                    self.assertIn(
                        index,
                        span,
                        f"line {index + 1} restates the data posture outside the channel "
                        f"table: {line.strip()!r}",
                    )


class TheCloneRouteIsDescribedTruthfully(unittest.TestCase):
    """The route that actually works at moment one, and who it is for."""

    def setUp(self):
        self.readme = APP_README.read_text(encoding="utf-8")
        self.plain = _plain(self.readme)
        self.table = channel_table(self.readme)

    def test_the_table_names_clone_the_linux_and_contributor_route(self):
        audience = _plain(_row(self.table, "who", "audience", "for")["clone"])
        self.assertIn(
            "linux", audience, "the clone row must name Linux — no installer targets it"
        )
        self.assertTrue(
            "contribut" in audience or "extend" in audience or "code" in audience,
            f"the clone row must name the contributor audience, got: {audience!r}",
        )

    def test_the_table_says_the_clone_route_works_today(self):
        status = _plain(_row(self.table, "status", "available")["clone"])
        self.assertIn(
            "today",
            status,
            f"the clone row must say the route works today, got: {status!r}",
        )

    def test_the_prerequisite_floor_is_stated(self):
        """git, a python3, a POSIX shell; Linux and macOS; Windows only via a bash.

        Read against the prerequisites block, and — for the two that name a program — as
        a **code span**, because two weaker drafts of this test could not fail:

        * searching the whole README for ``"git"`` matched the ``github.com`` in three
          clone URLs, and ``"python3"`` matched the commands in the code fences, so it
          passed with the prerequisites deleted outright;
        * narrowing to the block and matching ``\\bgit\\b`` still passed with ``git``
          deleted, because the same block ends *"run them from a bash: Git Bash or WSL"*
          and ``Git Bash`` is a word-boundary match for git.

        Both were caught by deleting the requirement and watching the test stay green.
        Requiring the backticks is what separates the tool being named as a prerequisite
        from the word appearing in a sentence about something else, and it costs nothing
        this page was not already doing: it writes every program it demands as a code
        span.
        """
        block = _prerequisites_block(self.readme)
        for pattern, requirement in (
            (r"`git`", "`git`"),
            (r"`python3`", "a `python3`"),
            (r"posix shell", "a POSIX shell"),
        ):
            self.assertTrue(
                re.search(pattern, block, re.IGNORECASE),
                f"the prerequisites block must state {requirement} — it reads: {block!r}",
            )
        self.assertTrue(
            "windows" in block.lower() and "bash" in block.lower(),
            f"the prerequisites block must say Windows needs a bash — it reads: {block!r}",
        )


class TheContributingNoteStatesTheInboundFlow(unittest.TestCase):
    """The export is one-way, and that costs a contributor something specific.

    A squashed whole-tree export means a public commit does not survive into the private
    tree, so credit has to be given by hand. #239 wants that documented rather than
    discovered by the first person it happens to.
    """

    def setUp(self):
        self.assertTrue(
            CONTRIBUTING.is_file(),
            "a contributing note must exist at the repo root, where GitHub links it from "
            "the issue and pull-request forms",
        )
        self.text = CONTRIBUTING.read_text(encoding="utf-8")
        self.plain = _plain(self.text)

    def test_it_says_public_issues_are_triaged_into_the_private_tracker(self):
        self.assertIn("issue", self.plain)
        self.assertTrue(
            "triag" in self.plain or "re-filed" in self.plain or "refiled" in self.plain,
            "the note must say what happens to a public issue after it is filed",
        )
        self.assertTrue(
            "stay open" in self.plain or "stays open" in self.plain or "remain open" in self.plain,
            "the note must say public issues stay open, so a reporter is not left guessing "
            "why their issue was never closed",
        )

    def test_it_says_pull_requests_are_reapplied_by_hand(self):
        self.assertIn("pull request", self.plain)
        self.assertTrue(
            "by hand" in self.plain or "manually" in self.plain,
            "the note must say a public pull request is reapplied by hand rather than merged",
        )

    def test_it_says_credit_is_given_manually_and_why(self):
        self.assertTrue(
            "credit" in self.plain,
            "the note must state that contributor credit is given manually",
        )
        self.assertTrue(
            "squash" in self.plain or "does not survive" in self.plain,
            "the note must say why credit is manual — the squash export means a "
            "contributor's commit does not survive",
        )

    def test_it_points_at_the_channel_table_rather_than_restating_it(self):
        self.assertIn(
            "streamlit_app/README.md",
            self.text,
            "the note must link to the app README, which holds the canonical channel table",
        )
        for marker in DATA_POSTURE_MARKERS:
            self.assertNotIn(
                marker,
                self.plain,
                f"the note restates the data posture ({marker!r}); it must point at the table",
            )


class NoNewCopyNamesThePrivateRepo(unittest.TestCase):
    """The two files this ticket writes may not name the repo the public cannot read.

    Issue #243 sweeps the whole tree for the same literal. Two guards asserting the same
    zero is deliberate: this one fails on the files most likely to reintroduce it, and
    keeps failing if the tree-wide sweep is ever narrowed.
    """

    def test_the_app_readme_names_only_the_public_repo(self):
        self.assertNotIn(PRIVATE_REPO, APP_README.read_text(encoding="utf-8"))

    def test_the_contributing_note_names_only_the_public_repo(self):
        self.assertNotIn(PRIVATE_REPO, CONTRIBUTING.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
