"""The pipeline's version identity, and the three ways it can drift.

Issue #375, on map #371. ``nextflow.config``'s ``manifest.version`` held one value across six
published releases — the banner said ``1.0.0`` at public ``v1.0`` in March 2023 and still says
``1.0.0`` at public HEAD — and nothing failed. This module is what makes that impossible to
repeat, and ``.github/workflows/version-contract.yml`` is where it runs: on every push and
every pull request, which is the placement half of the ticket. A check in an unrun file cannot
fail a pull request, which is how issue #355 merged red and sat there.

**The pipeline's number has nothing behind it.** MAFigate stamps every artifact from
``APP_VERSION`` through one reader (#260), so its guard can ask whether two things that should
be equal are. The manifest string *is* the pipeline's source of truth — there is no second
value to compare it against — so this module holds it against the three things that can
disagree with it instead:

1. **A pushed tag naming another version.** ``variantalker-v`` + ``manifest.version``, exact
   string, no normalising and no suffix allowance (#372, #374). The comparison lives in the
   workflow, as a self-contained script over ``$PUSHED_TAG``; this module extracts that script
   and *runs* it against the tags it must refuse. Asserting it by grepping the YAML for
   ``exit 1`` would pass on an inverted comparison.

2. **A version that has already been released.** Set membership, not ordering (#375's ruling):
   a ``variantalker-v<manifest.version>`` tag that exists and does not point at this commit
   means the tree is sitting on a number it has already shipped. This is the assertion that
   would have fired on every pull request for two years, and the cost the ruling accepted is
   that the release act has to end with a bump. It deliberately says nothing about *direction*
   — going backwards is uncaught, because ordering version strings needs a comparison rule the
   contract has not written down (#377).

3. **A second place stating a pipeline version.** #373 measured that no document in this tree
   states one independently, and the sweep below keeps that true. It is worth more than it
   sounds: the number is read in exactly one place in the codebase,
   ``NfcoreTemplate.version()``, so a version written anywhere *else* is a claim that nothing
   derives and nothing can keep honest. Both statements that existed when this landed were
   removed rather than tolerated: a commented-out ``docker://`` container tag in each of the
   two profiles, naming a version older than the manifest's, and ``bin/report.Rmd``'s methods
   sentence, which told every reader of every patient report that the variants were annotated
   with a release this tree has never held.

**No worked example in this module names a version beside the pipeline's name** — not in a
docstring, not in a comment, not in a fixture. The same rule
``streamlit_app/build/version.py``'s docstring states for itself, and here it is enforced
rather than remembered: the sweep reads every tracked file, this file included. Every example
below is either a placeholder or composed from :data:`PIPELINE_NAME` at import time. That was
not foresight. The first version of this module spelled its fixtures out, passed while it was
still untracked, and went red on itself the moment it was committed.

**What this module cannot see.** A version literal with no pipeline name beside it —
``params.pipeline_version = "1.0.0"`` — is indistinguishable from the Annovar, COSMIC or DRAGEN
versions this tree legitimately carries, so the sweep does not look for one. That is the whole
reason the number being read in a single place matters: the sweep guards the *stated* copies,
and the codebase's own shape is what prevents derived ones.

Nor does it compare anything to MAFigate's ``APP_VERSION``. #374 ruled the two numbers
independent and ruled that the contract must say so out loud; an assertion here would be the
opposite claim, and there is no assertion that expresses "these are deliberately unrelated".
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

TESTS_DIR = Path(__file__).resolve().parent
REPO_ROOT = TESTS_DIR.parents[1]

TOOL = REPO_ROOT / "tools" / "manifest_version.py"
CONFIG = REPO_ROOT / "nextflow.config"
WORKFLOW = REPO_ROOT / ".github" / "workflows" / "version-contract.yml"

#: Written out again here rather than imported from ``tools/manifest_version.py``. The
#: deliberate kind of duplication, on the ground ``test_installer_version`` states: a test that
#: imports the value it is checking asserts only that the module equals itself, which is how a
#: renamed namespace stays green.
TAG_PREFIX = "variantalker-v"

#: The pipeline's name as it appears beside a version — ``manifest.name``'s value, spelled out
#: for the same reason as the prefix.
PIPELINE_NAME = "variantalker"

#: The workflow's two jobs and the step this module executes. Named so that renaming or
#: deleting either fails loudly here instead of quietly removing every assertion that reaches
#: through it.
CONTRACT_JOB = "contract"
TAG_JOB = "tag"
TAG_STEP = "Check the pushed tag against manifest.version"

#: A pipeline version written out beside the pipeline's name — ``<name> v<version>``,
#: ``<name>:v<version>``, ``<name>) v<version>``, ``<name>-v<version>``. Matched by *shape*
#: rather than against the number this tree happens to hold, because a sweep pinned to today's
#: literal goes blind the moment the literal is bumped — the hole issue #71's version guard had
#: on the app side.
#:
#: The separator class is what makes this catch the real defects rather than a tidy subset: the
#: two container tags used ``:``, and ``report.Rmd`` reached its version through ``) ``.
VERSION_IN_PROSE = re.compile(
    rf"{PIPELINE_NAME}[^A-Za-z0-9]{{0,3}}v?(\d+\.\d+(?:\.\d+)?)", re.IGNORECASE
)

#: A release tag written out. Held to the derived tag rather than banned outright, so that the
#: export's own release procedure and the contract (#377) can name the tag scheme in prose
#: without this guard and that document contradicting each other — the arrangement
#: ``test_installer_version`` arrived at for ``mafigate-v``. Prose that wants to name the shape
#: rather than a release writes a placeholder, as ``tools/manifest_version.py``'s own docstring
#: does.
TAG_SHAPED = re.compile(rf"{TAG_PREFIX}(\d+\.\d+(?:\.\d+)?)", re.IGNORECASE)

#: The one path prefix the sweep does not read, and why. The prefix below is the private
#: notes, which hold the map's own records — issue #373's surface measurements quote the
#: banner's own output a dozen times because quoting what the pipeline printed *was* the
#: measurement. Those records state nothing about this tree, they are one of
#: ``.publicignore``'s six deny-list rules so they never reach a public reader, and pinning
#: them to today's manifest would rewrite a record of what was true on the day it was taken.
SWEEP_EXCLUSIONS = ("docs/wayfinder/",)

#: Files the sweep must be reading. A guard whose file list quietly emptied — a broken
#: ``git ls-files``, an exclusion that grew a leading ``*`` — would pass over nothing at all,
#: so the three files that have actually carried a pipeline version are named. Paths and
#: counts, never line numbers.
SWEPT_FILES_THAT_HAVE_CARRIED_A_VERSION = (
    "nextflow.config",
    "bin/report.Rmd",
    "README.md",
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _tracked_files() -> list[str]:
    """Every file git tracks, as repo-relative paths.

    Read through git rather than :meth:`Path.rglob`, for the reason
    ``test_installer_version`` and ``test_public_export`` both give: ``.claude/worktrees``
    holds full checkouts of this same tree, so a filesystem walk answers questions about a
    branch this test is not running on. It also keeps ``.nextflow.log``, ``work/`` and a local
    ``annovar/`` download out — untracked run artifacts do carry the banner's version, and none
    of them is a statement this repository makes.
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


def _swept_files() -> list[str]:
    """The tracked files this module reads: everything but :data:`SWEEP_EXCLUSIONS`."""
    return [
        path
        for path in _tracked_files()
        if not any(path.startswith(prefix) for prefix in SWEEP_EXCLUSIONS)
    ]


def _statements(sources: dict[str, str]) -> list[tuple[str, int, str]]:
    """Every pipeline-name-plus-version pair in *sources*, as ``(path, line number, line)``.

    A pure function over path-to-text, so the sweep can be pointed at planted content and
    shown to fail. Tag-shaped matches are left to
    :func:`test_every_release_tag_named_in_prose_is_the_derived_one`, which pins them instead
    of banning them; two guards contradicting each other would be decided by whichever fired
    first.
    """
    found: list[tuple[str, int, str]] = []
    for path, text in sources.items():
        for number, line in enumerate(text.splitlines(), start=1):
            for match in VERSION_IN_PROSE.finditer(line):
                if TAG_SHAPED.fullmatch(match.group(0)):
                    continue
                found.append((path, number, line.strip()))
                break
    return found


def _read_swept_sources() -> dict[str, str]:
    """The swept files' text. Undecodable files are skipped, not guessed at."""
    sources: dict[str, str] = {}
    for path in _swept_files():
        absolute = REPO_ROOT / path
        if not absolute.is_file():  # a submodule entry, or a file deleted but still staged
            continue
        try:
            sources[path] = absolute.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            continue
    return sources


def _tool(*args: str, repo_root: Path | None = None) -> subprocess.CompletedProcess[str]:
    """``tools/manifest_version.py``, run through this interpreter."""
    command = [sys.executable, str(TOOL), *args]
    if repo_root is not None:
        command += ["--repo-root", str(repo_root)]
    return subprocess.run(command, capture_output=True, text=True)


def _tool_output(*args: str, repo_root: Path | None = None) -> str:
    completed = _tool(*args, repo_root=repo_root)
    assert completed.returncode == 0, (
        f"tools/manifest_version.py {' '.join(args)} exited {completed.returncode}:\n"
        f"{completed.stderr}"
    )
    return completed.stdout.strip()


def _throwaway_tree(tmp_path: Path, config: str) -> Path:
    """A repository root holding *config* as its ``nextflow.config``, and nothing else.

    The reader has to read a real config file, so the fake carries the same one name in the
    same place. Everything a pipeline puts around it is irrelevant to a tool reading one
    assignment.
    """
    root = tmp_path / "tree"
    root.mkdir()
    (root / "nextflow.config").write_text(config, encoding="utf-8")
    return root


def _manifest(version_line: str, *, name: str = '    name = "variantalker"\n') -> str:
    """A minimal config whose manifest block carries *version_line*."""
    return "manifest {\n" + name + version_line + "}\n"


def _workflow() -> dict:
    return yaml.safe_load(WORKFLOW.read_text(encoding="utf-8"))


def _tag_step_script() -> str:
    """The tag job's comparison step, extracted from the workflow by name."""
    jobs = _workflow()["jobs"]
    assert TAG_JOB in jobs, (
        f"{WORKFLOW.name} has no {TAG_JOB!r} job. Every assertion about the tag comparison "
        "reaches through it, so its absence must fail here rather than pass silently."
    )
    for step in jobs[TAG_JOB]["steps"]:
        if step.get("name") == TAG_STEP:
            return step["run"]
    raise AssertionError(
        f"{WORKFLOW.name}'s {TAG_JOB} job has no step named {TAG_STEP!r}; its steps are "
        f"{[step.get('name') or step.get('uses') for step in jobs[TAG_JOB]['steps']]}"
    )


def _run_tag_step(pushed_tag: str) -> subprocess.CompletedProcess[str]:
    """The workflow's tag step, executed against *pushed_tag*."""
    return subprocess.run(
        ["bash", "-c", _tag_step_script()],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        env={"PATH": "/usr/bin:/bin:/usr/local/bin", "PUSHED_TAG": pushed_tag},
    )


def _git(*args: str) -> str:
    completed = subprocess.run(
        ["git", *args], cwd=REPO_ROOT, capture_output=True, text=True, check=True
    )
    return completed.stdout.strip()


def _already_released(
    version: str, tags: list[str], tags_at_head: list[str]
) -> str | None:
    """The tag that says *version* has already shipped, or ``None``.

    Set membership on the tag namespace, and nothing else. No ordering, no normalising: the
    contract compares versions as exact strings (#372) and has ruled nothing about which of two
    version strings is later, so this asks only whether this exact version has a tag.

    The exemption for a tag pointing at this commit is what keeps the release commit itself
    green — in both repositories. #374's release act tags private ``HEAD``, exports, then tags
    the resulting public commit with the same name, so on each of those two commits a tag for
    the manifest's own version exists and is correct. It is the *next* commit that has to bump.
    """
    expected = f"{TAG_PREFIX}{version}"
    if expected in tags and expected not in tags_at_head:
        return expected
    return None


# ---------------------------------------------------------------------------
# The reader: one place holds the number, and it really reads it
# ---------------------------------------------------------------------------


def test_the_reader_is_present():
    """Every other test here runs it, so its absence must say so once."""
    assert TOOL.is_file(), (
        f"{TOOL} is missing. It is the only thing that turns manifest.version into a tag, so "
        "without it the workflow's comparison has nothing to compare against."
    )


def test_the_reader_reports_the_manifest_version(tmp_path):
    """And reports it by *reading* it, shown by moving the manifest out from under it.

    The first assertion alone would be satisfied by a tool printing a literal of its own —
    which is the drift this guard exists for, reintroduced one layer down. The second is what
    makes it a derivation.
    """
    # Read independently here rather than through the tool: a different expression than the
    # tool's own, and no brace matching at all.
    #
    # The quote is what makes this work, and it is the trap #373 measured: `nextflow.config`
    # holds a *second* line beginning `version` — `version = false` in the params block, the
    # boolean behind `--version`. A line-level read that ignored quoting would find two
    # assignments and have to guess. Every caller of this pipeline's version wants the quoted
    # literal; the unquoted boolean is a flag that happens to share a name.
    in_config = [
        line.split("=", 1)[1].strip().strip("\"'")
        for line in CONFIG.read_text(encoding="utf-8").splitlines()
        if line.strip().startswith("version") and ('"' in line or "'" in line)
    ]
    assert len(in_config) == 1, (
        f"{CONFIG.name} has {len(in_config)} version assignments at line level; this test "
        "reads the manifest's own without parsing the block, so it needs exactly one"
    )
    assert _tool_output() == in_config[0]

    moved = _throwaway_tree(tmp_path, _manifest('    version = "9.9.9"\n'))
    assert _tool_output(repo_root=moved) == "9.9.9", (
        "the reader did not follow manifest.version to a tree that says 9.9.9, so it is not "
        "reading the manifest — it is repeating a number of its own"
    )


def test_the_reader_derives_the_release_tag(tmp_path):
    assert _tool_output("--tag") == f"{TAG_PREFIX}{_tool_output()}"

    moved = _throwaway_tree(tmp_path, _manifest('    version = "9.9.9"\n'))
    assert _tool_output("--tag", repo_root=moved) == f"{TAG_PREFIX}9.9.9"


#: Every way a config can fail to yield a version this guard could compare to a tag. All of
#: them must be fatal: a reader that returned a plausible string for any of these would let
#: the workflow compare a tag to something the pipeline does not actually announce.
BROKEN_CONFIGS = (
    ("no manifest block", 'params {\n    version = "1.2.0"\n}\n'),
    (
        "two manifest blocks",
        _manifest('    version = "1.2.0"\n') + _manifest('    version = "9.9.9"\n'),
    ),
    ("no version assignment", _manifest("")),
    (
        "two version assignments",
        _manifest('    version = "1.2.0"\n    version = "9.9.9"\n'),
    ),
    ("not a literal", _manifest("    version = params.release\n")),
    ("empty", _manifest('    version = ""\n')),
    ("not a version", _manifest('    version = "unreleased"\n')),
    ("block never closed", 'manifest {\n    version = "1.2.0"\n'),
)


@pytest.mark.parametrize("case,config", BROKEN_CONFIGS, ids=[c for c, _ in BROKEN_CONFIGS])
def test_the_reader_refuses_a_config_it_cannot_read(tmp_path, case, config):
    root = _throwaway_tree(tmp_path, config)
    completed = _tool(repo_root=root)
    assert completed.returncode != 0, (
        f"the reader accepted a config with {case} and printed "
        f"{completed.stdout.strip()!r}. Every caller treats its stdout as the pipeline's "
        "version."
    )
    assert "manifest.version" in completed.stderr


def test_the_reader_refuses_a_tree_with_no_config(tmp_path):
    empty = tmp_path / "empty"
    empty.mkdir()
    assert _tool(repo_root=empty).returncode != 0


# ---------------------------------------------------------------------------
# Policy 1: a pushed tag names this tree's version
# ---------------------------------------------------------------------------


def test_the_workflow_fires_where_it_can_fail_something(tmp_path):
    """Placement, which is half of what issue #375 asked.

    ``on: push`` plus ``on: pull_request`` is the trigger the four sibling contract jobs use,
    and the reason a check belongs here rather than in a test file no workflow names: nothing
    in this repository runs the full suite, so an unrun assertion cannot fail a pull request.
    """
    workflow = _workflow()
    # PyYAML reads an unquoted `on:` key as the boolean True; this is the trap
    # `verifying GitHub workflows offline` names, and the reason for the fallback.
    triggers = workflow.get("on", workflow.get(True))
    assert set(triggers) == {"push", "pull_request"}, (
        f"{WORKFLOW.name} fires on {sorted(triggers)}. A tag-only trigger is the shape this "
        "ticket rejected: it would say nothing about the tree between releases."
    )

    jobs = workflow["jobs"]
    assert CONTRACT_JOB in jobs and TAG_JOB in jobs

    contract = jobs[CONTRACT_JOB]
    runs = " ".join(step.get("run", "") for step in contract["steps"])
    assert "pytest tools/tests" in runs, (
        f"{WORKFLOW.name}'s {CONTRACT_JOB} job does not run this directory, so every "
        "assertion in this module would be a file nothing executes"
    )

    checkout = [step for step in contract["steps"] if "checkout" in str(step.get("uses", ""))]
    assert len(checkout) == 1
    fetched = checkout[0].get("with", {})
    # Without both of these the already-released check reads an empty tag set on every run and
    # passes over it — the vacuous shape this repository has been bitten by. Asserted here
    # rather than trusted to stay in the YAML.
    assert fetched.get("fetch-depth") == 0, (
        "the checkout is shallow, so `git tag` sees nothing and the already-released check "
        "would pass over an empty set on every run"
    )
    assert fetched.get("fetch-tags") is True


def test_the_tag_job_is_gated_on_the_pipeline_namespace():
    """``on: push`` fires for ``mafigate-v*`` too, and that tag says nothing about the manifest."""
    condition = _workflow()["jobs"][TAG_JOB]["if"]
    assert f"refs/tags/{TAG_PREFIX}" in condition, (
        f"the {TAG_JOB} job's condition is {condition!r}, which does not name the pipeline's "
        "tag namespace. Ungated, a MAFigate release tag would fail the pipeline's guard — the "
        "two release lines #374 ruled independent would be failing each other."
    )


def test_the_tag_step_accepts_the_derived_tag():
    completed = _run_tag_step(_tool_output("--tag"))
    assert completed.returncode == 0, (
        f"the tag step refused this tree's own derived tag:\n{completed.stdout}\n"
        f"{completed.stderr}"
    )


def _refused_tags() -> list[tuple[str, str]]:
    """Tags the step must refuse, built from this tree's own version."""
    version = _tool_output()
    tag = f"{TAG_PREFIX}{version}"
    return [
        ("another version", f"{TAG_PREFIX}9.9.9"),
        ("no namespace", f"v{version}"),
        ("bare version", version),
        ("the app's namespace", f"mafigate-v{version}"),
        # #374 fixed the prefix as text: if manifest.name is ever reworded the namespace must
        # not follow it, so a tag built from a reworded label has to be refused.
        ("a reworded label", f"variant-alker-v{version}"),
        # mafigate-release.yml accepts a suffix as a rehearsal because it builds installers.
        # Nothing is built here, so a suffix is only a way to tag a disagreeing tree.
        ("a rehearsal suffix", f"{tag}-rc1"),
        ("empty", ""),
    ]


@pytest.mark.parametrize("case,pushed", _refused_tags(), ids=[c for c, _ in _refused_tags()])
def test_the_tag_step_refuses_a_tag_naming_anything_else(case, pushed):
    completed = _run_tag_step(pushed)
    assert completed.returncode != 0, (
        f"the tag step accepted {pushed!r} ({case}). A release tag naming a version this tree "
        "does not hold points everyone who pins by revision at a pipeline that announces a "
        "different number."
    )
    assert "::error::" in completed.stdout


# ---------------------------------------------------------------------------
# Policy 2: the tree does not sit on an already-released version
# ---------------------------------------------------------------------------


def test_the_clone_can_see_tags():
    """The non-vacuity check for the policy below.

    ``git tag`` on a shallow clone with no tags fetched returns nothing, and a membership test
    over nothing passes. There is no way to tell that state apart from "nothing has been
    released yet" — which is the true state of this namespace today — so the clone itself is
    what gets asserted.
    """
    assert _git("rev-parse", "--is-shallow-repository") == "false", (
        "this is a shallow clone, so `git tag` cannot be trusted to list the release "
        "namespace and the already-released check below would pass over an empty set"
    )


def test_the_manifest_version_is_not_already_released():
    version = _tool_output()
    tags = _git("tag", "--list", f"{TAG_PREFIX}*").splitlines()
    at_head = _git("tag", "--points-at", "HEAD", "--list", f"{TAG_PREFIX}*").splitlines()

    offender = _already_released(version, tags, at_head)
    assert offender is None, (
        f"manifest.version is {version}, and {offender} already exists on a different commit. "
        "This tree is sitting on a version it has already shipped — which is exactly the state "
        "the pipeline sat in for six releases. The release act ends by bumping the manifest to "
        "the version the next release will carry (#372)."
    )


#: A version this repository does not hold, for the fabricated cases below. Every tag literal
#: in this module is *composed* from it and from a namespace constant rather than spelled out,
#: and that is not a style preference: the moment this file became tracked, its own fixtures
#: became statements, and both this module's sweep and MAFigate's went red on them. The same
#: trick ``.publicignore`` uses to name the private repository without naming it.
FABRICATED = "1.2.0"

#: The states the membership rule has to separate, and today none of them can be produced by
#: this repository: the ``variantalker-v`` namespace is empty until the first release is cut.
#: A guard whose only evidence is an empty set has proved nothing, so the rule is exercised
#: against fabricated tag lists — the whole reason it is a pure function.
MEMBERSHIP_CASES = (
    ("nothing released", FABRICATED, [], [], None),
    ("a different version released", FABRICATED, [f"{TAG_PREFIX}1.1.0"], [], None),
    (
        "this version already released",
        FABRICATED,
        [f"{TAG_PREFIX}{FABRICATED}"],
        [],
        f"{TAG_PREFIX}{FABRICATED}",
    ),
    (
        "the release commit itself",
        FABRICATED,
        [f"{TAG_PREFIX}{FABRICATED}"],
        [f"{TAG_PREFIX}{FABRICATED}"],
        None,
    ),
    # The three tag shapes that already exist across the two repositories, none of which is a
    # statement about the pipeline's manifest: this repo's bare line, the retired public line,
    # and the app's namespace.
    ("the bare private line", FABRICATED, ["1.0", "1.1"], [], None),
    ("the retired public line", FABRICATED, ["v1.0", "v1.5"], [], None),
    ("the app's namespace", FABRICATED, [f"mafigate-v{FABRICATED}"], [], None),
)


@pytest.mark.parametrize(
    "case,version,tags,at_head,expected",
    MEMBERSHIP_CASES,
    ids=[case for case, *_ in MEMBERSHIP_CASES],
)
def test_the_membership_rule_separates_the_states(case, version, tags, at_head, expected):
    assert _already_released(version, tags, at_head) == expected, case


# ---------------------------------------------------------------------------
# Policy 3: no second place states a pipeline version
# ---------------------------------------------------------------------------


def test_the_sweep_reads_the_files_that_have_carried_a_version():
    """A file list that quietly emptied would pass over everything.

    Named files and a count, never line numbers: the lesson from issue #370, where a guard's
    own README claim went stale minutes after it was written.
    """
    swept = _swept_files()
    assert len(swept) > 100, f"the sweep is reading only {len(swept)} files"
    for path in SWEPT_FILES_THAT_HAVE_CARRIED_A_VERSION:
        assert path in swept, (
            f"{path} is not in the swept set. It is one of the files that has actually carried "
            "a pipeline version, so a sweep that cannot see it guards nothing."
        )


def test_the_exclusion_is_only_the_maps_own_records():
    """An allow-list is where a guard goes quietly vacuous, so this one is held to its size.

    The emptiness anchor below is unaskable in the exported tree: the excluded prefix is one
    of ``.publicignore``'s own deny rules, so a public clone contains nothing under it by
    design rather than by staleness — and this workflow travels, so without the skip the
    guard fails in the tree it most needs to pass in. Recognised by the deny-list's own
    absence, since the list strips itself; the shape ``test_self_reference.py``'s rule 3 set,
    and it was found the same way — red in a rehearsed export, not by reading.
    """
    assert SWEEP_EXCLUSIONS == ("docs/wayfinder/",)
    tracked = _tracked_files()
    excluded = [path for path in tracked if path.startswith(SWEEP_EXCLUSIONS[0])]
    if not (REPO_ROOT / ".publicignore").exists():
        pytest.skip(
            "no .publicignore in this tree, which is exactly the exported tree — the "
            "excluded prefix is one of its deny rules, so here it matches nothing by design. "
            "Skipping and saying so rather than passing over nothing."
        )
    assert excluded, (
        "the one exclusion matches no tracked file, so it is either stale or misspelt — and a "
        "misspelt exclusion is indistinguishable from a correct one until it stops excluding"
    )


def test_the_sweep_reads_this_module():
    """The guard is inside its own remit, and that is load-bearing.

    A guard module is the most likely place in the tree for a version literal to appear —
    fixtures are *made* of the defect being guarded against — so exempting it would put the
    hole exactly where the material is. Every example here is composed from
    :data:`PIPELINE_NAME` instead, which is only sustainable because this assertion is what
    catches the next person spelling one out.
    """
    assert "tools/tests/test_version_contract.py" in _swept_files()


def test_no_second_place_states_a_pipeline_version():
    found = _statements(_read_swept_sources())
    assert not found, "a pipeline version is stated outside nextflow.config's manifest:\n" + "\n".join(
        f"  {path}:{number}: {line}" for path, number, line in found
    )


def test_the_manifests_own_version_does_not_trip_the_sweep():
    """The sweep reads ``nextflow.config``, and must find the manifest's own line clean.

    A sweep that flagged the source of truth would have to be given an exception for it, and
    that exception is where the two commented container tags would have hidden.
    """
    assert "nextflow.config" in _swept_files()
    assert not _statements({"nextflow.config": CONFIG.read_text(encoding="utf-8")})


#: The statements this repository actually carried, and the shapes they took. Every one is a
#: real line from before this guard landed, which is what makes them a test rather than an
#: imagination of what a defect looks like — and every one is composed from
#: :data:`PIPELINE_NAME` rather than written out, or this table would be the seventh statement
#: and the sweep would fail on the file that defines it. It did, once, which is how the rule
#: got written down.
PLANTED_STATEMENTS = (
    ("a container tag", f'process.container = "docker://yinxiu/{PIPELINE_NAME}:v1.0"'),
    (
        "a commented container tag",
        f'//process.container = "docker://yinxiu/{PIPELINE_NAME}:v1.0"',
    ),
    (
        "a report's methods sentence",
        f"prioritized using [{PIPELINE_NAME}](https://github.com/zhanyinx/{PIPELINE_NAME}) v1.4",
    ),
    (
        "the banner quoted in a doc",
        f"the pipeline prints {PIPELINE_NAME} v1.0.0 on startup",
    ),
    ("a bare pair", f"{PIPELINE_NAME} {FABRICATED}"),
    ("a wrong tag in prose", f"tag it {TAG_PREFIX}9.9.9"),
)


@pytest.mark.parametrize(
    "case,line", PLANTED_STATEMENTS, ids=[case for case, _ in PLANTED_STATEMENTS]
)
def test_the_sweep_finds_a_planted_statement(case, line):
    """The sweep, shown red before it is trusted green.

    A tag-shaped literal is the one case that must reach the *other* guard rather than this
    one — it is pinned to the derived tag, not banned — so it is checked there.
    """
    planted = {"docs/somewhere.md": f"# a heading\n{line}\n"}
    if TAG_SHAPED.search(line):
        assert not _statements(planted)
        assert TAG_SHAPED.search(line).group(0) != f"{TAG_PREFIX}{_tool_output()}"
        return
    assert _statements(planted), (
        f"the sweep did not see {case}: {line!r}. Every line in this table is one this "
        "repository actually shipped."
    )


def test_every_release_tag_named_in_prose_is_the_derived_one():
    """Tag-shaped literals are pinned rather than banned.

    Banning them would put this guard in contradiction with the documents that have to explain
    the tag scheme — the export's own release procedure, and the contract #377 writes — and the
    guard that fired first would decide. Pinning lets prose name the tag and turns red when the
    tree moves
    past it. Prose that means the *shape* rather than a release writes a placeholder, the way
    ``tools/manifest_version.py``'s docstring does.
    """
    expected = f"{TAG_PREFIX}{_tool_output()}"
    wrong: list[str] = []
    for path, text in _read_swept_sources().items():
        for number, line in enumerate(text.splitlines(), start=1):
            for match in TAG_SHAPED.finditer(line):
                if match.group(0) != expected:
                    wrong.append(f"  {path}:{number}: {match.group(0)}")
    assert not wrong, (
        f"a release tag is named in prose that is not this tree's own ({expected}):\n"
        + "\n".join(wrong)
    )
