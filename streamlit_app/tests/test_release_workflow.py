"""The workflow that cuts a release: one tag, two runners, one commit.

Issue #264. Cutting a release used to need a Mac *and* a physical Windows box, both driven
by hand, and the failure mode was asymmetric and silent — with the Windows machine not to
hand, Mac users got the new version and Windows users kept the old one without being told.
``.github/workflows/mafigate-release.yml`` replaces both machines with a ``macos-latest`` +
``windows-latest`` matrix on a ``mafigate-v*`` tag.

**Why this module runs the workflow's shell instead of reading its YAML.** Every claim worth
making here is about behaviour: that a tag naming the wrong version is refused, that two
artifacts built from two different commits are refused, that an empty commit report is
refused rather than passing as *three equal strings*. Asserting those by grepping the YAML
for ``exit 1`` would pass on a workflow whose comparison had been inverted. So the two steps
that decide anything are written as self-contained scripts reading environment variables, and
this module extracts them by step name and executes them — the same discipline
``test_installer_version`` uses when it proves ``build/version.py`` derives by running it
against a throwaway tree.

What is left as a text assertion, and why each is honest that way:

* The Windows build step's ``< NUL``. ``build_installer.bat`` ends in ``pause``, so without
  it the job waits for a keypress until the six-hour runner timeout. The failure is a hang,
  which cannot be provoked here — but its cause is one token.
* The unversioned-compile check. ``installer.iss``'s ``#error`` is the one assertion
  ``test_installer_version`` says it cannot make, "the compile that would prove it belongs to
  the release workflow (#264), which builds on Windows". This module holds the workflow to
  making it; ISCC is what makes it.
* ``--draft``. This workflow never publishes. It attaches both artifacts to a draft and
  stops, so #265 can author the page and a human decides when a download exists — which is
  also what keeps ``build/RELEASES.md`` a record of a human act rather than of CI's.

Not asserted here: that the builds work. That is what the rehearsal tag is for, and no test
in this suite can run ISCC or ``hdiutil``.
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

TESTS_DIR = Path(__file__).resolve().parent
STREAMLIT_APP = TESTS_DIR.parent
REPO_ROOT = STREAMLIT_APP.parent

WORKFLOW = REPO_ROOT / ".github" / "workflows" / "mafigate-release.yml"

#: The tag namespace, spelled out rather than imported from ``build/version.py`` for the
#: reason ``test_installer_version`` gives: a test that imports the value it checks asserts
#: only that the module equals itself.
TAG_PREFIX = "mafigate-v"

#: The jobs this module has opinions about. Named so a renamed or deleted job fails loudly
#: here instead of quietly removing every assertion that reached through it.
JOBS = ("version", "macos", "windows", "release")

#: Steps this module extracts and runs, or reads. Same reason: a renamed step would
#: otherwise take its assertions with it.
DERIVE_STEP = "Derive the version and check the tag against it"
COMMIT_STEP = "Check both artifacts came from one commit"
WINDOWS_BUILD_STEP = "Build the installer"
UNVERSIONED_COMPILE_STEP = "A compile with no version supplied still refuses"


def _document():
    """The workflow, parsed.

    ``on:`` is read through both spellings on purpose: YAML 1.1 reads a bare ``on`` as the
    boolean ``True``, so PyYAML hands back ``{True: {...}}`` and a lookup of ``"on"`` finds
    nothing. A guard that quietly found nothing would assert nothing about the trigger, which
    is the one part of this workflow that decides how often macOS runners are spent.
    """
    assert WORKFLOW.is_file(), (
        f"{WORKFLOW.relative_to(REPO_ROOT).as_posix()} is missing. It is the workflow that "
        "builds both installers on a tag; without it a release needs two machines driven by "
        "hand."
    )
    return yaml.safe_load(WORKFLOW.read_text(encoding="utf-8"))


DOCUMENT = _document()


def _triggers():
    triggers = DOCUMENT.get("on", DOCUMENT.get(True))
    assert isinstance(triggers, dict), (
        f"the workflow's triggers did not parse as a mapping: {triggers!r}. See _document "
        "for the `on:` / `True` trap."
    )
    return triggers


def _job(name):
    jobs = DOCUMENT.get("jobs") or {}
    assert name in jobs, (
        f"the workflow has no {name!r} job, so every assertion this module makes through it "
        f"is asserting nothing. Its jobs are {sorted(jobs)}."
    )
    return jobs[name]


def _steps(job_name):
    return _job(job_name)["steps"]


def _step(job_name, step_name):
    """One step of a job, by its ``name``, or a failure listing the names there are."""
    for step in _steps(job_name):
        if step.get("name") == step_name:
            return step
    raise AssertionError(
        f"the {job_name!r} job has no step named {step_name!r}; it has "
        f"{[step.get('name') or step.get('uses') for step in _steps(job_name)]}"
    )


def _step_index(job_name, step_name):
    for index, step in enumerate(_steps(job_name)):
        if step.get("name") == step_name:
            return index
    raise AssertionError(f"no step named {step_name!r} in {job_name!r}")


def _run(script, env, cwd=STREAMLIT_APP):
    """Execute one of the workflow's own ``run:`` scripts under bash."""
    return subprocess.run(
        ["bash", "-c", script],
        cwd=cwd,
        env={"PATH": "/usr/bin:/bin:/usr/local/bin:/opt/homebrew/bin", **env},
        capture_output=True,
        text=True,
    )


def _outputs(path):
    """A ``$GITHUB_OUTPUT`` file, as a dict."""
    values = {}
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            values[key] = value
    return values


def _derive(pushed_tag, tmp_path):
    """Run the version job's derive step against *pushed_tag*."""
    output = tmp_path / "output"
    output.write_text("", encoding="utf-8")
    summary = tmp_path / "summary"
    completed = _run(
        _step("version", DERIVE_STEP)["run"],
        {
            "PUSHED_TAG": pushed_tag,
            "GITHUB_OUTPUT": str(output),
            "GITHUB_STEP_SUMMARY": str(summary),
        },
    )
    return completed, _outputs(output)


def _check_commits(**env):
    return _run(_step("release", COMMIT_STEP)["run"], env)


# ---------------------------------------------------------------------------
# Proving the parse
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", JOBS)
def test_the_workflow_has_the_job_this_module_reaches_through(name):
    """Vacuity check: every assertion below reaches a job through :func:`_job`."""
    assert _job(name)


def test_the_two_decisive_steps_are_scripts_this_module_can_run():
    """The other half of the parse proof: both scripts are present and non-trivial.

    Length rather than content, because content is what the behavioural tests assert. A
    step reduced to ``true`` would satisfy every "exits 0" case below and nothing else.
    """
    for job, step_name in (("version", DERIVE_STEP), ("release", COMMIT_STEP)):
        script = _step(job, step_name).get("run", "")
        assert len(script.splitlines()) >= 5, (
            f"{job}/{step_name!r} is {script!r}, which is too short to be deciding anything"
        )


# ---------------------------------------------------------------------------
# When it runs, and on what
# ---------------------------------------------------------------------------


def test_only_a_release_tag_starts_it():
    """A tag push and nothing else — which is the runner-cost answer as well as the trigger.

    macOS runners bill at a multiplier on private repositories, so this workflow must not be
    reachable from an ordinary push or pull request; the two existing guards are the jobs
    that run on those. It also means the rehearsal an acceptance criterion asks for happens
    where it is free: pushing a throwaway ``<tag>-rc1`` to the public repository.
    """
    triggers = _triggers()
    assert set(triggers) == {"push"}, (
        f"the release workflow triggers on {sorted(triggers)}. Only a tag push may cut a "
        "release — everything else spends macOS minutes on every commit."
    )
    push = triggers["push"] or {}
    assert push.get("tags") == [f"{TAG_PREFIX}*"], (
        f"the tag filter is {push.get('tags')!r}; it must be exactly ['{TAG_PREFIX}*'] — the "
        "namespace that keeps MAFigate releases out of the Nextflow pipeline's live v1.x line"
    )
    assert "branches" not in push, "a branch push must not cut a release"


def test_one_job_builds_on_macos_and_one_on_windows():
    """The matrix, which is the whole point: two artifacts, no physical machines."""
    assert _job("macos")["runs-on"].startswith("macos-"), _job("macos")["runs-on"]
    assert _job("windows")["runs-on"].startswith("windows-"), _job("windows")["runs-on"]


def test_the_release_job_waits_for_all_three():
    needs = _job("release")["needs"]
    assert set(needs) == {"version", "macos", "windows"}, (
        f"the release job needs {needs}; it must wait for both builds and the version check, "
        "or it will draft a release with one artifact or none"
    )


def test_only_the_release_job_may_write_to_the_repository():
    """Least privilege, and it is the release job's ability to attach that is being granted."""
    assert DOCUMENT.get("permissions") == {"contents": "read"}, DOCUMENT.get("permissions")
    assert _job("release").get("permissions") == {"contents": "write"}, (
        "the release job must ask for contents: write — it creates the release and uploads "
        "both artifacts"
    )
    for name in ("version", "macos", "windows"):
        assert "permissions" not in _job(name), (
            f"the {name} job asks for permissions of its own; it only reads the tree"
        )


# ---------------------------------------------------------------------------
# The version and the tag come from the constant, and disagreement is fatal
# ---------------------------------------------------------------------------


def test_the_derive_step_reads_the_one_reader_of_the_constant():
    """Positive half of "no version is typed into the workflow".

    The negative half — that no file under ``.github/workflows`` spells a MAFigate version
    out — is swept by ``test_installer_version``, along with every other build input, so it
    covers a second workflow nobody has written yet.
    """
    script = _step("version", DERIVE_STEP)["run"]
    assert "build/version.py" in script, (
        "the derive step must read the version through build/version.py, the one route from "
        f"APP_VERSION to a build. It reads: {script!r}"
    )
    assert "--tag" in script, "the release tag must be derived too, not assembled here"


def test_the_derive_step_accepts_the_tag_the_constant_names(tmp_path):
    from config.constants import APP_VERSION

    tag = f"{TAG_PREFIX}{APP_VERSION}"
    completed, outputs = _derive(tag, tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert outputs["version"] == APP_VERSION, outputs
    assert outputs["tag"] == tag, outputs
    assert outputs["prerelease"] == "false", outputs
    assert re.fullmatch(r"[0-9a-f]{40}", outputs["commit"]), (
        f"the step must report the commit it is building, got {outputs.get('commit')!r}"
    )


def test_the_derive_step_marks_a_suffixed_tag_a_prerelease(tmp_path):
    """The rehearsal route, and the reason it is safe to rehearse on the real workflow.

    An acceptance criterion asks that this workflow's first run not be the one that has to
    produce a real artifact. ``<derived tag>-rc1`` is the same version deliberately not
    claiming to be the release, so a rehearsal cannot be mistaken for one.
    """
    from config.constants import APP_VERSION

    completed, outputs = _derive(f"{TAG_PREFIX}{APP_VERSION}-rc1", tmp_path)

    assert completed.returncode == 0, completed.stderr
    assert outputs["prerelease"] == "true", outputs
    assert outputs["version"] == APP_VERSION, outputs
    assert outputs["tag"] == f"{TAG_PREFIX}{APP_VERSION}-rc1", (
        "the release must be cut from the tag that was pushed, not from the derived one"
    )


@pytest.mark.parametrize(
    "pushed",
    [
        "mafigate-v9.9.9",  # a version this tree does not hold
        "mafigate-v0.1",  # nor this one
        "v1.0.0",  # the Nextflow pipeline's namespace
        "mafigate-vsomething",
        "mafigate-v",
    ],
)
def test_the_derive_step_refuses_a_tag_the_constant_does_not_name(pushed, tmp_path):
    """A tag is a human typing a version, so it is checked against the constant.

    Without this the workflow would happily build ``mafigate-v9.9.9`` out of a tree whose app
    reports something else — an artifact whose filename, About dialog and release page
    disagree, which is the exact class of failure #260 removed everywhere else.
    """
    completed, outputs = _derive(pushed, tmp_path)

    assert completed.returncode != 0, (
        f"{pushed!r} was accepted; it names no version this tree holds.\n{completed.stdout}"
    )
    assert not outputs, f"a refused tag still wrote outputs: {outputs}"
    assert TAG_PREFIX in completed.stdout + completed.stderr, (
        "the refusal must say what the tag should have been — the next thing whoever pushed "
        "it needs is the derived tag"
    )


def test_a_suffix_is_not_a_way_past_the_version_check(tmp_path):
    """``mafigate-v9.9.9-rc1`` is a wrong version with a suffix, not a rehearsal."""
    completed, _ = _derive("mafigate-v9.9.9-rc1", tmp_path)
    assert completed.returncode != 0, completed.stdout


# ---------------------------------------------------------------------------
# One commit, proved rather than assumed
# ---------------------------------------------------------------------------


def test_each_build_job_reports_the_commit_it_built():
    for name in ("macos", "windows"):
        outputs = _job(name).get("outputs") or {}
        assert "commit" in outputs, (
            f"the {name} job publishes no commit output, so nothing downstream can compare "
            f"what the two artifacts were built from. Its outputs: {sorted(outputs)}"
        )


def test_the_one_commit_check_passes_when_all_three_agree():
    completed = _check_commits(
        TAG_COMMIT="a" * 40, MAC_COMMIT="a" * 40, WINDOWS_COMMIT="a" * 40
    )
    assert completed.returncode == 0, completed.stderr


@pytest.mark.parametrize(
    "env",
    [
        {"TAG_COMMIT": "a" * 40, "MAC_COMMIT": "b" * 40, "WINDOWS_COMMIT": "a" * 40},
        {"TAG_COMMIT": "a" * 40, "MAC_COMMIT": "a" * 40, "WINDOWS_COMMIT": "b" * 40},
        {"TAG_COMMIT": "a" * 40, "MAC_COMMIT": "b" * 40, "WINDOWS_COMMIT": "b" * 40},
    ],
)
def test_the_one_commit_check_refuses_two_commits(env):
    """"Both artifacts are built from one commit" is the claim the release page makes."""
    completed = _check_commits(**env)
    assert completed.returncode != 0, (
        f"{env} was accepted, so the release would attach two artifacts built from two "
        f"commits and say they were one.\n{completed.stdout}"
    )


@pytest.mark.parametrize(
    "missing", ["TAG_COMMIT", "MAC_COMMIT", "WINDOWS_COMMIT"]
)
def test_the_one_commit_check_refuses_an_empty_report(missing):
    """The vacuity case, and the realistic one: three empty strings are all equal.

    A build job that stopped publishing its ``commit`` output — renamed step, renamed id —
    hands this check an empty string, and an equality test alone would call that agreement
    and pass. This is the shape ``test_delivery_channels_copy`` and the vendor-drift suite
    both refuse: a guard that reads nothing must not report success.
    """
    env = {"TAG_COMMIT": "a" * 40, "MAC_COMMIT": "a" * 40, "WINDOWS_COMMIT": "a" * 40}
    env[missing] = ""
    completed = _check_commits(**env)
    assert completed.returncode != 0, (
        f"{missing} was empty and the check passed anyway, so it is comparing nothing"
    )
    assert missing in completed.stdout + completed.stderr, (
        "the refusal must name which report was missing"
    )


# ---------------------------------------------------------------------------
# The Windows job, where the traps are
# ---------------------------------------------------------------------------


def test_the_windows_build_cannot_hang_on_pause():
    """``build_installer.bat`` ends in ``pause``, on both its success and failure paths.

    On a developer's machine that is a courtesy — the window stays open long enough to read
    the output. On a runner with no keyboard it is a six-hour job that fails with no diagnosis
    at all. Redirecting stdin is the fix, and it belongs here rather than in the batch file:
    the pause is right for the humans the file was written for.
    """
    step = _step("windows", WINDOWS_BUILD_STEP)
    assert "build_installer.bat" in step["run"], step
    assert re.search(r"<\s*nul", step["run"], re.IGNORECASE), (
        f"the Windows build step does not redirect stdin, so `pause` will hang it: {step['run']!r}"
    )


def test_the_windows_job_proves_the_unversioned_compile_refuses():
    """The one assertion ``test_installer_version`` says it cannot make.

    ``installer.iss`` stops with ``#error "No AppVersion. …"`` when nothing generated
    ``version.iss``, because Inno Setup expands an undefined symbol to nothing and would
    otherwise ship an unversioned installer in silence. That guard can only read the text of
    the ``#error``; ISCC is on this runner, so the workflow compiles without the include and
    requires the refusal.
    """
    script = _step("windows", UNVERSIONED_COMPILE_STEP)["run"]
    assert "version.iss" in script, script
    assert "No AppVersion" in script, (
        "the check must assert ISCC refused for the *stated* reason; a non-zero exit alone "
        f"could be the missing python\\ tree. It reads: {script!r}"
    )
    assert _step_index("windows", UNVERSIONED_COMPILE_STEP) < _step_index(
        "windows", WINDOWS_BUILD_STEP
    ), "the unversioned compile must run before the build regenerates version.iss"


def test_no_job_hands_iscc_a_version():
    """``installer.iss`` refuses ``/DAppVersion=`` by design: the include is the authority.

    Its own comment names this workflow as the reason an earlier draft honoured such an
    override, and says why it does not: an externally-typed number is a way to stamp an
    installer with something the app does not agree with.
    """
    text = WORKFLOW.read_text(encoding="utf-8")
    assert not re.search(r"/D\s*AppVersion", text, re.IGNORECASE), (
        "the workflow passes a version to ISCC on the command line; it must run "
        "build_installer.bat, which generates version.iss from APP_VERSION"
    )


# ---------------------------------------------------------------------------
# What the release says, and that it is not published
# ---------------------------------------------------------------------------


def _release_scripts():
    return "\n".join(step.get("run", "") for step in _steps("release"))


NOTES_STEP = "Write the release notes"

#: A stand-in for the hashing tool, so these cases run the same on any machine. The workflow
#: calls ``sha256sum``, which ubuntu-latest has and macOS spells differently; what is under
#: test here is the step's own handling of what the tool returns, not the tool.
SHA256SUM_SHIM = (
    "#!/usr/bin/env python3\n"
    "import hashlib, sys\n"
    "for path in sys.argv[1:]:\n"
    "    print(hashlib.sha256(open(path, 'rb').read()).hexdigest(), '', path)\n"
)

#: The same shim, broken the way a renamed or missing tool breaks: nothing on stdout.
BROKEN_SHIM = "#!/bin/sh\nexit 1\n"

#: And broken the other way, which is the one a non-zero exit does not catch: a tool that
#: succeeds while printing something that is not a digest.
GARBAGE_SHIM = "#!/bin/sh\necho 'sha256 of something  file'\n"


def _notes_step(tmp_path, artifacts, shim=SHA256SUM_SHIM, page_text=None):
    """Run the notes step in a throwaway tree holding *artifacts*.

    The step is what turns two files into what the release page says about them, so it is run
    rather than read. The tree carries only what the script touches: the artifacts, the page's
    copy, and the first-open note that copy quotes.

    *page_text* replaces ``build/RELEASE_PAGE.md`` for the run. It exists for the cases that have
    to watch the assembly *refuse* something — a lost placeholder, a renamed one — which cannot
    be staged from the real page, since the real page is the thing that must always be right.
    ``tests/test_release_page.py`` is the caller; sharing this rather than re-implementing it is
    the point, because a second copy of the assembly would drift and go green on a page the
    workflow rejects.
    """
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    tool = bin_dir / "sha256sum"
    tool.write_text(shim, encoding="utf-8")
    tool.chmod(0o755)

    work = tmp_path / "work"
    (work / "streamlit_app" / "build").mkdir(parents=True)
    # The two files the step assembles the page out of: its copy, and the first-open note that
    # copy quotes. Real ones rather than stand-ins — what these tests check is the assembly, and
    # a seeded page would let the step keep passing while the real one had lost a placeholder.
    for name in ("OPENING_MAFIGATE.md", "RELEASE_PAGE.md"):
        (work / "streamlit_app" / "build" / name).write_text(
            (STREAMLIT_APP / "build" / name).read_text(encoding="utf-8"),
            encoding="utf-8",
        )
    if page_text is not None:
        (work / "streamlit_app" / "build" / "RELEASE_PAGE.md").write_text(
            page_text, encoding="utf-8"
        )
    for relative in artifacts:
        path = work / "artifacts" / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(relative.encode())

    step_output = work / "output.txt"
    step_output.write_text("", encoding="utf-8")
    completed = subprocess.run(
        ["bash", "-c", _step("release", NOTES_STEP)["run"]],
        cwd=work,
        env={
            "PATH": f"{bin_dir}:/usr/bin:/bin:/usr/local/bin:/opt/homebrew/bin",
            "VERSION": "1.0.0",
            "TAG": "mafigate-v1.0.0",
            "COMMIT": "a" * 40,
            "RUN_URL": "https://example.invalid/run/1",
            "GITHUB_OUTPUT": str(step_output),
        },
        capture_output=True,
        text=True,
    )
    notes = work / "release-notes.md"
    return (
        completed,
        (notes.read_text(encoding="utf-8") if notes.is_file() else ""),
        _outputs(step_output),
    )


ONE_OF_EACH = (
    "macos-dmg/MAFigate-1.0.0-macOS-universal.dmg",
    "windows-installer/MAFigate-1.0.0-Windows-Setup.exe",
)


def test_the_notes_name_both_files_their_digests_and_the_commit(tmp_path):
    completed, notes, outputs = _notes_step(tmp_path, ONE_OF_EACH)

    assert completed.returncode == 0, completed.stdout + completed.stderr
    for name in ("MAFigate-1.0.0-macOS-universal.dmg", "MAFigate-1.0.0-Windows-Setup.exe"):
        assert name in notes, f"the notes never name {name}"
    assert notes.count("a" * 40) >= 1, "the notes must name the commit both were built from"
    digests = re.findall(r"`([0-9a-f]{64})`", notes)
    assert len(digests) == 2, f"expected two sha256 digests in the notes, found {digests}"
    assert outputs.get("dmg", "").endswith(".dmg"), outputs
    assert outputs.get("exe", "").endswith(".exe"), outputs


def test_the_notes_quote_the_note_rather_than_summarising_it(tmp_path):
    """The whole file, so a Windows recipient meets it before SmartScreen has the last word."""
    _, notes, _ = _notes_step(tmp_path, ONE_OF_EACH)
    source = (STREAMLIT_APP / "build" / "OPENING_MAFIGATE.md").read_text(encoding="utf-8")
    assert source.strip() in notes, (
        "the release notes must carry build/OPENING_MAFIGATE.md verbatim; quoting part of it "
        "is how the tree ends up with two versions of the wording"
    )


@pytest.mark.parametrize(
    "artifacts",
    [
        pytest.param(ONE_OF_EACH[:1], id="no-installer"),
        pytest.param(ONE_OF_EACH[1:], id="no-dmg"),
        pytest.param(
            ONE_OF_EACH + ("macos-dmg/MAFigate-1.0.0-macOS-arm64.dmg",), id="two-dmgs"
        ),
    ],
)
def test_the_notes_step_refuses_anything_but_one_of_each(artifacts, tmp_path):
    """A release with one artifact is the asymmetric failure this workflow exists to remove.

    The *reason* is asserted, not only the refusal. Without the count check these cases still
    fail — a two-line filename or an empty one makes the hashing fall over a step later — so a
    test satisfied by any non-zero exit reported green with the check deleted, and what was
    lost was the diagnosis: "expected one of each, found …" against "sha256sum: No such file".
    Measured by deleting it and watching this file stay green.
    """
    completed, _, _ = _notes_step(tmp_path, artifacts)
    assert completed.returncode != 0, (
        f"{list(artifacts)} was accepted, so the release page would describe files it does "
        f"not have.\n{completed.stdout}"
    )
    assert "one .dmg and one .exe" in completed.stdout + completed.stderr, (
        "the refusal must say what was expected and what was found; whoever reads this log is "
        f"looking for which artifact is missing.\n{completed.stdout}{completed.stderr}"
    )


@pytest.mark.parametrize(
    "shim, description",
    [
        pytest.param(BROKEN_SHIM, "produces nothing", id="silent"),
        pytest.param(GARBAGE_SHIM, "produces something that is not a digest", id="garbage"),
    ],
)
def test_the_notes_step_refuses_a_hashing_tool_that_does_not_hash(shim, description, tmp_path):
    """The regression this was written after, found by running the step rather than reading it.

    A command substitution that fails inside a ``printf`` argument leaves an empty cell and
    exits 0 — ``set -e`` never sees it — so an earlier draft of this step wrote a release page
    offering two blank checksums: a page saying *verify this download* with nothing to verify
    it against.

    Two shims, because the two failures are different. A tool that is missing exits non-zero and
    a tool that has changed its output format does not, so only the second reaches the digest
    check; without that check the *garbage* case writes a release page with a plausible-looking
    wrong hash.

    Asserted on the notes file **not existing**, not on its contents. Its contents are ``""``
    either way — the step exits before the block that creates the file — so an assertion about
    what the string does not contain could not fail, which an earlier draft of this test did
    make and a reviewer caught.
    """
    completed, notes, outputs = _notes_step(tmp_path, ONE_OF_EACH, shim=shim)
    assert completed.returncode != 0, (
        f"a hashing tool that {description} was accepted; the notes say: {notes!r}"
    )
    assert notes == "", (
        f"the step wrote release notes anyway, which is the blank-or-wrong-checksum page this "
        f"refusal exists to prevent: {notes!r}"
    )
    assert not outputs, f"a refused step still published outputs: {outputs}"


def test_the_release_is_drafted_and_never_published_by_ci():
    """CI attaches artifacts; a human publishes.

    Which keeps two things true. #265 gets a page to author rather than one already public,
    and ``build/RELEASES.md`` — the switch the README's installer copy is held against — stays
    a record of a human act. A workflow that published would flip that switch on its own,
    from a machine, in a run nobody was reading.
    """
    scripts = _release_scripts()
    assert "--draft" in scripts, (
        "the release must be created as a draft; this workflow does not decide when a "
        "download exists"
    )


DRAFT_STEP = "Draft the release"

#: A ``gh`` that records what it was asked to do instead of doing it. Two knobs stand in for
#: the state of the repository: whether a release already exists for the tag, and whether that
#: release is still a draft.
GH_SHIM = """#!/usr/bin/env python3
import os, sys, pathlib
argv = sys.argv[1:]
log = pathlib.Path(os.environ["GH_CALLS"])
if argv[:3] == ["release", "view", os.environ.get("TAG", "")] and "--json" in argv:
    print("true" if os.environ.get("FAKE_IS_DRAFT") == "true" else "false")
    sys.exit(0)
if argv[:2] == ["release", "view"]:
    sys.exit(0 if os.environ.get("FAKE_RELEASE_EXISTS") == "true" else 1)
with log.open("a") as handle:
    handle.write(" ".join(argv) + "\\n")
"""


def _draft_step(tmp_path, prerelease="false", exists=False, is_draft=True):
    """Run the draft step against a recording ``gh``, and return what it asked ``gh`` to do."""
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    shim = bin_dir / "gh"
    shim.write_text(GH_SHIM, encoding="utf-8")
    shim.chmod(0o755)
    calls = tmp_path / "gh-calls.txt"
    calls.write_text("", encoding="utf-8")

    completed = subprocess.run(
        ["bash", "-c", _step("release", DRAFT_STEP)["run"]],
        cwd=tmp_path,
        env={
            "PATH": f"{bin_dir}:/usr/bin:/bin:/usr/local/bin:/opt/homebrew/bin",
            "GH_TOKEN": "not-a-token",
            "TAG": "mafigate-v1.0.0",
            "VERSION": "1.0.0",
            "PRERELEASE": prerelease,
            "DMG": "artifacts/macos-dmg/MAFigate-1.0.0-macOS-universal.dmg",
            "EXE": "artifacts/windows-installer/MAFigate-1.0.0-Windows-Setup.exe",
            "GH_CALLS": str(calls),
            "FAKE_RELEASE_EXISTS": "true" if exists else "false",
            "FAKE_IS_DRAFT": "true" if is_draft else "false",
        },
        capture_output=True,
        text=True,
    )
    return completed, calls.read_text(encoding="utf-8")


def test_the_release_is_created_as_a_draft_with_both_artifacts(tmp_path):
    """The step that touches the outside world, run rather than read.

    Both files by path, ``--draft``, and no ``--prerelease`` for a release tag. The earlier
    draft of this test searched the release job's concatenated scripts for ``.dmg`` and
    ``.exe``, which the notes step's own ``find`` satisfied — it stayed green with the upload
    deleted. A reviewer caught that; running the step is the fix.
    """
    completed, calls = _draft_step(tmp_path)

    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert "release create" in calls, f"nothing created a release: {calls!r}"
    assert "--draft" in calls, f"the release was not created as a draft: {calls!r}"
    assert "--prerelease" not in calls, f"a release tag was marked a pre-release: {calls!r}"
    for artifact in (".dmg", ".exe"):
        assert artifact in calls, f"the {artifact} was never attached: {calls!r}"


def test_a_rehearsal_tag_is_marked_a_prerelease(tmp_path):
    completed, calls = _draft_step(tmp_path, prerelease="true")

    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert "--prerelease" in calls, f"a rehearsal was drafted as a release: {calls!r}"
    assert "--draft" in calls, calls


def test_a_re_run_updates_an_existing_draft_rather_than_failing(tmp_path):
    """A rehearsal is worth re-running after a fix, so a tag that already has a draft is edited."""
    completed, calls = _draft_step(tmp_path, exists=True, is_draft=True)

    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert "release upload" in calls and "--clobber" in calls, calls
    assert "release edit" in calls, calls
    assert "release create" not in calls, f"it tried to create a release that exists: {calls!r}"


def test_a_published_release_is_never_touched(tmp_path):
    """The reverse of "never publishes", and the one that costs something if it is wrong.

    ``gh release edit --draft`` on a *published* release un-publishes it, so a re-run of a
    shipped tag would retract the page a clinician is downloading from and report success. A
    rebuild of a shipped tag is either a mistake or a decision to withdraw it, and the second is
    not one a re-run gets to take.
    """
    completed, calls = _draft_step(tmp_path, exists=True, is_draft=False)

    assert completed.returncode != 0, (
        f"a published release was edited; --draft would have retracted it: {calls!r}"
    )
    assert calls.strip() == "", f"it modified a published release before stopping: {calls!r}"
    assert "already published" in completed.stdout + completed.stderr, (
        "the refusal must say why it stopped and what to do instead"
    )


@pytest.mark.parametrize("missing", ["TAG", "DMG", "EXE"])
def test_the_draft_step_refuses_an_input_that_arrived_empty(missing, tmp_path):
    """Its three inputs come from earlier steps, so each is checked rather than trusted."""
    completed, calls = _draft_step(tmp_path)
    assert completed.returncode == 0, "the control case must pass, or this proves nothing"

    bin_dir = tmp_path / "bin"
    environment = {
        "PATH": f"{bin_dir}:/usr/bin:/bin:/usr/local/bin:/opt/homebrew/bin",
        "GH_TOKEN": "not-a-token",
        "TAG": "mafigate-v1.0.0",
        "VERSION": "1.0.0",
        "PRERELEASE": "false",
        "DMG": "artifacts/macos-dmg/MAFigate-1.0.0-macOS-universal.dmg",
        "EXE": "artifacts/windows-installer/MAFigate-1.0.0-Windows-Setup.exe",
        "GH_CALLS": str(tmp_path / "gh-calls.txt"),
        "FAKE_RELEASE_EXISTS": "false",
        "FAKE_IS_DRAFT": "true",
    }
    environment[missing] = ""
    completed = subprocess.run(
        ["bash", "-c", _step("release", DRAFT_STEP)["run"]],
        cwd=tmp_path,
        env=environment,
        capture_output=True,
        text=True,
    )
    assert completed.returncode != 0, f"an empty {missing} was accepted"
    assert missing in completed.stdout + completed.stderr, (
        f"the refusal must name {missing} as the empty one"
    )


def test_the_release_notes_quote_the_first_open_note_rather_than_restating_it():
    """SmartScreen fires before anything we ship can be drawn, so the release page is the
    only surface that reaches a Windows recipient in time — and ``BUILD_INSTRUCTIONS.md`` says
    it must quote ``build/OPENING_MAFIGATE.md``. Read at release time, from the one home that
    file has; ``test_unsigned_artifact_copy`` fails on a second copy of the wording.
    """
    assert "OPENING_MAFIGATE.md" in _release_scripts(), (
        "the release notes must quote build/OPENING_MAFIGATE.md — a Windows recipient meets "
        "SmartScreen on the download, before any surface of ours exists"
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
