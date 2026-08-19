"""How a MAFigate build knows what it is, and what it says when nothing told it.

Issue #263. The About dialog names the version, the channel and a short commit, so that a bug
report identifies a *build* rather than a release — necessary because there is **no update
check ever** (#229: it would be this app's first outbound call and would falsify its own
privacy claim), which leaves About as the only place the fact is available to a user at all.

Three files carry that between them, and this module holds them together:

* ``build/build_stamp.py`` writes ``config/build_stamp.py`` — gitignored, one channel, one
  short commit — and each build script calls it;
* ``config/build_identity.py`` reads it, and answers *source checkout* when it is not there;
* ``tests/test_app_identity.py`` owns the other half of the question — that the two new
  strings render in the About dialog and on no other surface.

**The default path is the one under test here.** A stamp exists only inside a built .dmg or
.exe; every clone, every ``setup.sh`` install, every CI run and every run of this suite has
none. So the assertions below are written so that they hold either way: what a stamp means is
tested through ``sys.modules`` rather than by writing a generated file into the tree, which
also means a developer who has just built an installer locally does not change their outcome.

Filed as a **guard**, and the label was argued rather than assumed: review read the reader
tests — the field-by-field fallbacks, the unfamiliar channel — as plain unit behaviour, which
would make this ``unit``. What decides it is that the module's reason to exist is a set of
copies held against their sources: the channel list, which is written twice because the writer
may not import the app; the stamp's path, derived here from the module name the app imports
rather than typed; the channel tokens in the build documentation; and the removal line in the
Makefile. The reader tests are what make those comparisons mean something, in the same way
``test_windows_installer.py`` — also a guard — asserts launcher behaviour to give its
``Source:`` comparisons their point.

Proving the derivation, not the agreement
-----------------------------------------
``test_installer_version.py``'s rule, applied to a second pair of facts: two literals that
agree today are exactly the state that drifted before. So the writer's channel list is held
equal to the app's *in both directions* rather than one importing the other, the stamp's
location is derived from the module name the app imports rather than typed twice, and the
writer is run — in a subprocess, isolated — rather than read.
"""

from __future__ import annotations

import importlib.util
import re
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import pytest

from config.build_identity import (
    INSTALLER_CHANNELS,
    SOURCE_CHANNEL,
    STAMP_MODULE,
    UNRECORDED_BUILD,
    BuildIdentity,
    build_identity,
    identify,
)

STREAMLIT_APP = Path(__file__).resolve().parents[1]
REPO_ROOT = STREAMLIT_APP.parent

WRITER = STREAMLIT_APP / "build" / "build_stamp.py"
READER = STREAMLIT_APP / "config" / "build_identity.py"
BUILD_DOCS = STREAMLIT_APP / "build" / "BUILD_INSTRUCTIONS.md"
MAKEFILE = STREAMLIT_APP / "Makefile"
DMG_SCRIPT = STREAMLIT_APP / "build" / "mac" / "build_dmg.sh"
ISS_BAT = STREAMLIT_APP / "build" / "windows" / "build_installer.bat"
ISS_SCRIPT = STREAMLIT_APP / "build" / "windows" / "installer.iss"

#: Where the stamp lives, **derived** from the module name the app imports rather than written
#: out. Typed twice, the writer could be pointed at one path while the reader looked at
#: another, and the app would report a source checkout from inside an installer — the exact
#: failure that has no symptom until a user's bug report is unattributable.
STAMP_RELATIVE = Path(*STAMP_MODULE.split(".")).with_suffix(".py")
STAMP = STREAMLIT_APP / STAMP_RELATIVE

#: A short commit as the writer produces one, with the dirty flag optional. Matched by shape:
#: pinning it to this checkout's own HEAD would make the test drift with every commit.
SHORT_COMMIT = re.compile(r"^[0-9a-f]{7,}(-dirty)?$")


def _writer():
    """``build/build_stamp.py`` as a module, loaded by path.

    Loaded rather than imported, because it deliberately is not part of any package: it runs
    on a bare build machine before this app is importable, so it lives in ``build/`` beside
    ``version.py`` and off every import path. Under a name of its own, so that nothing here
    can confuse the *writer* with the module it writes.
    """
    spec = importlib.util.spec_from_file_location("mafigate_build_stamp_writer", WRITER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _run_writer(*arguments: str) -> subprocess.CompletedProcess:
    """The writer as a build script runs it: a subprocess, and an isolated interpreter.

    ``-I`` is the load-bearing flag. It drops the user site directory and ``PYTHONPATH``, so a
    writer that had quietly grown an import of the app — or of anything installed in this
    developer's environment — fails here rather than on the first build machine that has
    neither.
    """
    return subprocess.run(
        [sys.executable, "-I", str(WRITER), *arguments],
        capture_output=True,
        text=True,
        check=False,
    )


def _stamp_values(path: Path) -> dict[str, str]:
    """The two literals in a written stamp, read the way the app reads them — by importing it.

    Executed rather than pattern-matched: the claim is that ``config/build_identity.py`` can
    read this file, and a regex over it would pass on something Python cannot import.
    """
    namespace: dict[str, object] = {}
    exec(compile(path.read_text(encoding="utf-8"), str(path), "exec"), namespace)  # noqa: S102
    return {
        "BUILD_CHANNEL": namespace.get("BUILD_CHANNEL"),
        "BUILD_COMMIT": namespace.get("BUILD_COMMIT"),
    }


def _tracked_files() -> list[str]:
    """Every file git tracks, as repo-relative paths.

    Read through git rather than :meth:`Path.rglob` for the reason four other guard modules
    give: ``.claude/worktrees`` holds full checkouts of this same tree, so a filesystem walk
    answers questions about a branch this test is not running on. Not shared with them, on the
    ground they record — hoisting it into ``conftest`` would edit five guard modules for one
    refactor, and these modules are worked by parallel branches.
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


# ---------------------------------------------------------------------------
# The writer
# ---------------------------------------------------------------------------


def test_the_writer_and_the_reader_are_both_present():
    """The anchor. Every assertion below reads one of these two files."""
    assert WRITER.is_file(), f"{WRITER} is gone — nothing writes the build stamp any more"
    assert READER.is_file(), f"{READER} is gone — nothing reads it"


def _throwaway_repo(root: Path) -> str:
    """A one-commit git repository at *root*, and the full hash of that commit.

    The writer is pointed at a repository built here rather than at this checkout, and that is
    not fussiness about isolation. Asserting against this tree's own HEAD makes the guard
    assert something about the *checkout* — it fails on a shallow clone, on an archive
    download with no ``.git``, and it says nothing about whether the writer read the right
    tree, since one wrong path would still find this repository by walking up. A repository
    with a known commit in a known place answers both.

    ``-c`` for identity and ``commit.gpgsign=false`` so a developer's global config — a signing
    key, a hooks path, a template directory — cannot decide whether this test passes.
    """
    subprocess.run(["git", "init", "--quiet", str(root)], check=True)
    (root / "a-file").write_text("first\n", encoding="utf-8")
    git = [
        "git",
        "-C",
        str(root),
        "-c",
        "user.name=MAFigate test",
        "-c",
        "user.email=test@example.invalid",
        "-c",
        "commit.gpgsign=false",
    ]
    subprocess.run([*git, "add", "a-file"], check=True)
    subprocess.run([*git, "commit", "--quiet", "--no-verify", "-m", "one commit"], check=True)
    return subprocess.run(
        ["git", "-C", str(root), "rev-parse", "HEAD"],
        capture_output=True,
        text=True,
        check=True,
    ).stdout.strip()


def test_the_writer_stamps_a_channel_and_a_commit_the_app_can_read(tmp_path):
    """End to end, in the direction a build runs it, for both installer channels.

    The commit is checked against the one this test made — a *prefix* of it, since the writer
    shortens — so the assertion is that the writer read the tree it was given rather than that
    it produced something hash-shaped.
    """
    for channel in INSTALLER_CHANNELS:
        root = tmp_path / channel
        head = _throwaway_repo(root)

        result = _run_writer("--channel", channel, "--app-root", str(root))
        assert result.returncode == 0, f"the writer failed:\n{result.stderr}"

        values = _stamp_values(root / STAMP_RELATIVE)
        assert values["BUILD_CHANNEL"] == channel
        assert SHORT_COMMIT.match(values["BUILD_COMMIT"] or ""), (
            f"the writer stamped {values['BUILD_COMMIT']!r}, which is not a short commit. "
            "A build identifier a maintainer cannot check out identifies nothing."
        )
        assert head.startswith(values["BUILD_COMMIT"]), (
            f"the writer stamped {values['BUILD_COMMIT']!r}, which is not a prefix of "
            f"{head!r}, the commit of the tree it was pointed at. It read some other tree."
        )
        assert "-dirty" not in values["BUILD_COMMIT"], (
            "a tree with nothing modified was stamped dirty, so the flag says nothing"
        )

        identity = identify(values["BUILD_CHANNEL"], values["BUILD_COMMIT"])
        assert identity == BuildIdentity(channel, values["BUILD_COMMIT"]), (
            "the app does not read back what the writer wrote, so an installer would report "
            "something other than the build it is"
        )


def test_a_build_from_a_modified_tree_is_stamped_dirty(tmp_path):
    """The flag that means *this artifact matches no commit you can check out*.

    Asserted on a **tracked** modification, and separately on an untracked file, because the
    writer counts only the first: this repository is routinely checked out beside untracked
    reference data, and a flag that fires on that would be on every build and mean nothing.
    """
    head = _throwaway_repo(tmp_path)

    (tmp_path / "untracked").write_text("ignore me\n", encoding="utf-8")
    assert _run_writer("--channel", "macos-dmg", "--app-root", str(tmp_path)).returncode == 0
    untracked_only = _stamp_values(tmp_path / STAMP_RELATIVE)["BUILD_COMMIT"]
    assert untracked_only and "-dirty" not in untracked_only, (
        f"an untracked file made the build dirty ({untracked_only!r}). Every build beside "
        "untracked data would carry the flag, which is how a flag stops meaning anything."
    )

    (tmp_path / "a-file").write_text("modified\n", encoding="utf-8")
    assert _run_writer("--channel", "macos-dmg", "--app-root", str(tmp_path)).returncode == 0
    modified = _stamp_values(tmp_path / STAMP_RELATIVE)["BUILD_COMMIT"]
    assert modified == f"{head[: len(untracked_only)]}-dirty", (
        f"a tracked file was modified and the stamp says {modified!r}. A build from an edited "
        "tree that names a bare commit sends a maintainer to code that was never built."
    )


def test_the_writer_defaults_to_the_module_the_app_imports(tmp_path):
    """Where a build script's bare invocation lands, held to where the app looks.

    Neither build script passes ``--write``: they call the writer with a channel and nothing
    else, so the default path *is* the contract between the two halves. Checked against
    :data:`STAMP_RELATIVE`, which is derived from ``STAMP_MODULE``, so renaming the module the
    app imports moves this assertion with it instead of leaving it pointing at the old name.
    """
    result = _run_writer("--channel", "macos-dmg", "--app-root", str(tmp_path))
    assert result.returncode == 0, f"the writer failed:\n{result.stderr}"

    written = tmp_path / STAMP_RELATIVE
    assert written.is_file(), (
        f"the writer put nothing at {STAMP_RELATIVE.as_posix()} under the app root it was "
        f"given. The app imports {STAMP_MODULE}, so a stamp anywhere else is a stamp nothing "
        f"reads. Tree written: {sorted(p.name for p in tmp_path.rglob('*'))}"
    )
    assert _stamp_values(written)["BUILD_CHANNEL"] == "macos-dmg"


def test_a_build_machine_without_git_still_produces_an_installer(tmp_path):
    """The asymmetry with ``build/version.py``, asserted because it is a deliberate choice.

    That tool refuses to guess a version: the version names the release, the artifact filename
    and the tag. A commit is diagnostic metadata, and someone building from a downloaded
    archive has no ``.git`` to read — refusing to build for want of a bug-report field would
    trade a working installer for a tidier one. So the build succeeds, the commit is written
    empty, a warning says so, and the *app* supplies the word for a field it was not given.

    ``tmp_path`` stands in for that machine. If it ever sat inside a repository this test would
    fail rather than skip, which is the failure worth having: a silent skip here would leave
    the whole fallback unmeasured.
    """
    result = _run_writer("--channel", "windows-installer", "--app-root", str(tmp_path))
    assert result.returncode == 0, (
        f"the writer refused to stamp a tree git cannot describe:\n{result.stderr}"
    )
    assert "warning" in result.stderr.lower(), (
        f"the writer stamped an empty commit silently. stderr was: {result.stderr!r}"
    )

    values = _stamp_values(tmp_path / STAMP_RELATIVE)
    assert values["BUILD_COMMIT"] == "", (
        f"expected an empty commit from a tree with no git in it, got "
        f"{values['BUILD_COMMIT']!r} — is {tmp_path} inside a git repository?"
    )

    identity = identify(values["BUILD_CHANNEL"], values["BUILD_COMMIT"])
    assert identity.channel == "windows-installer", (
        "an unknown commit cost the build its channel too; About would call this .exe a "
        "source checkout"
    )
    assert identity.build == UNRECORDED_BUILD


@pytest.mark.parametrize(
    "channel, why",
    [
        ("hosted", "there is no hosted channel — that deployment was removed and #229 ruled "
                   "it out of scope permanently"),
        (SOURCE_CHANNEL, "a checkout is what an *absent* stamp means; stamping one would give "
                         "that channel two representations and leave the default path — the "
                         "one every clone and every CI run takes — unexercised"),
        ("", "a build with no channel identifies nothing"),
    ],
)
def test_the_writer_refuses_a_channel_that_is_not_a_route(tmp_path, channel, why):
    """Closed rather than free-form, and refused at the writer rather than at the reader.

    The reader passes an unfamiliar channel through to the screen on purpose — a maintainer
    told ``channel hosted`` learns the build is broken, where one told ``source-checkout`` by a
    silent rewrite is told something false. That only works if the writer is the place a wrong
    channel is stopped, while a build can still be failed.
    """
    result = _run_writer("--channel", channel, "--app-root", str(tmp_path))
    assert result.returncode != 0, f"the writer accepted --channel {channel!r}: {why}"
    assert not (tmp_path / STAMP_RELATIVE).exists(), (
        f"the writer refused --channel {channel!r} and wrote a stamp anyway"
    )


def test_the_writer_and_the_app_name_the_same_installer_channels():
    """Two lists, held equal in both directions.

    The writer is stdlib-only and never imports the app — it runs before this app is
    importable, the rule ``build/version.py`` and ``vendor/_sync.py`` both follow — so the
    channel list exists twice. Both directions, because each failure is real and they are not
    the same failure: a channel the writer can stamp but the app does not know about would
    reach About unrecognised, and one the app knows about but no writer can stamp is a route
    with no way to identify itself.
    """
    assert set(_writer().INSTALLER_CHANNELS) == set(INSTALLER_CHANNELS), (
        f"build/build_stamp.py stamps {_writer().INSTALLER_CHANNELS} while the app knows "
        f"{INSTALLER_CHANNELS}. The two lists are separate on purpose and equal by test."
    )
    assert SOURCE_CHANNEL not in _writer().INSTALLER_CHANNELS, (
        "build/build_stamp.py can stamp the source channel. A checkout is what an absent "
        "stamp means, so a second representation of it would leave the default path unused."
    )


def test_the_stamp_is_generated_and_never_committed():
    """Two claims, because either alone is weak — and a third the version include does not owe.

    git must not be tracking the stamp *today*, and the ignore rule must exist so tomorrow's
    ``git add -A`` after a local build cannot commit it. The third: a tracked stamp would not
    merely be a generated file in git, it would make every clone and every CI run *claim to be
    an installer build*, which deletes the default path rather than duplicating a value.
    """
    relative = STAMP.relative_to(REPO_ROOT).as_posix()

    assert relative not in _tracked_files(), (
        f"{relative} is committed. It is written by a build, so in git it is both the copy "
        "nobody updates and a clone claiming to be an installer."
    )

    ignored = subprocess.run(
        ["git", "check-ignore", "-q", relative],
        cwd=REPO_ROOT,
        capture_output=True,
    )
    assert ignored.returncode == 0, (
        f"{relative} is not ignored, so a `git add -A` after a build commits it. "
        "streamlit_app/config/.gitignore is where that rule lives."
    )


# ---------------------------------------------------------------------------
# The reader
# ---------------------------------------------------------------------------


def test_the_app_reports_a_source_checkout_when_there_is_no_stamp():
    """The path every contributor and every CI run takes, tested explicitly.

    ``None`` in ``sys.modules`` is how the absent case is produced rather than by deleting a
    file: the stamp is absent in this checkout already, but it is *present* on a developer
    machine that has just built an installer, and a guard whose outcome depends on whether
    someone ran ``make build-mac`` this morning is not a guard.
    """
    with patch.dict(sys.modules, {STAMP_MODULE: None}):
        identity = build_identity()

    assert identity == BuildIdentity(SOURCE_CHANNEL, UNRECORDED_BUILD)
    assert str(identity) == f"channel {SOURCE_CHANNEL} · build {UNRECORDED_BUILD}", (
        "the line About prints has changed shape; tests/test_app_identity.py sweeps the "
        "surfaces for these two values and this is where their spelling is settled"
    )


def test_the_app_reads_the_stamp_when_there_is_one():
    """A stamp stood up in ``sys.modules``, for the same reason as above."""
    stamp = SimpleNamespace(BUILD_CHANNEL="macos-dmg", BUILD_COMMIT="a1b2c3d")

    with patch.dict(sys.modules, {STAMP_MODULE: stamp}):
        identity = build_identity()

    assert identity == BuildIdentity("macos-dmg", "a1b2c3d")
    assert str(identity) == "channel macos-dmg · build a1b2c3d"


@pytest.mark.parametrize(
    "stamp, why",
    [
        (SimpleNamespace(), "a stamp written by an older writer, carrying neither key"),
        (
            SimpleNamespace(BUILD_CHANNEL="macos-dmg"),
            "a stamp carrying a channel and no commit — a key the writer stopped writing",
        ),
        (
            SimpleNamespace(BUILD_CHANNEL="  ", BUILD_COMMIT="\n"),
            "whitespace, which is what a template with an unfilled hole leaves behind",
        ),
        (
            SimpleNamespace(BUILD_CHANNEL=None, BUILD_COMMIT=("a1b2c3d",)),
            "values that are not strings at all",
        ),
    ],
)
def test_a_damaged_stamp_falls_back_field_by_field(stamp, why):
    """Each field decided on its own, and nothing unstated reaching the screen.

    A stamp is a generated file, so the failures worth surviving are the generated ones. Field
    by field rather than all-or-nothing, because the case that actually happens — a build
    machine with no git — is a real channel beside an empty commit, and discarding the channel
    over the commit would tell a bug reporter their installer is a clone.
    """
    with patch.dict(sys.modules, {STAMP_MODULE: stamp}):
        identity = build_identity()

    assert identity.channel in {SOURCE_CHANNEL, *INSTALLER_CHANNELS}, (
        f"{why}: channel read as {identity.channel!r}"
    )
    assert identity.build, f"{why}: the build identifier read as empty"
    assert "None" not in str(identity), f"{why}: {str(identity)!r} reached the About line"


def test_an_unfamiliar_channel_is_passed_through_rather_than_hidden():
    """A hand-edited stamp says something visibly odd instead of something quietly false.

    The alternative — mapping anything unrecognised to :data:`SOURCE_CHANNEL` — would tell a
    maintainer that a shipped installer is a clone, which is a wrong answer where this is
    merely a strange one. The writer is where a wrong channel is refused, and it is refused
    while the build can still be stopped.
    """
    stamp = SimpleNamespace(BUILD_CHANNEL="hosted", BUILD_COMMIT="a1b2c3d")

    with patch.dict(sys.modules, {STAMP_MODULE: stamp}):
        assert build_identity().channel == "hosted"


def test_only_the_reader_imports_the_stamp():
    """One reader, so the source-checkout default cannot be bypassed.

    A page that imported ``config.build_stamp`` directly would raise ``ImportError`` in every
    clone — the majority case — so the fallback is not a nicety that a second reader may skip.
    Asserted as an exact set rather than a maximum, so a third file naming the stamp fails
    here and someone has to decide, which is the point at which the decision is cheap.
    """
    named = {
        relative
        for relative in _tracked_files()
        if relative.endswith(".py")
        and not relative.startswith("streamlit_app/tests/")
        and STAMP_MODULE.split(".")[-1] in (REPO_ROOT / relative).read_text(
            encoding="utf-8", errors="replace"
        )
    }

    assert named == {
        "streamlit_app/build/build_stamp.py",
        "streamlit_app/config/build_identity.py",
    }, (
        f"the build stamp is named by {sorted(named)}. Exactly two files may: the one that "
        "writes it, and the one reader that knows what its absence means."
    )


# ---------------------------------------------------------------------------
# The two build scripts
# ---------------------------------------------------------------------------


def _line_index(text: str, needle: str, description: str) -> int:
    """The first line that *runs* *needle*, or an assertion naming what was being looked for.

    Comments are skipped, and so is ``echo`` — which is not fussiness. Both scripts print the
    writer's own name in the message they fail with, so counting an ``echo`` as the step would
    have let a build that deleted the invocation and kept the error message pass the ordering
    assertion below. Found while proving that assertion fails when it should.
    """
    for index, line in enumerate(text.splitlines()):
        stripped = line.strip().lower()
        if stripped.startswith(("#", "rem ", "echo ", 'echo "')):
            continue
        if needle in line:
            return index
    raise AssertionError(
        f"no line in this script {description} — it contains no {needle!r} outside its "
        "comments, so an assertion about where that step happens is about nothing"
    )


@pytest.mark.parametrize(
    "script, channel, before, what_follows",
    [
        (
            DMG_SCRIPT,
            "macos-dmg",
            "for item in MAFigate.py",
            "copies the app source into the bundle",
        ),
        (
            ISS_BAT,
            "windows-installer",
            "installer.iss",
            "hands the script to the Inno Setup compiler",
        ),
    ],
    ids=["macOS .dmg", "Windows .exe"],
)
def test_each_build_script_stamps_its_own_channel_before_it_packages(
    script, channel, before, what_follows
):
    """Both builds, and the ordering that makes the stamp reach the artifact.

    Ordering, not merely presence: the .dmg copies ``config`` into the bundle and the .exe
    compiles ``config\\*`` into the installer, so a stamp written *after* that step is a stamp
    in the developer's tree and nowhere else — an installer reporting itself as a source
    checkout, with the build having printed nothing but success.

    Read from the scripts because neither can be run here: one needs macOS, the other needs
    Windows and Inno Setup. The first real build is #264/#265's, and this is what stands in
    for it until then.
    """
    text = script.read_text(encoding="utf-8")

    stamp_at = _line_index(text, "build_stamp.py", "runs the build-stamp writer")
    package_at = _line_index(text, before, what_follows)

    assert channel in text, (
        f"{script.name} runs the stamp writer without naming --channel {channel}. Each build "
        "declares its own route; a build that let the writer guess would tell a bug reporter "
        "something unverifiable."
    )
    assert stamp_at < package_at, (
        f"{script.name} writes the build stamp at line {stamp_at + 1}, after the step that "
        f"{what_follows} at line {package_at + 1}. The stamp would stay in the build tree and "
        "the artifact would call itself a source checkout."
    )


def test_the_build_instructions_name_the_real_channels():
    """The documented commands are copies, so they are held to the channel list.

    ``BUILD_INSTRUCTIONS.md`` prints ``build/build_stamp.py --channel macos-dmg`` for a builder
    to read, which is a second place the channel tokens are written. The version's guard sweeps
    that same prose for exactly this reason — and *banned* is not available here, because the
    document has to name the commands. So they are pinned instead: rename a channel and the
    document goes red until someone rewrites it.

    Both directions. A token the document names that is not a channel is a command that fails
    when typed; a channel the document never names is a route nobody is told how to build.
    """
    text = BUILD_DOCS.read_text(encoding="utf-8")

    documented = set(re.findall(r"--channel\s+([A-Za-z0-9-]+)", text))
    assert documented, (
        f"{BUILD_DOCS.name} shows no `--channel` invocation at all, so this guard is reading "
        "nothing — and a builder is not told how the stamp is written."
    )
    assert documented == set(INSTALLER_CHANNELS), (
        f"{BUILD_DOCS.name} documents `--channel` values {sorted(documented)} while the "
        f"writer accepts {sorted(INSTALLER_CHANNELS)}."
    )


def test_make_clean_removes_a_stamp_a_local_build_left_behind():
    """The remedy for the one thing a build leaves in the developer's tree.

    A build writes the stamp into the working tree, and nothing removes it afterwards — so
    running the app from source after ``make build-mac`` reports ``channel macos-dmg``, which
    is false of what is running. Cleanup *inside* the build scripts would be worse than this:
    ``set -e`` aborts before any trailing cleanup, so the failure path would leave the stamp
    behind anyway while the tidy path pretended otherwise.

    So the remedy is ``make clean``, and this holds the recipe to it — by the path derived from
    the module the app imports, not a typed one. Asserted as text rather than by running
    ``make clean``, which would delete this run's own ``.pytest_cache`` from under it.
    """
    recipe = MAKEFILE.read_text(encoding="utf-8")
    wanted = STAMP_RELATIVE.as_posix()

    removals = [
        line.strip()
        for line in recipe.splitlines()
        if wanted in line and not line.strip().startswith("@#")
    ]
    assert removals, (
        f"nothing in the Makefile removes {wanted}. A stamp left by a local build outlives it, "
        "and the app then reports a source run as an installer build."
    )
    assert any(line.startswith(("rm ", "rm -f", "-rm")) for line in removals), (
        f"the Makefile names {wanted} but does not remove it: {removals}"
    )


def test_the_stamp_travels_inside_both_installers():
    """The stamp reaches a user because ``config`` is shipped whole, and that is asserted.

    Neither installer names the stamp file. That is the design — a build product listed one by
    one is the line someone forgets — but it means both installers ship it only for as long as
    they keep shipping the whole package, so the two lines that do it are held here.
    """
    dmg = DMG_SCRIPT.read_text(encoding="utf-8")
    copy_lists = [line for line in dmg.splitlines() if line.startswith("for item in ")]
    assert copy_lists, "build_dmg.sh no longer has a copy list this test can read"
    assert any(
        re.search(r"(?:^|\s)config(?:\s|$)", line) for line in copy_lists
    ), (
        f"build_dmg.sh's copy list no longer copies `config` whole: {copy_lists}. The build "
        "stamp is written inside it, so the .dmg would ship an app with no stamp in it."
    )

    iss = ISS_SCRIPT.read_text(encoding="utf-8")
    shipped = [line for line in iss.splitlines() if r"config\*" in line]
    assert shipped, (
        r"installer.iss no longer ships config\* — the build stamp is written there and would "
        "not reach the installed app"
    )
    assert all("recursesubdirs" in line for line in shipped), (
        f"installer.iss ships config without recursesubdirs: {shipped}"
    )
