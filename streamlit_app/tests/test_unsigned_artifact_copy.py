"""What we tell someone whose computer has just refused to open MAFigate — said once.

Neither installer is signed with a paid certificate, and neither will be, so macOS
Gatekeeper and Windows SmartScreen **will** fire on every recipient's first open. That
makes the wording of the reassurance a shipped feature rather than a footnote, and issue
#262 settled two things about it.

**One home.** Before this ticket the tree told recipients two different `xattr`
incantations — ``-cr`` in ``build/BUILD_INSTRUCTIONS.md``'s troubleshooting table, ``-dr
com.apple.quarantine`` in the DMG build script's closing message — each denying the other
was the answer. The fix is not to correct both copies but to have one, so the guards here
pin the canonical wording to a **single file** and let every other surface link to it. A
maintainer who pastes the instruction onto a second surface gets a red test rather than a
second copy to keep true; that is the whole mechanism, and it is aimed squarely at the
failure this repo has already had.

**GUI-first, and the ordering is load-bearing.** These recipients are clinicians and
molecular tumour board members who will not open a terminal, so a command line is the
wrong first instruction in any spelling. :func:`test_the_gui_route_comes_first` reads the
note's own structure — where each phrase falls in the file, and which heading the
command sits under — rather than trusting that it is *mentioned*, because a fallback
mentioned above the primary route is not a fallback.

**The timing that decides where the text can go.** SmartScreen fires on the downloaded
executable before any surface of ours can be drawn, so no file we ship reaches a Windows
user in time — only the release page does, which is why the note is written to be quotable
there. Gatekeeper's alert fires when the app is *opened*, after the DMG window is already
on screen, so macOS gets a second chance and takes it:
:func:`test_the_dmg_ships_the_note_beside_the_app` holds the build script to putting the
note in the staging directory it assembles, with a position of its own in the mounted
window so it does not land under the app icon.

**And the publisher URL, which has now been found three times.** ``installer.iss``'s
``AppPublisherURL`` named a scaffolding placeholder where the account belongs — the token
:data:`PLACEHOLDER_ORG` assembles: broken on the day it was written,
rendered in Windows' Add/Remove Programs, and permanently invisible to
``test_public_repo_name.py``'s sweep because it names *neither* repository correctly.
Issue #234 considered and declined the widened rule — every github.com URL naming a
variantalker repo must name the public one — because that needs an exemption list for
genuine third-party links. The cheap covering rule this ticket found instead needs no
exemptions at all: **the placeholder token itself may not appear in the tree.** It is a
placeholder, so there is no legitimate use to exempt.

**Why the literals below are assembled from halves.** Every claim here is *"this phrase
lives in exactly one file"* or *"this phrase lives in none"*, swept over ``git
ls-files``. Spelled out, this module would be the second home of every phrase it pins and
could never be green. The precedent is ``test_public_repo_name.py``, which composes for
the same reason and states the same cost: a human grepping for one of these phrases will
not find the guard that keeps it single. That is what this paragraph is for. Composing
was chosen over exempting this file from its own sweep, deliberately — the sibling guard's
note that *a guard which can be argued into narrowing is a guard that will be* applies to
an exemption list of one as much as to a long one.
"""

from __future__ import annotations

import base64
import re
import subprocess
from pathlib import Path

import pytest

# The public repository's URL, imported rather than spelled a second time. It is the one
# literal here nobody had to compose, so a second copy would be free to drift from the
# sibling guard that owns which repository this tree names out loud — and the publisher URL
# Windows renders should be *that* repository by construction, not by coincidence.
from tests.test_public_repo_name import PUBLIC_REPO

STREAMLIT_APP = Path(__file__).resolve().parents[1]
REPO_ROOT = STREAMLIT_APP.parent

BUILD = STREAMLIT_APP / "build"

#: The single home. Named as a path rather than discovered, because "which file is
#: canonical" is the decision this ticket made and a guard that searched for it would
#: accept whichever file happened to hold the phrases.
CANONICAL_NOTE = BUILD / "OPENING_MAFIGATE.md"

#: The build script that assembles the mounted DMG window, and the instructions that used
#: to restate the incantation.
DMG_BUILD_SCRIPT = BUILD / "mac" / "build_dmg.sh"
BUILD_INSTRUCTIONS = BUILD / "BUILD_INSTRUCTIONS.md"
WINDOWS_INSTALLER = BUILD / "windows" / "installer.iss"

#: What the note is called once it is sitting in the mounted DMG beside the app. A ``.txt``
#: rather than the source's ``.md`` on purpose: a Mac with no developer tooling has no
#: default application for ``.md``, and a note that will not open when double-clicked is
#: not a note. See the build script's own comment.
STAGING_NOTE_NAME = "Opening MAFigate.txt"

#: The phrases that may live in exactly one file, assembled from halves — see the module
#: docstring. Keyed by what each one is, because the failure message has to tell a
#: maintainer which of the three they duplicated.
CANONICAL_PHRASES = {
    "the macOS button to click": "Open " + "Anyway",
    "the Windows button to click": "Run " + "anyway",
    "the Terminal fallback command": "xattr -dr com.apple." + "quarantine",
    # The other half of the ticket's sentence, and the half a pin over buttons alone
    # misses: *"one canonical wording ... covering what each alert says **and** what to
    # click."* A tracked file holding a second copy of every alert quotation while the
    # button names stayed single would satisfy a buttons-only guard completely — and the
    # alert text is the more likely thing to be copied, because it is what a release page
    # or a support email quotes to be recognised.
    #
    # Each is a fragment rather than the whole sentence, because the sentences wrap: the
    # confirmation's full wording is one line in the note precisely so that it can be
    # pinned, and it went unpinnable once already by sitting across a line break.
    "the macOS alert's own words": "free of " + "malware",
    "the two buttons that alert offers": "Move to " + "Trash",
    "the confirmation sheet's sentence": "always allow it to run " + "on this Mac",
    "the SmartScreen headline": "Windows protected " + "your PC",
}

#: Which entries above are *what the alert says* rather than *what to click*. Named so the
#: presence check below reads one list instead of keeping a second copy of the quotations.
ALERT_KEYS = (
    "the macOS alert's own words",
    "the two buttons that alert offers",
    "the confirmation sheet's sentence",
    "the SmartScreen headline",
)

#: The Control-click route, named so the ordering guard can hold the deviation in place.
#: Not pinned to one home — ``BUILD_INSTRUCTIONS.md`` names it too, when citing Apple for
#: why it is no longer the headline, and that is a statement about the OS rather than an
#: instruction to a recipient.
SUPERSEDED_MAC_ROUTE = "Control-click"

#: The spelling this ticket retired. It clears *every* extended attribute rather than the
#: one quarantine flag, so it is not a rougher phrasing of the same instruction but a
#: different instruction — and it was the one the troubleshooting table gave. Zero
#: occurrences, not one: the point is that neither is left saying something the other
#: denies. Assembled from halves for the same reason as everything else here.
SUPERSEDED_XATTR = "xattr " + "-cr"

#: The placeholder that stood where an account name belongs. Assembled for the same reason
#: as everything else here.
PLACEHOLDER_ORG = "your" + "-org"

#: An independent spelling of each literal swept to **zero**, base64-encoded, checked
#: against the composition above by :func:`test_the_literals_swept_to_zero_are_the_right
#: _strings`.
#:
#: This exists because the obvious vacuity check does not work here, and it is worth being
#: precise about why. The three phrases swept to *one* witness themselves: mistype one of
#: those compositions and its count drops to zero, which is not one, so the guard goes red.
#: The two swept to *zero* have no such luck — a mistyped composition finds nothing, zero is
#: the answer the guard wants, and it passes for ever. Nor can a planted-hit check catch it:
#: planting the composed value and then searching for that same value succeeds no matter
#: what the value is. It proves the matcher reads files; it cannot prove the string is still
#: the forbidden one.
#:
#: So the witness has to be a spelling this module did not build from the same halves.
#: Encoded rather than written out for the reason the whole module is: in plain text these
#: would be the tracked occurrences their own guards forbid. Decode them to read them.
ZERO_SWEPT_ENCODED = {
    "the retired xattr spelling": "eGF0dHIgLWNy",
    "the placeholder account name": "eW91ci1vcmc=",
}

#: Every literal any sweep in this module looks for, keyed by what it is.
#:
#: The planted-hit check below runs over *this*, not over :data:`CANONICAL_PHRASES` alone,
#: and the two extra entries are the reason. :data:`SUPERSEDED_XATTR` and
#: :data:`PLACEHOLDER_ORG` are claims of **zero**, and they are assembled from halves — so
#: a one-character slip in either composition (a space, an underscore for the hyphen) makes
#: its guard permanently green while reading exactly like a clean tree. A matcher that has
#: stopped seeing its literal reports zero occurrences, which is the answer the guard wants
#: to hear. The convention is ``test_public_export.py``'s, and ``tests/README.md`` states
#: it in that module's row: *a planted hit of every rule, since a matcher that had stopped
#: seeing one rule would also report zero.*
SWEPT_LITERALS = {
    **CANONICAL_PHRASES,
    "the retired xattr spelling": SUPERSEDED_XATTR,
    "the placeholder account name": PLACEHOLDER_ORG,
}

#: The one excluded tree, as a ``git ls-files`` path prefix, and the only exclusion any
#: sweep here takes. That tree is the measurement record of the effort itself —
#: reports written *to* the private tracker, quoting a ticket's acceptance criteria as the
#: subject of the work rather than as an instruction to a user. Nothing in it is packaged
#: or rendered. It is also stripped by the public export, so in a public clone this prefix
#: matches nothing and the sweep is simply total.
EXCLUDED_PREFIX = "docs/wayfinder/"

#: One screen, asserted as a number so that "one screen" cannot drift into four. Generous
#: enough for the note to name both alerts and both routes, tight enough that a section on
#: something else does not fit.
MAX_NOTE_LINES = 80


def _tracked_files() -> list[str]:
    """Every path git tracks, relative to the repository root.

    Read through git rather than :meth:`Path.rglob` for a reason that has bitten this repo
    before: ``.claude/worktrees/`` holds whole checkouts of this same tree, so a
    filesystem walk from the root finds another session's files and answers questions
    about a branch this test is not running on.

    A failure to run git is reported as a failure of the guard rather than a skip. Every
    claim below is a count, and a skip in a CI log is indistinguishable from a pass.
    """
    try:
        completed = subprocess.run(
            ["git", "-C", str(REPO_ROOT), "ls-files", "-z"],
            capture_output=True,
            check=False,
        )
    except OSError as exc:  # pragma: no cover - depends on the machine, not the code
        pytest.fail(
            f"could not run git to list the tracked tree ({exc}). These guards count "
            "occurrences across what the repository ships, so they need the repository; "
            "they do not excuse themselves when they cannot look."
        )

    assert completed.returncode == 0, (
        f"`git ls-files` failed in {REPO_ROOT} (exit {completed.returncode}): "
        f"{completed.stderr.decode('utf-8', 'replace').strip()}"
    )
    out = completed.stdout.decode("utf-8", "replace")
    return [path for path in out.split("\0") if path]


def _files_containing(phrase: str, swept: list[str]) -> list[str]:
    """The swept paths whose text contains ``phrase``, in sorted order.

    Read as bytes and decoded permissively: the tracked tree holds MAFs, a spreadsheet, a
    jar and several PNGs, and a decode error must not be able to hide a hit. Binary files
    are swept for the same reason — every phrase here is ASCII and would survive into one.
    """
    found = []
    for path in swept:
        full = REPO_ROOT / path
        try:
            text = full.read_bytes().decode("utf-8", "replace")
        except OSError as exc:  # pragma: no cover - depends on the checkout
            pytest.fail(
                f"cannot read the tracked file {path} ({exc}), so this guard cannot say "
                "whether it holds a second copy. An unreadable file is a hole in a count, "
                "not a file with nothing in it."
            )
        if phrase in text:
            found.append(path)
    return sorted(found)


def _headings(text: str) -> list[tuple[int, str]]:
    """The note's headings as ``(offset, text)``, read as **setext** headings.

    Setext — a line of prose underlined with ``===`` or ``---`` — rather than ATX's
    leading ``#``, because this file is read twice: rendered as Markdown in the repository,
    and as literal text by a recipient who double-clicks it in the mounted DMG. Setext is
    the one heading syntax that is a heading in the first reading and still looks like one
    in the second. See :func:`test_the_note_is_readable_as_plain_text`.
    """
    lines = text.splitlines()
    found = []
    offset = 0
    for index, line in enumerate(lines):
        following = lines[index + 1] if index + 1 < len(lines) else ""
        underlined = len(following) >= 3 and set(following) in ({"="}, {"-"})
        if underlined and line.strip():
            found.append((offset, line.strip()))
        offset += len(line) + 1
    return found


@pytest.fixture(scope="module")
def swept() -> list[str]:
    """The paths these guards sweep: everything tracked, outside the one excluded tree."""
    return [path for path in _tracked_files() if not path.startswith(EXCLUDED_PREFIX)]


@pytest.fixture(scope="module")
def note_text() -> str:
    assert CANONICAL_NOTE.is_file(), (
        f"the canonical note is missing from {CANONICAL_NOTE.relative_to(REPO_ROOT)}. "
        "It is the single home issue #262 gave this wording and the file the DMG ships "
        "beside the app; without it every surface below has nowhere to point."
    )
    return CANONICAL_NOTE.read_text(encoding="utf-8")


# --------------------------------------------------------------------------------------
# Vacuity first, because every claim below is a count.
# --------------------------------------------------------------------------------------


def test_the_sweep_reads_a_real_tree(swept):
    """A truncated or detached ``git ls-files`` satisfies a count of one as easily as a
    count of a thousand, and satisfies a count of zero completely. Two anchors: the list
    is large, and it names files these guards are actually about — named rather than
    counted alone, so a truncated list cannot pass under a threshold.
    """
    assert len(swept) > 100, (
        f"the tracked-file sweep found only {len(swept)} paths outside "
        f"{EXCLUDED_PREFIX!r}, which is too few for this repository. A broken sweep "
        "reports every phrase as living in zero files."
    )
    for expected in ("streamlit_app/MAFigate.py", "streamlit_app/build/mac/build_dmg.sh"):
        assert expected in swept, (
            f"the sweep does not include {expected}, so whatever it is reading is not "
            f"this repository. {len(swept)} paths were listed."
        )


def test_the_literals_swept_to_zero_are_the_right_strings():
    """That the two claims of zero are still claims about the strings they were written for.

    The other half of vacuity, and the half a planted-hit check cannot reach — see
    :data:`ZERO_SWEPT_ENCODED`. Without this, a slip of one character in either composition
    leaves a guard that is green because it is looking for a string that occurs nowhere,
    which is indistinguishable in a CI log from a guard that is green because the tree is
    clean.
    """
    for what, encoded in ZERO_SWEPT_ENCODED.items():
        expected = base64.b64decode(encoded).decode("ascii")
        assert SWEPT_LITERALS[what] == expected, (
            f"{what} is composed as {SWEPT_LITERALS[what]!r}, but the independent spelling "
            f"says it should be {expected!r}. Its guard sweeps for a string that is not the "
            "one being forbidden, so it reports a clean tree whatever the tree holds."
        )

    assert set(ZERO_SWEPT_ENCODED) < set(SWEPT_LITERALS), (
        "the witnesses name literals this module does not sweep, so at least one of them "
        f"guards nothing: {sorted(set(ZERO_SWEPT_ENCODED) - set(SWEPT_LITERALS))}"
    )


def test_the_matcher_can_actually_find_a_planted_copy(tmp_path, monkeypatch):
    """The other half: that :func:`_files_containing` reports a hit it is shown.

    Over **every** swept literal, including the two whose guards claim a zero — see
    :data:`SWEPT_LITERALS` for why those are the ones that would go quiet unnoticed.
    Asserted against fabricated files rather than against the real tree, so the machinery
    is proved by something this test controls and stays proved if the note is rewritten.
    """
    monkeypatch.setattr("tests.test_unsigned_artifact_copy.REPO_ROOT", tmp_path)

    for what, phrase in SWEPT_LITERALS.items():
        (tmp_path / "planted.md").write_text(f"a first line\nclick {phrase} to continue\n")
        (tmp_path / "clean.md").write_text("a note that says none of it\n")
        assert _files_containing(phrase, ["planted.md", "clean.md"]) == ["planted.md"], (
            f"the matcher does not find {phrase!r} ({what}) in a file that contains it. "
            "Its count over the real tree therefore means nothing — and for the two "
            "literals swept to zero, that count is what the guard wants to hear."
        )


# --------------------------------------------------------------------------------------
# One canonical wording, in one place.
# --------------------------------------------------------------------------------------


def test_each_canonical_phrase_lives_in_exactly_one_file(swept):
    """The pin itself, and the reason this module exists.

    Reported per phrase with every offending path, because the fix differs by surface: a
    build script's closing message becomes a pointer to the note, a troubleshooting table
    row becomes a link to it.
    """
    canonical = str(CANONICAL_NOTE.relative_to(REPO_ROOT))
    wrong = {}
    for what, phrase in CANONICAL_PHRASES.items():
        homes = _files_containing(phrase, swept)
        if homes != [canonical]:
            wrong[what] = (phrase, homes)

    report = "\n".join(
        f"  {what} ({phrase!r}) lives in {homes or 'no tracked file'}"
        for what, (phrase, homes) in sorted(wrong.items())
    )
    assert not wrong, (
        f"{len(wrong)} of {len(CANONICAL_PHRASES)} canonical phrases are not in exactly "
        f"one file (issue #262):\n{report}\n"
        f"The one home is {canonical}. Another surface that needs to say this links to "
        "that file instead of restating it — the two contradictory `xattr` lines this "
        "ticket removed are what restating costs."
    )


def test_the_superseded_xattr_spelling_is_gone(swept):
    """Neither incantation may be left saying something the other denies.

    A zero rather than a one: ``-cr`` clears every extended attribute instead of the
    single quarantine flag, so it is not a second-best phrasing of the same instruction —
    it is a different instruction, and it was the one the troubleshooting table gave.
    """
    homes = _files_containing(SUPERSEDED_XATTR, swept)
    assert homes == [], (
        f"the retired `{SUPERSEDED_XATTR}` spelling still appears in {homes}. Issue #262 "
        "reduced the two contradictory incantations to one — the targeted one, which "
        f"clears only the quarantine flag and lives only in "
        f"{CANONICAL_NOTE.relative_to(REPO_ROOT)}."
    )


def test_the_note_fits_one_screen(note_text):
    """It is read at the moment something has already gone wrong, in a TextEdit window on
    a Mac that has just refused to open an app. A second screen is a second chance to
    stop reading.
    """
    lines = note_text.splitlines()
    assert len(lines) <= MAX_NOTE_LINES, (
        f"{CANONICAL_NOTE.relative_to(REPO_ROOT)} is {len(lines)} lines, over the "
        f"{MAX_NOTE_LINES}-line budget. It ships in the mounted DMG window as a "
        "one-screen note; anything that does not fit belongs on a surface the reader "
        "chose to open."
    )
    overlong = [(number, line) for number, line in enumerate(lines, 1) if len(line) > 100]
    assert not overlong, (
        "these lines are over 100 characters, so they wrap unreadably in the plain-text "
        "copy the DMG ships:\n"
        + "\n".join(f"  line {number}: {len(line)} chars" for number, line in overlong)
    )


def test_the_note_is_readable_as_plain_text(note_text):
    """Because the copy the recipient reads is this file, renamed.

    ``build_dmg.sh`` stages it as ``.txt`` for a good reason — a Mac with no developer
    tooling has no default application for ``.md``, and a note that will not open when
    double-clicked is not a note. But that rename does nothing to the *content*, so any
    Markdown syntax in here is delivered raw to the one surface that reaches a macOS user
    in time: a frightened reader meeting ``## macOS`` and ``> Apple could not verify`` in
    TextEdit. The script's comment argued for the rename precisely because of the audience;
    shipping markup through it would contradict the argument.

    So the note is written in the syntax that survives both readings: setext headings,
    indented blocks for quoted material and for the command, and no inline emphasis at all.
    Asserted rather than remembered, since the natural thing to type is ``#``.
    """
    offenders = []
    for number, line in enumerate(note_text.splitlines(), start=1):
        if line.startswith("#"):
            offenders.append((number, "an ATX heading marker", line))
        elif line.startswith(">"):
            offenders.append((number, "a blockquote marker", line))
        elif "**" in line or "`" in line:
            offenders.append((number, "inline emphasis or code markup", line))

    report = "\n".join(f"  line {n}: {why} — {line!r}" for n, why, line in offenders)
    assert not offenders, (
        f"{CANONICAL_NOTE.name} carries Markdown syntax that a recipient reads literally, "
        f"because the DMG ships this file as {STAGING_NOTE_NAME!r} without converting it:\n"
        f"{report}\n"
        "Use a setext heading (underline the line with === or ---), indent quoted material "
        "and commands by four spaces, and leave emphasis out."
    )

    assert _headings(note_text), (
        f"{CANONICAL_NOTE.name} has no setext headings, so either it has no sections or "
        "they are marked in a syntax this guard reads as prose — and one of those means "
        "the fallback below is no longer under a heading anyone can see."
    )


def test_the_note_says_what_each_alert_says(note_text):
    """Both alerts, not just both fixes.

    The alert text is what a frightened reader matches against their screen before they
    will follow any instruction, and it is the half a summary drops first. Held to Apple's
    own strings, read from ``CoreServicesUIAgent``'s ``Quarantine.loctable`` on macOS 26 —
    ``Q_DETAIL_CASPIAN_UNVERIFIED`` for the alert, ``Q_BUTTON_MOVE_TO_TRASH`` for the pair
    of buttons it offers, and ``Q_DETAIL_CASPIAN_ALLOW`` for the sheet raised by the button
    in Settings — and to SmartScreen's headline. That sheet is pinned because it is where
    the route is easiest to abandon: its button is ``Q_BUTTON_OPEN``, plain *Open*, not a
    repeat of the one just clicked in Settings, so a reader who was told to expect a repeat
    starts hunting for a button that is not on screen. Apple's own user guide
    (``support.apple.com/guide/mac-help/mh40616``) documents what follows as entering your
    login password — not Touch ID, which an earlier draft of the note claimed.
    """
    quotations = {what: CANONICAL_PHRASES[what] for what in ALERT_KEYS}
    # The one alert word not pinned to a single home: two generic words any surface may
    # legitimately use, so it is asserted present rather than unique.
    quotations["the SmartScreen link that reveals the button"] = "More info"

    for whose, fragment in quotations.items():
        assert fragment in note_text, (
            f"the note does not carry {fragment!r} — {whose}. A reader compares the words "
            "on their screen to the words in the note before trusting either."
        )


def test_the_gui_route_comes_first(note_text):
    """GUI-first, read off the note's structure rather than its intent.

    Two separate claims, because they fail separately. The command must come *after* both
    GUI instructions — a fallback printed above the primary route is the headline
    whatever it is called — and it must sit under a heading that says it needs the
    Terminal, so a reader who will not open one knows to stop before the command rather
    than after it.
    """
    command = CANONICAL_PHRASES["the Terminal fallback command"]
    where_command = note_text.index(command)
    for what in ("the macOS button to click", "the Windows button to click"):
        phrase = CANONICAL_PHRASES[what]
        assert note_text.index(phrase) < where_command, (
            f"{phrase!r} ({what}) appears after the Terminal command in "
            f"{CANONICAL_NOTE.name}. For this audience the command line is the wrong "
            "first instruction in any spelling: it needs a Terminal they will not open."
        )

    headings = _headings(note_text)
    assert headings, f"{CANONICAL_NOTE.name} has no headings, so nothing marks the fallback"
    above = [text for start, text in headings if start < where_command]
    assert above, f"the Terminal command in {CANONICAL_NOTE.name} sits under no heading"
    assert "Terminal" in above[-1], (
        f"the Terminal command sits under the heading {above[-1]!r}, which does not warn "
        "that it needs a Terminal. The `xattr` line survives only as a clearly-marked "
        "fallback for technical users, never as part of the route everyone else takes."
    )


def test_the_settings_route_leads_and_the_control_click_route_follows(note_text):
    """The deviation from the ticket's named macOS route, held in place.

    Issue #262's second criterion names *right-click → Open → Open* as the primary macOS
    instruction. **macOS 15 removed that override** for software that is not correctly
    signed or notarized — which is this DMG, ad-hoc signed and never notarized — so
    following the criterion literally would ship GUI-first advice that fails on every
    current Mac. Apple's own announcement is the source: users "will no longer be able to
    Control-click to override Gatekeeper" and "will need to visit System Settings >
    Privacy & Security" (*Updates to runtime protection in macOS Sequoia*,
    ``developer.apple.com/news/?id=saqachfa``, 6 August 2024).

    So the note leads with Settings and keeps Control-click as the macOS 14-and-older
    case, and this guard is why that cannot quietly revert. Without it the note could be
    reordered to lead with the shortcut and every other test here would stay green: the
    ordering check above only asks that both GUI routes precede the Terminal command, and
    a headline Control-click route satisfies that perfectly.
    """
    assert SUPERSEDED_MAC_ROUTE in note_text, (
        f"the note no longer mentions {SUPERSEDED_MAC_ROUTE} at all. It is still the route "
        "that works on macOS 14 and older, so it belongs here — second."
    )
    settings_route = CANONICAL_PHRASES["the macOS button to click"]
    assert note_text.index(settings_route) < note_text.index(SUPERSEDED_MAC_ROUTE), (
        f"{SUPERSEDED_MAC_ROUTE} appears before {settings_route!r} in {CANONICAL_NOTE.name}, "
        "so the note leads with the route macOS 15 removed. Whichever route comes first is "
        "the one most readers will try, and on any current Mac that one produces the same "
        "dead-end alert they are already looking at."
    )


# --------------------------------------------------------------------------------------
# The surfaces that used to restate it.
# --------------------------------------------------------------------------------------


def test_the_build_instructions_point_at_the_note(swept):
    """The instructions keep their troubleshooting table; the rows that used to restate
    the wording become links. Asserted as a link rather than a mention, because a
    maintainer reading a table needs somewhere to click.
    """
    text = BUILD_INSTRUCTIONS.read_text(encoding="utf-8")
    assert CANONICAL_NOTE.name in text, (
        f"{BUILD_INSTRUCTIONS.relative_to(REPO_ROOT)} does not name "
        f"{CANONICAL_NOTE.name}. Its Gatekeeper and SmartScreen rows previously carried "
        "their own copy of the instructions — one of which contradicted the build "
        "script — so with the copy gone the row has to point somewhere."
    )
    assert f"({CANONICAL_NOTE.name})" in text, (
        f"{BUILD_INSTRUCTIONS.relative_to(REPO_ROOT)} names {CANONICAL_NOTE.name} but "
        "not as a Markdown link target, so a reader of the rendered table cannot reach it."
    )


def test_the_dmg_ships_the_note_beside_the_app():
    """The one surface that reaches a macOS user in time.

    Gatekeeper's alert fires when the app is opened, which is *after* the DMG window is
    already on screen — so the note has to be in the staging directory the script
    assembles, and it has to have a position of its own in that window or it lands under
    the app icon.

    Parsed from the script rather than executed: a real run downloads two Python
    distributions. The precedent is the vendor-drift guard, which parses this same
    script's copy list, and which fails when it cannot find one rather than passing over
    an empty parse.
    """
    script = DMG_BUILD_SCRIPT.read_text(encoding="utf-8")

    staged = [
        line
        for line in script.splitlines()
        if "STAGING_DIR" in line and STAGING_NOTE_NAME in line and CANONICAL_NOTE.name in line
    ]
    assert staged, (
        f"{DMG_BUILD_SCRIPT.relative_to(REPO_ROOT)} never copies {CANONICAL_NOTE.name} "
        f'into the staging directory as "{STAGING_NOTE_NAME}". The DMG window is the '
        "only surface a Mac recipient sees before the alert, so the note has to be in it."
    )

    assert f'--icon "{STAGING_NOTE_NAME}"' in script, (
        f'the script stages the note but gives "{STAGING_NOTE_NAME}" no --icon position '
        "in the create-dmg layout. create-dmg places unlisted items itself, which can "
        "put the note under the app icon in the very window it exists to be seen in."
    )


# --------------------------------------------------------------------------------------
# The publisher URL, and the cheap rule that covers it.
# --------------------------------------------------------------------------------------


def test_the_windows_installer_publisher_url_is_the_public_repository():
    """Read off the ``[Setup]`` key rather than swept for, so that a URL moved into a
    comment cannot answer for the one Windows renders in Add/Remove Programs.
    """
    text = WINDOWS_INSTALLER.read_text(encoding="utf-8")
    match = re.search(r"^AppPublisherURL=(.*)$", text, flags=re.MULTILINE)
    assert match, (
        f"{WINDOWS_INSTALLER.relative_to(REPO_ROOT)} has no AppPublisherURL. Windows "
        "renders it in Add/Remove Programs; issue #262 owns it because no sweep can see "
        "a URL that names neither repository correctly."
    )
    assert match.group(1).strip() == PUBLIC_REPO, (
        f"AppPublisherURL is {match.group(1).strip()!r}, which is not the public "
        f"repository {PUBLIC_REPO}. It is rendered to every Windows user who opens "
        "Add/Remove Programs, and it has now been found broken three times."
    )


def test_no_tracked_file_carries_the_placeholder_account_name(swept):
    """The covering rule, and the reason it needs no exemptions.

    Issue #234 declined the widened version — *every* github.com URL naming a variantalker
    repo must name the public one — because it would need an exemption list for genuine
    third-party links, and an exemption list is the thing these sweeps exist without. This
    rule is narrower and free: :data:`PLACEHOLDER_ORG` is a scaffolding placeholder, so
    there is no legitimate occurrence to exempt, and any occurrence is a URL that was
    broken when it was written.
    """
    homes = _files_containing(PLACEHOLDER_ORG, swept)
    assert homes == [], (
        f"the {PLACEHOLDER_ORG!r} placeholder still appears in {homes}. It is a broken "
        "URL wherever it occurs, and it is invisible to the private-repository sweep "
        "because it names neither repository. Issue #262 fixed the one in "
        f"{WINDOWS_INSTALLER.relative_to(REPO_ROOT)}; this claim is what keeps the next "
        "one from going unnoticed."
    )
