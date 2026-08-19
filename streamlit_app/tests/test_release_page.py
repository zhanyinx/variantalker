"""What the release page must say, checked on the page a reader actually gets.

Issue #265. ``build/RELEASE_PAGE.md`` is the body of every MAFigate release, and it is the
only surface that reaches a Windows recipient at all: SmartScreen fires on the downloaded
installer before any screen of ours can be drawn. It is also where someone who has read the
paper and has no MAF of their own finds that out — before spending the download rather than
after, which is the decision the delivery matrix (#229) made and called copy rather than code.

Copy that load-bearing is guarded like code. Four things the matrix decided are asserted here,
and each one is a claim a maintainer could quietly drop while the page still read well:

* **The MAF prerequisite comes first.** Not merely present — *above* the downloads. A
  prerequisite below the file you have already downloaded is a note about what went wrong.
* **The three failure modes are named**, including the two that do not look like failures: an
  install that reports success and then will not start (no internet on first launch), and a
  process the operating system stops with no message at all (a MAF too large for free memory).
* **The memory rule is given as something to check**, with the numbers #225 measured.
* **No update check is promised**, because there is none and there will never be one — it would
  be the app's first outbound call and would falsify the data-posture claim two paragraphs
  above it (#229).

**Asserted on the assembled page, not on the source file.** The page has two placeholders that
CI fills, so the file on disk is not what anyone reads: an ordering rule checked against the
source would be checking a document with a hole where the downloads go. The assembly is run
through ``test_release_workflow._notes_step`` — imported rather than re-implemented, which is
this suite's standing rule for a guard that has to agree with a script (#247: a sweep that
copies the scanner's loop goes green on a file the real one refuses).

**What is deliberately *not* here.** The first-open wording is not checked for presence in the
page's source, and must not be: ``tests/test_unsigned_artifact_copy.py`` sweeps every tracked
file for those phrases and fails a second home, so this page carries a placeholder and would
fail that sibling the moment it grew a paraphrase. Two guards, one rule, no overlap. Version
literals are likewise already covered — ``tests/test_installer_version.py`` sweeps
``streamlit_app/build/`` by prefix, which is why the page names no version and lets the
generated filenames carry it.
"""

from __future__ import annotations

import re

import pytest

# The assembly, shared rather than repeated. See the module docstring.
from tests.test_release_workflow import ONE_OF_EACH, STREAMLIT_APP, _notes_step

# The data-posture sentence and the fact that earns it, imported from the guard that owns
# both. Respelling either here would let this page keep promising something after the README's
# switch had flipped the other way — two surfaces making one promise, drifting apart, which is
# the failure #229 chose a single pinned wording to prevent.
from tests.test_delivery_channels_copy import (
    UNQUALIFIED_CLAIM,
    telemetry_is_turned_off_in_the_tree,
)

PAGE = STREAMLIT_APP / "build" / "RELEASE_PAGE.md"

#: The two placeholders CI fills, written as they must appear: on a line of their own.
PLACEHOLDERS = ("{{ downloads }}", "{{ opening_the_app }}")

#: The header row of the generated downloads block. What "above the download" is measured
#: against — the table is where the filenames appear, so it is where a reader stops reading
#: prose and starts clicking.
DOWNLOADS_MARKER = "| file | platform | sha256 |"

#: Phrases that would promise an update mechanism. ``check for update`` covers the app doing
#: it; the other two cover a page telling someone to let it.
UPDATE_PROMISES = (
    "check for updates",
    "checks for updates",
    "automatically update",
    "update automatically",
    "auto-update",
)


@pytest.fixture(scope="module")
def page_source():
    return PAGE.read_text(encoding="utf-8")


@pytest.fixture(scope="module")
def assembled(tmp_path_factory):
    """The release page as CI writes it, placeholders filled."""
    completed, notes, _ = _notes_step(tmp_path_factory.mktemp("notes"), ONE_OF_EACH)
    assert completed.returncode == 0, (
        "the release page did not assemble, so nothing below is checking a page anyone would "
        f"read:\n{completed.stdout}\n{completed.stderr}"
    )
    return notes


def _plain(text):
    """Markdown emphasis and code ticks stripped, lowercased, for phrase matching."""
    return re.sub(r"[*`_]", "", text).lower()


class TestThePageIsAssembledFromItsPlaceholders:
    """The mechanism, before anything is asserted through it."""

    def test_each_placeholder_appears_exactly_once_on_a_line_of_its_own(self, page_source):
        lines = page_source.splitlines()
        for placeholder in PLACEHOLDERS:
            own_line = [line for line in lines if line.strip() == placeholder]
            assert len(own_line) == 1, (
                f"{placeholder} must appear exactly once, on a line of its own — the workflow "
                f"matches whole lines, and found {len(own_line)}"
            )

    def test_the_assembled_page_holds_no_placeholder(self, assembled):
        assert "{{" not in assembled, (
            "the assembled page still carries a placeholder, so a release would ship braces "
            "where the downloads or the first-open note belong"
        )

    @pytest.mark.parametrize(
        "page_text, description",
        [
            pytest.param(
                "# Page\n\nProse only.\n", "no placeholders at all", id="none"
            ),
            pytest.param(
                "# Page\n\n{{ downloads }}\n", "no first-open note", id="note-missing"
            ),
            pytest.param(
                "# Page\n\n{{ opening_the_app }}\n", "no downloads", id="downloads-missing"
            ),
            pytest.param(
                "# Page\n\n{{ downloads }}\n\n{{ downloads }}\n\n{{ opening_the_app }}\n",
                "the downloads twice",
                id="downloads-twice",
            ),
            pytest.param(
                "# Page\n\n{{ download }}\n\n{{ opening_the_app }}\n",
                "a renamed placeholder",
                id="renamed",
            ),
        ],
    )
    def test_the_assembly_refuses_a_page_it_cannot_fill(self, page_text, description, tmp_path):
        """The failure this exists for is silent: a page that assembles without its downloads
        block is a release page with no files named on it, and nothing about it looks wrong.
        """
        completed, notes, _ = _notes_step(tmp_path, ONE_OF_EACH, page_text=page_text)
        assert completed.returncode != 0, (
            f"the assembly accepted a page with {description}; it wrote:\n{notes}"
        )


class TestTheReaderIsToldWhatTheyNeedBeforeTheyDownload:
    """The matrix's answer to the visitor who arrives with no MAF (#229, #231)."""

    def test_the_prerequisite_is_stated_above_the_downloads(self, assembled):
        prerequisite = _plain(assembled).find("annotated")
        downloads = assembled.find(DOWNLOADS_MARKER)
        assert prerequisite != -1, "the page never says an annotated MAF is needed"
        assert downloads != -1, (
            "the downloads block is not on the assembled page, so its position cannot be "
            "compared with anything"
        )
        assert prerequisite < downloads, (
            "the MAF prerequisite must come before the downloads. Below them it is a note "
            "about why the download was wasted."
        )

    def test_the_prerequisite_names_the_three_kinds_of_file_the_app_opens(self, assembled):
        before = _plain(assembled.split(DOWNLOADS_MARKER)[0])
        for extension in (".maf", ".tsv", ".txt"):
            assert extension in before, (
                f"the prerequisite does not name {extension}, which the app's own chooser "
                "accepts — a reader with a .tsv would read this as not applying to them"
            )

    def test_the_page_says_where_such_a_file_comes_from(self, assembled):
        before = _plain(assembled.split(DOWNLOADS_MARKER)[0])
        assert "pipeline" in before, (
            "the prerequisite never names the pipeline the file comes from, which leaves a "
            "reader who has none with nothing to do about it"
        )

    def test_the_page_does_not_offer_a_file_to_try_it_with(self, assembled):
        """#231 parked the bundled demo MAF, so promising one here would be a 404 in prose."""
        assert "example file" not in _plain(assembled) or "no example file" in _plain(assembled)


class TestTheThingsThatGoWrongAreNamed:
    """Three failure modes, two of which do not present as failures at all."""

    def test_the_first_open_is_said_to_be_blocked(self, assembled):
        assert "unsigned" in _plain(assembled) or "not signed" in _plain(assembled)
        assert "blocked" in _plain(assembled)

    def test_an_offline_first_launch_is_said_to_install_and_then_not_start(self, assembled):
        plain = _plain(assembled)
        assert "internet" in plain, "the page never mentions needing the network at all"
        assert re.search(r"(offline|no internet|without internet)", plain), (
            "the page never names the offline case"
        )
        assert re.search(r"(fails to start|does not start|will not start)", plain), (
            "the page does not say what an offline first launch actually does — the install "
            "completes and reports success, and only the launch fails, which is why this reads "
            "as a broken app rather than a missing download"
        )

    def test_a_file_too_large_is_said_to_stop_the_app_with_no_message(self, assembled):
        plain = _plain(assembled)
        assert re.search(r"(no message|no error|no diagnostic|without a message)", plain), (
            "the page does not say that running out of memory kills the app silently. Someone "
            "who is not told this reports it as a crash, and the answer — a smaller file or a "
            "bigger machine — is not one they can guess."
        )


class TestTheMemoryRuleIsUsable:
    """#225's law, given as something to check rather than as a benchmark."""

    def test_the_rule_is_stated_with_both_of_its_terms(self, assembled):
        plain = _plain(assembled)
        assert "140" in plain, "the rule's fixed term (about 140 MB) is not stated"
        assert re.search(r"(ten times|10 times|10×|10x)", plain), (
            "the rule's per-file term — ten times the file's size on disk — is not stated, so "
            "the number cannot be applied to the file the reader actually has"
        )

    def test_the_largest_realistic_file_is_costed(self, assembled):
        assert re.search(r"1\.[34]\s*GB", assembled, re.IGNORECASE), (
            "the page does not say what the largest MAFs need in practice (~1.4 GB free). The "
            "rule alone leaves a clinician multiplying, at the moment they are least likely to."
        )


class TestNoUpdateCheckIsPromised:
    """There is none, and there never will be one — see #229."""

    @pytest.mark.parametrize("promise", UPDATE_PROMISES)
    def test_the_page_promises_no_update_mechanism(self, assembled, promise):
        plain = _plain(assembled)
        occurrences = [
            line for line in plain.splitlines() if promise in line and "never" not in line
        ]
        assert not occurrences, (
            f"the page implies MAFigate {promise!r}: {occurrences}. It makes no outbound call "
            "of any kind, and an update check would be the first — which is exactly what the "
            "data-posture claim on this page says it does not do."
        )

    def test_the_page_says_outright_that_it_never_checks(self, assembled):
        assert re.search(r"never checks", _plain(assembled)), (
            "the absence of an update check is left to be inferred. It is a privacy property "
            "worth stating, and stating it is what stops a reader waiting for a prompt that is "
            "never coming."
        )


class TestTheDataPostureClaimIsTheOneTheMatrixChose:
    """One wording across every surface, and it has to stay earned."""

    def test_the_page_makes_the_pinned_claim(self, assembled):
        assert UNQUALIFIED_CLAIM in _plain(assembled), (
            f"the release page does not carry {UNQUALIFIED_CLAIM!r}. #229 chose one wording so "
            "that the README, the installers and this page cannot drift apart; a paraphrase "
            "here is a second wording."
        )

    def test_the_claim_is_still_earned_by_the_tree(self):
        """The same switch the README's copy hangs on, read from the same function.

        If someone removes the Streamlit telemetry config, Streamlit reports usage again and
        this page's promise becomes false. That must fail *here*, on the surface making the
        promise, rather than only in the README's guard.
        """
        assert telemetry_is_turned_off_in_the_tree(), (
            "nothing in the tree turns Streamlit's usage reporting off any more, so the release "
            f"page may not say {UNQUALIFIED_CLAIM!r}. Restore the config, or rewrite the page's "
            "posture section — but the two have to agree."
        )


class TestThePageNamesNoChannelThatDoesNotExist:
    """#229 removed the hosted channel permanently; #231 parked the demo."""

    @pytest.mark.parametrize(
        "absent",
        ["streamlit.io", "share.streamlit", "community cloud", "hosted demo", "try it online"],
    )
    def test_the_page_does_not_mention_it(self, assembled, absent):
        assert absent not in _plain(assembled), (
            f"the release page mentions {absent!r}. That channel was ruled out permanently, so "
            "naming it here sends a reader somewhere that does not exist."
        )
