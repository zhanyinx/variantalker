"""The check that refuses to publish a draft the tree has moved on from.

Issue #328, and it exists because of a specific near-miss. `mafigate-v1.0.0` was drafted from
what was then the tip of public ``main``. Twenty-two commits later it was not, and four of them
changed code that ships — one being the fix for a variant dialog that opened the **wrong
variant** (#310). The draft sat there naming its own commit accurately, correctly built,
correctly hashed, with the whole suite green; it read exactly as a fresh one reads, and the
handover comment on #265 said to publish it.

So the property under test is not "does the script run". It is **does it refuse**, in the state
that actually occurred — which is why :func:`verdict` is a pure function taking the six facts
its caller had to go and find out. Every case below is otherwise unreachable: reproducing them
against a live repository would mean drafting a release and then moving a branch under it.

**The one that must not drift.** ``test_it_refuses_even_when_nothing_that_ships_changed``.
Whether the intervening commits touched an installer's copy list decides how the refusal is
*worded* and nothing else. A version number is a claim about a tree, and "the artifacts would
probably be identical" is not that claim — so the temptation to pass a docs-only drift is the
one change here that would quietly restore the hole this ticket closed.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

STREAMLIT_APP = Path(__file__).resolve().parents[1]
SCRIPT = STREAMLIT_APP / "build" / "release_preflight.py"


def _preflight():
    """``build/release_preflight.py`` as a module, loaded by path.

    The same route ``test_build_identity.py`` takes to ``build/build_stamp.py``: these scripts
    deliberately live outside the app's import graph — a release machine should need nothing
    installed — so a test reaches them by file rather than by package.
    """
    spec = importlib.util.spec_from_file_location("mafigate_release_preflight", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


#: A healthy release: the tag names the tip of the published branch, and that branch was
#: exported from where the private tree actually is.
FRESH = dict(
    tag="mafigate-v1.0.0",
    tag_commit="a" * 40,
    branch="main",
    branch_commit="a" * 40,
    behind=0,
    is_draft=True,
    shipped_changes=[],
    exported_from="c" * 40,
    private_head="c" * 40,
)

#: **The first incident, with every commit real** — issue #328, and the reason issue #405
#: exists. `mafigate-v1.0.0` was cut at public `0b81ec4`; `dev_marco` then took 22 commits,
#: four of them shipped code including the wrong-variant fix; and **nobody exported in
#: between**.
#:
#: So the published branch never moved: the tag and the tip are the *same commit*, `behind` is
#: 0, and the compare between them is empty. That is what the first version of this guard
#: looked at, and why it would have passed. The staleness is visible only in the last two
#: fields — `0b81ec4` was exported from `139fe9c`, while the private tree had reached
#: `b4537d8`.
#:
#: The previous fixture set `branch_commit` to `"b" * 40`, a commit that never existed, and so
#: described a world in which the published branch had moved on. It had not. That single
#: invented value is what made a guard that could not catch this look like it had.
THE_FIRST_INCIDENT = dict(
    tag="mafigate-v1.0.0",
    tag_commit="0b81ec4db99fd28ccad46b355004726e29246b5b",
    branch="main",
    branch_commit="0b81ec4db99fd28ccad46b355004726e29246b5b",
    behind=0,
    is_draft=True,
    shipped_changes=[],
    exported_from="139fe9cf8fceb9f015a791e66d8897f0782af0db",
    private_head="b4537d8",
)

#: **The second, found by rechecking** — the same shape four days later, and the state the
#: live tree was in when issue #405 was written: the re-cut draft at `40d28a0`, exported from
#: `7d89a52`, with `dev_marco` 74 commits further on at `4203535`.
THE_SECOND_INCIDENT = dict(
    tag="mafigate-v1.0.0",
    tag_commit="40d28a042b256b3aa17d687a2f39b688bb6a094a",
    branch="main",
    branch_commit="40d28a042b256b3aa17d687a2f39b688bb6a094a",
    behind=0,
    is_draft=True,
    shipped_changes=[],
    exported_from="7d89a5212189bbb8e6dcef6bf2fdf970994cc559",
    private_head="4203535",
)

#: A tag pushed at an *older* public commit — a different mistake, and one this guard has
#: always caught. Named hypothetical because it is: it has not happened here, so its commits
#: are placeholders rather than history.
HYPOTHETICAL_TAG_BEHIND_THE_BRANCH = dict(
    tag="mafigate-v1.0.0",
    tag_commit="a" * 40,
    branch="main",
    branch_commit="d" * 40,
    behind=3,
    is_draft=True,
    shipped_changes=["streamlit_app/components/variant_table.py"],
    exported_from="c" * 40,
    private_head="c" * 40,
)


def test_the_script_is_there_and_loads():
    """Loaded rather than assumed: everything below reads through this module."""
    module = _preflight()
    assert callable(module.verdict)
    assert (module.OK, module.REFUSED, module.USAGE, module.UNKNOWN) == (0, 1, 2, 3)


def test_a_draft_built_from_the_tip_is_safe_to_publish():
    """The pass has to be reachable, or the refusals below prove nothing."""
    module = _preflight()
    code, lines = module.verdict(**FRESH)
    assert code == module.OK, lines
    assert any("Safe to publish" in line for line in lines), lines


def test_it_refuses_both_states_that_actually_occurred():
    """The two real ones, read off the public repository's own history.

    Neither has the published branch moving. That is the whole point: the tag and the tip are
    the same commit in both, so everything the first version of this guard looked at agreed,
    and it said *safe to publish* — on the incident it was written for, and again four days
    later. What is wrong is one level further out, and only the last two fields show it.
    """
    module = _preflight()
    for name, facts in (
        ("the first", THE_FIRST_INCIDENT),
        ("the second", THE_SECOND_INCIDENT),
    ):
        code, lines = module.verdict(**facts)
        body = "\n".join(lines)
        assert code == module.REFUSED, f"{name} incident was not refused:\n{body}"
        assert "an export is pending" in body, body
        assert facts["exported_from"][:7] in body, body
        assert facts["private_head"][:7] in body, body


def test_the_tag_and_the_tip_agreeing_is_not_enough():
    """Stated as its own property, because it is the assumption the first version rested on."""
    module = _preflight()
    facts = THE_FIRST_INCIDENT
    assert facts["tag_commit"] == facts["branch_commit"], (
        "this fixture is only interesting while the tag and the tip are the same commit; if "
        "they differ it is testing the other rule"
    )
    assert facts["behind"] == 0 and not facts["shipped_changes"], (
        "and while there is nothing between them to notice"
    )
    assert module.verdict(**facts)[0] == module.REFUSED


def test_a_tag_behind_the_published_branch_is_still_refused():
    """The rule the guard already had, kept — a different mistake needing a different fix."""
    module = _preflight()
    code, lines = module.verdict(**HYPOTHETICAL_TAG_BEHIND_THE_BRANCH)
    body = "\n".join(lines)
    assert code == module.REFUSED, body
    assert "commit(s) behind" in body, body
    assert "an export is pending" not in body, (
        "a stale export and a mis-placed tag need different messages: they have different fixes"
    )


def test_it_refuses_even_when_nothing_that_ships_changed():
    """The load-bearing one for the *tag behind the branch* rule — see the module docstring.

    A drift of documentation only is still a draft built from a tree that is not the published
    one. Passing it would make the check advisory, and an advisory check is what #265 already
    had: a paragraph in the build instructions.
    """
    module = _preflight()
    code, lines = module.verdict(
        **{**HYPOTHETICAL_TAG_BEHIND_THE_BRANCH, "shipped_changes": []}
    )
    body = "\n".join(lines)
    assert code == module.REFUSED, body
    assert "probably be identical" in body, body
    assert "'probably' is not what a version number means" in body, body


def test_a_public_clone_cannot_answer_and_says_so():
    """The private tree is the only thing that knows what is waiting to be exported.

    Absent, the honest answer is *unknown*, not *fine* — the same rule as an unreadable
    repository below.
    """
    module = _preflight()
    code, lines = module.verdict(**{**FRESH, "private_head": None})
    body = "\n".join(lines)
    assert code == module.UNKNOWN, body
    assert "could not be checked" in body, body


def test_a_tip_that_did_not_come_from_the_export_is_not_a_pass():
    """Someone committed to the published branch by hand, so *exported from* is unanswerable."""
    module = _preflight()
    code, lines = module.verdict(**{**FRESH, "exported_from": None})
    assert code == module.UNKNOWN, lines


@pytest.mark.parametrize(
    "message, expected",
    [
        ("Sync from variantalker" + "_ieo" + " @ 4203535\n\nbody", "4203535"),
        ("Merge pull request #7 from zhanyinx/dev_marco", None),
        ("", None),
    ],
)
def test_the_private_commit_is_read_off_the_published_tip(message, expected):
    """The export writes that line on every commit it makes (#227); this reads it back."""
    module = _preflight()
    assert module.exported_from(message) == expected


def test_the_sync_prefix_is_not_spelled_out_in_this_tree():
    """It is composed from halves, as `test_public_repo_name.py` requires of every file here."""
    module = _preflight()
    source = SCRIPT.read_text(encoding="utf-8")
    assert module.SYNC_PREFIX.strip().endswith("@")
    assert "variantalker" + "_ieo" not in source, (
        "release_preflight.py spells the private repository's name, which the tree-wide sweep "
        "asserts occurs zero times — compose it from halves, as the export tooling does"
    )


def test_an_already_published_release_is_not_refused():
    """There is nothing left to protect, and a refusal here would teach people to ignore it."""
    module = _preflight()
    code, lines = module.verdict(**{**THE_FIRST_INCIDENT, "is_draft": False})
    assert code == module.OK, lines
    assert any("already published" in line for line in lines), lines


def test_no_release_at_all_is_not_a_pass():
    """`I could not check` must never read like `I checked`."""
    module = _preflight()
    code, lines = module.verdict(**{**FRESH, "is_draft": None})
    assert code == module.UNKNOWN, lines
    assert any("no release exists" in line for line in lines), lines


def test_a_repository_it_cannot_read_is_not_a_pass(monkeypatch):
    """The other half of the same rule, at the level that talks to the network."""
    module = _preflight()
    monkeypatch.setattr(module, "gather", lambda tag, repo: None)
    assert module.main(["--tag", "mafigate-v9.9.9"]) == module.UNKNOWN


@pytest.mark.parametrize(
    "path, ships",
    [
        ("streamlit_app/components/variant_detail.py", True),
        ("streamlit_app/MAFigate.py", True),
        ("streamlit_app/vendor/pipeline_utils.py", True),
        ("streamlit_app/.streamlit/config.toml", True),
        # Neither installer copies these, and test_vendor_drift.py fails one that does.
        ("streamlit_app/tests/test_components.py", False),
        ("streamlit_app/tests/fixtures/parity/somatic.maf", False),
        ("streamlit_app/build/BUILD_INSTRUCTIONS.md", False),
        # Outside the app altogether.
        ("bin/filter_variants.py", False),
        (".github/workflows/mafigate-release.yml", False),
        ("README.md", False),
    ],
)
def test_which_paths_count_as_shipping(path, ships):
    """Only the *wording* of a refusal depends on this, which is why it may be approximate."""
    module = _preflight()
    assert (module.shipped([path]) == [path]) is ships


def test_a_public_clone_is_told_apart_by_the_export_tooling(tmp_path):
    """The discriminator, exercised — a mutation removing it was otherwise invisible.

    Getting this wrong fails *safe*: a public clone's HEAD is never a private SHA, so the
    comparison refuses rather than passes. But it refuses saying *an export is pending*, which
    is the wrong diagnosis handed to someone who cannot act on it.
    """
    module = _preflight()
    (tmp_path / "tools").mkdir()
    assert module._private_head(tmp_path) is None, (
        "a tree with no tools/export_public.py is not the private tree, whatever git says"
    )


def test_the_default_target_is_the_public_repository():
    """A preflight that read the target from the local git config would check a fork."""
    module = _preflight()
    assert module.DEFAULT_REPO == "zhanyinx/variantalker"
