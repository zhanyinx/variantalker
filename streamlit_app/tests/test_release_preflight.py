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


#: The facts of a healthy release: drafted from the tip.
FRESH = dict(
    tag="mafigate-v1.0.0",
    tag_commit="a" * 40,
    branch="main",
    branch_commit="a" * 40,
    behind=0,
    is_draft=True,
    shipped_changes=[],
)

#: The facts as they actually were, the evening of 2026-08-19. The commits are the real ones.
THE_INCIDENT = dict(
    tag="mafigate-v1.0.0",
    tag_commit="0b81ec4db99fd28ccad46b355004726e29246b5b",
    branch="main",
    branch_commit="b" * 40,
    behind=22,
    is_draft=True,
    shipped_changes=[
        "streamlit_app/components/variant_detail.py",
        "streamlit_app/page_modules/results.py",
    ],
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


def test_it_refuses_the_state_that_actually_occurred():
    module = _preflight()
    code, lines = module.verdict(**THE_INCIDENT)
    assert code == module.REFUSED, lines
    body = "\n".join(lines)
    assert "22 commit(s) behind" in body, body
    assert "variant_detail.py" in body, (
        "the refusal must name what changed; a bare 'this is stale' leaves the reader deciding"
        " whether to believe it"
    )
    assert "older behaviour under the newer version number" in body, body


def test_it_refuses_even_when_nothing_that_ships_changed():
    """The load-bearing one — see the module docstring.

    A drift of documentation only is still a draft built from a tree that is not the published
    one. Passing it would make the check advisory, and an advisory check is what #265 already
    had: a paragraph in the build instructions.
    """
    module = _preflight()
    code, lines = module.verdict(**{**THE_INCIDENT, "shipped_changes": []})
    assert code == module.REFUSED, lines
    body = "\n".join(lines)
    assert "probably be identical" in body, body
    assert "'probably' is not what a version number means" in body, body


def test_an_already_published_release_is_not_refused():
    """There is nothing left to protect, and a refusal here would teach people to ignore it."""
    module = _preflight()
    code, lines = module.verdict(**{**THE_INCIDENT, "is_draft": False})
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


def test_the_default_target_is_the_public_repository():
    """A preflight that read the target from the local git config would check a fork."""
    module = _preflight()
    assert module.DEFAULT_REPO == "zhanyinx/variantalker"
