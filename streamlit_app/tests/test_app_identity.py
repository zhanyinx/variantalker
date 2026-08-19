"""What the app says it is, and where it says its version.

The app described itself in three rendered voices and printed its version in six places
(issue #71). The voices agreed in meaning and disagreed in words — a header tagline, an
About blurb and a sidebar footer, two of them differing only by a truncation — and two of
the About blurb's five feature bullets were **untrue of the app a user can actually drive**:
it advertised exports "in multiple formats" when every reachable variant download is CSV,
and "clinical database integration" when the app integrates with nothing and reads the
annotations already present in your file.

Copy has no compiler, so this file is the only thing standing between the app and a fourth
voice. What it holds is the *settings* rather than the rendering: that both surviving
surfaces read one constant, that the version reaches exactly one of them, and that the two
retired claims are gone by name. The last of those is the one with teeth — a false claim is
easy to re-add in a hurry, and reads as documentation while it sits there.

The header's own sentence is asserted as an identity against ``APP_TAGLINE``, never quoted:
a test that repeats the sentence is a fourth copy of it, which is the thing being fixed.

Since issue #263 the About dialog names two more facts — the **channel** an artifact arrived
by and the **build** it was cut from — and the one-surface sweep covers all three rather than
only the version. That was the constraint the ticket was written under: the rule exists
because review found a hand-typed version in the Help FAQ, and a rule that stops at the
oldest of the three strings it should cover is a rule half-applied. The mechanism behind the
two new strings — the gitignored stamp, the writer, the source-checkout default — belongs to
``tests/test_build_identity.py``; what is asserted here is only *where they render*.
"""

from __future__ import annotations

import os
import sys
from unittest.mock import MagicMock, patch

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import MAFigate  # noqa: E402
from config.build_identity import build_identity  # noqa: E402
from config.constants import APP_NAME, APP_TAGLINE, APP_VERSION  # noqa: E402

# Borrowed rather than copied. Its own docstring explains why a bare `Mock()` will not
# stand in for `session_state` here — `.get()` answering every key with one value, and `in`
# raising rather than failing — and a second copy of that reasoning would be a second thing
# to keep true. It lives in `tests/fakes.py` rather than in the UI-component tests that
# declared it first, so that borrowing it does not make a collected module into a library.
from tests.fakes import FakeSessionState, page_config_kwargs  # noqa: E402


def _strings_passed_to(mock) -> list[str]:
    """Every string argument this mock was called with, at any depth.

    ``mock_calls`` reaches nested attributes, so one sweep of ``st`` catches
    ``st.markdown``, ``st.sidebar.markdown`` and anything a later edit adds without this
    helper having to name it. Only positional and keyword *strings* are collected; the
    frames and containers that also flow through these calls are not copy.
    """
    found = []
    for call in mock.mock_calls:
        found.extend(a for a in call.args if isinstance(a, str))
        found.extend(v for v in call.kwargs.values() if isinstance(v, str))
    return found


def _surfaces_naming(token: str, surfaces) -> dict[str, list[str]]:
    """Which of *surfaces* drew a line containing *token*, as ``name -> the offending lines``.

    One sweep for all three identity strings, rather than the version's own loop plus a copy
    of it per string added later. *surfaces* is ``(name, drawn)`` pairs so a failure names the
    screen a reader would see it on, and so
    :func:`test_the_identity_sweep_catches_a_planted_string` can hand it a planted one.
    """
    found = {}
    for name, drawn in surfaces:
        hits = [line for line in drawn if token in line]
        if hits:
            found[name] = hits
    return found


def _rendered_surfaces(header, sidebar, help_page):
    """The three surfaces the one-surface rule is enforced against, named for a failure message.

    The browser tab is not among them because it is not a list of drawn strings — it is one
    keyword, asserted directly by each caller.
    """
    return (
        ("The header", header),
        ("The sidebar", sidebar),
        ("The Help page", help_page),
    )


#: This checkout's channel and build, read the way the app reads them. In the suite's ordinary
#: state — no stamp, because it is gitignored and no build has run here — that is the
#: source-checkout default, which is the state every contributor and every CI run is in. Read
#: rather than written out, so a developer who *has* just built an installer locally runs these
#: sweeps against that stamp instead of a red herring: what is asserted below is that these two
#: strings render in one place, whatever they happen to say.
IDENTITY = build_identity()


@pytest.fixture
def header(monkeypatch):
    """Render the header with Streamlit mocked out, and return what it drew.

    ``show_discarded_cache_banner`` is stubbed because it belongs to issue #40 and reads
    its *own* module's ``st``, so patching this module's would not contain it. The session
    state is empty on purpose: ``cache_loaded`` absent is the ordinary case, and it keeps
    the cache-restored banner — which is not a description of the app — out of the sweep.
    """
    with patch("MAFigate.st") as mock_st, patch("MAFigate.show_discarded_cache_banner"):
        mock_st.session_state = FakeSessionState()
        MAFigate.render_header()
        yield _strings_passed_to(mock_st)


@pytest.fixture
def page_config():
    """The keyword arguments the app hands ``st.set_page_config``.

    The driving moved to ``tests/fakes.py`` when ``test_public_repo_name.py`` needed the same
    dict for a different question — the same move ``FakeSessionState`` made above, and for the
    same reason.
    """
    return page_config_kwargs()


@pytest.fixture
def help_page():
    """Every string the Help page draws, by running the page itself.

    The first version of this fixture did not exist: the Help page was asserted at its
    module namespace instead — ``not hasattr(help_page, "APP_VERSION")`` — which is a check
    on the *import* and therefore blind to the one case that was actually there. The
    citation FAQ had the version written out by hand, as a literal, so nothing about the
    module's imports could reveal it and the guard passed over a live violation. Review
    caught it; this fixture is the repair, and it is why every surface here is now swept by
    running it rather than by reading it.

    The mock is sized where the page unpacks: ``st.columns`` is called with 4 and with 2 on
    different rows, and ``st.tabs`` with five names, so a fixed-length return value fails
    on the second caller rather than the first.
    """
    from page_modules import help as help_module

    with patch("page_modules.help.st") as mock_st:
        mock_st.session_state = FakeSessionState()
        mock_st.columns.side_effect = lambda spec, **kw: [
            MagicMock() for _ in range(spec if isinstance(spec, int) else len(spec))
        ]
        mock_st.tabs.return_value = [MagicMock() for _ in range(5)]
        mock_st.button.return_value = False
        mock_st.text_input.return_value = ""
        mock_st.selectbox.return_value = "All"
        help_module.show_help_page()
        yield _strings_passed_to(mock_st)


@pytest.fixture
def sidebar():
    """Every string the sidebar draws on an ordinary render, with nothing loaded.

    The radio is answered with the page already current so no rerun is provoked, and the
    buttons with ``False`` so no navigation fires: this is the resting state of the column,
    which is the state a footer would have sat in on every page.
    """
    from components import sidebar as sidebar_module

    with patch("components.sidebar.st") as mock_st:
        mock_st.session_state = FakeSessionState(current_page="home", maf_data=None)
        mock_st.sidebar.button.return_value = False
        mock_st.sidebar.radio.return_value = "🏠 Home"
        sidebar_module.create_sidebar_navigation()
        yield _strings_passed_to(mock_st)


def test_the_header_draws_the_one_shared_sentence(header):
    """Not *a* sentence saying roughly this — the constant itself, unedited."""
    assert APP_TAGLINE in header, (
        "render_header no longer draws APP_TAGLINE. The header and the About dialog say "
        "what the app is by reading one constant; a literal here is the third voice "
        "returning."
    )


def test_the_about_dialog_reads_the_same_sentence_as_the_header(page_config, header):
    """One sentence, two surfaces. They cannot drift apart if neither owns the words."""
    about = page_config["menu_items"]["About"]
    assert APP_TAGLINE in about
    assert APP_TAGLINE in header


@pytest.mark.parametrize(
    "retired",
    [
        # Advertised CSV/TSV/Excel. The one function offering all three,
        # `components.results_view.create_export_buttons`, had no caller and was deleted
        # in issue #83 — so the claim is now false by construction: every export a user
        # can reach is CSV, the run report is plain text, and the parameter files are JSON
        # or YAML. Nothing in the repo writes an .xlsx at all.
        "multiple formats",
        "Export capabilities",
        # The app queries no database. It reads the annotation columns your file already
        # carries — which is what the Home page's own copy says.
        "database integration",
        # A sales word, and this copy is for clinicians.
        "Advanced",
    ],
)
def test_the_about_dialog_makes_none_of_the_retired_claims(page_config, retired):
    about = page_config["menu_items"]["About"]
    assert retired not in about, (
        f"The About dialog says {retired!r} again. Every bullet in that dialog was checked "
        "against the code in issue #71 and this one did not survive the check."
    )


def test_the_about_dialog_is_the_only_place_the_version_renders(
    page_config, header, sidebar, help_page
):
    """Six places, then one.

    The browser tab, the header line, the sidebar footer, the Help page's subtitle and the
    citation answer in the Help page's FAQ all printed it. A release number is a fact a bug
    report needs and a clinician does not, and the About menu is where a bug reporter
    already looks — so it keeps the number and the other five surfaces do not.
    """
    assert APP_VERSION in page_config["menu_items"]["About"]

    assert page_config["page_title"] == APP_NAME, (
        "The browser tab is naming a version again. The tab is a label the user reads at a "
        "glance, not a release note."
    )

    offenders = _surfaces_naming(APP_VERSION, _rendered_surfaces(header, sidebar, help_page))
    assert not offenders, (
        f"a surface other than About is printing the version again: {offenders}. The About "
        "dialog is the one place that number belongs."
    )


def test_the_about_dialog_names_the_channel_and_the_build(page_config):
    """All three facts, on one line, in the dialog a bug reporter already opens (issue #263).

    Asserted as *one line* rather than three separate ``in`` checks, because the point of the
    two new fields is that a reporter can quote a single line and have it fully identify the
    build. Three facts scattered over a dialog would satisfy three ``in`` checks and still
    arrive in an issue as a version and a shrug.

    The values are read from ``config.build_identity`` rather than written out here, for the
    reason the tagline is never quoted in this file: a test that repeats the string it checks
    is another copy of it.
    """
    about = page_config["menu_items"]["About"]

    lines = [line for line in about.splitlines() if APP_VERSION in line]
    assert lines, f"the About dialog names no version at all:\n{about}"

    for fact in (IDENTITY.channel, IDENTITY.build):
        assert any(fact in line for line in lines), (
            f"the About dialog does not name {fact!r} beside the version. A version alone "
            "cannot identify a build: the same release reaches a user as a .dmg, as a "
            f"Windows .exe or as a clone.\n{about}"
        )


def test_the_about_dialog_is_the_only_place_the_channel_and_build_render(
    page_config, header, sidebar, help_page
):
    """The one-surface rule, extended to the two strings issue #263 added.

    The version's guard was written after review found a hand-typed copy in the Help FAQ, and
    the two new fields are more tempting to repeat than the version ever was — a build string
    reads like debug output, and debug output ends up in a sidebar caption. So they are held
    to the same rule from the start rather than added to the sweep after the first copy
    appears.

    The browser tab is checked here too, and not only in the version's test: a tab title is
    the one surface whose whole job is to be short, and ``MAFigate (macos-dmg)`` is exactly
    the kind of well-meant edit this rule exists to refuse.
    """
    surfaces = _rendered_surfaces(header, sidebar, help_page)

    for fact in (IDENTITY.channel, IDENTITY.build):
        assert fact in page_config["menu_items"]["About"], (
            f"About no longer names {fact!r}, so this sweep is proving that a string nothing "
            "renders renders nowhere."
        )
        offenders = _surfaces_naming(fact, surfaces)
        assert not offenders, (
            f"a surface other than About is naming {fact!r}: {offenders}. The channel and "
            "the build belong in the About dialog, with the version, and nowhere else."
        )
        assert fact not in page_config["page_title"], (
            f"the browser tab names {fact!r}. The tab is a label the user reads at a glance, "
            "not a build record."
        )


def test_the_sweep_reads_all_three_surfaces(header, sidebar, help_page):
    """Vacuity check one, and it guards a hole this file's own refactor opened.

    The three surfaces used to be written out inside the version's test. Issue #263 moved them
    into :func:`_rendered_surfaces` so that one sweep serves all three identity strings — which
    means deleting a line from that tuple now narrows *every* sweep at once, silently, and the
    planted-string check below would still pass. That is a worse failure than the duplication
    the refactor removed, so the tuple is pinned here.

    The three names are written out rather than read from the helper. A test that took them
    from the thing it is checking would assert only that the helper equals itself, and a
    surface dropped from it would take this assertion with it.

    Each surface must also have *drawn something*: a fixture that silently rendered nothing —
    an exception swallowed by a mock, a page that returned early — makes "the version appears
    nowhere but About" true of a screen nobody rendered.
    """
    surfaces = dict(_rendered_surfaces(header, sidebar, help_page))

    assert set(surfaces) == {"The header", "The sidebar", "The Help page"}, (
        f"_rendered_surfaces now names {sorted(surfaces)}. Every identity sweep in this file "
        "reads that tuple, so a surface missing from it is a surface nothing checks."
    )
    for name, drawn in surfaces.items():
        assert drawn, (
            f"{name} drew no strings at all, so every sweep over it passes for the wrong "
            "reason. The fixture is rendering nothing."
        )


@pytest.mark.parametrize("position", (0, 1, 2), ids=("header", "sidebar", "help page"))
def test_the_identity_sweep_catches_a_planted_string(header, sidebar, help_page, position):
    """Vacuity check two: the sweep finds a violation when there is one to find.

    Planted on **each** surface in turn, not just the first. One position proves the sweep
    reads that position; the failure being guarded against is a *particular* surface growing a
    build caption, and the sweep is only as good as the weakest slot it reads.

    The plant is a real surface's real strings with one line added — what a sidebar caption
    printing the build would actually look like — and the assertion is the exact dict, so a
    sweep that noticed the string but attributed it to the wrong screen fails too.
    """
    surfaces = list(_rendered_surfaces(header, sidebar, help_page))
    name, drawn = surfaces[position]

    for fact in (IDENTITY.channel, IDENTITY.build):
        line = f"Build: {fact}"
        planted = list(surfaces)
        planted[position] = (name, list(drawn) + [line])

        assert _surfaces_naming(fact, planted) == {name: [line]}, (
            f"the sweep did not notice {fact!r} planted on {name.lower()}, so it would not "
            "notice that surface printing it for real"
        )

    assert _surfaces_naming(IDENTITY.build, [("A planted surface", [IDENTITY.build])]), (
        f"the sweep does not match {IDENTITY.build!r} on its own; it is checking for whole "
        "lines rather than for the string appearing in one"
    )


def test_the_citation_answer_points_at_about_rather_than_naming_a_version(help_page):
    """A version a user is told to quote must not be a literal in the prose.

    The FAQ answered "How do I cite MAFigate?" with a worked example — ``MAFigate v2.0.0
    (VariantTalker IEO, 2025)`` — which is the failure ``tests/README.md`` names: derive it
    or guard it, never copy it. Hand-copied, it had already drifted twice, carrying a frozen
    year and a misspelling of the pipeline's own name. Asserted as *points at About* rather
    than as *lacks the current version*, because the sweep above only catches a stale
    literal while it happens to match ``APP_VERSION``.
    """
    citation = [line for line in help_page if "cite MAFigate" in line]
    assert citation, "the citation FAQ is gone; this guard now watches nothing"

    answer = citation[0]
    assert "About" in answer, (
        "The citation answer no longer sends the reader to About for the version. If it "
        "names a version itself, that number is a copy that will go stale."
    )


def test_the_sidebar_draws_no_footer(sidebar):
    """The three-line footer is gone, not merely trimmed.

    It restated the tagline drawn above it on the same screen, printed the version a second
    time, and credited Streamlit — three lines on every page, competing for the column that
    issue #58 had just filled with the file you have open and the way back to your results.
    """
    for gone in ("Made with", "Advanced MAF File Analysis"):
        hits = [line for line in sidebar if gone in line]
        assert not hits, f"{gone!r} is back at the bottom of the sidebar: {hits}"


def test_the_help_page_has_no_reason_to_know_the_version():
    """Belt to the sweep's braces, and a different failure: the *import* coming back.

    The Help page's subtitle read ``Comprehensive guide for **MAFigate v2.0.0**``; once the
    number came off, what remained said what the line directly beneath it already said, so
    the line went and with it the module's reason to know either constant. This is not a
    substitute for driving the page — asserting only this is what let the citation FAQ's
    hand-written version through — but a module that imports ``APP_VERSION`` again is a
    surface about to print it.
    """
    from page_modules import help as help_module

    assert not hasattr(help_module, "APP_VERSION"), (
        "page_modules/help.py imports APP_VERSION again. Nothing on that page needs it: the "
        "About dialog is where the number renders."
    )
