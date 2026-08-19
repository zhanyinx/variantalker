"""What the Help page may and may not tell a user, asserted by rendering it.

Issue #79 checked all five tabs of ``page_modules/help.py`` against the code and found
nineteen claims that did not hold. Most were corrected by editing prose, and prose is
exactly what drifts back — so the ones with a checkable counterpart in the code are pinned
here, against the strings the page actually draws.

**Rendering, not reading.** Almost every assertion below runs ``show_help_page`` and inspects
what it passed to Streamlit. Issue #71 is why: the version guard on this page asserted
``not hasattr(help, "APP_VERSION")``, a property of the module's *imports*, and passed over
a version number written out by hand in the FAQ. A claim is made by what renders.

Two tests deliberately do not render, and both are about a claim this page shares with
another surface rather than one it makes alone: ``test_the_parameter_tooltip_defers_to_the_
help_page_about_escat`` reads the tooltip beside the control, and
``test_every_escat_level_the_control_offers_is_defined`` reads ``ESCAT_DEFINITIONS`` against
``ESCAT_OPTIONS``. The second is here rather than in ``test_config.py``, whose territory
constants otherwise are, because it is what makes the rendering test above it *exhaustive*:
"every definition reaches the page" only covers every offered level while the keys are the
options. Split apart, each half would pass over a level that had quietly lost its gloss.

**Three kinds of test here**, and the difference matters:

* *No such number* — the page may not say 500MB anywhere. A literal, checkable forever.
* *Derived, not transcribed* — the core-column list and the ESCAT levels must equal the
  objects the app filters and validates by, so the page cannot come to describe a different
  MAF from the one the code refuses.
* *One story per fact* — where two surfaces glossed the same value differently, neither may
  reintroduce the gloss that lost.

What is deliberately **not** here: the accuracy of prose with no counterpart in the code.
Whether "ESCAT ranks how actionable a target is" is a good sentence is not a thing a test
can hold, and a test asserting the sentence verbatim would pin the wording rather than the
claim. #79's finding was that the checkable half is what the check catches.

Issue #89 moved one claim across that line without erasing it. The ESCAT definitions now
have a counterpart — ``ESCAT_DEFINITIONS``, condensed from the paper that defines the scale
— so what is checkable is that every offered level is defined, that both halves of each
definition reach the page, that the citation is drawn beside them, and that no surface
acquires a second copy. Whether the condensation is *faithful* to Table 2 still is not, and
no test here claims otherwise: that was read against the paper once, by a human-readable
citation left in the constant so it can be read again.
"""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

from config.columns import REQUIRED_COLUMNS
from config.vocabularies import (
    CIVIC_DEFINITIONS,
    CIVIC_OPTIONS,
    CIVIC_SOURCE,
    CIVIC_SOURCE_URL,
    CLINVAR_OPTIONS,
    ESCAT_DEFINITIONS,
    ESCAT_OPTIONS,
    ESCAT_SOURCE,
    ESCAT_SOURCE_DOI,
    ESCAT_STRONGEST,
)
from tests.fakes import FakeSessionState


def _strings_drawn_to(mock) -> list:
    """Every string handed to a Streamlit call on ``mock``, including its sub-objects.

    The page draws into tab and column context managers, which are ``MagicMock`` children,
    so walking ``mock_calls`` on the root is what reaches them: a claim made inside
    ``with tab5:`` is as rendered as one made at the top level, and the FAQ — where four of
    #79's findings lived — is drawn entirely inside one.
    """
    strings = []
    for call in mock.mock_calls:
        for argument in call.args:
            if isinstance(argument, str):
                strings.append(argument)
        for argument in call.kwargs.values():
            if isinstance(argument, str):
                strings.append(argument)
    return strings


@pytest.fixture
def help_page_calls():
    """The mock the Help page drew into, by running the page itself.

    Sized where the page unpacks: ``st.columns`` is called with 3 and with 2 on different
    rows, and ``st.tabs`` with five names, so a fixed-length return value would fail on the
    second caller rather than the first.

    Most claims about this page are about the strings it drew and take :func:`help_page`
    instead. This one exists for the claims that are about *how* something was drawn — the
    reference table reaches Streamlit as a frame plus a ``column_config``, and neither half
    is a string (issue #124).
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
        yield mock_st


@pytest.fixture
def help_page(help_page_calls) -> list:
    """Every string the Help page draws."""
    return _strings_drawn_to(help_page_calls)


@pytest.fixture
def page_text(help_page) -> str:
    """The page as one blob, for the claims that are about a phrase appearing at all."""
    return "\n".join(help_page)


@pytest.fixture
def page_prose(page_text) -> str:
    """The same blob with every run of whitespace collapsed to one space.

    For the tests that assert a sentence is **not** on the page. This prose is written as
    wrapped f-strings, so where a forbidden phrase falls across a wrap it reaches the page with
    a newline and the block's indentation inside it — and ``"would normally be" not in
    page_text`` then passes with the claim still on screen. Issue #150 hit the same wrap from
    the other side: its first positive assertion looked for a phrase that spanned one, and
    failed against a page that said exactly what it was asked to say.

    A negative assertion that can be defeated by reflowing a paragraph is not a guard on the
    claim, only on one typesetting of it — and reflowing is what editing prose does.
    """
    return " ".join(page_text.split())


# ---------------------------------------------------------------------------
# No such number
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("forbidden", ["500MB", "500 MB", "500Mb"])
def test_the_page_never_promises_a_500mb_upload(page_text, forbidden):
    """The number the app cannot honour, in any spelling.

    It stood in three places — the troubleshooting list, the data-preparation checklist and
    the FAQ — and nothing in this repository sets ``server.maxUploadSize``, so Streamlit's
    default of 200 MB is the real cap and the app refused a 300 MB MAF on a screen that had
    just invited it. Three copies of one wrong number is the argument for asserting the
    number rather than any one site: correcting two of three is the failure mode.
    """
    assert forbidden not in page_text, (
        f"the Help page says {forbidden!r} again. Nothing in this repo raises Streamlit's "
        "200 MB default — no .streamlit/config.toml, no flag in run_mafigate.sh or the "
        "Makefile — so a file above it is refused by the server (issue #79)."
    )


def test_the_page_states_the_limit_the_app_actually_enforces(page_text):
    """And says the true one, so the fix cannot be to fall silent instead."""
    assert "200 MB" in page_text


def test_the_page_never_tells_a_user_a_column_is_undescribed(page_text):
    """"Description not available for X" is a sentence, so it reads like an answer.

    That is why six columns went undescribed until the reference table was measured (issue
    #120), and `config.columns` keeps the fallback rather than raising precisely because
    this page calls `get_column_description` with lists of its own — `_KEY_CLINICAL_COLUMNS`
    and the core-required eight — where a `KeyError` would take the whole page down over one
    missing sentence.

    The guard in `test_column_resolver.py` covers the reference table, and it cannot cover
    these: the table is handed to `st.dataframe` as a frame, while these are drawn as
    markdown, and two of the names on this page's own shortlist — `gnomAD_exome_AF` and
    `gnomAD_genome_AF` — are deliberately *not* rows of that table. So the two guards are
    complementary rather than duplicated, and this one is the reason the fallback is allowed
    to stay a string here.

    Matched on the prefix taken from the constant, so renaming the fallback cannot leave a
    guard looking for wording nothing emits.
    """
    from config.columns import _MISSING_DESCRIPTION

    prefix = _MISSING_DESCRIPTION.split("{")[0]
    undescribed = [line for line in page_text.splitlines() if prefix in line]

    assert not undescribed, (
        f"the Help page tells the user a column has no description: {undescribed}"
    )


def test_the_reference_table_configures_only_columns_it_actually_has(help_page_calls):
    """A ``column_config`` key that matches no column is ignored in silence.

    So the failure mode of renaming a header is not an error — it is a column that quietly
    loses its width and its label and nothing says so. Issue #124 renamed one ("Required" →
    "In your MAF?"), which is what made this worth holding: both ends now read
    ``COLUMN_SOURCE_HEADER``, and this is what proves the frame and the config agree rather
    than merely being written from the same constant on the day.

    Asserted over *every* key, not just the renamed one — the next rename should fail here
    without anyone having to add a line for it.
    """
    from config.columns import COLUMN_SOURCE_HEADER

    configured = [
        call
        for call in help_page_calls.dataframe.call_args_list
        if "column_config" in call.kwargs
    ]
    assert configured, "the reference table is no longer drawn with a column_config"

    for call in configured:
        frame = call.args[0] if call.args else call.kwargs["data"]
        unmatched = sorted(set(call.kwargs["column_config"]) - set(frame.columns))
        assert not unmatched, (
            f"st.dataframe is configured for {unmatched}, which the frame does not have — "
            f"Streamlit drops those silently. The frame has {sorted(frame.columns)}"
        )

    assert any(
        COLUMN_SOURCE_HEADER in call.kwargs["column_config"] for call in configured
    ), f"nothing configures the {COLUMN_SOURCE_HEADER!r} column"


@pytest.mark.parametrize("format_name", ["TSV", "Excel", "xlsx"])
def test_the_page_offers_no_export_format_the_app_cannot_write(page_text, format_name):
    """Every download a user can reach writes CSV, plus the summary as text.

    The exporter that offered TSV and Excel had no caller anywhere in the app — the same
    false claim #55 removed from both READMEs and #71 from the About dialog — and #83 has
    since deleted it, so there is now no code in the app that could write either format.
    ``Excel`` is forbidden as a *format* here but the word may appear in prose ("opens in
    Excel, R and Python"), which is why this matches the offer-shaped spellings rather than
    the bare word.
    """
    offers = [
        line
        for line in page_text.splitlines()
        if format_name in line and ("**" + format_name in line or "as " + format_name in line)
    ]
    assert not offers, (
        f"the Help page offers {format_name} as an export format again: {offers}. Every "
        "reachable download is CSV; the only exporter that writes anything else has no "
        "caller (issue #79)."
    )


def test_the_page_describes_the_downloads_the_way_83_settled_them(page_text):
    """Two claims about the downloads that this page got wrong from opposite directions.

    "Timestamped filenames" was true when #79 checked it and false by the time #79 merged:
    #83 landed in between and made every download lead with the loaded MAF's name, because
    four files called ``passed_variants_all.csv`` are four files a clinician cannot tell
    apart. And "the columns on screen" is the wording **#83 deliberately removed** — the grid
    carries *Show all columns* and *Add columns*, so what is on screen is the user's to
    change while the narrow CSV is always the resolver's list. #79's own merge resolution
    reintroduced it and this test is why it did not ship.

    Asserted here rather than left to prose because it is the one claim on this page that a
    *sibling ticket* can falsify without touching Help — which is exactly the drift the map
    keeps finding.
    """
    assert "Timestamped filenames" not in page_text
    assert "columns on screen" not in page_text
    assert "name of the MAF you loaded" in page_text


def test_the_page_does_not_send_a_maf_through_excel(page_text):
    """The one Quick Fix that could destroy the file it claimed to repair.

    Excel truncates silently past 1,048,576 rows and rewrites gene symbols as dates —
    ``SEPT9``, ``MARCH1``. Advising a round-trip through it for clinical variant data was
    worse than saying nothing, so the advice is gone rather than qualified.
    """
    assert "Tab delimited" not in page_text
    assert "re-saving" not in page_text


# ---------------------------------------------------------------------------
# Derived, not transcribed
# ---------------------------------------------------------------------------


def test_the_core_required_columns_are_the_ones_the_app_refuses_without(help_page):
    """The list the page shows is ``REQUIRED_COLUMNS["core"]``, not a copy of it.

    ``validate_required_columns`` refuses a MAF missing any of these and names them, so a
    page listing a different set would send a user to check the wrong headers. The eight
    matched when #79 checked; this is what stops them coming apart.
    """
    core_section = [line for line in help_page if line.startswith("• **")]
    for column in REQUIRED_COLUMNS["core"]:
        assert any(
            line.startswith(f"• **{column}**") for line in core_section
        ), f"{column} is core-required but the Help page no longer lists it"


def test_the_frequency_columns_the_filter_reads_are_described_somewhere(page_text):
    """Including the two the column table cannot show.

    ``create_column_info_table`` is built from ``resolve_visible_columns``, so it holds the
    columns MAFigate *displays* — and ``gnomAD_exome_AF``/``gnomAD_genome_AF``, which the
    frequency filter reads first, are not among them, while their ``_raw`` variants are. Tab
    1 used to claim the table described "all columns used in MAFigate analysis", which sent a
    reader searching it for the very columns the next tab recommends. The table is right to be
    the resolver's set; the page has to cover the filter's inputs anyway.
    """
    from filters.variant_filters import FREQUENCY_COLUMNS

    for column in FREQUENCY_COLUMNS:
        assert column in page_text, (
            f"{column} is read by the population-frequency filter but named nowhere on the "
            "Help page (issue #79)"
        )


def test_every_escat_level_the_control_offers_is_named(page_text):
    """Rendered from ``ESCAT_OPTIONS``, so the page cannot offer a level the filter lacks.

    What this replaces was a hand-written meaning per level, gated on membership — which
    silently *dropped* any level the vocabulary gained and, worse, was wrong about two of
    the eight it did explain.
    """
    for level in ESCAT_OPTIONS:
        assert f"`{level}`" in page_text, f"ESCAT level {level} is selectable but unlisted"


def test_every_escat_level_the_control_offers_is_defined():
    """The keys are the options, in the options' order — the source's order (issue #89).

    #79 could name the levels and not define them, there being no authority in the checkout
    to define them from. Now there is one, and the risk inverts: a definition per level is
    exactly the shape that drifted before, so the two lists are pinned to each other rather
    than maintained in parallel.
    """
    assert list(ESCAT_DEFINITIONS) == ESCAT_OPTIONS, (
        "ESCAT_DEFINITIONS and ESCAT_OPTIONS have come apart; a level the control offers "
        "with no definition is what #79 found, and a definition with no level is dead prose"
    )


def test_every_escat_definition_reaches_the_page(page_text):
    """Both halves of each level, drawn — the evidence that earns it and what it implies.

    The implication is the half a clinician reads for, and the half easiest to drop while
    still looking complete: eight levels, each with a plausible sentence, is what the page
    had before.
    """
    for level, meaning in ESCAT_DEFINITIONS.items():
        assert meaning.evidence in page_text, f"ESCAT {level} is listed without its evidence"
        assert meaning.implication in page_text, (
            f"ESCAT {level} says what earns it but not what it implies"
        )


def test_the_page_says_where_the_escat_definitions_come_from(page_text):
    """A definition the app cannot cite is the one #79 deleted.

    The page carries the paper and its DOI, so a reader who doubts a level can check it —
    which is the whole difference between this block and the prose it replaces.
    """
    assert ESCAT_SOURCE in page_text, "the ESCAT definitions are rendered without their source"
    assert ESCAT_SOURCE_DOI in page_text


def test_every_civic_level_the_control_offers_is_defined():
    """The same pinning for CIViC's scale, and for the same reason (issue #109).

    The keys are the options in the options' order, which is CIViC's own — strongest evidence
    first — so a level cannot be offered undefined or defined unoffered. ESCAT's pair was found
    to be in its paper's order *by accident* and asserted nowhere (#89); this one is asserted.
    """
    assert list(CIVIC_DEFINITIONS) == CIVIC_OPTIONS, (
        "CIVIC_DEFINITIONS and CIVIC_OPTIONS have come apart; a level the control offers with "
        "no definition is unexplained on the page, and one without a level is dead prose"
    )


def test_every_civic_definition_reaches_the_page(page_text):
    """Both halves of each level: CIViC's name for it, and the study type that earns it.

    This is what replaces the transcribed dict that stood in the page — and the guard it never
    had. That dict glossed ``A`` as "Validated - Strong clinical significance", which names the
    *other* axis: clinical significance is what ClinVar reports, while a CIViC level says how
    strong the evidence behind an assertion is. Read off the constant, so the page cannot drift
    from what the circles are graded on.
    """
    for level, meaning in CIVIC_DEFINITIONS.items():
        assert meaning.name in page_text, f"CIViC {level} is listed without CIViC's name for it"
        assert meaning.evidence in page_text, (
            f"CIViC {level} is named without the study type that earns it"
        )


def test_the_page_says_where_the_civic_definitions_come_from(page_text):
    """A definition the app cannot cite is one it should not make — #89's line, second scale.

    The page carries CIViC's documentation and the URL, so a reader who doubts a level can check
    it. The dict this replaces cited nothing, which is why "roughly right" was as good as it got.
    """
    assert CIVIC_SOURCE in page_text, "the CIViC definitions are rendered without their source"
    assert CIVIC_SOURCE_URL in page_text


def test_the_page_grades_civic_evidence_rather_than_the_variant(page_text):
    """The one claim the whole change rests on, said where a user reads the levels.

    Not a style point: the app draws a coloured circle per CIViC level beside five circles that
    *are* pathogenicity calls, so a page that lists A-to-E without saying what the scale grades
    leaves the reader to assume it is the same axis — which is precisely the assumption the code
    used to make, mapping ``E`` to Benign.
    """
    assert "grades the evidence, not the variant" in page_text, (
        "the CIViC block lists the levels without saying the scale grades the evidence; the "
        "circles beside it are graded on this scale and the key says so"
    )


def test_every_clinvar_term_the_control_offers_is_named(page_text):
    """Every term the control offers, not the five the old expander had a gloss for.

    The gate ran the other way round: it walked five hand-written meanings and showed a term
    only if the vocabulary also had it, so the selectable values without one — among them
    ``Conflicting_classifications_of_pathogenicity``, ``drug_response`` and
    ``not_provided`` — were absent from a list presenting itself as the ClinVar options.

    The counts this docstring used to carry ("all eleven", "six selectable values") are
    gone: #88 took the two dead composite terms out of ``CLINVAR_OPTIONS``, and a number
    written beside a list it does not read is the drift this file exists to catch.
    """
    for term in CLINVAR_OPTIONS:
        assert term in page_text, f"ClinVar term {term} is selectable but unlisted"


#: The offered ClinVar terms whose *names* describe a record rather than a finding.
#:
#: A clinician can read ``Likely_benign`` off its name. These three cannot be read off theirs
#: at all — ``-`` is a bare dash in a list — so being *listed* is not being explained, which is
#: the distinction issue #103 was asked about: the ticket's words were that one of these
#: "is not a clinical distinction a dropdown can carry unexplained".
#:
#: Three and not more because each is quoted from a source this repository holds, both cited in
#: ``components/clinical_summary.py`` where #98 recorded them: the institute's term table and
#: ClinVar's own definition of ``other``. ``no_classifications_from_unflagged_records`` is
#: offered and deliberately unglossed — nothing here defines it, and #79 deleted invented
#: clinical prose rather than write it from memory.
_RECORD_STATE_TERMS = ("-", "other", "no_classification_for_the_single_variant")


@pytest.mark.parametrize("term", _RECORD_STATE_TERMS)
def test_the_record_state_clinvar_terms_carry_a_gloss(help_page, term):
    """Listed is not explained, for the terms whose names say nothing (issue #103).

    #103 tripled this list, and most of what it added is readable: ``risk_factor`` and
    ``protective`` say what they are. These three do not, and one of them is a dash.

    Asserted on the *line* rather than on the page blob, because every one of these terms is
    guaranteed to appear in the page text anyway — it is in the list — so a page-contains
    check would pass with the gloss deleted. The line is the claim: the term, then a colon,
    then words. That is also what makes this fail when the gloss goes, which a check for the
    gloss's own wording would do too but would re-break on any rewording.
    """
    prefix = f"• **{term}**"
    lines = [line for line in help_page if line.startswith(prefix)]
    assert lines, (
        f"the Help page does not list {term!r} as a ClinVar term at all; it is offered by the "
        "control, so it must appear in this list"
    )

    for line in lines:
        remainder = line[len(prefix):]
        assert remainder.startswith(":") and remainder[1:].strip(), (
            f"the Help page lists {term!r} with no explanation — the line drawn was "
            f"{line!r}. Its name says nothing a clinician can read, so listing it is not "
            "explaining it. Gloss it from a source this repo holds, as `clinvar_meanings` "
            "does, or say here why it cannot be glossed"
        )


def test_the_two_conflicting_spellings_are_explained_as_one_call(page_text):
    """The page lists both, so it has to say why (issues #88, #79).

    ClinVar renamed ``Conflicting_interpretations_of_pathogenicity`` to
    ``Conflicting_classifications_of_pathogenicity`` in 2023, and this repo pins no ClinVar
    release, so both are live spellings depending on when a file was annotated — #88 put
    both on offer for that reason, and the test above then requires this page to list both.

    Which creates a duplicate on a clinician's screen: two lines differing by one word,
    reading as two different classifications. The parameter page's tooltip explains it and
    this page did not, so the page was the less honest of the two surfaces about a control
    it exists to document.

    Said in the caption rather than as per-term glosses, deliberately: a gloss here is a
    *clinical* claim, and the ESCAT levels are what happens when this page makes one from
    memory. Naming the rename is a fact about ClinVar's vocabulary, which is checkable.
    """
    assert "renamed" in page_text.lower(), (
        "the Help page lists both spellings of the conflicting-classification call and "
        "never says one is a rename of the other, so they read as two classifications"
    )
    assert "2023" in page_text, (
        "the rename is unattributed to a date, which is the only thing that tells a user "
        "which spelling their own file should carry"
    )


def test_the_underscored_spellings_are_explained_as_the_same_terms(page_text):
    """The same obligation as the rename, for the pairs #99 put on offer.

    ``other`` and ``_other`` are one ClinVar call written two ways — a modifier following a
    call keeps its leading underscore in some files and loses it in others — so the test above
    requires this page to list both, and a clinician reading two lines that differ by a
    punctuation mark will take one of them for a typo and unselect it. On the two files in the
    measured corpus that spell it with the underscore, unselecting ``_other`` loses **every**
    variant carrying that call.

    Checked by what renders, and against the vocabulary rather than against a hardcoded pair:
    the caption reads its members out of ``CLINVAR_OPTIONS`` precisely so that #103 cannot add
    a term and leave this page naming the wrong ones, and a test spelling them out again would
    be the second place to keep in step that the caption avoids being.
    """
    underscored = [term for term in CLINVAR_OPTIONS if term.startswith("_")]
    assert underscored, (
        "no underscore-prefixed spelling is offered any more, so this test is asserting "
        "nothing — delete it with the caption it guards, or restore the terms"
    )

    for term in underscored:
        assert term in page_text, f"the underscored spelling {term} is offered but unlisted"

    assert "underscore" in page_text.lower(), (
        "the page lists spellings differing only by a leading underscore and never says they "
        "are the same terms, so each pair reads as two classifications and a user drops one"
    )
    assert "not separate classifications" in page_text, (
        "nothing on the page states the one fact that stops a clinician unselecting half of "
        "each pair: the two spellings are not different classifications"
    )


# ---------------------------------------------------------------------------
# One story per fact
# ---------------------------------------------------------------------------


def test_escat_tiers_are_not_described_as_resistance_mechanisms(page_text):
    """ESCAT has no resistance tier; ``IIIA``/``IIIB`` were described as one.

    Resistance levels are OncoKB's ``R1``/``R2``, which this app does not read — Table 1 of
    the ESCAT paper attributes resistance-biomarker grading to OncoKB and to nothing else.
    #79 deleted the gloss for want of an authority; #89 supplied one, so the check now holds
    the definitions as well as the page: the word must not come back through either door.
    """
    assert "Resistance mechanism" not in page_text
    assert "esistance" not in " ".join(
        meaning.evidence + meaning.implication for meaning in ESCAT_DEFINITIONS.values()
    ), "an ESCAT definition describes resistance; the scale has no resistance level"


def test_no_escat_definition_repeats_a_claim_that_lost(page_text):
    """The two stories about level ``V`` that this page and the tooltip used to tell.

    ``V`` read "Not actionable" here and "case reports" in the tooltip; ESCAT's ``V`` is
    neither — it is an alteration whose matched drug produces objective responses without
    improving outcome. "Lack of evidence for actionability" is the scale's ``X``, which this
    control does not offer at all, so "not actionable" was a gloss of the wrong level.

    Two of three, not three: the Pathogenicity Overview still maps ``V`` to ``Unknown`` by
    string prefix, and issue #100 owns it. Said here because a check that quietly covers two
    surfaces while its ticket named three is how "covered elsewhere" becomes false.

    **The "case reports" half is now measured against CIViC's own definition** (issue #109).
    It was a page-wide search for the phrase, which was sound while ESCAT was the only level
    scale this tab described — and CIViC's ``C`` is *defined* by CIViC as "Individual case
    reports from clinical journals", so the page now says it truthfully, once.

    So the page-wide reach is kept and the *budget* is derived: every occurrence of the phrase
    must be one a CIViC definition accounts for. A hand-written "case reports" returning
    anywhere — including the ESCAT expander's own intro and caption, which are prose and not
    rendered from :data:`ESCAT_DEFINITIONS` — pushes the count past what the constant explains
    and fails. Asserting only over the definitions would have left exactly those two strings
    uncovered, which is the half of this page ``V``'s losing gloss lived in.
    """
    assert "Not actionable" not in page_text

    phrase = "case report"
    explained = sum(
        meaning.evidence.lower().count(phrase) for meaning in CIVIC_DEFINITIONS.values()
    )
    assert page_text.lower().count(phrase) == explained, (
        "the page mentions case reports more often than CIViC's own level definitions "
        "explain — an ESCAT level is glossed as case reports again, or some other surface "
        "has grown a second copy of CIViC's"
    )
    assert phrase not in " ".join(
        meaning.evidence + meaning.implication for meaning in ESCAT_DEFINITIONS.values()
    ).lower(), "an ESCAT definition glosses a level as case reports; that was V's losing story"


def test_the_parameter_tooltip_defers_to_the_help_page_about_escat():
    """One definition per level, in one place, and the tooltip is not a second copy.

    The tooltip glossed the scale as "IA=strongest evidence to V=case reports" while Help
    said "V: Not actionable" — a contradiction that existed *because* two surfaces each held
    their own copy, so a third copy is how it recurs (issue #89). The tooltip therefore states
    the direction only, reads the strongest level off ``ESCAT_STRONGEST`` rather than spelling
    it out, and points at the page that defines them.
    """
    from page_modules.parameter_config import _GUIDELINE_CONTROLS

    (escat,) = [
        control for control in _GUIDELINE_CONTROLS["somatic"] if control[0] == "filter_escat"
    ]
    tooltip = escat[3]
    assert "case reports" not in tooltip
    assert "Not actionable" not in tooltip
    assert ESCAT_STRONGEST in tooltip
    for level, meaning in ESCAT_DEFINITIONS.items():
        assert meaning.evidence not in tooltip, (
            f"the tooltip has acquired its own copy of ESCAT {level}'s definition; that is "
            "the two-copy shape #79 found contradicting itself"
        )


def test_the_faq_does_not_say_a_missing_column_switches_a_filter_off(page_text):
    """The FAQ's pre-#39 answer, which its own Required Columns tab contradicted.

    "MAFigate gracefully handles missing columns by skipping related filters" is the
    opposite of what the app does — it fills the column with a stand-in and reports that the
    result is incomplete — and it is wrong in the expensive direction: filling ``CancerVar``
    neutrally drops about 70% of the somatic report. Two accounts of one behaviour on one
    page, and the reassuring one was false.
    """
    assert "skipping related filters" not in page_text
    assert "gracefully handles missing columns" not in page_text


@pytest.mark.parametrize(
    "forbidden",
    ["would normally be", "refused outright", "would normally be refused"],
)
def test_the_page_does_not_say_a_missing_filter_input_would_be_refused(page_prose, forbidden):
    """A sentence issue #136 deleted from the warnings, still drawn here until issue #150.

    The block read *"A MAF missing a filter input would normally be refused outright; MAFigate
    fills the column with a stand-in value"* — and **MAFigate never refuses on an absent
    column**. Absence fills and warns; only an unreadable *value* refuses, in
    ``filters/numeric_columns.py``. So "normally" could only have meant *in the pipeline*,
    which is the comparison decision 2 of the map retired.

    Asserted as a phrase rather than at a site, for the reason the 500MB test gives: this claim
    has now been found in three places — the two fill warnings (#136), ``_gene_notes``
    (#150) and this page (#150) — and correcting all but one is the demonstrated failure mode.

    Deliberately **not** ``"refused" not in page_prose``, which the ``filters/`` guard can
    afford and this page cannot: it says, truthfully, that MAFigate refuses a file without the
    core required columns, and ``validate_required_columns`` is what makes that true.

    Read off :func:`page_prose` rather than :func:`page_text`, and the third parameter is the
    whole sentence: in the deleted copy the wrap fell between *"would normally be"* and
    *"refused outright"*, so each half was findable by luck of the typesetting and the claim
    they make together was not findable at all.
    """
    assert forbidden not in page_prose, (
        f"the Help page says {forbidden!r} again. MAFigate does not refuse a file for a "
        "missing filter input — it fills the column and says the report is incomplete "
        "(issues #136, #150)."
    )


def test_the_page_still_says_what_a_missing_filter_input_costs(page_prose):
    """So falling silent cannot be the fix for the test above.

    The honest half of that sentence was always the second half, and it is the sentence
    ``filters.absent_columns._not_a_complete_result`` spells for the warnings.

    Pinned on the clause **only this block draws**, not on "stand-in value" alone: that phrase
    reaches the page five times, and the FAQ gives its own account of the same behaviour a
    thousand lines down — *"MAFigate fills the missing column with a stand-in value"*. A guard
    satisfied by either copy would let the Handling Missing Columns block fall silent and still
    pass, which is the shape of vacuity this suite keeps finding.

    Long enough to be that clause and no other, which :func:`page_prose` is what makes possible:
    the phrase spans a wrap in the source, and read off :func:`page_text` it failed against a
    page that was saying exactly what it had been asked to say.
    """
    assert "fills a missing filter input with a stand-in value" in page_prose
    assert "not a complete result" in page_prose


def test_the_page_does_not_call_pathogenic_retention_unconditional(page_prose):
    """The word #136 took out of the escalated warning for over-claiming.

    The rescue clears the *filter criteria*, but ``apply_filters`` masks it with the app's own
    population-frequency term — ``(criteria | rescue) & (frequency_ok | pathogenic)`` — and the
    only exemption is the narrower ClinVar-only :func:`~filters.variant_filters.
    pathogenic_exemption`. So a retained variant can still be dropped for being common, and
    "unconditional" promises a safety net with no holes in it.

    Found by issue #150 in the same three lines as the refusal claim above: the column *lists*
    in this block are derived from the vendored source and so cannot drift, which is exactly
    why the hand-written sentences beside them were the ones that had.
    """
    assert "unconditional" not in page_prose


def test_the_page_does_not_claim_the_app_runs_in_the_browser(page_text):
    """The page's most consequential false claim, on an app whose input is patient data.

    MAFigate is a server application reached *through* a browser: the file is uploaded to
    the Python process, spilled to a temporary file to be read, and held in memory for the
    session. Whether that is "local" depends on where the process runs, which is a property
    of how somebody launched it and not of this app — so neither claim is one the page can
    make. Both assertions survived issue #182 binding ``run_mafigate.sh`` to loopback: a
    server on ``127.0.0.1`` is still a server you reach through a browser, and the page's own
    answer already handles the shared-server case in the reader's terms.
    """
    assert "runs in your browser" not in page_text
    assert "processed locally" not in page_text


def test_the_page_promises_no_filter_editing_it_does_not_offer(page_text):
    """A whole feature that never existed.

    "Quick Parameter Adjustment: Modify key filters without leaving the page" described no
    control in the app: the Load section offers Re-apply Filters and a button that *leaves*
    for Configure Parameters, and there is nowhere else a threshold can be touched.
    """
    assert "without leaving the page" not in page_text
    assert "Quick Parameter Adjustment" not in page_text


def test_the_page_claims_no_comparison_feature(page_text):
    """"Compare different filtering strategies" named a feature the app does not have."""
    assert "Compare different filtering strategies" not in page_text


# ---------------------------------------------------------------------------
# The page's own controls
# ---------------------------------------------------------------------------


def test_the_column_search_survives_a_regex_metacharacter():
    """Typing ``(`` into the column search used to take the whole page down.

    ``str.contains`` treats its argument as a pattern by default, so ``(`` raised
    ``re.error`` out of the search, ``route_to_page`` caught it, and the user landed on Home
    under "Error loading page 'help'". Asserted through the page rather than against
    ``str.contains`` directly, because the fix is an argument the page passes.
    """
    from page_modules import help as help_module

    with patch("page_modules.help.st") as mock_st:
        mock_st.session_state = FakeSessionState()
        mock_st.columns.side_effect = lambda spec, **kw: [
            MagicMock() for _ in range(spec if isinstance(spec, int) else len(spec))
        ]
        mock_st.tabs.return_value = [MagicMock() for _ in range(5)]
        mock_st.button.return_value = False
        mock_st.selectbox.return_value = "All"
        mock_st.text_input.return_value = "("

        help_module.show_help_page()


def test_the_required_columns_tab_claims_an_arm_only_where_there_is_one(page_text):
    """The tab's per-arm blocks used to claim an arm for columns that have none (issue #127).

    `Tumor_Seq_Allele1`/`2` sat under **Somatic** and `n_alt_count`/`n_ref_count` under
    **Germline** while the reference table one tab over reported all four as **Both**, so the
    two surfaces disagreed on one screen. `tumor_f` and the tumour read counts were in the
    somatic block on the same false footing: both arms' filters read them.

    The block now says "annotations only a {arm} report **reads**", and that is the claim
    checked here — two ways, both sourced from the derivation over the vendored filters rather
    than from the set-difference the page itself computes, which would agree with the page by
    construction:

    * nothing both arms read may appear in an arm block (`tumor_f`, the tumour read counts);
    * nothing *no* filter reads may appear either (the allele columns, the normal read counts).

    **Not** asserted: that every column in a block is reported per-arm by the reference table.
    `ESCAT` is emitted on both arms and read by the somatic filter alone, and those are
    different questions — which arms carry the column, and which arm's filter consults it. An
    earlier draft of this test conflated them and failed on `ESCAT`; the tab is right there.

    Mutation: put `n_alt_count` back into the germline block → the second assertion names it.
    """
    from filters.absent_columns import PATHOGENIC_INPUTS, REQUIRED_INPUTS

    reads = {
        arm: set(REQUIRED_INPUTS[arm]) | set(PATHOGENIC_INPUTS[arm])
        for arm in ("somatic", "germline")
    }
    from page_modules.help import _DISPLAY_ONLY

    read_by_both = reads["somatic"] & reads["germline"]
    read_by_neither = _DISPLAY_ONLY
    assert read_by_both and read_by_neither, "this guard has nothing to check"

    # The arm blocks are the lines between the somatic heading and the leftover-columns block
    # after them; taking the whole page would sweep in the reference table's own rows.
    #
    # The end sentinel is read from the constant rather than written out (issue #219, which
    # renamed that heading from "Clinical Annotation Columns"). A literal here failed in the
    # worst available way: `next()` with no default raises `StopIteration`, so a renamed heading
    # was reported as an error in this guard rather than as a false claim about the page — and
    # the assertions below would never have run to say so. The `assert` keeps that distinction
    # if the block is ever deleted outright.
    from page_modules.help import _OTHER_COLUMNS_HEADING

    lines = page_text.splitlines()
    start = next(i for i, line in enumerate(lines) if "Somatic Analysis" in line)
    ends = [
        i
        for i, line in enumerate(lines[start:], start)
        if _OTHER_COLUMNS_HEADING in line
    ]
    assert ends, (
        f"the tab no longer renders {_OTHER_COLUMNS_HEADING!r}, so this guard cannot tell "
        "where the per-arm blocks stop and is checking nothing"
    )
    arm_prose = "\n".join(lines[start : ends[0]])

    shared = sorted(column for column in read_by_both if f"**{column}**" in arm_prose)
    assert not shared, (
        f"the tab lists {shared} as annotations only one arm reads, and both arms' filters "
        "read every one of them"
    )

    unread = sorted(column for column in read_by_neither if f"**{column}**" in arm_prose)
    assert not unread, (
        f"the tab lists {unread} under an arm, and no filter on either arm reads them — the "
        "reference table reports them as Both, which is the disagreement this fixed"
    )


# ---------------------------------------------------------------------------
# Issue #219 — the category words, and the surface they point at
# ---------------------------------------------------------------------------


def _tab_text(tab_function_name: str) -> str:
    """Every string one Help tab draws, by rendering just that tab.

    The whole-page fixtures cannot separate the tabs — they flatten five tab bodies into one
    blob — and two of the claims below are about *which* tab something is on. Issue #219's
    finding is exactly that kind: the FAQ sent a reader to Column Information for the
    classification systems, and every one of them renders in Filter Options.
    """
    from page_modules import help as help_module

    with patch("page_modules.help.st") as mock_st:
        mock_st.session_state = FakeSessionState()
        mock_st.columns.side_effect = lambda spec, **kw: [
            MagicMock() for _ in range(spec if isinstance(spec, int) else len(spec))
        ]
        mock_st.button.return_value = False
        mock_st.text_input.return_value = ""
        mock_st.selectbox.return_value = "All"
        getattr(help_module, tab_function_name)()
        return "\n".join(_strings_drawn_to(mock_st))


def test_the_page_calls_none_of_these_a_clinical_database(page_prose):
    """*Database* is false of most of what MAFigate reads (issue #219).

    Map #199 established it three times over: ``RENOVO`` is an in-silico predictor this pipeline
    ran (issue #201 — passing *computed here* is exactly the leg a third-party database fails),
    ``InterVar`` and ``CancerVar`` are guideline verdicts computed here (issue #187), and
    AlphaMissense is a predictor on its own scale (issue #203). ClinVar and CIViC are the only
    two the word fits, and neither is ever described collectively.

    The phrase, not the word: ``COLUMN_DESCRIPTIONS`` says *ClinVar database*, *CIViC database*
    and *COSMIC database* on the columns where that is simply true, and this page renders every
    one of them. A guard on ``"database"`` would forbid the accurate uses along with the
    inaccurate one, and would have to be weakened the first time it fired.

    Two places said it and both are fixed: the FAQ's interpretation answer, and the Filter
    Options tab's own lead sentence — which is why this is asserted over the page rather than
    over one tab. Map #199 rule 5 protects what the filters *read* and how that is documented;
    the category word naming a group of controls is neither.

    Mutation: restore either sentence's *database* → this names it.
    """
    assert "clinical database" not in page_prose.lower(), (
        "the Help page calls these sources clinical databases, and most of them are not — "
        "RENOVO is a predictor this pipeline ran, InterVar and CancerVar are guideline "
        "verdicts computed here"
    )


def test_the_faq_sends_the_reader_to_the_tab_that_holds_the_scales():
    """The classification systems are on Filter Options, and were never on Column Information.

    Issue #219's second finding, and the one the ticket did not anticipate: the FAQ answered
    *"How do I interpret the clinical annotations?"* with *"refer to the Column Information tab
    for detailed descriptions of each clinical database and their classification systems"*, and
    that tab holds no classification system at all. It is one line per column from
    ``COLUMN_DESCRIPTIONS``; every scale renders in ``show_filter_options_tab``.

    Asserted as a *contrast* rather than as a phrase, so it cannot be satisfied by naming the
    right tab in a sentence that has stopped being true. The probes are ESCAT's and CIViC's own
    glossed levels, which only the scale expanders draw.

    Mutation: point the answer back at Column Information → the first assertion names it.
    Mutation: move the ESCAT expander onto the Column Information tab → the third fires.
    """
    faq = _tab_text("show_faq_tab")
    interpretation = [
        line for line in faq.splitlines() if "interpret the clinical annotations" in line
    ]
    assert interpretation, "the FAQ no longer asks how to interpret the clinical annotations"

    answer = faq.split("interpret the clinical annotations", 1)[1]
    answer = answer.split("**Q:", 1)[0]
    assert "Filter Options" in answer, (
        f"the interpretation answer does not send the reader to the tab that holds the "
        f"scales: {answer.strip()!r}"
    )

    scales = _tab_text("show_filter_options_tab")
    columns_tab = _tab_text("show_column_information_tab")
    for probe, name in ((ESCAT_SOURCE_DOI, "the ESCAT citation"), ("Level A", "CIViC's levels")):
        assert probe in scales, f"{name} is not on the Filter Options tab any more"
        assert probe not in columns_tab, (
            f"{name} now renders on the Column Information tab too — the FAQ names one tab, "
            "and a second copy is the drift this map keeps finding"
        )


def test_the_clinical_shortlist_names_every_column_the_panel_draws():
    """The shortlist claims to hold what decides a clinical reading; #204 made it false of five.

    ``components/clinical_badges.py`` reads seven columns and ``components/alphamissense.py``
    one, and the variant panel draws all eight — ``ClinVar_VCF_CLNREVSTAT`` as stars,
    ``RENOVO_pls`` inline, ``ESCAT_TISSUE``/``ESCAT_CANCER`` inside the tier,
    ``am_pathogenicity`` as the band. Four of those eight were not on this page's shortlist
    (issue #219), and two of them are not rows of the reference table either, so a reader who
    saw the stars could look them up nowhere on the page.

    The panel's columns are read out of its **source** rather than listed here, so a ninth
    column the panel starts drawing fails this without anyone remembering to add it. The AST
    walk is over ``_text(row, "...")``, which is the one accessor ``clinical_badges`` reads a
    cell through.

    Mutation: drop ``RENOVO_pls`` from ``_KEY_CLINICAL_COLUMNS`` → this names it.
    Mutation: have a badge read a new column through ``_text`` → this names that one.
    """
    import ast
    import pathlib

    from components.alphamissense import AM_COLUMN
    from page_modules.help import _KEY_CLINICAL_COLUMNS

    source = pathlib.Path("components/clinical_badges.py").read_text()
    drawn = set()
    for node in ast.walk(ast.parse(source)):
        if (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Name)
            and node.func.id == "_text"
            and len(node.args) == 2
            and isinstance(node.args[1], ast.Constant)
            and isinstance(node.args[1].value, str)
        ):
            drawn.add(node.args[1].value)

    assert len(drawn) >= 7, (
        f"the AST walk found only {sorted(drawn)} — `clinical_badges` reads cells some other "
        "way now, so this guard is measuring the wrong thing"
    )
    drawn.add(AM_COLUMN)

    missing = sorted(drawn - set(_KEY_CLINICAL_COLUMNS))
    assert not missing, (
        f"the variant panel draws {missing} and the Key Clinical Columns shortlist does not "
        "name them, so its heading is false of them"
    )


def test_every_column_the_panel_draws_is_described(page_text):
    """And described by name, so the shortlist above cannot render the fallback sentence.

    ``ClinVar_VCF_CLNREVSTAT`` had **no** ``COLUMN_DESCRIPTIONS`` entry at all until issue #219,
    and it is not in ``resolve_visible_columns`` either — so the reference table's own
    completeness guard could not reach it, and adding it to the shortlist would have drawn
    *"Description not available"* on the page. The two guards are complementary in exactly the
    way ``test_the_page_never_tells_a_user_a_column_is_undescribed`` records.

    Each description must also name its tool, because *"Pathogenicity prediction from annotation
    tools"* is what ``am_pathogenicity`` said while the panel above it drew a heading reading
    *AlphaMissense (in silico missense predictor)* — a description true of nothing in particular
    is how a reader fails to connect the two.

    Mutation: delete the ``ClinVar_VCF_CLNREVSTAT`` description → the first assertion fires.
    Mutation: restore ``am_pathogenicity``'s old wording → the second names it.
    """
    from config.columns import COLUMN_DESCRIPTIONS, _MISSING_DESCRIPTION
    from page_modules.help import _KEY_CLINICAL_COLUMNS

    undescribed = [col for col in _KEY_CLINICAL_COLUMNS if col not in COLUMN_DESCRIPTIONS]
    assert not undescribed, (
        f"{undescribed} are on the Key Clinical Columns shortlist with no description, so the "
        f"page renders {_MISSING_DESCRIPTION.split('{')[0]!r} where the answer belongs"
    )

    for column, tool in (
        ("am_class", "AlphaMissense"),
        ("am_pathogenicity", "AlphaMissense"),
        ("ClinVar_VCF_CLNREVSTAT", "ClinVar"),
        ("RENOVO_pls", "RENOVO"),
    ):
        assert tool in COLUMN_DESCRIPTIONS[column], (
            f"{column}'s description does not name {tool}, so a reader who saw it on the "
            f"variant panel cannot match the two: {COLUMN_DESCRIPTIONS[column]!r}"
        )

    assert "AlphaMissense_score" in COLUMN_DESCRIPTIONS["am_pathogenicity"], (
        "am_pathogenicity's description does not warn that AlphaMissense_score is a different "
        "annotation — issue #203 measured 346 rows across 84 files falling in different classes"
    )
    assert _MISSING_DESCRIPTION.split("{")[0] not in page_text


def test_the_leftover_columns_block_duplicates_no_list_beside_it():
    """What survived the #219 deletion is the columns no other list on the tab carries.

    Ten of the block's twelve entries were a second copy: ``ESCAT``, ``CancerVar`` and
    ``CIViC_Evidence_Level`` are rendered by the per-arm derivation five lines above, and all
    seven ``FREQUENCY_COLUMNS`` fifty lines below. A duplicated list is free to drift from the
    one beside it, which is the failure issue #126 found in this very block.

    Held against the derivations rather than against the surviving names, so re-adding any of
    the ten fails here whichever bucket it is put in.

    Mutation: put ``list(FREQUENCY_COLUMNS)`` back into the block → this names every one.
    """
    import re

    from filters.absent_columns import PATHOGENIC_INPUTS, REQUIRED_INPUTS
    from filters.variant_filters import FREQUENCY_COLUMNS
    from page_modules.help import _OTHER_COLUMNS_HEADING

    tab = _tab_text("show_required_columns_tab")
    assert _OTHER_COLUMNS_HEADING in tab, "the leftover-columns block no longer renders"
    block = tab.split(_OTHER_COLUMNS_HEADING, 1)[1].split("Handling Missing Columns", 1)[0]

    reads = {
        arm: set(REQUIRED_INPUTS[arm]) | set(PATHOGENIC_INPUTS[arm])
        for arm in ("somatic", "germline")
    }
    arm_only = (reads["somatic"] ^ reads["germline"]) | set(FREQUENCY_COLUMNS)
    repeated = sorted(column for column in arm_only if f"**{column}**" in block)
    assert not repeated, (
        f"{repeated} are listed in the leftover block and again elsewhere on the same tab — "
        "the duplication issue #219 deleted"
    )

    named = set(re.findall(r"\*\*([A-Za-z0-9_.]+)\*\*", block))
    assert named, f"the block names no column at all: {block!r}"


def test_the_frequency_list_is_derived_rather_than_typed_out(page_text):
    """Both copies of it were hand-written or derived, and the derived one is the survivor.

    Issue #126 replaced a hand-written frequency list with ``FREQUENCY_COLUMNS`` after finding
    it had omitted ``Freq_esp6500siv2_all`` and named ``MAX_AF``, which the filter does not
    read. Issue #219 deleted that derived copy as a duplicate — and the list it duplicated was
    the *hand-written* one, so deleting it alone would have undone #126. The remaining list is
    derived instead.

    **Derivation is a claim about the source, not about the page**, and this guard was wrong
    about that once. Typing today's seven names back into the ``st.info`` block renders a list
    byte-identical to the helper's, so no assertion over the rendered text can tell the two
    apart — the copy only becomes false on the day ``FREQUENCY_COLUMNS`` changes, which is
    precisely the drift issue #126 found after the fact. The call site is asserted by walking the
    AST instead: the tab must *call* the helper. Measured, not assumed — the rendered-text
    version of this assertion was written first and the mutation went uncaught.

    Mutation: replace ``{_frequency_groups()}`` with today's seven names → the AST assertion
    fires, and nothing else does, which is the whole point of it being here.
    Mutation: drop the ``remaining`` branch from ``_frequency_groups`` → the unknown-source
    assertion fires.
    """
    import ast
    import pathlib

    from filters.variant_filters import FREQUENCY_COLUMNS
    from page_modules import help as help_module
    from page_modules.help import _frequency_groups

    tab = next(
        node
        for node in ast.walk(ast.parse(pathlib.Path("page_modules/help.py").read_text()))
        if isinstance(node, ast.FunctionDef) and node.name == "show_required_columns_tab"
    )
    calls = {
        node.func.id
        for node in ast.walk(tab)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
    }
    assert "_frequency_groups" in calls, (
        "show_required_columns_tab no longer calls _frequency_groups, so its frequency list is "
        "transcribed again — and a transcription that is correct today is exactly what issue "
        "#126 found had drifted"
    )

    grouped = _frequency_groups()
    for column in FREQUENCY_COLUMNS:
        assert f"`{column}`" in grouped, (
            f"{column} is read by the frequency filter and _frequency_groups does not name it"
        )
        assert f"`{column}`" in page_text, f"{column} does not reach the page"

    assert grouped in page_text, (
        "the page's frequency list is not the one _frequency_groups builds, so it is a "
        "hand-written copy again (issue #126)"
    )

    # A column matching none of the three source prefixes is listed rather than dropped. This is
    # the half of #126 a derived list does not fix on its own: silently omitting a column the
    # filter reads is exactly what the hand-written copy did.
    with patch.object(help_module, "FREQUENCY_COLUMNS", list(FREQUENCY_COLUMNS) + ["AF_new"]):
        assert "`AF_new`" in _frequency_groups(), (
            "a frequency column matching no known source is dropped from the page in silence"
        )


def test_the_pandas_fallback_is_gone_rather_than_unreachable():
    """It could never fire, and it rendered a second copy of the same table if it did.

    ``page_modules/help.py`` imported pandas at module scope and ``config.columns`` imports
    it too, so a checkout without pandas failed the page long before the ``except
    ImportError`` that claimed to handle it — which made ``show_basic_column_info``, sixty
    lines of hand-grouped column categories, unreachable code maintaining a second column
    vocabulary.
    """
    from page_modules import help as help_module

    assert not hasattr(help_module, "show_basic_column_info")
    assert not hasattr(help_module, "pd"), (
        "page_modules/help.py imports pandas again. Nothing in the file uses it; it existed "
        "only for the fallback that could not fire (issue #79)."
    )
