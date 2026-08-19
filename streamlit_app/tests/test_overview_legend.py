"""What the Pathogenicity Overview draws, and what its key says it drew (issue #95).

The column reduces each variant's annotation sources to one coloured circle apiece, and a
key above the grid explains the glyphs and names the sources in circle order. The key used
to name **all six** sources on both arms, from a hand-written parenthetical beside
abbreviations it read from :data:`components.clinical_summary.CLINICAL_CIRCLE_SOURCES`.

That was arm-blind in a way the data is not. ``compute_keep``'s germline branch removes
``CancerVar`` and every ``CIViC_*`` column; the somatic ``KEEP`` list has never held
``InterVar`` or ``RENOVO_Class``. So a germline clinician read two tool names that could
say nothing about their variants and a somatic one read two others, in a key explaining a
code neither could fully use.

Why the column narrowed rather than the key
-------------------------------------------
A key cannot drop a name without breaking the positional count it exists to serve: the
circle string was always six characters, so "the third position is CancerVar" had to stay
true even where the third position was blank. Annotating the dead positions instead —
greying them, or ``(not in germline)`` — keeps that count but prints a claim about the
user's file derived from a list of names, which is the shape issue #90's line was drawn
against. So the **column** stopped drawing them, and the key names exactly what is there.

That dissolves the ticket's second question rather than answering it. ``⬜`` used to mean
*this source had nothing to say about this variant* **and** *this arm never consults this
source* — the second being a whole column of ``⬜`` a reader could only identify by
scrolling. A position now exists only where its column does, so ``⬜`` has one job.

The intersection, and why both halves are needed
------------------------------------------------
A source is drawn when the arm's ``compute_keep`` emits its column **and** the file
carries it. Measured across 167 distinct annotated MAFs on this machine, each half rules
out a different false claim:

* **arm**: 1 of 64 germline MAFs carries ``CancerVar``. Reading the file alone would draw
  a circle for a column ``compute_keep`` removes, so it would appear in neither the table
  nor the export — a summary of something the reader cannot open and check. No fixture has
  that shape, which is why :func:`test_circle_sources_is_an_intersection` exists: mutating
  the arm half away is caught by nothing else in this module.
* **file**: ``compute_keep`` names ``InterVar``, ``RENOVO_Class``, ``ESCAT`` and
  ``ClinVar_VCF_CLNSIG`` unconditionally for their arm, and 2 of 64 germline MAFs do not
  carry ``InterVar``. Reading the arm alone would name a source over a column of ``⬜``.

CIViC is the case that *looks* like the second and is not: ``compute_keep`` guards those
columns itself, so they never reach the presence test. It still shows the intersection
working — ``somatic_reference.maf`` carries no CIViC and ``somatic_civic.maf`` does, and
the key differs between them on the same arm.

Counting the strip
------------------
The positional claims here are counts of the strip, and a strip is counted in **glyphs**
(:func:`_glyphs`) rather than in code points. ``⚠️`` and ``🛡️`` are a symbol plus U+FE0F, so
the two counts differ the moment either is drawn, and ``len`` was measuring the narrower
one — passing only because no committed fixture happened to draw those two glyphs. The
fixtures now plant a row that does (:data:`TWO_CODE_POINT_CLINVAR`), and each width
assertion says so, so the shape cannot go missing without a test going red.

Claims are read off the **rendered** page, not off the module — issue #79's rule, after
issue #71's ``not hasattr`` guard passed over a version written out by hand two lines
away. The one exception is :func:`test_circle_sources_is_an_intersection`, which tests the
derivation itself rather than a claim made to a user.
"""

from __future__ import annotations

import sys
import unicodedata
from pathlib import Path

import pytest

STREAMLIT_APP = Path(__file__).resolve().parent.parent
if str(STREAMLIT_APP) not in sys.path:
    sys.path.insert(0, str(STREAMLIT_APP))

FIXTURES = STREAMLIT_APP / "tests" / "fixtures" / "parity"

GERMLINE = FIXTURES / "germline_reference.maf"
SOMATIC = FIXTURES / "somatic_reference.maf"
SOMATIC_CIVIC = FIXTURES / "somatic_civic.maf"

#: Every source column the germline arm can draw. Dropping all of them leaves an arm with
#: nothing to summarise, which is the empty case: no column, no key.
GERMLINE_SOURCE_COLUMNS = ["ClinVar_VCF_CLNSIG", "InterVar", "RENOVO_Class", "ESCAT"]

#: A ClinVar call whose circle is **two code points** — ``⚠️``, the symbol plus U+FE0F.
#:
#: ``CLINICAL_VALUE_MAPPING`` reads ``Established_risk_allele`` as ``Disease_Risk``, and
#: ``_CLASS_GLYPHS`` draws that class ``⚠️``. It is also in the pipeline's ``CLINVAR_PATHO``
#: list and matches no keep term, so the row it is written onto passes by pathogenic rescue
#: whatever else it carries — which is what lets any row of any fixture here be turned into
#: one, and why the constructed fixture set carries rows of this shape.
#:
#: Every glyph in :data:`components.clinical_summary._CLASS_GLYPHS` is one code point except
#: ``⚠️`` and ``🛡️``, so a strip drawn from these fixtures without this had one character per
#: source by accident, and the width assertions below could not tell the two counts apart.
TWO_CODE_POINT_CLINVAR = "Established_risk_allele"

_SCRIPT = """
import os
import sys
sys.path.insert(0, {app!r})

import streamlit as st

from config.pipeline_params import pipeline_params
from page_modules import data_loading
from utils import read_maf

maf = read_maf({maf!r})
drop = {drop!r}
if drop:
    for column in drop:
        assert column in maf.columns, column
    maf = maf.drop(columns=list(drop))

two_code_point = {two_code_point!r}
if two_code_point is not None:
    assert "ClinVar_VCF_CLNSIG" in maf.columns, maf.columns.tolist()
    maf.loc[maf.index[0], "ClinVar_VCF_CLNSIG"] = two_code_point

st.session_state.filter_params = pipeline_params({arm!r})
st.session_state.maf_source_name = os.path.basename({maf!r})
st.session_state.maf_data = maf
data_loading.apply_filters_to_data(show_messages=False)

data_loading.show_results_section()
"""


def _render(maf, arm, drop=None, two_code_point=None):
    """The results section over ``maf``, filtered on ``arm``, rendered once.

    ``filter_params`` is seeded rather than left to the app: the real app opens on
    whichever arm this machine's ``~/.mafigate`` cache last held, so a harness that did not
    say would be testing that cache.

    ``two_code_point`` writes that ClinVar call onto the first row before filtering, so the
    rendered strip contains a glyph wider than one code point. Pass
    :data:`TWO_CODE_POINT_CLINVAR`; the shape is described there.
    """
    pytest.importorskip("streamlit", reason="streamlit not installed")
    from streamlit.testing.v1 import AppTest

    app = AppTest.from_string(
        _SCRIPT.format(
            app=str(STREAMLIT_APP),
            maf=str(maf),
            arm=arm,
            drop=drop,
            two_code_point=two_code_point,
        )
    )
    app.run(timeout=180)
    assert not app.exception, [str(e.value) for e in app.exception]

    # The page turns its own render failures into `st.error`, so an exception-free run is
    # not by itself a rendered results section.
    assert not app.error, [element.value for element in app.error]
    return app


@pytest.fixture(scope="module")
def germline():
    """A germline MAF carrying all four of its arm's sources, one row drawing ``⚠️``."""
    return _render(GERMLINE, "germline", two_code_point=TWO_CODE_POINT_CLINVAR)


@pytest.fixture(scope="module")
def somatic():
    """A somatic MAF with no CIViC annotation — the common somatic shape.

    Carries the same two-code-point row as :func:`germline`, because the width claim is
    made on both arms and each arm draws a different number of positions.
    """
    return _render(SOMATIC, "somatic", two_code_point=TWO_CODE_POINT_CLINVAR)


@pytest.fixture(scope="module")
def somatic_with_civic():
    """A somatic MAF that does carry CIViC, so the arm's fourth source can speak."""
    return _render(SOMATIC_CIVIC, "somatic")


@pytest.fixture(scope="module")
def germline_clinvar_only():
    """A germline MAF stripped down to one source.

    The companion that keeps the empty case below from passing vacuously: it goes through
    the same stripping, and it still draws a key. ClinVar is the survivor because it is a
    guideline source, so variants still pass and the *passed* grid is still drawn.
    """
    return _render(
        GERMLINE,
        "germline",
        drop=("InterVar", "RENOVO_Class", "ESCAT"),
        two_code_point=TWO_CODE_POINT_CLINVAR,
    )


@pytest.fixture(scope="module")
def germline_without_sources():
    """A germline MAF stripped of every source the overview could summarise.

    Nothing passes such a file, and that is not an artefact of the fixture: the germline
    guideline mask is ``InterVar | ClinVar | RENOVO``, and pathogenic retention reads
    ``InterVar`` and ``ClinVar``, so a MAF with no annotation source has nothing to be kept
    by. The rejects are where this state is actually met, which is why the assertions below
    read ``failed_data``.
    """
    return _render(GERMLINE, "germline", drop=tuple(GERMLINE_SOURCE_COLUMNS))


#: Code points that extend the glyph before them instead of drawing one of their own.
#:
#: ``Mn``/``Me``/``Mc`` are the combining categories, and U+FE0F — the variation selector
#: that makes ``⚠️`` and ``🛡️`` two code points — is ``Mn``, so it is covered by the category
#: rather than named. The emoji skin-tone modifiers are ``Sk``, a category that also holds
#: standalone characters like ``^``, so they are listed by code point instead of by class.
_EMOJI_MODIFIERS = {chr(cp) for cp in range(0x1F3FB, 0x1F400)}

#: The halves of a flag. Two of these are one glyph, four are two, and an odd one left over
#: is a glyph on its own \u2014 the only shape here that cannot be decided by looking at the code
#: point alone, which is why :func:`_glyphs` carries a second piece of state for it.
_REGIONAL_INDICATORS = {chr(cp) for cp in range(0x1F1E6, 0x1F200)}

#: Escaped rather than written, because the character itself is invisible in this file.
_ZWJ = "\u200d"


def _glyphs(strip):
    """``strip`` split into glyphs — grapheme clusters, not code points.

    ``len`` counts code points, and a glyph is not always one: this is the fault
    :func:`components.clinical_summary.pathogenicity_circle_glyphs` was written against, and
    measuring the strip with ``len`` reintroduces it on the other side. A strip of N sources
    is N glyphs however many code points each of them costs.

    Segmentation, not a lookup of the glyphs the app draws: a table of known glyphs would
    pass on a strip made of anything else, which is the opposite of what a width guard is
    for. Deliberately no dependency either — ``regex``'s ``\\X`` is the general answer and
    this app installs neither ``regex`` nor ``grapheme``.

    Not the whole of UAX #29 — Hangul jamo and the Indic conjunct rules are absent, and no
    glyph reachable from ``_CLASS_GLYPHS`` or writable into a MAF cell needs them. What is
    here is every shape emoji are built from: combining marks, U+FE0F, skin-tone modifiers,
    ZWJ sequences and flags. That list is what the tests below hold it to.
    """
    clusters = []
    after_zwj = False
    open_flag = False
    for char in strip:
        # A regional indicator joins the one before it only if that one is still unpaired,
        # so ``🇮🇹🇫🇷`` is two flags rather than one run of four halves.
        closes_flag = open_flag and char in _REGIONAL_INDICATORS
        extends = (
            after_zwj
            or char == _ZWJ
            or char in _EMOJI_MODIFIERS
            or closes_flag
            or unicodedata.category(char) in {"Mn", "Me", "Mc"}
        )
        if clusters and extends:
            clusters[-1] += char
        else:
            clusters.append(char)
        open_flag = char in _REGIONAL_INDICATORS and not closes_flag
        after_zwj = char == _ZWJ
    return clusters


def _glyph_widths(strips):
    """Every distinct glyph count down a ``Pathogenicity_Overview`` column.

    **Refuses a column no wider than its code points.** The planted row
    (:data:`TWO_CODE_POINT_CLINVAR`) is what separates counting glyphs from counting
    characters, and without it every width claim in this module passes for the reason it
    used to pass — one code point per source, by accident of which glyphs the fixtures
    happened to draw. So the shape is asserted here, once, rather than beside each claim.
    """
    assert any(len(strip) > len(_glyphs(strip)) for strip in strips), (
        "no strip costs more code points than it has glyphs, so this column cannot tell "
        "the two counts apart and the planted two-code-point row is gone"
    )
    return {len(_glyphs(strip)) for strip in strips}


def _keys(app):
    """Every Pathogenicity Overview key the section drew.

    There is one per grid, so a passed-and-failed report draws two; they are built from one
    list and must agree, which :func:`test_both_grids_draw_the_same_key` pins.
    """
    return [
        str(element.value)
        for element in app.markdown
        if "Pathogenicity Overview" in str(element.value)
    ]


def _order_line(key):
    """The ``Order:`` half of a key — the abbreviations and the names they expand to."""
    marker = "<b>Order:</b>"
    assert marker in key, key
    return key.split(marker, 1)[1]


def test_the_glyph_counter_is_not_just_len():
    """The instrument the width claims are made with, measured before they lean on it.

    ``_glyphs`` is the only reason those claims mean "one source per position" rather than
    "one code point per position", and a splitter that quietly degraded to ``len`` would
    make them pass again for the old wrong reason. So the cases are asserted here, together
    with the ``len`` each would have given.
    """
    assert _glyphs("⚠️") == ["⚠️"] and len("⚠️") == 2
    assert _glyphs("🛡️") == ["🛡️"] and len("🛡️") == 2

    strip = "⚠️🔵🟢⬜"
    assert _glyphs(strip) == ["⚠️", "🔵", "🟢", "⬜"]
    assert len(strip) == 5, "the fixture strip stopped being the case this test is about"

    # Single-code-point glyphs are left alone, so nothing merges that should not: the
    # failure a too-eager splitter would make is undercounting, not overcounting.
    assert _glyphs("🔴🟠🟡") == ["🔴", "🟠", "🟡"]
    assert _glyphs("") == []

    # The emoji shapes the app does not draw today, one per construction the docstring
    # claims: a skin-tone modifier, a ZWJ sequence, a keycap, and a flag. None is reachable
    # from `_CLASS_GLYPHS`, and each is what a glyph added later could be — a splitter that
    # broke one of them would fail a width claim for a counting quirk, which is the failure
    # this whole change exists to stop.
    assert _glyphs("👍🏽") == ["👍🏽"]
    assert _glyphs("🏳️‍🌈") == ["🏳️‍🌈"]
    assert _glyphs("1️⃣") == ["1️⃣"]
    assert _glyphs("🇮🇹") == ["🇮🇹"]

    # Flags pair off left to right rather than running together, and an odd half is its own
    # glyph. This is the one case a per-code-point rule cannot decide.
    assert _glyphs("🇮🇹🇫🇷") == ["🇮🇹", "🇫🇷"]
    assert _glyphs("🇮🇹🇫") == ["🇮🇹", "🇫"]
    assert _glyphs("🔵🇮🇹") == ["🔵", "🇮🇹"]


def test_the_grids_actually_rendered(germline):
    """The anchor under every absence assertion below, which is worthless without it.

    "``CancerVar`` is not in the key" passes just as well when there is no key — a section
    that returned early, a fixture that filtered to nothing, a harness that drew a blank
    page. So the key's own presence is asserted first, and separately.
    """
    assert _keys(germline), "the results section drew no Pathogenicity Overview key"


def test_the_germline_key_names_only_what_the_germline_arm_emits(germline):
    """The ticket's complaint, in the state it complained about.

    ``CancerVar`` and ``CIViC`` are removed from the germline output by ``compute_keep``,
    so before this change a germline reader was handed both names against positions that
    could only ever be ``⬜``.
    """
    order = _order_line(_keys(germline)[0])

    assert "CV IV RN ES" in order, order
    assert "ClinVar, InterVar, RENOVO, ESCAT" in order, order

    for absent in ("CancerVar", "CIViC", " CA ", " CI "):
        assert absent not in order, (
            f"{absent!r} is named in the germline key, and the germline pipeline does not "
            f"emit it: {order}"
        )


def test_the_somatic_key_names_only_what_the_somatic_arm_emits(somatic):
    """The mirror case, which the ticket noted and no surface had ever acted on."""
    order = _order_line(_keys(somatic)[0])

    assert "CV CA ES" in order, order
    assert "ClinVar, CancerVar, ESCAT" in order, order

    for absent in ("InterVar", "RENOVO", " IV ", " RN "):
        assert absent not in order, (
            f"{absent!r} is named in the somatic key, and the somatic pipeline does not "
            f"emit it: {order}"
        )


def test_the_file_decides_civic_on_the_somatic_arm(somatic, somatic_with_civic):
    """The half of the intersection the arm alone cannot supply.

    Both files are filtered on the same arm, and ``compute_keep`` emits the CIViC columns
    for it. Only one of them carries the annotation, and only that one gets the position —
    which is what stops the key naming ``CIViC`` over a blank column in the two thirds of
    somatic MAFs measured that do not carry it.
    """
    without = _order_line(_keys(somatic)[0])
    with_civic = _order_line(_keys(somatic_with_civic)[0])

    assert "CIViC" not in without, without
    assert "CV CA CI ES" in with_civic, with_civic
    assert "ClinVar, CancerVar, CIViC, ESCAT" in with_civic, with_civic


def test_the_key_names_one_source_per_drawn_glyph_not_per_code_point(germline, somatic):
    """The key and the cells agree on how many positions there are.

    This is the claim the old key could not make: six names over six glyphs, two or three
    of which were structurally blank. Read off the frame the app stored rather than the
    grid, because AgGrid does not render under ``AppTest`` — what is checked is the value
    the export would carry too.

    Counted with :func:`_glyphs`, and that is the whole of the test's claim rather than a
    detail of how it reads the column. This used to count code points, which is the same
    number only while every glyph is one code point; ``⚠️`` and ``🛡️`` are two, as
    :func:`components.clinical_summary.pathogenicity_circle_glyphs` says in the module the
    strip comes from. So the fixtures now carry a row drawing ``⚠️``
    (:data:`TWO_CODE_POINT_CLINVAR`) and the count is of glyphs, which is what "one source
    per position" ever meant.
    """
    for app in (germline, somatic):
        named = _order_line(_keys(app)[0]).split("(", 1)[0].split()
        frame = app.session_state["filtered_data"]
        assert "Pathogenicity_Overview" in frame.columns, frame.columns.tolist()

        widths = _glyph_widths(frame["Pathogenicity_Overview"])
        assert widths == {len(named)}, (
            f"the key names {len(named)} positions ({named}) and the column draws "
            f"{sorted(widths)} glyphs"
        )


def test_every_drawn_position_has_a_column_in_the_frame(germline, somatic_with_civic):
    """``⬜`` has one job now, and this is what makes that true.

    A position is drawn only where its column is present, so a blank can no longer mean
    *this arm never consults this source*. If a source were ever drawn without its column,
    the glyph would be structurally blank down every row and the key would be back to
    explaining a code the reader cannot use.
    """
    for app in (germline, somatic_with_civic):
        sources = app.session_state["overview_sources"]
        assert sources, "no sources were recorded for a MAF that has them"

        present = set(app.session_state["maf_data"].columns)
        for _abbreviation, name, column in sources:
            assert column in present, (
                f"{name} is drawn a circle from {column}, which this MAF does not carry"
            )


def test_one_surviving_source_still_gets_a_key(germline_clinvar_only):
    """A stripped MAF still draws a key, down to a single glyph.

    This is what makes the empty case below a real assertion rather than an observation
    that a degraded fixture drew nothing at all.
    """
    order = _order_line(_keys(germline_clinvar_only)[0])
    assert "CV (ClinVar)" in order, order

    frame = germline_clinvar_only.session_state["filtered_data"]

    # One source, so one glyph — including on the planted ``⚠️`` row, whose single glyph is
    # two code points. The starkest form of the count this module makes: a strip of `len` 2
    # that is still one position.
    assert _glyph_widths(frame["Pathogenicity_Overview"]) == {1}


def test_a_maf_with_none_of_its_arms_sources_gets_no_overview(germline_without_sources):
    """The floor: an overview of nothing is a pinned blank column, so there is none.

    ``Clinical_Summary`` still reports — it says "🔍 No Clinical Data", which is the honest
    sentence for this file and needs no empty column beside it. The resolver drops absent
    app extras silently, so nothing downstream has to be told.

    Read off ``failed_data``: this MAF passes no variants, for the reason the fixture's
    docstring gives, so the rejects are the only frame with rows in them.
    """
    app = germline_without_sources

    assert app.session_state["overview_sources"] == []

    frame = app.session_state["failed_data"]
    assert len(frame) > 0, "the rejects are empty, so this asserts nothing"
    assert "Pathogenicity_Overview" not in frame.columns, (
        "an overview column was built for a MAF carrying none of its arm's sources"
    )
    assert "Clinical_Summary" in frame.columns, frame.columns.tolist()

    assert not _keys(app), (
        f"a key was drawn with no sources to explain: {_keys(app)}"
    )


def test_both_grids_draw_the_same_key(germline):
    """Passed and failed are two grids reading one list, so they cannot disagree.

    The two verbatim copies issue #62 found in the Table Interaction Guide had already
    drifted apart; this key is drawn from ``overview_sources`` in one function, and this
    holds it that way.
    """
    keys = _keys(germline)
    assert len(keys) >= 2, f"expected a key above each grid, found {len(keys)}"
    assert len(set(keys)) == 1, keys


def test_circle_sources_is_an_intersection():
    """The derivation itself, on the arm and file combinations the fixtures cannot reach.

    A germline MAF *carrying* ``CancerVar`` is the case that separates "follows the file"
    from "follows the arm", and 1 of the 64 germline MAFs measured for issue #95 was one.
    ``compute_keep`` drops the column from the germline output, so drawing its circle would
    summarise a column absent from both the table and the export.
    """
    from components.clinical_summary import circle_sources

    def names(arm, columns):
        return [name for _abbreviation, name, _column in circle_sources(arm, columns)]

    germline_columns = ["ClinVar_VCF_CLNSIG", "InterVar", "RENOVO_Class", "ESCAT"]
    somatic_columns = ["ClinVar_VCF_CLNSIG", "CancerVar", "ESCAT"]

    drawn = names(
        "germline", germline_columns + ["CancerVar", "CIViC_Evidence_Level"]
    )
    assert drawn == ["ClinVar", "InterVar", "RENOVO", "ESCAT"], drawn

    # And the mirror: a somatic MAF carrying the germline pair.
    drawn = names("somatic", somatic_columns + ["InterVar", "RENOVO_Class"])
    assert drawn == ["ClinVar", "CancerVar", "ESCAT"], drawn

    assert circle_sources("germline", ["Hugo_Symbol"]) == []


def test_circle_sources_answers_what_the_grid_would():
    """The third input, which is easy to drop and makes the key disagree with the table.

    ``resolve_visible_columns`` is a pure function of the arm, ``skip_civic`` and the
    columns present, and the grid passes all three (``variant_table._visible_columns``).
    A derivation here that took only two would draw a CIViC circle beside a table with no
    CIViC column in it — the one outcome this whole mechanism exists to prevent.

    ``skip_civic`` is ``False`` on every live path today, so this is the contract rather
    than a reachable bug; it is asserted because the flag being dormant is exactly what
    would let the two drift apart unnoticed.
    """
    from components.clinical_summary import circle_sources
    from config.columns import MissingColumnsWarning, resolve_visible_columns

    columns = ["ClinVar_VCF_CLNSIG", "CancerVar", "ESCAT", "CIViC_Evidence_Level"]

    for skip in (False, True):
        drawn = {column for _abbreviation, _name, column in circle_sources(
            "somatic", columns, skip_civic=skip
        )}
        # The stand-in above is four columns, so the resolver rightly warns about the 41
        # the somatic pipeline emits and it does not carry — and `pytest.ini` promotes
        # `MissingColumnsWarning` to an error on purpose, which is what made this test fail
        # from the day it merged. Expected here rather than silenced: `pytest.warns` also
        # asserts the warning still fires, so the resolver going quiet fails this too. The
        # alternative, listing all 45 columns, would test the resolver's happy path and stop
        # exercising the absence this test is about.
        with pytest.warns(MissingColumnsWarning):
            shown = set(
                resolve_visible_columns(
                    "somatic", skip_civic=skip, available_columns=columns
                )
            )
        assert drawn <= shown, (
            f"with skip_civic={skip} the overview draws {sorted(drawn - shown)}, which "
            "the table does not show"
        )

    assert "CIViC_Evidence_Level" not in {
        column for _abbreviation, _name, column in circle_sources(
            "somatic", columns, skip_civic=True
        )
    }
