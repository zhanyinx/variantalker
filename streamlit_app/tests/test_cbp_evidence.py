"""CancerVar's AMP/ASCO/CAP evidence: the parser, the section, and the arm gate (issue #189).

Three groups of guards, and they fail for three different reasons:

* **The parser** — pure, so these are plain value assertions. The interesting ones are the
  *refusals*: every score feeds a sum and the sum is the tier, so a parser that patched over an
  unreadable entry would badge a tier CancerVar never assigned.
* **The section** — what a clinician sees. Driven through a fake ``st`` and asserted against the
  HTML, because the table is one ``st.markdown`` and element counts cannot see inside it (#188
  measured that: a 1-row and a 12-row table are one element each).
* **The arm gate and the column name** — driven through the *real fixture MAFs*, so the padded
  and dot-mangled spellings of the evidence column are exercised as bytes on disk rather than as
  strings retyped in a test. See ``tests/fixtures/cancervar/README.md``.

Each guard here was made to fail before being trusted, per this repo's standing rule.
"""

import os
import sys

import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from components.cbp_evidence import (  # noqa: E402
    CBP_COUNT,
    _CBP_CRITERIA,
    _TIER_COLORS,
    parse_cancervar,
    tier_color,
    tier_label,
)

FIXTURES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fixtures", "cancervar")

# Byte-for-byte from a real somatic MAF, padding included — see the fixture README.
REAL = " CancerVar: 8#Tier_II_potential EVS=[1, 0, 2, 1, 0, 0, 0, 1, 2, -1, 1, 1] "


# ---------------------------------------------------------------------------
# The twelve criteria, as a table
# ---------------------------------------------------------------------------


def test_there_are_exactly_twelve_criteria_in_cancervars_own_order():
    """The vector is positional, so the table's order *is* the mapping.

    A criterion inserted or reordered here would silently relabel every score in every file:
    ``EVS=[...]`` carries no codes, only twelve numbers in CancerVar's order.
    """
    codes = [c["code"] for c in _CBP_CRITERIA]

    assert len(_CBP_CRITERIA) == CBP_COUNT
    assert codes == [f"CBP{n}" for n in range(1, 13)]


def test_only_cbp5_and_cbp6_are_marked_never_evaluated():
    """#185 measured them 0 on all 109,416 rows: they ``return 0`` unconditionally.

    They are marked rather than dropped because a clinician who reads ``CBP5: 0`` as "no mosaic
    evidence in this variant" has been misled — the zero is CancerVar declining to evaluate.
    """
    never = [c["code"] for c in _CBP_CRITERIA if c.get("never_evaluated")]

    assert never == ["CBP5", "CBP6"]


def test_every_value_a_criterion_can_take_has_a_sentence_of_its_own():
    """The substantive departure from ``_ACMG_CRITERIA``: the sentence knows the score.

    ``CBP8 = 1`` and ``CBP8 = -1`` are opposite claims about the same evidence source, so a
    criterion whose ``says`` map is missing one of its own values would fall back to the
    value-independent phrasing and describe neither.
    """
    for criterion in _CBP_CRITERIA:
        for value in criterion["values"]:
            assert value in criterion["says"], (
                f"{criterion['code']} can be {value} and has no sentence for it"
            )


def test_the_criteria_that_can_go_negative_are_exactly_cbp8_and_cbp10():
    """#188 rejected the germline two-column mirror on this fact.

    A right-hand "evidence against" column whose complete possible content is "ClinVar
    disagrees" and "the predictors disagree" is not a peer of a column that can hold ten.
    """
    negative = [c["code"] for c in _CBP_CRITERIA if min(c["values"]) < 0]

    assert negative == ["CBP8", "CBP10"]


# ---------------------------------------------------------------------------
# parse_cancervar — what it reads
# ---------------------------------------------------------------------------


def test_a_real_evidence_string_parses_to_its_twelve_scores_its_sum_and_its_tier():
    parsed = parse_cancervar(REAL)

    assert parsed == {
        "scores": [1, 0, 2, 1, 0, 0, 0, 1, 2, -1, 1, 1],
        "sum": 8,
        "printed_sum": 8,
        "tier": "Tier_II_potential",
    }


def test_the_leading_space_and_the_label_are_both_tolerated():
    """The real cell carries both: a leading space *and* a ``CancerVar:`` label."""
    bare = "8#Tier_II_potential EVS=[1, 0, 2, 1, 0, 0, 0, 1, 2, -1, 1, 1]"

    assert parse_cancervar(REAL) == parse_cancervar(bare)


def test_a_negative_printed_sum_keeps_its_sign():
    """28 real rows have both negative criteria negative, and print ``-2#``.

    A sum regex of ``(\\d+)`` would read that as ``2`` — the same magnitude, the wrong side of
    zero — and the disagreement note in the caption would then fire on every such row.
    """
    parsed = parse_cancervar(
        " CancerVar: -2#Tier_IV_benign EVS=[0, 0, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0] "
    )

    assert parsed["printed_sum"] == -2
    assert parsed["sum"] == -2


def test_an_already_parsed_dict_is_returned_unchanged():
    """Mirrors ``parse_intervar``: the renderer may be handed either."""
    parsed = parse_cancervar(REAL)

    assert parse_cancervar(parsed) is parsed


# ---------------------------------------------------------------------------
# parse_cancervar — what it refuses, and why None beats a half-answer
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "value,why",
    [
        (" CancerVar: 1#Tier_IV_benign ", "a tier and a sum with no vector — 73 real cells"),
        ("EVS=[1, 0, 2, 1, 0, 0, 0, 1, 2, -1, 1]", "eleven scores: a truncated line"),
        ("EVS=[1, 0, 2, 1, 0, 0, 0, 1, 2, -1, 1, 1, 1]", "thirteen scores"),
        ("EVS=[1, 0, 2, 1, 0, 0, 0, 1, 2, x, 1, 1]", "a non-integer entry"),
        ("EVS=[]", "an empty vector"),
        (" CancerVar: 8#Tier_II_potential ", "the tier alone"),
        ("no evidence here at all", "nothing resembling a vector"),
    ],
)
def test_an_unreadable_vector_is_refused_rather_than_patched(value, why):
    """``None`` rather than a half-answer, and stricter than ``parse_intervar`` on purpose.

    ``parse_intervar`` zeroes a list entry it cannot read, because an ACMG criterion it cannot
    read is one it cannot claim was met — the classification comes from elsewhere in the string.
    Here every score feeds the **sum**, and the sum *is* the tier, so zeroing one unreadable
    entry would move the variant down CancerVar's own scale and put a tier on screen that the
    tool never assigned.
    """
    assert parse_cancervar(value) is None, why


@pytest.mark.parametrize(
    "value", [None, float("nan"), "", "   ", ".", "__UNKNOWN__", "UNKNOWN"]
)
def test_a_cell_that_says_nothing_parses_to_nothing(value):
    """The annotators' sentinels, from ``config/missing_values``.

    ``.`` alone occurs 39,423,699 times across the corpus; a parser that treated it as text
    would reach the ``EVS`` search and return ``None`` anyway, so this guard is about the
    predicate staying the shared one rather than a set retyped here.
    """
    assert parse_cancervar(value) is None


def test_a_vector_with_no_tier_prefix_still_parses():
    """The vector is the evidence; the tier is a value printed beside it.

    Never seen in 119,194 real vectors, so this is a not-fabricating guard: the scores are
    readable and must be shown, and the *renderer* is what declines to invent a tier.
    """
    parsed = parse_cancervar("EVS=[1, 0, 2, 1, 0, 0, 0, 1, 2, -1, 1, 1]")

    assert parsed["scores"] == [1, 0, 2, 1, 0, 0, 0, 1, 2, -1, 1, 1]
    assert parsed["sum"] == 8
    assert parsed["printed_sum"] is None
    assert parsed["tier"] == ""


# ---------------------------------------------------------------------------
# The tier: names, colours, and the one definition
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "tier,color,label",
    [
        ("Tier_I_strong", "#d62728", "Tier I — strong"),
        ("Tier_II_potential", "#ff7f0e", "Tier II — potential"),
        ("Tier_III_Uncertain", "#f0c420", "Tier III — uncertain"),
        ("Tier_IV_benign", "#2ca02c", "Tier IV — benign"),
    ],
)
def test_each_tier_has_its_own_colour_and_its_own_words(tier, color, label):
    """``"Tier_I" in val`` matched all four names — issue #187's defect, kept fixed here.

    The table moved into this module in #189 so that the guideline-row badge and this section's
    badge read one definition; the shadowing rule it has to satisfy is guarded in
    ``tests/test_components.py``.
    """
    assert tier_color(tier) == color
    assert tier_label(tier) == label


def test_a_tier_outside_cancervars_vocabulary_is_unknown_not_benign():
    """Not green: green is the benign colour, and an unknown value is not reassuring.

    And the name is echoed rather than smoothed into one of the four — the file said something
    this module does not know, and inventing a tier for it is the failure being guarded.
    """
    assert tier_color("Tier_V_invented") == "#7f7f7f"
    assert tier_label("Tier_V_invented") == "Tier_V_invented"


def test_the_tier_table_is_the_only_definition_of_these_colours():
    """One definition, so the two badges on one panel cannot disagree.

    ``components/variant_detail`` had its own copy until #189. If a second literal table of
    CancerVar tier colours appears anywhere in the app, the same tier can be painted two colours
    on one screen — which is the class of defect #187 was filed for.
    """
    import pathlib

    app_root = pathlib.Path(__file__).resolve().parent.parent
    hits = []
    for path in sorted(app_root.rglob("*.py")):
        if "tests" in path.parts or "vendor" in path.parts:
            continue
        source = path.read_text(encoding="utf-8")
        for tier, color, _label in _TIER_COLORS:
            # A file that names *both* a tier and its hex colour is defining the pairing.
            if f'"{tier}' in source and color in source:
                hits.append((path.name, tier))

    assert {name for name, _ in hits} <= {"cbp_evidence.py"}, (
        f"a second definition of CancerVar's tier colours: {hits}"
    )


# ---------------------------------------------------------------------------
# The section, as a clinician sees it
# ---------------------------------------------------------------------------


def _render(value):
    """Everything ``render_cbp_evidence`` draws, as (kind, text) pairs.

    Expanders are entered rather than skipped, and labelled, so a guard can tell "behind the
    *scored zero* expander" from "in the table" — the whole point of #188's split remainder.
    """
    from components import cbp_evidence

    drawn = []

    class _Expander:
        def __init__(self, label):
            self.label = label

        def __enter__(self):
            drawn.append(("expander", self.label))
            return self

        def __exit__(self, *exc):
            return False

    class _Fake:
        def markdown(self, text, *a, **k):
            drawn.append(("markdown", str(text)))

        def caption(self, text, *a, **k):
            drawn.append(("caption", str(text)))

        def info(self, text, *a, **k):
            drawn.append(("info", str(text)))

        def expander(self, label, *a, **k):
            return _Expander(label)

    original = cbp_evidence.st
    cbp_evidence.st = _Fake()
    try:
        cbp_evidence.render_cbp_evidence(value)
    finally:
        cbp_evidence.st = original
    return drawn


def _all_text(drawn):
    return "\n".join(text for _kind, text in drawn)


def _table_rows(drawn):
    """The ``<tr>`` count of the *first* table drawn — the one not behind an expander.

    #188 measured that element counts cannot answer "is this too much panel": every candidate
    shape draws its table as a single ``st.markdown``, so a 1-row and a 12-row table count
    identically. What a clinician scrolls is rows.
    """
    for kind, text in drawn:
        if kind == "expander":
            break
        if kind == "markdown" and "<table" in text:
            return text.count("<tr") - 1  # minus the header row
    return 0


def test_the_table_holds_only_the_criteria_that_fired():
    """"Fired" means non-zero. The commonest real vector is the reason.

    ``[0,0,0,0,0,0,0,-1,2,0,1,0]`` — 15,313 of 114,053 rows — is in COSMIC and ICGC, in a cancer
    gene, and ClinVar calls it **benign**. A positive-only table renders that as one row and
    silently discards the disagreement, which is the most clinically interesting thing in it.
    """
    drawn = _render(" CancerVar: 2#Tier_IV_benign EVS=[0, 0, 0, 0, 0, 0, 0, -1, 2, 0, 1, 0] ")
    table = next(text for kind, text in drawn if kind == "markdown" and "<table" in text)

    assert _table_rows(drawn) == 3
    for code in ("CBP8", "CBP9", "CBP11"):
        assert code in table
    assert "ClinVar calls this benign or likely benign" in table, (
        "the negative criterion was dropped, or described without knowing its value"
    )


def test_a_negative_criterion_reads_as_the_opposite_of_its_positive_self():
    """One sentence per code cannot serve both signs — #188's value-specific decision."""
    pathogenic = _all_text(
        _render("1#Tier_IV_benign EVS=[0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0]")
    )
    benign = _all_text(
        _render("-1#Tier_IV_benign EVS=[0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0]")
    )

    assert "ClinVar calls this pathogenic" in pathogenic
    assert "ClinVar calls this benign or likely benign" in benign
    assert "ClinVar calls this pathogenic" not in benign


def test_the_strongest_claim_comes_first_and_a_dissent_follows_what_it_qualifies():
    """Sorted by magnitude, negatives after positives of equal magnitude.

    So a dissent reads as a qualification of what precedes it rather than as the headline.
    """
    drawn = _render(" CancerVar: 8#Tier_II_potential EVS=[1, 0, 2, 1, 0, 0, 0, 1, 2, -1, 1, 1] ")
    table = next(text for kind, text in drawn if kind == "markdown" and "<table" in text)
    order = [
        code
        for code in ("CBP3", "CBP9", "CBP1", "CBP4", "CBP8", "CBP11", "CBP12", "CBP10")
        if code in table
    ]
    positions = sorted(order, key=table.index)

    # The two +2s, then the five +1s, then the one -1.
    assert positions[:2] == ["CBP3", "CBP9"]
    assert positions[-1] == "CBP10"


def test_all_twelve_at_zero_is_named_as_a_result_not_left_silent():
    """7.4% of real rows. Silence reads as a missing annotation — #187's exact defect."""
    drawn = _render(" CancerVar: 0#Tier_IV_benign EVS=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] ")
    kinds = {kind for kind, _ in drawn}
    text = _all_text(drawn)

    assert "info" in kinds, "the empty state drew no message at all"
    assert "not a missing annotation" in text
    assert _table_rows(drawn) == 0


def test_the_two_kinds_of_zero_are_not_one_list():
    """*Scored zero* and *never evaluated* are different facts about the file.

    CancerVar evaluated ten categories and found nothing either way; it does not evaluate CBP5
    or CBP6 at all. Collapsing them invites reading ``CBP5: 0`` as a finding about the variant.
    """
    drawn = _render(" CancerVar: 3#Tier_III_Uncertain EVS=[0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0] ")
    labels = [text for kind, text in drawn if kind == "expander"]

    assert len(labels) == 2
    assert "Scored zero" in labels[0] and "(7)" in labels[0]
    assert "Never evaluated" in labels[1] and "(2)" in labels[1]


def test_cbp5_and_cbp6_never_appear_in_the_table_of_what_fired():
    """They are 0 on every row of every file, so they can only ever be remainder."""
    drawn = _render(" CancerVar: 13#Tier_I_strong EVS=[2, 0, 2, 1, 0, 0, 1, 1, 2, 1, 1, 2] ")
    table = next(text for kind, text in drawn if kind == "markdown" and "<table" in text)

    assert "CBP5" not in table
    assert "CBP6" not in table


def test_the_tier_badge_carries_the_tier_from_the_string_in_the_panels_own_colour():
    drawn = _render(" CancerVar: 8#Tier_II_potential EVS=[1, 0, 2, 1, 0, 0, 0, 1, 2, -1, 1, 1] ")
    text = _all_text(drawn)

    assert "CancerVar AMP/ASCO/CAP tier" in text
    assert "Tier II — potential" in text
    assert tier_color("Tier_II_potential") in text


def test_a_vector_with_no_tier_is_said_to_have_none_rather_than_banded_into_one():
    """MAFigate does not add the twelve scores up — map #184's standing decision.

    The sum is right there and banding it would be easy, which is exactly why this is guarded:
    a tier on screen has to be one CancerVar wrote.
    """
    drawn = _render("EVS=[1, 0, 2, 1, 0, 0, 0, 1, 2, -1, 1, 1]")
    text = _all_text(drawn)

    assert "recorded no tier" in text
    assert "Tier II" not in text, "a tier was invented from the sum"


def test_the_arithmetic_is_one_caption_and_names_whose_thresholds_they_are():
    """#185 §6: the categories and tier names are the guideline's, the numbers are CancerVar's.

    And **no distance-to-next-tier** (#188): "2 points from Tier II" implies the app could move
    the variant, and invites reading CancerVar's thresholds as AMP/ASCO/CAP's.
    """
    drawn = _render(" CancerVar: 8#Tier_II_potential EVS=[1, 0, 2, 1, 0, 0, 0, 1, 2, -1, 1, 1] ")
    captions = [text for kind, text in drawn if kind == "caption"]
    text = "\n".join(captions)

    assert "+8" in text
    assert "the scores and the thresholds are\nCancerVar's" in text.replace("  ", " ") or (
        "thresholds are" in text and "CancerVar's" in text
    )
    assert "AMP/ASCO/CAP score" not in _all_text(drawn), "a false attribution"
    for forbidden in ("points from", "away from", "short of"):
        assert forbidden not in text, f"a distance-to-next-tier crept in: {forbidden!r}"


def test_a_printed_sum_that_disagrees_with_the_vector_is_named_not_hidden():
    """Never seen in 119,194 real vectors, and structural — but not asserted as impossible.

    The badge shows what CancerVar printed and the caption shows what the vector adds to, so if
    the two ever part company the panel says so rather than showing both as if they agreed.
    """
    text = _all_text(
        _render(" CancerVar: 99#Tier_I_strong EVS=[1, 0, 2, 1, 0, 0, 0, 1, 2, -1, 1, 1] ")
    )

    assert "+99" in text and "+8" in text
    assert "printed a sum" in text


def test_nothing_parseable_says_so_rather_than_drawing_an_empty_section():
    drawn = _render(" CancerVar: 1#Tier_IV_benign ")

    assert drawn == [("caption", "No CancerVar AMP/ASCO/CAP evidence available for this variant.")]


# ---------------------------------------------------------------------------
# The arm gate and the column name, driven through the real fixtures
# ---------------------------------------------------------------------------


def _fixture(name):
    return pd.read_csv(os.path.join(FIXTURES, name), sep="\t", low_memory=False)


def _panel(monkeypatch, row):
    """Every string the whole detail panel draws for one variant.

    The real panel, not the section alone: the questions here are *whether* the section attaches
    and *what the captions above it now say*, both of which live in ``variant_detail``.

    **Both modules' ``st`` are replaced.** ``render_cbp_evidence`` draws through
    ``cbp_evidence.st``, so patching only ``variant_detail.st`` leaves the section writing to the
    real Streamlit and returns a transcript with the section silently missing from it — which is
    indistinguishable from the section not being drawn at all, and would have made every guard
    below pass for the wrong reason had they been written the other way round.
    """
    from unittest.mock import MagicMock

    from components import cbp_evidence, variant_detail

    drawn = []

    class _Ctx:
        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

    fake = MagicMock()
    fake.columns.side_effect = lambda n, *a, **k: [
        _Ctx() for _ in range(n if isinstance(n, int) else len(n))
    ]
    fake.expander.side_effect = lambda *a, **k: _Ctx()
    fake.markdown.side_effect = lambda text, *a, **k: drawn.append(str(text))
    fake.caption.side_effect = lambda text, *a, **k: drawn.append(str(text))
    fake.info.side_effect = lambda text, *a, **k: drawn.append(str(text))
    fake.metric.side_effect = lambda label, value, *a, **k: drawn.append(f"{label}: {value}")
    monkeypatch.setattr(variant_detail, "st", fake)
    monkeypatch.setattr(cbp_evidence, "st", fake)
    monkeypatch.setattr(variant_detail, "render_acmg_evidence", lambda *a, **k: None)

    variant_detail.render_variant_detail_panel(row)
    return "\n".join(drawn)


def test_the_padded_column_name_real_files_carry_reaches_the_section(monkeypatch):
    """``' CancerVar: CancerVar and Evidence '`` — the spelling on 121 of 124 real files.

    Read off a fixture on disk rather than retyped, so the padding is bytes and not a literal a
    later edit could tidy away.
    """
    frame = _fixture("somatic_cancervar_evidence.maf")
    assert " CancerVar: CancerVar and Evidence " in frame.columns, "the fixture lost its padding"

    text = _panel(monkeypatch, frame.iloc[5])  # PIK3CA p.E545G, a real Tier I

    assert "CancerVar AMP/ASCO/CAP tier" in text
    assert "Tier I — strong" in text


def test_the_dot_mangled_column_name_reaches_the_section_too(monkeypatch):
    """``'CancerVar..CancerVar.and.Evidence'`` — 3 real files, R-mangled punctuation.

    ``_evidence_column``'s substring match survives both spellings *by luck rather than by
    design* (#188), which is why it is pinned.
    """
    frame = _fixture("somatic_cancervar_dotted_no_column.maf")
    assert "CancerVar..CancerVar.and.Evidence" in frame.columns

    text = _panel(monkeypatch, frame.iloc[0])  # TP53 p.T304A

    assert "CancerVar AMP/ASCO/CAP tier" in text
    assert "Tier III — uncertain" in text


def test_a_file_with_no_cancervar_column_is_not_told_it_holds_no_tier(monkeypatch):
    """The caption issue #188 obliged this ticket to reword.

    PR #195 told the reader of exactly these files — 16 of 124 — that the file "holds no
    AMP/ASCO/CAP tier". The tier is printed *inside* the evidence string, so with the section
    below badging it that sentence was false on screen while the tier sat beneath it.
    """
    frame = _fixture("somatic_cancervar_dotted_no_column.maf")
    assert "CancerVar" not in frame.columns, "the fixture grew the column it exists to lack"

    text = _panel(monkeypatch, frame.iloc[1])  # PIK3CA p.E545G

    assert "holds no AMP/ASCO/CAP tier" not in text
    assert "Tier I — strong" in text, "the tier the caption points at was not drawn"
    assert "printed inside its" in text and "own evidence string" in text


def test_a_tier_with_no_vector_and_no_column_is_not_promised_a_tier_below(monkeypatch):
    """The 73-row state, and the caption issue #210 found false on every one of them.

    All 73 real cells holding a tier and a sum with no ``EVS=[...]`` vector are in **one file,
    and that file has no ``CancerVar`` column** — which issue #189 assumed it did when it decided
    the section should draw nothing here ("the guideline row still badges it from the ``CancerVar``
    column"). With no column and no parse, the tier CancerVar printed is on screen **nowhere**,
    while the caption above said it was below.

    The dev kept the section's refusal and corrected the caption, so this test pins both halves:
    no tier is drawn, and nothing claims one is.
    """
    frame = _fixture("somatic_cancervar_dotted_no_column.maf")
    assert "CancerVar" not in frame.columns, "the fixture grew the column it exists to lack"

    text = _panel(monkeypatch, frame.iloc[2])  # KRAS, ' CancerVar: 1#Tier_IV_benign '

    assert "so it shows none" in text, "the caption did not say the tier is unavailable"
    assert "tier below is the one" not in text, "the caption promised a tier that is not drawn"
    assert "Tier IV" not in text, "a tier was drawn from a cell the parser refused"
    assert "No CancerVar AMP/ASCO/CAP evidence available for this variant." in text


def test_an_empty_evidence_cell_with_no_column_says_no_tier_is_shown(monkeypatch):
    """No ``CancerVar`` column and an evidence cell that says nothing.

    **0 rows in the corpus** are in this state, and the branch is reachable all the same:
    ``_evidence_value`` returns ``None`` for a blank cell, so the section below is not drawn at
    all and a caption pointing at it would point at nothing. Pinned rather than left to chance,
    because "no real row does this" is what the caption above was relying on when it was false.
    """
    frame = _fixture("somatic_cancervar_dotted_no_column.maf")

    text = _panel(monkeypatch, frame.iloc[3])  # EGFR, evidence cell '.'

    assert "so it shows none" in text
    assert "tier below is the one" not in text
    assert "CancerVar AMP/ASCO/CAP tier" not in text, "a section was drawn for a blank cell"


def test_a_vector_summing_to_zero_still_has_a_tier_to_point_at(monkeypatch):
    """**2,115 real rows in 13 files**: no ``CancerVar`` column, and a vector summing to zero.

    The tier is drawn — ``Tier_IV_benign``, from the string — so the caption must point at it. This
    exists because a mutation survived: reading the verdict from the parsed **sum** rather than the
    **tier** passed every other guard, and would say *"it shows none"* here, because
    ``0 or ""`` is ``""``. A tier of Tier IV is not the absence of a tier, and 7.4% of real rows
    score zero on all twelve criteria (#185), so the two are not interchangeable at the margin —
    they are interchangeable on 2,115 rows.
    """
    frame = _fixture("somatic_cancervar_dotted_no_column.maf")

    text = _panel(monkeypatch, frame.iloc[4])  # BRAF, EVS all zero, '0#Tier_IV_benign'

    assert "Tier IV — benign" in text, "the tier the caption points at was not drawn"
    assert "tier below is the one" in text, "a drawn tier was described as absent"
    assert "shows none" not in text


def test_the_marker_sentence_does_not_claim_an_arm_the_tier_below_contradicts(monkeypatch):
    """6,524 real rows: a germline marker above, and CancerVar's somatic tier below it.

    The one file shape that makes the arm word wrong — ``RENOVO_Class`` present, no ``InterVar``,
    no ``CancerVar`` **column**, and CancerVar's evidence vector in the file all the same. The
    caption used to open *"This file was annotated for germline analysis"* directly above a badged
    AMP/ASCO/CAP tier.

    Driven through the whole panel with **both** modules' ``st`` replaced, because the
    contradiction is between a caption ``variant_detail`` draws and a badge ``cbp_evidence`` draws:
    patch one and the test can only see half of it.
    """
    row = pd.Series(
        {
            "Hugo_Symbol": "TP53",
            "Chromosome": "17",
            "RENOVO_Class": "Pathogenic",
            " CancerVar: CancerVar and Evidence ": (
                " CancerVar: 10#Tier_II_potential EVS=[2, 2, 0, 1, 0, 0, 0, 1, 1, 0, 1, 2] "
            ),
        }
    )

    text = _panel(monkeypatch, row)

    assert "Tier II — potential" in text, "the somatic tier this is about was not drawn"
    assert "annotated for" not in text, (
        "the caption claims an arm while the tier below contradicts it"
    )
    assert (
        "This file carries `RENOVO_Class` but no `InterVar` column, so MAFigate has no ACMG/AMP "
        "classification to show for this variant." in text
    )


def test_a_variant_whose_evidence_cell_says_nothing_draws_no_section(monkeypatch):
    """Both columns present, both cells ``.`` — a *variant*-level absence, not a file-level one.

    #187 counted 53 of 24,330 somatic rows in this state. The tool ran and reached no verdict
    here, which is a different sentence from "this file has no tier".
    """
    frame = _fixture("somatic_cancervar_evidence.maf")
    text = _panel(monkeypatch, frame.iloc[7])  # BRAF, both cells '.'

    assert "CancerVar AMP/ASCO/CAP tier" not in text
    assert "CancerVar recorded no tier for this variant." in text


def test_a_tier_with_no_vector_keeps_the_guideline_badge_and_draws_no_evidence(monkeypatch):
    """73 real cells hold ``'1#Tier_IV_benign'`` with no vector, all in one file.

    The tier is real and the guideline row still badges it from the ``CancerVar`` column; what
    is missing is the reasoning this section exists to show, so the section says so rather than
    badging a tier with nothing behind it.
    """
    frame = _fixture("somatic_cancervar_evidence.maf")
    text = _panel(monkeypatch, frame.iloc[6])  # KRAS p.G12D

    assert "AMP/ASCO/CAP (CancerVar)" in text, "the guideline-row badge went missing"
    assert "No CancerVar AMP/ASCO/CAP evidence available for this variant." in text
    assert "CancerVar AMP/ASCO/CAP tier" not in text


def test_a_germline_row_draws_no_cancervar_section(monkeypatch):
    """Strict arm-gating, and the gate is the row's own columns (#187).

    ``CancerVar`` appears on 0 of 57 real germline files, so column presence gates without
    reading ``filter_params`` — which would hide a present verdict from a user on the wrong arm,
    against issue #135's *warn, never override*.
    """
    germline = pd.Series(
        {
            "Hugo_Symbol": "BRCA1",
            "Chromosome": "17",
            "Start_Position": 43000001,
            "Reference_Allele": "C",
            "Tumor_Seq_Allele2": "T",
            "InterVar": "Pathogenic",
            "InterVar: InterVar and Evidence": "InterVar: Pathogenic PVS1=1 PS=[0,0,0,0,0]",
            "RENOVO_Class": "High",
        }
    )

    text = _panel(monkeypatch, germline)

    assert "CancerVar" not in text
