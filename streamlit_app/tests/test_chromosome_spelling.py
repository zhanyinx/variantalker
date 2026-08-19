"""Where the ``chr`` prefix comes from, and that no surface adds a second one (issue #149).

Before this, the app answered *does this chromosome need a ``chr``?* at four render sites
and got it right at one. The other three wrote the prefix into the string unconditionally,
so on **174 of the 187 real MAFs** measured for #149 — every file that has been through
the pipeline — they printed ``chrchr17``: on the detail panel's ``Position`` line, on the
button a user presses to open a variant, and in the fallback row selector.

**The whole suite was green through all of it.** Nothing asserted what any of those four
sites renders, which is why this file exists and why every guard here is mutation-checked
against the old behaviour: a test that passes on both spellings is not a guard.

The three kinds of claim here:

* the spelling is settled **once**, by ``config.chromosome_spelling``, on the pipeline's own
  rule (whole-column, all-or-nothing, only when nothing already says ``chr``);
* every surface that renders a chromosome renders the value **as it is**, adding nothing;
* the load path runs it before the frame becomes ``maf_data``, and says so when it fires,
  because the re-spelling reaches the user's download.
"""

from __future__ import annotations

import sys
from pathlib import Path
from unittest.mock import MagicMock

import numpy as np
import pandas as pd

STREAMLIT_APP = Path(__file__).resolve().parent.parent
if str(STREAMLIT_APP) not in sys.path:
    sys.path.insert(0, str(STREAMLIT_APP))

from config.chromosome_spelling import (  # noqa: E402
    CHROMOSOME,
    normalise_chromosome_spelling,
)
from config.missing_values import says_nothing  # noqa: E402


def _maf(chromosomes, **extra):
    """A frame with a ``Chromosome`` column and enough beside it to be a MAF."""
    frame = pd.DataFrame(
        {
            CHROMOSOME: chromosomes,
            "Hugo_Symbol": ["TP53"] * len(chromosomes),
            "Start_Position": list(range(7577120, 7577120 + len(chromosomes))),
        }
    )
    for name, values in extra.items():
        frame[name] = values
    return frame


# ---------------------------------------------------------------------------
# The rule itself
# ---------------------------------------------------------------------------


def test_a_bare_column_gets_the_prefix():
    """The two real IEO germline MAFs that arrive this way — bare `1`, `2`, `X`."""
    maf = _maf(["1", "2", "X"])

    assert normalise_chromosome_spelling(maf) is True
    assert list(maf[CHROMOSOME]) == ["chr1", "chr2", "chrX"]


def test_a_prefixed_column_is_left_exactly_alone():
    """174 of 187 real MAFs. Mutation: drop the `.any()` guard and this reads `chrchr17`."""
    maf = _maf(["chr1", "chr17", "chrX"])

    assert normalise_chromosome_spelling(maf) is False
    assert list(maf[CHROMOSOME]) == ["chr1", "chr17", "chrX"]


def test_the_rule_is_whole_column_and_not_per_row():
    """The pipeline's decision, copied: one value saying `chr` settles the whole column.

    No real MAF mixes the two spellings — 0 of 184 carrying the column — precisely because
    the rule upstream is whole-column. A per-row rule would be a *different* rule that
    happens to agree on every file anyone has, and would quietly disagree on the first one
    that mixes them. Mutation: make this per-row and the bare `2` below becomes `chr2`.
    """
    maf = _maf(["chr1", "2"])

    assert normalise_chromosome_spelling(maf) is False
    assert list(maf[CHROMOSOME]) == ["chr1", "2"]


def test_an_integer_column_survives_the_assignment():
    """A bare all-numeric column reads as `int64`, not object.

    Writing strings into a subset of an int64 column through `.loc` leaves pandas to upcast
    mid-assignment. Building the replacement avoids the question entirely — and this is the
    shape the two real bare files would have if neither carried an `X`.
    """
    maf = _maf([1, 2, 3])
    assert maf[CHROMOSOME].dtype == np.dtype("int64")

    assert normalise_chromosome_spelling(maf) is True
    assert list(maf[CHROMOSOME]) == ["chr1", "chr2", "chr3"]


def test_a_blank_cell_stays_blank_rather_than_becoming_chrnan():
    """The one deliberate deviation from the pipeline's line.

    `"chr" + column.astype(str)` renders a missing cell as the string `chrnan`, which
    `says_nothing` does not recognise — so the detail panel would print `chrnan` at a
    clinician where it prints an em dash. Mutation: assign `rendered` wholesale instead of
    through the mask, and this row becomes `chrnan`.
    """
    maf = _maf(["1", None, "3"])

    assert normalise_chromosome_spelling(maf) is True
    assert list(maf[CHROMOSOME])[0] == "chr1"
    assert list(maf[CHROMOSOME])[2] == "chr3"
    assert says_nothing(list(maf[CHROMOSOME])[1]), (
        "a blank chromosome was re-spelled into a value — the panel draws this, and "
        "`chrnan` is a string no file contains and no reader can interpret"
    )


def test_a_column_of_nothing_is_not_touched_and_is_not_announced():
    """Nothing to prefix, so nothing happened, so the user is told nothing.

    The values are asserted as well as the return value: *not touched* is the half of this
    that a return value cannot show, and it is the half that matters — the pipeline's own
    expression would fill every cell here with ``chrnan``. This is the one column shape on
    which the app and the pipeline deliberately disagree; see
    ``tests/test_chromosome_rule_contract.py``, which asserts the disagreement rather than
    leaving it described.
    """
    maf = _maf([None, float("nan"), ""])
    before = list(maf[CHROMOSOME])

    assert normalise_chromosome_spelling(maf) is False
    assert all(says_nothing(value) for value in maf[CHROMOSOME])
    assert str(list(maf[CHROMOSOME])) == str(before)


def test_an_absent_column_is_not_this_function_s_problem():
    """`validate_required_columns` refuses such a file a moment later; refusing is its job."""
    maf = pd.DataFrame({"Hugo_Symbol": ["TP53"], "Start_Position": [7577120]})

    assert normalise_chromosome_spelling(maf) is False
    assert CHROMOSOME not in maf.columns


# That the app's rule *agrees with the pipeline's* is not asserted here. It cannot be —
# it needs the pipeline checkout, and everything above has to keep holding in one without
# it, which is the whole point of the `bin/`-free net. It lives in
# `tests/test_chromosome_rule_contract.py`, which is filed as a fifth module allowed to
# need `bin/` and argues there for why.


# ---------------------------------------------------------------------------
# What the surfaces render
# ---------------------------------------------------------------------------


def _prefixed_variant(**overrides):
    """A variant spelled the way 174 of 187 real MAFs spell it."""
    row = {
        "Hugo_Symbol": "TP53",
        "Variant_Classification": "Missense_Mutation",
        "Variant_Type": "SNP",
        CHROMOSOME: "chr17",
        "Start_Position": 7577120,
        "End_Position": 7577120,
        "Reference_Allele": "C",
        "Tumor_Seq_Allele2": "T",
    }
    row.update(overrides)
    return pd.Series(row)


def _capture_detail_writes(monkeypatch):
    """Stub ``variant_detail.st`` and hand back the list it writes into.

    One scaffold for the three readers below rather than three copies of it: the panel and
    the two link checks all need the same stub, and a `markdown` side effect that only two
    of three copies carried is exactly how a guard stops seeing the field it is about.

    Reads all three writers, not just ``markdown``. A claim about what a clinician sees has
    to cover the lot — the fields these guards are about are drawn by different ones, and a
    guard watching only ``markdown`` would pass over a metric that had stopped being drawn.
    """
    from components import variant_detail

    drawn = []

    class _Column:
        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

    fake = MagicMock()
    fake.columns.side_effect = lambda n, *a, **k: [
        _Column() for _ in range(n if isinstance(n, int) else len(n))
    ]
    fake.markdown.side_effect = lambda text, *a, **k: drawn.append(str(text))
    fake.caption.side_effect = lambda text, *a, **k: drawn.append(str(text))
    fake.metric.side_effect = lambda label, value, *a, **k: drawn.append(
        f"{label}: {value}"
    )
    monkeypatch.setattr(variant_detail, "st", fake)
    return variant_detail, drawn


def _detail_panel_text(monkeypatch, row):
    """Every string the detail panel draws, across all three of the writers it uses."""
    variant_detail, drawn = _capture_detail_writes(monkeypatch)
    monkeypatch.setattr(variant_detail, "render_acmg_evidence", lambda *a, **k: None)

    variant_detail.render_variant_detail_panel(row)
    return drawn


def _external_links(monkeypatch, row):
    """Every link the panel offers out, as the markdown it writes them in."""
    variant_detail, drawn = _capture_detail_writes(monkeypatch)

    variant_detail._render_external_links(row)
    return drawn


def test_the_detail_panel_names_the_position_once(monkeypatch):
    """The line a clinician reads to identify the variant in front of them.

    Mutation: put the `chr` literal back in `variant_detail.py` and this reads
    `**Position:** chrchr17:7577120-7577120`.
    """
    drawn = _detail_panel_text(monkeypatch, _prefixed_variant())

    position = [text for text in drawn if text.startswith("**Position:**")]
    assert position == ["**Position:** chr17:7577120-7577120"], (
        f"the detail panel's position line reads {position!r}"
    )


def test_no_surface_of_the_detail_panel_doubles_the_prefix(monkeypatch):
    """Every string the panel draws, not just the one field a guard was written for."""
    drawn = _detail_panel_text(monkeypatch, _prefixed_variant())

    doubled = [text for text in drawn if "chrchr" in text]
    assert not doubled, f"the panel drew a doubled chromosome prefix: {doubled!r}"


def test_the_ucsc_link_asks_for_the_chromosome_the_row_carries(monkeypatch):
    """UCSC wants `chr17`, and the value already is one — so nothing is added to it.

    This site is where the app's only prefix guard used to live. Deleting it is only safe
    because the column arrives settled; this reads the rendered URL rather than trusting
    that. Mutation: prefix here as well and the link asks UCSC for `chrchr17`.
    """
    drawn = _external_links(monkeypatch, _prefixed_variant())

    ucsc = [text for text in drawn if "genome.ucsc.edu" in text]
    assert ucsc, "the UCSC link was not drawn at all"
    assert "position=chr17:7577120-7577120" in "\n".join(ucsc), (
        f"the UCSC link does not name chr17: {ucsc!r}"
    )
    assert "chrchr" not in "\n".join(ucsc)


def test_the_chromosome_axis_still_reads_bare_numbers(monkeypatch):
    """The other consumer that wants the number without its prefix.

    ``components/charts.py`` strips unconditionally, so it was right under either spelling
    and is untouched by this change — but it is a *sixth* site answering the same question
    and nothing asserted it, which is how the four that got it wrong went unnoticed. The
    axis is the visible consequence: a chart labelled ``chr1`` where every other one on the
    dashboard says ``1`` would be this change leaking somewhere it was not meant to.
    """
    from components import charts

    frame = pd.DataFrame(
        {CHROMOSOME: ["chr1", "chr1", "chr2", "chrX"], "Hugo_Symbol": ["TP53"] * 4}
    )
    counts = frame[CHROMOSOME].astype(str).str.replace("chr", "", regex=False)

    assert sorted(set(counts)) == ["1", "2", "X"], (
        f"the chromosome axis would be labelled {sorted(set(counts))}"
    )
    assert "chr" not in "".join(counts)
    assert 'str.replace("chr", "", regex=False)' in Path(
        charts.__file__
    ).read_text(), (
        "components/charts.py no longer strips the prefix the way this guard models — "
        "re-derive the assertion above against what it now does"
    )


def test_the_gnomad_coordinate_link_strips_the_prefix(monkeypatch):
    """gnomAD wants the number bare, and still gets it now the column is always prefixed."""
    drawn = _external_links(monkeypatch, _prefixed_variant())

    gnomad = [text for text in drawn if "gnomad.broadinstitute.org" in text]
    assert gnomad, "the gnomAD link was not drawn at all"
    assert "/variant/17-7577120-C-T" in "\n".join(gnomad), (
        f"the gnomAD link does not carry a bare chromosome: {gnomad!r}"
    )


def test_the_view_details_button_names_the_variant_once(monkeypatch):
    """The label on the button a user presses to open a variant.

    Mutation: restore the `chr` literal in `variant_table.py` and this reads
    `🔍 View details: **TP53** chrchr17:7577120`.
    """
    from components import variant_table

    fake_st = MagicMock()
    fake_st.button.return_value = False
    monkeypatch.setattr(variant_table, "st", fake_st)
    monkeypatch.setattr(variant_table, "GridOptionsBuilder", MagicMock(), raising=False)

    grid_response = MagicMock()
    grid_response.selected_rows = [
        {"Hugo_Symbol": "TP53", CHROMOSOME: "chr17", "Start_Position": 7577120}
    ]
    monkeypatch.setattr(
        variant_table, "AgGrid", lambda *a, **k: grid_response, raising=False
    )

    frame = _maf(["chr17"])
    variant_table._render_aggrid_with_detail(
        frame, frame, 500, "passed_variants", list(frame.columns)
    )

    labels = [call.args[0] for call in fake_st.button.call_args_list if call.args]
    assert any("**TP53** chr17:7577120" in label for label in labels), (
        f"no button named the variant correctly: {labels!r}"
    )
    assert not [label for label in labels if "chrchr" in label], (
        f"a button label doubled the prefix: {labels!r}"
    )


def test_the_fallback_row_selector_names_the_variant_once(monkeypatch):
    """The selector shown where streamlit-aggrid is not installed.

    Reached by turning `AGGRID_AVAILABLE` off, and the label is read by calling the
    `format_func` the selectbox was handed — the string never exists anywhere else.
    """
    from components import variant_table

    fake_st = MagicMock()
    fake_st.checkbox.return_value = True
    fake_st.selectbox.return_value = None
    monkeypatch.setattr(variant_table, "st", fake_st)
    monkeypatch.setattr(variant_table, "AGGRID_AVAILABLE", False)

    frame = _maf(["chr17"])
    variant_table.create_data_table(frame, key_suffix="passed_variants")

    assert fake_st.selectbox.called, "the fallback selector was not drawn"
    format_func = fake_st.selectbox.call_args.kwargs["format_func"]
    label = format_func(0)

    assert "chr17:7577120" in label, f"the selector label reads {label!r}"
    assert "chrchr" not in label, f"the selector label doubled the prefix: {label!r}"


# ---------------------------------------------------------------------------
# The load path, which is the only place any of this runs
# ---------------------------------------------------------------------------


def _open(monkeypatch, maf, refused=False, **seed):
    """Drive the real shared load tail, not a helper only this test calls.

    ``_open_the_file_just_read`` is the one door all three load paths funnel through — the
    sidebar's chooser and the OS "Open With" handler's two readers — so a claim made here
    is a claim about the app rather than about one route into it. The same reasoning, and
    the same seam, as ``tests/test_file_chooser.py``'s.
    """
    from tests.fakes import FakeSessionState
    from page_modules import data_loading

    state = FakeSessionState(**seed)
    fake_st = MagicMock()
    fake_st.session_state = state
    monkeypatch.setattr(data_loading, "st", fake_st)
    monkeypatch.setattr(
        data_loading, "validate_required_columns", lambda data: not refused
    )
    monkeypatch.setattr(
        data_loading, "apply_filters_to_data", lambda show_messages=True: True
    )

    data_loading._open_the_file_just_read(maf, "sample.maf")
    return state


def test_the_stored_frame_and_the_note_key_are_both_settled(monkeypatch):
    """What `maf_data` holds, and what the note key makes of it.

    Named for what it can detect. It was first called *…before it becomes maf_data*, which
    asserts an **ordering** — and no assertion here can see one, because the frame is
    mutated in place and is the same object either side of the assignment. The claim that
    is worth making, and that this does make, is about the values the rest of the app then
    reads: the stored column, and the identity `_variant_key` builds out of it.
    """
    from components.variant_table import _rendered_key_components

    state = _open(monkeypatch, _maf(["1", "2"], Reference_Allele=["C", "G"],
                                    Tumor_Seq_Allele2=["T", "A"]))

    stored = state["maf_data"]
    assert list(stored[CHROMOSOME]) == ["chr1", "chr2"]
    assert _rendered_key_components(stored.iloc[0])[1] == "chr1", (
        "the note key was built from an unsettled chromosome, so the same variant would "
        "be filed under two identities depending on how its file spelled it (issue #67)"
    )


def test_opening_a_bare_file_parks_the_notice(monkeypatch):
    """It is stashed, not drawn: the load filters silently and jumps straight to Results."""
    from page_modules.data_loading import _CHROMOSOMES_RESPELLED

    state = _open(monkeypatch, _maf(["1", "2"]))

    assert state.get(_CHROMOSOMES_RESPELLED) is True


def test_opening_an_already_prefixed_file_says_nothing(monkeypatch):
    """174 of 187 real MAFs. Nothing was changed, so there is nothing to report.

    Mutation: stash unconditionally rather than on the return value → every user of every
    pipeline MAF is told their chromosomes were re-spelled, which is false.
    """
    from page_modules.data_loading import _CHROMOSOMES_RESPELLED

    state = _open(monkeypatch, _maf(["chr1", "chr2"]))

    assert _CHROMOSOMES_RESPELLED not in state


def test_a_refused_file_promises_nothing_about_a_download(monkeypatch):
    """The sentence ends *"your download will spell them that way too"*.

    A MAF the app refuses produces no report and offers no download, so stashed before the
    refusal the notice is drawn beside an error, promising something about a file the app
    has just declined to open. Mutation: stash above `validate_required_columns` — where
    this first sat — and a refused file carries the promise.
    """
    from page_modules.data_loading import _CHROMOSOMES_RESPELLED

    state = _open(monkeypatch, _maf(["1", "2"]), refused=True)

    assert state["maf_data"] is None, "the fixture did not actually refuse the file"
    assert _CHROMOSOMES_RESPELLED not in state


def test_a_file_that_loads_and_fails_to_filter_still_says_it(monkeypatch):
    """The other side of the same line, and the reason it is not simply moved to the end.

    A re-spelling that happened is a re-spelling that happened: the file opened, the
    columns it changed are reachable, and only the *filter* fell over. Mutation: put the
    stash below the filter run and this goes silent on a file it changed.
    """
    from page_modules import data_loading
    from page_modules.data_loading import _CHROMOSOMES_RESPELLED

    monkeypatch.setattr(
        data_loading, "apply_filters_to_data", lambda show_messages=True: False
    )
    state = _open(monkeypatch, _maf(["1", "2"]))

    assert state[_CHROMOSOMES_RESPELLED] is True


def test_a_notice_never_drawn_does_not_survive_the_next_file(monkeypatch):
    """Every sibling in that tail is reset; this is the one a *render* consumes.

    So a bare-spelled MAF loaded and replaced with no render of the data page between the
    two would otherwise leave its flag standing, and the next file — already prefixed,
    changed in no way — would draw the sentence. Mutation: drop the `pop` from the reset
    tail and the prefixed file inherits the bare one's notice.
    """
    from page_modules.data_loading import _CHROMOSOMES_RESPELLED

    state = _open(
        monkeypatch, _maf(["chr1", "chr2"]), **{_CHROMOSOMES_RESPELLED: True}
    )

    assert _CHROMOSOMES_RESPELLED not in state


def _banners(monkeypatch, state):
    """What the page-level slot draws this render, and what it leaves behind."""
    from page_modules import data_loading

    drawn = []
    fake_st = MagicMock()
    fake_st.session_state = state
    fake_st.info.side_effect = lambda text, *a, **k: drawn.append(str(text))
    fake_st.warning.side_effect = lambda text, *a, **k: drawn.append(str(text))
    fake_st.error.side_effect = lambda text, *a, **k: drawn.append(str(text))
    monkeypatch.setattr(data_loading, "st", fake_st)
    monkeypatch.setattr(data_loading, "show_parameter_notice", lambda: None)
    monkeypatch.setattr(data_loading, "_show_arm_notice", lambda: None)
    monkeypatch.setattr(data_loading, "drain_missing_column_reports", lambda: [])

    data_loading._show_stashed_banners()
    return drawn


def test_the_notice_reaches_the_page_and_names_what_changed(monkeypatch):
    """Stashed by the load, drawn by whichever render draws the page.

    Mutation: delete the drain from ``_show_stashed_banners`` → nothing is ever said, on a
    file whose export the app has just re-spelled.
    """
    from tests.fakes import FakeSessionState
    from page_modules.data_loading import _CHROMOSOMES_RESPELLED

    state = FakeSessionState(**{_CHROMOSOMES_RESPELLED: True})
    drawn = _banners(monkeypatch, state)

    said = [text for text in drawn if "chromosomes" in text]
    assert len(said) == 1, f"the chromosome notice was said {len(said)} times: {drawn!r}"
    assert "download" in said[0], (
        "the notice does not mention the download — which is the reason it is said at "
        "all, the re-spelling being invisible on screen and permanent in the export"
    )


def test_the_notice_is_drained_rather_than_redrawn(monkeypatch):
    """Nothing is left behind for the next render.

    Mutation: read the key without popping it → the sentence reappears on every render of
    the page for the rest of the session, long after the load it describes.
    """
    from tests.fakes import FakeSessionState
    from page_modules.data_loading import _CHROMOSOMES_RESPELLED

    state = FakeSessionState(**{_CHROMOSOMES_RESPELLED: True})
    _banners(monkeypatch, state)

    assert _CHROMOSOMES_RESPELLED not in state
    assert not [text for text in _banners(monkeypatch, state) if "chromosomes" in text]


def test_a_prefixed_file_draws_no_notice(monkeypatch):
    """The state the overwhelming majority of users are in."""
    from tests.fakes import FakeSessionState

    drawn = _banners(monkeypatch, FakeSessionState())

    assert not [text for text in drawn if "chromosomes" in text], drawn
