"""The gene-list contract: how a paste or an upload becomes the symbols we filter on.

This is the unit suite for ``filters/gene_lists.py`` — the tokeniser, the token
validation and the named panels (issue #34). It owns behaviour the parity harness cannot
judge, because the pipeline has no counterpart for it: ``bin/filter_variants.py`` is
handed a *path* by Nextflow and never sees a human typing symbols into a box.

The one claim that *is* a parity claim — that a one-symbol-per-line paste selects the
same rows the pipeline selects from the same file — is asserted in
``tests/parity/contract.py``'s ``paste`` cases, not here. Here we assert the tokenising,
the dropping, and what the user is told.

Why any of this is worth a suite of its own
-------------------------------------------
Two bugs closed here were silent in opposite directions:

* the multi-line box split on commas only, so a one-symbol-per-line paste became a
  single token matching nothing and the report came back empty. Not hypothetical: the
  project's own gene files are *named* for commas and *delimited* by newlines, which
  ``test_the_projects_own_gene_files_are_newline_delimited`` pins.
* matching was case-sensitive, so a lowercase paste silently emptied the report.

And one crash class: a gene file whose tokens are not symbols makes the vendored clause
raise, because ``pd.read_csv`` infers a non-string dtype and ``.str.upper()`` has
nothing to work on. Those are reproduced here as inputs to the *tokeniser*, so the
assertion is that they are unreachable rather than that they are handled.
"""

from __future__ import annotations

import io
import os
import sys
from pathlib import Path

import pandas as pd
import pytest

STREAMLIT_APP = Path(__file__).resolve().parent.parent
if str(STREAMLIT_APP) not in sys.path:
    sys.path.insert(0, str(STREAMLIT_APP))

from config.gene_panels import GENE_LISTS, PREDEFINED_GENE_SETS  # noqa: E402
from filters.gene_lists import (  # noqa: E402
    BOOL_TOKENS,
    GENE_LIST_EXTENSIONS,
    HEADER_TOKENS,
    NA_TOKENS,
    GeneSelection,
    missing_symbols,
    panel_symbols,
    parse_gene_list,
)
from filters.variant_filters import _gene_file, apply_filters  # noqa: E402
from tests.fakes import note_texts  # noqa: E402

FIXTURE_DIR = STREAMLIT_APP / "tests" / "fixtures" / "parity"


# ---------------------------------------------------------------------------
# One tokeniser, every separator a user actually produces
# ---------------------------------------------------------------------------


def test_a_one_symbol_per_line_paste_tokenises_to_its_symbols():
    """The headline bug: newline-delimited input used to become one token.

    ``"APC\\nTP53"`` split on commas alone is the single token ``"APC\\nTP53"``, which
    matches no ``Hugo_Symbol`` anywhere — so the report came back empty and nothing said
    why.
    """
    assert parse_gene_list("APC\nTP53\nBRCA1").symbols == ("APC", "TP53", "BRCA1")


@pytest.mark.parametrize(
    "raw",
    [
        "APC,TP53,BRCA1",
        "APC;TP53;BRCA1",
        "APC TP53 BRCA1",
        "APC\nTP53\nBRCA1",
        "APC\tTP53\tBRCA1",
        "APC, TP53 ; BRCA1",
        "APC,\nTP53,\r\nBRCA1\n",
        "  APC  \n\n  TP53,;BRCA1  ",
    ],
    ids=[
        "commas",
        "semicolons",
        "spaces",
        "newlines",
        "tabs",
        "mixed with padding",
        "crlf and trailing commas",
        "blank lines and doubled separators",
    ],
)
def test_every_separator_a_user_produces_yields_the_same_three_symbols(raw):
    """One tokeniser, and it does not care which of these the user used.

    Parametrised over the forms rather than asserted for one, because the bug was
    precisely that the code handled one form and silently mishandled the rest.
    """
    assert parse_gene_list(raw).symbols == ("APC", "TP53", "BRCA1")


def test_a_typed_paste_and_an_uploaded_file_go_through_the_same_tokeniser():
    """Same bytes, same symbols — whichever box they arrived in.

    The acceptance criterion is "one tokeniser serves both", and the way to assert *one*
    is to hand it the two shapes the two widgets produce: ``st.text_area`` gives a
    ``str``, ``st.file_uploader`` gives bytes that the caller decodes. Both land here.
    """
    text = (FIXTURE_DIR / "genes_somatic.txt").read_text()
    typed = parse_gene_list(text)
    uploaded = parse_gene_list(text.encode("utf-8").decode("utf-8"))
    assert typed == uploaded
    assert len(typed.symbols) == 12


def test_the_projects_own_gene_files_are_newline_delimited():
    """The premise behind the headline bug, pinned to the files themselves.

    The widget label said "comma-separated" and the parser obeyed it, while every gene
    file this project ships is one symbol per line with no comma in it. That mismatch is
    what made a comma-only parser a real bug rather than a theoretical one, so it is
    asserted rather than described.
    """
    for name in ("genes_somatic.txt", "genes_germline.txt", "genes_somatic_mixed_case.txt"):
        text = (FIXTURE_DIR / name).read_text()
        assert "," not in text, f"{name} gained a comma; this test's premise has moved"
        assert len(text.splitlines()) > 1, f"{name} is not multi-line"
        # The whole point: comma-splitting yields one token, the tokeniser yields many.
        assert len([t for t in text.split(",") if t.strip()]) == 1
        assert len(parse_gene_list(text).symbols) > 1


def test_a_single_symbol_stays_a_one_element_tuple():
    """Not a bare string — which would be iterated character by character downstream.

    ``_gene_file`` writes ``"\\n".join(symbols)``, so handing it ``"BRCA1"`` instead of
    ``["BRCA1"]`` writes five one-character gene tokens and filters on none of them.
    """
    selection = parse_gene_list("BRCA1")
    assert selection.symbols == ("BRCA1",)
    assert not isinstance(selection.symbols, str)


def test_duplicates_collapse_case_insensitively_keeping_the_first_spelling():
    """One restriction, written once.

    The vendored clause upper-cases both sides, so ``APC`` and ``apc`` are the same
    restriction; writing both would only make the temp file longer and report the same
    gene as absent twice.
    """
    assert parse_gene_list("APC, apc, Apc, TP53").symbols == ("APC", "TP53")


# ---------------------------------------------------------------------------
# The leading header token
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("header", sorted(HEADER_TOKENS))
def test_a_leading_header_token_is_dropped(header):
    """A pasted spreadsheet column brings its heading with it."""
    selection = parse_gene_list(f"{header}\nAPC\nTP53")
    assert selection.symbols == ("APC", "TP53")
    assert selection.header == header


@pytest.mark.parametrize("header", ["Gene", "HUGO_SYMBOL", "Hugo_Symbol"])
def test_the_header_match_ignores_case(header):
    assert parse_gene_list(f"{header}\nAPC").symbols == ("APC",)


def test_only_the_leading_token_is_treated_as_a_header():
    """A header word further down is data, not a heading — and is kept as such.

    Guards the obvious over-reach: dropping every occurrence would silently delete a
    symbol from the middle of a list, which is the narrowing direction and invisible.
    """
    selection = parse_gene_list("APC\ngene\nTP53")
    assert selection.symbols == ("APC", "gene", "TP53")
    assert selection.header is None


def test_a_real_symbol_in_first_position_is_not_mistaken_for_a_header():
    """The header drop must cost nothing on the overwhelmingly common input."""
    selection = parse_gene_list((FIXTURE_DIR / "genes_somatic.txt").read_text())
    assert selection.header is None
    assert selection.symbols[0] == "APC"


def test_a_header_only_paste_restricts_nothing():
    """Dropping the header must not leave a one-token list of nothing."""
    selection = parse_gene_list("Hugo_Symbol")
    assert selection.symbols == ()
    assert not selection.restricts


# ---------------------------------------------------------------------------
# Token validation, and the crash paths it closes
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "token", ["123", "0", "-", ".", "*", "1.5", "---", "  7  "]
)
def test_a_token_with_no_letter_is_dropped(token):
    """No gene symbol is letterless, and a letterless token is what crashes the clause."""
    selection = parse_gene_list(f"APC {token} TP53")
    assert selection.symbols == ("APC", "TP53")
    assert token.strip() in selection.rejected


def test_dropped_tokens_are_named_to_the_user():
    """"Dropped and named" — a silent drop is how a typo becomes an unexplained result."""
    selection = parse_gene_list("APC, 123, TP53, 456")
    assert selection.rejected == ("123", "456")
    message = " ".join(note.text for note in selection.messages())
    assert "123" in message and "456" in message


@pytest.mark.parametrize("token", sorted(NA_TOKENS))
def test_the_tokens_pandas_reads_as_missing_are_dropped_too(token):
    """The letter rule alone is not enough, and this is the measured reason why.

    ``NA`` contains letters, so a letter test keeps it — and then ``pd.read_csv`` turns
    it into ``NaN``. A file of nothing but ``NA`` therefore parses to a *float* column
    and ``.str.upper()`` raises ``AttributeError``, exactly like the all-numeric file.
    Measured, not argued: ``test_no_token_in_the_corpus_can_crash_the_vendored_clause``
    re-derives it from the vendored code for every token this list names.
    """
    assert token not in parse_gene_list(f"APC\n{token}").symbols


def test_our_na_token_list_still_matches_what_pandas_actually_does():
    """The premise behind ``NA_TOKENS``, asserted against pandas rather than trusted.

    ``NA_TOKENS`` is a copy of pandas' default missing-value vocabulary, which the
    vendored ``pd.read_csv(header=None)`` call uses and cannot be told not to. A copy is
    only worth anything while it is provably still a copy — the same house rule that
    gates ``vendor/`` — so this reads every token back through ``read_csv`` and asserts
    pandas really does call it missing. A pandas release that *adds* a token is caught by
    ``test_every_accepted_token_reads_back_as_a_string``; one that *removes* one is caught
    here.
    """
    frame = pd.read_csv(io.StringIO("\n".join(sorted(NA_TOKENS))), header=None)
    assert frame[0].isna().all(), (
        "pandas no longer reads these as missing: "
        f"{sorted(frame.loc[frame[0].notna(), 0])} — NA_TOKENS has drifted"
    )


@pytest.mark.parametrize("token", sorted(BOOL_TOKENS))
def test_every_bool_token_really_would_have_broken_the_clause(token):
    """The premise behind ``BOOL_TOKENS`` — the other copied list, and the other direction.

    ``NA_TOKENS`` is guarded by the test above; this is its sibling, and it exists because
    without it ``BOOL_TOKENS`` would be an unguarded copy of a pandas behaviour, which is
    exactly what ``vendor/README.md``'s "derive it or guard it, never copy it" forbids.
    The corpus test cannot cover these: it *skips* every token the tokeniser rejects, so a
    ``BOOL_TOKENS`` entry that pandas no longer coerces would be silently over-rejected
    forever.

    The assertion is the consequence rather than the dtype name: a one-token file of each
    of these must be something ``.str.upper()`` **cannot** work on. If pandas ever stops
    inferring bool here, this fails and the entry can be dropped — over-rejection is safe
    but it costs a user a gene they are entitled to filter on.
    """
    column = pd.read_csv(io.StringIO(token + "\n"), header=None)[0]
    with pytest.raises(AttributeError):
        column.str.upper()


#: Inputs that make the vendored gene clause raise when written to a file verbatim. Each
#: was reproduced against ``vendor.pipeline_filters.somatic_filters`` before the
#: tokeniser existed; the exception each one raised is named so a change of failure mode
#: is visible.
CRASHING_INPUTS = {
    "all-numeric tokens": ("123\n456\n", "AttributeError"),
    "one numeric token": ("7\n", "AttributeError"),
    "empty": ("", "EmptyDataError"),
    "whitespace only": ("\n\n  \n", "EmptyDataError"),
    "NA only": ("NA\nNA\n", "AttributeError"),
    # Not a dtype failure: read_csv's default quotechar makes an unbalanced quote run off
    # the end of the file looking for its partner. Found by review, after the four above.
    "unbalanced quote": ('"BRCA1\n', "ParserError"),
}


def _vendored_gene_mask(frame: pd.DataFrame, symbols) -> pd.Series:
    """The vendored somatic clause's gene mask, reached the way the app reaches it."""
    from vendor.pipeline_filters import somatic_filters

    with _gene_file(list(symbols)) as path:
        mask, _ = somatic_filters(
            frame,
            vaf=-1.0,
            somatic_genes=path,
            cancervar_keep=["Tier_I_strong"],
            civic_keep=[],
            escat_keep=[],
            clinvar_keep=[],
            skip_civic=True,
            skip_pathogenic=True,
        )
    return mask


@pytest.fixture(scope="module")
def tiny_maf() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "Hugo_Symbol": ["APC", "TP53", "BRCA1"],
            "CancerVar": ["Tier_I_strong"] * 3,
            "ClinVar_VCF_CLNSIG": ["Pathogenic"] * 3,
            "ESCAT": ["IA"] * 3,
            "tumor_f": [0.5, 0.5, 0.5],
        }
    )


@pytest.mark.parametrize("label", sorted(CRASHING_INPUTS))
def test_the_reproduced_crash_paths_are_unreachable_through_the_tokeniser(
    label, tiny_maf
):
    """Each input that crashed the clause now reaches it as "no gene filter".

    The assertion is *unreachability*, not error handling: after tokenising there is
    nothing left to write, so the adapter hands over the ``"null"`` sentinel and the
    clause is never asked to parse a file at all. The mask is therefore all-True — the
    widening direction, which is the safe one and which
    ``test_all_invalid_input_leaves_the_report_unfiltered_and_says_so`` pairs with a
    warning.
    """
    raw, _ = CRASHING_INPUTS[label]
    selection = parse_gene_list(raw)
    assert not selection.restricts, f"{label!r} still produced a gene restriction"
    assert bool(_vendored_gene_mask(tiny_maf, selection.symbols).all())


def test_the_crash_paths_really_do_crash_without_the_tokeniser(tiny_maf):
    """The premise: these inputs are not merely ugly, they raise.

    Without this the test above could pass against inputs that were never dangerous, and
    the criterion it claims to meet would be unfalsifiable. Written verbatim to a file
    the way the pre-#34 code would have, and each named exception asserted.
    """
    import tempfile

    for label, (raw, expected) in sorted(CRASHING_INPUTS.items()):
        handle = tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False)
        handle.write(raw)
        handle.close()
        try:
            with pytest.raises(Exception) as caught:
                from vendor.pipeline_filters import somatic_filters

                somatic_filters(
                    tiny_maf,
                    vaf=-1.0,
                    somatic_genes=handle.name,
                    cancervar_keep=["Tier_I_strong"],
                    civic_keep=[],
                    escat_keep=[],
                    clinvar_keep=[],
                    skip_civic=True,
                    skip_pathogenic=True,
                )
            assert type(caught.value).__name__ == expected, (
                f"{label}: raised {type(caught.value).__name__}, not {expected} — the "
                "vendored failure mode moved"
            )
        finally:
            os.unlink(handle.name)


#: A corpus spanning every token class pandas has an opinion about, plus real symbols
#: whose spelling flirts with those classes. The point is breadth: a user's paste is not
#: drawn from a list of five, and each of the four validation conditions in
#: ``_could_be_a_symbol`` has to be witnessed by something here or it is untested.
TOKEN_CORPUS = [
    # real symbols, including ones that look dangerous and are not
    "APC", "TP53", "HLA-A", "C1orf43", "MT-CO1", "NANS", "NANOS1", "AFF1", "TERT",
    "e", "E", "D", "J", "T", "F", "Y", "N", "AF", "NAN", "Nan", "yes", "no", "on",
    # numeric literals
    "123", "0", "-1", "1.5", "1e5", "1E5", "-2.5e-3", "1_000",
    # infinities, which have letters
    "inf", "-inf", "INF", "Infinity", "infinity",
    # pandas' missing-value vocabulary
    "NA", "N/A", "n/a", "nan", "NaN", "-nan", "null", "NULL", "None", "<NA>",
    "#N/A", "#NA", "1.#IND", "-1.#QNAN",
    # booleans
    "true", "True", "TRUE", "false", "False", "FALSE",
    # letterless junk
    "-", ".", "***", "/", "?", "()",
    # things that merely look numeric to a human
    "0x1F", "1j", "5'UTR", "p.V600E",
    # quote characters: read_csv's default quotechar, so an unbalanced one never finishes
    # the parse at all — the one crash path that is not a dtype problem
    '"BRCA1', 'BRCA1"', '"BRCA1"', 'BR"CA1', "'BRCA1",
]


@pytest.mark.parametrize("token", TOKEN_CORPUS)
def test_every_accepted_token_reads_back_as_a_string(token):
    """The validation rule, re-derived from pandas instead of trusted.

    ``_could_be_a_symbol`` is a hand-written statement of "pandas will give this back as
    a string", and a hand-written statement of a library's behaviour is exactly the kind
    of copy this project guards rather than trusts — the first draft of the rule was
    "contains a letter", which let ``1e5`` and ``true`` through to crash the clause.

    So the check runs the real inference: each accepted token is written *alone* to a
    one-column CSV, read the way the vendored clause reads its gene file, and required to
    come back as something ``.str.upper()`` can work on. Alone, because a bad token
    beside a good one is harmless — the column stays ``object`` — and it is the all-bad
    file that changes the inferred dtype and raises.

    A pandas release that widens its inference fails here, on the machine of whoever
    upgraded, rather than in a clinician's report.
    """
    if not parse_gene_list(token).restricts:
        pytest.skip(f"{token!r} is rejected by the tokeniser, so it is never written")

    column = pd.read_csv(io.StringIO(token + "\n"), header=None)[0]
    assert column.dtype == object, (
        f"the tokeniser accepted {token!r} but pandas reads it back as "
        f"{column.dtype} — writing it to the gene file makes the vendored "
        "`.str.upper()` raise. Add the class to _could_be_a_symbol."
    )
    column.str.upper()  # the vendored call itself; raises if the dtype is wrong


@pytest.mark.parametrize("token", TOKEN_CORPUS)
def test_no_token_in_the_corpus_can_crash_the_vendored_clause(token, tiny_maf):
    """And end to end: tokenise, write, call the real vendored function.

    The previous test checks the dtype premise; this one checks the consequence through
    the adapter and the vendored code, so a mistake in the adapter's *writing* — the
    one-character-tokens collapse, a stray quote, a missing newline — cannot pass by
    satisfying the premise alone.
    """
    selection = parse_gene_list(token)
    mask = _vendored_gene_mask(tiny_maf, selection.symbols)
    assert mask.dtype == bool, f"{token!r} produced a non-boolean mask"
    if not selection.restricts:
        assert bool(mask.all()), f"{token!r} was rejected but still restricted rows"


# ---------------------------------------------------------------------------
# All-invalid input: unfiltered, and said out loud
# ---------------------------------------------------------------------------


def test_all_invalid_input_leaves_the_report_unfiltered_and_says_so():
    """Extra rows are visible; missing rows are not.

    So when nothing usable was given, the gene filter is dropped rather than applied to
    an empty list — which would answer with an empty report and look like a clinical
    finding. The warning is what stops the widening from being silent.
    """
    selection = parse_gene_list("123, 456, ---")
    assert selection.symbols == ()
    assert not selection.restricts
    assert selection.rejected == ("123", "456", "---")

    message = " ".join(note.text for note in selection.messages())
    assert "no gene filter" in message.lower(), (
        f"all-invalid input produced no warning naming the consequence: {message!r}"
    )


def test_an_all_invalid_gene_list_warns_through_the_filter_seam():
    """The warning has to reach the report, not just the parameter page.

    The user can tokenise on one page and filter on another, so the message rides
    ``Diagnostics.warnings`` — the channel ``data_loading`` already renders — rather than
    living only where the text was typed.
    """
    from config.pipeline_params import pipeline_params
    from vendor.pipeline_utils import read_maf

    frame = read_maf(str(FIXTURE_DIR / "somatic_reference.maf"))
    params = pipeline_params("somatic")
    params["filter_genes"] = ["123", "456"]

    labelled, diagnostics = apply_filters(frame, params)
    unrestricted, _ = apply_filters(frame, pipeline_params("somatic"))

    from filters.variant_filters import MAFIGATE_FILTER

    assert list(labelled[MAFIGATE_FILTER]) == list(unrestricted[MAFIGATE_FILTER]), (
        "an all-invalid gene list changed the verdicts — it must be inert"
    )
    joined = " ".join(note_texts(diagnostics))
    assert "gene" in joined.lower() and "123" in joined


# ---------------------------------------------------------------------------
# Genes the MAF does not have
# ---------------------------------------------------------------------------


def test_symbols_absent_from_the_maf_are_named():
    """A typo'd symbol is indistinguishable from a gene with no variants — unless we say.

    Case-insensitively, because that is how the match itself works: a lowercase paste
    that the vendored clause happily matches must not be reported as absent.
    """
    present = ["APC", "TP53"]
    assert missing_symbols(["APC", "BRCA1"], present) == ("BRCA1",)
    assert missing_symbols(["apc", "tp53"], present) == ()
    assert missing_symbols([], present) == ()


def test_the_filter_seam_reports_genes_the_maf_does_not_carry():
    from config.pipeline_params import pipeline_params
    from vendor.pipeline_utils import read_maf

    frame = read_maf(str(FIXTURE_DIR / "somatic_reference.maf"))
    params = pipeline_params("somatic")
    params["filter_genes"] = ["APC", "NOTAREALGENE"]

    _, diagnostics = apply_filters(frame, params)
    joined = " ".join(note_texts(diagnostics))
    assert "NOTAREALGENE" in joined, (
        f"a requested gene absent from the MAF went unreported: {note_texts(diagnostics)}"
    )
    assert "APC" not in joined.replace("NOTAREALGENE", ""), (
        "a gene the MAF does carry was reported as absent"
    )


def test_a_gene_list_the_maf_matches_entirely_produces_no_absence_warning():
    from config.pipeline_params import pipeline_params
    from vendor.pipeline_utils import read_maf

    frame = read_maf(str(FIXTURE_DIR / "somatic_reference.maf"))
    params = pipeline_params("somatic")
    params["filter_genes"] = ["apc"]  # lowercase on purpose: the match is insensitive

    _, diagnostics = apply_filters(frame, params)
    assert not [t for t in note_texts(diagnostics) if "not present" in t.lower()], (
        f"a matching gene was reported absent: {note_texts(diagnostics)}"
    )


# ---------------------------------------------------------------------------
# Case-insensitivity, and where it comes from
# ---------------------------------------------------------------------------


def test_the_tokeniser_does_not_normalise_case_itself():
    """Case-insensitivity is achieved by *routing*, which the criterion is explicit about.

    Upper-casing in app code would work today and rot tomorrow: it would be a second
    normalisation to keep in step with the vendored ``.str.upper()`` on both sides. So
    the tokeniser preserves what the user typed, and the vendored clause is what makes
    the comparison case-blind. Asserting the symbols come back unchanged is how "we did
    not normalise here" is made checkable.
    """
    assert parse_gene_list("brca1, Tp53").symbols == ("brca1", "Tp53")


# ---------------------------------------------------------------------------
# The named panels
# ---------------------------------------------------------------------------


def test_every_offered_panel_resolves_to_a_gene_list():
    """Every option the selectbox shows must resolve, or the dropdown lies."""
    for option in GENE_LISTS:
        symbols = panel_symbols(option)
        assert isinstance(symbols, tuple), option
        assert all(isinstance(s, str) for s in symbols), option


@pytest.mark.parametrize("option", sorted(PREDEFINED_GENE_SETS))
def test_a_named_panel_resolves_to_its_own_symbols(option):
    resolved = panel_symbols(option)
    assert resolved == parse_gene_list(PREDEFINED_GENE_SETS[option]).symbols
    assert len(resolved) > 100, f"{option} resolved to only {len(resolved)} symbols"


@pytest.mark.parametrize("option", ["All", "Custom"])
def test_the_neutral_panel_options_restrict_nothing(option):
    """"All" is no filter, and "Custom" defers to the boxes — neither restricts by itself."""
    assert panel_symbols(option) == ()


def test_the_default_panel_is_neutral_at_parity():
    """The dropdown's default must contribute exactly the contract's own gene value.

    The parity cases run with ``filter_genes == []``, so if the panel's default resolved to
    anything else, every parity number in ``baseline.json`` would describe a configuration
    no user of the page can actually produce. Both arms, because the page renders per arm
    and a fix landing on one used to be how these two drifted.
    """
    from config.pipeline_params import pipeline_params

    default_panel = next(iter(GENE_LISTS))
    assert default_panel == "All", (
        f"the dropdown's first option is now {default_panel!r}; a fresh session would "
        "open with a gene restriction the parity baseline was not measured under"
    )
    for arm in ("somatic", "germline"):
        assert list(panel_symbols(default_panel)) == pipeline_params(arm)["filter_genes"]


def test_an_unknown_panel_restricts_nothing_rather_than_raising():
    """A stale saved parameter file must not take the page down.

    Widening, again, is the direction that is visible to the user; a traceback on the
    parameter page is not.
    """
    assert panel_symbols("A_PANEL_WE_RETIRED") == ()


def test_the_panel_choice_is_never_a_filter_parameter():
    """The selectbox is UI state, and the source is the only place that can be checked.

    ``config/pipeline_params.py`` says it: "Gene *set* selection (the panel dropdown) is
    UI state that resolves to this list, never a parameter of its own", and
    ``test_param_contract.py`` already asserts the contract carries no such key. What was
    unasserted is the *writer* — the page used to store ``somatic_gene_set`` straight
    into ``filter_params``, which made the panel name a filter parameter in every saved
    file and every cache entry. Read out of the source because the alternative is booting
    Streamlit to observe a dict.
    """
    source = (STREAMLIT_APP / "page_modules" / "parameter_config.py").read_text()
    for stale in (
        '"somatic_gene_set"',
        '"germline_gene_set"',
        '"somatic_genes"',
        '"germline_genes"',
    ):
        assert stale not in source, (
            f"page_modules/parameter_config.py still writes {stale} — the panel name or "
            "a per-arm gene string is being stored as a filter parameter. The panel is "
            "UI state; the one filter parameter is `filter_genes`."
        )
    assert '"filter_genes"' in source, (
        "the parameter page no longer writes `filter_genes` — the gene boxes have "
        "stopped reaching the filter altogether"
    )


# ---------------------------------------------------------------------------
# The page itself, driven
# ---------------------------------------------------------------------------

#: A script that renders just the Basic Filters tab, which is where the gene controls
#: live. Rendering the tab rather than the whole page keeps the run to the widgets under
#: test and avoids the parameter cache, which otherwise decides the arm from whatever the
#: developer last configured (see ``tests/run_app_check.py``).
_GENE_TAB_SCRIPT = """
import sys
sys.path.insert(0, {app!r})
import streamlit as st
from config.pipeline_params import pipeline_params
from page_modules.parameter_config import show_basic_filters_tab

if "filter_params" not in st.session_state:
    st.session_state.filter_params = pipeline_params("somatic")
{seed}
show_basic_filters_tab("somatic")
"""


def _render_gene_tab(seed=""):
    """The Basic Filters tab under ``AppTest``, over an optionally pre-seeded parameter dict.

    ``seed`` is a line of script run before the tab renders, which is how the parameter cache
    and an imported file both behave: something else decides ``filter_params`` and the page
    is rendered afterwards. Inserted into the script rather than set through
    ``AppTest.session_state`` so the seeding happens *inside* the run, as it really does.
    """
    pytest.importorskip("streamlit", reason="streamlit not installed")
    from streamlit.testing.v1 import AppTest

    app = AppTest.from_string(
        _GENE_TAB_SCRIPT.format(app=str(STREAMLIT_APP), seed=seed)
    )
    app.run()
    assert not app.exception, [str(e.value) for e in app.exception]
    return app


@pytest.fixture
def gene_tab():
    """The Basic Filters tab, rendered by ``streamlit.testing.v1.AppTest``."""
    return _render_gene_tab()


def test_the_page_writes_filter_genes_and_nothing_else_gene_shaped(gene_tab):
    """The wiring, observed rather than grepped.

    ``test_the_panel_choice_is_never_a_filter_parameter`` reads the source and can only
    prove an absence. This runs the real widgets and proves the presence: after a render,
    the parameter dict carries a ``filter_genes`` *list* and no panel name.
    """
    params = gene_tab.session_state["filter_params"]
    assert params["filter_genes"] == [], "the default panel must restrict nothing"
    assert isinstance(params["filter_genes"], list)
    assert not [key for key in params if key.endswith("_gene_set")], (
        f"the page stored a panel name as a filter parameter: {sorted(params)}"
    )


def test_a_paste_into_the_page_reaches_filter_genes_tokenised(gene_tab):
    """The end-to-end form of the headline fix: type it, and the filter gets symbols.

    Lowercase and newline-delimited on purpose — the two bugs this ticket closes — plus a
    letterless token, so the drop and the warning are observed on the page too.
    """
    gene_tab.selectbox(key="param_somatic_gene_set").select("Custom").run()
    gene_tab.text_area(key="gene_text_somatic").input("apc\ntp53\n123").run()
    assert not gene_tab.exception, [str(e.value) for e in gene_tab.exception]

    assert gene_tab.session_state["filter_params"]["filter_genes"] == ["apc", "tp53"]
    warnings = " ".join(w.value for w in gene_tab.warning)
    assert "123" in warnings, f"the dropped token was not named on the page: {warnings}"


def test_choosing_a_named_panel_on_the_page_resolves_it_to_symbols(gene_tab):
    """The selectbox still works, and resolves to a list rather than to its own name."""
    gene_tab.selectbox(key="param_somatic_gene_set").select("MSK_IMPACT").run()
    assert not gene_tab.exception, [str(e.value) for e in gene_tab.exception]

    resolved = gene_tab.session_state["filter_params"]["filter_genes"]
    assert resolved == list(panel_symbols("MSK_IMPACT"))
    assert len(resolved) > 100
    assert not [key for key in gene_tab.session_state["filter_params"] if key.endswith("_gene_set")]


def _tab_with_incoming_genes(genes):
    """The tab rendered over a ``filter_params`` that already carries a gene list.

    Which is what the ``~/.mafigate`` cache and an imported JSON/YAML both produce — the
    page is rendered *after* someone else has decided what ``filter_genes`` is.
    """
    return _render_gene_tab(
        f"st.session_state.filter_params['filter_genes'] = {genes!r}"
    )


def test_an_incoming_gene_list_is_adopted_rather_than_discarded():
    """A saved gene list must survive being rendered. It used to be silently thrown away.

    The regression this pins is the widening kind issue #28 says must never be silent. The
    widgets own their state by key, so a ``filter_genes`` restored from the parameter cache
    reached a dropdown that had never heard of it, and the page then overwrote the real list
    with ``[]``: measured at ``['TP53', 'BRCA1']`` in, ``[]`` out, and the page reading
    "Using all genes" while the user believed two genes were selected.
    """
    app = _tab_with_incoming_genes(["TP53", "BRCA1"])

    assert app.session_state["filter_params"]["filter_genes"] == ["TP53", "BRCA1"], (
        "the page discarded an incoming gene list — silent widening"
    )
    assert app.selectbox(key="param_somatic_gene_set").value == "Custom"
    assert set(app.text_area(key="gene_text_somatic").value.split()) == {
        "TP53",
        "BRCA1",
    }


def test_an_incoming_list_that_is_a_named_panel_is_adopted_as_that_panel():
    """A saved panel choice round-trips, though only its symbols were ever stored.

    The panel name is deliberately not a parameter, so the only way back to the dropdown is
    to recognise the symbols. Without this the user's MSK-IMPACT selection would reappear as
    505 symbols pasted into a Custom box — the same filter, but not the same page.
    """
    app = _tab_with_incoming_genes(list(panel_symbols("MSK_IMPACT")))

    assert app.selectbox(key="param_somatic_gene_set").value == "MSK_IMPACT"
    assert app.session_state["filter_params"]["filter_genes"] == list(
        panel_symbols("MSK_IMPACT")
    )


def test_an_incoming_legacy_per_arm_gene_string_is_adopted_too():
    """The shape a pre-#34 saved file holds: a comma string under a per-arm key.

    ``_gene_selection`` reads it at the filter seam; the page has to show it as well, or the
    user sees "all genes" while the filter restricts.
    """
    app = _render_gene_tab(
        "st.session_state.filter_params.pop('filter_genes', None)\n"
        "st.session_state.filter_params['somatic_genes'] = 'TP53, BRCA1'"
    )
    assert app.session_state["filter_params"]["filter_genes"] == ["TP53", "BRCA1"]


def test_adopting_does_not_overwrite_what_the_user_is_typing():
    """The other half of the fix: adoption must not fight the typist.

    A page that re-adopted on every run would rewrite the box from its own parsed output,
    so a half-typed symbol — or a rejected token the user is about to correct — would vanish
    mid-keystroke. The stamp is what distinguishes a list the page just produced from one
    that arrived.
    """
    app = _tab_with_incoming_genes(["TP53"])
    app.text_area(key="gene_text_somatic").input("TP53\nBRCA1\n12").run()
    assert not app.exception, [str(e.value) for e in app.exception]

    assert app.text_area(key="gene_text_somatic").value == "TP53\nBRCA1\n12", (
        "the page rewrote the box out from under the user"
    )
    assert app.session_state["filter_params"]["filter_genes"] == ["TP53", "BRCA1"]


def test_switching_back_to_all_clears_the_gene_restriction(gene_tab):
    """A user who narrows and then widens must actually get the wider report back.

    The old page left the previous panel's symbols in a per-arm key that the filter still
    read, so "All" was not always all.
    """
    gene_tab.selectbox(key="param_somatic_gene_set").select("MSK_IMPACT").run()
    assert gene_tab.session_state["filter_params"]["filter_genes"]

    gene_tab.selectbox(key="param_somatic_gene_set").select("All").run()
    assert gene_tab.session_state["filter_params"]["filter_genes"] == []


# ---------------------------------------------------------------------------
# The upload widget's accepted extensions
# ---------------------------------------------------------------------------


def test_the_usual_list_file_extensions_are_accepted():
    """What a user's gene list is actually called on disk.

    Extensions without the dot, which is the shape ``st.file_uploader(type=...)`` wants;
    asserted so a reformat cannot quietly turn them into a list Streamlit rejects.
    """
    assert set(GENE_LIST_EXTENSIONS) >= {"txt", "csv", "tsv", "list"}
    assert not any(ext.startswith(".") for ext in GENE_LIST_EXTENSIONS)


# ---------------------------------------------------------------------------
# The adapter's own contract, from this ticket's side
# ---------------------------------------------------------------------------


def test_the_adapter_refuses_a_bare_string():
    """The collapse this ticket forbids, made loud instead of silent.

    ``"\\n".join("BRCA1")`` is ``"B\\nR\\nC\\nA\\n1"`` — five one-character gene tokens,
    a filter that matches nothing, and no error. A ``TypeError`` at the seam is strictly
    better than an empty report.
    """
    with pytest.raises(TypeError, match="list of symbols"):
        with _gene_file("BRCA1"):
            pass


def test_the_adapter_writes_one_line_per_symbol():
    with _gene_file(["BRCA1"]) as path:
        assert Path(path).read_text().splitlines() == ["BRCA1"]
    with _gene_file(["BRCA1", "TP53"]) as path:
        assert Path(path).read_text().splitlines() == ["BRCA1", "TP53"]


def test_the_adapter_removes_its_temp_file_on_the_error_path():
    """The leak that matters most is the one on the path nobody tests.

    ``test_filter_entry_point.py`` covers the happy path. A vendored call that raises —
    which is the pipeline's own behaviour on an unreadable depth column (#38) — must not
    also leak a file per click.
    """
    escaped: list[str] = []
    with pytest.raises(RuntimeError):
        with _gene_file(["BRCA1"]) as path:
            escaped.append(path)
            raise RuntimeError("as the vendored clause would")
    assert not os.path.exists(escaped[0]), f"the adapter leaked {escaped[0]} on error"


def test_a_selection_is_comparable_and_hashable():
    """``GeneSelection`` is a value: two equal pastes must be one thing.

    Relied on by ``test_a_typed_paste_and_an_uploaded_file_go_through_the_same_tokeniser``,
    which asserts the typed and uploaded paths agree by comparing two selections outright.
    """
    assert parse_gene_list("APC TP53") == parse_gene_list("APC,TP53")
    assert isinstance(hash(parse_gene_list("APC")), int)
    assert GeneSelection(symbols=("APC",)) == parse_gene_list("APC")
