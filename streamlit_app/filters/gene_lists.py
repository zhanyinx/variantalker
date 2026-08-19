"""How a pasted or uploaded gene list becomes the symbols the filter restricts to.

One tokeniser, used by every path a gene symbol can arrive on: the multi-line text box,
the uploaded list file, and the named panel dropdown. That is the whole point of the
module — the bug it closes was two parsers disagreeing about what a separator is.

Where this sits
---------------
``filters/variant_filters.py`` owns the *decision* and delegates all of it to
``vendor/``. This module owns the step before: turning what a human produced into the
list of symbols that adapter writes to a file. It is deliberately free of Streamlit and
of pandas, so the parameter page, the filter seam and the unit suite can all reach it.

Three things the tokeniser deliberately does not do
--------------------------------------------------
* **It does not normalise case.** Matching is case-insensitive, but that is achieved by
  *routing*: the vendored clause is
  ``maf["Hugo_Symbol"].str.upper().isin(genes[0].str.upper().values)``, which upper-cases
  both sides itself. Upper-casing here as well would be a second normalisation to keep in
  step with that one, and the failure mode of a stale copy is silent. So symbols come back
  spelled the way the user typed them.
* **It does not resolve panels into the parameter dict as a name.** ``panel_symbols``
  turns a panel choice into symbols; the choice itself is UI state and never a filter
  parameter. ``config/pipeline_params.py`` states this, and
  ``tests/test_gene_lists.py::test_the_panel_choice_is_never_a_filter_parameter`` checks
  the writer.
* **It does not guess.** A token that cannot be a gene symbol is dropped and *named*,
  never repaired.

Why letterless tokens are dropped rather than passed through
-----------------------------------------------------------
The vendored clause reads its file with ``pd.read_csv(path, header=None)`` and then calls
``genes[0].str.upper()``. ``read_csv`` infers the column's dtype, so a file whose tokens
are all numeric parses to ``int64`` and ``.str`` raises ``AttributeError``; an empty or
whitespace-only file raises ``EmptyDataError`` before that. Five such inputs were
reproduced against the vendored code, and ``tests/test_gene_lists.py`` keeps both halves
of the claim honest — that they really do raise, and that nothing surviving this
tokeniser can reach them.

The letter rule alone is not sufficient, and each shortfall was measured rather than
guessed — which is why the guard is four conditions and not one:

======================  =========  ==================================================
token                   inferred   why the letter rule misses it
======================  =========  ==================================================
``123``, ``1.5``        int/float  caught by the letter rule
``NA``, ``<NA>``        float64    has letters; pandas reads it as *missing*
``1e5``, ``inf``        float64    has letters; pandas reads it as a *number*
``true``, ``FALSE``     bool       has letters; pandas reads it as a *boolean*
======================  =========  ==================================================

Hence :data:`NA_TOKENS`, :data:`BOOL_TOKENS` and the ``float()`` probe alongside it.
``tests/test_gene_lists.py::test_every_accepted_token_reads_back_as_a_string`` derives
the whole table from pandas again on every run, so a pandas release that widens its
inference fails there rather than in a clinician's report.

The cost of the rule is stated plainly: a gene whose symbol is spelled exactly like a
number, a boolean or a missing-value marker cannot be filtered on. No HGNC symbol is —
and if one ever were, the vendored reader could not represent it as a string either, so
this is the pipeline's limit rather than the app's. The token is named to the user
either way.

The safe direction
------------------
When nothing usable is left, the answer is **no gene filter**, not an empty gene filter.
An empty list would produce an empty report, which reads like a clinical finding; a
dropped filter produces extra rows, which the user can see. The warning is what keeps
the widening from being silent.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Iterable

from config.gene_panels import PREDEFINED_GENE_SETS
from filters.notes import Note

#: Every separator a real paste contains: the commas and semicolons a spreadsheet export
#: uses, and the newlines and tabs everything else uses. Whitespace is in here because
#: the project's own gene files are *named* for commas and *delimited* by newlines — the
#: mismatch that made a comma-only parser a real bug and not a theoretical one.
_SEPARATORS = re.compile(r"[,;\s]+")

#: File extensions the upload widget accepts, without the leading dot — the shape
#: ``st.file_uploader(type=...)`` wants.
GENE_LIST_EXTENSIONS = ("txt", "csv", "tsv", "list", "genes")

#: Words that are a column heading rather than a gene symbol, matched case-insensitively
#: and **only in first position**. No HGNC symbol collides with any of these, which is
#: what makes dropping them safe; dropping every occurrence would not be, because that
#: would delete a symbol from the middle of a list in the narrowing direction.
#:
#: Single words only, deliberately: the tokeniser splits on whitespace, so a two-word
#: heading such as ``Gene Symbol`` arrives as ``gene`` followed by ``symbol`` and the
#: first of those is matched here. The second is then a letterful token that matches no
#: MAF row — reported as absent rather than dropped, which is the visible failure.
HEADER_TOKENS = frozenset(
    {
        "gene",
        "genes",
        "gene_id",
        "gene_name",
        "gene_symbol",
        "genename",
        "genesymbol",
        "hgnc",
        "hgnc_symbol",
        "hugo",
        "hugo_symbol",
        "symbol",
        "symbols",
    }
)

#: pandas' default missing-value vocabulary, as of pandas 2.x, minus the empty string
#: (which this module drops as an empty token before ever getting here).
#:
#: A copy, and copies rot — so it is guarded rather than trusted:
#: ``test_our_na_token_list_still_matches_what_pandas_actually_does`` reads every token
#: below back through ``pd.read_csv`` and asserts pandas really does call it missing. It
#: cannot be derived instead, because the vendored call takes no ``na_values`` argument
#: for us to read it off, and ``pandas._libs.parsers.STR_NA_VALUES`` is private.
#:
#: Matched exactly, not case-insensitively, because pandas matches exactly: this set
#: already carries the case variants pandas honours, and widening the match would drop
#: tokens pandas would have kept.
NA_TOKENS = frozenset(
    {
        "#N/A",
        "#N/A N/A",
        "#NA",
        "-1.#IND",
        "-1.#QNAN",
        "-NaN",
        "-nan",
        "1.#IND",
        "1.#QNAN",
        "<NA>",
        "N/A",
        "NA",
        "NULL",
        "NaN",
        "None",
        "n/a",
        "nan",
        "null",
    }
)

#: The literals ``pd.read_csv`` infers a **bool** column from. A one-token file of any of
#: these parses to ``dtype=bool``, and ``.str.upper()`` raises on that exactly as it does
#: on a numeric column. Measured, and guarded like :data:`NA_TOKENS` — by its own sibling
#: test, ``test_every_bool_token_really_would_have_broken_the_clause``, because the corpus
#: test skips whatever the tokeniser rejects and so cannot see these at all.
BOOL_TOKENS = frozenset({"TRUE", "True", "true", "FALSE", "False", "false"})

#: ``pd.read_csv``'s default ``quotechar``, which the vendored call does not override. A
#: token carrying one is the fourth way to break the read and the only one that is not a
#: dtype problem: ``"BRCA1`` — a symbol picked up from a quoted CSV with its opening quote
#: attached — makes the parser look for the closing quote and raise
#: ``ParserError: EOF inside string``. Measured; balanced quotes survive, but no HGNC symbol
#: contains one at all, so the whole character is refused rather than counted.
QUOTE_CHAR = '"'

#: How many tokens a warning names before it summarises the rest. A paste of a whole
#: spreadsheet column can reject hundreds, and a warning nobody reads is a warning that
#: does not exist.
_MAX_NAMED = 10


def name_tokens(tokens, limit: int = _MAX_NAMED) -> str:
    """Tokens as a quoted, comma-joined phrase for a warning, truncated with a count.

    Shared with ``variant_filters._gene_notes``, which names the symbols a MAF does not
    carry. Both lists come from the same paste and are read by the same user in the same
    place, so a second copy of the truncation would let one warning name ten items and the
    other all four hundred.
    """
    named = ", ".join(f"`{token}`" for token in tokens[:limit])
    if len(tokens) > limit:
        named += f" and {len(tokens) - limit} more"
    return named


@dataclass(frozen=True)
class GeneSelection:
    """The symbols to restrict to, what was thrown away, and what to tell the user.

    A frozen value rather than a bare list, because three separate facts come out of one
    parse and the caller needs all three: the symbols go to the filter, the rejects go to
    a warning, and ``restricts`` is what decides whether a gene filter happens at all.
    Returning only the list is what let the "nothing usable" case look identical to the
    "no gene filter requested" case.

    Frozen, so it is a value: two equal pastes are one thing regardless of which widget
    produced them, which is what lets the typed and uploaded paths be compared for equality
    rather than symbol by symbol.
    """

    #: The usable symbols, in the order given, deduplicated case-insensitively, spelled
    #: as the user spelled them. Never a bare string.
    symbols: tuple[str, ...] = ()
    #: Tokens that cannot be a gene symbol, in the order given. Named to the user.
    rejected: tuple[str, ...] = ()
    #: The leading token dropped as a column heading, if there was one.
    header: str | None = None

    @property
    def restricts(self) -> bool:
        """Whether this selection asks for any gene restriction at all.

        The question the adapter actually needs answered, and not the same as "was
        anything typed": an all-invalid paste asked for a restriction and gets none.
        """
        return bool(self.symbols)

    def messages(self) -> tuple[Note, ...]:
        """What the user needs to be told about this parse, in display order.

        Returned rather than rendered so this module stays Streamlit-free and so the same
        notes can ride ``Diagnostics.notes`` to the report page — the user can type a gene
        list on one page and filter on another, and a warning shown only where the typing
        happened is a warning they will not see next to the result.

        The levels, and where issue #151 put the line
        --------------------------------------------
        Dropping a heading is ``INFO``: the user pasted a column of symbols with its header
        still attached, MAFigate did the obvious thing, and the restriction they get is the one
        they asked for. Nothing is wrong, so nothing warns.

        **Both rejections warn**, and the partial one is the case the dev was asked about
        directly. An unusable token could never have matched a gene, so on the report's own
        terms nothing changed and it could have been drawn as another note about a run that
        did what was asked. It warns because the token the user typed is the evidence of what
        they *meant*: a mistyped ``TP-53`` is a gene the clinician believes they restricted to
        and did not, and that belief is not visible anywhere in the table. Which makes this the
        one note here whose level is about the user's intent rather than about the frame.
        """
        notes: list[Note] = []

        if self.header is not None:
            notes.append(
                Note.info(
                    f"ℹ️ Gene list: dropped the leading `{self.header}` as a column heading."
                )
            )

        if self.rejected:
            named = name_tokens(self.rejected)
            if self.restricts:
                notes.append(
                    Note.warning(
                        f"⚠️ Gene list: ignored {len(self.rejected)} entr"
                        f"{'y' if len(self.rejected) == 1 else 'ies'} that cannot be a gene "
                        f"symbol — {named}. "
                        f"Filtering on the remaining {len(self.symbols)}."
                    )
                )
            else:
                notes.append(
                    Note.warning(
                        f"⚠️ Gene list: none of the {len(self.rejected)} entries given "
                        f"could be a gene symbol — {named}. **No gene filter was "
                        "applied**, so the report is wider than you asked for rather "
                        "than narrower."
                    )
                )

        return tuple(notes)


def parse_gene_list(raw) -> GeneSelection:
    """The one tokeniser: whatever the user gave us, as the symbols we can filter on.

    Args:
        raw: ``None``, a string, or any iterable of values. A string is what
            ``st.text_area`` yields and what a decoded upload yields; an iterable is what
            a saved parameter file or a named panel yields. Every element is stringified
            and then tokenised, so a list of one comma-separated string and a list of
            symbols both work — the caller does not have to know which shape it holds.

    Returns:
        A :class:`GeneSelection`. Symbols keep their given order and spelling;
        duplicates collapse case-insensitively onto the first spelling seen, because the
        vendored clause upper-cases both sides and so a repeat is the same restriction
        written twice.
    """
    tokens = _tokens(raw)

    header: str | None = None
    if tokens and tokens[0].casefold() in HEADER_TOKENS:
        header, tokens = tokens[0], tokens[1:]

    symbols: list[str] = []
    rejected: list[str] = []
    seen: set[str] = set()
    for token in tokens:
        if not _could_be_a_symbol(token):
            if token not in rejected:
                rejected.append(token)
            continue
        if token.upper() in seen:
            continue
        seen.add(token.upper())
        symbols.append(token)

    return GeneSelection(
        symbols=tuple(symbols), rejected=tuple(rejected), header=header
    )


def panel_symbols(panel: str) -> tuple[str, ...]:
    """The symbols a named panel denotes, or ``()`` for the options that restrict nothing.

    An unrecognised name resolves to no restriction rather than raising: the name can
    come out of a saved parameter file written by an older version, and a traceback on
    the parameter page is less visible to a clinician than a wider report is.

    Routed through :func:`parse_gene_list` rather than returning the constant directly,
    so a panel and a paste of the same panel cannot disagree — and so a placeholder panel
    that is still empty comes back as ``()`` by the same route as ``All``.
    """
    if panel in PREDEFINED_GENE_SETS:
        return parse_gene_list(PREDEFINED_GENE_SETS[panel]).symbols
    return ()


def missing_symbols(requested: Iterable, present: Iterable) -> tuple[str, ...]:
    """Requested symbols that the MAF does not carry, in the order requested.

    Case-insensitive, because the match itself is: a lowercase paste that the vendored
    clause happily matches must not then be reported as absent. Upper-casing here is for
    the *report* and not for the decision — the decision is made inside ``vendor/``,
    which does its own upper-casing on both sides.

    Worth reporting because a symbol with no rows is indistinguishable from a typo: both
    just make the report smaller, and only one of them is what the user meant.
    """
    have = {value.upper() for value in present if isinstance(value, str)}
    return tuple(
        symbol for symbol in requested if str(symbol).upper() not in have
    )


# ---------------------------------------------------------------------------
# Internals
# ---------------------------------------------------------------------------


def _tokens(raw) -> list[str]:
    """``raw`` flattened to non-empty tokens, whatever shape it arrived in."""
    if raw is None:
        return []
    pieces = [raw] if isinstance(raw, str) else list(raw)
    tokens: list[str] = []
    for piece in pieces:
        tokens += [t for t in _SEPARATORS.split(str(piece)) if t]
    return tokens


def _could_be_a_symbol(token: str) -> bool:
    """Whether ``token`` could name a gene — the guard that closes the crash paths.

    One question underneath, asked five ways: *can* ``pd.read_csv`` hand the vendored clause
    this token back as a **string**? Four of the five are about dtype — it infers the dtype
    of the whole column, so anything it can read as a number, a boolean or a missing value
    takes the column with it and ``.str.upper()`` has nothing to work on. The fifth,
    :data:`QUOTE_CHAR`, is about the parse never finishing at all. The module docstring
    tabulates which condition catches which token.

    None of this claims the symbol *exists* — that is :func:`missing_symbols`'s job,
    answered against the MAF in hand. This is only the claim that writing the token to
    the vendored code's file cannot break the call.
    """
    return (
        any(char.isalpha() for char in token)
        and token not in NA_TOKENS
        and token not in BOOL_TOKENS
        and not _reads_as_a_number(token)
        and QUOTE_CHAR not in token
    )


def _reads_as_a_number(token: str) -> bool:
    """Whether ``token`` is a numeric literal — ``1e5``, ``inf``, ``-Infinity``.

    ``float()`` rather than a regex: it is the same grammar the C parser accepts, and a
    hand-written pattern for "number" is how ``1e5`` got through in the first place. It is
    marginally *stricter* than pandas — pandas keeps ``1_000`` and ``NAN`` as strings
    where ``float()`` takes them — and stricter is the harmless direction here: the token
    is dropped and named, which widens the report visibly rather than emptying it
    silently.
    """
    try:
        float(token)
    except ValueError:
        return False
    return True
