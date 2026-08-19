"""What the last filter run was told to do, in the words the controls use.

The results view says *"{n} variants passed all applied filters"* three times over and
never says what those filters were. This builds the answer: one line per setting the
filter read, named as the parameter page names it and valued as the filter resolved it.

**Read off the resolved settings, never the dict.** ``_Settings.from_params`` is what the
filter itself runs on, and it is not a formality — the germline VAF threshold is written
by the widget as ``vaf_threshold_germline`` and by the contract as ``vaf_threshold``, so a
recap reading the dict by key reports ``None`` for it on the germline arm. Every other
two-spelling case (``keep_pathogenic``/``skip_pathogenic``, the gene list's parse) has the
same shape. Describing the object the filter was handed is what makes a wrong description
unreachable rather than merely unlikely.

That is also why this is not the surface issue #28 deleted. The app used to carry eight
on-screen parameter echoes and they were **wrong**: they read a catch-all sentinel that no
longer meant anything, the superseded germline VAF key, and the classification list as
though it were an include list. Each was a hand-written reading of the dict, made
somewhere the filter was not. This one cannot drift from the run that way, because it does
not do its own reading.

**A snapshot, not a live view.** :func:`describe_run` is called once, by the filter run,
and the result is carried with the report. Reading the current parameters instead would
quietly misdescribe the table beside it the moment a user changes a setting without
re-filtering — the state ``data_params_hash`` exists to detect, and the state this recap is
most needed in. ``components/clinical_summary.circle_sources`` is settled at filter time
for exactly this reason (issue #95).

Streamlit-free, like the rest of ``filters/``: it returns values and the results view
renders them.
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass, field
from typing import Any, Mapping

from config.param_labels import ParamLabel, labels_for

from .absent_columns import civic_clause_applies, filled_input_note
from .variant_filters import FREQUENCY_COLUMNS, Diagnostics, _Settings


def _has_frequency_column(available_columns) -> bool:
    """Whether the app's frequency layer had anything to judge a variant by.

    ``frequency_mask`` returns ``None`` where no member of :data:`FREQUENCY_COLUMNS` is
    present and the filter then skips, so the threshold was set and never applied. Asked
    through the same constant the mask loops over rather than a list written here.
    """
    return bool(set(FREQUENCY_COLUMNS) & set(available_columns or ()))

#: What an empty guideline selection does. Deliberately borrows the wording of the
#: parameter page's own caption — *"An empty selection means that source places no
#: restriction — it drops out of the comparison rather than widening it"* — rather than
#: sharing the string with it, because that caption is a sentence introducing a tab and
#: this is a value standing after a label. They are two renderings of one fact and the
#: phrase they share is the load-bearing half: **no restriction**, not "nothing".
#:
#: What stops them drifting is therefore not a shared constant but that neither may say
#: something else: an empty source is not a filter set to reject everything, which is what
#: rendering the empty list as ``""`` or ``"none"`` would read as, and
#: ``test_an_empty_source_says_it_places_no_restriction`` fails on either.
NO_RESTRICTION = "nothing selected — this source places no restriction"


@dataclass(frozen=True)
class RecapLine:
    """One setting of the run: what it is called, and what it was set to."""

    label: str
    value: str


@dataclass(frozen=True)
class RunRecap:
    """Everything the results view needs to say what produced the report beside it.

    Attributes:
        arm: the arm the run was made on. Named separately from the lines because it is
            the setting that decides which of the others exist at all.
        lines: one per setting live on this arm, in control order.
        note: what a reader needs to know about columns the file did not carry, or
            ``None``. This is the only account of a filled column that outlives the render
            it was produced in — the filter's own warning is drawn once and then gone.
        params: the parameters verbatim, carried for the downloadable run report. Not for
            display: the report is a file that outlives the session, and a dump cannot
            misdescribe what it prints, whereas prose in a downloaded text file has no
            control beside it to be checked against.
    """

    arm: str
    lines: tuple[RecapLine, ...] = ()
    note: str | None = None
    params: dict = field(default_factory=dict)


def _spell_depth(value) -> str:
    # `>=` in the vendored clause, so "or more" is exact. Zero is not "0 reads or more",
    # which reads as a setting when it is the absence of one.
    return f"{int(value):,} reads or more" if value else "no minimum"


def _spell_vaf(value) -> str:
    # `>` in the vendored clause, not `>=` — a variant sitting exactly on the threshold is
    # dropped, so "or above" would be wrong by one boundary.
    return f"above {value:.4g}"


def _spell_frequency(value) -> str:
    # The caller short-circuits at 1.0 and the mask is all-True there by algebra, so this
    # is not a rounding of "very permissive" — the filter genuinely does not run.
    if value >= 1.0:
        return "no limit — this filter is off"
    return f"{value:.4g} or below"


def _spell_genes(selection) -> str:
    # `restricts` rather than a truth test on the symbols: an all-invalid paste asked for a
    # restriction and got none, and the two states are not the same thing to say.
    if not selection.restricts:
        return "every gene — no gene restriction"
    count = len(selection.symbols)
    return f"{count:,} gene{'s' if count != 1 else ''}"


def _spell_terms(terms) -> str:
    return ", ".join(terms) if terms else NO_RESTRICTION


#: Appended to a setting the user made that the run then did not apply. The setting is
#: still printed, because the recap answers *what did I ask for* as well as *what happened*
#: and a user checking a selection against their memory needs to see it.
NOT_APPLIED = "{value} — not applied: {because}"

#: Two clauses can silently not run, both because the MAF lacks a column, and neither is
#: reported by :func:`~filters.absent_columns.filled_input_note` — nothing is filled in
#: either case, because the pipeline and the app respectively guard the clause instead.
#: They are the reason :func:`describe_run` needs the file's columns and not just its
#: parameters: without them the recap states a restriction that never applied.
NO_CIVIC = "this MAF carries no CIViC annotation"
CIVIC_SKIPPED = "CIViC annotation was skipped"
NO_FREQUENCY_COLUMNS = "this MAF carries no population-frequency columns"


def spell(row: ParamLabel, settings: _Settings, available_columns) -> str:
    """One setting's value, as a reader should meet it.

    The ClinVar rows are the reason a row carries ``terms``: one parameter is drawn by two
    controls, so the recap partitions the selection the same way the page does. A term the
    vocabulary has never heard of therefore appears under neither heading — which is the
    honest outcome, since neither control offers it and no widget could have produced it.

    ``available_columns`` is what lets a line say a setting **did not apply**. Two clauses
    can be dropped whole by a column the file does not carry, and in neither case is
    anything filled or warned about by ``absent_columns``: the pipeline guards CIViC for
    itself, and the app's frequency layer skips when no frequency column is present. On the
    reference MAFs the CIViC case is 100 files out of 100 — so a recap that could not say
    this would be stating a restriction that never ran, on almost every file the app has
    ever been given.
    """
    if row.id == "civic_keep" and not civic_clause_applies(
        settings.skip_civic, available_columns
    ):
        return NOT_APPLIED.format(
            value=_spell_terms(settings.civic_keep),
            because=CIVIC_SKIPPED if settings.skip_civic else NO_CIVIC,
        )

    if (
        row.kind == "frequency"
        and settings.max_freq_population < 1.0
        and not _has_frequency_column(available_columns)
    ):
        return NOT_APPLIED.format(
            value=_spell_frequency(settings.max_freq_population),
            because=NO_FREQUENCY_COLUMNS,
        )

    if row.kind == "depth":
        return _spell_depth(settings.min_depth)
    if row.kind == "vaf":
        return _spell_vaf(settings.vaf_threshold)
    if row.kind == "frequency":
        return _spell_frequency(settings.max_freq_population)
    if row.kind == "genes":
        return _spell_genes(settings.gene_selection)
    if row.kind == "retention":
        # The control is "Auto-retain pathogenic variants" and the parameter is
        # `skip_pathogenic`, so the polarity inverts here rather than in the label.
        return "off" if settings.skip_pathogenic else "on"
    if row.kind == "exclude":
        classifications = settings.exclude_classifications
        return ", ".join(classifications) if classifications else "nothing excluded"
    if row.kind == "keep":
        if row.terms:
            return _spell_terms(
                [term for term in settings.clinvar_keep if term in row.terms]
            )
        return _spell_terms(getattr(settings, row.id))
    raise ValueError(f"no spelling for {row.kind!r}")


def describe_run(
    params: Mapping[str, Any],
    available_columns,
    diagnostics: Diagnostics | None = None,
) -> RunRecap:
    """Describe the run these parameters just made over a MAF carrying these columns.

    Called by the filter run rather than by the page that draws it, so that what is carried
    is a statement about the cut on screen and not about the controls as they stand now.

    ``available_columns`` has no default, deliberately. A recap built from parameters alone
    cannot tell whether the CIViC or frequency clause ran, and the failure is silent and
    confident: it prints the user's selection as though it had been applied. Making the
    caller supply the file's columns is what keeps that state unreachable rather than
    merely unlikely — there is no "columns unknown" branch to fall into.

    ``diagnostics`` is optional only so the recap can be built for a run that produced
    none; every live call site has one, and without it the report loses its account of the
    columns the file did not carry.
    """
    settings = _Settings.from_params(params)
    lines = tuple(
        RecapLine(row.label, spell(row, settings, available_columns))
        for row in labels_for(settings.arm)
    )
    note = None
    if diagnostics is not None:
        note = filled_input_note(diagnostics.filled_columns, diagnostics.degraded_columns)
    # Deep, because the caller keeps mutating the live dict — every widget on the parameter
    # page writes into it in place, so a shallow copy would leave the lists shared and this
    # snapshot would follow the controls it exists to be independent of.
    return RunRecap(
        arm=settings.arm, lines=lines, note=note, params=deepcopy(dict(params))
    )
