"""The results grid: its columns, its widths, and the notes written onto it.

The AgGrid table plus everything that decides what it shows — which columns are
visible (through ``config.columns.resolve_visible_columns``), which are frozen at the
left edge, how wide each is — and the variant dialog a row opens.

It also owns the columns the *user* writes. ``Notes`` and any annotation column they
invent live in session state, are built onto a copy of the frame at display time by
:func:`_add_derived_columns`, and never touch ``st.session_state.filtered_data``.
:func:`with_user_columns` is the seam for a download outside this module, and
:func:`shown_columns` the set the grid *opens* with — one definition each, because a
download that quietly drops a note is the single failure this app cannot afford
(issue #67).

Since issue #92 the download beneath the grid delivers what the grid is *currently*
showing, not what it opened with, so :func:`create_data_table` returns the list it
resolved and the download reads it. That closes the last gap of the shape issue #60
opened: a column the user adds through *Add columns* is now in the file offered
directly beneath it, as a column they invent in the variant dialog already was.

**Why the column tables in here did not follow ``PRIORITY_COLUMNS`` into ``config/``.**
Issue #76 asked whether any pure column decision found here belongs in a module a bare
``pytest`` can read, and one did. The rest — ``_COLUMN_WIDTH_MAP``, ``_PINNED_LEFT``,
``_NUMERIC_COLUMNS``, ``_leading_columns`` and ``_display_order`` — are *viewport*
decisions: pixel widths, which columns freeze at the left edge, which get a numeric
filter box. They are pure too, and that is not the test. They configure one widget, they
are meaningless to anything that is not drawing it, and ``_PINNED_LEFT``'s own comment
gives the reason they must stay downstream of ``resolve_visible_columns`` rather than
beside it: the resolver's output has to open with the pipeline's list as an exact
prefix, so a module that could reorder it is a module that could break the export.
"""

# For `_variant_key`'s `str | None`, which a Python below 3.10 would otherwise evaluate at
# definition time and reject. The same import for the same reason as
# `filters/absent_columns.py` and `filters/numeric_columns.py`.
#
# The floor this once cited — `requirements.txt`'s 3.9 — is gone: the pins now require 3.11
# (issue #256), which supports the syntax outright. The import stays because it costs nothing
# and these three modules should not each depend on that floor never dropping again.
from __future__ import annotations

import warnings

import streamlit as st
import pandas as pd

# Try to import AgGrid - only show warning when actually using the component
try:
    from st_aggrid import AgGrid, GridOptionsBuilder, JsCode

    AGGRID_AVAILABLE = True
except ImportError:
    # Fallback for when st_aggrid is not available
    AGGRID_AVAILABLE = False

from config.columns import MissingColumnsWarning, resolve_visible_columns
from config.constants import TABLE_HEIGHT
from config.missing_values import is_blank, says_nothing_over

from .clinical_summary import circle_axis_notes, pathogenicity_circle_legend
from .variant_detail import render_variant_detail_panel


# A note is **not** cleared when a different file is opened, and there is deliberately no
# function here that does it. Issue #64 added one — notes are keyed by variant identity and
# not by file, so it read the cross-file surfacing as one sample's note leaking onto
# another's report — and issue #67 then settled the question #64's own docstring referred to
# it: a note is a statement *about the variant*, meant to be read by everyone at the
# institute once there is a server to hold it. The surfacing is the design.
#
# Which makes the protection a matter of what a note may contain rather than how long it
# lives: the variant dialog says a note describes the variant and not the patient, and asks
# for no patient identifiers. Restoring a clearing rule means reopening #67, not fixing an
# oversight — see `_show_variant_dialog` below for the copy that carries this.


#: The five columns a note is stored against, in the order they are spelled into the key.
#:
#: Four of them are *core*, so a MAF that does not have the **column** is refused before any of
#: this runs — but a column is not a value, and this module's question is whether the row says
#: anything. ``Tumor_Seq_Allele2`` is the fifth and nothing requires it at all, which is why it
#: is the case issue #127 measured; the other four are guarded on the same footing because
#: :func:`_variant_key` cannot tell the two failures apart and must not name the wrong one.
#:
#: Measured across the 165 real MAFs on this machine, 326,204 rows: **no component is blank on
#: any row of any file**, this column included. So every hole here is a hole nobody has fallen
#: into — the reason to handle it is that the app's own reference table says one of the five is
#: a column your MAF need not carry.
_KEY_COLUMNS = (
    "Hugo_Symbol",
    "Chromosome",
    "Start_Position",
    "Reference_Allele",
    "Tumor_Seq_Allele2",
)


def _missing_key_components(row) -> tuple[str, ...]:
    """Which of :data:`_KEY_COLUMNS` this row does not say. Empty when it identifies a variant.

    Separate from :func:`_variant_key` so the dialog can *name* the hole rather than guess at
    it. The first version of this told every keyless row to check ``Tumor_Seq_Allele2``, which
    is the likely case and not the only one — a blank ``Hugo_Symbol`` would have sent the user
    to look at a column that was fine.

    Asks ``config.missing_values.is_blank`` and **not** its sibling ``says_nothing``, which is
    the app's one answer for everything a user reads. The two differ on exactly the annotator's
    sentinels, and an identity is the one place where that difference has to be honoured: a
    variant whose ``Hugo_Symbol`` is ANNOVAR's ``Unknown`` — 2,605 rows across 136 of the 303
    MAFs measured for issue #131 — is still identified by its position and alleles, so it can
    carry a note, while a variant with an *empty* key column cannot, because two rows differing
    only in that component would collapse onto one key (issue #127). Reaching for the display
    predicate here would have quietly withdrawn notes from nearly half of real files.

    ``str(row.get(name, ""))`` spells an absent column and a missing value identically, which is
    why the check is on the value rather than on the frame's columns: both cases are a hole in
    the identity, and only one of them is visible from the header row.
    """
    return tuple(
        name
        for name, value in zip(_KEY_COLUMNS, _rendered_key_components(row))
        if is_blank(value)
    )


def _rendered_key_components(row) -> list[str]:
    """The five values as the key spells them, so one rendering serves both readers above."""
    return [str(row.get(name, "")).strip() for name in _KEY_COLUMNS]


def _spelled(columns) -> str:
    """Column names for the dialog, backticked and never truncated.

    At most five, so an "and N more" would be shortening a list that fits — and a user who has
    to check their header row needs all of them at once. The same reasoning, and the same
    shape, as ``filters/absent_columns._names``; not imported from there because that module
    is the filter's and this string is the grid's.
    """
    return ", ".join(f"`{column}`" for column in columns)


def _variant_key(row) -> str | None:
    """The identity a note is stored against, or ``None`` if this row does not carry one.

    ``gene:chrom:pos:ref>alt``, which is issue #67's decision about what a note is *about*: the
    variant, not the file and not the case, so it surfaces on any MAF carrying that variant.

    Returns ``None`` rather than a key with a hole in it. A missing component used to be
    interpolated as the empty string, so every variant that differed only in that component
    shared one key — on a MAF with no ``Tumor_Seq_Allele2``, every alternate allele at one
    position collapsed onto ``GENE:chr:pos:REF>`` and one note was shown against all of them.
    That is the failure #67's design makes worst: a note is written for everyone at the institute
    to read, so attaching it to the wrong variant is worse than not being able to attach it at
    all. It also silently reused one Streamlit widget key across those rows.

    Nothing is lost by tightening it: a note lives in session state and nothing writes one to
    disk (issue #67), so there is no stored key under the old spelling to migrate.

    Callers must handle ``None``: :func:`_add_derived_columns` leaves ``Notes`` empty for such a
    row and :func:`_show_variant_dialog` names the missing component instead of drawing fields
    whose contents it cannot store honestly.
    """
    if _missing_key_components(row):
        return None
    gene, chrom, pos, ref, alt = _rendered_key_components(row)
    return f"{gene}:{chrom}:{pos}:{ref}>{alt}"


def _components_of_key(key: str) -> list[str]:
    """Take :func:`_variant_key`'s spelling apart again, in :data:`_KEY_COLUMNS` order.

    The stores hold the joined string and nothing else, so :func:`written_notes` — which
    deliberately consults no report frame, that being what lets it rescue a note written on a
    variant the filters rejected — has to reverse the join rather than look a row up.

    Both splits are **forward**, and the direction is measured rather than reasoned. Across
    340,972 rows of the 198 real MAFs on this machine carrying all five columns,
    ``Hugo_Symbol``, ``Chromosome``, ``Start_Position`` and ``Reference_Allele`` never contain
    ``:`` or ``>`` — but ``Tumor_Seq_Allele2`` contains ``>`` on two of them, where the allele
    is the symbolic ``<DEL>`` or ``<DUP>``. So ``rsplit(">", 1)``, the obvious reading of "the
    separator is the last one", takes ``C><DEL>`` apart into ``C><DEL`` and ``""``, losing the
    reference allele and the alternate together. The separator that can repeat is on the side
    that never carries it, so splitting from the front is what is correct here.

    Total by construction. ``_variant_key`` writes exactly three colons and refuses a key with
    a blank component, so the padding is unreachable through it — but this runs inside a crash
    handler, where raising a second exception would replace the app's account of the first one
    with Streamlit's.
    """
    parts = key.split(":", 3)
    parts += [""] * (4 - len(parts))
    gene, chrom, pos, alleles = parts
    ref, _, alt = alleles.partition(">")
    return [gene, chrom, pos, ref, alt]


def _add_derived_columns(data: pd.DataFrame) -> pd.DataFrame:
    """Add derived columns: Sample_Name and Notes."""
    # Sample_Name: prefer Tumor_Sample_Barcode, fall back to Matched_Norm_Sample_Barcode.
    #
    # This one asks the *display* predicate, unlike the note key above: the question is whether
    # the barcode is worth showing as a name, not whether it identifies anything. Its old set
    # was the only one in the app that knew ``__UNKNOWN__``, and that member is load-bearing —
    # Funcotator writes it for a value it could not determine, and ``Tumor_Sample_Barcode`` is
    # entirely ``__UNKNOWN__`` in 70 of the 303 MAFs measured for issue #131.
    if "Tumor_Sample_Barcode" in data.columns:
        tsb = data["Tumor_Sample_Barcode"].astype(str)
        nsb = (
            data["Matched_Norm_Sample_Barcode"].astype(str)
            if "Matched_Norm_Sample_Barcode" in data.columns
            else pd.Series("", index=data.index)
        )
        data["Sample_Name"] = tsb.where(~says_nothing_over(tsb), nsb)
    # Notes from session state
    if "variant_notes" not in st.session_state:
        st.session_state.variant_notes = {}
    notes_dict = st.session_state.variant_notes
    data["Notes"] = data.apply(lambda row: _stored(notes_dict, row), axis=1)
    # Custom annotation columns from session state
    if "custom_annotations" not in st.session_state:
        st.session_state.custom_annotations = {}
    for col_name, values_dict in st.session_state.custom_annotations.items():
        data[col_name] = data.apply(
            lambda row, vd=values_dict: _stored(vd, row), axis=1
        )
    return data


def _stored(store: dict, row) -> str:
    """What the user wrote against this row, or the empty string if it has no identity.

    Both stores are keyed by :func:`_variant_key`, and both read through here rather than
    calling ``store.get(_variant_key(row), "")`` — which would have worked, since ``None`` is
    never written as a key, but only by accident. Saying it once means the ``None`` case is
    handled where it is decided instead of relying on two call sites to keep getting the
    default right.
    """
    key = _variant_key(row)
    return store.get(key, "") if key is not None else ""


def shown_columns(data: pd.DataFrame) -> list:
    """The columns the user is looking at: the resolver's list, then their own annotations.

    One list behind the grid and the "shown columns" downloads. They used to build it
    separately — the grid called the resolver and appended ``_custom_annotation_columns()``
    (issue #60), the download called the resolver alone — so a column the user invented was
    on screen and missing from the CSV offered directly beneath it. ``Notes`` escaped that
    only because it is in ``APP_EXTRA_COLUMNS`` and the resolver returns it; the columns the
    user names themselves cannot be there and must not be, since the resolver is
    streamlit-free and those names live in session state. So the append is the app's job,
    and with a download reading it too it belongs in exactly one place.

    Public, unlike its neighbours here, because that download is in another module now:
    ``results_view`` offers a *Download shown columns* button beside each tab's grid and
    has to ask the grid what "shown" means rather than answer it a second time. The
    underscore it used to carry described where it sat, not who reads it.
    """
    columns = _visible_columns(data)
    # Appended, not inserted: the resolver's order is the pipeline's and must stay that way
    # for the download. Pinning is what puts these at the left edge.
    return columns + [
        col
        for col in _custom_annotation_columns()
        if col in data.columns and col not in columns
    ]


def with_user_columns(data: pd.DataFrame) -> pd.DataFrame:
    """``data`` with the app's derived columns on it, including the ones the user wrote.

    The seam for a download button outside this module. ``Notes`` and the user's own
    annotation columns are built onto a *copy* at display time, so they are not on
    ``st.session_state.filtered_data``: a download that serialises that frame directly
    carries no ``Notes`` column at all, which is what the two buttons under *Download
    Results* did. Nothing stores a note (issue #67), so a download is the only way one
    leaves the session, and a download that quietly drops it is the single failure this
    app cannot afford.
    """
    return _add_derived_columns(data.fillna(""))


def written_notes() -> pd.DataFrame:
    """Everything the user has written, against the five columns their own MAF names.

    The second download seam out of this module, and it is the mirror image of
    :func:`with_user_columns`: that one puts the writing onto a report, this one takes the
    writing *without* a report. Issue #144 needed it for ``MAFigate.py``'s outermost handler,
    where the app's advice of last resort is to refresh — which starts the session over, and a
    note lives in session state and nowhere else (issue #67), so the advice as it stood was
    *discard your writing*, offered by a surface that could not tell whether there was any.

    Consulting no frame is the point, three times over. ``st.session_state.filtered_data`` holds
    only the variants that **passed**, and both tabs draw the note dialog, so serialising it
    would silently drop every note written on a rejected variant — the one failure
    :func:`with_user_columns` exists to prevent, reintroduced by the rescue meant to prevent it.
    It is also ``None`` whenever a file has been opened and not yet filtered, while the notes
    survive that (issue #67 reversed #64's clearing rule), so a report-shaped rescue would have
    had nothing to offer exactly when there was something to lose. And a full derived-column
    build over both frames is real work inside a handler that only runs because something has
    already gone wrong.

    Empty — with its columns, so a caller can read them off it — when nothing has been written.
    That is the whole gate: there is no state where the user has written something and this has
    nothing to serialise, so the handler asks this one question rather than two.

    A key with no value against it anywhere is dropped rather than written as a blank row. The
    save path already pops an emptied note instead of storing ``""``, so this is defensive, but
    a rescue file listing variants the user wrote nothing about would misreport what it holds.
    """
    notes = st.session_state.get("variant_notes", {}) or {}
    annotations = st.session_state.get("custom_annotations", {}) or {}
    # Sorted so the column order is the user's own naming, stably, rather than dict insertion
    # order — a rescue file is read next to the MAF it came from, possibly much later.
    annotation_columns = sorted(annotations)

    keys = {key for key, value in notes.items() if str(value).strip()}
    for values in annotations.values():
        keys |= {key for key, value in values.items() if str(value).strip()}

    rows = []
    for key in sorted(keys):
        row = dict(zip(_KEY_COLUMNS, _components_of_key(key)))
        row["Notes"] = notes.get(key, "")
        for name in annotation_columns:
            row[name] = annotations[name].get(key, "")
        rows.append(row)

    return pd.DataFrame(rows, columns=[*_KEY_COLUMNS, "Notes", *annotation_columns])


def _view_details_label(row: "pd.Series") -> str:
    """The label on the button that opens the variant dialog, wherever it is drawn.

    There are two such buttons since issue #160 — the grid's, below the selection, and the
    grid-less path's, below the selector — and one spelling of the label, because #149's
    ``chrchr17`` was a `chr` literal written into two places that had drifted apart. No
    literal here either: the column arrives prefixed (``config/chromosome_spelling.py``).
    """
    return (
        f"🔍 View details: **{row.get('Hugo_Symbol', '?')}** "
        f"{row.get('Chromosome', '?')}:{row.get('Start_Position', '?')}"
    )


def create_data_table(
    data: pd.DataFrame, height: int = TABLE_HEIGHT, key_suffix: str = None
) -> list:
    """Draw the interactive table, and return the columns it decided to show.

    The return value is the seam issue #92 needed: the download beneath this grid
    delivers what the grid is showing, so it has to ask rather than answer a second
    time. The alternative was for ``results_view`` to read ``show_all_{key_suffix}``
    and ``extra_cols_{key_suffix}`` back out of session state and re-derive the list —
    which would put this function's widget-key spelling in another module and give the
    app two answers to *what is shown*, the one thing :func:`shown_columns` exists to
    prevent.

    Both exits return it, including the AgGrid-less fallback, because that path draws a
    table too and a user on it downloads from the same button. A caller that never gets
    here — ``create_data_table`` raised and the error fallback rendered instead — has no
    grid to follow and falls back to :func:`shown_columns`.

    The list is in the order the grid resolved: the pipeline's, then the user's
    annotations, then anything added through *Add columns*. It is **not** display order.
    ``_display_order`` and ``_PINNED_LEFT`` move ``Notes`` and the clinical columns to
    the left edge of the screen only, and issue #60 put that step downstream of the
    resolver precisely so the export keeps the pipeline's order.
    """

    # Clean data for display (fillna returns a new DataFrame, no need to copy)
    cleaned_data = data.fillna("")

    # Convert problematic columns to strings to avoid display issues
    for col in cleaned_data.columns:
        if cleaned_data[col].dtype == "object":
            cleaned_data[col] = cleaned_data[col].astype(str)

    # Add derived columns (Sample_Name, Notes)
    cleaned_data = _add_derived_columns(cleaned_data)

    # --- Column selector (outside fragment so it's always visible) ---
    all_columns = list(cleaned_data.columns)
    show_all = st.checkbox(
        "Show all columns", value=False, key=f"show_all_{key_suffix}"
    )
    # Resolved on both branches, and unused on one of them. Ticking *Show all columns*
    # used to skip this call, which skipped the resolver's report with it — so the app
    # went silent about a MAF missing three pipeline columns exactly when the user had
    # asked to see everything (issue #93, measured at zero warnings). The absence is a
    # fact about the file: showing every column the file *has* cannot show one it has
    # not, so widening the view is the last state in which the sentence should vanish.
    default_cols = shown_columns(cleaned_data)
    if show_all:
        visible_columns = all_columns
    else:
        extra_pool = [c for c in all_columns if c not in default_cols]
        extra_cols = st.multiselect(
            "Add columns:",
            options=extra_pool,
            default=[],
            key=f"extra_cols_{key_suffix}",
            placeholder="Search and add columns...",
        )
        visible_columns = default_cols + extra_cols

    display_data = cleaned_data[visible_columns]

    if not AGGRID_AVAILABLE:
        st.info(
            "📊 Using standard Streamlit dataframe display (streamlit-aggrid not available)"
        )
        st.dataframe(
            display_data,
            height=height,
            use_container_width=True,
            # No pinning here, so the lead columns are brought forward instead. The
            # frame itself is untouched — column_order is a display argument, so the
            # export still reads display_data in the resolver's order.
            column_order=_display_order(visible_columns),
            key=(
                f"fallback_table_{key_suffix}_{hash(str(visible_columns))}_{len(display_data)}"
                if key_suffix
                else None
            ),
        )
        st.caption(f"📊 Showing {len(display_data):,} rows")

        # Row selection fallback for variant detail panel.
        #
        # The selector chooses; the button opens. Both halves are issue #160, and neither is
        # cosmetic: the selectbox was given no `index`, so it answered `0` rather than `None`
        # and the guard beneath it — `if row_idx is not None` — was true on every single
        # render. Measured before the fix, on three untouched renders of a three-variant
        # frame: three openings of the dialog on the first variant, nobody having touched
        # anything. This path is not inside a fragment either, so *any* interaction on the
        # page reran it and the dialog came back. That is the standing-condition failure
        # issue #159 designed the grid's double-click around, sitting here all along.
        #
        # `index=None` alone would only postpone it — once a variant is chosen the widget
        # keeps answering with it, and the dialog would re-open on every later rerun. So the
        # opening is moved behind a press, which is a thing that happens once. The dialog is
        # unreachable by double-click here because `st.dataframe` has no such event, and this
        # is what the two paths are allowed to differ on: the gesture, not whether a dialog
        # appears unbidden.
        #
        # It is not allowed to differ on *which* variant opens — and this route was never at
        # risk of that, which issue #312 checked here rather than assuming. The grid's defect
        # (#309/#310) needed a row identifier making a round trip through browser JSON; this
        # route has none. `range(len(cleaned_data))` makes the selectbox's value a position in
        # the frame in hand, `format_func` reads the label the reader sees off that same
        # position, and the row opened is `cleaned_data.iloc[row_idx]` — one frame, one
        # positional read, no `_full_row` and no three-column match to fall through to. The
        # variant named in the option is therefore the variant the dialog gets, by
        # construction.
        if len(cleaned_data) > 0:
            row_idx = st.selectbox(
                "Select a variant to view details:",
                range(len(cleaned_data)),
                index=None,
                placeholder="Choose a variant...",
                format_func=lambda i: (
                    f"{cleaned_data.iloc[i].get('Hugo_Symbol', '?')} "
                    f"{cleaned_data.iloc[i].get('Chromosome', '?')}:"
                    f"{cleaned_data.iloc[i].get('Start_Position', '?')}"
                ),
                key=f"row_select_{key_suffix}",
            )
            # No selection, no button — the same shape as the grid path, where the button is
            # drawn only under a selected row. Left on `index=0` it would sit there
            # permanently offering the first variant in the file, which nobody chose.
            if row_idx is not None:
                selected = cleaned_data.iloc[row_idx]
                if st.button(
                    _view_details_label(selected),
                    key=f"view_detail_{key_suffix}",
                ):
                    _show_variant_dialog(selected)
        return visible_columns

    # Wrap AgGrid + detail panel in a fragment so row selection
    # only re-runs this section, not the entire page.
    @st.fragment
    def _render_table_fragment():
        _render_aggrid_with_detail(
            display_data, cleaned_data, height, key_suffix, visible_columns
        )

    _render_table_fragment()
    return visible_columns


@st.dialog("Variant Details", width="large")
def _show_variant_dialog(row: "pd.Series") -> None:
    """Show variant detail panel in a modal dialog (closable with ESC)."""
    render_variant_detail_panel(row)

    # --- Notes & custom annotations ---
    st.markdown("---")
    st.markdown("**📝 Notes & Annotations**")
    # What a note is, and what becomes of it. Both halves are load-bearing (issue #67):
    # a note is keyed by the variant alone, so it appears on any file carrying that
    # variant and is a statement about the variant rather than about the case in front
    # of you — and nothing writes one to disk, so the download is the only way it
    # survives the tab being closed.
    st.caption(
        "A note describes the variant, not the patient — it appears on any file "
        "carrying this variant, so do not include patient identifiers. Notes live "
        "for this session only: download the table to keep them."
    )
    vkey = _variant_key(row)

    # No identity, no note. The caption above has just promised that a note appears on any file
    # carrying this variant, and a key with a hole in it cannot keep that promise: it would
    # attach one note to every alternate allele at this position, on this file and on the next
    # one. So the fields are not drawn, and the reason is given in the terms the caption used —
    # `_variant_key`'s docstring has why refusing beats guessing.
    #
    # The message names the components that are actually missing rather than the likely one. It
    # said "check that your MAF carries `Tumor_Seq_Allele2`" while the key refuses on any of
    # five, so a blank `Hugo_Symbol` would have sent the reader to a column that was fine.
    #
    # This returns before the annotation widgets as well as the note field, and says so: an
    # invented column is keyed exactly as a note is, so there is no value this row could hold in
    # one either. Offering *Add new annotation column* here would let the user create a column
    # they cannot then fill for the variant in front of them.
    missing = _missing_key_components(row)
    if missing:
        st.info(
            "**Notes are not available for this variant.** A note — and any annotation column "
            "you add — is stored against the variant itself: gene, position, and both alleles. "
            f"This row does not say what {_spelled(missing)} "
            f"{'is' if len(missing) == 1 else 'are'}"
            ", so anything written here could not be told apart from writing about a different "
            "variant at the same position. Every other variant in this file is unaffected."
        )
        return

    # Initialise session state stores
    if "variant_notes" not in st.session_state:
        st.session_state.variant_notes = {}
    if "custom_annotations" not in st.session_state:
        st.session_state.custom_annotations = {}

    # Notes
    current_note = st.session_state.variant_notes.get(vkey, "")
    new_note = st.text_area(
        "Notes",
        value=current_note,
        key=f"note_input_{vkey}",
        placeholder="Add your notes here...",
    )

    # Custom annotation columns — show existing ones
    ann_values = {}
    for col_name in sorted(st.session_state.custom_annotations.keys()):
        current_val = st.session_state.custom_annotations[col_name].get(vkey, "")
        ann_values[col_name] = st.text_input(
            f"🏷️ {col_name}",
            value=current_val,
            key=f"ann_{col_name}_{vkey}",
        )

    # Create new annotation column
    with st.expander("Add new annotation column"):
        new_col_name = st.text_input(
            "Column name:",
            key=f"new_ann_col_{vkey}",
            placeholder="e.g. Actionable, Review_Status, Panel...",
        )
        if st.button("Create column", key=f"create_ann_{vkey}"):
            if new_col_name and new_col_name.strip():
                col_name = new_col_name.strip()
                if col_name not in st.session_state.custom_annotations:
                    st.session_state.custom_annotations[col_name] = {}
                # Fragment scope, because this rerun exists to redraw *this dialog* with a
                # field for the column just created. `st.dialog` wraps its content in a
                # fragment and an app-scope rerun is Streamlit's documented way to close a
                # dialog, so the plain `st.rerun()` that stood here closed the dialog on
                # the user instead — they named a column, pressed Create, and had to find
                # the row again to type into it. Found while resolving issue #140, which
                # did not ask for it.
                st.rerun(scope="fragment")

    # Save button
    if st.button("💾 Save", key=f"save_note_{vkey}", type="primary"):
        # Save notes
        if new_note.strip():
            st.session_state.variant_notes[vkey] = new_note
        else:
            st.session_state.variant_notes.pop(vkey, None)
        # Save custom annotations
        for col_name, new_val in ann_values.items():
            if new_val.strip():
                st.session_state.custom_annotations[col_name][vkey] = new_val.strip()
            else:
                st.session_state.custom_annotations[col_name].pop(vkey, None)
        # "Saved!" was the app's strongest claim of durability and its least true one:
        # this writes to session state and nothing else. The scope goes in the notice
        # rather than only in the caption above, because this is the moment the user
        # believes their writing is safe.
        #
        # Parked, not drawn: `st.rerun` discards this frame, so the `st.success` that stood
        # here reached nobody (issue #140). The rerun stays and is app-scoped on purpose —
        # it is what closes the dialog *and* rebuilds the grid so `Notes` shows what was
        # just written. The alternative, dropping the rerun to say it here, was measured
        # and rejected: dismissing a dialog is client-side only and sends the server
        # nothing, so the table behind would have gone on showing no note while the app
        # claimed one was saved.
        park_note_confirmation("Saved for this session.")
        st.rerun()


#: Where the note dialog leaves its confirmation, for the run after the rerun to say.
#:
#: A survive-the-run stash, unlike ``_MISSING_COLUMN_REPORTS`` below: the dialog writes it
#: and then reruns, and it is the *next* render that draws it. That is the whole point —
#: ``st.rerun`` raises immediately and throws away everything the current run drew, so a
#: sentence written beside the save reaches nobody (issue #140, following issue #133 on the
#: parameter page). It is not that page's ``park_confirmation``: that one parks into the
#: parameter page's own notice slot and can carry an arm clause, and neither is true here.
_NOTE_CONFIRMATION = "note_confirmation"


def park_note_confirmation(confirmation: str) -> None:
    """Park one sentence for the render after the rerun, because this one is discarded.

    Args:
        confirmation: one sentence saying what happened, ending in a full stop.
    """
    st.session_state[_NOTE_CONFIRMATION] = confirmation


def draw_parked_note_confirmation() -> None:
    """Say what the dialog parked, once, on the run that the user actually sees.

    A ``st.toast`` rather than a banner, and that is a claim about *where* rather than a
    dislike of banners: by the time this runs the dialog has gone and the grid has been
    rebuilt with the note in it, so the table is already the evidence that the save
    happened. What it cannot show is how long the note lasts, which is the one thing this
    sentence is for (issue #67) — so it is said over the table, briefly, rather than given
    a permanent place in the page beside the filter's own warnings.

    Drains rather than reads: a toast redrawn on every later render would follow the user
    around the app.
    """
    confirmation = st.session_state.pop(_NOTE_CONFIRMATION, None)
    if confirmation:
        st.toast(confirmation, icon="📝")


#: Where a resolver call leaves what it found missing, for the page to say once.
#:
#: Not a survive-the-run stash like ``filter_notes`` next door, which is written by a
#: filter run and read by whichever later render draws a page. This one is filled and
#: drained inside a **single** render: :func:`clear_missing_column_reports` empties it
#: before the section draws, every resolver call on that render appends to it, and
#: :func:`drain_missing_column_reports` hands the collected text to the page-level slot at
#: the end. Carried through session state rather than returned up the call chain because
#: the collectors are three frames deep in a grid and the drain is in a page module, and
#: this is the seam the page already uses for the filter's own banners.
_MISSING_COLUMN_REPORTS = "missing_column_reports"


def clear_missing_column_reports() -> None:
    """Empty the per-render collection. Called by the page before its section draws."""
    st.session_state[_MISSING_COLUMN_REPORTS] = []


def drain_missing_column_reports() -> list:
    """What this render's resolver calls found missing — each distinct sentence once.

    Deduplicated by exact text, in first-said order. The two grids resolve over frames
    with the same columns, so on any real MAF they produce one identical sentence and the
    user reads it once (issue #93). Dedup is by *text* rather than by a flag saying
    "already said", so two calls that genuinely disagreed would both be reported: the
    collapse is of repetition, not of information.

    Drains rather than reads, so a render that draws no grid — the user is on the Load
    section, or the MAF was refused — cannot show the last render's sentence.
    """
    reports = st.session_state.get(_MISSING_COLUMN_REPORTS) or []
    st.session_state[_MISSING_COLUMN_REPORTS] = []
    seen = set()
    unique = []
    for text in reports:
        if text not in seen:
            seen.add(text)
            unique.append(text)
    return unique


def _visible_columns(data: pd.DataFrame) -> list:
    """The columns to show for ``data``, from the app's single resolver.

    Every display path goes through here — both grids, and a download that has no grid to
    follow — so what the user sees and what they download cannot disagree. This replaced a
    hand-maintained ``_DEFAULT_VISIBLE_COLUMNS`` list; ``config/columns.py`` records
    why. The order is the pipeline's, with the app's extras after it.

    The arm comes from the parameters the user filtered with, so the germline arm gets
    ``InterVar``/``RENOVO_*`` and the somatic arm ``CancerVar`` and ``cosmic``, rather
    than the old list's union of both, which showed empty columns on either arm.

    ``skip_civic`` is read for completeness and is ``False`` on this path for every set the
    app itself writes: it is a contract key, so the completion both parameter-carrying pages
    now do fills it from the contract, where it is ``False`` (issues #280, #289 — the older
    claim here, that ``parameter_config.py`` strips it as a deprecated name and nothing else
    writes it, stopped being true when it became a completed key). That changes no output —
    ``compute_keep`` drops
    the CIViC columns by their absence from the frame anyway — so this is passing the
    flag the resolver's contract asks for, not relying on it.

    Any column the MAF does not carry is dropped rather than raised on. **Reporting is
    still not optional per call site** — every call appends what it found — but since
    issue #93 it appends to :data:`_MISSING_COLUMN_REPORTS` instead of drawing an
    ``st.warning`` where it stands, and the page says the collected sentence once, above
    the section switch.

    The guarantee that clause was defending is unchanged and now reaches further. It was
    written when the two exports called this resolver themselves, so "a variant that
    stayed quiet would let an export silently lose columns" named a live call site; issue
    #92 made the exports inherit the grid's list, and both remaining callers are grids. A
    per-call-site ``st.warning`` was therefore no longer buying the export anything, while
    costing one sentence per grid — four copies before #92, two after, and **none at all**
    with *Show all columns* ticked, because that branch skipped the resolver entirely. The
    stash keeps the report unconditional and moves only where it is drawn.

    The one thing it needs in return: a caller that renders a grid outside
    ``show_data_loading_page`` would collect into a stash nothing drains, which is the
    silence the old clause forbade. There is no such caller, and
    ``tests/test_missing_column_warning.py`` holds the page to draining it.
    """
    # These two stay *defaulting* reads, deliberately, where the page's own equivalents now
    # index (issue #289). What they decide is which columns are shown, not what is kept, and
    # the page completes the parameters against their arm's contract before any grid draws —
    # so on every live path the fallbacks are unreachable. What is not unreachable is a
    # caller with no parameters at all: `create_data_table` is reached directly by the grid
    # and download seams' own tests, and making a column list raise for want of a
    # `sample_type` would ask each of those to state an arm it is not about.
    params = st.session_state.get("filter_params") or {}
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", MissingColumnsWarning)
        columns = resolve_visible_columns(
            sample_type=params.get("sample_type", "somatic"),
            skip_civic=bool(params.get("skip_civic", False)),
            available_columns=list(data.columns),
        )
    for entry in caught:
        if issubclass(entry.category, MissingColumnsWarning):
            st.session_state.setdefault(_MISSING_COLUMN_REPORTS, []).append(
                f"⚠️ {entry.message}"
            )
    return columns


# Column width overrides for known columns
_COLUMN_WIDTH_MAP = {
    "Pathogenicity_Overview": 180,
    "Clinical_Summary": 250,
    # Narrowed from 200 when Notes joined _PINNED_LEFT: three frozen columns at the old
    # widths take 630px, half a 13" laptop's table. 160 still previews a note's opening
    # words, the full text is a click away in the variant dialog, and the column is
    # resizable for anyone who wants more.
    "Notes": 160,
    "Sample_Name": 150,
    "project_id": 120,
    "Hugo_Symbol": 130,
    "Tumor_Sample_Barcode": 150,
    "Variant_Classification": 170,
    "Variant_Type": 100,
    "Chromosome": 110,
    "Start_Position": 120,
    "Protein_Change": 140,
    "ClinVar_VCF_CLNSIG": 160,
    "InterVar": 160,
    "CancerVar": 160,
    "RENOVO_Class": 130,
    "ESCAT": 80,
    "ESCAT_TISSUE": 120,
    "am_class": 120,
    "tumor_f": 90,
    "DP": 80,
    "t_alt_count": 100,
    "t_ref_count": 100,
    # The four spellings compute_keep can actually append. These replaced
    # gnomAD_exome_AF_grpmax and gnomad41_genome_AF, which sized columns that no MAF
    # carries and the resolver cannot produce.
    "gnomAD_exome_AF_raw": 160,
    "gnomAD_exome_AF": 160,
    "gnomAD_genome_AF_raw": 150,
    "gnomAD_genome_AF": 150,
    "dbSNP_ID": 130,
}

# Columns frozen to the left edge of the grid, in this order.
#
# Pinning is a *viewport* decision and lives only here, downstream of
# ``config.columns.resolve_visible_columns``. That matters: the resolver's output must
# open with the pipeline's column list as an exact prefix, so no app column can be moved
# up it. These three are already far to the right of the frame — all of them sit in
# ``APP_EXTRA_COLUMNS``, after every pipeline column — and AgGrid still renders them
# leftmost, in frame order, because they are pinned. So the table can lead with what the
# user needs while the column list, and therefore the export, stays the pipeline's.
#
# ``Notes`` earns its place by being the one column the user wrote themselves: unpinned
# it sat past roughly forty pipeline columns on the somatic arm, further right than every
# column they did not write. See issue #60.
_PINNED_LEFT = ("Clinical_Summary", "Pathogenicity_Overview", "Notes")

# Numeric filter columns
_NUMERIC_COLUMNS = {
    "Start_Position", "End_Position", "Chromosome",
    "tumor_f", "DP", "t_alt_count", "t_ref_count",
    "gnomAD_exome_AF_raw", "gnomAD_exome_AF",
    "gnomAD_genome_AF_raw", "gnomAD_genome_AF",
}


def _custom_annotation_columns() -> list:
    """The annotation columns this user created, in name order.

    ``_show_variant_dialog`` lets the user invent a named column — "Actionable",
    "Review_Status" — and ``_add_derived_columns`` writes it onto the frame. It is as
    user-authored as ``Notes``, so it gets the same treatment: shown by default and
    pinned left. Neither is optional. These names are not in ``APP_EXTRA_COLUMNS`` and
    cannot be — the resolver is streamlit-free and these live in session state — so
    without this the user creates a column, types into it, and it appears nowhere until
    they hunt it out of the "Add columns" list.
    """
    return sorted(st.session_state.get("custom_annotations", {}))


def _leading_columns(columns) -> list:
    """The entries of ``columns`` that belong at the left edge, in display order.

    ``_PINNED_LEFT`` first, then the user's own annotation columns after ``Notes``.
    Both the grid (via AgGrid pinning) and the AgGrid-less fallback (via
    ``st.dataframe(column_order=...)``) lead with exactly this, so the two paths agree
    on where a user's own writing sits.
    """
    present = set(columns)
    return [
        col
        for col in list(_PINNED_LEFT) + _custom_annotation_columns()
        if col in present
    ]


def _display_order(columns: list) -> list:
    """``columns`` with :func:`_leading_columns` moved to the front.

    For the fallback table only. AgGrid needs no such reordering — pinning already
    lifts those columns out of frame order — and the frame handed to AgGrid is left in
    resolver order so the export it feeds stays the pipeline's.
    """
    lead = _leading_columns(columns)
    return lead + [col for col in columns if col not in lead]


def render_fallback_table(
    data: pd.DataFrame, height: int, error: Exception, what: str
) -> None:
    """What a user sees when the grid could not be drawn: the frame, plainly.

    One definition, because there were three and they disagreed. This module's
    handler configured three columns for readability and stretched to the container;
    the two in ``results_view`` did neither, and named the tab in the error where this
    one named the table. A path a user only ever reaches on a bad day is exactly the
    path that drifts, so the differences were accidents of where each was written
    rather than decisions — issue #73 asked whether the handlers agreed about what a
    broken table looks like, and this is them agreeing.

    ``what`` names what failed, so the message stays specific: three call sites, three
    different things that can break.
    """
    st.error(f"❌ Error displaying {what}: {str(error)}")
    st.write("**Fallback: Simple table view**")

    # For fallback, use st.dataframe with column configuration for better display
    column_config = {}

    # Configure key columns for better display in fallback mode
    if "Clinical_Summary" in data.columns:
        column_config["Clinical_Summary"] = st.column_config.TextColumn(
            "Clinical Summary",
            width="medium",
            help="Consolidated clinical significance",
        )

    if "Hugo_Symbol" in data.columns:
        column_config["Hugo_Symbol"] = st.column_config.TextColumn(
            "Gene", width="medium"
        )

    if "Tumor_Sample_Barcode" in data.columns:
        column_config["Tumor_Sample_Barcode"] = st.column_config.TextColumn(
            "Sample ID", width="medium"
        )

    st.dataframe(
        data,
        height=height,
        use_container_width=True,
        column_config=column_config,
    )


#: The column a double-click writes into, hidden and kept out of the side bar's column
#: panel. It is not part of the report: it is added to the copy of the frame handed to the
#: grid, never to ``display_data`` itself, so the download beneath the grid — which reads the
#: list ``create_data_table`` returns — cannot pick it up.
_DOUBLE_CLICK_COLUMN = "_double_click_stamp"

#: ``st_aggrid``'s own name for the row identifier it sends alongside the cells, and which it
#: sets as the index of ``selected_rows``. The double-clicked row arrives as plain JSON, so
#: there it is a field like any other and has to be read — and then hidden — by hand.
#:
#: **It is a position, not a label**, and the name invites the opposite reading. ``st_aggrid``
#: builds it as ``[str(i) for i in range(len(frame))]`` before serialising the rows
#: (``st_aggrid/AgGrid.py:41``, streamlit-aggrid 1.1.6) and deletes it again afterwards, so it
#: is the row's **position** in the frame the grid was handed, as a string, and it carries no
#: trace of that frame's own index. Issue #310 is what reading it as a label cost: the frames
#: this app draws are the report's passed/failed splits, boolean-mask slices whose labels are
#: neither contiguous nor zero-based, and ``full_data.loc[2]`` on labels ``[2, 3, 5, 7]`` is
#: the *first* row where position 2 is the third. ``tests/test_row_recovery.py`` holds the
#: installed component to this, and issue #311 read every published release to find out what
#: it is holding against: the field is positional in **all nine**, 1.1.5 to 1.2.1.post2, and it
#: has to be — it feeds AG Grid's ``getRowId``, which requires a unique string, and a pandas
#: label guarantees neither. So the hazard is not an upgrade that starts sending labels. It is
#: the **rename** to ``::auto_unique_id::`` at 1.1.8, the release straight after the pin, which
#: makes the field *absent* rather than misleading: this app degrades to the match below, and
#: the contract test fails at import. The pin is the boundary, and ``requirements.txt`` is
#: where it is argued (issue #321).
#:
#: **The two routes received this value differently, and one of them never used it at all.**
#: Issue #309 measured the asymmetry, which matters here because it is invisible from the code
#: as it now stands. ``selected_rows`` sets its index from this same field and leaves the
#: ordinal a **string** (``st_aggrid/AgGridReturn.py:266``), so ``'5' in full_data.index`` was
#: ``False`` against a MAF's int64 labels — always — and the ``🔍 View details`` button had been
#: resolving variants by the three-column match all along. The double-click **cast** the string
#: to the frame's index dtype, added by #159 to make the lookup work, and that cast is what
#: turned a harmless miss into a silent wrong hit: in one Chromium run, on one row, the two
#: routes named different variants on screen at the same moment. The routes were never
#: symmetric, so do not "simplify" the two reads back together on the assumption that they
#: were. They agree now because :func:`_grid_position` reads both as positions — that agreement
#: is the fix, not a property of what arrives.
_AGGRID_INDEX_KEY = "__pandas_index"

#: Why a double-click writes a cell instead of just being reported.
#:
#: ``st_aggrid`` will report *any* AG Grid event named in ``update_on`` — its
#: ``attachStreamlitRerunToEvents`` calls ``api.addEventListener`` for each one — so
#: ``rowDoubleClicked`` reaches Python on its own, carrying the row and its index. That is
#: not enough, and the reason is measured rather than assumed, under issue #159:
#:
#: * A component's value **persists across reruns**, so a standing "the last event was a
#:   double-click" condition re-opens the dialog on every later interaction. Dismissing a
#:   dialog reruns nothing, so there is no close event to clear it on.
#: * The payload of a genuine second double-click on the same row is **byte-identical** to
#:   that persisted value — the whole 249KB of it, not merely the event — so no comparison
#:   can tell them apart. That also rules out ``AgGrid(callback=…)``: Streamlit only fires a
#:   widget callback when the value *changed*.
#: * A JsCode handler cannot stamp the outgoing event either. It does run — it can say so on
#:   the browser console — but ``st_aggrid`` has already cloned and sent the event by then.
#:
#: ``setDataValue`` dispatches a *fresh* ``cellValueChanged`` after this handler returns, and
#: that event carries the value written. So every double-click arrives with a stamp of its
#: own, and the stamp is the one-shot token: :func:`_fresh_double_click` consumes it on open.
#: Re-arming on a second subscribed event instead — ``rowClicked`` before the double-click —
#: was tried and rejected: Streamlit coalesces grid events under load (a subscribed
#: ``selectionChanged`` disappeared outright in that run), so the re-arm can be dropped and
#: the double-click silently do nothing.
_DOUBLE_CLICK_JS = """
function(event) {
    event.node.setDataValue(
        "%s", String(Date.now()) + "-" + String(event.rowIndex)
    );
}
""" % (
    _DOUBLE_CLICK_COLUMN,
)


def _fresh_double_click(grid_response, state, state_key: str):
    """Return the row of a double-click not yet acted on, or ``None``.

    The stamp is consumed here, at the moment the row is handed back, because that is the
    only moment there is: an ``st.dialog`` dismissed with Escape reruns nothing at all, so
    a token cleared on close would never be cleared.

    ``isinstance`` rather than truthiness on both reads, and not for tidiness: several tests
    stand ``AgGrid`` up as a ``MagicMock``, whose ``event_data`` is a truthy mock whose
    ``.get`` answers with another truthy mock. Read loosely, every one of those tests would
    look like a live double-click and open the variant dialog.
    """

    event = getattr(grid_response, "event_data", None)
    if not isinstance(event, dict):
        return None
    row = event.get("data")
    if not isinstance(row, dict):
        return None

    stamp = row.get(_DOUBLE_CLICK_COLUMN)
    if not stamp or state.get(state_key) == stamp:
        return None

    state[state_key] = stamp
    return row


def _grid_position(raw, full_data: pd.DataFrame) -> int | None:
    """The position in ``full_data`` that ``_AGGRID_INDEX_KEY`` names, or ``None``.

    One reader for both routes, because the two receive the same value by different means —
    the double-click as a field in the row's JSON, the button as the name of the row
    ``selected_rows`` hands back (``st_aggrid/AgGridReturn.py:266``) — and nothing else about
    them differs here.

    It arrives as a string and is answered against ``full_data``, which is sound because
    ``full_data`` is the frame the grid was drawn from: ``display_data`` is
    ``cleaned_data[visible_columns]``, the same rows in the same order, and the copy handed to
    ``AgGrid`` only adds a column. Sorting and filtering in the grid do not disturb it either
    — the value is a field on the row, written before the rows were sent, so it survives the
    reader reordering the view and is exactly *not* the row's position on screen.

    Out of range is refused rather than clamped, and negatives with it: ``iloc[-1]`` is a
    perfectly good pandas call that returns the last variant in the report, which is the kind
    of answer this function exists to stop giving.
    """

    if raw is None:
        return None

    try:
        position = int(raw)
    except (TypeError, ValueError):
        return None

    return position if 0 <= position < len(full_data) else None


def _full_row(
    partial: "pd.Series", full_data: pd.DataFrame, grid_position: int | None = None
) -> "pd.Series":
    """Recover the pristine row behind the partial one the grid handed back.

    Lifted out of the ``🔍 View details`` handler unchanged so the dialog is handed the same
    row whichever way it was opened, rather than the double-click growing a second opinion
    about which variant the user picked.

    ``grid_position`` is the component's own answer to *which row was that*, read by
    :func:`_grid_position`, and it is exact: it identifies one row of ``full_data`` and needs
    no columns to do it. Both call sites pass it, so on every live path the match below is not
    reached at all.

    The match is kept anyway, because issue #147's reason for it stands — the row makes a
    round trip through browser JSON and the identifier is the component's to send, not ours to
    promise — but it is no longer allowed to guess. Issue #310 is what guessing cost:

    * The mask began all-True and was narrowed only by the columns the partial row carried, so
      a partial carrying none of the three narrowed it by nothing and ``iloc[0]`` answered
      with **row 0 of the report**, dressed as the variant the reader picked.
    * With some of the three it narrowed to a *set* and took the first of it. That is the
      ordinary case rather than a corner: 37 of the 104 real MAFs on this machine repeat at
      least one ``(gene, chromosome, position)`` triple — one pooled file on 1711 of its 2361
      rows — because one variant called in several samples is one variant in all three of
      those columns. ``Tumor_Sample_Barcode`` is what tells those rows apart and it is not
      matched on.

    So a clause must have been applied and exactly one row may satisfy every clause; anything
    else returns ``partial``. That is a thinner dialog — the columns the grid was drawn with,
    and not the rest of the report's — but it is the row the reader actually picked, and this
    module's standing rule is that a wrong answer indistinguishable from a right one is worse
    than a short one. Nothing is refused *visibly*: saying so belongs to the dialog's surface,
    which issue #308 rules out redesigning here.
    """

    # Dropped here rather than at either call site, because *no* caller wants them and the
    # fallback below returns `partial` itself: a row the grid could not be matched back to
    # would otherwise show `_double_click_stamp` and `__pandas_index` as fields in the variant
    # dialog — both meaningless to a reader, and the kind of internal name that is not allowed
    # to reach the interface. The selection route carries the stamp too, and only the
    # double-click route carries the index, which ``selected_rows`` turns into an index
    # instead.
    partial = partial.drop(
        labels=[_DOUBLE_CLICK_COLUMN, _AGGRID_INDEX_KEY], errors="ignore"
    )

    # `is not None`, not truthiness: position 0 is the first row of every report and the one
    # falsy value in the range, so a truthiness test sends exactly that row down the match —
    # where a repeated triple answers it with a different sample's row.
    if grid_position is not None:
        return full_data.iloc[grid_position]

    mask = pd.Series(True, index=full_data.index)
    applied = 0
    for match_col in ["Hugo_Symbol", "Chromosome", "Start_Position"]:
        if match_col in partial.index and match_col in full_data.columns:
            mask &= full_data[match_col] == partial[match_col]
            applied += 1

    if applied:
        matches = full_data[mask]
        if len(matches) == 1:
            return matches.iloc[0]

    return partial


def _double_clicked_series(row: dict) -> "pd.Series":
    """Turn the double-clicked row's JSON into the partial row :func:`_full_row` takes.

    It used to do more, and issue #310 retired the rest of it. ``__pandas_index`` arrives as a
    string, so this cast it to the frame's own index dtype to stop it missing on a MAF's
    integer index — a conversion that made the lookup *hit*, on a label that names a different
    variant, and so turned a miss into a wrong answer. The value is a position; positions are
    read by :func:`_grid_position` and passed to :func:`_full_row` beside the row, and the
    series no longer carries an index label at all.

    It takes the row alone now — reconciling anything against the frame was the whole of what
    it did wrong — and stays a named step rather than being inlined so that what was retired
    has somewhere to be recorded.
    """

    return pd.Series(row)


def _render_aggrid_with_detail(
    display_data: pd.DataFrame,
    full_data: pd.DataFrame,
    height: int,
    key_suffix: str,
    visible_columns: list,
) -> None:
    """Render AgGrid table with variant detail panel. Designed to run inside @st.fragment."""

    # --- Configure AgGrid ---
    #
    # The grid is given a copy carrying one extra column, `_DOUBLE_CLICK_COLUMN`, which a
    # double-click stamps so that the click reaches Python exactly once. `display_data`
    # itself is left alone: it is what the fallback draws and what the row counts below are
    # taken from, and the download reads the column list this function returns, which the
    # stamp is not in.
    grid_data = display_data.assign(**{_DOUBLE_CLICK_COLUMN: ""})

    gb = GridOptionsBuilder.from_dataframe(grid_data)
    gb.configure_pagination(paginationAutoPageSize=True)
    gb.configure_side_bar()
    gb.configure_selection("single", use_checkbox=False)

    gb.configure_default_column(
        minWidth=80,
        resizable=True,
        sortable=True,
        filter=True,
        floatingFilter=True,
    )

    gb.configure_grid_options(
        enableRangeSelection=True,
        suppressHorizontalScroll=False,
        onRowDoubleClicked=JsCode(_DOUBLE_CLICK_JS),
    )

    # Hidden, and suppressed from the columns tool panel as well, so the side bar cannot
    # offer to show it. `editable=False` is the same statement the rest of the grid makes:
    # the stamp is written through the grid's API by the double-click handler, and no cell in
    # this table is typed into — the Notes column included, which is written in the dialog.
    gb.configure_column(
        _DOUBLE_CLICK_COLUMN,
        hide=True,
        editable=False,
        suppressColumnsToolPanel=True,
    )

    # Pinned in frame order, which is why APP_EXTRA_COLUMNS' order — Clinical_Summary,
    # Pathogenicity_Overview, Notes — is the order they appear in at the left edge, with
    # the user's own annotation columns after them.
    pinned_left = set(_leading_columns(visible_columns))

    for col in visible_columns:
        kwargs = {
            "width": _COLUMN_WIDTH_MAP.get(col, 120),
            "resizable": True,
            "sortable": True,
            "filter": "agNumberColumnFilter" if col in _NUMERIC_COLUMNS else True,
            "floatingFilter": True,
        }
        if col in pinned_left:
            kwargs["pinned"] = "left"
        gb.configure_column(col, **kwargs)

    gridOptions = gb.build()

    col_hash = hash(tuple(visible_columns))
    if key_suffix:
        table_key = f"aggrid_{key_suffix}_{col_hash}"
    else:
        table_key = f"aggrid_table_{col_hash}"

    # Key for the Pathogenicity_Overview circles. Nothing in it is written out by hand,
    # and the two halves were freed by different tickets for the same reason.
    #
    # The **glyphs** used to be one hand-written string naming six of the eight the circles
    # can draw, which is how ⚪ and 💊 came to be drawn beside a key that did not list them —
    # a user met an unexplained white circle immediately beside the ⬜ the key *does* define,
    # as "No data" (issue #100). They now read the constants they are drawn from, so a glyph
    # cannot outlive its explanation.
    #
    # The **sources** used to be all six on both arms, so a germline reader was given
    # `CancerVar` and `CIViC` — two tools that arm's `compute_keep` does not emit, so two
    # positions that could only ever be ⬜ — and a somatic reader got `InterVar` and
    # `RENOVO` on the same terms. Issue #95 narrowed the column rather than annotating the
    # dead positions: a key cannot drop a name without breaking the positional count, so
    # the column stopped drawing them instead. Both the abbreviations and their expansions
    # now come from the list the circles were actually drawn from.
    #
    # `overview_sources`, not a fresh `circle_sources` call: see this module's
    # `clinical_summary` counterpart for why the answer is carried rather than recomputed.
    #
    # The **axis clauses** follow the sources too. A position graded on something other than
    # pathogenicity needs a sentence saying so, and a MAF whose ESCAT column the arm does not
    # emit — or that simply does not carry one — draws no ES position for it to be about. There
    # are two such positions since issue #109, so which sources need a clause is asked of
    # `circle_axis_notes` rather than tested for here: this module draws no circles, and the
    # `any(source[0] == "ES" …)` it used to run was one hand-written exception per exceptional
    # source, sitting a long way from the mappings that make them exceptional.
    sources = st.session_state.get("overview_sources") or []
    if "Pathogenicity_Overview" in visible_columns and sources:
        abbreviations = " ".join(source[0] for source in sources)
        names = ", ".join(source[1] for source in sources)
        glyphs = "  ".join(
            f"{glyph} {label}" for glyph, label in pathogenicity_circle_legend()
        )
        axis_notes = "".join(f"<br>{note}" for note in circle_axis_notes(sources))
        st.markdown(
            '<div style="font-size:0.80em;color:#888;margin-bottom:4px;">'
            f"<b>Pathogenicity Overview</b> — {glyphs}"
            f"<br><b>Order:</b> {abbreviations} "
            f"({names})"
            f"{axis_notes}</div>",
            unsafe_allow_html=True,
        )

    # Only the render is guarded, and the guard ends where the render does.
    #
    # This `try` used to reach as far as `_show_variant_dialog`, so a failure while
    # opening one variant's details was reported as "Error displaying table" and
    # replaced a grid that had drawn perfectly well with a fallback for it. Catching
    # broadly is right *here* — AgGrid fails on data it cannot serialise, and there is
    # no narrower exception to name — but it is only right for the call it wraps.
    #
    # There is no grid response after a failed render, so there is no selection to
    # handle and nothing below this to run: the handler returns rather than falling
    # through to code that would read the name AgGrid never bound. That was the
    # original defect in a different form — the handler referred to `cleaned_data`, a
    # local of `create_data_table` several hundred lines away, so every AgGrid failure
    # raised `NameError` from inside the handler for it and lost the real error.
    # Issue #73 found it with `flake8 --select=F821`; #54 carried it here.
    try:
        # Create a container for the table to isolate state
        with st.container():
            grid_response = AgGrid(
                grid_data,
                gridOptions=gridOptions,
                data_return_mode="FILTERED_AND_SORTED",
                update_mode="SELECTION_CHANGED",
                # The double-click's cell stamp. `update_mode` is kept as well and still
                # reports selection, which is what the button below reads.
                update_on=["cellValueChanged"],
                fit_columns_on_grid_load=False,
                theme="streamlit",
                enable_enterprise_modules=True,
                height=height,
                width="100%",
                reload_data=True,
                allow_unsafe_jscode=True,
                key=table_key,
                try_to_convert_back_to_original_types=False,
            )
    except Exception as e:
        # `display_data` is the honest substitute: it is the frame this function was
        # handed and just failed to draw, in the columns the user asked for. `full_data`
        # would answer a failure to show forty columns by showing four hundred.
        render_fallback_table(display_data, height, e, "table")
        return

    # Here, not in the two results tabs: only this path has filter boxes to explain (#62).
    st.caption(
        "Number columns take a comparison, not just an exact value — type "
        "`>50` or `<0.01` in the box under the heading."
    )

    # Display interaction information if grid response is available
    if (
        grid_response
        and hasattr(grid_response, "data")
        and len(grid_response.data) > 0
    ):
        filtered_count = len(grid_response.data)
        total_count = len(display_data)
        if filtered_count != total_count:
            st.info(
                f"📊 Table filtered: Showing {filtered_count} of {total_count} variants"
            )

    # Variant detail panel — use full row data (all columns) for detail view
    has_selection = (
        grid_response
        and grid_response.selected_rows is not None
        and len(grid_response.selected_rows) > 0
    )

    if has_selection:
        selected_partial = pd.DataFrame(grid_response.selected_rows).iloc[0]

        if st.button(
            # Read `_view_details_label` for why no `chr` literal is written here (#149), and
            # why the string is composed there rather than in the two places that draw it.
            _view_details_label(selected_partial),
            key=f"view_detail_{key_suffix}",
        ):
            # The button's position arrives as the row's *name*: `st_aggrid` has already made
            # `_AGGRID_INDEX_KEY` the index of `selected_rows`, so there is no such field to
            # read. It is the same string the double-click gets as a field, and both go
            # through `_grid_position` so the two routes cannot come to different answers
            # about which variant the reader picked (#310).
            _show_variant_dialog(
                _full_row(
                    selected_partial,
                    full_data,
                    _grid_position(selected_partial.name, full_data),
                )
            )

    # Last, and after the button, so a run that opens the dialog has already drawn the grid
    # and everything under it: `st.dialog` takes over the screen, and what is behind it is
    # what the user comes back to on dismissing it.
    #
    # Nothing here reads the *selection*. A row is selected by the first click of the
    # double-click, so a dialog opened on "a row is selected" would re-open on every later
    # interaction while it stayed selected; opening on the stamp instead means one dialog per
    # double-click, and dismissing it leaves nothing behind to re-open it.
    double_clicked = _fresh_double_click(
        grid_response, st.session_state, f"{_DOUBLE_CLICK_COLUMN}_{key_suffix}"
    )
    if double_clicked is not None:
        _show_variant_dialog(
            _full_row(
                _double_clicked_series(double_clicked),
                full_data,
                _grid_position(double_clicked.get(_AGGRID_INDEX_KEY), full_data),
            )
        )
