"""What each filter setting is called on screen, and which arm it applies to.

The vocabulary the app already uses for its filter settings, lifted out of the widget
calls that were its only home so that a second surface can name a setting without
inventing a second name for it.

**Why this exists at all.** Issue #137 puts a recap of the filters beside the report, and
the Style rule forbids internal parameter names reaching the interface — so the recap
cannot render ``filter_cancervar``. The parameter page's own control labels are the
existing honest vocabulary, but they were literals inside ``st.multiselect`` calls, which
made them unreadable from anywhere else. A second hand-written set of names is the thing
this repository has been punished for repeatedly: two copies of one fact drift, and the
copy beside the control is the one that stays right (issue #79). So the labels move here
and the controls read them, which leaves exactly one place a setting is named.

**Why not the whole control table.** Only the parts a second surface needs are here — the
identity, the label, the arms, and how a value is spelled. The option lists live in
:mod:`config.vocabularies` and the help prose stays beside the widget, because help text
is written *for* the control and reads wrongly anywhere else.

**Rows, not keys.** ``filter_clinvar`` is drawn by two controls since issue #103 — the
pathogenicity calls and the other assertions — so a table keyed by parameter name could
not describe the screen. Each row therefore has its own ``id``, and the ClinVar pair
carry the term lists that partition one parameter between them.

This module is importable without Streamlit, like the rest of ``config/``.
"""

from dataclasses import dataclass

from config.vocabularies import (
    CLINVAR_OTHER_ASSERTION_TERMS,
    CLINVAR_PATHOGENICITY_TERMS,
)


@dataclass(frozen=True)
class ParamLabel:
    """One filter setting, as the interface names it.

    Attributes:
        id: stable identity for this row. It is *usually* also the attribute the filter's
            resolved settings carry the value under, which is what lets the ``keep`` rows
            be read with one ``getattr``. Three rows are exceptions and each is handled by
            name where it is spelled: the two ClinVar halves, which partition one
            attribute between them, and ``vaf_threshold_germline``, which is a *control's*
            name — the resolved settings hold one ``vaf_threshold`` per run, the germline
            spelling having been resolved away by then.
        label: the words the parameter page's own control uses. Any surface naming this
            setting uses this string, so a user reading it elsewhere can find the control.
        arms: the arms this setting applies to. The somatic guideline sources are absent
            from the germline filter's signature entirely, so naming them on that arm
            would describe a clause that does not exist.
        kind: how a value of this setting is spelled for a reader — see
            :func:`filters.run_recap.spell`.
        terms: for a row that holds one half of a parameter, the terms belonging to that
            half. Empty for every other row.
    """

    id: str
    label: str
    arms: tuple[str, ...]
    kind: str
    terms: tuple[str, ...] = ()

    def applies_to(self, arm: str) -> bool:
        return arm in self.arms


_BOTH = ("somatic", "germline")

#: Every setting the filter reads, in the order of the controls that set them: the Basic
#: tab, then the Clinical Classification tab, then Population Frequency. That order is not
#: decoration — a user reading the recap is on their way to change one of these, and the
#: reading order is the order they will meet the controls in.
#:
#: ``skip_civic`` is deliberately absent. It is in the parameter contract, no control sets
#: it, and ``components/variant_table.py`` records that it is always ``False`` in this app —
#: so a row for it would name a setting the user has never made and cannot change.
PARAM_LABELS: tuple[ParamLabel, ...] = (
    ParamLabel("min_depth", "Minimum read depth", _BOTH, "depth"),
    # Two labels for one concept, because the two arms' controls are two widgets with two
    # labels — and the threshold itself is spelled two ways in live code, which is why the
    # recap reads it off the filter's resolved settings rather than the dict.
    ParamLabel("vaf_threshold", "VAF Threshold", ("somatic",), "vaf"),
    ParamLabel("vaf_threshold_germline", "VAF Threshold (Germline)", ("germline",), "vaf"),
    ParamLabel(
        "exclude_classifications", "Exclude Variant Classifications", _BOTH, "exclude"
    ),
    ParamLabel(
        "skip_pathogenic", "Auto-retain pathogenic variants", _BOTH, "retention"
    ),
    # The gene controls are three widgets — a panel dropdown, a paste box and an upload —
    # that all resolve to one list, so no single control label names what the filter
    # received. This is the section's own heading, minus its emoji.
    ParamLabel("gene_selection", "Gene Filters", _BOTH, "genes"),
    ParamLabel(
        "cancervar_keep", "Keep CancerVar Classifications", ("somatic",), "keep"
    ),
    ParamLabel("civic_keep", "Keep CIViC Evidence Levels", ("somatic",), "keep"),
    ParamLabel("escat_keep", "Keep ESCAT Levels", ("somatic",), "keep"),
    ParamLabel(
        "intervar_keep", "Keep InterVar Classifications", ("germline",), "keep"
    ),
    ParamLabel("renovo_keep", "Keep RENOVO Classifications", ("germline",), "keep"),
    ParamLabel(
        "clinvar_pathogenicity",
        "Keep ClinVar Classifications",
        _BOTH,
        "keep",
        tuple(CLINVAR_PATHOGENICITY_TERMS),
    ),
    ParamLabel(
        "clinvar_other",
        "Keep Other ClinVar Assertions",
        _BOTH,
        "keep",
        tuple(CLINVAR_OTHER_ASSERTION_TERMS),
    ),
    ParamLabel(
        "max_freq_population", "Maximum Population Frequency", _BOTH, "frequency"
    ),
)

#: By id, for a caller that knows which setting it wants — the parameter page reads its
#: guideline control labels through this.
LABELS_BY_ID = {row.id: row for row in PARAM_LABELS}


def label_of(row_id: str) -> str:
    """The words the interface uses for one setting.

    A function rather than ``LABELS_BY_ID[row_id].label`` at each call site: what a caller
    wants is the label, not the row, and the walk to reach it is this module's business.
    The lookup raises either way — and raises *early*, because the parameter page's control
    table is built at import, so an id that no longer exists stops the page from loading
    rather than waiting for a user to open the tab that draws it.
    """
    return LABELS_BY_ID[row_id].label


def labels_for(arm: str) -> tuple[ParamLabel, ...]:
    """Every setting that applies on ``arm``, in control order."""
    return tuple(row for row in PARAM_LABELS if row.applies_to(arm))
