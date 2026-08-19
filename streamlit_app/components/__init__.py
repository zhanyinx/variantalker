"""The app's UI components, one module per kind.

| module | holds |
| --- | --- |
| `clinical_summary` | the two derived columns: one classification, and one circle per annotation source this arm and this file carry |
| `sidebar` | what file is open, where you are, and how to open another |
| `variant_table` | the grid, its columns, and the notes written onto it |
| `results_view` | the passed/failed/summary tabs of the Results section |
| `charts` | the summary dashboard's six Plotly charts |
| `variant_detail` | one variant, read out in full |
| `clinical_badges` | the clinical row's three annotations, and what it says when none draws |
| `acmg_evidence` | the ACMG/AMP criteria, and InterVar's report of them |
| `cbp_evidence` | the twelve AMP/ASCO/CAP criteria, and CancerVar's report of them |
| `cancervar_markers` | the markers behind CBP1/CBP2/CBP3, from CancerVar's own table |
| `predictor_context` | the raw predictor scores behind InterVar's PP3/BP4 and CancerVar's CBP10 |
| `alphamissense` | one in-silico predictor's own band, on its own three-class scale |
| `reference_scales` | the ClinGen SVI PP3/BP4 scores — a separate calibration that fed no verdict |
| `arm_notice` | what the app says when the loaded MAF and the selected arm disagree |

They depend on each other in one direction only:

    results_view ──▶ variant_table ──▶ variant_detail ─┬─▶ acmg_evidence ◀────┐
         │                 │                           │        ▲             │
         │                 │                           ├─▶ cbp_evidence ◀─────┤
         │                 │                           │                      │
         │                 │                           ├─▶ predictor_context ─┘
         │                 │                           │        ▲
         │                 │                           ├─▶ alphamissense ─┐
         │                 │                           ├─▶ reference_scales ◀┘
         │                 │                           ├─▶ clinical_badges
         │                 │                           └─▶ cancervar_markers
         ├──▶ charts       └──▶ clinical_summary ◀───────────────┐
         └──────────────────────────────────────────────────────┘

``clinical_badges`` is a **leaf**, and that is what it is for. Issue #204 gave the clinical row
three annotations that each need a vocabulary — ClinVar's star hierarchy, ESCAT's tiering-table
vintages, RENOVO's score — plus five computed empty-state sentences, and none of that wants a
Streamlit call to be asked a question. It imports no component and no ``st``: it turns a row into
badge values and sentences, and ``variant_detail`` renders them. That is what lets
``tests/test_clinical_badges.py`` drive the vocabularies directly rather than through a fake
``st``, which matters most for the two states no MAF on this machine reaches.

``alphamissense`` imports ``predictor_context`` for issue #190's **table shape** — ``_Reading``
and ``_table_html`` — and ``reference_scales`` for ``CLASSIFIER_INPUTS``, which is how its
provenance caption asserts that neither classifier read AlphaMissense instead of claiming it. Both
are one-way: neither imports back. It is a separate module rather than a fourth table inside
``predictor_context`` because that module's whole subject is *the scores a classifier actually
read*, and this is a predictor neither one did — the distinction issue #191 moved
``reference_scales`` out to keep, paid once more.

``cancervar_markers`` is the only module here that reads a **file** — the vendored copy of
CancerVar's marker table, which a MAF's ``Therap_list`` / ``Diag_list`` / ``Prog_list`` cell
holds line offsets into. It is its own module for exactly that reason: issue #198 kept the
file read out of ``cbp_evidence`` so that module keeps its single Streamlit seam, and kept it
out of ``variant_detail`` so it can be tested with no session and no pandas. It imports
nothing from this package, so it is a leaf.

``cbp_evidence`` imports ``acmg_evidence`` for one thing: the tier-cell palette. Issue #188
settled that a somatic CBP score cell reuses the germline ACMG colours, so that a clinician
reading either arm reads one colour language — and deriving the somatic palette from the
germline one, rather than re-typing its hex values, is what stops the two drifting apart.

``predictor_context`` imports both evidence modules for their **parsers** — ``parse_intervar``
and ``parse_cancervar``. Issue #190's section names the criterion its scores sit under, and it
reads that criterion out of the same evidence string the section above it drew, so the two cannot
disagree about what the tool recorded. It does not import ``variant_detail``, which imports it:
the one thing they share — how a contaminated cell is turned into text — lives in
``config.contaminated_columns.held_value`` instead of being reached for across the cycle.

``reference_scales`` imports **neither** evidence module, and that absence is the point. Issue
#191 established that its four columns are disjoint from the nine either classifier reads, so it
has no criterion to name and no evidence string to agree with — unlike ``predictor_context``,
which exists to sit under one. It takes the tool's *name* as an argument purely so its label can
say whose verdict it did not contribute to. If it ever needs a parser, the claim its label makes
has stopped being true.

    sidebar — depends on nothing here
    arm_notice — depends on nothing here; it reads `filters.arm_detection` and is handed
    its switch as a callback, so it does not reach back into the page that draws it

``results_view`` is handed a ``filters.run_recap.RunRecap``, whether it is stale, and a
re-apply callback, and renders the filter recap from them (issue #137). It imports neither
``filters`` nor the page: the recap is a value the *filter run* built, staleness is a
comparison the page makes, and re-applying is the page's own function — so this package
still imports no page, the rule ``summary_report`` set when the run report moved here.

**This package re-exports nothing, deliberately.** It used to open with
``from .ui_components import *`` and carry an ``__all__`` of twelve names. **Seven** of
those named functions with no caller anywhere in the app — leaving five that were
really read — and the ``__all__`` entry was the whole reason nothing reported them: a
name re-exported from a package reads as public API, so a dead function listed here
looks like one kept for its callers rather than one nobody calls. Six were deleted with
that change; the seventh, ``results_view.create_export_buttons``, was left standing for
issue #83 to decide beside the downloads it would replace — and #83 deleted it too, on
the ground that it never routed through ``variant_table.with_user_columns`` and so could
only have re-introduced the note-dropping export issue #67 had just fixed. All seven are
gone, and the package still re-exports nothing.

``config/__init__.py`` lost the same wildcard for the same reason in issue #54 — it
is how a name added in one place becomes visible in three. Import from the module
that defines a thing, so the import says where it lives.
"""
