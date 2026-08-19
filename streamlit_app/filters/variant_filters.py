"""The app's one filter decision, made by the pipeline's own code.

This module used to hold a re-implementation of ``bin/filter_variants.py``'s filter, and
the two had drifted apart in eleven measurable ways — the largest of them a germline
ESCAT clause the pipeline does not have, worth 540 rows and 51% of all attributed
divergence on the reference. It now holds no filter logic at all. :func:`apply_filters`
loads the parameters, hands the frame to the vendored pipeline functions, and labels
every row with the verdict those functions reached and the reason they reached it.

Parity is therefore **structural**: there is no second implementation to keep in step.
What is left here is the impedance matching the vendored code cannot do for itself —
turning a list of gene symbols into the file path its body expects, composing the app's
own population-frequency layer, and turning three Boolean masks into two columns and a
diagnostics channel.

The seam
--------
``apply_filters(maf, params) -> (DataFrame, Diagnostics)``

**Every** input row comes back, in input order, with its index preserved, plus:

* ``MAFigate_filter`` — ``PASS`` or ``NOPASS``, the pipeline's verdict;
* ``MAFigate_reason`` — which of the four cells decided the row.

Callers derive whatever splits they need. The two previous shapes (``apply_filters``
returning only the survivors, and ``apply_filters_with_separation`` returning a pair)
are gone: one entry point means one place the decision is made and one seam to test.

How the verdict is composed
---------------------------
Exactly as ``bin/filter_variants.py:445-447`` composes it::

    PASS = filter_patho | (filter_common & filter_specific)

The union is **unconditional**. The app does not branch on ``skip_pathogenic``; it passes
the flag down, and the vendored functions return an all-``False`` rescue mask when it is
set. Branching in the app was one of the divergences.

``filter_specific`` already carries the guideline block, the VAF gate and the gene list —
the vendored functions return it that way, and nothing here reaches inside it.

This module is deliberately free of Streamlit. It is the one seam the parity harness and
the ``bin/``-free net can both reach by calling it.

One polarity, stated plainly
----------------------------
``filter_variant_classification`` is an **exclude** list, because that is what the
pipeline's parameter of that name means and this module does not reinterpret it. Issue
#36 converted the widget and the SOFT and CLINICAL presets to match, so every producer of
this key now agrees with this consumer.

The interim state is worth recording, because it is what the diagnostics decomposition was
built to make visible. While the presets still held include lists, the app's default cut
was *degraded* rather than merely different: measured on ``somatic_reference.maf``, the
criteria path was **2 rows of 82** against 19 under the contract — and yet the union PASS
total was **20 either way**, because pathogenic retention re-admitted 18 of the rows the
inverted list dropped. Exactly the coincidence issue #28 warns about ("counts that
coincide"): the totals agreed while the reason every row was in the report had changed
completely.

There is deliberately **no compatibility shim** for this key, unlike the three below whose
*spellings* changed. Nothing can distinguish a legacy include list from a deliberate
exclude-list deviation — issue #28 measured that the parity default is legal as both — so a
shim would have to guess, and guessing wrong here is silent in the widening direction. The
fallbacks below are safe because a *name* is unambiguous; a semantic inversion is not.

The refusal, since #38
----------------------
The first thing :func:`apply_filters` does is refuse a MAF whose depth or VAF columns hold
values that are not numbers, naming every offending column and value.
``filters/numeric_columns.py`` owns it, including the ``ast`` derivation of *which* columns
those are — read out of the vendored source rather than written down here. It is a
**refusal, not a verdict**: a typed exception rather than an entry in
:attr:`Diagnostics.notes`, because the app has no report to attach a warning to. What it
replaces is the ``TypeError`` the vendored comparison used to throw straight through this
function, which was already parity with the pipeline and was already the right answer — it
just arrived as a stack trace naming neither the column nor the value.

The fill, since #39
-------------------
The other half of that boundary. A filter-input column that is **absent** rather than
unreadable is filled with a neutral value and the frame is filtered anyway, so an incomplete
MAF is still usable — ``filters/absent_columns.py`` owns which columns those are, what stands
in for each, and what to say about it. Two properties of the wiring are this module's, and
both are load-bearing:

* the fill goes onto a **throwaway copy** that the vendored calls and the app's own frequency
  layer receive, while :func:`_label` labels the frame the caller passed. A filled ``tumor_f``
  of 1.0 reaching the grid would read as a genuine 100% VAF — the app lying to the user about
  their own data — so the fill reaches the masks and nothing else;
* filling is **off parity by construction**. The pipeline raises ``KeyError`` on these files,
  so a filled run has no pipeline verdict to match, and the warning says so in those words.

A missing ``Hugo_Symbol`` is not in that set and is handled here instead: it is the one
filter input the *app* guards, because the app writes the gene file whose existence decides
whether the vendored body reads the column at all. See :func:`_gene_file`.

The frequency extra, since #37
------------------------------
The population-frequency filter is the app's own — the pipeline has none at any value, and
suppressing common polymorphisms is a real clinical need the report does not serve. It is
also the only filter logic left in this module, and it is composed so that having it
available cannot cost parity::

    final_pass = pipeline_pass & (frequency acceptable | ClinVar-pathogenic)

Three properties, each of which was a measured decision rather than a preference:

* It sits **on top of** the pipeline's verdict. It is not passed into the vendored
  functions and their signatures are untouched. Joining it into the common filters instead
  was measured to leave 115 of the 136 disputed somatic rows in the report anyway, because
  the pathogenic union re-admits them — so that composition costs criteria-path rows while
  barely changing what the user sees.
* Pathogenic calls are **exempt**, reusing the vendored ClinVar function *and* the vendored
  constant that pathogenic retention already uses, so there is no second list to drift.
  Masking the union without the exemption drops a genuinely pathogenic low-penetrance
  allele for being common: on the reference the exemption keeps exactly one variant,
  ``SERPINA1`` chr14:94847262 ``Pathogenic/Pathogenic,_low_penetrance|other`` at 0.0233,
  and drops 132 rows that are ClinVar ``Benign``, ``drug_response`` or ``Conflicting``. The
  115 rows the joined composition would have kept are five variants — one
  ``drug_response`` and four ``Benign``, at frequencies up to 0.9949 — every one of them
  surviving on a potential-significance CancerVar tier.

  Issue #28 states these as "111 of 115 disputed rows". Re-running issue #16's own
  instrument at the app's actual rule shows 115 is the *disagreement* count and 136 the
  disputed count; the spec transposes the two and the conclusion is unaffected. Measured
  against this code over all 100 reference MAFs in ``docs/wayfinder/issue-37/``.
* It is **neutral at the parity default** of 1.0, algebraically and not by the guard:
  frequencies lie in [0, 1], so the comparison is a tautology there. Asserted that way with
  the short-circuit bypassed, and separately end-to-end.

A **missing** frequency cannot sink a row, and a row no panel measured passes. That is the
opposite of the rule for depth and VAF, and deliberately so — but a missing frequency does
not *rescue* a row a populated column calls common either, which is #122's change. See
:func:`frequency_mask`.

The gene list, since #34
------------------------
The tokeniser, the token validation and the named panels live next door in
``filters/gene_lists.py``, which is Streamlit-free and pandas-free so that the parameter
page and this seam can share one parser. Everything gene-shaped that arrives here — a
list, a legacy comma string, a one-symbol-per-line paste — goes through
``parse_gene_list`` before the adapter writes anything, so the running app and the parity
harness tokenise identically. What is owned *here* is the part that needs the MAF: which
requested symbols the frame does not carry, reported through
:attr:`Diagnostics.notes`.
"""

from __future__ import annotations

import os
import tempfile
from contextlib import contextmanager
from dataclasses import dataclass, field
from typing import Any, Iterator, Mapping

import pandas as pd

from filters.absent_columns import FillPlan, plan_fills
from filters.gene_lists import (
    GeneSelection,
    missing_symbols,
    name_tokens,
    parse_gene_list,
)
from filters.notes import Note
from filters.numeric_columns import require_readable_numerics
from vendor.pipeline_filters import common_filters as pipeline_common_filters
from vendor.pipeline_filters import germline_filters as pipeline_germline_filters
from vendor.pipeline_filters import somatic_filters as pipeline_somatic_filters
from vendor.pipeline_utils import CLINVAR_PATHO, has_clinvar_term

#: The verdict column, and its two values. Named for the app rather than ``filter``,
#: which is what the pipeline calls its own column: a MAF that has already been through
#: the pipeline carries ``filter``, and a user comparing the two needs to see both.
MAFIGATE_FILTER = "MAFigate_filter"
PASS = "PASS"
NOPASS = "NOPASS"

#: The reason column, and the four cells it can name. These are the four cells of the
#: union's decomposition, so the reason is a *rearrangement* of the masks that made the
#: verdict rather than a second opinion about it.
MAFIGATE_REASON = "MAFigate_reason"
REASON_CRITERIA = "criteria"
REASON_RESCUE = "pathogenic_rescue"
REASON_BOTH = "both"
REASON_REJECTED = "rejected"

#: The population-frequency columns the app's own extra reads, in no significant order.
#: The pipeline has no frequency filter at any value, so this list has no pipeline
#: counterpart to be guarded against — unlike every other list in this module, which is
#: imported from the vendored code rather than written down. It is, however, the *only*
#: statement of what counts as a frequency column: ``page_modules/data_loading.py`` reads
#: this list to decide whether a MAF has any, so the load-time account and the filter
#: cannot come to disagree about it.
#:
#: **Every entry is a panel's opinion about one variant, and the two pre-QC ``_raw``
#: columns are opinions of the same kind** — issue #126's decision, over 164 real MAFs
#: (104 somatic, 60 germline) through the whole filter, by
#: ``docs/wayfinder/issue-126/measure_shapes.py``. Three things settled it.
#:
#: That harness reproduces #122's per-arm deltas exactly — 435, 91, 0, 0 — over 164 files
#: where #122 reported 171, and the difference is corpus selection rather than a changed
#: answer: this one requires the arm's own annotation markers and drops a MAF the app itself
#: refuses. The absolute report totals therefore sit a little below #122's, and the deltas,
#: which are what a shape is judged on, are identical. **Shape A is asserted against the
#: shipped ``apply_filters`` on every file, 0 mismatches**, which is why the other three
#: columns of that table are evidence at all.
#:
#: The minimum rule already makes ``_raw`` deferential. Where both gnomAD exome columns
#: are populated the ``_raw`` value calls a variant common that the adjusted value does
#: not on **1,052 germline and 1,549 somatic rows** at 1% — and :func:`frequency_mask`
#: keeps every one of them, because the lower value decides and the adjusted one is the
#: lower in 1,052 of those 1,087 germline disagreements. The reverse case, the adjusted
#: column being the stricter of the two, is 35 germline and 32 somatic rows. So counting
#: ``_raw`` does not let the weaker column veto anything: it speaks only when nothing
#: better does.
#:
#: When nothing better does, it is the *only* gnomAD value there is. Six of the 164 files
#: carry a ``_raw`` column and no adjusted one — an annotation vintage that emits
#: ``gnomAD_exome_AF_raw`` and **not** ``gnomAD_exome_AF``, so the adjusted column is
#: *absent*, not blank (present-but-empty is 2 rows in the whole corpus). Three of those six
#: hold a row this decides, and they hold **all** of them: **435** Broad and **91**
#: Stringent germline removals resting on a ``_raw`` value alone. Dropping the ``_raw`` names
#: would leave the app with no gnomAD opinion at all on any of the six.
#:
#: **On the somatic arm the question is inert: 0 such removals at either preset.** Somatic
#: files in this corpus carry both columns, so the adjusted value always speaks and always
#: wins — the 1,549 disagreements above are a different claim, and none of them moves a row.
#: This decision is therefore paid for entirely on the germline arm.
#:
#: And that shape is the pipeline's own report. ``bin/filter_variants.py:384-392`` appends
#: ``gnomAD_exome_AF_raw`` **in preference to** ``gnomAD_exome_AF``, and the same for the
#: genome pair — so a MAF that has been through the pipeline carries the pre-QC column and
#: not the adjusted one. An app that refused to read ``_raw`` would be blind to gnomAD on
#: the output of the tool it re-implements.
#:
#: The third shape issue #126 weighed — ``_raw`` consulted only where the adjusted column
#: of the same build is blank or absent — was measured and is **not** a middle ground: it
#: keeps all 435 and 91 removals, since those rest on absence, and additionally drops the
#: 27 rows ``_raw`` rescues. It is today's rule minus 27 rows, not a softer one.
#:
#: **``MAX_AF`` is gone, and it is the one entry that was never a panel.** It is a maximum
#: across subpopulations, so under the minimum rule it could only ever rescue a row, never
#: remove one, unless it was the sole populated column — and issue #122 refused the max
#: rule partly because ``MAX_AF``'s presence would have made the cut mean "no subpopulation
#: above 1%". It is also absent from **all 164** real MAFs measured, so nothing it could
#: have done has ever been exercised. Removing it moved **0** report rows on either arm at
#: either preset, measured rather than argued.
#:
#: The two ``gnomAD_genome_*`` names stay on the other side of that line: they are panels,
#: the pipeline's keep list already anticipates them, and a whole-genome run would carry
#: them. They are unobserved over a 2022–2024 corpus, which is not the same as impossible.
FREQUENCY_COLUMNS = (
    "gnomAD_exome_AF",
    "gnomAD_genome_AF",
    "gnomAD_exome_AF_raw",
    "gnomAD_genome_AF_raw",
    "Freq_ExAC_ALL",
    "Freq_esp6500siv2_all",
    "Freq_1000g2015aug_all",
)

#: The ClinVar significance column the pathogenic exemption reads. The vendored filters read
#: this same column on **both** arms and without a guard, so it is one of the columns
#: :mod:`filters.absent_columns` fills — and the exemption is computed on the *filled* frame,
#: where it is therefore always present. Its fill is NaN, on which ``has_clinvar_term`` is
#: False, so an absent ClinVar column exempts nothing rather than exempting everything.
#:
#: This *is* a name the vendored code also spells, which normally earns a guard here rather
#: than a copy — ``vendor/README.md`` records a hand-copied ``KEEP`` that drifted to 39
#: entries against the pipeline's 45. What makes this one different is that a drift cannot
#: be silent: the vendored calls subscript the same literal six times and run *first*, so a
#: renamed column raises ``KeyError`` out of them before this constant is ever consulted.
#: The failure mode a guard protects against — the app quietly filtering on something the
#: pipeline no longer reads — is unreachable. It is guarded anyway, cheaply, by
#: ``tests/test_filter_app_extras.py``, which asserts the vendored source still subscripts
#: this exact name.
CLINVAR_SIGNIFICANCE = "ClinVar_VCF_CLNSIG"

#: The gene-symbol column the vendored gene clause matches against. Spelled here for the same
#: reason as :data:`CLINVAR_SIGNIFICANCE` and guarded the same way, by
#: ``tests/test_filter_app_extras.py`` asserting the vendored source still subscripts this
#: exact name — the app has to test for the column's *presence* before deciding whether to
#: write a gene file at all, and a presence test is a place a copied name can rot silently
#: rather than raising.
GENE_SYMBOL = "Hugo_Symbol"

#: What the vendored functions take to mean "no gene filter". Their bodies compare the
#: argument against this literal and then call ``os.path.exists`` on it.
NO_GENE_FILE = "null"


@dataclass(frozen=True)
class Diagnostics:
    """The four-cell decomposition of the union, plus anything worth telling the user.

    The cells come from the same three masks that produced the verdict, so they cannot
    disagree with it, and they **partition** the frame: every row lands in exactly one,
    and the four sum to the row count. That is asserted as a partition rather than
    against recorded numbers, because the property is what matters and a second set of
    committed numbers would be a second thing to keep in step.

    What this replaces is worth naming: 28 diagnostic strings, of which the five funnel
    counts were measured wrong on their own terms — the classification count read 68,593
    where the pipeline removes 20,386, and the depth count read 41 where the truth is
    261, because it was measuring the wrong column. A decomposition derived from the
    deciding masks cannot be wrong in that way.

    The decomposition also reports something nothing reported before: how much of the
    report exists *only* by pathogenic rescue. On the reference that is 157 of 818
    germline rows — 19.2% — against 3 of 411 somatic.
    """

    rows: int
    criteria_only: int
    rescue_only: int
    both: int
    rejected: int
    #: Everything the filter has to say about this run, in display order, each carrying the
    #: level it should be drawn at. Renamed from ``warnings`` by issue #151, when the channel
    #: stopped being all-warnings: it now carries notes about a run that did exactly what was
    #: asked, which on a complete MAF on its own arm is the **only** thing in here.
    notes: tuple[Note, ...] = field(default_factory=tuple)
    #: Filter-input columns this MAF did not carry, which were filled neutrally so a report
    #: could be produced (#39). Non-empty means **this report is off parity by construction**:
    #: the pipeline raises on the same file, so there is no verdict for it to match.
    #:
    #: Structured as well as narrated, for the same reason
    #: :attr:`~filters.numeric_columns.UnreadableNumericColumns.offenders` is: a caller that
    #: needs to render these differently — the UI escalates them, the parity harness records
    #: them — should not have to parse prose to find out which columns they were.
    filled_columns: tuple[str, ...] = field(default_factory=tuple)
    #: The subset of :attr:`filled_columns` that feeds the pipeline's pathogenic retention.
    #: Non-empty means rows are **missing** from the report rather than merely re-cut, because
    #: the fill has emptied the unconditional safety net as well as one guideline term.
    degraded_columns: tuple[str, ...] = field(default_factory=tuple)

    @property
    def passed(self) -> int:
        """The **union** PASS count — not the criteria path, which is a different cell.

        Spelled out because these two numbers are routinely confused: on the somatic
        reference the criteria path is 408 and the union is 411, and quoting "the
        baseline" without saying which left a real discrepancy unreconciled for several
        tickets.
        """
        return self.criteria_only + self.rescue_only + self.both

    @property
    def criteria_path(self) -> int:
        """Rows meeting the filter criteria, rescued or not. The *other* cell."""
        return self.criteria_only + self.both

    def cells(self) -> dict[str, int]:
        """The partition, as a mapping. Sums to :attr:`rows`."""
        return {
            "criteria_only": self.criteria_only,
            "rescue_only": self.rescue_only,
            "both": self.both,
            "rejected": self.rejected,
        }


@dataclass(frozen=True)
class _Settings:
    """The parameters the vendored functions need, resolved from the app's dict.

    Resolution is separated from application so that the awkwardness of reading a
    parameter dict — which spellings are accepted, what a missing key means — is not
    tangled up with the decision itself.

    Two spellings are accepted per concept while the rest of the parity work lands.
    ``config/pipeline_params.py`` is the contract: one ``vaf_threshold``, one
    ``filter_genes``, ``skip_pathogenic`` with the pipeline's polarity. The app now
    *opens* on the contract (issue #36) and a saved parameter file is migrated onto it on
    the way in (issue #40), so the file case is closed. Two live sources still write the
    older spellings, and they are what keeps the fallbacks here:

    * the SOFT and CLINICAL presets carry ``keep_pathogenic``, and the germline pair carry
      ``vaf_threshold_germline`` with no ``vaf_threshold`` at all — so dropping that
      fallback would read their threshold as ``0.0`` and widen the report;
    * the germline VAF widget still writes ``vaf_threshold_germline``, which is why the
      migration sets both keys on that arm rather than collapsing to one.

    Delete the fallbacks with the ticket that moves the widget and the presets onto the
    contract's single spelling; until then this is the transitional code, and it is all of
    it.
    """

    arm: str
    min_depth: int
    exclude_classifications: list
    vaf_threshold: float
    clinvar_keep: list
    cancervar_keep: list
    civic_keep: list
    escat_keep: list
    intervar_keep: list
    renovo_keep: list
    #: The parsed gene list, not a bare list of symbols: the parse also produced the
    #: tokens it threw away and whether any restriction survived, and both have to reach
    #: the user. Carrying only ``symbols`` is what made "nothing usable was given" look
    #: identical to "no gene filter was asked for".
    gene_selection: GeneSelection
    skip_civic: bool
    skip_pathogenic: bool
    max_freq_population: float

    @classmethod
    def from_params(cls, params: Mapping[str, Any]) -> "_Settings":
        arm = params.get("sample_type") or "somatic"
        if arm not in ("somatic", "germline"):
            raise ValueError(
                f"sample_type must be somatic or germline; provided {arm!r}"
            )
        return cls(
            arm=arm,
            min_depth=params.get("min_depth", 0),
            # An *exclude* list, which is what the pipeline's parameter of this name
            # means: these classifications are dropped and everything else is kept. The
            # app's control says so too since issue #36 — "Exclude Variant
            # Classifications", defaulting to the config's three values.
            exclude_classifications=_terms(
                params.get("filter_variant_classification")
            ),
            vaf_threshold=_vaf_threshold(params, arm),
            clinvar_keep=_terms(params.get("filter_clinvar")),
            cancervar_keep=_terms(params.get("filter_cancervar")),
            civic_keep=_terms(params.get("filter_civic")),
            escat_keep=_terms(params.get("filter_escat")),
            intervar_keep=_terms(params.get("filter_intervar")),
            renovo_keep=_terms(params.get("filter_renovo")),
            gene_selection=_gene_selection(params, arm),
            skip_civic=bool(params.get("skip_civic", False)),
            skip_pathogenic=_skip_pathogenic(params),
            max_freq_population=params.get("max_freq_population", 1.0),
        )


def apply_filters(
    maf: pd.DataFrame, params: Mapping[str, Any]
) -> tuple[pd.DataFrame, Diagnostics]:
    """Label every row of ``maf`` with the pipeline's verdict and its reason.

    Args:
        maf: an annotated MAF, loaded by ``utils.read_maf`` — which is the pipeline's own
            reader. A frame loaded any other way may carry string dtypes, and the
            vendored comparisons raise on those rather than deciding anything.
        params: the filter parameters. See :class:`_Settings` for the spellings accepted.

    Returns:
        ``(labelled, diagnostics)``. ``labelled`` is ``maf`` with :data:`MAFIGATE_FILTER`
        and :data:`MAFIGATE_REASON` appended — every input row, in input order, index
        preserved, so a caller derives the passed and failed frames by masking on the
        verdict column. The two new columns come *last*, leaving the pipeline's own
        column order recognisable at a glance.

    Raises:
        UnreadableNumericColumns: if a depth or VAF column is present but holds a value
            that is not a number. The pipeline raises ``TypeError`` on exactly these
            files; this is the same non-verdict, naming every offending column and value.

    An **absent** filter-input column does not raise. It is filled neutrally on a throwaway
    copy and named in :attr:`Diagnostics.filled_columns`, which puts that run off parity by
    construction — see the module docstring.
    """
    # First, before the parameters are even read. The pipeline's raise does not depend on
    # any threshold, so neither can the refusal — and a caller with both an unreadable MAF
    # and a bad parameter dict needs to hear about the file, which is the thing they can
    # see. Refusing at *load* time instead would be worse still: the user would never reach
    # the screen that could show them which values stopped it.
    require_readable_numerics(maf)

    settings = _Settings.from_params(params)

    # Which columns are missing depends on the arm, so this cannot precede the parameters the
    # way the refusal does. Everything downstream reads `masked`; only `_label` and the gene
    # report read `maf`, and that split is what keeps a filled value out of the user's sight.
    plan = plan_fills(maf, settings.arm, skip_pathogenic=settings.skip_pathogenic)
    masked = plan.frame_for_masks(maf)

    notes: list[Note] = list(plan.warnings())
    notes.extend(_gene_notes(maf, settings.gene_selection))

    common = pipeline_common_filters(
        masked, settings.min_depth, settings.exclude_classifications
    )

    # The gene adapter must wrap the *call*, not just the write: the vendored body checks
    # for the file itself, so a file removed before the call is a file that was never
    # there, and the gene filter silently widens to everything.
    with _gene_file(_gene_symbols(maf, settings.gene_selection)) as genes_path:
        if settings.arm == "somatic":
            specific, rescue = pipeline_somatic_filters(
                masked,
                vaf=settings.vaf_threshold,
                somatic_genes=genes_path,
                cancervar_keep=settings.cancervar_keep,
                civic_keep=settings.civic_keep,
                escat_keep=settings.escat_keep,
                clinvar_keep=settings.clinvar_keep,
                skip_civic=settings.skip_civic,
                skip_pathogenic=settings.skip_pathogenic,
            )
        else:
            # No ESCAT argument, and that is the point: the germline guideline block is
            # InterVar | ClinVar | RENOVO. The clause the app used to OR in here was the
            # largest divergence on real data.
            specific, rescue = pipeline_germline_filters(
                masked,
                vaf=settings.vaf_threshold,
                germline_genes=genes_path,
                intervar_keep=settings.intervar_keep,
                renovo_keep=settings.renovo_keep,
                clinvar_keep=settings.clinvar_keep,
                skip_pathogenic=settings.skip_pathogenic,
            )

    criteria = common & specific

    # The app's own extra, composed on top of the verdict the pipeline just reached:
    #
    #     final_pass = (criteria | rescue) & (frequency_ok | pathogenic)
    #
    # It is applied to the two deciding masks rather than to their union, and that is the
    # *same* mask on the verdict — `(a & k) | (b & k)` is `(a | b) & k`. What it buys is
    # that the four cells below are still derived from the masks that decided the row: mask
    # the union afterwards and a dropped row reads NOPASS while its reason still claims it
    # met the criteria, which is a report contradicting itself.
    #
    # The guard is a short-circuit, not the neutrality argument. At 1.0 the mask is all-True
    # by algebra because frequencies lie in [0, 1], which is what the tests assert with this
    # branch bypassed — the guard only saves two passes over the frame.
    if settings.max_freq_population < 1.0:
        frequency_ok = frequency_mask(maf, settings.max_freq_population)
        if frequency_ok is None:
            # A ``WARNING``, and ``⚠️`` to match it since issue #151: the filter the user chose
            # did not run, so this report is not the one they asked for. See `filters.notes` for
            # why the glyph moved off #150's `ℹ️`.
            notes.append(
                Note.warning(
                    "⚠️ Population frequency filter: no frequency columns in this MAF, so "
                    "the filter was skipped"
                )
            )
        else:
            # `masked`, not `maf`: the exemption reads ClinVar_VCF_CLNSIG, which is one of the
            # columns that can have been filled. On the *frequency* mask above the distinction
            # does not arise — no frequency column is a filter input, so none is ever filled.
            exempt = pathogenic_exemption(masked)
            keep = frequency_ok | exempt
            passed = criteria | rescue
            removed = int((passed & ~keep).sum())
            exempted = int((passed & ~frequency_ok & exempt).sum())

            criteria, rescue = criteria & keep, rescue & keep

            if removed or exempted:
                notes.append(
                    _frequency_note(settings.max_freq_population, removed, exempted)
                )

    cells = _cells(criteria, rescue)
    return _label(maf, cells), _diagnose(cells, notes, plan)


# ---------------------------------------------------------------------------
# Turning masks into columns
# ---------------------------------------------------------------------------

#: Which reason each cell of the decomposition writes into :data:`MAFIGATE_REASON`.
#: One mapping, consulted by the labeller *and* by :func:`decomposition`, so the column
#: and the channel cannot come to name the same cell differently.
_CELL_REASONS = {
    "criteria_only": REASON_CRITERIA,
    "rescue_only": REASON_RESCUE,
    "both": REASON_BOTH,
    "rejected": REASON_REJECTED,
}


def _cells(criteria: pd.Series, rescue: pd.Series) -> dict[str, pd.Series]:
    """The union's four-cell decomposition, as four disjoint masks.

    Computed **once** and then used for both the row labels and the counts, so the verdict
    column, the reason column and the diagnostics channel are three readings of one
    Boolean rearrangement rather than three re-derivations that could drift apart.

    Disjoint and exhaustive by construction: the four are the truth table of two
    predicates, so every row falls in exactly one and the counts sum to the row count.
    """
    return {
        "criteria_only": criteria & ~rescue,
        "rescue_only": rescue & ~criteria,
        "both": criteria & rescue,
        "rejected": ~criteria & ~rescue,
    }


def _label(maf: pd.DataFrame, cells: dict[str, pd.Series]) -> pd.DataFrame:
    """``maf`` with the verdict and reason columns appended.

    ``assign`` rather than an explicit ``copy()`` because it does not mutate the caller's
    frame, which is Streamlit session state. It is *not* free: measured on the 427-column
    somatic fixture under pandas 2.3.1 with copy-on-write off, the returned frame shares
    no blocks with the input, so the data is duplicated once. That is still strictly less
    work than the implementation this replaces, which called ``data.copy()`` on every
    invocation *and* mutated its argument anyway by adding a ``CancerVar`` column in
    place. Enabling copy-on-write would make this genuinely cheap; that is a project-wide
    decision, not this seam's to take.
    """
    verdict = pd.Series(NOPASS, index=maf.index, dtype=object)
    verdict[~cells["rejected"]] = PASS

    reason = pd.Series(index=maf.index, dtype=object)
    for cell, label in _CELL_REASONS.items():
        reason[cells[cell]] = label

    return maf.assign(**{MAFIGATE_FILTER: verdict, MAFIGATE_REASON: reason})


def decomposition(reasons: pd.Series) -> dict[str, int]:
    """The four cells, re-derived from a :data:`MAFIGATE_REASON` column.

    For callers holding the labelled rows but not the :class:`Diagnostics` that came with
    them — the UI, which has already split the frame in two and stored the halves.
    Counting the reason column is not a second opinion about the decision: the column was
    written from the deciding masks, so this recovers the same partition rather than
    recomputing it.

    Takes the *column*, not the frame, because that is all it needs and the difference is
    not academic: the caller has two frames to combine, and combining them whole would
    copy 427 columns of every row on each rerender to read one of them.
    """
    return {
        cell: int((reasons == label).sum()) for cell, label in _CELL_REASONS.items()
    }


def _diagnose(
    cells: dict[str, pd.Series], notes: list[Note], plan: FillPlan
) -> Diagnostics:
    """The four cells counted, straight off the masks that made the decision.

    The plan travels as structured fields *and* as prose in ``notes``, and the two come
    from the same object rather than being assembled twice — a caller that renders the columns
    and a caller that renders the text cannot end up naming different ones.
    """
    counts = {cell: int(mask.sum()) for cell, mask in cells.items()}
    return Diagnostics(
        rows=len(cells["rejected"]),
        notes=tuple(notes),
        filled_columns=plan.filled,
        degraded_columns=plan.degraded,
        **counts,
    )


# ---------------------------------------------------------------------------
# The app's own extra
# ---------------------------------------------------------------------------


def frequency_mask(maf: pd.DataFrame, threshold: float) -> pd.Series | None:
    """Rows whose population frequency is acceptable, or ``None`` if the MAF has no
    frequency columns to judge them by.

    Answers one question, so the caller decides what "no columns" should mean rather than
    receiving a pre-formatted warning alongside a mask, and composes the exemption itself.

    Public, deliberately, though only :func:`apply_filters` calls it in the app. The
    neutrality claim this whole layer rests on is a claim about *this expression* — that it
    is all-True at 1.0 with the caller's ``< 1.0`` short-circuit bypassed — and a test that
    restates the expression to assert it is a test of the restatement. The parity suite
    previously kept its own copy of both this loop and :data:`FREQUENCY_COLUMNS`; now it
    calls this.

    The rule, in one line: **a row passes if the minimum over the columns that actually
    measured it is at or below the threshold, and a row no column measured passes.** A
    blank is not a measurement, so it contributes nothing either way.

    That minimum is also what makes :data:`FREQUENCY_COLUMNS`'s membership safe to be
    generous with, and issue #126 is the measurement of it: **the weakest column cannot
    veto a row.** Where two columns disagree about the threshold the more permissive one
    decides, so a pre-QC value that calls a variant common is overruled by an adjusted value
    that does not. A column earns its place in that list by being able to *rescue* a variant
    nothing else speaks for; it takes the absence of every other column for it to remove one.
    The measurement, and what taking the two ``_raw`` names out would cost, is recorded on
    :data:`FREQUENCY_COLUMNS` and is deliberately not restated here — one set of numbers, in
    the one place the decision they support is made.

    A blank still cannot *sink* a row, which is #28's and #37's decision and is untouched
    here. Missing rates run 13–35% per column, and a panel that never saw a variant is not
    evidence the variant is common; dropping on a blank would invert the filter and
    discard rare pathogenic variants for being rare. That is the opposite of the
    missing-value rule for depth and VAF, where a blank means *we cannot tell whether this
    call is sound* and the pipeline drops the row. Here it means *not seen in this panel*,
    and the whole-row fallback keeps a variant no panel has an opinion on.

    What issue #122 changed is the other half: a blank no longer *rescues* a row that a
    populated column calls common. Until then the loop read ``(v <= threshold) | v.isna()``,
    which is a blank counted as zero — so one empty column carried a row whatever the
    populated ones said, and the app reported variants at 100% allele frequency.

    **Measured before the change**, over 171 real MAFs (108 somatic, 63 germline) through
    this whole filter, the blank rescue was holding **4,167 rows across 534 distinct
    variants** in Broad reports — a third of the somatic Broad report. The top of that list
    was not marginal: ``ATM`` chr11:108183167 at gnomAD 1.0, ExAC 1 and 1000G 1, in every
    report only because ``Freq_esp6500siv2_all`` was ``.``; ``PTEN`` chr10:89623861 at
    1000G 1.0; ``DHFR`` at 0.86; ``HLA-A`` at 0.57. #114's ``HMGCR`` and ``NOS3`` — the two
    variants that opened #122 — were ranks 3 and 1 of those 534, so they were the visible
    edge of the mechanism rather than the mechanism.

    Whole-filter cost, Broad then Stringent: somatic 9,034 → 6,162 and 1,940 → 1,884;
    germline 12,873 → 11,578 and 1,287 → 1,142. **Every removed row is one where every
    panel that measured the variant put it above the threshold** — by construction, so no
    variant any panel calls rare is affected, and 2,567 of the 4,167 sat above 25%.

    No ClinVar-pathogenic row moves: :func:`pathogenic_exemption` still composes over this.
    303 removed rows did carry a strong call that exemption deliberately does not cover —
    252 RENOVO ``HP Pathogenic``, 50 CancerVar Tier I/II, 1 CIViC A — and those were
    measured rather than assumed, because they are the cost of the change. They are the
    same pathology and not a counterweight to it: the largest is ``PTEN`` at 1000G 1.0
    called ``HP Pathogenic`` by RENOVO while ClinVar calls it ``Benign``, and all 50
    CancerVar rows are one ``KMT2C`` variant at gnomAD 0.494. This is why the exemption is
    still ClinVar-only — widening it to the pipeline's whole rescue mask would re-admit
    exactly these, which is the trade #37 already refused.

    **The absent-versus-empty asymmetry dissolved rather than being fixed.** Under the old
    rule a column present-but-empty could rescue a row while a column missing entirely
    could not, because the missing one is not in ``present`` to contribute a blank — an
    asymmetry that fell out of the comprehension below and that nothing ever argued for. It
    decided ``HMGCR``/``NOS3`` one way and ``FUT2`` the other. Now neither contributes, so
    the two cases agree without a second mechanism to keep them agreeing.

    Neutrality at 1.0 survives, which is what the composition in :func:`apply_filters`
    rests on: every populated value is ``<= 1.0`` and every all-blank row is carried by the
    fallback, so the mask is still all-True by algebra rather than by the caller's guard.
    """
    present = [column for column in FREQUENCY_COLUMNS if column in maf.columns]
    if not present:
        return None

    acceptable = pd.Series(False, index=maf.index)
    measured = pd.Series(False, index=maf.index)
    for column in present:
        values = pd.to_numeric(maf[column], errors="coerce")
        acceptable |= (values <= threshold) & values.notna()
        measured |= values.notna()
    return acceptable | ~measured


def pathogenic_exemption(maf: pd.DataFrame) -> pd.Series:
    """Rows the frequency layer must not drop: ClinVar calls them pathogenic.

    Public for the same reason :func:`frequency_mask` is, and only that reason: the two are
    the halves of the composition ``apply_filters`` performs, and the app-owned suite states
    the composition as the equation #37 specifies rather than re-deriving either half. A
    private name would push the test into restating what it is checking.

    A common variant is not necessarily a harmless one. The reference's witness is
    ``SERPINA1`` chr14:94847262, ``Pathogenic/Pathogenic,_low_penetrance|other`` at a
    population frequency of 2.3% — the only ClinVar-pathogenic call in 181,566 rows above
    0.01, and precisely the kind of low-penetrance allele a frequency cut-off is blind to.
    Without this carve-out the layer discards it for being common.

    The test is the **vendored** one, over the **vendored** constant: the same
    ``has_clinvar_term`` and the same :data:`~vendor.pipeline_utils.CLINVAR_PATHO` that the
    pipeline's own pathogenic retention uses. Nothing here restates the list — a fourth
    drifted copy of a shared list is the failure mode the vendoring exists to prevent, and
    ``tests/test_filter_app_extras.py`` fails if a literal collection in this module ever
    names one of its terms.

    Narrower than the pipeline's rescue mask, which is the point rather than an oversight.
    That mask also admits CancerVar Tier I/II, InterVar pathogenic and CIViC A/B calls, and
    the rows actually in dispute here are the ones surviving on a *potential*-significance
    CancerVar tier while ClinVar calls them ``Benign`` or ``drug_response`` — four of the
    five disputed variants, at frequencies up to 0.995. Exempting the whole rescue mask
    would keep exactly those.

    Independent of ``skip_pathogenic``, which is not an oversight. That flag empties the
    *pipeline's* rescue mask — a statement about which variants the pipeline reports — while
    this exemption is a bound on what the **app's own** filter is allowed to throw away. A
    user who has turned pathogenic rescue off is asking for a stricter report, not for the
    frequency cut-off to start deciding clinical significance, so a ClinVar-pathogenic
    variant that reaches the report on its own criteria keeps its exemption. #37 specifies
    the composition unconditionally.

    Reads :data:`CLINVAR_SIGNIFICANCE` without a presence guard, and the caller is what makes
    that safe: it passes the frame the fills were applied to, where the column is present
    whether or not the MAF carried it. A guard here would be a second, quieter answer to a
    question :mod:`filters.absent_columns` already answers loudly.
    """
    return maf[CLINVAR_SIGNIFICANCE].apply(
        lambda value: has_clinvar_term(value, CLINVAR_PATHO)
    )


def _frequency_note(threshold: float, removed: int, exempted: int) -> Note:
    """What to tell the user the layer did, including what it spared.

    ``INFO`` since issue #151, and the note that made the case for a third level existing at
    all. Nothing here is wrong: the filter the user chose read the columns it needed and
    removed what they asked it to remove. It was nevertheless a yellow warning box, and
    measurably the *only* box the slot drew on either reference MAF on its own arm under all
    four shipped presets — so the app's ordinary output was a warning about success, and yellow
    had nothing left to mean by the time a user met a real one.

    Both counts, always — the exemption clause is not decoration. A user who reads
    "5 variants removed" on a report where six passing rows are above the threshold has been
    told something that looks wrong, and the sixth row is missing from the count for a
    clinical reason they need to know. Stating the rule with a count of zero is the honest
    reading of "nothing was common *and* pathogenic here"; omitting it would leave the
    arithmetic unexplained on exactly the reports where it matters.
    """
    return Note.info(
        f"🌍 Population frequency filter: {removed} variant(s) above {threshold} in every "
        f"available frequency column were removed. ClinVar-pathogenic variants are exempt: "
        f"{exempted} above the threshold were kept."
    )


# ---------------------------------------------------------------------------
# The gene-list adapter
# ---------------------------------------------------------------------------


def _gene_symbols(maf: pd.DataFrame, selection: GeneSelection) -> list:
    """The symbols to restrict to, or none at all if this MAF cannot be restricted.

    ``Hugo_Symbol`` is the one filter input :mod:`filters.absent_columns` does not fill, and
    the reason is that the *app* guards it rather than the pipeline: the vendored body only
    reads the column when a gene file exists, and the app is what writes that file. So the
    neutral response to an absent ``Hugo_Symbol`` is not a value — it is to write no file, at
    which point the vendored gene clause falls back to its own all-``True`` default and the
    restriction is simply not applied.

    Filling the column instead would be the dangerous direction, and this is the one place
    where the choice is visible: any value stood in for a gene symbol matches *no* requested
    gene, so the clause would go all-``False`` and empty the report. Issue #28's principle —
    extra rows are visible, missing rows are not — settles it. The user is told, through
    :func:`_gene_notes`, that the restriction they asked for was not applied.
    """
    if _gene_filter_is_unapplicable(maf, selection):
        return []
    return list(selection.symbols)


def _gene_filter_is_unapplicable(maf: pd.DataFrame, selection: GeneSelection) -> bool:
    """The user asked for a gene restriction and this MAF cannot carry one.

    One predicate rather than the same two-clause test written at both sites, because the two
    sites are the *action* and the *telling* — drop the restriction, and say it was dropped —
    and those going out of step is precisely the silence this ticket exists to remove: a report
    quietly widened with no warning, or a warning about a restriction that was in fact applied.
    """
    return selection.restricts and GENE_SYMBOL not in maf.columns


@contextmanager
def _gene_file(symbols: list) -> Iterator[str]:
    """A path the vendored gene clause will read, for as long as the call needs it.

    The vendored body takes a *file* — it compares the argument against ``"null"``, then
    calls ``os.path.exists``, then ``pd.read_csv(header=None)`` — so a list of symbols
    has to become a one-per-line file. Written here rather than adapted inside the
    vendored code, which stays byte-identical to ``bin/``.

    A context manager, and one that must wrap the vendored *call*: the existence check
    happens inside that call, so closing early makes the gene filter apply to nothing and
    widen the report silently. The file is removed on the way out, including on the error
    path — the naive form leaks one file per re-filter click, on a platform that does not
    reap its temp directory.

    Case-insensitivity comes free from routing: the vendored clause upper-cases both the
    MAF's symbols and the file's. The app used to compare them verbatim, so a lowercase
    paste silently emptied the report.

    Raises:
        TypeError: if handed a bare string rather than a sequence of symbols. A string is
            iterable, so ``"\\n".join("BRCA1")`` is ``"B\\nR\\nC\\nA\\n1"`` — five
            one-character gene tokens, a filter matching nothing, and no error anywhere.
            That collapse is what a single-element gene list degrades into whenever some
            caller joins a list back into a scalar, so it is refused loudly here rather
            than diagnosed later from an empty report.
    """
    if isinstance(symbols, str):
        raise TypeError(
            "_gene_file takes a list of symbols, not a string: a bare "
            f"{symbols!r} would be written one character per line. Tokenise it with "
            "filters.gene_lists.parse_gene_list first."
        )
    if not symbols:
        yield NO_GENE_FILE
        return

    handle = tempfile.NamedTemporaryFile(
        "w", suffix=".txt", prefix="mafigate_genes_", delete=False, encoding="utf-8"
    )
    try:
        handle.write("\n".join(str(symbol) for symbol in symbols) + "\n")
        handle.close()
        yield handle.name
    finally:
        os.unlink(handle.name)


# ---------------------------------------------------------------------------
# Parameter resolution
# ---------------------------------------------------------------------------


def _terms(value) -> list:
    """The list of terms a filter parameter denotes, as the pipeline builds one.

    Named for *terms* rather than for a keep-list because it serves both polarities: five
    of its callers pass a keep-list, and ``filter_variant_classification`` passes an
    exclude list. What it does is the parsing, which is the same either way.

    ``bin/filter_variants.py`` does ``[x.strip() for x in arg.split(",") if x.strip()]``,
    so blanks are dropped and an empty list is a real, expressible value. An empty
    keep-list makes its source contribute nothing to the guideline union — neutral, not
    narrowing — which is exactly what the vendored ``isin([])`` does.

    There is no catch-all sentinel here, and that is the fix for divergence #5. The
    deleted code skipped a whole clause when a source was set to ``All``, which left the
    guideline block unapplied and widened the somatic report 167-fold. The vendored code
    has no such concept: a sentinel value is simply a value nothing matches, so passing
    one through is neutral. Issue #36 removed it from the widgets too, and it is that
    neutrality here which makes dropping it from a stale parameter file safe — the two
    reach this function as the same restriction.
    """
    if value is None:
        return []
    if isinstance(value, str):
        value = [value]
    return [str(item).strip() for item in value if str(item).strip()]


def _vaf_threshold(params: Mapping[str, Any], arm: str) -> float:
    """The one VAF threshold for this arm.

    The contract carries a single ``vaf_threshold`` whose value differs per arm, because
    the app knows its arm and the pipeline does not. The germline-specific spelling is
    read first on the germline arm so that a legacy dict carrying *both* keys — where
    ``vaf_threshold`` holds the somatic value — still gets the germline number. Three
    sources once disagreed here: the config and both presets said 0.2 while a widget
    fallback said 0.3.
    """
    if arm == "germline" and params.get("vaf_threshold_germline") is not None:
        return params["vaf_threshold_germline"]
    return params.get("vaf_threshold", 0.0)


def _gene_selection(params: Mapping[str, Any], arm: str) -> GeneSelection:
    """The gene list to restrict to, tokenised.

    ``filter_genes`` is the contract's one unprefixed key and holds a list. The legacy
    per-arm keys hold a comma-separated string. Both go through the *same* tokeniser, and
    that is the fix: the previous reader split the contract's list on nothing and the
    legacy string on commas only, so a one-symbol-per-line paste — the shape every gene
    file in this project actually has — became a single token matching no row anywhere,
    and the report came back empty with nothing said about why.

    A present-but-empty ``filter_genes`` wins over a legacy key rather than falling
    through to it. That is deliberate: ``[]`` is the parameter page saying "no gene
    filter", and treating it as absent would resurrect a stale preset string and narrow
    the report behind the user's back.
    """
    declared = params.get("filter_genes")
    if declared is None:
        declared = params.get(f"{arm}_genes", "")
    return parse_gene_list(declared)


def _gene_notes(maf: pd.DataFrame, selection: GeneSelection) -> tuple[Note, ...]:
    """What to tell the user about their gene list, including what this frame lacks.

    Two sources, in order. The tokeniser's own messages — tokens dropped, a heading
    dropped, or nothing usable left at all — travel with the report rather than staying on
    the parameter page, because those are different pages and the warning is only useful
    beside the result.

    Then the part that needs the MAF, and so cannot live in the tokeniser: a requested
    symbol this frame does not carry. Worth saying because a gene with no variants and a
    mistyped gene look identical from the report — both simply make it smaller — and only
    one of them is what the user meant.

    An absent ``Hugo_Symbol`` is the case above that one: not "which of your genes are
    missing" but "this MAF cannot be restricted by gene at all". It widens the report rather
    than narrowing it (see :func:`_gene_symbols`), so it is said plainly and not escalated —
    the rows it lets through are on screen to be judged.

    Issue #150 took a sentence off the end of it. The note used to close *"a file without that
    column would normally be refused rather than filtered"* — the same clause
    :func:`~filters.absent_columns._not_a_complete_result` had already deleted at #136, still
    shipping here because #136 was reading that module and this copy was out of its view.
    **MAFigate never refuses on an absent column**, so "normally" could only have meant *in the
    pipeline*, which is the comparison decision 2 of the map retired. What survived the deletion
    is the clause that was doing the work anyway: *wider than the one you asked for* measures
    this report against **the report the user asked for**, which is a comparison the user can
    make, not one about what some other tool would have done with their file.

    This is the one member of the missing-column family whose direction runs the other way —
    #136 settled that a *fill* "can only take rows out", and nothing is filled here — so it is
    the one that says so, and it says it as a property of what MAFigate did.

    ``⚠️`` and a ``WARNING`` since issue #151, which is the ticket #150 raised and declined. Why
    the glyph moved off #150's ``ℹ️`` is told once, in :mod:`filters.notes`, and not retold here.

    What belongs here is why this note's *level* is the one it is, which is not about loudness.
    *We could not apply the restriction you asked for* means the report is not the one the user
    requested — :data:`~filters.notes.WARNING`. The third level #151 added is for an account of a
    run that did what was asked, and this is not one. Issue #28's *extra rows are visible,
    missing rows are not* is still why it is said plainly rather than escalated: the rows it lets
    through are on screen to be judged.
    """
    notes = list(selection.messages())

    if _gene_filter_is_unapplicable(maf, selection):
        notes.append(
            Note.warning(
                f"⚠️ Gene list: this MAF carries no `{GENE_SYMBOL}` column, so the gene "
                f"restriction you asked for ({len(selection.symbols)} symbol(s)) **was not "
                "applied** and this report covers every gene in the file, so it is wider than "
                "the one you asked for."
            )
        )

    if selection.restricts and GENE_SYMBOL in maf.columns:
        absent = missing_symbols(selection.symbols, maf[GENE_SYMBOL])
        if absent:
            named = name_tokens(absent)
            # ``INFO``, and the glyph stays ``🧬`` rather than becoming ``ℹ️``: in the info tier
            # the glyph marks the *topic* — this note's sibling there is the frequency filter's
            # ``🌍`` — because the box already carries the level and a topic is the thing the
            # box cannot say. The restriction the user asked for did apply; a requested symbol
            # simply has no variants in this file, which is a fact about the MAF, not about the
            # report being wrong. Contrast the tokeniser's rejections next door, which warn
            # because there the symbol was never filtered on at all.
            notes.append(
                Note.info(
                    f"🧬 Gene list: {len(absent)} of the {len(selection.symbols)} requested "
                    f"gene(s) not present in this MAF — {named}. Those contribute no rows; "
                    "the filter still applies to the rest."
                )
            )

    return tuple(notes)


def _skip_pathogenic(params: Mapping[str, Any]) -> bool:
    """The pipeline's polarity, from either spelling.

    The flag is only resolved here and then passed down; nothing in this module branches
    on it. Getting the polarity wrong is worth +2 rows on each arm, silently.
    """
    if "skip_pathogenic" in params:
        return bool(params["skip_pathogenic"])
    return not params.get("keep_pathogenic", True)
