"""Which predictor columns hold the chromosome instead of a score, settled once per file.

A predictor column that ANNOVAR did not supply does **not** arrive as ``.``. Both classifiers
resolve their columns positionally, over a dict whose values are initialised to ``0``
(``resources/CancerVar/CancerVar.py:1933``, ``resources/InterVar/Intervar.py:2179``)::

    def search_key_index(line, dict):
        cls = line.split("\\t")
        for key in dict.keys():
            for i in range(1, len(cls)):
                ii = i - 1
                if key == cls[ii]:
                    dict[key] = ii
                    break

If the name is not in the ANNOVAR header the index **stays 0**, the writer emits ``cls[0]``, and
column 0 of a ``_multianno.txt`` is ``Chr``. No error, no ``.`` — the chromosome is written into
the score column, silently. Issue #186 charted this; issue #194 is about what the app does with it.

Why it has to be caught here rather than by ``says_nothing``
------------------------------------------------------------
``config.missing_values.says_nothing`` recognises a blank, a ``.`` and a dash. It cannot recognise
this, because the cell is not blank: it is a plausible-looking number. A chromosome ``22`` in
``CADD_phred`` reads as a CADD phred of 22, and ``components/variant_detail.py`` scores that
against ClinGen's PP3/BP4 calibration and reports evidence from it. That is the on-screen failure
this module exists to stop.

Whole column, not one cell
--------------------------
The rule asks of a *column in a file*, not of a cell: do (almost) all of its non-missing values
equal their own row's chromosome? Measured over 185 real MAFs on this machine (188 byte-distinct
``.maf`` files found under ``~``, 185 carrying a readable ``Hugo_Symbol``/``Chromosome`` header),
the two populations do not overlap and are nowhere near overlapping:

===========================  ============  =====================  ==================  ====================
column                       clean files   max clean fraction     contaminated files  min bad fraction
===========================  ============  =====================  ==================  ====================
``CADD_phred``               16            0.0004                 38                  **1.0000**
``CADD_raw``                 16            0.0004                 38                  **1.0000**
``phastCons20way_mammalian`` 15            0.0345                 35                  **1.0000**
``Polyphen2_HDIV_score``     120           0.1509                 35                  **1.0000**
``SIFT_score``               161           0.0323                 0                   —
``GERP++_RS``                159           0.0017                 0                   —
``REVEL_score``              105           0.0000                 0                   —
``VEST4_score``              105           0.0000                 0                   —
``BayesDel_noAF_score``      146           0.0000                 0                   —
``MetaSVM_score``            161           0.0000                 0                   —
``am_pathogenicity``         132           0.0000                 0                   —
===========================  ============  =====================  ==================  ====================

The last row is issue #204's, measured when issue #203's AlphaMissense section wired that column:
34,279 drawable cells across 180 byte-distinct MAFs, and **not one** equals its own row's
chromosome. It is listed because a clean column measures clean — the rule is that a table scoring
a column asks, not that only dirty columns are named.

Contamination is **all-or-nothing**: in a contaminated file *every* non-missing cell equals the
chromosome — 15,636 of 15,636 for ``CADD_phred``. So a column-level verdict separates the two
populations with a 6.6x margin at :data:`_CONTAMINATED_FRACTION`, and costs one vectorised
comparison per column per file rather than one per rendered cell.

The per-cell alternative was measured and rejected
--------------------------------------------------
Asking the same question of a single cell — *does this value equal this row's chromosome?* — is
equally sensitive and considerably more dangerous. It flags **1,622 real** ``SIFT_score`` cells and
**581 real** ``Polyphen2_HDIV_score`` cells: a legitimate score of exactly ``1.0`` on a ``chr1``
variant. Under issue #194's decision to draw the value with a warning rather than suppress it,
those become alarming warnings over correct numbers. The column-level rule cannot make that
mistake, because one coincidence in a column of thousands moves the fraction by nothing.

The documented-range rule was measured and rejected too
-------------------------------------------------------
Issue #186 offered a second candidate: a value outside the column's documented range means
*absent*. It does not work where it is needed most. ``CADD_phred``'s range is 0–99 and every
chromosome number lies inside it, so the range rule catches **0 of 15,636** contaminated CADD
cells — on the one row the app draws today, it is blind. It does better on the 0–1 columns (8,161
of 9,930 ``phastCons20way_mammalian`` cells) but misses ``chr1`` there as well, since ``1.0`` is
in range. It also has false positives of its own: 6 ``BayesDel_noAF_score`` cells in clean files
sit outside the documented range for reasons that are not contamination.

No arm gate, deliberately
-------------------------
Issue #186 established that provenance differs per arm for the same column name: ``.intervar``
never emits ``Polyphen2_HDIV_score``, so germline receives it *by name* from the ``.grl_p`` merge
and cannot be contaminated, while somatic receives it positionally from ``.cancervar`` and can.
It would therefore be possible to check only the columns that are positional on this file's arm.

This module does not, and the reason is that it **measures** rather than predicts. Arm-gating a
measurement can only subtract coverage, and can only be wrong: a future dbNSFP release or a change
to the ``.grl_p`` merge would move a column from by-name to positional, and an ungated check
notices while a gated one reads ``clean`` forever. On the corpus as it stands the two give
identical answers — which is exactly when the cheaper-to-be-right choice is free.

What this module does not decide
--------------------------------
Whether a *present and absent* column should say that the file's ANNOVAR release predates it. Every
MAF written since ``bin/utils.py:283`` stamps ``## Annovar database_names …`` into its header, but
``vendor.pipeline_utils.read_maf`` counts those lines only to pass ``header=c`` to ``read_csv``, so
the stamp reaches no frame and nothing in ``streamlit_app/`` reads it. Issue #194 ruled that out of
its own scope: for *this* case the measurement already says everything the warning needs, and #186
measured contamination as confined to un-stamped and ``dbnsfp42c`` files anyway, so the stamp would
be redundant with the frame. It stays the map's open "annotation vintage" question.
"""

from __future__ import annotations

import pandas as pd

from config.columns import spelled_in
from config.missing_values import says_nothing, says_nothing_over

#: The column every check compares against. Spelled here rather than imported from
#: ``config.chromosome_spelling`` for the reason that module gives for not offering its own
#: constant around: this module reads the column, it does not decide its spelling.
CHROMOSOME = "Chromosome"

#: Where ``page_modules.data_loading`` parks this file's verdict, for
#: ``components.variant_detail`` to read on every render of the panel.
#:
#: Named here rather than in either of those two modules because both touch it and neither owns
#: the question. The panel importing the key from the load page would say a component depends on
#: a page, which ``tests/test_component_seams.py`` is about not being true.
#:
#: **Assigned, never popped.** Unlike the load page's banners, this is not a message waiting to
#: be drawn once; it is a standing fact about the file on screen, read again every time a user
#: opens a variant. It is replaced wholesale in the same tail that sets ``maf_data``, so it
#: cannot outlive the frame it describes.
SESSION_KEY = "contaminated_predictor_columns"

#: This module is deliberately free of Streamlit — ``config/`` is imported by the two
#: streamlit-free CI jobs — so the session read lives with the panel that needs it. What is
#: shared is the key and the empty answer, so a reader with no verdict stashed and a file with
#: no contamination take the same path.
NOTHING_CONTAMINATED = frozenset()

#: Above this fraction of a column's non-missing cells equalling their row's chromosome, the
#: column is holding chromosomes rather than scores.
#:
#: Issue #186 published ``>0.30``; this is ``>0.50``. **The two give identical verdicts on all
#: 185 files** — contaminated columns sit at 1.0000 and the highest clean one at 0.1509 — so the
#: difference is margin against a false warning, not a different answer. The higher number was
#: chosen because a warning over a correct score is the costlier of the two mistakes here.
_CONTAMINATED_FRACTION = 0.50

#: The columns asked about. Every predictor column the variant panel draws or is about to draw:
#: the five of ``components.reference_scales.CLINGEN_SVI_TABLE`` that have a wired dbNSFP
#: column, plus the five CancerVar's CBP10 and InterVar's PP3/BP4 read, plus the two legacy
#: conservation columns #186 measured as contaminated.
#:
#: Not gated by arm — see the module docstring. Not derived from
#: ``CLINGEN_SVI_TABLE`` either, because that table is the *ClinGen* scale's five and this
#: question is about any positionally-resolved column, which is a larger set — as issue #190's
#: score-context section then proved by drawing more of them.
#:
#: **The three ``_pred`` columns are asked about even though no file holds a chromosome in one.**
#: Issue #190 measured them for the first time — 0 contaminated of 94 somatic and 57 germline files
#: for ``FATHMM_pred`` and ``MutationAssessor_pred``, 0 of 57 and 54 for ``Polyphen2_HDIV_pred``,
#: with a maximum clean fraction of exactly 0.0 — because a file missing one of them turns out to
#: drop the column rather than fill it positionally. They are here anyway for the reason the module
#: docstring gives for not arm-gating: this is a measurement, and a measurement that is currently
#: negative costs one vectorised comparison and notices if that ever changes. What would otherwise
#: happen is the #194 failure in a new place — ``components.predictor_context`` reporting
#: ``chr7`` as PolyPhen-2's call.
PREDICTOR_COLUMNS = (
    "CADD_phred",
    "CADD_raw",
    "phastCons20way_mammalian",
    "phyloP46way_placental",
    "Polyphen2_HDIV_score",
    "Polyphen2_HDIV_pred",
    "REVEL_score",
    "VEST4_score",
    "BayesDel_noAF_score",
    "SIFT_score",
    "GERP++_RS",
    "MetaSVM_score",
    "MetaLR_score",
    "FATHMM_score",
    "FATHMM_pred",
    "MutationAssessor_score",
    "MutationAssessor_pred",
    "dbscSNV_ADA_SCORE",
    "dbscSNV_RF_SCORE",
    # Added with issue #203's AlphaMissense section, which scores this column. Measured before
    # being wired: **clean on every file**, maximum chromosome fraction exactly 0.0000. The check
    # is here anyway, because the rule is that a table scoring a column asks rather than concludes
    # — a predictor column ANNOVAR did not supply arrives holding the chromosome, not a `.`, so a
    # score can look plausible and be nothing.
    "am_pathogenicity",
)


def held_value(row: "pd.Series", column: str) -> "str | None":
    """What a contaminated column holds for this variant, as a table should print it.

    Args:
        row: the variant.
        column: a column the load-time verdict named.

    Returns:
        str: the cell, formatted as a score where it is a number and as the file spells it where it
        is not — ``22`` and ``X`` respectively, which is the difference between an autosome and a
        sex chromosome landing in a score column. ``None`` when the column is not in this row or
        says nothing, where there is no row to draw and drawing nothing is right.

    Separate from a score parser rather than a flag on one, because the two want opposite things
    from an unparseable cell. A parser returns ``None``, correctly: a value it cannot read is one it
    cannot score. This hands ``X`` back, because a value that cannot be read is exactly what needs
    saying.

    It lives here rather than in either renderer because **two** tables now print such a cell —
    the ClinGen PP3/BP4 table (issue #194) and issue #190's score-context tables — and the rule
    for turning a contaminated cell into text is one rule. Streamlit-free like the rest of this
    module, so the two streamlit-free CI jobs still import it.

    ``column`` is the canonical spelling, matching what :func:`contaminated_columns` reports; the
    cell is found under the file's own spelling of it.
    """
    spelling = spelled_in(row.index, column)
    if spelling is None:
        return None
    value = row[spelling]
    if says_nothing(value):
        return None

    held = str(value).strip()
    try:
        number = float(held)
    except (ValueError, TypeError):
        return held
    return f"{number:.3g}" if abs(number) < 100 else f"{number:.1f}"


def holds_the_chromosome(values: pd.Series, chromosomes: pd.Series) -> bool:
    """Whether this column holds chromosomes rather than the scores it is named for.

    Args:
        values: the predictor column, of any dtype.
        chromosomes: the ``Chromosome`` column, index-aligned with ``values``.

    Returns:
        bool: ``True`` when more than :data:`_CONTAMINATED_FRACTION` of the column's non-missing
        cells equal their own row's chromosome. ``False`` when the column says nothing at all,
        which is a missing column's shape and not this module's business.

    Compared as **rendered strings with any ``chr`` prefix stripped**, and then again as floats
    where both sides parse. Two spellings have to agree for one comparison to be right:
    ``normalise_chromosome_spelling`` has already given every chromosome the ``chr`` prefix by the
    time this runs, so the score cell holds ``22`` while the chromosome column holds ``chr22``;
    and a column pandas has typed as float renders the same cell as ``22.0``. Comparing only one
    way misses the contamination this exists to find.
    """
    said_nothing = says_nothing_over(values)
    said_something = ~said_nothing
    n = int(said_something.sum())
    if n == 0:
        return False

    scores = values[said_something].astype(str).str.strip()
    chroms = chromosomes[said_something].astype(str).str.strip()
    bare = chroms.str.replace(r"^chr", "", case=False, regex=True)

    same = scores.str.lower() == bare.str.lower()

    # And numerically, for the float-typed columns: '22.0' is the same cell as chromosome '22',
    # and a string comparison alone calls it a score.
    as_float = pd.to_numeric(scores, errors="coerce")
    chrom_float = pd.to_numeric(bare, errors="coerce")
    both_parse = as_float.notna() & chrom_float.notna()
    same = same | (both_parse & (as_float == chrom_float))

    return bool(int(same.sum()) / n > _CONTAMINATED_FRACTION)


def contaminated_columns(maf: pd.DataFrame) -> frozenset:
    """The predictor columns in this MAF that hold the chromosome instead of a score.

    Args:
        maf: a MAF, read and with its chromosome spelling already settled. Not modified.

    Returns:
        frozenset: the names of the contaminated columns, of :data:`PREDICTOR_COLUMNS`. Empty
        when ``Chromosome`` is absent — the comparison has nothing to be against, and
        ``validate_required_columns`` refuses such a file a moment later anyway.

    Runs once per file, in ``page_modules.data_loading``, *after*
    ``config.chromosome_spelling.normalise_chromosome_spelling``: the answer is a comparison
    against ``Chromosome``, so it must be taken once that column's spelling is a settled value
    rather than one of two.
    """
    if CHROMOSOME not in maf.columns:
        return frozenset()

    chromosomes = maf[CHROMOSOME]
    verdict = set()
    for column in PREDICTOR_COLUMNS:
        # The file's own spelling, because 2 of 167 real MAFs write ``GERP++_RS`` as
        # ``GERP.._RS`` — see :func:`config.columns.spelled_in`. Looked up under the file's
        # spelling and **reported under the canonical one**, so every reader keeps asking this
        # verdict the question it already asks (``"CADD_phred" in untrustworthy``) and only the
        # measurement has to know a column can arrive named twice.
        spelling = spelled_in(maf.columns, column)
        if spelling is None:
            continue
        if holds_the_chromosome(maf[spelling], chromosomes):
            verdict.add(column)
    return frozenset(verdict)
