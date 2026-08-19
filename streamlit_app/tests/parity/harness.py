"""The parity harness: run both filters over one MAF, diff, and attribute the diff.

The measurement instrument issue #17 defines the destination in terms of. Two sides:

* **pipeline** — ``bin/filter_variants.py`` in a fresh subprocess. Fresh is not
  optional: ``keep = KEEP`` at ``filter_variants.py:352`` is an alias, so one
  germline run mutates the module-level list and the next somatic run in the same
  process emits germline columns or raises (issue #15 predicted it, issue #18
  reproduced it). ``psycopg`` is stubbed because ``filter_variants.py:10`` imports
  ``database.database_utils``, which imports it; it is not a dependency of the app.
* **app** — ``filters.variant_filters.apply_filters`` in process. Since issue #33 that
  function makes no filter decision of its own: it calls the same functions the pipeline
  side calls, vendored. So this harness has changed what it measures — from "do two
  implementations agree?" to "does the app's plumbing around one implementation preserve
  its verdict?". The pipeline subprocess stays, because that plumbing includes reading
  parameters, adapting a gene list to a file, and labelling rows, and any of those can
  still lose a verdict.

Both sides are handed **the same DataFrame**, loaded by the vendored
``pipeline_utils.read_maf``. That is deliberate: feeding one frame to both sides is what
makes every divergence this harness reports a divergence in *filter logic* rather than
in I/O.

It began as a workaround. The app used to have a loader of its own that split numeric
columns into strings, which made ``common_filters`` raise ``TypeError: '>=' not
supported between str and int`` before any comparison was possible (issue #10), so the
harness could not have used it. Issue #32 replaced that loader with a delegation to the
vendored reader, and the workaround became an identity —
``test_parity.py::test_harness_loads_the_frame_the_app_would_load`` asserts per fixture
that what this loads is what the app loads.

What is measured, per case:

* PASS/NOPASS set equality on the four-column variant key;
* the output **column set** of ``filtered.*.pass.tsv`` against the app's two
  competing column lists;
* **attribution**: for every variant the two sides disagree on, which of the map's
  known divergences can explain it.

Run standalone for the baseline table::

    python streamlit_app/tests/parity/harness.py --report
    python streamlit_app/tests/parity/harness.py --update-baseline
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
STREAMLIT_APP = HERE.parent.parent
REPO_ROOT = STREAMLIT_APP.parent
BIN_DIR = REPO_ROOT / "bin"
BASELINE_PATH = HERE / "baseline.json"

if str(STREAMLIT_APP) not in sys.path:
    sys.path.insert(0, str(STREAMLIT_APP))

from tests.parity.contract import (  # noqa: E402
    APP_VARIANT_CLASSIFICATIONS,
    CASES,
    CASES_BY_NAME,
    FIXTURE_DIR,
    KEY_COLUMNS,
    Case,
    app_params,
    app_variant_classifications,
    pipeline_args,
    read_gene_tokens,
)
from vendor.pipeline_utils import CLINVAR_PATHO, has_clinvar_term, read_maf  # noqa: E402

PIPELINE_AVAILABLE = BIN_DIR.is_dir()


# ---------------------------------------------------------------------------
# Pipeline side
# ---------------------------------------------------------------------------


@lru_cache(maxsize=1)
def _psycopg_stub_dir() -> str:
    """A directory holding a stub ``psycopg`` package, prepended to PYTHONPATH.

    ``bin/filter_variants.py`` imports ``database.database_utils``, which imports
    psycopg at module scope. Nothing on the filtering path touches a database, so an
    import-satisfying stub is enough; installing psycopg would pull a Postgres client
    into the app's test environment for no benefit.
    """
    stub = Path(tempfile.mkdtemp(prefix="parity_psycopg_stub_"))
    pkg = stub / "psycopg"
    pkg.mkdir()
    (pkg / "__init__.py").write_text(
        '"""Import-satisfying stub. See harness._psycopg_stub_dir."""\n\n\n'
        "class Connection:\n    pass\n\n\n"
        "def connect(*args, **kwargs):\n"
        '    raise RuntimeError("parity harness: psycopg is stubbed")\n'
    )
    (pkg / "rows.py").write_text(
        "def dict_row(*args, **kwargs):\n"
        '    raise RuntimeError("parity harness: psycopg is stubbed")\n'
    )
    return str(stub)


@dataclass
class PipelineResult:
    verdicts: dict[tuple, str] = field(default_factory=dict)
    output_columns: list[str] = field(default_factory=list)
    error: str | None = None
    stdout: str = ""

    @property
    def passed(self) -> set[tuple]:
        return {k for k, v in self.verdicts.items() if v == "PASS"}


def run_pipeline(case: Case) -> PipelineResult:
    """Run ``bin/filter_variants.py`` over the case's fixture in a fresh process."""
    work = Path(tempfile.mkdtemp(prefix="parity_pipeline_"))
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join(
        [_psycopg_stub_dir(), str(BIN_DIR), env.get("PYTHONPATH", "")]
    ).rstrip(os.pathsep)

    cmd = [sys.executable, str(BIN_DIR / "filter_variants.py")]
    cmd += pipeline_args(case)
    cmd += ["--output", "out.maf"]
    proc = subprocess.run(cmd, cwd=work, env=env, capture_output=True, text=True)

    if proc.returncode != 0:
        return PipelineResult(error=_last_exception_line(proc.stderr), stdout=proc.stdout)

    out = read_maf(work / "out.maf")
    verdicts = dict(zip(variant_keys(out), out["filter"]))

    # Column-set parity is about the filtered outputs, not out.maf. Both branches
    # end in the same `out[keep]` reselection, so either file answers the question;
    # a branch is absent only when that verdict has no rows.
    columns: list[str] = []
    for name in ("filtered.out.maf.pass.tsv", "filtered.out.maf.nopass.tsv"):
        path = work / name
        if path.exists():
            columns = list(pd.read_csv(path, sep="\t", nrows=0, low_memory=False).columns)
            break

    return PipelineResult(verdicts=verdicts, output_columns=columns, stdout=proc.stdout)


def _last_exception_line(stderr: str) -> str:
    """The ``TypeError: ...`` line from a traceback, or the whole thing if unclear."""
    lines = [line for line in stderr.strip().splitlines() if line.strip()]
    return lines[-1] if lines else "(no stderr)"


# ---------------------------------------------------------------------------
# App side
# ---------------------------------------------------------------------------


@dataclass
class AppResult:
    verdicts: dict[tuple, str] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)
    error: str | None = None
    #: The four-cell decomposition the entry point reports, or ``None`` if it raised.
    #: Recorded so the harness can check the app's own account of *why* each row got its
    #: verdict against the PASS set it actually produced.
    cells: dict[str, int] | None = None
    #: The columns the app **refused** the file over (issue #38), or ``None`` if it did not
    #: refuse. A refusal is not a crash: it is the pipeline's own non-verdict, delivered as
    #: a typed exception naming what stopped the file. Distinguished from any other
    #: exception because "the app refused where the pipeline raised" is parity, while "the
    #: app blew up for some unrelated reason" is not, and one string field cannot say which.
    refused_columns: list[str] | None = None
    #: The filter-input columns the app **filled** because the MAF did not carry them (issue
    #: #39). Distinct from :attr:`refused_columns` in every respect that matters here: a
    #: refusal is the pipeline's own non-verdict and therefore parity, while a fill is a
    #: verdict the pipeline never reaches and therefore off parity *by construction*. Recorded
    #: rather than inferred from the warnings, so the harness never has to read prose.
    filled_columns: list[str] = field(default_factory=list)

    @property
    def passed(self) -> set[tuple]:
        return {k for k, v in self.verdicts.items() if v == "PASS"}


def run_app(case: Case, frame: pd.DataFrame | None = None) -> AppResult:
    """Run the app's filter path over the case's fixture.

    One call, one entry point. The app returns every row labelled with its verdict, so the
    harness reads the verdict column rather than diffing two frames — which is also what
    lets it record the diagnostics decomposition alongside.
    """
    from filters.numeric_columns import UnreadableNumericColumns
    from filters.variant_filters import MAFIGATE_FILTER, apply_filters

    frame = load_frame(case) if frame is None else frame
    try:
        labelled, diagnostics = apply_filters(frame, app_params(case))
    except UnreadableNumericColumns as refusal:
        # Recorded by *type and columns*, not by message. The prose is asserted in
        # tests/test_numeric_columns.py, where a wording change belongs; putting it in
        # baseline.json would make every rephrasing look like a parity movement.
        return AppResult(
            error=type(refusal).__name__,
            refused_columns=sorted(refusal.offenders),
        )
    except Exception as exc:  # noqa: BLE001 - the error itself is the measurement
        return AppResult(error=f"{type(exc).__name__}: {exc}")

    verdicts = dict(zip(variant_keys(labelled), labelled[MAFIGATE_FILTER]))
    return AppResult(
        verdicts=verdicts,
        # The prose only. ``AppResult.warnings`` exists to be diffed against a recorded
        # baseline, and the level a note is drawn at is a UI decision that has never been part
        # of the parity question — recording it would make a re-levelling look like a movement.
        warnings=[note.text for note in diagnostics.notes],
        cells=diagnostics.cells(),
        filled_columns=list(diagnostics.filled_columns),
    )


# ---------------------------------------------------------------------------
# The comparison key
# ---------------------------------------------------------------------------


def variant_keys(frame: pd.DataFrame) -> list[tuple]:
    """The four-column variant identity, as hashable tuples of strings.

    Issue #18 measured this key unique within any single reference MAF (0 duplicates
    across all 100 files) and the fixture generator deduplicates on it while pooling,
    so it is unique on every committed fixture — asserted by
    ``test_parity.py::test_comparison_key_is_unique``. It cannot be extended with
    ``Tumor_Sample_Barcode``: that column is ``__UNKNOWN__`` in all 92,547 germline
    reference rows.
    """
    missing = [c for c in KEY_COLUMNS if c not in frame.columns]
    if missing:
        raise KeyError(f"comparison key columns absent: {missing}")
    return list(zip(*(frame[c].astype(str) for c in KEY_COLUMNS)))


@lru_cache(maxsize=None)
def _read_maf_cached(path: str) -> pd.DataFrame:
    return read_maf(path)


def load_frame(case: Case) -> pd.DataFrame:
    """The fixture, loaded once by the vendored reader and shared by both sides.

    Cached on the path rather than the ``Case`` (which is unhashable: contract
    overrides are a dict) and copied out, because the app's ``apply_filters*`` adds a
    ``CancerVar`` column in place when a CancerVar-Evidence column is present.
    """
    return _read_maf_cached(str(case.maf_path)).copy()


# ---------------------------------------------------------------------------
# Attribution
# ---------------------------------------------------------------------------
#
# A count of diverging rows tells you parity is not reached. It does not tell you
# which of the twelve known divergences is responsible, which is what makes a
# baseline actionable. Each probe below is a predicate over one row that is true when
# that divergence *could* explain the row's disagreement. Probes are not mutually
# exclusive — a row can witness several, and the reference subsample was built so
# that they overlap (issue #18) — so the counts do not sum to the diff size.


def _numeric(row, column):
    if column not in row.index:
        return None
    return pd.to_numeric(pd.Series([row[column]]), errors="coerce").iloc[0]


def _probe_vc_semantics(row, c, case):
    """Divergence #1: the pipeline *excludes* the list, the app *includes* it.

    Stated as the decision itself rather than as list membership, because the two
    modes need opposite membership tests: under ``literal`` the app keeps exactly what
    the pipeline drops (so the decisions differ for every row), while under
    ``complement`` they differ only where the value is outside the app's vocabulary.
    """
    vc = row.get("Variant_Classification")
    pipeline_keeps = vc not in c["filter_variant_classification"]
    app_keeps = vc in app_variant_classifications(case)
    return pipeline_keeps != app_keeps


def _probe_vc_outside_app_vocabulary(row, c, case):
    """Divergence #1's floor: a value the app's hardcoded vocabulary cannot express.

    ``parameter_config`` runs before any MAF is loaded, so the vocabulary can never be
    data-driven; an include-list silently drops these rows however long the list gets.
    Two of the values in the reference are minted by the pipeline itself at
    ``bin/utils.py:628-630``. Issue #12.
    """
    vc = row.get("Variant_Classification")
    return pd.notna(vc) and vc not in APP_VARIANT_CLASSIFICATIONS


def _probe_depth_dp_vs_sum(row, c, case):
    """Divergence #2: the pipeline sums t_alt+t_ref, the app reads DP."""
    dp = _numeric(row, "DP")
    alt, ref = _numeric(row, "t_alt_count"), _numeric(row, "t_ref_count")
    if pd.isna(dp) or pd.isna(alt) or pd.isna(ref):
        return False
    return (dp >= c["min_depth"]) != ((alt + ref) >= c["min_depth"])


def _probe_depth_nan(row, c, case):
    """Divergence #2: the app keeps a missing depth, the pipeline drops it."""
    return any(
        pd.isna(_numeric(row, col)) for col in ("DP", "t_alt_count", "t_ref_count")
    )


def _probe_vaf_boundary(row, c, case):
    """Divergence #3: `tumor_f > vaf` is strict; the app used `>=`."""
    threshold = (
        c["vaf_threshold"] if case.arm == "somatic" else c["vaf_threshold_germline"]
    )
    tf = _numeric(row, "tumor_f")
    return pd.notna(tf) and tf == threshold


def _probe_vaf_nan(row, c, case):
    """Divergence #3: the app keeps a missing VAF, the pipeline drops it."""
    return pd.isna(_numeric(row, "tumor_f"))


def _probe_clnsig_split(row, c, case):
    """Divergence #4: has_clinvar_term splits on [|/;,]; the app matches exactly."""
    value = row.get("ClinVar_VCF_CLNSIG")
    keep = c["filter_clinvar"]
    return has_clinvar_term(value, keep) != (value in keep)


def _probe_clnsig_patho_split(row, c, case):
    """Divergence #4 inside the pathogenic rescue, against CLINVAR_PATHO."""
    value = row.get("ClinVar_VCF_CLNSIG")
    return has_clinvar_term(value, CLINVAR_PATHO) != (value in CLINVAR_PATHO)


def _probe_escat_dead_clause(row, c, case):
    """Divergence #7: the app rescues ESCAT in [I, II, III] — values that never occur.

    Confirmed dead on 182K reference rows (real values are IA, IB, IC, ...). The
    synthetic fixtures construct them so the clause is observable rather than merely
    argued.
    """
    return case.arm == "somatic" and row.get("ESCAT") in ("I", "II", "III")


def _probe_germline_escat_unmirrored(row, c, case):
    """Divergence #6: the app ORs an ESCAT clause into the *germline* guidelines.

    ``variant_filters.py:414-426`` adds ``ESCAT.isin(filter_escat)`` to the germline
    clinical OR. The pipeline's ``germline_filters`` guideline block is
    ``InterVar | ClinVar | RENOVO`` — there is no ESCAT term. So a germline row that
    is Benign on all three pipeline sources but carries an in-keep ESCAT grade passes
    the app and fails the pipeline.

    Written by issue #25 as a *candidate* probe in ``reference.py``, because the
    committed fixtures held no row that witnessed it: the one #6-shaped row in the
    set had ``DP`` 46 against a ``min_depth`` of 50, so divergence #2 silenced it and
    both sides answered NOPASS. Issue #27 added the ``germline_escat_only`` cell,
    which pins the row's whole path to the verdict, and moved the probe here — where
    an empty ``unattributed`` bucket now means the probes explain the diff rather
    than that no row reaches them.
    """
    if case.arm != "germline":
        return False
    return row.get("ESCAT") in c["filter_escat"]


def _probe_somatic_patho_no_clinvar(row, c, case):
    """Divergence #7's other half: the app's somatic rescue has no ClinVar clause at all.

    The pipeline rescues on ``CancerVar | ClinVar | CIViC A/B``; the app
    (``variant_filters.py:287-297``) rescues on ``CancerVar | ESCAT | CIViC A/B``. The
    ESCAT clause is dead, so the net effect is that a ClinVar-pathogenic call whose
    CancerVar is Tier III/IV is rescued by the pipeline and dropped by the app. This
    is what the reference's ``Tier_III_Uncertain`` + ``Pathogenic/Likely_pathogenic``
    rows witness, and it is invisible to the split-vs-exact probe because those values
    match ``CLINVAR_PATHO`` exactly.
    """
    if case.arm != "somatic" or c["skip_pathogenic"]:
        return False
    clinvar_patho = has_clinvar_term(row.get("ClinVar_VCF_CLNSIG"), CLINVAR_PATHO)
    cancervar_patho = row.get("CancerVar") in ("Tier_II_potential", "Tier_I_strong")
    civic_patho = _civic_matches(row, ["A", "B"])
    return clinvar_patho and not cancervar_patho and not civic_patho


def _probe_germline_patho_renovo(row, c, case):
    """Divergence #8: the app's germline rescue uses RENOVO where the pipeline uses ClinVar.

    Both rescue on ``InterVar``, so a row only diverges here when InterVar does *not*
    already rescue it and exactly one of RENOVO / ClinVar does.
    """
    if case.arm != "germline" or c["skip_pathogenic"]:
        return False
    renovo_patho = row.get("RENOVO_Class") in (
        "LP Pathogenic",
        "IP Pathogenic",
        "HP Pathogenic",
    )
    intervar_patho = row.get("InterVar") in ("Pathogenic", "Likely pathogenic")
    clinvar_patho = has_clinvar_term(row.get("ClinVar_VCF_CLNSIG"), CLINVAR_PATHO)
    return (renovo_patho != clinvar_patho) and not intervar_patho


def _civic_matches(row, keep) -> bool:
    """Substring CIViC match — the test *both* sides use. See PROBES' note on #12."""
    value = row.get("CIViC_Evidence_Level")
    if value is None or pd.isna(value):
        return False
    return any(elem in str(value) for elem in keep)


def _probe_gene_case(row, c, case):
    """Divergence #9: the pipeline upper-cases both sides; the app compares verbatim."""
    if c["genes"] is None:
        return False
    tokens = read_gene_tokens(case.genes_path)
    symbol = row.get("Hugo_Symbol")
    if not isinstance(symbol, str):
        return False
    upper_match = symbol.upper() in {t.upper() for t in tokens}
    exact_match = symbol in tokens
    return upper_match != exact_match


def _probe_skip_civic_unmirrored(row, c, case):
    """The app *had* no ``skip_civic``, so it filtered on CIViC when the pipeline did not.

    Closed by issue #33: the flag is now a contract parameter and is passed to the
    vendored ``somatic_filters``, which takes it. Retained as a regression probe — if the
    app ever stops forwarding it, ``civic_skipped`` diverges and this names why rather
    than leaving the rows unattributed.
    """
    return bool(c["skip_civic"]) and "CIViC_Evidence_Level" in row.index


#: **Every probe below describes a divergence issue #33 closed.** They are retained, and
#: the line numbers in their docstrings point into an implementation that no longer
#: exists, because a probe is now a *regression* detector rather than a measurement: it
#: names which known divergence has come back rather than leaving a row in the
#: ``unattributed`` bucket. The consequence for the suite is that no probe has a witness
#: any more — attribution is empty for every case — and ``test_baseline_integrity.py``
#: pins that condition, rather than reporting fourteen unwitnessed probes as a fault.
#:
#: **A probe is no longer the instrument that measures fixture coverage**, and issue #242
#: is where that changed. Requiring each probe to explain a diverging row was that
#: instrument, and it inverted when the diff went empty: the assertion relaxed into
#: asserting the *absence* of attribution and passed on any fixture set reaching row
#: parity. ``mutation_coverage.py`` measures coverage by re-injecting each divergence into
#: the app instead, which needs no divergence to already exist.
#:
#: Note on divergence #12. The map records the app as matching CIViC with an exact
#: ``.isin(["A", "B"])`` that "cannot match list-repr at all". Measured, that is not
#: what the code does: both the guideline clause (``variant_filters.py:238``) and the
#: pathogenic clause (``:295``) call the app's own ``has_element_from_list``, which is
#: a substring test behaviourally identical to the pipeline's. There is therefore no
#: #12 probe — the divergence does not exist, which the ``civic_present`` case
#: confirms end-to-end and ``test_civic_matching_already_agrees`` pins directly.
PROBES = {
    "#1 vc_exclude_vs_include": _probe_vc_semantics,
    "#1 vc_outside_app_vocabulary": _probe_vc_outside_app_vocabulary,
    "#2 depth_dp_vs_sum": _probe_depth_dp_vs_sum,
    "#2 depth_nan_kept_by_app": _probe_depth_nan,
    "#3 vaf_at_threshold": _probe_vaf_boundary,
    "#3 vaf_nan_kept_by_app": _probe_vaf_nan,
    "#4 clnsig_split_vs_exact": _probe_clnsig_split,
    "#4 clnsig_patho_split_vs_exact": _probe_clnsig_patho_split,
    "#6 germline_escat_unmirrored": _probe_germline_escat_unmirrored,
    "#7 escat_dead_retain_clause": _probe_escat_dead_clause,
    "#7 somatic_patho_no_clinvar_clause": _probe_somatic_patho_no_clinvar,
    "#8 germline_patho_renovo_vs_clinvar": _probe_germline_patho_renovo,
    "#9 gene_case_sensitivity": _probe_gene_case,
    "skip_civic_unmirrored": _probe_skip_civic_unmirrored,
}


def attribute(case: Case, frame: pd.DataFrame, keys: set[tuple]) -> dict[str, int]:
    """How many of ``keys`` each known divergence could explain."""
    if not keys:
        return {}
    c = case.contract
    all_keys = variant_keys(frame)
    positions = [i for i, k in enumerate(all_keys) if k in keys]
    subset = frame.iloc[positions]

    counts: dict[str, int] = {}
    unexplained = set()
    for pos, (_, row) in zip(positions, subset.iterrows()):
        hit = False
        for name, probe in PROBES.items():
            if probe(row, c, case):
                counts[name] = counts.get(name, 0) + 1
                hit = True
        if not hit:
            unexplained.add(all_keys[pos])
    if unexplained:
        counts["unattributed"] = len(unexplained)
    return dict(sorted(counts.items()))


# ---------------------------------------------------------------------------
# The comparison
# ---------------------------------------------------------------------------


@dataclass
class ParityResult:
    case: str
    arm: str
    fixture: str
    rows: int
    pipeline_pass: int | None = None
    app_pass: int | None = None
    only_pipeline: int = 0
    only_app: int = 0
    pipeline_error: str | None = None
    app_error: str | None = None
    pipeline_columns: list[str] = field(default_factory=list)
    app_columns: list[str] = field(default_factory=list)
    attribution: dict[str, int] = field(default_factory=dict)
    only_pipeline_keys: list[tuple] = field(default_factory=list)
    only_app_keys: list[tuple] = field(default_factory=list)
    #: The app's own four-cell account of how it reached the verdict. Deliberately kept
    #: out of ``to_json`` — the cells are asserted as a *partition* against the PASS set
    #: (``test_diagnostics_decomposition_reproduces_the_verdict``) rather than recorded as
    #: four more numbers that would become a second baseline to maintain.
    app_cells: dict[str, int] | None = None
    #: Which columns the app refused the file over, or ``None`` if it did not refuse.
    #: See :attr:`AppResult.refused_columns`.
    app_refused_columns: list[str] | None = None
    #: Which filter-input columns the app filled. Like :attr:`app_cells`, deliberately kept out
    #: of ``to_json``: no committed case has an absent filter input, so recording it would add
    #: an always-empty field to every baseline record and force a re-recording to say nothing.
    #: It is read by :attr:`off_parity_by_construction` and by
    #: ``tests/parity/test_absent_columns.py``, which is where the deviation is asserted.
    app_filled_columns: list[str] = field(default_factory=list)

    @property
    def in_parity(self) -> bool:
        return (
            self.errors_in_parity and self.only_pipeline == 0 and self.only_app == 0
        )

    @property
    def off_parity_by_construction(self) -> bool:
        """The app reported on a file the pipeline refused outright, because it filled.

        Not a divergence to be closed, and not parity either — a third state, and the whole
        reason it has a name. Issue #39 chose to keep the app usable on an incomplete MAF, so
        the app reaches a verdict where the pipeline raises ``KeyError``. There is no pipeline
        answer for that verdict to be right or wrong about, which means :attr:`in_parity` is
        false and must stay false: a case in this state has to be asserted against its own
        recorded expectations, never counted towards parity.

        Kept separate from :attr:`errors_in_parity` so the ratchet cannot be loosened by it. If
        a filled run were ever folded into "the errors agree", then *any* app verdict against a
        pipeline ``KeyError`` would read as parity — including a fill the app performed by
        mistake on a file it should have refused.
        """
        return bool(self.app_filled_columns) and self.pipeline_error is not None

    @property
    def errors_in_parity(self) -> bool:
        """Whether the two sides agree about *whether* this file can be filtered.

        Three shapes count as agreement, and the third is issue #38's:

        * neither side failed;
        * both failed identically — which is how this read before the app had any failure
          of its own;
        * the app **refused** and the pipeline raised ``TypeError``. Not the same exception,
          deliberately not: the pipeline's answer to an unreadable MAF is a stack trace, and
          the app's is a typed refusal naming every offending column and value. Those are the
          *same non-verdict*, so calling them a divergence would mean the only way to reach
          parity on ``somatic_dot_numeric.maf`` was to reproduce a traceback — which is
          exactly the answer the ticket set out to improve on.

        Two things it deliberately does **not** accept, and both matter:

        * any *other* app-side exception against a pipeline raise. "Both sides failed" is not
          parity on its own; a refusal is parity because it is the pipeline's own verdict,
          and an unrelated crash is not.
        * a refusal against a pipeline error that is **not** a ``TypeError`` — in practice
          the ``KeyError`` an absent column raises. The refusal explicitly *skips* absent
          columns, so a refusal standing in for a ``KeyError`` would mean the app stopped the
          file for a reason that has nothing to do with why the pipeline stopped it.
          Accepting it would let a real divergence read as parity, which is the one failure
          this whole harness exists to prevent.

        Since issue #39 the app usually does not fail at all on an absent column: it fills and
        reports. That is not a third kind of agreement, and this property does not treat it as
        one — it is :attr:`off_parity_by_construction`, a state with no pipeline answer to
        agree with.
        """
        if self.pipeline_error == self.app_error:
            return True
        if self.app_refused_columns is None or self.pipeline_error is None:
            return False
        return self.pipeline_error.startswith("TypeError")

    @property
    def expected_prefix(self) -> list:
        """The pipeline's columns the app undertakes to show, in the pipeline's order.

        The pipeline's own list, less ``PIPELINE_COLUMNS_THE_APP_REPLACES`` — the columns
        the app answers itself and deliberately keeps out of the default view (issue
        #117). Read from the app's constant rather than written out again here: this file
        checks that the app *keeps* its undertaking, and a second copy of the list would
        mean the harness could agree with a resolver that had quietly dropped something
        else.

        Everything the pipeline emits and this list does not name must still appear, in
        order. The subtraction is a named exception, so it makes the check narrower by
        exactly those names and by nothing else.
        """
        from config.columns import PIPELINE_COLUMNS_THE_APP_REPLACES

        return [c for c in self.pipeline_columns if c not in PIPELINE_COLUMNS_THE_APP_REPLACES]

    @property
    def pipeline_columns_not_shown(self) -> list:
        """The pipeline columns this case's output carries and the app does not show.

        Recorded by name in the baseline beside the extras, and for the same reason: the
        subtraction is as much a claim as the addition, and a count would let one name be
        swapped for another silently. **This is what makes an addition to that list a
        deliberate act**: it lands in the committed baseline, so it fails
        ``test_baseline_matches`` until someone regenerates it on purpose.

        Read straight off the constant rather than as ``pipeline_columns`` minus
        :attr:`expected_prefix` — the same answer by a double complement, routed through
        two derived properties to say something simple.
        """
        from config.columns import PIPELINE_COLUMNS_THE_APP_REPLACES

        return [c for c in self.pipeline_columns if c in PIPELINE_COLUMNS_THE_APP_REPLACES]

    @property
    def columns_in_parity(self) -> bool:
        """The app's list opens with the pipeline's, element for element.

        A **prefix**, not equality: the app is a deliberate superset, carrying its own
        derived columns after the pipeline's. And not a length or set comparison —
        before issue #35 the two lists measured 40 against 40 while differing by a
        substitution, so either of those would have called that parity.

        Since issue #117 the prefix compared against is :attr:`expected_prefix` rather
        than ``pipeline_columns``. That is a weaker claim by exactly the named columns,
        which is why they are recorded in the baseline: this property says the app shows
        what the pipeline emits *in the pipeline's order*, and the baseline says which
        names were let off.
        """
        expected = self.expected_prefix
        return self.app_columns[: len(expected)] == expected

    def to_json(self) -> dict:
        """The committed baseline record: counts and attribution, not row contents."""
        return {
            "arm": self.arm,
            "fixture": self.fixture,
            "rows": self.rows,
            "pipeline_pass": self.pipeline_pass,
            "app_pass": self.app_pass,
            "only_pipeline": self.only_pipeline,
            "only_app": self.only_app,
            "pipeline_error": self.pipeline_error,
            "app_error": self.app_error,
            "app_refused_columns": self.app_refused_columns,
            "pipeline_column_count": len(self.pipeline_columns),
            "app_column_count": len(self.app_columns),
            # The app's extras — the tail past the pipeline's prefix. Recorded by name,
            # not by count: this is the one place the superset is visible in the
            # baseline, and a count would hide a substitution in it.
            #
            # Measured past `expected_prefix`, not past `pipeline_columns`: since issue
            # #117 those differ, and slicing at the wrong one would report a pipeline
            # column as an app extra.
            "app_extra_columns": self.app_columns[len(self.expected_prefix) :],
            # And the other direction, for the same reason — see
            # `pipeline_columns_not_shown`.
            "pipeline_columns_not_shown": self.pipeline_columns_not_shown,
            "columns_in_parity": self.columns_in_parity,
            "in_parity": self.in_parity,
            "attribution": self.attribution,
        }


def app_columns(case: Case, frame: pd.DataFrame) -> list[str]:
    """The columns the app shows and exports, by calling the app's own resolver.

    One call, because there is now one answer — this used to compute two competing
    lists; ``config/columns.py`` records why.

    This *imports and calls* the resolver. The grid's list used to be read back out of
    the grid's own module with ``ast`` to avoid importing streamlit; the resolver lives
    in a streamlit-free module precisely so that workaround could go.
    """
    from config.columns import resolve_visible_columns

    available = list(frame.columns) + ["project_id", "filter"]
    return resolve_visible_columns(
        sample_type=case.arm,
        skip_civic=case.contract["skip_civic"],
        available_columns=available,
    )


def compare(case: Case) -> ParityResult:
    """Run both sides over one case and diff them."""
    frame = load_frame(case)
    result = ParityResult(
        case=case.name, arm=case.arm, fixture=case.fixture, rows=len(frame)
    )

    pipeline = run_pipeline(case)
    app = run_app(case, frame)
    result.pipeline_error = pipeline.error
    result.app_error = app.error
    result.app_refused_columns = app.refused_columns
    result.app_filled_columns = app.filled_columns

    if pipeline.error is None:
        result.pipeline_pass = len(pipeline.passed)
        result.pipeline_columns = pipeline.output_columns
    if app.error is None:
        result.app_pass = len(app.passed)
        result.app_cells = app.cells
        result.app_columns = app_columns(case, frame)

    if pipeline.error is None and app.error is None:
        only_pipeline = pipeline.passed - app.passed
        only_app = app.passed - pipeline.passed
        result.only_pipeline = len(only_pipeline)
        result.only_app = len(only_app)
        result.only_pipeline_keys = sorted(only_pipeline)
        result.only_app_keys = sorted(only_app)
        result.attribution = attribute(case, frame, only_pipeline | only_app)

    return result


def compare_all(cases=None) -> list[ParityResult]:
    return [compare(case) for case in (CASES if cases is None else cases)]


# ---------------------------------------------------------------------------
# Baseline + reporting
# ---------------------------------------------------------------------------


def load_baseline() -> dict:
    return json.loads(BASELINE_PATH.read_text())


def write_baseline(results: list[ParityResult]) -> None:
    payload = {
        "note": (
            "Measured baseline: what the parity harness catches on the app as it "
            "stands, before any fix. Regenerate with "
            "`python streamlit_app/tests/parity/harness.py --update-baseline`. "
            "Changing a number here is a deliberate act -- test_parity.py asserts "
            "against it, so an unexplained change is a regression."
        ),
        "cases": {r.case: r.to_json() for r in results},
    }
    BASELINE_PATH.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def report(results: list[ParityResult]) -> str:
    lines = [
        "| case | arm | rows | pipeline PASS | app PASS | pipeline-only | app-only | cols (pipe+extras) |",
        "|---|---|---:|---:|---:|---:|---:|---|",
    ]
    for r in results:
        if r.pipeline_error or r.app_error:
            pipe = r.pipeline_error or str(r.pipeline_pass)
            app = r.app_error or str(r.app_pass)
            lines.append(
                f"| `{r.case}` | {r.arm} | {r.rows} | {pipe} | {app} | — | — | — |"
            )
            continue
        n = len(r.pipeline_columns)
        cols = f"{n}+{len(r.app_columns) - n}" if r.columns_in_parity else f"{n}/DIVERGED"
        lines.append(
            f"| `{r.case}` | {r.arm} | {r.rows} | {r.pipeline_pass} | {r.app_pass} "
            f"| {r.only_pipeline} | {r.only_app} | {cols} |"
        )

    lines.append("")
    lines.append("### Attribution")
    lines.append("")
    lines.append("| case | divergence | rows |")
    lines.append("|---|---|---:|")
    for r in results:
        for name, count in r.attribution.items():
            lines.append(f"| `{r.case}` | {name} | {count} |")

    totals: dict[str, int] = {}
    for r in results:
        for name, count in r.attribution.items():
            totals[name] = totals.get(name, 0) + count
    lines.append("")
    lines.append("### Divergence witnessed across all cases")
    lines.append("")
    lines.append("| divergence | total diverging rows |")
    lines.append("|---|---:|")
    for name, count in sorted(totals.items(), key=lambda kv: -kv[1]):
        lines.append(f"| {name} | {count} |")

    in_parity = [r.case for r in results if r.in_parity]
    lines.append("")
    lines.append(
        f"**{len(in_parity)} of {len(results)} cases in row parity**"
        + (f": {', '.join(in_parity)}" if in_parity else "")
    )
    cols_parity = [r.case for r in results if r.columns_in_parity]
    lines.append(
        f"**{len(cols_parity)} of {len(results)} cases in column parity**"
        + (f": {', '.join(cols_parity)}" if cols_parity else "")
    )
    return "\n".join(lines)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--report", action="store_true", help="print the baseline table")
    parser.add_argument(
        "--update-baseline", action="store_true", help="rewrite baseline.json"
    )
    parser.add_argument("--case", action="append", help="limit to named case(s)")
    args = parser.parse_args(argv)

    if not PIPELINE_AVAILABLE:
        print(f"pipeline bin/ not found at {BIN_DIR}", file=sys.stderr)
        return 2

    cases = CASES
    if args.case:
        wanted = set(args.case)
        cases = [c for c in CASES if c.name in wanted]
        missing = wanted - {c.name for c in cases}
        if missing:
            print(f"unknown case(s): {sorted(missing)}", file=sys.stderr)
            return 2

    results = compare_all(cases)
    if args.update_baseline:
        if args.case:
            print("--update-baseline requires the full case set", file=sys.stderr)
            return 2
        write_baseline(results)
        print(f"wrote {BASELINE_PATH}")
    if args.report or not args.update_baseline:
        print(report(results))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
