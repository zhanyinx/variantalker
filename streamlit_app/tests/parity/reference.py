"""Run the parity harness over the full 50-sample GERSOM reference.

Issue #17 built the harness against ten committed fixtures because the reference
was not mounted on the machine that built it. The fixtures are a *purposive*
subsample — one named predicate per divergence, first-N matching rows — which is
what gives them power and also what makes their per-divergence row counts
estimates of nothing. This module runs the same contract and the same probes over
real data so those counts become rates.

Three decisions are baked in here, each forced by something measured:

* **Per file, never pooled.** The four-column comparison key is unique within any
  single MAF (issue #18: 0 duplicates across all 100 files) but collides on
  78,728 of 84,130 rows once the arm is pooled, and germline cannot be rescued by
  ``Tumor_Sample_Barcode`` (``__UNKNOWN__`` in all 92,547 rows). Comparing per
  file sidesteps the problem entirely rather than fabricating a key. It also
  keeps every ``iterrows`` attribution pass bounded by one sample's ~1,800 rows.

* **The contract and the probes are imported, not restated.** Only two values
  move: ``Case.fixture`` becomes an absolute reference path, and the gene cases'
  12-symbol fixture list becomes the project's real 36/52-symbol GERSOM list. A
  gene list drawn from a purposive subsample would not produce a rate either.

* **Seven cases do not map.** ``somatic_synthetic``, ``germline_synthetic``,
  ``somatic_synthetic_depth_0``, ``civic_present``, ``civic_skipped``,
  ``gnomad_genome_present`` and ``dot_numeric`` name shapes the reference does not
  contain — no CIViC column in 100/100 files, no missing numerics, no
  gnomAD_genome. They stay fixture-only, and that asymmetry is a result, not an
  omission.

Usage — ``PARITY_REFERENCE_ROOT`` says where the reference is mounted, and has no default
(see :data:`REFERENCE_ROOT_ENV`)::

    export PARITY_REFERENCE_ROOT=/path/to/the/reference/root

    python streamlit_app/tests/parity/reference.py --report
    python streamlit_app/tests/parity/reference.py --out results.json --jobs 8
    python streamlit_app/tests/parity/reference.py --samples 3   # smoke run

    # or, for one run, without the variable:
    python streamlit_app/tests/parity/reference.py --reference-root /path/to/root --report
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field, replace
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
STREAMLIT_APP = HERE.parent.parent

if str(STREAMLIT_APP) not in sys.path:
    sys.path.insert(0, str(STREAMLIT_APP))

from tests.parity import harness  # noqa: E402
from tests.parity.contract import CASES, Case  # noqa: E402

#: The environment variables that say where the reference data is mounted. There is no
#: default, and that is the change issue #247 made rather than a gap.
#:
#: This used to hold the mount path itself, spelled out: a cloud-drive path naming the
#: institute, the unit, the project and the manuscript the clinical data belongs to. It
#: was the last thing in the tree that the publication scan gate refused to export, and
#: the two ways out were to strip this file from the public tree or to stop keeping the
#: path in it. Stripping loses the module — the sweep is the evidence behind every rate
#: the parity README quotes, and a public reader who cannot see how it was produced
#: cannot check it — so the path went instead. The mount is a fact about one machine,
#: which is what an environment variable is for.
#:
#: Nothing else here needed changing: :func:`main` already reported only ``root.name`` in
#: its results, never the path.
REFERENCE_ROOT_ENV = "PARITY_REFERENCE_ROOT"
RESOURCE_ROOT_ENV = "PARITY_RESOURCE_ROOT"

GENE_FILES = {
    "somatic": "gersom_gene_somatic_V1.csv",
    "germline": "gersom_gene_germline_V1.csv",
}

#: Fixtures with a reference counterpart. Everything else is a constructed shape.
REFERENCE_FIXTURES = {
    "somatic_reference.maf": "somatic",
    "germline_reference.maf": "germline",
}


def reference_root() -> Path:
    """Where the 50-sample reference is mounted, from the environment. No default.

    Raises rather than returning a placeholder that would fail one frame later as "not
    found: ." — the message has to say what to set, because the person meeting it is
    typically running this module for the first time on a machine that has the drive.
    """
    value = os.environ.get(REFERENCE_ROOT_ENV)
    if not value:
        raise SystemExit(
            f"{REFERENCE_ROOT_ENV} is not set, and there is no default: the reference "
            "lives on a clinical shared drive whose path names the institute, the unit "
            "and the project, so this tree does not carry it.\n"
            f"  export {REFERENCE_ROOT_ENV}=/path/to/the/reference/root\n"
            "or pass --reference-root. The root is the directory holding the per-arm "
            f"subdirectories; {RESOURCE_ROOT_ENV} defaults to its sibling `resources`."
        )
    return Path(os.path.expanduser(value))


def resource_root(reference: Path | None = None) -> Path:
    """The project's own gene lists — the ones the reference run itself was filtered with.

    ``.csv`` by name, newline-delimited in fact (issue #14). Defaults to the reference
    root's sibling, which is where they sit on the drive.

    `reference` is how ``--reference-root`` reaches this: without it, a run that passed
    the flag and set no environment variable would fall through to :func:`reference_root`
    and be told to set ``PARITY_REFERENCE_ROOT`` — a variable it had just been given a
    perfectly good alternative to.

    It is optional rather than required because the no-argument call has real callers:
    the measurement scripts in the private notes reach for ``resource_root()``
    directly, having no argparse of their own. Read as an unused parameter and deleted,
    that breaks them.
    """
    value = os.environ.get(RESOURCE_ROOT_ENV)
    if value:
        return Path(os.path.expanduser(value))
    return (reference or reference_root()).parent / "resources"


# ---------------------------------------------------------------------------
# Discovery
# ---------------------------------------------------------------------------


def discover_samples(root: Path, arm: str) -> list[tuple[str, Path]]:
    """``(sample_name, maf_path)`` for one arm, in stable sorted order."""
    arm_dir = root / arm
    if not arm_dir.is_dir():
        return []
    found = []
    for sample_dir in sorted(p for p in arm_dir.iterdir() if p.is_dir()):
        maf = sample_dir / f"{sample_dir.name}.maf"
        if maf.exists():
            found.append((sample_dir.name, maf))
    return found


def mappable_cases() -> list[Case]:
    """Cases this sweep can retarget at a reference MAF, and can actually distinguish.

    Two exclusions, for different reasons. A case on a constructed fixture has no
    reference counterpart to be retargeted at. And a case whose ``genes_as`` is a paste
    mode has nothing to add *here*: this module configures the vendored functions from the
    gene **file path** directly, so it never runs the app's tokeniser and a pasted case is
    byte-for-byte its pre-split sibling. Including them would inflate the measurement count
    by cases that cannot witness anything, which is the exact fault issue #27 fixed for the
    unwitnessed probes. Their power lives in ``harness.py``, where the app side goes through
    ``apply_filters``.
    """
    return [
        c
        for c in CASES
        if c.fixture in REFERENCE_FIXTURES and c.genes_as == "tokens"
    ]


def fixture_only_cases() -> list[Case]:
    """Cases pinned to a constructed fixture, which the reference cannot reproduce.

    Not the complement of :func:`mappable_cases` — a paste-mode case on a reference
    fixture is in neither, because it *is* reproducible here and simply adds nothing.
    """
    return [c for c in CASES if c.fixture not in REFERENCE_FIXTURES]


def mixed_case_gene_file(source: Path, dest_dir: Path) -> Path:
    """The real gene list with its case mangled the way the fixture list is.

    ``genes_somatic_mixed_case.txt`` alternates upper and title case so that a
    case-sensitive comparison misses every symbol while a case-insensitive one
    matches all of them. Reproduced here against the 36-symbol GERSOM list rather
    than reused, because the fixture list holds only 12 symbols.
    """
    tokens = [t.strip() for t in source.read_text().splitlines() if t.strip()]
    mangled = [t.title() if i % 2 else t.lower() for i, t in enumerate(tokens)]
    dest = dest_dir / f"{source.stem}_mixed_case.txt"
    dest.write_text("\n".join(mangled) + "\n")
    return dest


def build_reference_case(case: Case, maf: Path, gene_paths: dict[str, Path]) -> Case:
    """``case`` retargeted at one reference MAF, with the real gene list.

    ``Case.maf_path`` is ``FIXTURE_DIR / fixture``; pathlib's ``/`` yields the
    right-hand side unchanged when it is absolute, so substituting an absolute
    path is all that is needed. Same for ``genes_path``.
    """
    overrides = dict(case.overrides)
    genes = overrides.get("genes")
    if genes is not None:
        key = "mixed_case" if "mixed_case" in genes else "plain"
        overrides["genes"] = str(gene_paths[f"{case.arm}_{key}"])
    return replace(case, fixture=str(maf), overrides=overrides)


# ---------------------------------------------------------------------------
# Candidate probes
# ---------------------------------------------------------------------------
#
# ``harness.PROBES`` is imported and run unchanged — that is what makes the
# reference numbers comparable to the fixture baseline. Anything below is
# measured *separately*, against the rows the unchanged probes leave in the
# ``unattributed`` bucket, so a candidate can be shown to explain the gap without
# retroactively editing the instrument that found it.


#: Empty, and deliberately kept rather than deleted. Issue #25 measured divergence
#: #6 here, as a candidate, because the fixtures held no row that witnessed it and
#: editing ``harness.PROBES`` to explain a bucket found by ``harness.PROBES`` would
#: have been circular. Issue #27 added the ``germline_escat_only`` fixture cell and
#: promoted the probe into the harness, so there is nothing left to hold — but the
#: separation is the reason the blind spot was findable at all, and the next
#: divergence the reference turns up needs somewhere to be measured before it earns
#: a place in the instrument.
CANDIDATE_PROBES: dict = {}


# ---------------------------------------------------------------------------
# The four-cell decomposition (issue #20)
# ---------------------------------------------------------------------------


def decompose(case: Case, frame: pd.DataFrame) -> dict[str, int]:
    """Criteria-only / patho-only / both / NOPASS, from the vendored masks.

    Issue #20 made this the app's diagnostic output precisely because a PASS
    total cannot see the difference between the two paths into the report — 157
    of 818 germline report rows exist solely by pathogenic rescue. It is also
    what reconciles issue #16's 411 against issue #13's 408: they are the union
    and the criteria path of the same run.

    Cross-validated against the pipeline subprocess by the caller: ``union`` here
    must equal the subprocess's PASS set, or the reconstruction is wrong and
    nothing downstream of it can be trusted.
    """
    import contextlib
    import io

    from vendor.pipeline_filters import common_filters, germline_filters, somatic_filters

    c = case.contract
    genes = "null" if case.genes_path is None else str(case.genes_path)

    # `check_civic_column_exists` prints to stdout on every call, and CIViC is absent
    # in 100/100 reference files. Unsuppressed that is 1,100 lines of warning on top
    # of the report; the vendored body is verbatim, so it is muzzled here, not edited.
    stdout = io.StringIO()
    with contextlib.redirect_stdout(stdout):
        common = common_filters(
            frame, c["min_depth"], c["filter_variant_classification"]
        )
        specific, patho = _arm_filters(
            case, frame, c, genes, somatic_filters, germline_filters
        )

    criteria = common & specific
    union = patho | criteria
    return {
        "criteria_only": int((criteria & ~patho).sum()),
        "both": int((criteria & patho).sum()),
        "patho_only": int((patho & ~criteria).sum()),
        "criteria_total": int(criteria.sum()),
        "patho_total": int(patho.sum()),
        "union_pass": int(union.sum()),
        "nopass": int((~union).sum()),
    }, set(harness.variant_keys(frame[union]))


def _arm_filters(case, frame, c, genes, somatic_filters, germline_filters):
    if case.arm == "somatic":
        return somatic_filters(
            frame,
            c["vaf_threshold"],
            genes,
            c["filter_cancervar"],
            c["filter_civic"],
            c["filter_escat"],
            c["filter_clinvar"],
            c["skip_civic"],
            c["skip_pathogenic"],
        )
    return germline_filters(
        frame,
        c["vaf_threshold_germline"],
        genes,
        c["filter_intervar"],
        c["filter_renovo"],
        c["filter_clinvar"],
        c["skip_pathogenic"],
    )


# ---------------------------------------------------------------------------
# One (case, sample) measurement
# ---------------------------------------------------------------------------


@dataclass
class Measurement:
    case: str
    arm: str
    sample: str
    rows: int
    pipeline_pass: int | None = None
    app_pass: int | None = None
    only_pipeline: int = 0
    only_app: int = 0
    pipeline_error: str | None = None
    app_error: str | None = None
    attribution: dict[str, int] = field(default_factory=dict)
    #: What CANDIDATE_PROBES explain of the rows PROBES left unattributed.
    candidate_attribution: dict[str, int] = field(default_factory=dict)
    still_unattributed: int = 0
    cells: dict[str, int] = field(default_factory=dict)
    decomposition_matches_pipeline: bool | None = None
    unattributed_keys: list = field(default_factory=list)

    def to_json(self, sample_alias=None, include_row_keys: bool = False) -> dict:
        """The committed record. Counts by default, never raw coordinates.

        ``unattributed_keys`` holds ``(Chromosome, Start_Position,
        Reference_Allele, Tumor_Seq_Allele2)`` for real **germline** calls, keyed
        to a real sample id. That is patient genetic data, and issue #18
        pseudonymised barcodes and scrubbed HPC paths from the fixtures for the
        same reason. So the identity stays local: sample names are aliased and
        the coordinates are reduced to a count unless explicitly asked for.
        """
        payload = self.__dict__.copy()
        payload["sample"] = (
            sample_alias(self.sample) if sample_alias else self.sample
        )
        if include_row_keys:
            payload["unattributed_keys"] = [list(k) for k in self.unattributed_keys]
        else:
            payload["unattributed_keys"] = len(self.unattributed_keys)
        return payload


def measure(case: Case, sample: str) -> Measurement:
    """Both sides over one reference MAF, diffed, attributed and decomposed."""
    frame = harness.load_frame(case)
    m = Measurement(case=case.name, arm=case.arm, sample=sample, rows=len(frame))

    pipeline = harness.run_pipeline(case)
    app = harness.run_app(case, frame.copy())
    m.pipeline_error = pipeline.error
    m.app_error = app.error

    if pipeline.error is None:
        m.pipeline_pass = len(pipeline.passed)
        cells, union_keys = decompose(case, frame)
        m.cells = cells
        m.decomposition_matches_pipeline = union_keys == pipeline.passed
    if app.error is None:
        m.app_pass = len(app.passed)

    if pipeline.error is None and app.error is None:
        only_pipeline = pipeline.passed - app.passed
        only_app = app.passed - pipeline.passed
        m.only_pipeline = len(only_pipeline)
        m.only_app = len(only_app)
        diff = only_pipeline | only_app
        m.attribution = harness.attribute(case, frame, diff)
        if m.attribution.get("unattributed"):
            m.unattributed_keys, m.candidate_attribution, m.still_unattributed = (
                _examine_unattributed(case, frame, diff)
            )
    return m


def _examine_unattributed(case: Case, frame: pd.DataFrame, diff: set):
    """The diverging rows no probe explains — the fixture set's blind spot.

    ``harness.attribute`` counts them but discards their identity. On the
    fixtures the bucket was empty for every case, so it never mattered; on the
    reference it is the single most informative thing the sweep produces. Each
    such row is then offered to ``CANDIDATE_PROBES``, so "the bucket is explained
    by divergence #6" is a measurement rather than a reading of six rows.
    """
    c = case.contract
    all_keys = harness.variant_keys(frame)
    positions = [i for i, k in enumerate(all_keys) if k in diff]

    keys: list = []
    candidates: dict[str, int] = {}
    still = 0
    for pos in positions:
        row = frame.iloc[pos]
        if any(probe(row, c, case) for probe in harness.PROBES.values()):
            continue
        keys.append(all_keys[pos])
        hit = False
        for name, probe in CANDIDATE_PROBES.items():
            if probe(row, c, case):
                candidates[name] = candidates.get(name, 0) + 1
                hit = True
        if not hit:
            still += 1
    return keys, candidates, still


def _worker(payload):
    case, sample = payload
    try:
        return measure(case, sample)
    except Exception as exc:  # noqa: BLE001 - a crashed cell is a result too
        return Measurement(
            case=case.name,
            arm=case.arm,
            sample=sample,
            rows=-1,
            pipeline_error=f"harness crash: {type(exc).__name__}: {exc}",
        )


# ---------------------------------------------------------------------------
# The sweep
# ---------------------------------------------------------------------------


def build_workload(root: Path, gene_paths: dict[str, Path], limit: int | None):
    work = []
    samples = {arm: discover_samples(root, arm) for arm in ("somatic", "germline")}
    if limit:
        samples = {arm: s[:limit] for arm, s in samples.items()}
    for case in mappable_cases():
        for sample, maf in samples[case.arm]:
            work.append((build_reference_case(case, maf, gene_paths), sample))
    return work, samples


def run_sweep(root: Path, gene_paths: dict, limit: int | None, jobs: int):
    work, samples = build_workload(root, gene_paths, limit)
    total = len(work)
    print(
        f"reference: {len(samples['somatic'])} somatic + {len(samples['germline'])} "
        f"germline samples, {len(mappable_cases())} mappable cases -> {total} measurements",
        file=sys.stderr,
    )
    results = []
    if jobs == 1:
        for i, payload in enumerate(work, 1):
            results.append(_worker(payload))
            print(f"\r  {i}/{total}", end="", file=sys.stderr, flush=True)
    else:
        with ProcessPoolExecutor(max_workers=jobs) as pool:
            futures = [pool.submit(_worker, payload) for payload in work]
            for i, fut in enumerate(as_completed(futures), 1):
                results.append(fut.result())
                print(f"\r  {i}/{total}", end="", file=sys.stderr, flush=True)
    print(file=sys.stderr)
    return results


# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------


def aggregate(results: list[Measurement]) -> dict:
    """Per-case totals across samples. Sums, not means: rates need numerators."""
    by_case: dict[str, dict] = {}
    for m in results:
        agg = by_case.setdefault(
            m.case,
            {
                "arm": m.arm,
                "samples": 0,
                "rows": 0,
                "pipeline_pass": 0,
                "app_pass": 0,
                "only_pipeline": 0,
                "only_app": 0,
                "samples_in_parity": 0,
                "pipeline_errors": defaultdict(int),
                "app_errors": defaultdict(int),
                "attribution": defaultdict(int),
                "candidate_attribution": defaultdict(int),
                "still_unattributed": 0,
                "cells": defaultdict(int),
                "decomposition_mismatches": 0,
                "unattributed_samples": 0,
            },
        )
        agg["samples"] += 1
        agg["rows"] += max(m.rows, 0)
        if m.pipeline_error:
            agg["pipeline_errors"][m.pipeline_error] += 1
        else:
            agg["pipeline_pass"] += m.pipeline_pass or 0
            for k, v in m.cells.items():
                agg["cells"][k] += v
            if m.decomposition_matches_pipeline is False:
                agg["decomposition_mismatches"] += 1
        if m.app_error:
            agg["app_errors"][m.app_error] += 1
        else:
            agg["app_pass"] += m.app_pass or 0
        agg["only_pipeline"] += m.only_pipeline
        agg["only_app"] += m.only_app
        if not m.pipeline_error and not m.app_error and not m.only_pipeline and not m.only_app:
            agg["samples_in_parity"] += 1
        for k, v in m.attribution.items():
            agg["attribution"][k] += v
        for k, v in m.candidate_attribution.items():
            agg["candidate_attribution"][k] += v
        agg["still_unattributed"] += m.still_unattributed
        if m.unattributed_keys:
            agg["unattributed_samples"] += 1

    for agg in by_case.values():
        for key in (
            "pipeline_errors",
            "app_errors",
            "attribution",
            "candidate_attribution",
            "cells",
        ):
            agg[key] = dict(sorted(agg[key].items()))
    return by_case


def report(by_case: dict) -> str:
    lines = [
        "| case | arm | rows | pipeline PASS | app PASS | pipe-only | app-only | samples at parity |",
        "|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    for name, a in by_case.items():
        err = ""
        if a["pipeline_errors"] or a["app_errors"]:
            err = " ⚠"
        lines.append(
            f"| `{name}`{err} | {a['arm']} | {a['rows']} | {a['pipeline_pass']} | "
            f"{a['app_pass']} | {a['only_pipeline']} | {a['only_app']} | "
            f"{a['samples_in_parity']}/{a['samples']} |"
        )

    lines += ["", "### Four-cell decomposition (issue #20 vocabulary)", ""]
    lines += [
        "| case | criteria-only | both | patho-only | criteria total | union PASS |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for name, a in by_case.items():
        c = a["cells"]
        if not c:
            continue
        lines.append(
            f"| `{name}` | {c.get('criteria_only', 0)} | {c.get('both', 0)} | "
            f"{c.get('patho_only', 0)} | {c.get('criteria_total', 0)} | "
            f"{c.get('union_pass', 0)} |"
        )

    totals: dict[str, int] = defaultdict(int)
    for a in by_case.values():
        for k, v in a["attribution"].items():
            totals[k] += v
    lines += ["", "### Divergence rates on the full reference", ""]
    lines += ["| divergence | diverging rows (summed over cases) |", "|---|---:|"]
    for name, count in sorted(totals.items(), key=lambda kv: -kv[1]):
        lines.append(f"| {name} | {count} |")

    cand: dict[str, int] = defaultdict(int)
    still = 0
    for a in by_case.values():
        for k, v in a["candidate_attribution"].items():
            cand[k] += v
        still += a["still_unattributed"]
    if totals.get("unattributed"):
        lines += [
            "",
            "### What explains the `unattributed` bucket",
            "",
            f"The unchanged probes leave **{totals['unattributed']}** diverging rows "
            "unexplained. Offered to the candidate probes:",
            "",
            "| candidate probe | rows explained |",
            "|---|---:|",
        ]
        for name, count in sorted(cand.items(), key=lambda kv: -kv[1]):
            lines.append(f"| {name} | {count} |")
        lines.append(f"| *still unexplained* | {still} |")

    mismatches = sum(a["decomposition_mismatches"] for a in by_case.values())
    lines += [
        "",
        f"**decomposition vs pipeline subprocess: {mismatches} mismatching sample-cases**"
        " (0 means the vendored mask reconstruction reproduces the pipeline exactly)",
    ]
    return "\n".join(lines)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--reference-root", type=Path, default=None)
    parser.add_argument("--out", type=Path, default=None, help="write full JSON results")
    parser.add_argument("--samples", type=int, default=None, help="limit samples per arm")
    parser.add_argument("--jobs", type=int, default=max(os.cpu_count() - 2, 1))
    parser.add_argument("--case", action="append", help="limit to named case(s)")
    parser.add_argument("--report", action="store_true")
    parser.add_argument(
        "--include-row-keys",
        action="store_true",
        help="write raw variant coordinates for unattributed rows. Real germline "
        "calls -- keep the output local, do not commit it.",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="exit non-zero if any diverging row is unattributed. This is what "
        "`make parity-reference` runs: the invariant needs no committed numbers, "
        "so it cannot go stale, and it is exactly the property that was false "
        "(3,687 rows) while the fixture baseline read as complete. Leave it off "
        "while investigating -- an unattributed bucket is the finding, and "
        "CANDIDATE_PROBES is where the next explanation gets measured first.",
    )
    args = parser.parse_args(argv)

    root = args.reference_root or reference_root()
    if not root.is_dir():
        print(f"reference root not found: {root}", file=sys.stderr)
        return 2

    import tempfile

    gene_dir = Path(tempfile.mkdtemp(prefix="parity_reference_genes_"))
    gene_paths = {}
    for arm, filename in GENE_FILES.items():
        source = resource_root(root) / filename
        if not source.exists():
            print(f"gene list not found: {source}", file=sys.stderr)
            return 2
        plain = gene_dir / filename
        plain.write_text(source.read_text())
        gene_paths[f"{arm}_plain"] = plain
        gene_paths[f"{arm}_mixed_case"] = mixed_case_gene_file(plain, gene_dir)

    global CASES_FILTER  # noqa: PLW0603
    results = run_sweep(root, gene_paths, args.samples, args.jobs)
    if args.case:
        wanted = set(args.case)
        results = [r for r in results if r.case in wanted]

    by_case = aggregate(results)
    # Stable order: contract order, not completion order.
    order = [c.name for c in mappable_cases()]
    by_case = {k: by_case[k] for k in order if k in by_case}

    if args.out:
        aliases: dict[str, str] = {}

        def alias(sample: str) -> str:
            return aliases.setdefault(sample, f"S{len(aliases):03d}")

        for arm in ("somatic", "germline"):
            for sample, _ in discover_samples(root, arm):
                alias(sample)

        args.out.write_text(
            json.dumps(
                {
                    # The absolute path names a shared clinical drive and the
                    # operator's home directory; neither belongs in git.
                    "reference_root": root.name,
                    "samples": {"aliased": len(aliases)},
                    "row_keys_included": args.include_row_keys,
                    "by_case": by_case,
                    "measurements": [
                        m.to_json(alias, args.include_row_keys) for m in results
                    ],
                },
                indent=2,
                sort_keys=True,
                default=str,
            )
            + "\n"
        )
        print(f"wrote {args.out}", file=sys.stderr)

    print(report(by_case))

    if args.strict:
        failed = False
        unattributed = sum(a["attribution"].get("unattributed", 0) for a in by_case.values())
        if unattributed:
            print(
                f"\nFAIL: {unattributed} diverging rows no probe explains. The fixture "
                f"set cannot see them and baseline.json will not go red for them. "
                f"Measure the cause as a CANDIDATE_PROBE here first, then add the "
                f"fixture cell that witnesses it and promote the probe to "
                f"harness.PROBES.",
                file=sys.stderr,
            )
            failed = True

        # The four-cell decomposition, asserted rather than merely reported. It is a
        # strictly stronger signal than any PASS count: the cells are compared to
        # `bin/filter_variants.py` in a fresh subprocess as a *set of variant keys*,
        # per sample, so two runs cannot agree by cancelling errors the way equal
        # totals can. It costs nothing extra -- the masks are already computed -- and
        # the expectation is the constant 0, so there is no number to keep current.
        # If this fires, the vendored masks no longer reconstruct the pipeline's
        # verdict and every diagnostic built on them (issue #20's cells, issue #21's
        # `MAFigate_reason`) is reporting on a decision it did not make.
        mismatches = sum(a["decomposition_mismatches"] for a in by_case.values())
        if mismatches:
            print(
                f"\nFAIL: the four-cell decomposition disagrees with "
                f"bin/filter_variants.py on {mismatches} sample-cases. The vendored "
                f"masks no longer reconstruct the pipeline's PASS set, so issue #20's "
                f"diagnostics describe a decision the pipeline did not make.",
                file=sys.stderr,
            )
            failed = True
        if failed:
            return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
