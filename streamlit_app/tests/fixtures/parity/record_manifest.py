"""Write ``MANIFEST.json`` for a fixture directory by *running the pipeline* over it.

``test_filter_entry_point.py`` uses the manifest as its oracle: for each named row of the
constructed fixtures it holds the verdict ``bin/filter_variants.py`` reached at contract
defaults, and the test asserts the app reaches the same one. The rows are addressed
**positionally** — the Nth key of ``expected`` is the Nth row of the file — so the manifest
must be written in row order, from the same build that produced the MAF.

Recorded rather than declared, for the reason the fixture README gives: a verdict written
by hand is a second implementation of the filter, and it will be wrong. That is also why
this is a second step rather than something ``build_parity_fixtures.py`` does — the
generator can say what it *built*, but only a pipeline run can say what the pipeline makes
of it.

The manifest carries one claim the tests do not read and the export does:
``"provenance": "generator-constructed"`` per MAF, with the sha256 of the bytes it was
recorded from. The public export refuses to export a MAF that is not recorded that
way in a manifest **in its own directory**, because the bytes that must not travel are, by
content, indistinguishable from the bytes that must — five of the seven fixtures this set
replaced passed every content pattern while carrying real calls. A claim about where a file
came from is the only check that separates them, and the checksum is what stops the claim
being transferred to different bytes by a rename.

Usage::

    python streamlit_app/tests/fixtures/parity/record_manifest.py
    python streamlit_app/tests/fixtures/parity/record_manifest.py --fixture-dir DIR
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
STREAMLIT_APP = HERE.parent.parent.parent
sys.path.insert(0, str(STREAMLIT_APP))

from tests.parity import contract as C  # noqa: E402
from tests.parity import harness as H  # noqa: E402

#: The one provenance value the public export accepts for a MAF. Spelled here as
#: a literal rather than imported: the export tool is outside the app's import path, and a
#: constant that travels by copy is the point — if the two ever disagree the export fails
#: loudly rather than passing something through.
CONSTRUCTED = "generator-constructed"

#: ``fixture -> [(manifest verdict key, case overrides)]``. Mirrors
#: ``test_filter_entry_point.VERDICT_FIXTURES``: the CIViC fixture records two verdicts per
#: row, one per ``skip_civic`` setting.
VERDICT_KEYS = {
    "somatic_synthetic.maf": [("filter", {})],
    "germline_synthetic.maf": [("filter", {})],
    "somatic_civic.maf": [
        ("filter_skip_civic_false", {"skip_civic": False}),
        ("filter_skip_civic_true", {"skip_civic": True}),
    ],
}

#: ``kind`` is descriptive documentation — what role the fixture plays — and is deliberately
#: **not** what the export gate reads. Three of these words read as "constructed" and every
#: file the old set shipped under them was in fact built from real rows, so a gate trusting
#: them would have passed five of seven. ``provenance`` carries the claim; this carries the
#: description. No entry says "real subsample" any more, because nothing here is one.
KINDS = {
    "somatic_reference.maf": "constructed reference replacement",
    "germline_reference.maf": "constructed reference replacement",
    "somatic_synthetic.maf": "synthetic",
    "germline_synthetic.maf": "synthetic",
    "somatic_civic.maf": "grafted",
    "somatic_gnomad_genome.maf": "grafted",
    "somatic_dot_numeric.maf": "robustness",
}

ARMS = {
    "somatic_reference.maf": "somatic",
    "germline_reference.maf": "germline",
    "somatic_synthetic.maf": "somatic",
    "germline_synthetic.maf": "germline",
    "somatic_civic.maf": "somatic",
    "somatic_gnomad_genome.maf": "somatic",
    "somatic_dot_numeric.maf": "somatic",
}


def cell_names(path: Path) -> list[str]:
    """The constructed-cell name of every row, in row order.

    The generator parks it in ``Otherinfo`` — a ``KEEP`` column the filters never read, so
    it survives into the pipeline's output and stays readable in a failure message.

    NOT ``comment="#"``, and the difference is silent rather than loud. That option strips
    everything after a ``#`` **anywhere on a line**, not only at the start of one, and this
    fixture set carries ``#`` inside cell values on purpose: CancerVar writes its tier into
    the evidence string as ``-2#Tier_IV_benign EVS=[...]``. Read that way, the three
    ``cancervar_evidence_*`` rows lose their ``Otherinfo`` cell and the manifest records
    ``"nan"`` for them — a positional oracle that has quietly stopped naming three of its
    rows. Issue #187 recorded the same trap truncating ``CancerVar`` to
    ``1``, and ``read_maf`` avoids it the way this does, by counting the ``#``-prefixed
    lines and passing that as ``header=``.
    """
    import pandas as pd

    lines = path.read_text().splitlines()
    frame = pd.read_csv(
        path,
        sep="\t",
        header=sum(1 for line in lines if line.startswith("#")),
        low_memory=False,
        usecols=["Otherinfo"],
    )
    names = [str(v).replace("constructed_cell=", "") for v in frame["Otherinfo"]]
    unnamed = [index for index, name in enumerate(names) if name in ("nan", "")]
    if unnamed:
        raise SystemExit(
            f"{path.name}: rows {unnamed} carry no cell name in Otherinfo, so the "
            "positional oracle cannot name them"
        )
    return names


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixture-dir", type=Path, default=HERE)
    args = parser.parse_args(argv)

    C.FIXTURE_DIR = args.fixture_dir
    H.FIXTURE_DIR = args.fixture_dir
    H._read_maf_cached.cache_clear()

    fixtures: dict = {}
    for name, arm in ARMS.items():
        path = args.fixture_dir / name
        lines = path.read_text().splitlines()
        record = {
            "arm": arm,
            "bytes": path.stat().st_size,
            "columns": len(
                lines[sum(1 for line in lines if line.startswith("#"))].split("\t")
            ),
            "provenance": CONSTRUCTED,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        }
        cells = cell_names(path)
        record["rows"] = len(cells)
        record["kind"] = KINDS[name]
        # ``rows_named`` is the positional index the per-row oracle is addressed through:
        # the Nth entry is the Nth row of the file. Written for every fixture, not only the
        # ones carrying verdicts, so a diff can name the row that moved.
        record["rows_named"] = cells

        if name in VERDICT_KEYS:
            expected: dict = {cell: {} for cell in cells}
            for verdict_key, overrides in VERDICT_KEYS[name]:
                case = C.Case(f"record_{name}_{verdict_key}", name, arm, overrides=overrides)
                result = H.run_pipeline(case)
                if result.error:
                    raise SystemExit(f"{name}: pipeline failed: {result.error}")
                verdicts = list(result.verdicts.values())
                if len(verdicts) != len(cells):
                    raise SystemExit(
                        f"{name}: {len(verdicts)} verdicts for {len(cells)} rows — the "
                        "positional addressing the manifest relies on is broken"
                    )
                for cell, verdict in zip(cells, verdicts):
                    expected[cell][verdict_key] = verdict
            record["expected"] = expected
        fixtures[name] = record

    for name in ["genes_somatic.txt", "genes_germline.txt", "genes_somatic_mixed_case.txt"]:
        path = args.fixture_dir / name
        fixtures[name] = {
            "bytes": path.stat().st_size,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "genes": len(path.read_text().split()),
        }

    manifest = {
        "fixtures": dict(sorted(fixtures.items())),
        "parameters": {
            key: C.CONTRACT[key]
            for key in sorted(C.CONTRACT)
            if key not in {"genes", "max_freq_population"}
        },
    }
    (args.fixture_dir / "MANIFEST.json").write_text(json.dumps(manifest, indent=1) + "\n")
    named = sum(len(r.get("expected", {})) for r in fixtures.values())
    verdicts = sum(
        len(v) for r in fixtures.values() for v in r.get("expected", {}).values()
    )
    print(f"wrote {args.fixture_dir / 'MANIFEST.json'}: {named} named rows, {verdicts} verdicts")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
