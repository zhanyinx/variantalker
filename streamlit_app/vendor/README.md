# `vendor/` — the pipeline's code and data, copied

This directory is **not** app code. It is a copy of parts of the VarianTalker
pipeline, held here so that MAFigate's PASS/NOPASS verdict *is* the pipeline's
verdict rather than a re-implementation of it — and, since issue 198, so that the
markers CancerVar cites *are* the ones it cited.

| File | What it is | Origin |
|---|---|---|
| `pipeline_utils.py` | the entire file, **byte-for-byte** | `bin/utils.py` |
| `pipeline_filters.py` | five filter functions verbatim, plus `compute_keep` — see below | `bin/filter_variants.py` |
| `cancervar_markers.txt` | the entire file, **byte-for-byte** (1.2 MB of data, not code) | `resources/CancerVar/cancervardb/cancervar.out.txt` |
| `_sync.py` | the developer-side copy/repair tool (not shipped in installers) | — |

**Two source trees, so absence is per source.** The first two units come from `bin/`;
the marker table comes from `resources/CancerVar/`. A checkout can have one without
the other, so `_sync.py`'s `GUARDED_SOURCES` skips them independently rather than
skipping everything the moment `bin/` is missing.

Five of the six are plain functions. The sixth, `compute_keep`, is the output-column
computation — which is not a function in the pipeline at all, but eight statements
strewn through `main()`. `_sync.py` lifts five of them into a `def` so the app can
call it, and the pipeline's own parameter names (`args`, `out`) survive as the
signature: callers supply an argparse-like namespace and the frame as shims.

### `compute_keep` is the one place the copy is deliberately unfaithful

`main()` opens with `keep = KEEP`, which **aliases** the module-level list from
`utils.py` and then mutates it in place. A pipeline process runs `main()` once and
exits, so it never sees the consequence. The app keeps the module loaded for the whole
session, where the mutations accumulate: a germline run leaves the shared list three
entries short, so a following **somatic** run silently returns germline columns, and a
second germline run raises `ValueError: list.remove(x): x not in list`.

`compute_keep` therefore opens with `keep = list(KEEP)` instead. The copy is taken
*before* the verbatim statements rather than by editing them, so every statement the
guard compares is still an unedited copy. `tests/test_vendor_compute_keep.py`
reproduces both symptoms against the fixed function.

## Do not edit these files

Change `bin/` and re-sync:

```bash
python streamlit_app/vendor/_sync.py          # repair the copies
python streamlit_app/vendor/_sync.py --check   # report drift, change nothing
make -C streamlit_app check-vendor             # same check, via make
```

### Which `bin/` is the authority

`bin/` is authoritative for `vendor/`, and the **private tree** is authoritative for
`bin/`. This file is exported to the public repository unchanged, so it is read from both
trees — and only in the private one does "change `bin/` and re-sync" end with a change
that survives. A public clone receives each release as one squashed commit that replaces
the working tree wholesale, so an edit to `bin/` there is overwritten by the next export
rather than merged; see [`CONTRIBUTING.md`](../../CONTRIBUTING.md), which spells out what
happens to a pull request. Propose the change there and it is reapplied by hand into the
private tree, where `bin/` moves and `_sync.py` follows it.

## Why copies and not an import

The app ships as a standalone `.dmg` / `.exe` that carries no `bin/`, so importing the
pipeline at runtime is not an option. And `bin/filter_variants.py` cannot be imported
even in a full checkout: its line 10 pulls in `database.database_utils`, and with it
**`psycopg`**, a Postgres driver the packaged app will never have. `bin/utils.py` has
no such problem — it imports only `os`, `numpy`, `re` and `pandas` — so it is copied
whole, and only `filter_variants.py` needs symbol extraction.

The marker table is here for the same reason and it is the plainer case: it is a data
file under `resources/CancerVar/`, which the `.dmg` / `.exe` does not carry at all. The
alternative issue 198 weighed — read `resources/` where it exists, show nothing where
it does not — would have put the marker detail in developer checkouts only, which is to
say not in the product.

Copying is a real risk, taken deliberately and with a guard. It has already gone wrong
here once without one: `streamlit_app/utils/main_utils.py` **held** a hand-copied `KEEP`
that drifted to **39 entries against the pipeline's 45**, with different membership and
different order, while nothing read it. That copy is gone — issue 32 deleted it, and
`main_utils.py`'s docstring is where the history is now recorded. The live list is
`vendor.pipeline_utils.KEEP`, which this guard holds to the pipeline's 45.

## How drift is caught

Four mechanisms, one per unit:

- **`pipeline_utils.py`** — whole-file sha256 against `bin/utils.py`. Intolerant of even
  cosmetic edits, which is fine because the repair is one command.
- **`cancervar_markers.txt`** — whole-file sha256 against the pipeline's copy, plus one
  test that row 0 is still the header. Here intolerance is the *point*, not a tolerable
  side effect. A MAF's `Therap_list` / `Diag_list` / `Prog_list` cell holds **0-based
  line offsets** into this file — `CancerVar.py:193` reads it with
  `list(csv.reader(...))`, header included, and `:1091` indexes straight into the
  result — so a copy of a different vintage does not fail to resolve. Every index after
  an inserted row resolves to the *neighbouring* marker, and the app names a drug the
  file never associated with the variant. There is no cosmetic edit to a data file.

  The header test is not redundant with the hash: a header stripped from *both* copies
  leaves the sha256 matching while shifting every index by one.
  `test_cancervar_markers_still_has_its_header_at_row_zero` is what notices, and the
  redundancy was checked by mutation rather than assumed.
- **the five filter functions** — per-symbol `ast.dump()` equality against the matching
  function in `bin/filter_variants.py`. Tolerates comment and formatting churn in
  `bin/`; catches anything semantic.
- **`compute_keep`** — an eight-statement AST slice of `main()`. Five statements are
  vendored and compared one by one. The other three are guarded **without** being
  vendored, because nothing else would notice them changing:
  - `keep = KEEP`, the aliasing statement `compute_keep` replaces. If the pipeline ever
    changes how it seeds the list, the replacement stops being equivalent.
  - `out_nopass = out_nopass[keep]` and `out = out[keep]`, the two reselections that
    apply the list. They are what make it *the pipeline's output columns* rather than a
    list the pipeline happens to build. Delete one and the pipeline stops narrowing
    that frame, while the app carries on narrowing it to a set the pipeline no longer
    emits.

  `_sync.py` locates the eight by source-line anchor; `tests/test_vendor_drift.py`
  finds them a different way, by searching `main()` for an AST match. A missing or
  duplicated statement is a hard error in both, never a silently shorter slice.

The guard has been shown to fire, not just assumed to: `test_column_guard_fires_on`
plants a throwaway repo, renames a column, mistypes one, reseeds the list and alters
and deletes a reselection, and requires both forms to go red on each.
`test_column_guard_tolerates` requires both to stay green on comments, re-wrapped
lines, and edits elsewhere in `main()`.

Four enforcement points, because a test nobody runs is what produced the drift above:

1. `streamlit_app/tests/test_vendor_drift.py` — runs with the normal suite.
2. `.github/workflows/vendor-drift.yml` — on every push and pull request. Runs both
   forms: the stdlib-only gate, then the pytest guard. It also runs
   `tests/test_vendor_compute_keep.py`, the one behaviour suite in this job: the guard
   proves `compute_keep` is still the pipeline's, and that proves the one place it
   deliberately is not. It needs pandas and nothing else, so the job stays narrow.
3. `make build-mac` / `make build-win` — both depend on `check-vendor`, so an installer
   cannot be cut from a drifted copy.
4. `.pre-commit-config.yaml` — fires when `bin/` or `vendor/` changes.

Points 3 and 4 run `_sync.py --check` and nothing else — the Makefile through
`$(PYTHON)` (`$MAFIGATE_PYTHON`, else this directory's `.venv`, else a bare `python3`),
the pre-commit hook through plain `python`. Either interpreter will do, because the
script is stdlib-only — no pytest, no pandas — so a bare build machine needs nothing
installed. Because the rule is therefore implemented twice,
`test_sync_check_agrees_with_this_guard` pins the two together: if they ever disagree, one
of them is lying about whether the app still matches the pipeline.

**Where a source tree is absent, every enforcement point skips rather than fails** —
the packaged app and app-only checkouts stay green, because there the comparison is
impossible to make rather than failed. The skip is **per source**: a tree with `bin/`
but no `resources/CancerVar/` still checks the filter code. The one exception is
`_sync.py` in *repair* mode, which still exits non-zero: with no `bin/` to copy from,
succeeding would claim work it never did.

### Telling a skip from a check

That skip is the guard's most dangerous state, because a guard that compared nothing is
still green. So green is not the thing to read. Both forms report **how many sources they
compared**, and that is:

- **`_sync.py --check`** ends in `vendored copies are in sync (2 of 2 sources compared)`.
  The count is the message: `1 of 2` means one tree was absent, and a `SKIP` line above
  names which. When *every* source is missing it prints `The vendored copies were not
  compared against the pipeline.` and never the words "in sync" —
  `test_sync_check_skips_cleanly_without_bin` pins that wording, because a build log must
  not read an impossible comparison as a passing one. Both shapes exit 0.
- **the pytest guard** says it in the **skip count** of its summary line, and nowhere
  else. Measured on this tree, with the two files CI runs
  (`tests/test_vendor_drift.py tests/test_vendor_compute_keep.py`): a full checkout is
  `50 passed` and **nothing skipped**; `bin/` present but no `resources/CancerVar/` is
  `48 passed, 2 skipped`, the two marker-table tests; neither tree present is
  `20 passed, 30 skipped`. So `50 passed` with a bare zero skips is the only shape that
  means the whole copy was compared, and any skip count at all means less than everything
  was. `-rs` names which; the workflow's `-v` shows each one.

## Two things that would quietly break the copy

Both are handled; this is a note for whoever changes the tooling next.

- **The formatter.** `bin/utils.py` is not black-formatted at width 100 — `black
  --line-length 100` rewrites 296 of its lines — and its imports are not isort-ordered.
  `streamlit_app/pyproject.toml` therefore excludes `vendor/` from both. Without that,
  the first `make format` after vendoring turns the sha256 guard red for a reason that
  has nothing to do with drift.
- **Coverage.** These ~1,100 lines are exercised by the app, not by unit tests, and
  would dominate any coverage measurement of code we actually own.
  `streamlit_app/.coveragerc` omits them. That file, not `pytest.ini`, is the only place
  pytest-cov reads coverage settings from.
