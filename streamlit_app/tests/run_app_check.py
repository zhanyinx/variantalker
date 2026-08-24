"""Does the running app still load and browse a MAF, on both load paths?

    cd streamlit_app && make app-load-check

Issue #32's one acceptance criterion that pytest cannot cover, because it is a claim
about the running Streamlit app rather than about a function. Both real load call sites
are driven through `streamlit.testing.v1.AppTest`:

1. **"Open With" / auto-load** — sets `MAFIGATE_OPEN_FILE` and runs `MAFigate.py`
   unmodified, which routes through `_auto_load_file_from_path`.
2. **Upload** — runs `upload_driver.py`, which hands the page an object shaped like
   Streamlit's `UploadedFile` (a `name` plus `getvalue()` bytes) through the session key the
   sidebar's chooser writes, and then renders the data-loading page.

Each is checked for: nothing raised, no error banner, every fixture row loaded, a
filtered result set produced, **both results grids drawn and holding the report behind
them**, and `DP` arriving as `int64` — the last being what the whole swap was for. Exits
non-zero with a list of failures.

The grid assertion (`check_grid`) replaced a line that counted `st.dataframe` elements and
called a count of zero *"nothing rendered — the MAF cannot be browsed"*. The results table
is a **custom component**, not an `st.dataframe`, so that line was not merely blind but
**inverted**: zero is what a healthy app gives, and the only thing that draws an
`st.dataframe` on this page is `render_fallback_table` — the thing that runs when the grid
has *failed*. It reported the opposite of the truth on both load paths from 2026-08-17
until issue #326, which is long enough for two sessions to file the false red as a release
blocker and for issue #163 to rule the repair impossible ("expect this red indefinitely")
on the belief that `AppTest` could not see the grid. It can: issue #324 measured it.

**Four of the six scenarios take the grid assertion, and two must not.** The rule is
whether the scenario draws a report at all. *Auto-load*, *upload*, *degraded* (#39) and
*sidebar load* (#64) each do, and each now asserts on the grid; for the last two this is a
sharper claim than what stood before, because both used to read the report off **session
state** — and a report that exists in state and renders as an empty table is exactly the
shape their own bug would take. The *refusal* scenario (#38) produces no report by design
and the *arm mismatch* scenario (#135) deliberately withholds the jump to `results`, so
both legitimately draw **no grid**; asserting one there would assert the bug. Their
sections are measured, not assumed — the arm-mismatch run lands on `filter_run`.

A third scenario was added by issue #38 and asks the opposite question: on
`somatic_dot_numeric.maf`, whose depth and VAF columns hold `.`, does the app **refuse**
with a message rather than a stack trace? The refusal itself is asserted thoroughly in
`tests/test_numeric_columns.py`; what only the running app can show is that the message
survives the trip to the screen and arrives framed as a refusal rather than as a filter
error. It reaches the user through the *silent* load path — filtering happens before the
results tab exists — so the banner is stored in session state and rendered a rerender
later, which is exactly the seam where the framing used to be lost.

Both runs are seeded with `pipeline_params("somatic")` — the same expression `MAFigate.py`
itself falls back to — which is what makes this check a property of the code rather than of
the machine it runs on. (It was `DEFAULT_PARAMS.copy()` until issue #77 deleted that
constant as a second, app-unread statement of the default; the shallow `.copy()` went with
it, the app's own seed being a deep one.) `MAFigate.py` otherwise initialises
`filter_params` from `~/.mafigate/cached_parameters.json`, so the auto-load result
depended on whatever cut the developer last configured — and issue #33 found exactly that
biting: a cache left at `sample_type: germline` made the auto-load path filter the
*somatic* fixture with germline parameters, and the vendored `germline_filters` raised
`KeyError: 'InterVar'` on a file that has no germline annotation columns.

That failure is worth understanding rather than only avoiding, because it is the app
getting *better*: before the routing, the app's own germline filter guarded every source
with `if column in maf.columns`, so an arm mismatch silently skipped the entire guideline
block and reported a confident, wrong cut. It then refused; since #39 it fills the absent
columns and reports with an escalated warning. The trigger disappeared with #40 (the
parameter cache is discarded, not migrated) — the upload driver has always seeded the
defaults for the same reason.

A fourth scenario is issue #39's, and it is the mirror of the third: on a MAF missing
`CancerVar` the app must *not* refuse. It fills, reports, and escalates — and only the
running app can show the three things that matter about that. The warning has to arrive as
an **error** rather than as one more yellow box; the report has to still be there
underneath it, since being usable is the entire reason for filling; and the filled column
must **not** be in `maf_data`, because a fill the user can see is the app presenting
invented data as their own. This scenario found a real defect: every load call site filters
with `show_messages=False`, so the warnings were computed and discarded, and a user opening
an incomplete MAF got a report a quarter smaller with nothing on screen saying why.

A fifth scenario belongs to issue #64, which moved the file chooser into the sidebar so a
sample can be opened from any page. Two of its claims are only true of the whole app rather
than of the page: that a file chosen while a **different** file is open is loaded and lands on
its own report — the sidebar is built before the page beside it, and the nav radio rebuilds
itself when a page is selected from outside it, so a chosen file crosses two seams before
anything reads it — and that what the user wrote is **still there** when that report is on
screen.

The second claim was the opposite of this when #64 wrote it: notes are keyed by variant
identity rather than by file, and #64 read that as one sample's note leaking onto another's
report, so it cleared them at the load. Issue #67 settled the question the other way with the
dev — a note is a statement about the *variant*, meant to be read across the institute once
there is a server to hold it — so the surfacing is the design and the clearing rule is gone.
The check is inverted rather than deleted, because "notes survive a change of file" is exactly
as load-bearing as its opposite was, and it is the kind of behaviour a later refactor of the
shared load tail could drop without anyone noticing.
"""

import json
import os
import sys
from pathlib import Path

APP = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
MAF = (APP / "tests" / "fixtures" / "parity" / "somatic_reference.maf").resolve()
DRIVER = Path(__file__).with_name("upload_driver.py")
EXPECTED_ROWS = 82

os.environ["MAFIGATE_APP_DIR"] = str(APP)
os.environ["MAFIGATE_UPLOAD_FILE"] = str(MAF)
sys.path.insert(0, str(APP))

from streamlit.testing.v1 import AppTest  # noqa: E402

# Imported rather than spelled as a literal, the way `tests/test_row_recovery.py` and
# `tests/test_double_click_dialog.py` already read it: the assertion below asks whether the
# field *the app itself* resolves through reached the browser, so a rename in the app must
# move this check with it rather than leave it passing against a string nothing reads.
from components.variant_table import _AGGRID_INDEX_KEY  # noqa: E402

failures = []

# What `st_aggrid` calls itself in the element tree, and the two roles `components/
# results_view.py` draws (`:394` and `:444`). Grids are matched on the `aggrid_<role>_
# variants_` **prefix** of the element id and never on the whole widget key: the key ends in
# `hash(tuple(visible_columns))` (`variant_table.py:1198-1199`), which Python salts per
# process — so it differs on every run of this file, **and it can be negative**, which puts a
# second `-` in the id and breaks any `id.split("-")[-1]` parse.
GRID_COMPONENT = "st_aggrid.AgGrid.agGrid"
GRID_ROLES = {"passed": "aggrid_passed_variants_", "failed": "aggrid_failed_variants_"}


def state(app, key):
    """AppTest's session_state has no `.get()`."""
    try:
        return app.session_state[key]
    except (KeyError, AttributeError):
        return None


def grids(app):
    """The row payloads of the results grids `app` drew, keyed by role.

    `AppTest` can see the grid — `app.get("component_instance")` returns one node per
    `AgGrid` call, each carrying the whole `gridOptions` as ~100 KB of `json_args`,
    `rowData` included. This was believed impossible, and issue #163 wrote the belief down as
    a ruling; issue #324 measured it and it is false. The line that used to stand where
    `check_grid` is now counted `st.dataframe` elements, which the results table is not.

    Returns `(by_role, unmatched)`. A grid whose id matches neither role is *reported*
    rather than dropped: an unrecognised role means the two call sites moved, and a checker
    that silently found nothing to assert about would read as a pass.
    """
    by_role = {}
    unmatched = []
    for node in app.get("component_instance"):
        proto = node.proto
        if proto.component_name != GRID_COMPONENT:
            continue
        role = next(
            (name for name, prefix in GRID_ROLES.items() if prefix in proto.id), None
        )
        if role is None:
            unmatched.append(proto.id)
            continue
        payload = json.loads(proto.json_args)
        by_role[role] = payload.get("gridOptions", {}).get("rowData")
    return by_role, unmatched


def check_grid(label, app, expected):
    """Record every way the two results grids fell short of the report behind them.

    `expected` maps role to the row count that role's grid must carry, derived by the
    caller from session state — never written down, because the passing count has already
    drifted (59 at issue #309, 60 now) and a literal would be a guard with a shelf life.

    Three assertions, each naming a different failure. The second is the one this exists
    for: a grid that renders with **nothing in it** raises nothing, draws no banner, and
    is invisible to every other assertion in this file.
    """
    by_role, unmatched = grids(app)

    counts = {role: (None if rows is None else len(rows)) for role, rows in by_role.items()}
    print(f"  grids rendered:      {sorted(by_role) or 'NONE'}")
    print(f"  rows in each grid:   {counts or 'none'} (expected {expected})")

    if unmatched:
        failures.append(
            f"{label}: an st_aggrid grid rendered under an id matching neither role — "
            f"{unmatched}. The two call sites in components/results_view.py are how this "
            "check finds the Passed and Failed tables; a grid it cannot name is a grid it "
            "cannot assert anything about"
        )

    for role, want in expected.items():
        rows = by_role.get(role)

        if rows is None:
            # A fallback table is not a tolerable substitute here. `streamlit-aggrid==1.1.6`
            # is pinned exactly, `tests/test_launcher_dependencies.py` fails if any line
            # stops pinning, and both installers consume the file — so a fallback in *this*
            # interpreter means this interpreter is not the shipping one, which is a fact
            # worth reporting rather than absorbing. Same position
            # `test_an_interpreter_that_can_run_the_app_has_plotly_too` takes about plotly.
            drew = len(app.dataframe)
            instead = (
                f"{drew} st.dataframe fallback table(s) drew instead, so st_aggrid is "
                "absent or raised in this interpreter — requirements.txt pins "
                "streamlit-aggrid==1.1.6 exactly and both installers consume it, so this "
                "interpreter is not the one that ships"
                if drew
                else "and nothing drew in its place"
            )
            empty_frame = (
                " The report has no rows for this tab, so there was nothing to draw a grid "
                "from — which is a failure of the cut above it, not of the grid."
                if want == 0
                else ""
            )
            failures.append(
                f"{label}: the {role} variants grid did not render — {instead}."
                f"{empty_frame} A reader on that tab has no table to browse"
            )
            continue

        if len(rows) != want:
            blind = (
                " An empty grid raises nothing and draws no banner, so nothing else in "
                "this file can tell it from a healthy one."
                if not rows
                else ""
            )
            failures.append(
                f"{label}: the {role} variants grid holds {len(rows)} rows, not {want} — "
                f"the table on screen is not the cut the app made.{blind}"
            )

        if rows and _AGGRID_INDEX_KEY not in rows[0]:
            failures.append(
                f"{label}: the {role} variants grid's rows carry no {_AGGRID_INDEX_KEY!r}, "
                "so the row identifier never reached the browser and both dialog routes "
                "fall back to matching three columns (issue #310) — which refuses an "
                "ambiguous match, and 37 of 104 real MAFs repeat a triple, so readers get "
                "a partial row where they asked for a variant. See the pin comment beside "
                "streamlit-aggrid in requirements.txt"
            )


def check(label, app):
    """Report what one run did, and record every way it fell short."""
    exceptions = [str(e.value) for e in app.exception]
    errors = [e.value for e in app.error]
    rows = state(app, "maf_data")
    filtered = state(app, "filtered_data")

    print(f"\n--- {label} ---")
    print(f"  exceptions:          {exceptions or 'none'}")
    print(f"  st.error banners:    {errors or 'none'}")
    print(f"  rows loaded:         {None if rows is None else len(rows)}")
    print(f"  rows passing:        {None if filtered is None else len(filtered)}")
    print(f"  DP dtype:            {None if rows is None else rows['DP'].dtype}")
    print(f"  tumor_f dtype:       {None if rows is None else rows['tumor_f'].dtype}")
    print(f"  fallback tables:     {len(app.dataframe)}")

    if exceptions:
        failures.append(f"{label}: uncaught exception {exceptions}")
    blocking = [e for e in errors if "Missing Core Required Columns" in e or "Error" in e]
    if blocking:
        failures.append(f"{label}: page reported {blocking}")
    if rows is None or len(rows) != EXPECTED_ROWS:
        failures.append(
            f"{label}: expected {EXPECTED_ROWS} rows, got "
            f"{None if rows is None else len(rows)}"
        )
    if filtered is None:
        failures.append(f"{label}: filters never produced a result set")
    if filtered is not None and rows is not None:
        # Derived, and exact in both halves. The split at `page_modules/data_loading.py:1061`
        # is a boolean mask (`labelled[verdict]` / `labelled[~verdict]`), so the two frames
        # are complements of `maf_data`; `add_clinical_summary_column` only adds columns and
        # sorts, and nothing between session state and the grid caps the frame —
        # `components/results_view.py:394` and `:444` hand `prioritize_columns(...)` straight
        # to `create_data_table`, which reorders columns without dropping rows.
        #
        # `failed_data` is in session state too and could be read instead, but the
        # complement is the stronger guard: it makes the two grids account for **every row
        # of the file**, where grid-agrees-with-state would survive a split that lost rows.
        #
        # The total is `len(rows)` rather than `EXPECTED_ROWS` — the two are the same number
        # here, and the assertion above already pins `len(rows)` to `EXPECTED_ROWS`. Taking
        # it from the loaded frame keeps a short load reported once, at its cause, instead of
        # also as two grids that disagree with a constant.
        check_grid(
            label, app, {"passed": len(filtered), "failed": len(rows) - len(filtered)}
        )
    if rows is not None and str(rows["DP"].dtype) != "int64":
        failures.append(
            f"{label}: DP is {rows['DP'].dtype}, not int64 — the loader is not the "
            "pipeline's"
        )


from config.pipeline_params import pipeline_params  # noqa: E402

os.environ["MAFIGATE_OPEN_FILE"] = str(MAF)
auto = AppTest.from_file(str(APP / "MAFigate.py"), default_timeout=120)
# Seeded so the run cannot pick up this machine's parameter cache: MAFigate.py only
# consults the cache when `filter_params` is absent. See the module docstring.
auto.session_state["filter_params"] = pipeline_params("somatic")
auto.run()
check("auto-load via MAFigate.py", auto)

del os.environ["MAFIGATE_OPEN_FILE"]
upload = AppTest.from_file(str(DRIVER), default_timeout=120)
upload.run()
check("upload path", upload)


# --- issue #38: the refusal reaches the screen as a refusal -------------------
# No `check_grid` here, and that is a ruling rather than an omission: a refused MAF yields no
# report, so there is nothing for a grid to draw and demanding one would assert the bug this
# scenario exists to catch. See the module docstring's paragraph on the four-and-two split.
UNREADABLE = (APP / "tests" / "fixtures" / "parity" / "somatic_dot_numeric.maf").resolve()

os.environ["MAFIGATE_OPEN_FILE"] = str(UNREADABLE)
refused = AppTest.from_file(str(APP / "MAFigate.py"), default_timeout=120)
refused.session_state["filter_params"] = pipeline_params("somatic")
refused.run()
del os.environ["MAFIGATE_OPEN_FILE"]

exceptions = [str(e.value) for e in refused.exception]
banners = [e.value for e in refused.error]
named = [b for b in banners if "cannot be filtered" in b]

print("\n--- refusal on an unreadable MAF ---")
print(f"  exceptions:          {exceptions or 'none'}")
print(f"  refusal banner:      {(named[0][:110] + '…') if named else 'NONE'}")

if exceptions:
    failures.append(f"refusal: the app raised instead of refusing — {exceptions}")
if not named:
    failures.append(
        "refusal: no banner explained why the MAF was rejected — the user sees either "
        "nothing or a stack trace"
    )
else:
    missing = [c for c in ("t_alt_count", "t_ref_count", "tumor_f") if c not in named[0]]
    if missing:
        failures.append(f"refusal: the banner does not name {missing}")
    if "'.'" not in named[0]:
        failures.append("refusal: the banner names no offending value")
    if not named[0].startswith("🛑"):
        failures.append(
            f"refusal: framed as {named[0][:40]!r}, not as a refusal — a refused MAF must "
            "not read as a filter error, which sends the user to change a threshold"
        )

# --- issue #39: the escalated warning reaches the screen as an error ----------
# The one part of the fill that pytest cannot see. `tests/test_absent_columns.py` asserts
# which warnings are escalated; only the running app can show that an escalated one arrives
# as an **error box** rather than as the eighth yellow box in a column of them — and that the
# report is still produced and still browsable underneath it, which is the whole point of
# filling rather than refusing.
DEGRADED = Path(os.environ.get("TMPDIR", "/tmp")) / "mafigate_no_cancervar.maf"

import pandas as pd  # noqa: E402

from vendor.pipeline_utils import read_maf  # noqa: E402

with MAF.open(encoding="utf-8") as handle:
    comments = [line for line in handle if line.startswith("#")]
with DEGRADED.open("w", encoding="utf-8") as handle:
    handle.writelines(comments)
    read_maf(str(MAF)).drop(columns=["CancerVar"]).to_csv(handle, sep="\t", index=False)

# Driven through the **upload** path rather than the auto-load one. Not a preference: the
# auto-load hook fires once per session (`auto_load_checked`), and AppTest instances in one
# process share the Streamlit session, so a third auto-load run in this file silently reuses
# the first run's frame. The upload driver stubs the uploader and therefore always loads what
# it is pointed at, which is what makes it usable as the third scenario here.
os.environ["MAFIGATE_UPLOAD_FILE"] = str(DEGRADED)
degraded = AppTest.from_file(str(DRIVER), default_timeout=120)
degraded.run()
os.environ["MAFIGATE_UPLOAD_FILE"] = str(MAF)

exceptions = [str(e.value) for e in degraded.exception]
banners = [e.value for e in degraded.error]
escalated = [b for b in banners if "PATHOGENIC RETENTION DEGRADED" in b]
rows = state(degraded, "maf_data")
filtered = state(degraded, "filtered_data")

print("\n--- absent CancerVar: filled, reported, escalated ---")
print(f"  exceptions:          {exceptions or 'none'}")
print(f"  escalated banner:    {(escalated[0][:110] + '…') if escalated else 'NONE'}")
print(f"  rows loaded:         {None if rows is None else len(rows)}")
print(f"  rows passing:        {None if filtered is None else len(filtered)}")

if exceptions:
    failures.append(f"degraded: the app raised instead of filling — {exceptions}")
if rows is None or len(rows) != EXPECTED_ROWS:
    failures.append("degraded: the MAF did not load, so the fill kept nothing usable")
if filtered is None:
    failures.append(
        "degraded: no report was produced — filling exists so an incomplete MAF is still "
        "usable, and a MAF that yields nothing has not been filled, it has been refused"
    )
if not escalated:
    failures.append(
        "degraded: the pathogenic-retention warning did not reach the screen as an error — "
        "a user reading a report 25% smaller than the pipeline's is not being told why"
    )
elif "CancerVar" not in escalated[0]:
    failures.append("degraded: the escalated banner does not name the missing column")
if rows is not None and "CancerVar" in rows.columns:
    failures.append(
        "degraded: a filled CancerVar reached session state, so the fill is visible to the "
        "user as if it were their own data"
    )
# The grid assertion belongs here more than anywhere: "the report is still there underneath
# the banner" is this scenario's whole claim, and until now it was read off session state.
# A report exists in state and renders as an empty table is exactly the shape a fill that
# had gone wrong would take, and it would have passed. The counts are derived here too —
# the degraded cut is a different size from the reference one, so nothing may be reused.
if rows is not None and filtered is not None:
    check_grid("degraded", degraded, {"passed": len(filtered), "failed": len(rows) - len(filtered)})

DEGRADED.unlink(missing_ok=True)

# --- issue #64: a file chosen in the sidebar loads, and keeps what the user wrote --------
# Driven by seeding what the chooser's own callback writes — `PENDING_UPLOAD_KEY` plus the page
# to route to — because `AppTest` has no file-uploader accessor to click. Everything after
# that seam is the real app: `MAFigate.py` renders the sidebar, then the page, then fills the
# status slot, and the load happens in the page body between them.
#
# A second file, under a name of its own: the load is guarded by an upload token, so copying
# the fixture is what makes this scenario a *second* file rather than the same one arriving
# again — and a genuine change of file is what the notes have to survive.
from components.sidebar import PENDING_UPLOAD_KEY  # noqa: E402

SECOND = Path(os.environ.get("TMPDIR", "/tmp")) / "mafigate_second_sample.maf"
SECOND.write_bytes(MAF.read_bytes())

NOTED_VARIANT = "TP53:17:7577120:C>T"


class SecondSample:
    """What the sidebar's chooser hands over: a name plus bytes via getvalue()."""

    name = SECOND.name
    type = "maf"

    def getvalue(self):
        return SECOND.read_bytes()


switched = AppTest.from_file(str(APP / "MAFigate.py"), default_timeout=120)
switched.session_state["filter_params"] = pipeline_params("somatic")
switched.session_state["current_page"] = "data_loading"
switched.session_state["maf_source_name"] = "first_sample.maf"
switched.session_state["variant_notes"] = {NOTED_VARIANT: "discussed at MTB"}
switched.session_state["custom_annotations"] = {"Reviewed by": {NOTED_VARIANT: "MR"}}
switched.session_state[PENDING_UPLOAD_KEY] = SecondSample()
switched.run()

exceptions = [str(e.value) for e in switched.exception]
rows = state(switched, "maf_data")
filtered = state(switched, "filtered_data")
section = state(switched, "data_page_section")
notes = state(switched, "variant_notes")
annotations = state(switched, "custom_annotations") or {}
named = state(switched, "maf_source_name")
cleared_banner = [i.value for i in switched.info if "cleared" in i.value.lower()]
sidebar_says = [m.value for m in switched.sidebar.markdown if SECOND.name in m.value]

print("\n--- a file chosen in the sidebar ---")
print(f"  exceptions:          {exceptions or 'none'}")
print(f"  rows loaded:         {None if rows is None else len(rows)}")
print(f"  rows passing:        {None if filtered is None else len(filtered)}")
print(f"  section on screen:   {section}")
print(f"  file the app names:  {named}")
print(f"  notes carried over:  {notes}")
print(f"  annotation columns:  {sorted(annotations)}")
print(f"  values carried over: {annotations.get('Reviewed by')}")
print(f"  a 'cleared' banner:  {(cleared_banner[0][:80] + '…') if cleared_banner else 'none'}")

if exceptions:
    failures.append(f"sidebar load: the app raised — {exceptions}")
if rows is None or len(rows) != EXPECTED_ROWS:
    failures.append(
        "sidebar load: the chosen file was never read — a chooser that reaches no loader is "
        "the friction this ticket exists to remove, with a widget on top of it"
    )
if filtered is None:
    failures.append("sidebar load: the file loaded but was never filtered")
if section != "results":
    failures.append(
        f"sidebar load: landed on {section!r} rather than on the report of the file just "
        "chosen — the load happened somewhere the user is not looking"
    )
if named != SECOND.name:
    failures.append(f"sidebar load: the app still names {named!r} as the open file")
if (notes or {}).get(NOTED_VARIANT) != "discussed at MTB":
    failures.append(
        "sidebar load: the note written before this file was opened did not survive it — "
        f"got {notes!r}. A note is about the variant, not the sample in front of you, and "
        "nothing stores one, so a load that drops it destroys the user's only copy (issue #67)"
    )
if "Reviewed by" not in annotations:
    failures.append(
        "sidebar load: the user's own annotation column was dropped by the load — "
        f"columns present: {sorted(annotations)}"
    )
if annotations.get("Reviewed by", {}).get(NOTED_VARIANT) != "MR":
    failures.append(
        "sidebar load: the column survived but what the user wrote in it did not — "
        f"{annotations!r}. Keeping the column and emptying it is what issue #64 did; #67 "
        "reversed it, and a half-applied reversal is the likeliest way back to it"
    )
if cleared_banner:
    failures.append(
        "sidebar load: the app told the user their notes were cleared, over a rule that no "
        f"longer runs — {cleared_banner[0]!r}"
    )
if not sidebar_says:
    failures.append("sidebar load: the sidebar does not name the file that is now open")
# "Landed on its own report" is asserted above by reading `data_page_section`, which is a
# claim about routing rather than about what the user can see. This is the other half: the
# report the load routed to is drawn, and holds this file's rows.
if rows is not None and filtered is not None:
    check_grid(
        "sidebar load", switched, {"passed": len(filtered), "failed": len(rows) - len(filtered)}
    )

SECOND.unlink(missing_ok=True)


# --- issue #135: a MAF annotated for the other arm says so, and does not open its report --
# Same seam as the scenario above, with one thing changed: the file is the *germline*
# reference and the app is set to somatic. That is the reported session, and it is the one
# case pytest cannot fully speak for — the withheld jump is a decision the *page* makes about
# where to put the user, so it only shows up when the whole app runs and the section is read
# off the session afterwards.
#
# It is also the scenario that catches the likeliest regression: the jump is withheld in
# `_load_pending_upload`, and every other load path grants it, so a future ticket adding a
# fourth way in would get a green suite and a user dropped straight into a wrong-arm report.
GERMLINE = (APP / "tests" / "fixtures" / "parity" / "germline_reference.maf").resolve()

WRONG_ARM = Path(os.environ.get("TMPDIR", "/tmp")) / "mafigate_germline_sample.maf"
WRONG_ARM.write_bytes(GERMLINE.read_bytes())


class GermlineSample:
    """The same shape the chooser hands over, carrying a germline-annotated MAF."""

    name = WRONG_ARM.name
    type = "maf"

    def getvalue(self):
        return WRONG_ARM.read_bytes()


mismatched = AppTest.from_file(str(APP / "MAFigate.py"), default_timeout=120)
mismatched.session_state["filter_params"] = pipeline_params("somatic")
mismatched.session_state["current_page"] = "data_loading"
mismatched.session_state[PENDING_UPLOAD_KEY] = GermlineSample()
mismatched.run()

exceptions = [str(e.value) for e in mismatched.exception]
section = state(mismatched, "data_page_section")
filtered = state(mismatched, "filtered_data")
notice = [w.value for w in mismatched.warning if "and MAFigate is set to" in w.value]
switch = [b.label for b in mismatched.button if b.label.startswith("🔄 Switch")]

print("\n--- a germline MAF opened on the somatic arm ---")
print(f"  exceptions:          {exceptions or 'none'}")
print(f"  rows passing:        {None if filtered is None else len(filtered)}")
print(f"  section on screen:   {section}")
print(f"  the notice:          {(notice[0][:100] + '…') if notice else 'none'}")
print(f"  the way across:      {switch or 'none'}")

if exceptions:
    failures.append(f"arm mismatch: the app raised — {exceptions}")
if not notice:
    failures.append(
        "arm mismatch: a germline MAF was filtered on the somatic arm and the app said "
        "nothing about it — which is the whole of issue #135"
    )
if section == "results":
    failures.append(
        "arm mismatch: the app jumped to the report of a file cut on the wrong arm. The "
        "jump is withheld so the user meets the notice before the variants (issue #135)"
    )
if not switch:
    failures.append(
        "arm mismatch: the notice was drawn but offered no way across, so the user is told "
        "their report is wrong and left to find the arm control themselves"
    )
# No `check_grid` here either, for the reason above and one of its own: the withheld jump is
# the *point* of issue #135, so this run lands on `filter_run` and draws no grid on purpose.
# The `section == "results"` assertion above is what guards that, and adding a grid assertion
# here would contradict it.

WRONG_ARM.unlink(missing_ok=True)

# Issue #163 read this line as the part of the output not worth having — "the diagnostics
# are informative and the `RESULT` line is not" — and it was right while the verdict was a
# false red that no one could act on. What made it uninformative was the inverted assertion
# above it, not the line, so the fix is the assertion; what the line owed on top of that is a
# count, so a run that fails in eight places does not read like a run that fails in one.
print("\n=== RESULT ===")
if failures:
    print(f"{len(failures)} check(s) failed:")
    for failure in failures:
        print(f"FAIL {failure}")
    sys.exit(1)
print(
    "OK: both load paths loaded, filtered and rendered a MAF with pipeline dtypes; the "
    "Passed and Failed grids each drew, carrying exactly the rows the report behind them "
    "holds, each row identified as the dialog needs; an unreadable MAF was refused with a "
    "message naming every offending column and value; a MAF missing CancerVar was filled, "
    "reported, escalated to an error banner, and still browsable; a file chosen in the "
    "sidebar loaded, landed on its own drawn report, and carried the notes and annotation "
    "values written before it into it; and a germline MAF opened on the somatic arm was "
    "told so, offered the way across, and kept off its own wrong-arm report"
)
