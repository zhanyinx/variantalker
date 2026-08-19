# The MAFigate test suite

What is asserted where, and why it is asserted *there*. Written as the closing ticket of #28 — the
effort that took the app's filter from twelve divergences against `bin/filter_variants.py` to none —
because a suite that reached parity and could not say how is a suite the next person will erode.

```sh
cd streamlit_app

make test          # everything
make test-cov      # everything, with bin/ absent, against the coverage floor
make parity-reference   # the full 50-sample GERSOM reference (needs PARITY_REFERENCE_ROOT)
make app-load-check     # boot the real app and load a MAF through both load paths
```

## The governing rule

> **The harness owns every behaviour with a pipeline counterpart; the unit suite owns
> every behaviour without one.**

The rule exists because the two kinds of claim fail differently. "The app agrees with the
pipeline" can only be settled by running the pipeline, and an app-side test asserting it
against written-down numbers decays into a test of the numbers. "The app warns the user
about a filled column" has no pipeline counterpart at all — the pipeline refuses that file
— so the harness has nothing to compare and would have to invent an expectation.

There is a **third instrument**, because those two leave a gap. `tests/parity/` needs
`bin/`, and the packaged .dmg/.exe ships neither `bin/` nor the harness: on every
pipeline-less checkout and every CI job without the pipeline, the whole filter path was
unasserted. The `bin/`-free net closes it by asserting against verdicts **already frozen
in git** — `tests/fixtures/parity/MANIFEST.json` records, per named fixture row, the
verdict `bin/filter_variants.py` reached at the contract's defaults. Those verdicts came
from the pipeline, so asserting against them is a real parity claim that needs no
pipeline.

And a fourth kind, which is not about behaviour at all: the **guards**, which assert that
a thing this app restates still matches what it restates. Two of them hold copies against
their source — `vendor/` against `bin/`, and the parameter contract against
`nextflow.config` plus the ESCAT vocabulary against `resources/escat_tiering.csv` — and the
third holds the table below against the files on disk.

| module | instrument | needs `bin/` | what it owns |
|---|---|---|---|
| `parity/test_parity.py` | harness | yes | per-row PASS-set equality, column parity, the baseline, the refusal against the pipeline's raise |
| `parity/test_absent_columns.py` | harness | yes | the one deliberate deviation: that the pipeline really does refuse each filter-input column the app fills |
| `parity/test_mutation_coverage.py` | harness | yes | which divergences the fixtures still **catch** — each one re-injected into the app, 13 of 13 noticed, and the instrument shown reporting a gap when a witness is withheld |
| `parity/test_baseline_integrity.py` | net | no | what the committed baseline must say: no case diverges, no diverging row is unattributed, no attribution key outlives its probe |
| `parity/test_constructed_fixtures.py` | net | no | the **shape** of the constructed set, which no verdict can see: per-arm rows and columns, both blanking conventions, no column present-but-entirely-blank, one arm's markers per file, the two AlphaMissense columns still disagreeing, the gene lists verbatim, and the manifest's provenance claim against the committed bytes |
| `parity/test_loader_premise.py` | net | no | the harness's premise: its frame is the frame the app loads |
| `test_filter_entry_point.py` | net | no | the filter entry point against the recorded pipeline verdicts — three layers, below |
| `test_column_resolver.py` | net | no | the visible-column list opening with the vendored `compute_keep`'s, element for element |
| `test_clinical_summary.py` | net | no | `Clinical_Summary` against `variantalker_naive`, the same verdict written by the pipeline's own process and carried in the MAF — row for row on the committed germline reference; also the only pin on RENOVO's six tiers, on the `RN` glyphs they draw, and on no render table naming a filter parameter where it means a column |
| `test_vendor_compute_keep.py` | net | no | the behaviour of the vendored column computation, including the alias that leaks germline columns into a later somatic call |
| `test_absent_columns.py` | unit | no | the neutral fill, the escalation, and that a fill never reaches the returned frame |
| `test_filter_notes.py` | unit | no | which of three levels each of the filter's notes is drawn at, that a complete MAF on its own arm warns about nothing, and that the renderer cannot fall behind the levels |
| `test_stale_banners.py` | unit | no | the filter's slot describes the open file's latest run or nothing: a silent run clears the last one's notes, a refusal clears the notes it cannot be about, a change of file clears both, and every stash the slot drains is cleared by the load tail |
| `test_filter_app_extras.py` | unit | no | the app-only frequency layer and its ClinVar-pathogenic exemption |
| `test_attribution.py` | unit | no | why a particular variant is not in the report: that isolating one setting isolates only that one, that a variant outside several is told about all of them, and that the per-row and whole-report answers cannot disagree about the same row |
| `test_gene_lists.py` | unit | no | how a paste or an upload becomes the symbols we filter on |
| `test_numeric_columns.py` | unit | no | the refusal, asserted as a biconditional with the vendored functions as oracle |
| `test_param_migration.py` | unit | no | a saved parameter file still means what it meant |
| `test_parameter_adoption.py` | unit | no | no path changes the germline/somatic arm without saying so: the three controls that replace the parameters whole, the sentence they share, and that a fourth one cannot be added in silence |
| `test_partial_parameter_dict.py` | unit | no | a parameter the page cannot find is completed from the arm's contract rather than invented, and a keep-list that is present and empty is left the empty choice it is |
| `test_one_resolver_on_the_filter_path.py` | unit | no | the data page completes its parameters against their arm's contract before anything filters, so the engine's neutral literals stop being a second opinion about an absent key — and the one key completion cannot supply, `sample_type`, is refused rather than defaulted to somatic. Carries the premise as its own test: the partial set and the completed one must still filter differently |
| `test_incomplete_cache_is_retired.py` | unit | no | a stamped cache that is not a complete set for its own arm is one this app could only have written by inventing values, so it is set aside and the banner says why — while a *complete* dict with every guideline source emptied on purpose is left untouched on both arms. Emptiness is never consulted. Carries the live fire for its own render harness: a completion made to miss `skip_civic` must be seen |
| `test_app_defaults.py` | unit | no | the app opening at parity, with the guideline sentinel gone |
| `test_data_page_sections.py` | unit | no | the Load/Results switch: a fixed order, and a jump to Results that fires once per file |
| `test_download_contents.py` | unit | no | every download carrying the notes and annotation columns the user typed — the only place they survive |
| `test_missing_column_warning.py` | unit | no | the missing-column sentence said once per render at page level, in every state including *Show all columns*, and drained rather than left behind |
| `test_results_summary.py` | unit | no | the Summary tab inventories no columns, a complete MAF is reported complete, and a column the file really lacks is still named against the arm that emits it |
| `test_overview_legend.py` | unit | no | the Pathogenicity Overview drawing one circle per source this arm emits and this file carries, and its key naming exactly those |
| `test_arm_mismatch.py` | unit | no | the app saying so when the loaded MAF and the selected arm disagree: the detector set held to a narrowing of the filter contract, silence where the file cannot be placed, a consequence sentence chosen by this file's own counts, and the jump to Results withheld |
| `test_file_chooser.py` | unit | no | opening a MAF from the sidebar, and that a change of file leaves what the user wrote in place |
| `test_utils.py` | unit | no | the MAF loader, now a delegation to the vendored reader |
| `test_missing_values.py` | unit | no | the app's two answers to *this row does not say*: that they differ by exactly the annotator's sentinels, that the display one catches every spelling this pandas calls missing, and that `-` is a real value in every column |
| `test_run_recap.py` | unit | no | the recap of the filters above the report: that it names settings in the interface's own words, reads their values off the object the filter ran on rather than the dict, and describes the run rather than the controls as they stand now |
| `test_components.py` | unit | no | UI components |
| `test_config.py` | unit | no | configuration and constants |
| `test_app_identity.py` | unit | no | the app's one self-description, the claims it may not make about itself, and the single surface that prints its version — and, since issue #263, its channel and its build, all three swept the same way and the sweep shown catching a planted copy |
| `test_build_identity.py` | guard | no | what a build knows about itself: the writer run isolated for both installer channels, its channel list held equal to the app's in both directions, the stamp's location derived from the module the app imports, a wrong channel refused where a build can still be stopped, an unfamiliar one passed through rather than quietly rewritten, a damaged stamp falling back field by field, one reader so the source-checkout default cannot be bypassed, and both build scripts stamping *before* they package — plus the deliberate asymmetry with `build/version.py`: a build machine with no `git` gets a warning and an installer, not a failure. The commit is proved read from the tree it was pointed at, by stamping a throwaway one-commit repository rather than this checkout, and `-dirty` is held to a tracked edit and to *not* firing on an untracked file. The two copies a builder reads are pinned too: the documented `--channel` commands, and `make clean`'s removal of the stamp a local build leaves behind — because a build stamp that outlives its build makes a source run call itself an installer |
| `test_help_claims.py` | unit | no | the Help page's checkable claims: the upload limit, the export formats, the column and vocabulary lists it must derive rather than transcribe, and the ESCAT and CIViC definitions it must render, cite, and not keep a second copy of |
| `test_component_seams.py` | guard | no | the shape of `components/` against the account `__init__.py` gives of it: one-way imports, no re-exports, a module table that still names what exists |
| `test_vendor_drift.py` | guard | yes | `vendor/` against `bin/`, by hash and per-symbol AST equality |
| `test_param_contract.py` | guard | yes | `config/pipeline_params.py` against `nextflow.config`, and the ESCAT levels the control offers against the tiers `resources/escat_tiering.csv` assigns |
| `test_messages_reach_the_user.py` | unit | no | the two sentences issue #140 moved, read back off a running app: the note confirmation drawn and drained on the run after it is parked, and a failing page reported in place under the name the sidebar gives it |
| `test_discarded_frames.py` | guard | no | no message is drawn before an `st.rerun` in the same block, app-wide — the frame is thrown away, so the sentence reaches nobody |
| `test_suite_organisation.py` | guard | no | this table: every module has a home, and no home names a module that is gone |
| `test_parameters_page_stays_put.py` | unit | no | changing any control on the parameters page leaves you on it: issue #277's reproducer, now issue #283's regression guard with its strict `xfail` deleted. A filter multiselect and a number input, so the defect is held to be neither multiselect- nor page-specific; a control case with no file held, which is what pins the cause on the sidebar chooser's `on_change` rather than on anything the page itself does; and the mechanism counted — one file chosen, one hand-over, no matter how many interactions follow. Plus the two things the fix had to keep while removing three reruns: the nav radio still naming the page the body renders after a programmatic navigation, which is what the deleted reconcile rerun was for, and the route back to your file still reaching it in one click from the parked rerun that replaced its own. The browser's held file is supplied at the `WidgetStates` layer because `AppTest` has no `file_uploader` |
| `test_chromosome_spelling.py` | unit | no | where the `chr` prefix comes from and that no surface adds a second one: the rule settled once on load, four render sites that add nothing to it, and the notice shown only on the file it changed |
| `test_contaminated_columns.py` | unit | no | a predictor column holding the chromosome instead of a score never reaches the ClinGen scale: the whole-column rule and the measured coincidence it must stay quiet about, a warned row drawn but unscored and out of the tally, and the verdict taken after the chromosome spelling is settled |
| `test_double_click_dialog.py` | unit | no | the variant dialog opened by double-clicking a row: one opening per click and none on the reruns after it, the same row re-openable, the stamp that makes a second click distinguishable, and no machinery of it reaching the dialog or the download |
| `test_fallback_dialog.py` | unit | no | the route into the variant dialog where `streamlit-aggrid` is absent: a selector that chooses and a button that opens, neither offering nor opening anything on a render nobody touched, and one spelling of the button's label across both table paths |
| `test_chromosome_rule_contract.py` | guard | yes | `config/chromosome_spelling.py`'s rule against the one `bin/add_guidelines_escat_to_funcotator.py` applies — the decision and the value written, compared behaviourally because the two are written in different idioms |
| `test_sidebar_doors.py` | guard | no | the sidebar offers no second door onto a page its nav radio already lists: a button there must reach what the radio cannot, checked over all three syntaxes that could route one |
| `test_sidebar_reruns_last.py` | guard | no | the sidebar takes its reruns after the file chooser has been drawn, or not at all: `st.rerun` is not a premature stop, so a run that ends earlier prunes the chooser's state and its `on_change` re-fires on the user's next interaction anywhere in the app (issues #277, #283). Two rules — every rerun in the module inside `render_into_status_slot`, and inside it below the call that draws the column — with non-emptiness asserted by the first so the second is never quantified over nothing |
| `test_installer_version.py` | guard | no | that `APP_VERSION` is the only place a MAFigate version is written: every installer input held to carrying *no* version literal and to naming its derivation — the Inno script's two fields, the DMG script's filename, the app bundle's two `Info.plist` keys, the Makefile, and the build prose — plus the `mafigate-v` tag namespace, kept out of the Nextflow pipeline's live `v1.x` line. The derivation is proved by *moving the constant*: `build/version.py` run against a throwaway tree that says `9.9.9` must say `9.9.9`, since a tool that had gone back to a literal would agree with the real constant and only disagree with that one |
| `test_launcher_dependencies.py` | guard | no | installing and launching land in the same interpreter: one variable through both scripts and the Makefile, pip/streamlit/pytest reached as `-m` modules rather than as PATH entries of their own, a launcher that checks before it boots, and a checker that really reads `requirements.txt` — plus what that file may declare: every line pinned to one exact version with the reason for the number beside it, the development extras still only comments, and every sweep in the class proving its own parse against a named anchor before it asserts, since a sweep that read nothing would run zero iterations and pass. Also the pin pattern itself, tabled in both directions — a floor or a marker must not read as pinned, and a legal exact pin with extras or spaces must not read as unpinned — and the Python floor the pins impose (`numpy==2.3.1` needs 3.11), held against the prerequisite the README promises. Since #258 that interpreter is a virtual environment in the checkout, so the resolution is asserted by *running* it against stub environments — the venv preferred, the Windows `Scripts/` layout found, an interpreter that cannot be run refused, `MAFIGATE_PYTHON` short-circuiting — plus the rule that nothing but `mafigate_python.sh` may spell the path, which the Makefile is held to as well |
| `test_delivery_channels_copy.py` | guard | no | what the README may promise about each delivery channel, tied to whether the tree can keep the promise yet: no link to a release artifact while `build/RELEASES.md` records no published release, and — the other direction, added when #264 re-keyed this switch off the mere *presence* of a release workflow — no copy still calling the installers unavailable once one is recorded, since the workflow's own arrival used to make both tests skip rather than flip. Both halves skip in one of the two states, so the flip is additionally **driven** over seeded copy — each rule evaluated against the state it must accept and the state it must refuse — because a dormant assertion is one nobody notices has stopped working, and this ticket's burden was that the guard still refuses both wrong states. Plus the ledger's parse, which must state one of its two states out loud, may not name a version ahead of `APP_VERSION`, and does not read a worked example inside an HTML comment as an event; the data-posture claim carried in one canonical table and qualified for exactly as long as no config turns Streamlit's usage reporting off; and a contributing note stating the one-way inbound flow and its manual credit |
| `test_release_page.py` | guard | no | what the release page must say, asserted on the page a reader gets rather than on the file on disk — the two placeholders CI fills mean the source is a document with a hole where the downloads go, so the assembly is run through `test_release_workflow`'s own step, imported rather than re-implemented for the reason #247 recorded: a check that copies a script's loop goes green on input the real one refuses. Four decisions of the delivery matrix are held to the page — the MAF prerequisite stated **above** the downloads (below them it is a note about why the download was wasted) and naming all three extensions the chooser accepts; the three failure modes, including the two that do not present as failures, an install that reports success and then will not start and a process the operating system stops with no message at all; #225's memory rule with both of its terms, so it can be applied to the file the reader has, plus what the largest MAFs cost in practice; and no update mechanism promised, with the absence said out loud, since it would be the app's first outbound call and would falsify the posture claim two paragraphs above it. The pinned data-posture sentence and the tree fact that earns it are both **imported** from `test_delivery_channels_copy`, so this page cannot keep promising something after that switch has flipped the other way. Deliberately not here: the first-open wording, which `test_unsigned_artifact_copy` sweeps to one home — the page carries a placeholder and would fail that sibling the moment it grew a paraphrase — and version literals, already swept by prefix out of `streamlit_app/build/`. Plus the assembly's own refusals, driven over seeded pages: a placeholder lost, duplicated or renamed, the failure that would otherwise ship a release page with no files named on it and nothing about it looking wrong |
| `test_release_workflow.py` | guard | no | the workflow that builds both installers on a `mafigate-v*` tag, asserted by **running** the two steps that decide anything rather than by reading its YAML: the version check accepts the tag `APP_VERSION` derives, marks the same version with a suffix a rehearsal pre-release, and refuses five tags naming a version this tree does not hold; the one-commit check refuses two commits and refuses an *empty* commit report, which three equal empty strings would otherwise pass as agreement; and the notes step, run over a throwaway tree of fake artifacts, must name both files with a real digest each and quote `build/OPENING_MAFIGATE.md` **whole**, must refuse anything but one artifact of each kind *and say which was missing*, and must refuse a hashing tool that either produces nothing or produces something that is not a digest — the second being the one a non-zero exit does not catch, and the pair a regression this file found by running the step, since a command substitution failing inside a `printf` argument leaves an empty cell and exits 0, so `set -e` never sees it and the release page offers two blank checksums. The step that touches the outside world is run too, against a `gh` that records rather than acts: both artifacts attached, `--draft` always, `--prerelease` only for a rehearsal tag, an existing *draft* edited rather than re-created, and an already-**published** release refused outright — because `gh release edit --draft` would retract the page a clinician is downloading from and report success. Left as text, their failure modes being unprovokable here: the `< NUL` that stops `pause` hanging the Windows job for six hours, and the unversioned ISCC compile that proves `installer.iss`'s `#error` (the one assertion `test_installer_version` says it cannot make). Plus that only a tag push can start it, which is the runner-cost answer as much as the trigger |
| `test_telemetry_config.py` | guard | no | Streamlit's usage reporting is off on every one of the four launch routes, answered by Streamlit's own resolver in a subprocess standing where each route stands — with the environment scrubbed, so the tracked `.streamlit/config.toml` is the only thing that can turn it off; both installers' copy lists stood up as trees so a config no build ships is a red test; the default proved to be *on* and the check proved able to name an unreachable route; a completeness sweep so a fifth launcher cannot arrive unmeasured, and so a deleted one is reported rather than silently unmeasured; and the three neighbouring surfaces that promise silence — the installer's welcome screen, the Help page's *queries no external service*, the parameter cache's privacy line — each held to what the config makes true. Its sibling is `test_no_outbound_calls.py`: that one is about the app's own code, this one about the framework it runs inside |
| `test_no_outbound_calls.py` | guard | no | that the shipped app makes no outbound network call — the property issue #229's *no update check, ever* rests on, asserted rather than measured once, in all three languages that reach a user's machine: across the 52 tracked Python files, no network client imported, dynamically named, called, remote-URL-fed to a pandas reader, or shelled out to; across the shipped shell entry points — **derived** from `git ls-files`, not listed — no command reaching a host that is not this one, judged per command and against that command's own URLs, because a loopback address anywhere on a line was found vouching for every other command on it; and no network API in the native Swift launcher, which is what the double-click runs. Scope stated in its own docstring and held there: `tests/` and `update_db/` are out and `build/` deliberately is not, the first-run `pip install -r requirements.txt` is the one allowed fetch, and each rule is seeded and watched refusing, with a matching case per rule that must *not* trip it |
| `test_dead_session_keys.py` | guard | no | no branch waits on a session key nothing ever sets: every key the app reads held against the writes that put a value in — an assignment, a widget's `key=`, or mutation of the container it holds — and no key written and never read back |
| `test_predictor_context.py` | unit | no | the raw predictor scores behind InterVar's PP3/BP4 and CancerVar's CBP10: every letter in the measured vocabulary read the way its own tool reads it, the two tools' opposite treatment of a missing input said out loud, `GERP++_RS` read at `>= 2` on one arm and `> 2` on the other, no near-neighbour substituted for an absent column, and a necessary-condition test that names the 15 files whose cells cannot account for the verdict |
| `test_public_export.py` | guard | no | the tree `tools/export_public.py` publishes, held against what the deny-list and the scan gate claim about it — real exports into temporary repositories, asserting the notes and the private repository's name are absent, `.github/` and `bin/` present, the superseded public-only paths deleted, one commit naming the private SHA with no private history crossing; plus a seeded violation per scan rule to watch the gate refuse, and a case per rule that must *not* trip it. Also the one claim here about the **real** tree rather than a seeded one: a sweep of `git ls-files` minus the deny-list through the gate's own `scan_file`, so a case identifier committed on an ordinary day fails on that day instead of at release — with its own two vacuity checks, since a claim of zero is satisfied by a sweep that reads nothing, and a planted hit of *every* rule, since a matcher that had stopped seeing one rule would also report zero |
| `test_column_spelling.py` | guard | no | one rule for every column the variant panel reads by name: a name `make.names` would rewrite must reach the row through `config.columns.spelled_in` and never by equality — held both against the panel's own AST, which catches a literal lookup, and against a wholesale-mangled header rendered through the whole panel, which catches a table-driven one the AST cannot see; plus the resolver itself, including the padding 93 of 96 real files arrive with and the tie it must keep refusing |
| `test_predictor_cutoff_contract.py` | guard | yes | the five thresholds `components/predictor_context.py` prints against the assignments `resources/{CancerVar,InterVar}` actually make — read from each function's own AST, so a number in a comment or a neighbouring function cannot answer for it, plus the two comparison operators that differ where the numbers are equal |
| `test_clinical_badges.py` | guard | no | map #199's five re-homed annotations, each drawn exactly once — measured off a recording `st` planted in *every* module the panel draws through, since a text guard over `variant_detail.py` is satisfied by that module's own account of what moved where — plus the four vocabularies the row now carries: ClinVar's star hierarchy total over the corpus and over the level no MAF here holds, with an unmapped status falling back to the full string and an earned zero rendering unlike an absent one; RENOVO's badge degrading to class-only on the 6,597 rows with no score rather than to `(nan)`; ESCAT's three tiering-table vintages, the 2 files where no tumour type may be invented and no gloss borrowed, and a nine-pair cell that groups to three badges keeping every tumour type; and the empty state's invariant — exactly two captions on every empty row, one per member state, never naming RENOVO |
| `test_alphamissense.py` | unit | no | one predictor's own section: the column it reads (`am_pathogenicity`, not the dbNSFP score that disagrees with it on 346 rows), the band named at AlphaMissense's own cutoffs rather than echoing a class label 100 of 139 files overstate, a provenance sentence *derived* from what the classifiers actually read rather than transcribed — asserted both ways, including the case where a classifier does read it — an ACMG strength claimed nowhere, no row inside the ClinGen table, and the position that was the whole of the decision |
| `test_reference_scales.py` | unit | no | the ClinGen SVI PP3/BP4 section says what it is and does not claim to be the verdict: the disjointness the label's *"read none of these predictors"* rests on, asserted against what `predictor_context` actually reads rather than transcribed, the germline-calibration sentence that appears only on the somatic arm, no tool named on a row that carries neither, a collapsed expander whose label says whether it is empty, an empty state that names the variant instead of claiming no data, and the section's position *after* both evidence sections asserted as order |
| `test_cbp_evidence.py` | unit | no | CancerVar's twelve AMP/ASCO/CAP criteria and the somatic evidence section: a parser that refuses a vector it cannot read whole rather than moving the variant down CancerVar's scale, a table of what fired in *either* direction with sentences that know the score, two kinds of zero kept apart, a tier read from the string and never banded from the sum, and both real spellings of the evidence column driven off fixtures on disk — plus, since issue #210, the captions above that section: a file with no `CancerVar` column is promised a tier only when one is actually drawn below (false on 73 real rows before), a vector summing to zero still counts as having a tier, and the caption never claims a germline arm on a file whose somatic tier draws beneath it |
| `test_cancervar_markers.py` | guard | no | the markers behind CBP1/CBP2/CBP3, resolved as 0-based line offsets into the **vendored** copy of CancerVar's marker table — so the expected drugs and levels are read off real indices in that file and a re-vendored vintage turns them red rather than silently renaming drugs. Plus the tumour-type split the module computes rather than reads (only 10.6% of real therapy markers match the sample's own tissue), three named refusals that each keep the cited count, and the disclosure's HTML — including the wiring, which mutation testing showed the direct-call tests did not cover |
| `test_public_repo_name.py` | guard | no | which repository this tree names out loud: the private repo's name swept as a literal to a count of **zero** across every tracked file outside `docs/wayfinder/` — no test-tree exemption, and its own two vacuity checks, since a claim of zero is satisfied by a sweep that reads nothing — plus the two clickable menu items, Get Help reaching the app's own README rather than the pipeline's and Report a bug the public issue tracker |
| `test_unsigned_artifact_copy.py` | guard | no | what we tell a recipient whose computer has just refused to open MAFigate, held to **one** home: the macOS button, the Windows button and the Terminal fallback command each swept to exactly one tracked file, the retired second `xattr` spelling swept to zero, and both alerts quoted rather than only their fixes — the macOS one against Apple's own `Quarantine.loctable` string. Plus the ordering, read off the note's structure, since a fallback printed above the primary route is the headline whatever it is called; the DMG build script staging the note into the mounted window with a position of its own; and the publisher URL Windows renders in Add/Remove Programs, read off its `[Setup]` key so a URL moved into a comment cannot answer for it, with the placeholder token forbidden tree-wide as the covering rule issue #234 could not write without exemptions. Its vacuity check has two halves on purpose: the phrases swept to *one* witness themselves, but a literal swept to *zero* and assembled from halves cannot be proved by planting it — searching for the value you just planted succeeds whatever the value is — so those two are held against an independently encoded spelling as well. Three further claims about the note itself: that it stays readable as **plain text**, since the DMG ships this same file renamed and any `#` or `>` reaches the reader raw; that the route macOS 15 removed stays *second*, which no ordering check over the Terminal fallback can see; and that what each alert **says** is pinned as tightly as what to click, the quotations being the half a release page is likeliest to copy |
| `test_windows_installer.py` | guard | no | the two claims the Windows .exe makes that no import can reach: that it brings its own interpreter — one `Source:` line shipping it, a launcher pointing at it with no search of the environment and no python.org error left to reach, a setup hook with the soft *"continue anyway?"* prompt gone, the release pinned in **one** file both builds *read* so the .dmg and the .exe cannot ship different CPythons, and the two builds' `tar --exclude` lists held equal in both directions so one platform cannot ship the test suite and `tkinter` the other trims — and that uninstalling reclaims only the venv and the logs, never the `~/.mafigate` parent that also holds the parameter cache. Plus the two ways a `.bat` reads correctly and behaves otherwise, asserted over both of this repo's: `%errorlevel%` inside a parenthesised block, which `cmd` substitutes before the block runs, and an unquoted `%VAR%` echoed inside one, which a path containing a parenthesis closes early. Plus every path the **compiler** reads — `[Setup]`'s `LicenseFile`/`SetupIconFile` and each unflagged `[Files]` `Source:` — resolved against the tree, which is the one class of defect here that only a compiler had ever been able to find: `mafigate-v1.0.0-rc2`, the first Windows compile in this project's history to get past parsing, died on a `LicenseFile` two levels up where the licence is three, wrong since it was written and invisible to every test that read this script as text. A `[Setup]` directive has no `skipifsourcedoesntexist`, so tolerance is spelled `#ifexist` and the parser understands it as such; build products are exempted by asking **git** whether the tree ignores them, in both spellings, since a directory rule written with a trailing slash does not match the path without one and the staged interpreter would otherwise read as missing |

The **needs `bin/`** column is not a description, it is a constraint, and the five `yes` rows are
the complete list — `test_suite_organisation.py` names them and refuses a sixth without a deliberate
edit. It keys on the module's *name*, not on how it is filed, so nothing can exempt itself by
choosing a label. `make test-cov` then runs the whole suite in a tree with no `bin/` at all, so a
`net` or `unit` module that quietly acquired a pipeline dependency would either fail there or stop
asserting anything — which is #24's failure exactly: a module-level `skipif` removed the entire
parity suite from every pipeline-less checkout, and nothing said so.

The check reads `pytest.skip(...)` as well as `skipif(..., reason=...)`, and that is not
pedantry: this suite's own preferred shape puts the gate in a *fixture*, where a call site
cannot forget it, and a checker blind to that shape would wave through the better pattern
while catching only the worse one.

## The house pattern

> **Where the app must agree with the pipeline about a list, derive it or guard it —
> never copy it.**

Every drifted `KEEP` copy in this codebase is what happens when this rule is broken. There
were three of them — `bin/utils.py:KEEP`, the app's utility-module copy, and a third in
the column config — reported as 40 columns against 40 while their *contents* differed by a
substitution. Three instances of the rule, kept honest three different ways:

1. **The vendored constant behind the frequency exemption.** The app's frequency layer
   exempts ClinVar-pathogenic calls. It does not hold a pathogenic vocabulary to do it:
   `filters/variant_filters.py` imports `CLINVAR_PATHO` and `has_clinvar_term` from
   `vendor.pipeline_utils` — the pipeline's own constant and the pipeline's own test over
   it. `test_filter_app_extras.py::test_the_module_holds_no_second_pathogenic_list`
   asserts there is no second list to drift.
2. **The AST-guarded column computation.** `compute_keep` decides which columns the report
   carries. It is *vendored*, not restated, and `vendor/_sync.py` compares it to `bin/`
   symbol by symbol on the parsed tree — so a whitespace edit passes and a changed branch
   does not. It gates both installer builds, because a `.dmg` cut from a drifted copy
   carries no `bin/` to be caught against.
3. **The AST-derived validation column set.** Which columns the app refuses a MAF over is
   computed at import by parsing the vendored source for columns reaching an arithmetic or
   ordering operator (`filters/numeric_columns.py`). Nobody can make the app fussier or
   laxer than the pipeline by editing a list, because there is no list.

Where a value genuinely cannot be derived, it is **guarded behaviourally**, which is
stronger than comparing spellings because it asserts the consequence:

- `NO_GENE_FILE` is a bare string inside two vendored function bodies, not a named
  constant, so it cannot be imported. `test_filter_entry_point.py` instead asserts that
  the mask the vendored clause returns for the sentinel is identical to the mask for a
  gene list naming every symbol in the MAF. If `bin/` changes the sentinel, that fails.
- `config/pipeline_params.py` restates `nextflow.config`'s `params` block, because the
  packaged app ships no `nextflow.config` to read. `test_param_contract.py` and
  `make check-params` hold it against the real file — deliberately *not* a build gate: a
  moved pipeline default is a decision to make, not a broken build.

## Two failure modes this suite is designed against

**Counts that coincide.** Column parity read 40-against-40 while the contents differed by
one substitution. `somatic_depth_500` sat at pipeline 10, app 10 — equal totals — while one
row passed only in the pipeline and a *different* row passed only in the app. And the
fixture cell that was supposed to witness divergence #6 was silenced by divergence #2
pointing the other way, so a cell reported coverage it did not have for the whole effort.
**Compare sets and per-row values, never totals alone.**

**Numbers that are not self-describing.** The somatic reference baseline is **408** on the
criteria path and **411** on the union — `(common & specific)` against
`(common & specific) | filter_patho`. Both are true; neither is self-describing; calling
either "the somatic baseline" is how a real discrepancy stayed open across three tickets.
**Quote the cell, not the number.** `Diagnostics.passed` is always the union;
`Diagnostics.criteria_path` is always the other one; the two cells are held apart by name
everywhere they appear. `test_filter_entry_point.py` pins which fixtures can actually tell
them apart — the two reference MAFs, where the pathogenic rescue admits a row the criteria
path did not — and which cannot, so nobody writes a cell-sensitive assertion on a fixture
that would pass either way.

## The `bin/`-free net, layer by layer

Three layers, weakest **last**, and the order is the design:

1. **Per-row verdict equality** on the fixtures' named rows, addressed *positionally*.
   This is the strong layer, and the addressing is what makes S1's "every row returned, in
   input order, index preserved" a load-bearing assertion rather than a claim: if the entry
   point ever reorders, drops or re-indexes a row, the names stop lining up.
2. **Oracle-free contract checks** on the real rows — every row returned, order and index
   preserved, both verdict columns present and populated, the reason unable to disagree
   with the verdict, the four diagnostic cells partitioning the frame, and the visible
   column list opening with the pipeline's own (computed from the vendored `compute_keep`,
   since there is no header to read without `bin/`).
3. **Aggregate PASS counts**, which are **deliberately the weakest layer and say so in the
   file**, for the coincidence reason above. They are a coarse check that the union is not
   wildly off. If you are tempted to assert a new number there, assert a row in layer 1.

## Diagnostics are asserted as a partition

The four cells — `criteria_only`, `rescue_only`, `both`, `rejected` — are asserted to sum
to the row count and to agree cell-for-cell with the reason column, **not** against fixed
numbers. Fixed cell counts would be a second baseline to keep in step, and the property
that makes the channel trustworthy is the arithmetic. What this replaced is why: 28
diagnostic strings, of which the classification count read 68,593 against a true 20,386.
`tests/parity/reference.py` runs the stronger per-variant form over all 1,100 reference
measurements.

## The ratchet, and why it is gone

While divergences were open, a diverging case carried `xfail(strict=True)` so the suite
went **red the moment the case started passing** — whoever fixed it re-recorded
`baseline.json` and dropped the marker in the same commit, so progress could be neither
silently lost nor silently claimed. It was proven in both directions on a real one-line
fix.

With the last divergence closed the marker inverts: one that can still be handed out is
one a *regression* can be absorbed into. Issue #35 made this argument for column parity
when it became structural; issue #41 applied it to the rows and removed the machinery, so
a case that stops passing now fails. **Green here is not green-by-skipping**, and the
claim rests on checks rather than on this sentence: `parity/test_baseline_integrity.py`
asserts off the committed baseline, with no `bin/`, that no case diverges and that no
diverging row is unattributed.

Nor is it green-by-emptiness. Agreement is the cheapest property a fixture set can hold —
rows neither side does anything with agree trivially — so `parity/test_mutation_coverage.py`
re-injects each of the map's divergences into the app and requires some case to notice.
That is the check issue #242 built to replace `test_attribution_coverage.py`, whose
soundness assertion had inverted into asserting the absence of what it was written to
detect.

## Coverage

The floor is a ratchet against regression, not an aspiration, and it is **measured with
`bin/` absent** — see the Makefile, which builds a pipeline-less tree to run it in. That
is not a detail: the parity module calls the app's entry point in-process, so the same
command reported one number with `bin/` present and another without, and a floor that
moves with the checkout is not a floor. Raise `COV_MIN` as real coverage improves.

## No timing tests

There are none, and that is deliberate. The two that existed went with
`test_integration.py` in issue #33: one bounded a 1000-row filter at 5 seconds — roughly
1000x slack, green on any machine and on any implementation, including a quadratic one —
and the other built 100 rows and then asserted nothing about cost at all. What replaced
them asserts the costs that actually exist, deterministically: `frame_for_masks`
returning the caller's own object on the happy path, measured with `tracemalloc` against a
400,000-row frame, and the gene adapter's temp file existing inside the context manager
and gone after it. The `performance` marker is not declared in `pytest.ini` any more, so
the next timing test cannot land under a name that reads as sanctioned.

## Privacy

`tests/fixtures/parity/` **no longer carries patient data.** It did until issue #246: all
seven MAFs were cut from a 50-sample paired clinical run and pseudonymised, which is
weaker than it sounds — the five not named "reference" were real rows with two to four
cells edited, and the one that rewrote `Start_Position` to a round number still carried
the real coordinate verbatim in `Genome_Change`. They are now constructed end to end by
`fixtures/parity/build_parity_fixtures.py`, which needs nothing mounted; `--check`
rebuilds and compares byte for byte. `fixtures/cancervar/` is likewise emitted by its own
builder. Both installers exclude `tests/` anyway, asserted by
`parity/test_parity.py::test_harness_is_excluded_from_packaged_builds`.

Publication goes through `tools/export_public.py`, whose scan gate refuses to export any
MAF that its fixture set's `MANIFEST.json` does not record as `generator-constructed` with
a matching checksum. That is a claim about where a file came from rather than what is in
it, and it is deliberately not a content pattern: five of the seven files this set
replaced passed every identifier pattern the gate can apply while carrying real calls, so
by content the bytes that must not travel are indistinguishable from the bytes that must.
`test_public_export.py` is the proof, and
`parity/test_constructed_fixtures.py` holds the manifest's claim against the committed
bytes.
