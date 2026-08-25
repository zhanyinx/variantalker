# Published MAFigate releases

The tree's record of which releases have actually been **published** — not which ones could
be built, and not which ones are drafted. One line per release, newest first.

This file exists because a promise in the README is only true once a download exists, and
something checkable had to answer *does it?*. Issue 264 re-keyed that switch onto this
file: the guard in `tests/test_delivery_channels_copy.py` reads the entries below, and while
there are none the README may not link to a release artifact and must say the installers are
not available yet. Adding an entry here flips it, in both directions.

Before this, the switch keyed on whether a release workflow existed in `.github/workflows`.
That is a **capability**, not an event, and it was wrong by exactly the window between the
workflow landing (issue 264) and the first release being cut (issue 265): the workflow's own arrival
would have made the guard skip — not fail, *skip* — leaving both wrong states unguarded in
the one window where a README could link to a 404. A record of the event cannot be wrong
that way. The cost is that the record is written by hand, which is the same cost every
other prose guard in this suite carries.

## Releases

- mafigate-v1.0.0 — published 2026-08-25 — commit b086d27

     The first MAFigate release. The commit is the public repository's, which is what the
     tag names and what the release page prints; the private tree it was exported from is
     recorded in that commit's own message.

     Three rehearsal tags were needed to get a build at all, each finding a defect no test
     here could see, because all three lived in files no Python imports and one of them only
     existed for a compiler: a check that passed and exited non-zero anyway, a `LicenseFile`
     two levels up where the licence is three, and `{userprofile}`, which Inno Setup does not
     define. Each is guarded now. See `BUILD_INSTRUCTIONS.md` on rehearsing, and rehearse.

     It was cut four times, and the reason is the more useful thing recorded on this line.
     Every draft after the first was replaced because the tree moved on underneath it while
     the draft sat there naming its own commit correctly and reading exactly as a fresh one
     reads — first under twenty-two commits including the fix for a dialog that opened the
     *wrong variant*, then under seventy-four more. `make release-preflight` exists because
     of the first of those, and step 4 is not optional.

     It did not catch the second. The check compared the tag against the tip of what is
     published, and those are the same commit whenever nobody has exported since the tag was
     cut — the ordinary state. Issue 405 moved it one level out, onto the private commit the
     published tree says it was exported from. A release cut from a tree with a pending
     export is now refused.

     And publishing was not the end of it: this tag went live twice. The first publish,
     2026-08-24, lasted under an hour — the first external mac users hit a second launch
     that died, because the venv had pinned the bundled interpreter at the ephemeral path
     the app held on first launch, and an app that had moved left that pin dangling
     (issue 411). Every rehearsal, human and scripted alike, had only ever exercised first
     launches — and first launches all work. Every download of the retracted build was a
     test download, so the release was re-drafted, the tag deleted, and the number reused
     rather than bumped (issue 415). The relaunch is guarded twice now: the
     launcher-contract workflow runs the launcher's venv mechanism against a moved
     interpreter on every push (issue 414), and the release route gained a manual
     move-then-relaunch rehearsal step — launch, quit, move the `.app` by hand, launch
     again (issue 418).

<!-- One list item per published release, newest first, in this form:

- mafigate-v<version> — published <YYYY-MM-DD> — commit <short sha>

     The tag is what the guard reads, so it has to open the line. Nothing else on the line
     is parsed, and nothing else is required; the date and the commit are there because
     the questions asked of a shipped release are always "when" and "from what".

     Do NOT write an artifact filename here with the version spelled into it —
     `MAFigate-<version>-macOS-universal.dmg` and no other way. This file sits under
     `build/`, which `tests/test_installer_version.py` sweeps for version literals, and the
     release page is where the real filenames are already written for you.

     The release itself is not published by CI. The workflow attaches both artifacts to a
     **draft** and stops; a human authors the page, presses Publish, and then records it
     here. So this file cannot say something CI did behind anyone's back. -->

## What to do when you publish one

First the tag. `cd streamlit_app && make release-tag` prints it, derived from `APP_VERSION`;
push it, and the workflow builds the DMG on macOS and the installer on Windows from the one
commit the tag names, attaching both to a **draft** release.

Then these seven steps — the same seven, in the same order and under the same numbers, as
`build/BUILD_INSTRUCTIONS.md` lists under *Cutting a release*. They are deliberately
parallel, so that a step remembered by its number from one page is that step on the other.
That page carries the detail; if the two ever disagree, it is the one to trust.

1. Watch the run. `gh run watch <id> --exit-status`, or the Actions tab.
2. Check the draft. Both files attached, both named from `<version>`, and the notes naming
   the commit and each file's sha256.
3. Read the release page. Its copy is `build/RELEASE_PAGE.md` and the workflow has already
   assembled it into the draft (issue 265) — you read it here, you do not write it here.
4. **`make release-preflight`.** It refuses a draft built from a tree the repository has since
   moved on from. This step is not a formality: `mafigate-v1.0.0` was drafted from the tip of
   `main`, sat while twenty-two commits landed, and came within one click of publishing a known
   wrong-variant bug under the fixed version's number (issue 328). Nothing about the draft showed it.
5. Rehearse a second launch on a real Mac: install the draft's DMG by moving the app — drag
   to Applications or an explicit `mv`, never Gatekeeper translocation, which can silently
   not happen — launch, quit, move the `.app` again by hand, and launch again. First
   launches all pass; a moved app's second launch is what died in the field (issue 418).
6. Publish it. CI never does — it drafts and stops, so a human decides when a download exists.
7. Add the line here, and update the channel table in `streamlit_app/README.md` — the guard
   will hold you to it, and until step 6 it holds you to the opposite.

`build/BUILD_INSTRUCTIONS.md` has the rest of the release route, including how to rehearse
the workflow on a throwaway pre-release tag without publishing anything.
