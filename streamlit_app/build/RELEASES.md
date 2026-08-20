# Published MAFigate releases

The tree's record of which releases have actually been **published** — not which ones could
be built, and not which ones are drafted. One line per release, newest first.

This file exists because a promise in the README is only true once a download exists, and
something checkable had to answer *does it?*. Issue #264 re-keyed that switch onto this
file: the guard in `tests/test_delivery_channels_copy.py` reads the entries below, and while
there are none the README may not link to a release artifact and must say the installers are
not available yet. Adding an entry here flips it, in both directions.

Before this, the switch keyed on whether a release workflow existed in `.github/workflows`.
That is a **capability**, not an event, and it was wrong by exactly the window between the
workflow landing (#264) and the first release being cut (#265): the workflow's own arrival
would have made the guard skip — not fail, *skip* — leaving both wrong states unguarded in
the one window where a README could link to a 404. A record of the event cannot be wrong
that way. The cost is that the record is written by hand, which is the same cost every
other prose guard in this suite carries.

## Releases

**No release has been published yet.** The workflow that builds the installers exists
(`.github/workflows/mafigate-release.yml`) and drafts a release on a `mafigate-v*` tag, but a
draft is not a download. Issue #265 cuts the first one.

That sentence is load-bearing: the guard demands the ledger say one of its two states out
loud, so a file that lost its entries — or never grew any — cannot read as *not published
yet* by accident. Delete the sentence when the first entry lands.

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

1. `cd streamlit_app && make release-tag` — that is the tag, derived from `APP_VERSION`.
2. Push it. The workflow builds the DMG on macOS and the installer on Windows, from the one
   commit the tag names, and attaches both to a draft release.
3. Read the release page. Its copy is `build/RELEASE_PAGE.md` and the workflow has already
   assembled it into the draft (issue #265) — you read it here, you do not write it here.
4. **`make release-preflight`.** It refuses a draft built from a tree the repository has since
   moved on from. This step is not a formality: `mafigate-v1.0.0` was drafted from the tip of
   `main`, sat while twenty-two commits landed, and came within one click of publishing a known
   wrong-variant bug under the fixed version's number (#328). Nothing about the draft showed it.
5. Publish it.
6. Add the line here, and update the channel table in `streamlit_app/README.md` — the guard
   will hold you to it, and until step 5 it holds you to the opposite.

`build/BUILD_INSTRUCTIONS.md` has the rest of the release route, including how to rehearse
the workflow on a throwaway pre-release tag without publishing anything.
