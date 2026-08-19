<!--
The release page, written once and assembled by CI.

Issue #265. This is the body of every MAFigate release on GitHub, and it is the only
surface that reaches a Windows recipient before SmartScreen does — the installer is
blocked before any screen of ours can be drawn, so a person who is not told here is not
told at all. It is also where a reader of the paper who has no MAF of their own finds
that out, before spending a download rather than after.

Two lines below are placeholders, and `.github/workflows/mafigate-release.yml` fills
them in from the run that built the artifacts:

    {{ downloads }}         the two files, their platforms, their sha256 digests, and
                            the one commit both were built from
    {{ opening_the_app }}   build/OPENING_MAFIGATE.md, quoted whole

They are written in braces rather than as HTML comments for two reasons. Comments are
stripped from the published page — this block included, which is why it can be written
to whoever edits the file rather than to whoever reads the release — and a placeholder
that survived unfilled would then be invisible on the page instead of loud. The step
refuses to write a page still holding one either way.

The second placeholder is a placeholder rather than prose for a reason this repo has
already paid for once: the first-open wording lives in exactly one file and
`tests/test_unsigned_artifact_copy.py` fails a second copy of it. So this page quotes
that file; it does not summarise it, and it must not paraphrase it here.

That is also why it has no heading and no lead-in sentence of its own here. A draft of
this page had both, and the lead-in said what the note's own first paragraph says — a
second copy of the wording, arrived at in the one way the guard cannot see, by writing
an introduction to it. The quoted file carries its own title; let it.

No version number is written on this page. `tests/test_installer_version.py` sweeps
`build/` for version literals, and the filenames — which carry the version — come from
the downloads block, which is generated from what the build actually produced.
-->

## Before you download: MAFigate needs a MAF

MAFigate reads variant files; it does not produce them. To see anything at all you need an
**annotated `.maf`, `.tsv` or `.txt` written by the [VarianTalker
pipeline](https://github.com/zhanyinx/variantalker)** — or another MAF-format file carrying the
same annotation columns.

Without one, the app opens on a screen whose only offer is *choose a file*, and there is nothing
on this page or inside the download that will fill it. MAFigate ships with no example file.

If you have such a file, everything below applies to you. If you are here from the paper and do
not, this is the point at which to stop: read on only if you expect to have pipeline output.

## What it needs from your machine

- **macOS**, or **Windows on x64**. Neither installer targets Linux — on Linux, clone the
  repository and run `setup.sh`, which is a full route rather than a fallback.
- **Internet, once**, the first time you launch it. The installer brings its own Python, but the
  Python libraries are fetched on that first launch. Offline, the install still completes and the
  app then fails to start.
- **Free memory of roughly 140 MB plus ten times the size of the file you open.** That is the
  number to check before you open a large one: a 100 MB MAF needs about 1.1 GB free, and the
  largest we have measured — 118 MB — peaks near 1.4 GB. A typical file is a few megabytes and
  costs a few hundred.
- **Roughly half a gigabyte of disk** once it has finished. Most of that is the Python libraries
  fetched on the first launch, about 375 MB; the application itself accounts for about 110 MB on
  macOS, where it carries both Intel and Apple Silicon builds, and less on Windows.

## Downloads

{{ downloads }}

---

{{ opening_the_app }}

---

## When it does not work

Three things account for nearly all of it, and each looks different from a bug.

- **The first open is blocked.** Covered above; it happens once per machine, and neither the
  installer nor the app is damaged.
- **No internet on the first launch.** The install finishes, reports success, and then the app
  does not start — because the libraries it needs were never fetched. Connect and launch again.
- **The file is too big for the memory you have free.** Past the rule above, the app is stopped by
  the operating system, and it is stopped in a way that leaves no message and no error screen: the
  window simply goes away. Nothing is wrong with the MAF. Close other applications, or open the
  file on a machine with more memory.

Anything else: open an issue at
[zhanyinx/variantalker/issues](https://github.com/zhanyinx/variantalker/issues), and say which
file you have and what the screen showed.

## What happens to your file

**Your MAF: nowhere — nothing leaves your machine.** MAFigate reads it from your own disk and
makes no network request of its own. The ClinVar, gnomAD, dbSNP and UCSC links in the variant view
open in your browser; they are not lookups the app performs. Streamlit's own anonymous usage
reporting is turned off in what ships here.

That is also why **MAFigate never checks for new versions**. A check would be the one outbound
call the app makes, and it would make the paragraph above untrue. New releases are announced on
this page, and inside the app the About dialog names the version and the commit it was built
from — so you can always tell what you are running, and compare it against what is here.
