# Contributing

Contributions are welcome, and this repository receives them in an unusual way. Please read
this first — it is short, and it changes what you should expect after you press the button.

## Where the code comes from

Development happens in a private repository at the European Institute of Oncology, on
patient data that cannot be published. What you see here is that tree, exported: each export
arrives as **one squashed commit** that replaces the working tree wholesale, with a message
naming the internal revision it was cut from.

Exports are not releases, and there are far more of them — around two hundred commits here
against seven tagged releases. So this repository's history is a list of exports, not a list of
changes. `git log` will not show you who wrote a line or why, and `git blame` will point at the
export rather than at the author. That is a deliberate safety property — the internal history
holds clinical data in past commits, and exporting only the current tree makes it structurally
impossible for any of it to travel — but it has consequences for you, below.

## Issues

**Please open them here.** The issue tracker on this repository is read.

An issue you file is triaged into the internal tracker, where the work actually happens.
Your issue **stays open here** until the fix reaches you in an export — so an issue that
sits open is not an issue that was ignored, and you will not see it closed the moment
someone starts on it. If we need more from you, we will ask in the thread.

Useful things to include: the version from the app's **About** entry, in the menu at the
top right of the window, which route you installed by (see the table in
[`streamlit_app/README.md`](streamlit_app/README.md)), and what the MAF looked like —
the arm, roughly how many rows, and which annotation columns it carries. Please do not
attach patient data.

## Pull requests

**Please open them here too**, and know what happens next: a pull request is read and
**reapplied by hand** into the internal tree, then it comes back to you in the next export.
It is not merged in the usual sense, and GitHub will most likely show it as closed rather
than merged even when every line of it shipped.

The cost of this falls on you rather than on us, so we would rather state it than have you
discover it: **your commit does not survive the export, so GitHub will not credit you.**
Your name will not appear in the contributors graph, and the commit that lands will not
carry your authorship.

We therefore credit contributors **manually**, and we would like to get it right:

- Say in the pull request how you want to be named, if it is not your GitHub handle.
- Tell us if you would rather not be named at all.
- We will tell you in the pull request thread which export carried your change.

One thing we cannot do, so that you are not looking for it: the credit is not in the commit.
The export commit's message is written by a script from a fixed one-line template —
`Sync from <private repo> @ <revision>` and nothing else — with nowhere for a name to go.

If that trade is not one you want to make, please open an issue describing the change
instead of writing it. A well-described issue is genuinely more useful to us than a patch
we cannot merge, and it costs you less.

## Before you write code

- **Say what you are planning** in an issue first, especially for anything larger than a
  fix. The internal tree may already be moving on the same thing, and the export means we
  cannot show you that it is.
- **Run the tests.** From `streamlit_app/`: `make test`. Several guards there assert facts
  about the tree rather than about behaviour, and they are meant to fail loudly.
- **Do not edit `streamlit_app/vendor/`.** It is a byte-for-byte copy of the pipeline's own
  filtering code, kept honest by a drift guard; read `streamlit_app/vendor/README.md` before
  changing anything near filtering.
- **Format before you push:** `make format` (black + isort) — but read the developer section of
  [`streamlit_app/README.md`](streamlit_app/README.md) first. That target prints a warning, does
  nothing and **exits 0** when the tools are not on your `PATH`, so a clean run is not evidence
  the tree is formatted; and this tree predates the current black, which reformats most of it if
  you let it loose.

## The app itself

For what MAFigate does, how to install it, and what each install route asks of you, see
[`streamlit_app/README.md`](streamlit_app/README.md). For the pipeline, see the
[README](README.md) at the root of this repository.
