"""The names of the files you have opened, kept between sessions.

A recents list, and deliberately nothing more. It answers one question — *what have I
looked at?* — for a user who has opened four samples this morning and cannot tell from the
app which. It is not a way back into any of them: see :func:`record_open_file`.

What is at rest, and why that is the whole design
-------------------------------------------------
An entry is a **name and a time**. Not a path, not the file's bytes, not its variant
counts, not the arm it was read under.

The names are the concession, and they were measured before being accepted: of 205 distinct
real MAF basenames on the dev's own machine, **153 carry a sample-code-shaped token** — in
three shapes: a cohort code with a normal/tumour suffix, a two-part accession with ``_TUM``,
a bare number with ``_blood``. The shapes are described rather than quoted, and deliberately:
a real one of each used to stand here as the example, and they are clinical identifiers, which
the publication scan gate keeps out of the tree. The next person to reach for a concrete
example wants ``CASE01ABCD-1N.maf``, the fabricated one ``tests/test_file_history.py`` uses.
The incremental exposure is small, because those codes are already on that disk in those
filenames. What is *not* already true
is that this file **outlives the MAF**: delete the sample and ``~/.mafigate`` still names it,
in a directory nobody thinks to look in.

So three things follow, and they are requirements rather than polish:

* the list is **capped** (:data:`HISTORY_LIMIT`), so it cannot grow into an archive;
* the app **says what it keeps**, in the sidebar beside the list itself; and
* there is a **way to clear it** that does not also throw away the user's filter settings —
  wanting the file names gone is not wanting to be reconfigured.

Why it is stamped, like the parameter cache
-------------------------------------------
``page_modules/param_store.py`` owns ``~/.mafigate`` today and writes a version-stamped
envelope, because issue #40's lesson was that an unstamped cache silently restored stale
data. This file follows that precedent, with one deliberate difference in what happens next:
a parameter cache in an unknown format is **moved aside and announced**, because it is the
user's only copy of settings they spent time on, while a history in an unknown format is
simply **dropped in silence**. A recents list is reconstructible by using the app, so there
is nothing to rescue and nothing worth a banner about.

Free of Streamlit, like everything else in ``config/``: the store is a file and a list, and
:mod:`components.sidebar` owns every sentence drawn from it.
"""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path


#: Where the history lives. The same directory the parameter cache uses — one place the app
#: keeps things between sessions, not two.
HISTORY_DIR = Path.home() / ".mafigate"
HISTORY_FILE = HISTORY_DIR / "opened_files.json"

#: The history format's own version, bumped by the change that alters the format.
#:
#: Separate from the parameter format's stamp on purpose: the two documents are written by
#: different code for different reasons, and tying them together would mean a change to one
#: format invalidating the other's files.
HISTORY_SCHEMA_VERSION = 1

#: Where the stamp and the list live inside the envelope.
SCHEMA_VERSION_KEY = "schema_version"
FILES_KEY = "files"

#: The keys of one entry.
NAME_KEY = "name"
OPENED_KEY = "opened"

#: How many names the list keeps.
#:
#: Ten, and the number is doing work rather than being round. It is the cap that makes this
#: a *recents list* rather than a log: a clinician reading the sidebar is looking for a file
#: from this morning or yesterday, and an eleventh name is a name still on disk after it has
#: stopped being useful to anybody. It is also small enough to read without scrolling in a
#: sidebar column that already holds four things.
HISTORY_LIMIT = 10


def _read_document() -> dict | None:
    """The history file as a dict, or ``None`` if there is nothing readable there."""
    try:
        if not HISTORY_FILE.exists():
            return None
        with open(HISTORY_FILE, "r", encoding="utf-8") as handle:
            document = json.load(handle)
        return document if isinstance(document, dict) else None
    except Exception:
        # A recents list is optional. A truncated write or a hand-edit that does not parse
        # costs the user the list, not the session, and there is nothing here worth an
        # error banner beside their variants.
        return None


def _valid_entry(entry) -> bool:
    """Whether one item out of a stored list is an entry this version can draw.

    Checked per entry rather than trusted from the stamp, because the stamp only says which
    *writer* produced the file and a hand-edited or partially-written list can carry the
    right stamp over the wrong rows. An entry missing its name would render as a blank line
    in the sidebar, which reads as a file with no name rather than as a damaged store.
    """
    return (
        isinstance(entry, dict)
        and isinstance(entry.get(NAME_KEY), str)
        and bool(entry[NAME_KEY])
    )


def read_history() -> list[dict]:
    """The files opened, newest first — **only** if the file carries the current stamp.

    An unstamped or newer-stamped document reads as an empty history. Unlike the parameter
    cache it is not moved aside and not announced: the user loses a list they can rebuild by
    opening files, so there is nothing to rescue, and a banner about the internals of a
    recents list would be the app talking about itself.
    """
    document = _read_document()
    if document is None:
        return []
    if document.get(SCHEMA_VERSION_KEY) != HISTORY_SCHEMA_VERSION:
        return []
    files = document.get(FILES_KEY)
    if not isinstance(files, list):
        return []
    return [entry for entry in files if _valid_entry(entry)][:HISTORY_LIMIT]


def record_open_file(name: str, now: datetime | None = None) -> list[dict]:
    """Put ``name`` at the top of the history, and return the list as it now stands.

    Callers pass the name of a file the app has **already opened**. That is the one rule
    about *when* this is called and it is load-bearing: a MAF the app could not read, or
    read and then refused for a missing column, leaves no entry — so the list means *files
    you have looked at* rather than *files you have pointed at*. The account of a refusal is
    the error on the page, which has room for the reason.

    Re-opening a file **moves** its entry rather than adding a second one. A recents list
    with the same sample three times in it is a log of interactions, and the question this
    answers is which files, not how often.

    Only the basename is stored, even where the caller has a full path — the OS "Open With"
    route does. A path is a way back to the file, and this list is deliberately not one
    (see the module docstring); it is also a second piece of information about where a
    patient's data sits on this disk.

    Silent on failure, for the reason :func:`_read_document` gives: a home directory that
    cannot be written to costs the user their recents list and must not cost them the file
    they just opened.
    """
    entry_name = Path(str(name)).name
    if not entry_name:
        return read_history()

    stamped = (now or datetime.now()).isoformat(timespec="seconds")
    kept = [entry for entry in read_history() if entry[NAME_KEY] != entry_name]
    files = [{NAME_KEY: entry_name, OPENED_KEY: stamped}, *kept][:HISTORY_LIMIT]

    try:
        HISTORY_DIR.mkdir(exist_ok=True)
        with open(HISTORY_FILE, "w", encoding="utf-8") as handle:
            json.dump(
                {SCHEMA_VERSION_KEY: HISTORY_SCHEMA_VERSION, FILES_KEY: files},
                handle,
                indent=2,
            )
    except Exception:
        return files

    return files


def clear_history() -> bool:
    """Delete the history file. Returns whether there is now no history to read.

    Deleted rather than moved aside, which is the other half of the difference from
    :func:`page_modules.param_store.discard_stale_cache`: that function keeps a superseded
    parameter cache because it is the user's only copy of work they did, while a user
    clearing this list is asking for these names to stop being on their disk. Keeping a copy
    beside it would be the app declining to do the one thing the control promises.
    """
    try:
        HISTORY_FILE.unlink(missing_ok=True)
        return True
    except Exception:
        return False


def format_opened(stamp: str, now: datetime | None = None) -> str:
    """When a file was opened, as the sidebar says it.

    *today 14:02* for today, *16 Aug 09:20* otherwise. Relative for today because that is
    the distinction the list is actually for — telling this morning's sample from the one
    before it — and absolute beyond that because "3 days ago" makes a reader do arithmetic
    to compare it with anything else on their screen.

    Returns the empty string for a stamp this version cannot read, and the caller draws the
    name alone. A list of names with no times is still the list; a list with ``Unknown``
    beside three of them is the app reporting on its own storage.
    """
    if not isinstance(stamp, str) or not stamp:
        return ""
    try:
        opened = datetime.fromisoformat(stamp)
    except ValueError:
        return ""

    reference = now or datetime.now()
    if opened.date() == reference.date():
        return f"today {opened:%H:%M}"
    # `%-d` is a platform extension; built by hand so the string is the same everywhere.
    return f"{opened.day} {opened:%b %H:%M}"
