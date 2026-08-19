"""Where a user's parameters live between sessions: the cache, and the file they export.

One module behind five names the rest of the app calls — :func:`load_parameters_from_cache`,
:func:`save_parameters_to_cache`, :func:`get_cache_info`, :func:`clear_parameters_cache` and
:func:`show_discarded_cache_banner` — plus :func:`export_text`, the single serialiser behind
the page's four download buttons.

It was inside ``parameter_config.py``, and the give-away was ``MAFigate.py``: the app's boot
sequence restores the cache and renders the discard banner *before any page exists*, and to
do that it had to import three functions out of a page module. ``config/param_migration.py``
names the same tangle from the other side.

Why it is not in ``config/``, where its siblings :mod:`config.param_migration` and
:mod:`config.pipeline_params` are: ``config/`` and ``filters/`` are importable without
Streamlit, deliberately, so the parity harness and a bare ``pytest`` can read the app's
contract without booting a UI. This module reports its failures through ``st.error`` and
remembers the discard banner in ``st.session_state``. Rewriting that error path to raise
instead would be a behaviour change, so the module sits beside the page whose parameters it
keeps, and the streamlit-free line in ``config/`` stays intact.
"""

import json
from dataclasses import dataclass
from pathlib import Path

import streamlit as st
import yaml

from config.param_migration import (
    PARAM_SCHEMA_VERSION,
    SCHEMA_VERSION_KEY,
    param_document,
    unwrap_document,
)


# Cache file location
CACHE_DIR = Path.home() / ".mafigate"
PARAMS_CACHE_FILE = CACHE_DIR / "cached_parameters.json"

#: Where an unstamped cache is moved when it is discarded. Named beside the cache rather
#: than deleted: see :func:`discard_stale_cache`.
SUPERSEDED_CACHE_SUFFIX = ".superseded.json"


def save_parameters_to_cache(params):
    """Save parameters to the local cache file, in the stamped envelope.

    The stamp is the point. Without it the next version of this app faces exactly what
    this one faced — a file it cannot tell apart from one written five formats ago — and
    has to discard these caches too.
    """
    try:
        # Create cache directory if it doesn't exist
        CACHE_DIR.mkdir(exist_ok=True)

        with open(PARAMS_CACHE_FILE, "w") as f:
            json.dump(param_document(params), f, indent=2)

        return True
    except Exception as e:
        st.error(f"Failed to save parameters to cache: {str(e)}")
        return False


def _read_cache_document():
    """The cache file as a dict, or ``None`` if there is nothing readable there."""
    try:
        if not PARAMS_CACHE_FILE.exists():
            return None
        with open(PARAMS_CACHE_FILE, "r") as f:
            document = json.load(f)
        return document if isinstance(document, dict) else None
    except Exception:
        # Silently fail - cache is optional
        return None


def load_parameters_from_cache():
    """Parameters from the local cache — **only** if the cache carries the current stamp.

    An unstamped cache is not migrated and not read. It cannot be migrated safely: the
    migration of ``filter_variant_classification`` is a set complement, so applying it to
    a file that has already been converted returns a list of the same shape meaning the
    opposite, and no property of the value tells the two vintages apart. A parameter *file*
    escapes that because the user uploading it is the missing signal; a cache is restored
    before anything renders, with no one to ask.

    A cache from a *newer* format is likewise not read — this version cannot know what its
    values mean — but it is not set aside either; see :func:`discard_stale_cache`.

    :func:`discard_stale_cache` is what actually clears one away, so this function stays a
    pure read that a caller can stub.
    """
    document = _read_cache_document()
    if document is None:
        return None
    params, version = unwrap_document(document)
    if version != PARAM_SCHEMA_VERSION:
        return None
    return params


@dataclass(frozen=True)
class DiscardedCache:
    """A pre-stamp ``~/.mafigate`` cache that was moved aside, and what it held."""

    #: Where the file was moved to, so the user can upload it and migrate it deliberately.
    kept_at: str
    #: The arm it was configured for, as far as it can be read.
    arm: str
    #: When it was last saved.
    timestamp: str
    #: The parameter names it carried.
    keys: tuple = ()

    def summary(self):
        """One paragraph naming what was dropped. Shown as a banner, once per session."""
        when = f" saved {self.timestamp[:16].replace('T', ' ')}" if self.timestamp else ""
        held = f" It held {', '.join(self.keys)}." if self.keys else ""
        return (
            f"⚠️ Your saved parameter settings{when} were **discarded**, not migrated. "
            f"They were written for the **{self.arm}** arm by an older MAFigate, whose "
            "parameters do not mean the same thing here — and nothing in the file says "
            f"which version wrote it.{held} The app has opened at its own default "
            f"settings instead. The old file is kept at `{self.kept_at}`: upload it under "
            "*Presets → Upload Custom Preset* to migrate it deliberately."
        )


def discard_stale_cache():
    """Move a pre-stamp cache aside, and describe what it held. ``None`` if there is none.

    The acceptance criterion is that the cache is *discarded, not migrated*, and this is
    the whole of it. The reason is the one :mod:`config.param_migration` gives: every cache
    ever written carries the literal ``"app_version": "2.0.0"``, never bumped in any commit
    that touched it, so an old cache and a current one are the same state to any reader —
    and the app restored it before a single page rendered, which is how "the app opens at
    parity" was true only for a user who had never opened the app before.

    Moved rather than deleted. Discarding it is the app declining to guess; destroying the
    user's only copy of their settings is a different and unnecessary act, and the file is
    still perfectly migratable *by upload*, where the user supplies the one piece of
    information the stamp is missing.
    """
    document = _read_cache_document()
    if document is not None:
        params, version = unwrap_document(document)
        # Only an *older* format is set aside. A cache from a newer MAFigate is left
        # exactly where it is — this version declines to read it (see
        # `load_parameters_from_cache`) but has no business moving another version's file,
        # and the user is very likely to run that version again.
        if version is not None and version >= PARAM_SCHEMA_VERSION:
            return None
    elif not PARAMS_CACHE_FILE.exists():
        return None
    else:
        # Unreadable — a truncated write, or a hand-edit that did not parse. It cannot be
        # restored either, so it goes the same way rather than being read again next run.
        params = {}

    # Never clobber an earlier superseded copy. The whole point of moving rather than
    # deleting is that the user still has the file; a second discard overwriting the first
    # would take back that promise for the one user it matters to.
    kept_at = PARAMS_CACHE_FILE.with_suffix(SUPERSEDED_CACHE_SUFFIX)
    attempt = 1
    while kept_at.exists():
        attempt += 1
        kept_at = PARAMS_CACHE_FILE.with_suffix(f".superseded-{attempt}.json")

    try:
        PARAMS_CACHE_FILE.replace(kept_at)
    except Exception as e:
        st.error(f"Failed to set the old parameter cache aside: {str(e)}")
        return None

    return DiscardedCache(
        kept_at=str(kept_at),
        arm=params.get("sample_type", "somatic") if isinstance(params, dict) else "somatic",
        timestamp=(document or {}).get("timestamp", ""),
        keys=tuple(sorted(params)) if isinstance(params, dict) else (),
    )


def show_discarded_cache_banner():
    """Show the discard banner once per session, discarding the cache on first ask.

    Called from the app header and from the parameter page, so that whichever renders
    first tells the user — and only one of them does.
    """
    if "discarded_cache" not in st.session_state:
        st.session_state["discarded_cache"] = discard_stale_cache()

    discarded = st.session_state["discarded_cache"]
    if discarded is None or st.session_state.get("discarded_cache_shown"):
        return

    st.session_state["discarded_cache_shown"] = True
    st.warning(discarded.summary())


def get_cache_info():
    """Get information about the cached parameters."""
    document = _read_cache_document()
    if document is None:
        return None
    try:
        return {
            "timestamp": document.get("timestamp", "Unknown"),
            "app_version": document.get("app_version", "Unknown"),
            "schema_version": document.get(SCHEMA_VERSION_KEY, "Unknown"),
            "file_size": PARAMS_CACHE_FILE.stat().st_size,
        }
    except Exception:
        return None


def clear_parameters_cache():
    """Clear the parameters cache file."""
    try:
        if PARAMS_CACHE_FILE.exists():
            PARAMS_CACHE_FILE.unlink()
            return True
    except Exception as e:
        st.error(f"Failed to clear cache: {str(e)}")
        return False


def export_text(params, fmt):
    """Parameters as the text of a downloadable file, in the stamped envelope.

    One seam behind four download buttons — two formats, offered twice on the parameter
    page. They
    used to serialise ``st.session_state.filter_params`` directly at each site, which is
    how exports came to be the one artefact with no envelope: the cache recorded a
    timestamp and a version beside its parameters, and the file the user actually keeps
    recorded neither, so nothing could tell an export from a hand-written dict.
    """
    document = param_document(params)
    if fmt == "json":
        return json.dumps(document, indent=2)
    if fmt == "yaml":
        return yaml.dump(document, indent=2, sort_keys=False)
    raise ValueError(f"unknown export format {fmt!r}; expected 'json' or 'yaml'")
