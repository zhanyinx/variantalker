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
from config.pipeline_params import pipeline_params, stated_arm
from config.vocabularies import GUIDELINE_SOURCES


# Cache file location
CACHE_DIR = Path.home() / ".mafigate"
PARAMS_CACHE_FILE = CACHE_DIR / "cached_parameters.json"

#: Where an unstamped cache is moved when it is discarded. Named beside the cache rather
#: than deleted: see :func:`discard_stale_cache`.
SUPERSEDED_CACHE_SUFFIX = ".superseded.json"

#: Why a cache was set aside. The two reasons need two sentences, not one: almost nothing
#: the unstamped wording says is true of an incomplete one — see :meth:`DiscardedCache.summary`.
DISCARDED_UNSTAMPED = "unstamped"
DISCARDED_INCOMPLETE = "incomplete"


def missing_contract_keys(params):
    """Which of its own arm's contract keys a parameter dict does not carry.

    Read off :func:`config.pipeline_params.pipeline_params` rather than a key list restated
    here, so a parameter added to the contract cannot silently stop being required.

    The arm is the dict's **own** ``sample_type``, and a ``sample_type`` that is not an arm
    is reported as the missing key: there is no contract such a dict could be complete
    against, and the app could not have written it either, which is what the caller is
    really asking (issue #279). That reading is :func:`config.pipeline_params.stated_arm`,
    shared with the data page's refusal rather than spelled twice (issue #289).

    Emptiness is never consulted. An empty keep-list is a legal choice — the pipeline's own
    ``--filter_cancervar ""`` — so this asks only which keys are *absent*, which is the one
    thing the issue #276 defect produced and a deliberate choice cannot forge.
    """
    arm = stated_arm(params)
    if arm is None:
        return ("sample_type",)
    return tuple(sorted(set(pipeline_params(arm)) - set(params)))


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

    What this route deliberately does **not** do is complete a partial cache
    ---------------------------------------------------------------------
    It could: an incomplete set is exactly the asymmetry issue #276 found — the upload
    route's ``migrate_params`` is documented to always produce *a complete parameter set for
    one arm* while this one adopts whatever the file carries — and
    :func:`config.param_migration.complete_params` would fill it in a line.

    It is left alone because issue #279 decided what such a cache *means*, and the answer was
    not "quietly repair it": a stamped cache is restored only if its parameters are a complete
    set for their own arm, and one that is not is superseded and named to the user, because
    incompleteness is the #276 defect's signature and a deliberate choice cannot forge it.
    Completing here would make that condition unreachable — anything measured after this
    function would already be complete — so what happens to an incomplete cache is *refusal*,
    on the raw document, and :func:`discard_stale_cache` moves the file and says so
    (issue #286).

    Both halves are needed, because the boot sequence puts this read **first**:
    ``initialize_session_state`` restores the cache and ``render_header`` draws the banner,
    in that order, so a read that adopted an incomplete cache and left the moving to the
    banner would open the session on the very parameters it was about to set aside.

    What stops a partial dict reaching a widget by any *other* route is the completion at the
    top of ``show_parameter_config_page``, where every route's parameters pass through.
    """
    document = _read_cache_document()
    if document is None:
        return None
    params, version = unwrap_document(document)
    if version != PARAM_SCHEMA_VERSION:
        return None
    if missing_contract_keys(params):
        return None
    return params


@dataclass(frozen=True)
class DiscardedCache:
    """A ``~/.mafigate`` cache that was moved aside, and what it held."""

    #: Where the file was moved to, so the user can upload it and migrate it deliberately.
    kept_at: str
    #: The arm it was configured for, as far as it can be read.
    arm: str
    #: When it was last saved.
    timestamp: str
    #: The parameter names it carried.
    keys: tuple = ()
    #: Why it was set aside: :data:`DISCARDED_UNSTAMPED` or :data:`DISCARDED_INCOMPLETE`.
    reason: str = DISCARDED_UNSTAMPED
    #: For an incomplete cache, the contract keys its own arm names and it did not carry.
    missing: tuple = ()
    #: For an incomplete cache, whether *every* guideline source in it was empty — the
    #: state that quietly cut the report down to the pathogenic-rescue floor.
    no_guideline_source_kept_anything: bool = False

    def summary(self):
        """One paragraph naming what was dropped. Shown as a banner, once per session.

        Two reasons, two sentences. Reusing the unstamped wording for an incomplete cache
        would tell the user three things that are false of it — that an *older* MAFigate
        wrote it (this one did), that nothing in the file says which version wrote it (it is
        stamped), and that uploading it migrates it (for a stamped file ``migrate_params``
        is a verbatim no-op, so the remedy restores the damage). Issue #279 measured that
        last one; issue #286 is why there is a branch here rather than a second banner.
        """
        when = f" saved {self.timestamp[:16].replace('T', ' ')}" if self.timestamp else ""
        if self.reason == DISCARDED_INCOMPLETE:
            return self._incomplete_summary(when)
        held = f" It held {', '.join(self.keys)}." if self.keys else ""
        return (
            f"⚠️ Your saved parameter settings{when} were **discarded**, not migrated. "
            f"They were written for the **{self.arm}** arm by an older MAFigate, whose "
            "parameters do not mean the same thing here — and nothing in the file says "
            f"which version wrote it.{held} The app has opened at its own default "
            f"settings instead. The old file is kept at `{self.kept_at}`: upload it under "
            "*Presets → Upload Custom Preset* to migrate it deliberately."
        )

    def _incomplete_summary(self, when):
        """The wording for a cache this app could only have written by inventing values.

        Says what was actually measured and nothing more. The missing parameters are named
        because they are the evidence; the guideline clause is drawn only when every source
        really was empty, so the sentence cannot be false of a cache that lost one key and
        kept its keep-lists.
        """
        if "sample_type" in self.missing:
            # No arm, no contract, and nothing to name as missing but the arm itself —
            # which is what `missing_contract_keys` reports and, for a dict that *does*
            # state an arm, can never report. Not read off `self.arm`: that field carries
            # the same `.get("sample_type", "somatic")` fallback every reader in the app
            # applies, so it says "somatic" for a file that said nothing at all.
            #
            # The app reads such a file as somatic everywhere else; saying "the somatic arm
            # needs" would put words in it that it never said.
            return (
                f"⚠️ Your saved parameter settings{when} were **not the ones you chose**, "
                "so MAFigate has set them aside. The file does not say which analysis arm "
                "it was written for, so there is no way to tell what its settings mean. "
                "The app has opened at its own settings instead. Your old file is kept at "
                f"`{self.kept_at}`, so nothing of yours is destroyed."
            )

        quoted = [f"`{key}`" for key in self.missing]
        named = " or ".join([", ".join(quoted[:-1]), quoted[-1]] if len(quoted) > 1 else quoted)
        collapsed = (
            " Every guideline source in it was empty, so any report cut while it was in "
            "effect was reduced to the **pathogenic-rescue** floor — only the variants "
            "pathogenic retention keeps on their own."
            if self.no_guideline_source_kept_anything
            else ""
        )
        return (
            f"⚠️ Your saved parameter settings{when} were **not the ones you chose**, so "
            f"MAFigate has set them aside. The file did not carry {named}, which the "
            f"**{self.arm}** arm needs — so MAFigate filled those in for itself and saved "
            "them back as though you had chosen them. That was a defect in the app, not "
            f"anything you did.{collapsed} The app has opened at its own settings "
            f"instead. Your old file is kept at `{self.kept_at}`, so nothing of yours is "
            "destroyed — but there is nothing to upload back: these are the settings that "
            "were wrong."
        )


def discard_stale_cache():
    """Move a cache this app cannot honestly restore aside, and say what it held.

    ``None`` if there is none. Two kinds go:

    * **unstamped, or stamped lower than this format** — discarded, not migrated. The reason
      is the one :mod:`config.param_migration` gives: a pre-stamp cache carries an
      ``app_version`` and nothing else that dates it, and that number has only ever moved for
      a packaging decision (issue #260) rather than a format one, so an old cache and a
      current one are the same state to any reader — and the app restored it before a single
      page rendered, which is how "the app opens at parity" was true only for a user who had
      never opened the app before.
    * **stamped with this format, but not a complete set for its own arm** — retired. Such a
      dict is one this app could only have written by *inventing* the values it could not
      find (issue #276), so restoring it opens the app on settings nobody chose; issue #279
      measured that no contract, preset or deliberate choice produces one, on either arm.

    A cache from a **newer** format is left where it is whether it is complete or not: this
    version declines to read it, but a format it has never seen is not its to judge, and the
    user is very likely to run that version again.

    Moved rather than deleted, in both cases. Discarding it is the app declining to guess;
    destroying the user's only copy of their settings is a different and unnecessary act. For
    the unstamped kind the file is also still perfectly migratable *by upload*, where the user
    supplies the one piece of information the stamp is missing — which is **not** true of the
    incomplete kind, where ``migrate_params`` is a verbatim no-op and re-uploading would
    restore exactly the settings that were wrong. The banner says so per reason.
    """
    document = _read_cache_document()
    reason = DISCARDED_UNSTAMPED
    missing = ()
    if document is not None:
        params, version = unwrap_document(document)
        if version is not None and version > PARAM_SCHEMA_VERSION:
            return None
        if version == PARAM_SCHEMA_VERSION:
            # This version's own stamp. The only thing that retires it is incompleteness,
            # measured on the raw document — before the completion at the top of
            # `show_parameter_config_page`, which would otherwise repair the file out from
            # under the test that is supposed to catch it (issue #286).
            missing = missing_contract_keys(params)
            if not missing:
                return None
            reason = DISCARDED_INCOMPLETE
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
        reason=reason,
        missing=missing,
        no_guideline_source_kept_anything=_no_guideline_source_kept_anything(params),
    )


def _no_guideline_source_kept_anything(params):
    """Whether every guideline source in a dict was empty.

    Only the *banner* asks this, and only to say what a report cut under those settings
    actually was. Nothing in the retirement condition consults emptiness — issue #279's
    whole point — so this is a description of the file, never a reason to move it.

    The sources are :data:`config.vocabularies.GUIDELINE_SOURCES`, the same relation
    ``_warn_if_every_guideline_source_is_empty`` reads, so the banner and the parameters
    page cannot come to disagree about what "every source" means. It matters that they
    do not: that warning draws **only** on the parameters page, so this banner is the
    only thing a user who loaded a MAF and cut a report without opening that page ever
    saw about it.
    """
    arm = stated_arm(params)
    if arm is None:
        return False
    sources = GUIDELINE_SOURCES[arm]
    return not any(params.get(key) for key in sources)


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
