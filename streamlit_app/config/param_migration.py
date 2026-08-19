"""Bringing a parameter file written by an older MAFigate onto the contract.

``config/pipeline_params.py`` says what a parameter dict means *now*. This module says
what an older one meant, and how to say the same thing in the current vocabulary. It is
free of Streamlit, so the page, the tests and any future command line can all reach it.

The gate is a version stamp, and it has to be
---------------------------------------------
The tempting design is to sniff the dict — "this looks like a legacy include list, so
invert it". It cannot work. The contract's own default for
``filter_variant_classification`` is ``["Silent", "IGR", "RNA"]``, which is legal as the
post-parity *exclude* list it now is **and** as a legacy *include* list naming the three
classifications a user wanted reported. Membership, length and vocabulary are identical in
the two readings; nothing in the value discriminates them.

That would be merely imprecise if the migration were idempotent. It is not: the
classification migration is a set complement, an involution, so guessing wrong on a
post-parity file returns a list of the same shape meaning the opposite. Measured on
``somatic_reference.maf``, migrating the parity value a second time takes the criteria
path from 19 rows to 2.

So a file is migrated if and only if it does not carry :data:`PARAM_SCHEMA_VERSION`, and
every file this app writes carries it from now on.

Why the ``~/.mafigate`` cache is discarded rather than migrated
--------------------------------------------------------------
The stamp is powerless over files that already exist, and the cache is nothing but files
that already exist. Every pre-stamp cache carries an ``app_version`` and nothing else that
dates it, and the literal it carries had never been bumped in any commit that touched it —
so "an old cache" and "a cache written by an app that stamps its output" are the same state
to any reader.

Issue #260 has since given that number a reason to move: it is ``APP_VERSION``, and it is
now the one source both desktop installers and the release tag cut from. That does not make
it able to date a cache — it moved for a packaging decision, which is precisely the kind of
event a format version must not follow, and the argument below for keeping the two apart.

For an uploaded file that does not matter: the *user* supplies the file, and the act of
uploading is the missing information ("this one is old"). For the cache there is no such
act — it is restored before anything renders — so the only safe reading of an unstamped
cache is to not read it. ``page_modules/parameter_config.discard_stale_cache`` moves it
aside and names what it held, which leaves the user able to upload it if they want it
back.

Why every per-key rule below is a widening if you get it wrong
--------------------------------------------------------------
Each mistake is silent and they compound: an unrenamed gene key drops the gene
restriction entirely and is worth about sixteen extra rows per arm on the reference
fixtures; an uninverted pathogenic flag two more; a germline dict that keeps its ESCAT key
carries a value nothing reads while looking like a setting. None of them raises, and the
report is simply bigger than the one the user's file asked for.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass, field
from datetime import datetime

from config.constants import APP_VERSION
from config.pipeline_params import ARMS, PIPELINE_PARAMS, pipeline_params
from filters.gene_lists import parse_gene_list

#: The parameter-file format's own version, bumped by the change that alters the format.
#:
#: Deliberately **not** ``APP_VERSION``. They were the same field once — the cache's
#: ``app_version`` — and that is exactly why the cache cannot be migrated: a release
#: number moves for reasons that have nothing to do with the format, and for the whole life
#: of that cache it had not moved at all. Issue #260 then moved it, for a packaging
#: decision: the first public release is numbered from what the public README already
#: promised. No format changed that day, which is the argument in one event. Tying them
#: together would mean a patch release renumbering every saved file and a format change
#: waiting on a release.
PARAM_SCHEMA_VERSION = 1

#: Where the stamp lives in an exported document, and in the cache.
SCHEMA_VERSION_KEY = "schema_version"

#: Where the parameters live inside the envelope. The name the cache already used, so a
#: cache file and an exported file are the same document with the same reader.
PARAMETERS_KEY = "parameters"

#: The catch-all a pre-#36 file can carry in any multiselect parameter. Stripped from
#: keep-lists, where it and ``[]`` already decide every row identically; read as "every
#: classification" in the classification list, where it meant the opposite of an empty
#: list. See :func:`_migrate_classification`.
LEGACY_SENTINEL = "All"

#: The variant classifications an older MAFigate could offer, **frozen**.
#:
#: A copy of ``constants.VARIANT_CLASSIFICATIONS`` as it stood when the classification
#: parameter changed polarity, and it must stay a copy. The migration takes a set
#: complement, which is only well defined against the vocabulary the *writer* could choose
#: from: taken against the live list, a classification added next year would start being
#: excluded by every file written before it existed. Retroactive, silent, and invisible in
#: the file itself.
LEGACY_VARIANT_CLASSIFICATIONS = (
    "Missense_Mutation",
    "Nonsense_Mutation",
    "Silent",
    "Frame_Shift_Del",
    "Frame_Shift_Ins",
    "In_Frame_Del",
    "In_Frame_Ins",
    "Nonstop_Mutation",
    "Translation_Start_Site",
    "Splice_Site",
    "Splice_Region",
    "5'UTR",
    "3'UTR",
    "5'Flank",
    "3'Flank",
    "Intron",
    "IGR",
    "RNA",
)

#: The six guideline keep-lists, in the order a note names them.
#:
#: A second copy of a relation :data:`config.vocabularies.GUIDELINE_SOURCES` already states.
#:
#: The original reason was that the relation lived in a Streamlit page and this module must
#: stay importable without Streamlit. That is fixed — it is in ``config/`` now — but the copy
#: survives for a narrower reason: :mod:`config.vocabularies` imports
#: :data:`LEGACY_SENTINEL` from *here*, so importing the relation back would close an import
#: cycle. Dissolving the copy means moving the sentinel into the vocabulary module too,
#: which is a separate change with its own consequences for the test below.
#:
#: Until then the remedy is unchanged — restate, then **pin it with a test**.
#: ``test_param_migration`` checks this tuple against the union of the per-arm sources, so a
#: seventh guideline source cannot be added to the contract and silently dropped from every
#: uploaded file as an "unknown key".
_GUIDELINE_KEYS = (
    "filter_cancervar",
    "filter_civic",
    "filter_clinvar",
    "filter_intervar",
    "filter_renovo",
    "filter_escat",
)

#: Per-source spelling repairs, kept because the contract requires exactly this and no
#: more: **InterVar and ReNovo values carry spaces, CancerVar and ClinVar values carry
#: underscores**, because that is how the annotations are written. There is no
#: normalisation to apply here — a single rule either way would break two of the four
#: sources — so what a migration can do is repair the spellings an older app wrote.
_SPACING_REPAIRS = {
    "filter_cancervar": {
        "Tier_I_Strong": "Tier_I_strong",
        "Tier_II_Potential": "Tier_II_potential",
        "Tier_III_Unknown": "Tier_III_Uncertain",
        "Tier_IV_Benign": "Tier_IV_benign",
    },
    "filter_intervar": {
        "Likely_pathogenic": "Likely pathogenic",
        "Uncertain_significance": "Uncertain significance",
        "Likely_benign": "Likely benign",
    },
    # The two composite repairs this map dropped at issue #88 -- "Benign/Likely benign"
    # and "Pathogenic/Likely pathogenic" -- are gone rather than re-pointed. They repaired
    # a spelling into a keep-term that no longer exists, and re-pointing them at an atom
    # would be this module *choosing* which half of a two-call selection the user meant.
    # A saved file carrying either now passes through unrepaired and is dropped at the
    # widget by `filter_terms`, which is what already happened to the repaired spelling:
    # it never matched a row (see CLINVAR_OPTIONS in config/vocabularies.py), so the
    # migration's output is unchanged in the only sense that reaches a report.
    "filter_clinvar": {
        "Likely pathogenic": "Likely_pathogenic",
        "Uncertain significance": "Uncertain_significance",
        "Likely benign": "Likely_benign",
        "Conflicting classifications of pathogenicity": (
            "Conflicting_classifications_of_pathogenicity"
        ),
        "Conflicting interpretations of pathogenicity": (
            "Conflicting_interpretations_of_pathogenicity"
        ),
        "Drug response": "drug_response",
        "Not provided": "not_provided",
    },
    "filter_renovo": {
        "Likely_pathogenic": "LP Pathogenic",
        "Pathogenic": "HP Pathogenic",
        "Uncertain_significance": "LP Benign",
        "Likely_benign": "IP Benign",
        "Benign": "HP Benign",
    },
}

#: The gene-panel keys an older app wrote. Dropped rather than mapped: the panel *name* is
#: UI state that resolves to a gene list, never a filter parameter (issue #34), and the
#: page recovers the panel from the symbols themselves — so a saved panel choice survives
#: in the only form that ever filtered anything.
_PANEL_KEYS = ("somatic_gene_set", "germline_gene_set", "gene_set")

#: Scalars whose name and meaning never changed.
_UNCHANGED_SCALARS = ("min_depth", "max_freq_population", "skip_civic")

#: Straight renames, where only the app's name for the parameter moved.
_RENAMED_SCALARS = {"min_freq": "max_freq_population"}


@dataclass(frozen=True)
class Migration:
    """What a parameter document became, and what was left behind on the way."""

    #: The parameters, on the contract, for the document's own arm.
    params: dict
    #: Whether the legacy transform ran. ``False`` for an already-stamped document, which
    #: is the no-op the re-upload case needs.
    migrated: bool
    #: Keys from the document that have no home in the migrated dict, named so nothing
    #: disappears silently.
    dropped: tuple = ()
    #: Lines for a human, in the order the migration made the decisions.
    notes: tuple = ()


def param_document(params, *, app_version: str = APP_VERSION, timestamp=None) -> dict:
    """Wrap parameters in the stamped envelope every export and cache write carries.

    The envelope is the cache's own, with the stamp added: the cache already recorded a
    timestamp and a version beside its parameters while exports were a bare dict, which is
    why an exported file could never be told apart from a hand-written one.
    """
    return {
        SCHEMA_VERSION_KEY: PARAM_SCHEMA_VERSION,
        "app_version": app_version,
        "timestamp": timestamp or datetime.now().isoformat(),
        PARAMETERS_KEY: dict(params),
    }


def unwrap_document(document) -> tuple:
    """``(parameters, schema version or None)`` for an envelope *or* a bare dict.

    Every file written before this ticket is a bare parameter dict, so the unwrapping has
    to accept both. ``None`` means unstamped, which is the state that gets migrated —
    including the pre-stamp cache envelope, which carries ``parameters`` and an
    ``app_version`` but no schema version.
    """
    if not isinstance(document, dict):
        return {}, None
    if isinstance(document.get(PARAMETERS_KEY), dict):
        version = document.get(SCHEMA_VERSION_KEY)
        return document[PARAMETERS_KEY], version if isinstance(version, int) else None
    return document, None


def complete_params(params) -> dict:
    """One parameter dict, with every key its arm's contract names but it does not carry.

    Why this exists, and why it is **at the boundary** rather than in each control
    ---------------------------------------------------------------------------------
    A dict reaching the parameters page with a key missing used to have that key *invented*
    by whichever widget drew it: a list control read ``.get(key)`` with no fallback, so
    ``filter_terms(None, options)`` was ``[]`` — *keep nothing* — and a number control read
    ``.get(key, <literal>)`` and produced something that looked chosen. The page then wrote
    the invention back as the user's own selection and the auto-save persisted it, and since
    the page's written key set is a fixed point the file re-wrote itself identically for ever
    (issue #276).

    The repair completes the dict **once, where it enters**, instead of giving each ``.get``
    a better fallback, for a reason the code itself already demonstrated: the per-control
    shape is the one that drifted. Three of those literals had come to contradict the
    contract — ``max_freq_population`` invented ``0.01`` where both arms say ``1.0``, the
    somatic VAF widget invented ``0.05`` where the contract says ``0.01``, and the germline
    VAF widget invented its own ``0.2`` under a key the contract does not even spell. A
    single completion also covers keys *no* control draws (``skip_civic``) and any control
    added later, which no amount of per-widget care can promise.

    What it does **not** do: overrule a value the dict states
    --------------------------------------------------------
    **An absent key and a present-but-empty key are different states, and only the first is
    a defect.** An empty keep-list is a legal, expressible choice — the pipeline's own
    ``--filter_cancervar ""`` — which is why issue #36 deleted the empty→``["All"]``
    backfill and why the page *warns* about an all-empty guideline block instead of
    preventing it. So this fills what is **missing** and never touches what is **empty on
    purpose**; a cache already holding empty clinical filters is untouched here, and what to
    do about one is issue #279's question, not this function's.

    A key can also be stated under a second spelling, and that counts as stated. Two pairs
    exist and both are read exactly as :func:`migrate_params` reads them — deliberately, so
    that a dict arriving by upload and one arriving any other way cannot end up meaning
    different things by the same key:

    * ``keep_pathogenic`` is the app's polarity for the contract's ``skip_pathogenic`` — all
      four presets still carry only the former, so taking the contract's value would let
      this function silently overrule a preset's own retention setting;
    * on germline the widget and the filter both prefer ``vaf_threshold_germline`` while the
      contract spells the parameter ``vaf_threshold``, so the germline arm carries both keys
      set to one number, as ``_migrate_vaf`` already arranges.

    Why it is not :func:`migrate_params`
    ------------------------------------
    Because that function's legacy transform is **not idempotent** and this one is run on
    every render. ``_migrate_classification`` takes a set complement, so handing a
    current dict to the legacy path returns a ``filter_variant_classification`` of the same
    shape meaning the opposite — the module header's own argument, and the reason completion
    is a separate, additive function that never rewrites a value.

    Args:
        params: a parameter dict, complete or not. Filled **in place** — the app's pages
            edit ``st.session_state.filter_params`` in place and callers hold references to
            it — and returned for the convenience of a caller that wants an expression.

    Returns:
        The same dict, now carrying every key its own arm's contract names. Values taken
        from the contract are fresh deep copies (:func:`pipeline_params`), so a widget
        appending to a keep-list cannot redefine the contract for the rest of the process.
    """
    if not isinstance(params, dict):
        params = {}

    stated = params.get("sample_type")
    # The arm is read the way every other reader reads it: a dict that states no arm is
    # somatic to this app, and so is one that states something that is not an arm.
    arm = stated if stated in ARMS else "somatic"
    contract = pipeline_params(arm)

    # The pathogenic flag before the plain fill, so a dict carrying only the app's polarity
    # keeps its own answer rather than taking the contract's.
    if "skip_pathogenic" not in params and "keep_pathogenic" in params:
        params["skip_pathogenic"] = not bool(params["keep_pathogenic"])

    # The germline VAF threshold under both spellings, from whichever one the dict has.
    if arm == "germline":
        germline_vaf = params.get(
            "vaf_threshold_germline",
            params.get("vaf_threshold", contract["vaf_threshold"]),
        )
        params.setdefault("vaf_threshold_germline", germline_vaf)
        params.setdefault("vaf_threshold", germline_vaf)

    for key, value in contract.items():
        # `sample_type` is arm *identity*, not a filter setting — it selects which filter
        # runs rather than tuning one, which is why `pipeline_params` gives it no
        # `ParamOrigin`. This function reads the arm to know which contract to complete
        # against and has no business deciding it: a dict that states no arm goes on stating
        # none, and is read as somatic by every reader in the app exactly as before. Writing
        # it here would make completion a path that can move a user between arms without
        # saying so, which is the one thing `adopt_parameters` exists to prevent.
        if key == "sample_type":
            continue
        params.setdefault(key, value)

    return params


def migrate_params(document) -> Migration:
    """A parameter document, as parameters the current filter reads.

    Args:
        document: an uploaded file's contents — a stamped envelope, a bare parameter dict
            written by any version of this app, or the pre-stamp cache envelope.

    Returns:
        A :class:`Migration`. Its ``params`` are always a complete parameter set for one
        arm: the contract for that arm, with the document's own values layered on top.
    """
    params, version = unwrap_document(document)

    if version is not None and version >= PARAM_SCHEMA_VERSION:
        notes = ()
        if version > PARAM_SCHEMA_VERSION:
            notes = (
                f"This file was written by a newer MAFigate (parameter format {version}; "
                f"this one reads {PARAM_SCHEMA_VERSION}). It was loaded unchanged, "
                "because there is no way to know what a format we have never seen means.",
            )
        return Migration(params=copy.deepcopy(params), migrated=False, notes=notes)

    return _migrate_legacy(params)


# ---------------------------------------------------------------------------
# The legacy transform
# ---------------------------------------------------------------------------


@dataclass
class _Ledger:
    """The migration in progress: the dict being built, and what to tell the user."""

    arm: str
    params: dict
    seen: set = field(default_factory=set)
    dropped: list = field(default_factory=list)
    notes: list = field(default_factory=list)

    def drop(self, key, note=None):
        """Record that a key from the file has no home in the migrated dict."""
        self.seen.add(key)
        self.dropped.append(key)
        if note:
            self.notes.append(note)

    def take(self, key, target, value):
        """Record that a key was read, and store it under its contract name."""
        self.seen.add(key)
        self.params[target] = value

    def mark(self, key, note=None):
        """Record that a key was *considered* — read, but with nothing to store.

        Distinct from :meth:`drop`, which says the key's value did not survive, and from
        :meth:`take`, which stores one. Its whole job is to keep a key out of the unknown
        list at the end, and it is spelled out because ``take(key)`` with no target read
        like a typo.
        """
        self.seen.add(key)
        if note:
            self.notes.append(note)


def _migrate_legacy(raw) -> Migration:
    """Layer one legacy dict onto the pipeline parameter set for **its own** arm.

    Layering onto the contract is what replaces the block that back-filled every
    multiselect key with ``[]`` regardless of arm. That block was wrong twice over: ``[]``
    is not "unset" but "this source keeps nothing", so a file that simply never mentioned a
    source had that source's clause removed from the union; and doing it for all seven keys
    put three somatic keys in every germline export and two germline keys in every somatic
    one.

    A key the file does not carry now takes the pipeline's own default for that arm, which
    is what "the file said nothing about it" should mean in an app that opens at parity.
    """
    stated = raw.get("sample_type")
    arm = stated if stated in ARMS else "somatic"
    ledger = _Ledger(arm=arm, params=pipeline_params(arm))
    # Named rather than quietly corrected. Falling back to somatic is the app's behaviour
    # everywhere else, so a corrupt file opens rather than raising — but which arm filters
    # decides *which filter runs*, and a user whose file said something else has to be told
    # that this is not what they asked for.
    ledger.mark(
        "sample_type",
        None
        if stated in ARMS
        else f"sample_type was {stated!r}, which is not an arm; read as {arm}.",
    )

    _migrate_pathogenic_flag(raw, ledger)
    _migrate_genes(raw, ledger)
    _migrate_classification(raw, ledger)
    _migrate_guideline_lists(raw, ledger)
    _migrate_vaf(raw, ledger)
    _migrate_scalars(raw, ledger)

    unknown = sorted(set(raw) - ledger.seen)
    if unknown:
        ledger.dropped.extend(unknown)
        noun = "key" if len(unknown) == 1 else "keys"
        ledger.notes.append(
            f"Dropped {len(unknown)} {noun} this app has no parameter for: "
            + ", ".join(unknown)
            + "."
        )

    return Migration(
        params=ledger.params,
        migrated=True,
        dropped=tuple(ledger.dropped),
        notes=tuple(ledger.notes),
    )


def _migrate_pathogenic_flag(raw, ledger):
    """``keep_pathogenic`` is the app's polarity; ``skip_pathogenic`` is the pipeline's.

    Both spellings are read, and only the pipeline's is written. Leaving the old one
    behind is what made the retention checkbox dead for a while: the filter prefers
    ``skip_pathogenic`` whenever a dict carries it, so a stale ``keep_pathogenic`` sits in
    the file reading like a setting and deciding nothing.
    """
    if "keep_pathogenic" in raw:
        retained = bool(raw["keep_pathogenic"])
        ledger.take("keep_pathogenic", "skip_pathogenic", not retained)
        ledger.notes.append(
            f"keep_pathogenic={retained} became skip_pathogenic={not retained}: the "
            "pipeline states this flag the other way up."
        )
    if "skip_pathogenic" in raw:
        # The pipeline's own spelling wins where a file carries both, matching what the
        # filter itself does with such a dict.
        ledger.take("skip_pathogenic", "skip_pathogenic", bool(raw["skip_pathogenic"]))


def _migrate_genes(raw, ledger):
    """Two arm-prefixed string keys become one ``filter_genes`` list, picked by arm.

    A legacy dict can carry *both* prefixed keys, because the page wrote whichever arm the
    user was on and left the other behind — so picking without regard to arm applies one
    arm's panel to the other's report. Dropping both, which is what happens if the rename
    is missed, removes the gene restriction entirely: about sixteen extra rows per arm on
    the reference fixtures, silently.

    The value is *tokenised* rather than copied, through the same parser the filter seam
    uses. The legacy key held a string whose label promised commas while every gene file
    in the project is one symbol per line, and a single-element list flattened to a bare
    string — which the old preset processor did — is otherwise written one character per
    line at the gene adapter.
    """
    mine = f"{ledger.arm}_genes"
    theirs = "germline_genes" if ledger.arm == "somatic" else "somatic_genes"

    if raw.get("filter_genes") is not None:
        ledger.take(
            "filter_genes", "filter_genes", list(parse_gene_list(raw["filter_genes"]).symbols)
        )
        if mine in raw:
            ledger.drop(mine)
    elif raw.get(mine) is not None:
        symbols = list(parse_gene_list(raw[mine]).symbols)
        ledger.take(mine, "filter_genes", symbols)
        ledger.notes.append(
            f"{mine} became filter_genes ({len(symbols)} "
            f"{'symbol' if len(symbols) == 1 else 'symbols'})."
        )

    if theirs in raw:
        ledger.drop(
            theirs,
            f"{theirs} is the other arm's gene list; this file is {ledger.arm}.",
        )

    for key in _PANEL_KEYS:
        if key in raw:
            ledger.drop(key)
    if any(key in raw for key in _PANEL_KEYS):
        ledger.notes.append(
            "The gene-panel name is no longer a filter parameter; the panel is recognised "
            "from the gene list itself."
        )


def _migrate_classification(raw, ledger):
    """The include list becomes the exclude list the pipeline's parameter has always been.

    Two orderings are possible and only one is right. The sentinel means "every
    classification", so it has to be read *before* the complement is taken: stripping it
    first would leave an empty include list, whose complement is the whole vocabulary —
    every classification excluded, an empty report, from a file that asked for everything.

    An empty include list is inverted faithfully, to the whole vocabulary. It was
    unreachable in a saved file — the old validator back-filled every empty selection to
    the sentinel — and an empty report is at least loud, where the widening it would take
    to avoid it is not.
    """
    if "filter_variant_classification" not in raw:
        return

    included = _as_terms(raw["filter_variant_classification"], strip_sentinel=False)
    if LEGACY_SENTINEL in included:
        excluded = []
    else:
        excluded = [c for c in LEGACY_VARIANT_CLASSIFICATIONS if c not in included]

    ledger.take("filter_variant_classification", "filter_variant_classification", excluded)
    ledger.notes.append(
        f"filter_variant_classification listed {len(included)} classifications to report; "
        f"it is now the {len(excluded)} to exclude."
    )


def _migrate_guideline_lists(raw, ledger):
    """The six keep-lists: shape repaired, spelling repaired, other arms' keys dropped.

    Which keys belong to an arm is read off the contract rather than restated here, so a
    source added to one arm cannot be forgotten in this file. ESCAT is the one worth
    naming: ``germline_filters()`` takes no ESCAT argument, and the app-side germline ESCAT
    clause was the single largest divergence on real data — 540 rows, 51% of all attributed
    divergence — so a germline file that keeps the key carries a value that reads like a
    setting and reaches nothing.

    Why the sentinel is stripped here and *expanded* in the classification list
    --------------------------------------------------------------------------
    Because the two parameters have opposite polarity, and the sentinel's legacy meaning —
    "this source places no restriction" — lands on opposite values as a result.
    ``filter_variant_classification`` is an **exclude** list, so "no restriction" is the
    empty list, and reading ``All`` as "every classification" produces exactly that.
    These are **keep** lists ORed together, so "no restriction" can only be spelled by
    listing every term the source can carry.

    That spelling was measured and rejected. It is what the pre-#36 app did, and it
    saturates the union: on ``somatic_reference.maf``, expanding the sentinel on
    ``filter_cancervar`` alone takes the criteria path from the contract's 19 to **58**, and
    passing rows from 20 to **59 of 82**; on ``germline_reference.maf``, ``filter_intervar``
    alone takes 27 → **61** and 31 → **65 of 94**. That is the signature issue #36 deleted
    the sentinel over — one source reading as unrestricted making the OR true for nearly
    every row, 68,364 somatic rows against a parity baseline of 411 at full scale.

    So the sentinel has no faithful translation, and the choice is between a silent ~3x
    widening and dropping the source's clause. Dropping it is the loud option: the page
    warns when every source is empty, and the report visibly collapses to the
    pathogenic-rescue floor rather than quietly acquiring rows that look like findings.
    The parameter-migration prototype deleted in #52 recorded the expanding rule; it was
    written before #36 took that measurement, and this is where it is superseded.

    Loud, not silent, is the whole of the argument — so a stripped sentinel is *named*.
    """
    for key in _GUIDELINE_KEYS:
        if key not in raw:
            continue
        if key not in PIPELINE_PARAMS[ledger.arm]:
            ledger.drop(
                key,
                f"{key} is not a {ledger.arm} guideline source; the {ledger.arm} filter "
                "takes no such argument.",
            )
            continue
        terms = _as_terms(raw[key])
        repairs = _SPACING_REPAIRS.get(key, {})
        ledger.take(key, key, [repairs.get(term, term) for term in terms])
        if LEGACY_SENTINEL in _as_terms(raw[key], strip_sentinel=False):
            ledger.notes.append(
                f"{key} carried the old {LEGACY_SENTINEL!r} catch-all, which this app no "
                f"longer has: it is dropped, leaving {len(terms)} term(s). Honouring it "
                "would make this one source keep nearly every variant in the file."
            )


def _migrate_vaf(raw, ledger):
    """One threshold per arm, under the spelling that arm's widget and filter read.

    The contract carries a single ``vaf_threshold`` whose value depends on the arm, but
    the germline widget still writes ``vaf_threshold_germline`` and the filter still
    prefers it — so the germline arm gets both keys, set to the same number. Setting only
    one leaves the other to decide the report: a stale somatic 0.05 under the contract's
    name, or the widget's own fallback.
    """
    if ledger.arm == "somatic":
        if "vaf_threshold" in raw:
            ledger.take("vaf_threshold", "vaf_threshold", raw["vaf_threshold"])
        if "vaf_threshold_germline" in raw:
            ledger.drop(
                "vaf_threshold_germline",
                "vaf_threshold_germline is the germline threshold; this file is somatic.",
            )
        return

    if "vaf_threshold_germline" in raw:
        threshold = raw["vaf_threshold_germline"]
        ledger.take("vaf_threshold_germline", "vaf_threshold_germline", threshold)
        ledger.params["vaf_threshold"] = threshold
        if "vaf_threshold" in raw:
            ledger.mark(
                "vaf_threshold",
                None
                if raw["vaf_threshold"] == threshold
                else (
                    f"The germline arm uses vaf_threshold_germline={threshold}; the "
                    f"file's vaf_threshold={raw['vaf_threshold']} was the somatic tab's "
                    "leftover."
                ),
            )
    elif "vaf_threshold" in raw:
        ledger.take("vaf_threshold", "vaf_threshold", raw["vaf_threshold"])
        ledger.params["vaf_threshold_germline"] = raw["vaf_threshold"]


def _migrate_scalars(raw, ledger):
    """The keys whose meaning never moved, and the one whose name did."""
    for key in _UNCHANGED_SCALARS:
        if key in raw:
            value = bool(raw[key]) if key == "skip_civic" else raw[key]
            ledger.take(key, key, value)

    for old, new in _RENAMED_SCALARS.items():
        if old in raw:
            ledger.take(old, new, raw[old])
            ledger.notes.append(f"{old} is now {new}.")


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------


def _as_terms(value, strip_sentinel: bool = True) -> list:
    """One multiselect parameter's value, as a list of terms.

    A hand-edited file can hold a bare string where a list belongs, and anything else is
    unusable — repaired to empty rather than to a selection the user never made.
    """
    if isinstance(value, str):
        value = [value]
    if not isinstance(value, list):
        return []
    if not strip_sentinel:
        return list(value)
    return [term for term in value if term != LEGACY_SENTINEL]
