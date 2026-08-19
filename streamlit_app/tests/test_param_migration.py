"""A saved parameter file still means what it meant, and the cache is discarded (issue #40).

Two artefacts survive a MAFigate upgrade: a parameter file the user downloaded, and the
``~/.mafigate`` cache the app wrote for itself. Issues #31–#36 changed what several
parameter *keys* mean, so neither can simply be handed to the current filter. This module
is the test for what happens to each.

Why the gate is a version stamp and not a shape test
----------------------------------------------------
The obvious design — sniff the dict and migrate it if it "looks legacy" — cannot work,
and the reason is worth stating because it is not obvious. The parity default for
``filter_variant_classification`` is ``["Silent", "IGR", "RNA"]``. That is legal as a
post-parity *exclude* list (drop those three) and equally legal as a legacy *include*
list (report only those three). Nothing about membership, length or vocabulary
discriminates them.

And the migration of that key is a set complement, which is an *involution*: migrating
twice returns the original value while meaning the opposite. So a shape test that guesses
wrong is not merely imprecise, it silently produces a report built from the complement of
what the user asked for. Measured on ``somatic_reference.maf``: the parity value passes 19
rows on the criteria path, the twice-migrated value passes 2.

Hence :data:`PARAM_SCHEMA_VERSION`, and hence the second decision — the *cache* is
discarded rather than migrated. Every cache ever written carries the literal
``"app_version": "2.0.0"``, which the app has never bumped, so an existing cache and an
unstamped one are the same state and the stamp has no power over either. A parameter file
the user can re-upload is a different matter: they still have it, and migrating it is
their explicit act.

What is asserted here
---------------------
* the gate: a stamped document is a no-op, an unstamped one is migrated, and the two are
  told apart by the stamp alone;
* the key mapping, key by key, because every per-key mistake **widens** the report
  silently — an unrenamed gene key is worth about sixteen extra rows per arm and an
  uninverted pathogenic flag two more, and the two compound;
* that migration lands on the pipeline parameter set *for the uploaded dict's own arm*,
  so a germline file cannot arrive carrying somatic guideline sources;
* that the cache is discarded with a banner that names what was dropped;
* that the dead code this ticket removes is actually gone.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

STREAMLIT_APP = Path(__file__).resolve().parent.parent
if str(STREAMLIT_APP) not in sys.path:
    sys.path.insert(0, str(STREAMLIT_APP))

from config.constants import APP_VERSION  # noqa: E402
from config.vocabularies import VARIANT_CLASSIFICATIONS  # noqa: E402
from config.param_migration import (  # noqa: E402
    LEGACY_VARIANT_CLASSIFICATIONS,
    PARAM_SCHEMA_VERSION,
    SCHEMA_VERSION_KEY,
    migrate_params,
    param_document,
    unwrap_document,
)
from config.pipeline_params import ARMS, PIPELINE_PARAMS, pipeline_params  # noqa: E402
from filters.variant_filters import apply_filters  # noqa: E402
from vendor.pipeline_utils import read_maf  # noqa: E402

FIXTURE_DIR = STREAMLIT_APP / "tests" / "fixtures" / "parity"

#: The catch-all a pre-#36 file can still carry.
SENTINEL = "All"


class FakeSessionState(dict):
    """A dict that also answers to attribute access, as Streamlit's session state does.

    The page reaches its parameters both ways — ``st.session_state.filter_params`` in most
    places, ``st.session_state["discarded_cache"]`` where the name is computed — so a plain
    dict stubs only half of it and fails on the other half.
    """

    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError as missing:
            raise AttributeError(name) from missing

    def __setattr__(self, name, value):
        self[name] = value


@pytest.fixture(scope="module")
def reference_mafs():
    """Both reference fixtures, loaded by the pipeline's own reader."""
    return {
        "somatic": read_maf(str(FIXTURE_DIR / "somatic_reference.maf")),
        "germline": read_maf(str(FIXTURE_DIR / "germline_reference.maf")),
    }


def _legacy_somatic():
    """A parameter file the app would have written before the parity work.

    Deliberately carries the traps rather than a clean dict: an include-style
    classification list, ``keep_pathogenic``, both arms' gene keys, a panel name, the
    catch-all sentinel, and a key nothing has ever read.
    """
    return {
        "sample_type": "somatic",
        "min_depth": 30,
        "vaf_threshold": 0.05,
        "vaf_threshold_germline": 0.3,
        "somatic_genes": "TP53, KRAS\nEGFR",
        "germline_genes": "BRCA1, BRCA2",
        "somatic_gene_set": "Custom",
        "germline_gene_set": "All",
        "filter_variant_classification": ["Missense_Mutation", "Nonsense_Mutation"],
        "filter_cancervar": ["Tier_I_Strong", "Tier_II_Potential"],
        "filter_civic": ["A", "B"],
        "filter_clinvar": ["Pathogenic", "Likely pathogenic"],
        "filter_escat": [SENTINEL],
        "filter_intervar": ["Pathogenic"],
        "filter_renovo": ["HP Pathogenic"],
        "max_freq_population": 0.01,
        "keep_pathogenic": True,
        "skip_civic": True,
        "notes_for_myself": "please keep",
    }


def _legacy_germline():
    """The germline counterpart, carrying the somatic arm's keys as well."""
    return {
        "sample_type": "germline",
        "min_depth": 30,
        "vaf_threshold": 0.05,
        "vaf_threshold_germline": 0.25,
        "somatic_genes": "TP53",
        "germline_genes": "BRCA1, BRCA2, ATM",
        "filter_variant_classification": ["Missense_Mutation"],
        "filter_cancervar": ["Tier_I_Strong"],
        "filter_civic": ["A"],
        "filter_escat": ["IA", "IB"],
        "filter_intervar": ["Likely_pathogenic"],
        "filter_renovo": ["Pathogenic"],
        "filter_clinvar": ["Pathogenic"],
        "keep_pathogenic": False,
    }


# ---------------------------------------------------------------------------
# The gate is a version stamp
# ---------------------------------------------------------------------------


def test_the_schema_version_is_not_the_app_version():
    """A constant of its own, so the format can move without a release and vice versa.

    They were the same field once — the cache's ``app_version`` — and that is precisely
    why the cache cannot be migrated: the literal ``"2.0.0"`` has never been bumped in any
    commit that touched it, so it never distinguished one format from another. A schema
    version is bumped by the change that alters the format, which is a different event
    from a release.
    """
    assert isinstance(PARAM_SCHEMA_VERSION, int)
    assert PARAM_SCHEMA_VERSION >= 1
    assert PARAM_SCHEMA_VERSION != APP_VERSION

    source = (STREAMLIT_APP / "config" / "param_migration.py").read_text()
    stamp = next(
        line for line in source.splitlines() if line.startswith("PARAM_SCHEMA_VERSION")
    )
    assert "APP_VERSION" not in stamp, (
        f"the schema version is derived from the app version: {stamp!r}. Tying them "
        "together means a release renumbers the format and a format change needs a "
        "release."
    )


def test_a_document_carries_the_stamp_and_the_parameters():
    """The envelope the cache already had, plus the stamp, for exports too."""
    params = pipeline_params("somatic")
    document = param_document(params)

    assert document[SCHEMA_VERSION_KEY] == PARAM_SCHEMA_VERSION
    assert document["parameters"] == params
    assert document["app_version"] == APP_VERSION
    assert document["timestamp"]

    # It survives the trip through JSON that an export and an upload make.
    assert json.loads(json.dumps(document)) == document


def test_a_bare_dict_and_an_envelope_both_unwrap():
    """Every file written before this ticket is a bare parameter dict, not an envelope."""
    params = {"sample_type": "somatic", "min_depth": 50}

    unwrapped, version = unwrap_document(params)
    assert unwrapped == params
    assert version is None, "a bare dict must not be read as stamped"

    unwrapped, version = unwrap_document(param_document(params))
    assert unwrapped == params
    assert version == PARAM_SCHEMA_VERSION


@pytest.mark.parametrize("arm", ARMS)
def test_re_uploading_an_already_migrated_file_is_a_no_op(arm):
    """The criterion, and the thing a shape test cannot deliver.

    A user who exports their parameters and imports them again — or imports the same file
    into two sessions — must get back exactly what they exported. Under a shape test they
    would get the complement of their classification list, silently.
    """
    exported = param_document(pipeline_params(arm))

    result = migrate_params(exported)

    assert result.migrated is False
    assert result.params == PIPELINE_PARAMS[arm]
    assert result.dropped == ()

    # And again, because "idempotent" has to survive more than one round trip.
    again = migrate_params(param_document(result.params))
    assert again.params == result.params


def test_the_parity_classification_list_is_not_inverted_a_second_time(reference_mafs):
    """The involution, asserted on rows rather than on the list.

    ``["Silent", "IGR", "RNA"]`` is the contract's exclude list *and* a legal legacy
    include list. Migrating it a second time returns a 15-value list which, read as an
    exclude list, drops everything the report is made of. This is what the stamp buys, and
    it is quoted on the criteria path because the union count hides it: pathogenic
    retention re-admits most of what the wrong list dropped.
    """
    maf = reference_mafs["somatic"]

    stamped = migrate_params(param_document(pipeline_params("somatic"))).params
    _, correct = apply_filters(maf, stamped)

    unstamped = migrate_params(dict(pipeline_params("somatic"))).params
    _, twice = apply_filters(maf, unstamped)

    assert correct.criteria_path == 54
    assert twice.criteria_path < correct.criteria_path, (
        "migrating the parity value a second time changed nothing, so this fixture no "
        "longer demonstrates why the stamp is needed"
    )
    assert stamped["filter_variant_classification"] == ["Silent", "IGR", "RNA"]


def test_a_document_from_a_newer_schema_is_left_alone_and_named():
    """We cannot know what a future format means, so we do not guess at it."""
    document = param_document({"sample_type": "somatic", "min_depth": 7})
    document[SCHEMA_VERSION_KEY] = PARAM_SCHEMA_VERSION + 5

    result = migrate_params(document)

    assert result.migrated is False
    assert result.params["min_depth"] == 7
    assert any("newer" in note.lower() for note in result.notes), result.notes


# ---------------------------------------------------------------------------
# The key mapping
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("kept,expected", [(True, False), (False, True)])
def test_the_pathogenic_flag_is_inverted(kept, expected):
    """``keep_pathogenic`` is the app's polarity; ``skip_pathogenic`` is the pipeline's.

    Worth +2 rows on each arm if it survives unconverted, and worth the whole setting if
    it survives *alongside* the new key: the filter prefers ``skip_pathogenic`` whenever a
    dict carries it, so a stale ``keep_pathogenic`` is a value that reads like a setting
    and decides nothing.
    """
    result = migrate_params({"sample_type": "somatic", "keep_pathogenic": kept})

    assert result.params["skip_pathogenic"] is expected
    assert "keep_pathogenic" not in result.params
    assert any("pathogenic" in note.lower() for note in result.notes), result.notes


def test_an_absent_pathogenic_flag_takes_the_contracts_value():
    """Absent is not False. It is "the file does not say", which is the contract."""
    result = migrate_params({"sample_type": "somatic"})
    assert result.params["skip_pathogenic"] == PIPELINE_PARAMS["somatic"]["skip_pathogenic"]


@pytest.mark.parametrize(
    "arm,expected",
    [("somatic", ["TP53", "KRAS", "EGFR"]), ("germline", ["BRCA1", "BRCA2", "ATM"])],
)
def test_the_gene_keys_are_collapsed_and_picked_by_arm(arm, expected):
    """Two prefixed keys become one ``filter_genes``, chosen by the dict's own arm.

    A legacy dict can carry *both*, because the page wrote whichever arm the user was on
    and left the other behind. Picking the wrong one applies another arm's panel to this
    arm's report; keeping neither drops the restriction entirely, which is the silent
    widening — about sixteen extra rows per arm on the reference fixtures.

    The value is tokenised rather than copied: the legacy key held a string, the contract
    holds a list of symbols, and the separator is not only the comma the old label
    promised.
    """
    legacy = _legacy_somatic() if arm == "somatic" else _legacy_germline()

    result = migrate_params(legacy)

    assert result.params["filter_genes"] == expected
    assert "somatic_genes" not in result.params
    assert "germline_genes" not in result.params


def test_a_gene_list_already_collapsed_survives_migration():
    """A file written between #34 and this ticket carries ``filter_genes`` already."""
    result = migrate_params({"sample_type": "somatic", "filter_genes": ["TP53", "KRAS"]})
    assert result.params["filter_genes"] == ["TP53", "KRAS"]


def test_a_gene_list_flattened_to_a_string_is_read_as_one_symbol():
    """``["TP53"]`` collapsed to ``"TP53"`` must not become four one-letter genes.

    The old preset processor collapsed every single-element list to its element. At the
    gene adapter that string is written one character per line, so the filter restricts to
    genes ``T``, ``P``, ``5`` and ``3`` — which match nothing, emptying the report.
    """
    result = migrate_params({"sample_type": "somatic", "filter_genes": "TP53"})
    assert result.params["filter_genes"] == ["TP53"]


def test_the_panel_name_is_dropped_and_named():
    """The panel choice is UI state, recovered from the symbols (issue #34).

    Dropping it is safe because the page matches an incoming gene list back to the panel
    that denotes it, so a saved panel choice survives in the only form that ever filtered
    anything.
    """
    result = migrate_params(_legacy_somatic())

    assert "somatic_gene_set" not in result.params
    assert "germline_gene_set" not in result.params
    assert "gene_set" not in result.params
    assert "somatic_gene_set" in result.dropped


def test_the_escat_key_is_dropped_for_germline():
    """``germline_filters()`` takes no ESCAT argument.

    The app-side germline ESCAT clause was the single largest divergence on real data —
    540 rows, 51% of all attributed divergence — and issue #33 deleted it. A germline file
    that still carries the key would put a value in every export that nothing reads.
    """
    result = migrate_params(_legacy_germline())

    assert "filter_escat" not in result.params
    assert "filter_escat" in result.dropped
    assert any("escat" in note.lower() for note in result.notes), result.notes


def test_the_escat_key_survives_for_somatic():
    """It is a real guideline source on the arm whose filter takes it."""
    legacy = _legacy_somatic()
    legacy["filter_escat"] = ["IA", "IB"]

    result = migrate_params(legacy)

    assert result.params["filter_escat"] == ["IA", "IB"]


@pytest.mark.parametrize("skip", [True, False])
def test_the_civic_skip_flag_is_preserved(skip):
    """Preserved, not deleted — it was on a "deprecated parameters" list by mistake.

    ``skip_civic`` is a live pipeline parameter on both arms: it decides the output column
    set everywhere and, on somatic, which rows pass. Deleting it silently reset a user's
    choice to the contract's ``False``.
    """
    result = migrate_params({"sample_type": "somatic", "skip_civic": skip})

    assert result.params["skip_civic"] is skip
    assert "skip_civic" not in result.dropped


def test_unknown_keys_are_dropped_and_named():
    """Dropped so they cannot be mistaken for settings; named so nothing vanishes quietly."""
    result = migrate_params(
        {"sample_type": "somatic", "notes_for_myself": "please keep", "gene_filtering": {}}
    )

    assert "notes_for_myself" not in result.params
    assert "gene_filtering" not in result.params
    assert set(result.dropped) >= {"notes_for_myself", "gene_filtering"}
    joined = " ".join(result.notes)
    assert "notes_for_myself" in joined and "gene_filtering" in joined, result.notes


# ---------------------------------------------------------------------------
# The classification list, inverted over a frozen vocabulary
# ---------------------------------------------------------------------------


def test_the_classification_list_is_inverted_from_include_to_exclude():
    """The legacy key listed what to report; the contract's lists what to drop."""
    result = migrate_params(
        {
            "sample_type": "somatic",
            "filter_variant_classification": ["Missense_Mutation", "Nonsense_Mutation"],
        }
    )

    excluded = result.params["filter_variant_classification"]
    assert "Missense_Mutation" not in excluded
    assert "Nonsense_Mutation" not in excluded
    assert set(excluded) == set(LEGACY_VARIANT_CLASSIFICATIONS) - {
        "Missense_Mutation",
        "Nonsense_Mutation",
    }


def test_the_sentinel_classification_list_excludes_nothing():
    """``["All"]`` meant "report every classification", which is an empty exclude list.

    Stripping the sentinel first and *then* inverting would produce the complement of
    nothing — every classification excluded, an empty report. The order matters and this
    is what pins it.
    """
    for value in ([SENTINEL], [SENTINEL, "Missense_Mutation"], SENTINEL):
        result = migrate_params(
            {"sample_type": "somatic", "filter_variant_classification": value}
        )
        assert result.params["filter_variant_classification"] == [], value


def test_the_inversion_runs_over_a_frozen_vocabulary():
    """Adding an option to the live control must not change what an old file means.

    The complement is only well defined against the vocabulary the *writer* could choose
    from. Taken against today's list, a classification added next year would start being
    excluded by every file written before it existed — a silent narrowing, retroactive,
    and impossible to see in the file.

    Asserted on the module's own syntax tree rather than by monkeypatching
    ``VARIANT_CLASSIFICATIONS``: the migration does not read that constant, so a
    monkeypatch would prove nothing and the test could never fail. What can actually go
    wrong is someone wiring the live list in here — reasonably enough, to "keep them in
    step" — so that is what is checked.

    The tree, not the text: the module *discusses* the live constant at length, and a
    substring check over the source flags its own explanation of why it must not use it.
    """
    import ast

    tree = ast.parse((STREAMLIT_APP / "config" / "param_migration.py").read_text())
    referenced = {
        node.id if isinstance(node, ast.Name) else node.attr
        for node in ast.walk(tree)
        if isinstance(node, (ast.Name, ast.Attribute))
    }
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom):
            referenced.update(alias.name for alias in node.names)

    assert "VARIANT_CLASSIFICATIONS" not in referenced, (
        "config/param_migration.py now refers to the live VARIANT_CLASSIFICATIONS. The "
        "complement must be taken over the frozen copy, or a classification added later "
        "retroactively changes what every old file means."
    )

    excluded = migrate_params(
        {"sample_type": "somatic", "filter_variant_classification": ["Missense_Mutation"]}
    ).params["filter_variant_classification"]

    # The complement is over the frozen vocabulary exactly — no more, no fewer.
    assert set(excluded) == set(LEGACY_VARIANT_CLASSIFICATIONS) - {"Missense_Mutation"}
    assert len(LEGACY_VARIANT_CLASSIFICATIONS) == 18


def test_the_migration_knows_every_guideline_source_the_page_renders():
    """The second copy, pinned — the remedy ``pipeline_params.ARMS`` sets for its own.

    ``config/param_migration.py`` cannot import ``config.vocabularies.GUIDELINE_SOURCES``
    without closing an import cycle — the vocabulary module takes ``LEGACY_SENTINEL`` from
    the migration module — so it restates the six keys. Unpinned, that copy is how a
    seventh guideline source would be added to the contract and then silently dropped from
    every uploaded file as an "unknown key".
    """
    from config.param_migration import _GUIDELINE_KEYS
    from config.vocabularies import GUIDELINE_SOURCES

    rendered = {key for keys in GUIDELINE_SOURCES.values() for key in keys}
    assert set(_GUIDELINE_KEYS) == rendered, (
        "the migration's guideline keys and the page's guideline sources disagree: "
        f"{sorted(set(_GUIDELINE_KEYS) ^ rendered)}"
    )


def test_the_frozen_vocabulary_matches_the_live_one_today():
    """A record, not a rule — and the distinction is the point.

    The two coincide right now, and this says so, so a reader is not left wondering
    whether the frozen copy was already stale when it was taken. It is deliberately an
    *equality of contents today* rather than a constraint: when a classification is added
    to the control, this assertion is the one that should be updated to record the
    divergence, and :func:`test_the_inversion_runs_over_a_frozen_vocabulary` is the one
    that must keep passing untouched.
    """
    assert list(LEGACY_VARIANT_CLASSIFICATIONS) == VARIANT_CLASSIFICATIONS


def test_an_absent_classification_list_takes_the_contracts_three_exclusions():
    result = migrate_params({"sample_type": "somatic"})
    assert result.params["filter_variant_classification"] == ["Silent", "IGR", "RNA"]


# ---------------------------------------------------------------------------
# Per-source spacing, which is repaired rather than normalised
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "key,legacy,expected",
    [
        ("filter_cancervar", ["Tier_I_Strong"], ["Tier_I_strong"]),
        ("filter_cancervar", ["Tier_III_Unknown"], ["Tier_III_Uncertain"]),
        ("filter_intervar", ["Likely_pathogenic"], ["Likely pathogenic"]),
        ("filter_intervar", ["Uncertain_significance"], ["Uncertain significance"]),
        ("filter_clinvar", ["Likely pathogenic"], ["Likely_pathogenic"]),
        ("filter_clinvar", ["Drug response"], ["drug_response"]),
        (
            "filter_clinvar",
            ["Conflicting interpretations of pathogenicity"],
            ["Conflicting_interpretations_of_pathogenicity"],
        ),
        ("filter_renovo", ["Pathogenic"], ["HP Pathogenic"]),
    ],
)
def test_the_per_source_spacing_repair_survives(key, legacy, expected):
    """Two sources spaced, two underscored — no global normalisation.

    InterVar and ReNovo values carry spaces; CancerVar and ClinVar values carry
    underscores. That is how the annotations are written, so a repair table is exactly
    what the contract requires and a single normalisation rule would break two sources
    whichever way it went.
    """
    arm = "germline" if key in ("filter_intervar", "filter_renovo") else "somatic"
    result = migrate_params({"sample_type": arm, key: legacy})
    assert result.params[key] == expected


@pytest.mark.parametrize(
    "saved",
    [
        "Pathogenic/Likely_pathogenic",
        "Benign/Likely_benign",
        "Pathogenic/Likely pathogenic",
        "Benign/Likely benign",
    ],
)
def test_a_saved_composite_clinvar_term_survives_migration_without_being_reinterpreted(saved):
    """An older file's dead composite is carried, not repaired into an atom (issue #88).

    Two repairs used to map the spaced composites onto underscored ones. #88 took the
    underscored ones off the vocabulary, which left those repairs pointing at a keep-term
    the app no longer has — so both went, rather than being re-pointed at ``Pathogenic``
    or ``Benign``. Re-pointing is the tempting move and it is the wrong one: it would be
    this module deciding which half of a two-call selection the user meant, and it would
    *widen* their report, silently, during a format migration.

    What matters is that nothing is invented. The term passes through, and the widget
    drops it from its default because ``filter_terms`` keeps only what the control offers
    — which costs the user nothing, the term never having matched a row in its life.
    """
    result = migrate_params({"sample_type": "somatic", "filter_clinvar": [saved]})

    # Carried verbatim, asserted directly. An earlier version of this test paired the
    # composite with a plain "Pathogenic" and checked that every surviving term was one of
    # the two — which a repair re-pointing the composite *at* "Pathogenic" passes, since
    # the answer it produces is on the allowed list. Mutation caught it; the honest
    # property is that the term itself comes out the far side untouched.
    assert result.params["filter_clinvar"] == [saved], (
        f"a saved {saved!r} came back as {result.params['filter_clinvar']!r}; a composite "
        "is carried unchanged or dropped, never resolved into one of its halves — "
        "resolving it would silently widen a report during a format migration"
    )


def test_a_saved_composite_is_dropped_at_the_widget_rather_than_by_the_migration():
    """The two halves of the same guarantee, checked together.

    The migration carries the term (above) and the widget refuses it, so a saved file
    written by an older MAFigate opens without an exception and without a control showing
    a value that cannot do anything. Streamlit raises on a default outside its options, so
    this is the difference between a degraded parameter and a page that will not render.
    """
    from config.vocabularies import CLINVAR_OPTIONS, filter_terms

    carried = migrate_params(
        {
            "sample_type": "somatic",
            "filter_clinvar": ["Pathogenic/Likely_pathogenic", "Pathogenic"],
        }
    ).params["filter_clinvar"]

    assert filter_terms(carried, CLINVAR_OPTIONS) == ["Pathogenic"], (
        "the widget default must keep only what the control offers; a composite reaching "
        "st.multiselect as a default raises"
    )


def test_the_sentinel_is_stripped_from_a_keep_list_rather_than_expanded():
    """In a keep list ``All`` and ``[]`` already decide every row identically (issue #36)."""
    result = migrate_params({"sample_type": "somatic", "filter_clinvar": ["Pathogenic", SENTINEL]})
    assert result.params["filter_clinvar"] == ["Pathogenic"]


def test_a_stripped_sentinel_is_named():
    """Dropping it narrows the report, so it cannot also be silent."""
    result = migrate_params({"sample_type": "somatic", "filter_cancervar": [SENTINEL]})

    assert result.params["filter_cancervar"] == []
    assert any("filter_cancervar" in note and SENTINEL in note for note in result.notes), (
        result.notes
    )


@pytest.mark.parametrize(
    "arm,source,expanded_criteria,expanded_passed",
    [("somatic", "filter_cancervar", 67, 73), ("germline", "filter_intervar", 68, 76)],
)
def test_expanding_the_sentinel_would_saturate_the_guideline_union(
    arm, source, expanded_criteria, expanded_passed, reference_mafs
):
    """The rejected alternative, measured — so it is refused rather than merely not chosen.

    The parameter-migration prototype deleted in #52 recorded the other rule: read a legacy
    ``["All"]`` as "every option on offer", which is the only post-parity spelling of the
    "no restriction" the sentinel used to mean. It was written before issue #36 measured
    what that costs.

    One source is enough. Expanding it takes the report to most of the file, which is the
    exact signature #36 deleted the sentinel over — 68,364 somatic rows against a parity
    baseline of 411 at full scale. The strip narrows instead, and a narrowing is the loud
    direction: the page warns on the all-empty state and the report visibly collapses,
    where these rows would look like findings.

    **WHAT SATURATION IS ASSERTED AS, AND WHY THAT CHANGED AT ISSUE #246.** This used to
    read ``expanded.passed > 2 * contract.passed`` — "roughly three times the contract's" —
    which was a proxy for saturation that worked only because the fixture it was written
    against was a purposive subsample of real rows: 62 of its 82 rows carried no witness at
    all and failed everything, so the contract baseline was 20 and there was room to
    triple. The constructed set that replaced it packs witnesses into the same 82 rows, so
    its contract baseline is 60 and **no amount of saturation can be three times it** — 82
    rows is the ceiling. The multiplier was measuring witness density, not saturation.

    So the claim is now asserted as the *mechanism* it always meant: with every option
    selected the criteria path reaches **exactly** the rows that clear depth,
    classification and VAF, which is the guideline block excluding nothing at all. That is
    strictly what #36 objected to, and it cannot be satisfied by an accident of baseline
    density. Recorded plainly rather than quietly: this is a *different* assertion, not the
    old one with a smaller number in it.

    The second half is what keeps the first from going vacuous. "Everything reachable is
    reported" would also be true of a fixture whose guideline block never excluded anything
    in the first place, so the contract's own criteria path is asserted to fall short of
    the same bound — 54 of 67 somatic, 51 of 68 germline. An earlier draft of this checked
    the row *indices* of the returned frame instead, which is always empty: ``apply_filters``
    returns the annotated frame rather than the report, so the difference could not be
    anything else.

    **What is no longer demonstrable here, stated rather than lost.** The old assertion
    also showed the *report* ballooning — the visible symptom #36 deleted the sentinel over,
    68,364 somatic rows against a parity baseline of 411 at full scale. On this fixture the
    report goes 60 → 73, which is a widening but not a ballooning, and no fixture of 82 rows
    with this witness density can show one. `expanded_passed` is still recorded above so the
    movement is pinned; the *magnitude* of the real-scale symptom now lives only in #36's
    own measurement, which is where it was measured.
    """
    from config.vocabularies import CANCERVAR_OPTIONS, INTERVAR_OPTIONS
    from vendor.pipeline_filters import common_filters

    options = {"filter_cancervar": CANCERVAR_OPTIONS, "filter_intervar": INTERVAR_OPTIONS}[source]
    maf = reference_mafs[arm]
    params = pipeline_params(arm)

    _, contract = apply_filters(maf, params)
    _, expanded = apply_filters(maf, dict(params, **{source: list(options)}))
    _, stripped = apply_filters(
        maf, migrate_params({"sample_type": arm, source: [SENTINEL]}).params
    )

    assert (expanded.criteria_path, expanded.passed) == (expanded_criteria, expanded_passed)

    # The rows that clear depth, classification and VAF — everything the criteria path can
    # reach before any guideline block has an opinion. One ``vaf_threshold`` per arm:
    # ``pipeline_params`` has already resolved the config's two keys to this arm's.
    reachable = int(
        (
            common_filters(maf, params["min_depth"], params["filter_variant_classification"])
            & (maf["tumor_f"] > params["vaf_threshold"])
        ).sum()
    )
    assert expanded.criteria_path == reachable, (
        f"{arm}: with every {source} option selected the criteria path reaches "
        f"{expanded.criteria_path} of the {reachable} rows that clear depth, "
        "classification and VAF — so the guideline block is still excluding something and "
        "this no longer demonstrates why the sentinel is stripped"
    )
    assert contract.criteria_path < reachable, (
        f"{arm}: the contract's own guideline block already excludes nothing on this "
        "fixture, so the assertion above would hold however the sentinel were read"
    )
    assert stripped.passed <= contract.passed, (
        f"{arm}: stripping the sentinel widened the report past the contract's"
    )


# ---------------------------------------------------------------------------
# Layering onto the pipeline parameter set, for the file's own arm
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("arm", ARMS)
def test_an_absent_key_takes_the_contracts_value_rather_than_an_empty_list(arm):
    """The back-fill this replaces wrote ``[]``, which is not "unset" but "keeps nothing".

    A file that never mentioned a guideline source said nothing about it, and the app's
    answer to "nothing said" is the pipeline's own default. Back-filling ``[]`` instead
    removed that source's clause from the union — a narrowing the user never chose and the
    all-empty warning cannot see, because it only fires when *every* source is empty.
    """
    result = migrate_params({"sample_type": arm, "min_depth": 30})

    for key, value in PIPELINE_PARAMS[arm].items():
        if key == "min_depth":
            continue
        assert result.params[key] == value, key


@pytest.mark.parametrize(
    "arm,foreign",
    [
        ("germline", ("filter_cancervar", "filter_civic", "filter_escat")),
        ("somatic", ("filter_intervar", "filter_renovo")),
    ],
)
def test_the_other_arms_guideline_sources_do_not_survive(arm, foreign):
    """Layering is per arm, so a file cannot import the other arm's sources.

    The block this replaces back-filled all seven multiselect keys regardless of arm, so
    every germline export carried three somatic keys and every somatic export two germline
    ones — values that read as settings, that the arm's filter never passes on, and that
    the next reader has to know to ignore.
    """
    legacy = _legacy_somatic() if arm == "somatic" else _legacy_germline()

    result = migrate_params(legacy)

    for key in foreign:
        assert key not in result.params, key
        assert key in result.dropped, key


@pytest.mark.parametrize("arm", ARMS)
def test_migration_lands_on_exactly_the_arms_key_set(arm):
    """No key the arm has no use for, and none of the arm's keys missing."""
    legacy = _legacy_somatic() if arm == "somatic" else _legacy_germline()

    result = migrate_params(legacy)

    expected = set(PIPELINE_PARAMS[arm])
    if arm == "germline":
        # The page's germline VAF widget still writes this spelling; the contract's
        # single-key collapse is not at the widget yet.
        expected.add("vaf_threshold_germline")
    assert set(result.params) == expected


def test_an_unknown_arm_is_read_as_somatic_and_said_so():
    """The app's own fallback everywhere else, so a corrupt file opens rather than raises.

    Named rather than silently corrected. ``sample_type`` decides *which filter runs*, not
    how one is tuned, so a file that said something else has not been narrowed or widened —
    it has been answered with a different question, and the one rule this migration keeps
    everywhere is that nothing about it is silent.
    """
    result = migrate_params({"sample_type": "elephant", "min_depth": 30})

    assert result.params["sample_type"] == "somatic"
    assert any("elephant" in note for note in result.notes), result.notes


def test_a_stated_arm_is_not_remarked_on():
    """The note fires on the state it names, not on every file."""
    result = migrate_params({"sample_type": "germline", "min_depth": 30})
    assert not [note for note in result.notes if "sample_type" in note], result.notes


# ---------------------------------------------------------------------------
# The VAF thresholds, which the page still spells two ways
# ---------------------------------------------------------------------------


def test_the_somatic_arm_keeps_its_own_threshold_and_drops_the_germline_one():
    result = migrate_params(_legacy_somatic())

    assert result.params["vaf_threshold"] == 0.05
    assert "vaf_threshold_germline" not in result.params
    assert "vaf_threshold_germline" in result.dropped


def test_the_germline_arm_keeps_the_germline_threshold_under_both_spellings():
    """The page writes ``vaf_threshold_germline``; the contract carries ``vaf_threshold``.

    Both are set, to the same number, because the two live side by side until the widget
    is moved onto the contract's single key. Migrating to only one of them would let the
    other — a stale somatic 0.05, or the widget's own fallback — decide the germline arm.
    """
    result = migrate_params(_legacy_germline())

    assert result.params["vaf_threshold_germline"] == 0.25
    assert result.params["vaf_threshold"] == 0.25


def test_the_old_frequency_key_is_renamed():
    """``min_freq`` was the app's name for what is now ``max_freq_population``."""
    result = migrate_params({"sample_type": "somatic", "min_freq": 0.02})
    assert result.params["max_freq_population"] == 0.02
    assert "min_freq" not in result.params


# ---------------------------------------------------------------------------
# What the migration is worth, in rows
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("arm", ARMS)
def test_a_migrated_gene_restriction_actually_reaches_the_filter(arm, reference_mafs):
    """The per-key damage, measured rather than asserted from the dict.

    An unrenamed gene key leaves ``filter_genes`` at the contract's empty list, which
    means *no gene restriction at all* — the file asked for three genes and got the whole
    panel.

    Quoted on the criteria path rather than on the union, and this fixture is why: on
    ``somatic_reference.maf`` both settings pass **18** rows. The gene restriction does not
    shrink the union, it moves rows across it — the criteria path falls from 13 to 3 while
    pathogenic rescue re-admits exactly the ten rows it dropped, so ``rescue_only`` climbs
    from 5 to 15. A test written on the total would have passed with the migration deleted.

    The issue quotes about sixteen extra rows per arm for this key; that is measured on the
    GERSOM reference, and the exact cells below are this fixture's own numbers.

    **The coincidence is built rather than found**, and issue #246 is where that changed.
    On the constructed fixture set the three symbols the legacy file names — ``TP53``,
    ``KRAS``, ``EGFR`` — carry rows placed here on purpose (``legacy_gene_rows`` in
    ``tests/fixtures/parity/build_parity_fixtures.py``), and every criteria-path row that
    is *not* on that list is pathogenic, so the rescue re-admits precisely what the
    restriction drops. On the real subsample this set replaced the same property held by
    accident of which genes the sampled rows happened to carry, which is why it broke on
    the first swap and why it is now a documented obligation of the generator rather than
    an observation about the data. If these cells move again, the thing to re-check is
    whether the union still coincides — that is the argument, the numbers are only its
    evidence.
    """
    maf = reference_mafs[arm]
    legacy = _legacy_somatic() if arm == "somatic" else _legacy_germline()

    migrated = migrate_params(legacy).params
    _, restricted = apply_filters(maf, migrated)

    unrestricted = dict(migrated, filter_genes=[])
    _, widened = apply_filters(maf, unrestricted)

    assert restricted.criteria_path < widened.criteria_path, (
        f"{arm}: the migrated gene list reached no row the unrestricted one did not, so "
        "this fixture cannot show what dropping the key costs"
    )

    if arm == "somatic":
        assert (restricted.criteria_path, widened.criteria_path) == (3, 13)
        assert restricted.passed == widened.passed == 18, (
            "the union no longer coincides — the docstring's warning about totals is now "
            "describing something this fixture does not do"
        )


def test_a_migrated_file_filters_identically_after_a_round_trip(reference_mafs):
    """Export, re-import, filter: the same verdict on every row.

    The whole point of the envelope. A user comparing two sessions must not find that the
    act of saving and reloading their parameters moved the report.
    """
    maf = reference_mafs["somatic"]
    migrated = migrate_params(_legacy_somatic()).params

    reloaded = migrate_params(json.loads(json.dumps(param_document(migrated)))).params

    first, _ = apply_filters(maf, migrated)
    second, _ = apply_filters(maf, reloaded)
    assert first.equals(second)


# ---------------------------------------------------------------------------
# The cache, discarded rather than migrated
# ---------------------------------------------------------------------------


@pytest.fixture
def cache_at(tmp_path):
    """Point the cache constants at a temporary directory, with streamlit patched out.

    Patched at the module rather than through ``HOME`` because the path is a module-level
    constant evaluated at import, and because a developer's own ``~/.mafigate`` deciding a
    test's outcome has already cost this project one confusing failure.
    """
    from page_modules import param_store

    directory = tmp_path / ".mafigate"
    directory.mkdir()
    cache_file = directory / "cached_parameters.json"

    with patch.object(param_store, "CACHE_DIR", directory), patch.object(
        param_store, "PARAMS_CACHE_FILE", cache_file
    ), patch.object(param_store, "st", MagicMock()):
        yield cache_file


def test_an_unstamped_cache_is_discarded_and_never_restored(cache_at):
    """The state the app opened in: a cache that silently defeated opening at parity.

    Every cache ever written is unstamped, carries whatever arm the user last looked at,
    and is restored before any page renders — so "the app opens at parity" was true only
    for a user who had never opened the app before.
    """
    from page_modules.param_store import (
        discard_stale_cache,
        load_parameters_from_cache,
    )

    cache_at.write_text(
        json.dumps(
            {
                "parameters": _legacy_germline(),
                "timestamp": "2026-01-02T03:04:05",
                "app_version": "2.0.0",
            }
        )
    )

    discarded = discard_stale_cache()

    assert discarded is not None
    assert not cache_at.exists(), "the stale cache is still there to be read next time"
    assert load_parameters_from_cache() is None


def test_the_discard_banner_names_what_was_dropped(cache_at):
    """"Your settings are gone" is not enough to act on; the banner says which settings."""
    from page_modules.param_store import discard_stale_cache

    cache_at.write_text(
        json.dumps(
            {
                "parameters": _legacy_germline(),
                "timestamp": "2026-01-02T03:04:05",
                "app_version": "2.0.0",
            }
        )
    )

    summary = discard_stale_cache().summary()

    assert "germline" in summary, summary
    assert "2026-01-02" in summary, summary
    for named in ("min_depth", "filter_intervar"):
        assert named in summary, f"{named} was dropped without being named: {summary}"


def test_the_discarded_cache_is_kept_aside_so_the_user_can_migrate_it(cache_at):
    """Discarded by the app, not destroyed behind the user's back.

    The app refuses to migrate it — the stamp cannot tell an old cache from a current one,
    which is the whole reason this ticket discards rather than migrates. Uploading it is a
    different matter: that is the user saying "this file is old", which is exactly the
    information the stamp is missing.
    """
    from page_modules.param_store import discard_stale_cache

    cache_at.write_text(
        json.dumps({"parameters": _legacy_somatic(), "timestamp": "2026-01-02T03:04:05"})
    )

    discarded = discard_stale_cache()

    kept = Path(discarded.kept_at)
    assert kept.exists()
    assert json.loads(kept.read_text())["parameters"]["sample_type"] == "somatic"
    assert discarded.kept_at in discarded.summary()


def test_a_second_discard_does_not_overwrite_the_first_kept_copy(cache_at):
    """Moving rather than deleting is a promise that the file is still there."""
    from page_modules.param_store import discard_stale_cache

    cache_at.write_text(json.dumps({"parameters": {"sample_type": "somatic", "min_depth": 11}}))
    first = discard_stale_cache()

    cache_at.write_text(json.dumps({"parameters": {"sample_type": "germline", "min_depth": 22}}))
    second = discard_stale_cache()

    assert second.kept_at != first.kept_at
    assert json.loads(Path(first.kept_at).read_text())["parameters"]["min_depth"] == 11
    assert json.loads(Path(second.kept_at).read_text())["parameters"]["min_depth"] == 22


def test_a_cache_from_a_newer_format_is_left_where_it_is(cache_at):
    """Not read, and not moved either.

    This version cannot know what a newer format's values mean, so it declines to restore
    one — but moving another version's file is not its business, and the user is very
    likely to run that version again.
    """
    from page_modules.param_store import (
        discard_stale_cache,
        load_parameters_from_cache,
    )

    document = {
        SCHEMA_VERSION_KEY: PARAM_SCHEMA_VERSION + 1,
        "parameters": {"sample_type": "somatic", "min_depth": 7},
    }
    cache_at.write_text(json.dumps(document))

    assert discard_stale_cache() is None
    assert cache_at.exists()
    assert load_parameters_from_cache() is None


def test_a_stamped_cache_is_restored_untouched(cache_at):
    """This ticket discards the *pre-stamp* caches. The feature itself survives."""
    from page_modules.param_store import (
        discard_stale_cache,
        load_parameters_from_cache,
        save_parameters_to_cache,
    )

    params = pipeline_params("somatic")
    params["min_depth"] = 33
    assert save_parameters_to_cache(params)

    assert discard_stale_cache() is None
    assert cache_at.exists()
    assert load_parameters_from_cache() == params


def test_a_saved_cache_carries_the_stamp(cache_at):
    """Otherwise the next version's discard would throw away this version's caches too."""
    from page_modules.param_store import save_parameters_to_cache

    save_parameters_to_cache(pipeline_params("somatic"))

    document = json.loads(cache_at.read_text())
    assert document[SCHEMA_VERSION_KEY] == PARAM_SCHEMA_VERSION
    assert document["parameters"] == PIPELINE_PARAMS["somatic"]


def test_no_cache_is_not_a_discard(cache_at):
    """A first run must not announce that it dropped something."""
    from page_modules.param_store import discard_stale_cache

    assert discard_stale_cache() is None


def test_an_unreadable_cache_is_discarded_rather_than_raised_over(cache_at):
    """A truncated write must not turn the parameter page into a stack trace."""
    from page_modules.param_store import (
        discard_stale_cache,
        load_parameters_from_cache,
    )

    cache_at.write_text("{not json at all")

    assert load_parameters_from_cache() is None
    assert discard_stale_cache() is not None
    assert not cache_at.exists()


# ---------------------------------------------------------------------------
# Exports gain the envelope
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("fmt", ["json", "yaml"])
def test_an_export_carries_the_envelope_and_re_imports_as_a_no_op(fmt):
    """Both download buttons, through the same seam the page uses.

    The round trip is the assertion that matters: an envelope nothing can read back is
    just decoration.
    """
    import yaml

    from page_modules.param_store import export_text

    params = pipeline_params("germline")
    params["min_depth"] = 42

    text = export_text(params, fmt)
    document = json.loads(text) if fmt == "json" else yaml.safe_load(text)

    assert document[SCHEMA_VERSION_KEY] == PARAM_SCHEMA_VERSION
    assert document["parameters"] == params

    result = migrate_params(document)
    assert result.migrated is False
    assert result.params == params


@pytest.mark.parametrize("fmt", ["json", "yaml"])
def test_the_dash_clinvar_term_survives_an_export_and_re_import(fmt):
    """``-`` is a real ClinVar term, and the first one whose *spelling* can break a file.

    Issue #103 put it on offer and into both Broad presets: the institute's term table calls
    it *not submitted as a classification; the variant was submitted only as part of a
    haplotype*, which is a statement about the record and not an empty cell — #98 gave it a
    class on that basis.

    The test above round-trips ``pipeline_params("germline")``, whose ClinVar list is
    ``Pathogenic``/``Likely_pathogenic``, so it never carries this term and cannot see the
    risk. A bare ``-`` at the start of a YAML scalar is block-sequence syntax; a serialiser
    that emitted it unquoted would produce a parameter file that reads back as a null, or as a
    nested list, and the user's ClinVar selection would come back short by a term with nothing
    saying so. PyYAML quotes it, which is why this passes rather than why it is unnecessary —
    the point is that the download a user keeps is checked against the vocabulary that fills
    it, so the next term with awkward spelling fails here rather than in someone's export.
    """
    import yaml

    from config.presets import SOFT_SOMATIC_PARAMS
    from page_modules.param_store import export_text

    import copy

    params = copy.deepcopy(SOFT_SOMATIC_PARAMS)
    assert "-" in params["filter_clinvar"], (
        "the Broad somatic preset no longer keeps the `-` ClinVar term, so this test is "
        "guarding nothing -- check config/presets.py rather than relaxing this"
    )

    text = export_text(params, fmt)
    document = json.loads(text) if fmt == "json" else yaml.safe_load(text)

    assert document["parameters"]["filter_clinvar"] == params["filter_clinvar"], (
        f"the {fmt} export changed the ClinVar keep-list: "
        f"{document['parameters']['filter_clinvar']} came back for "
        f"{params['filter_clinvar']}"
    )

    result = migrate_params(document)
    assert result.params["filter_clinvar"] == params["filter_clinvar"], (
        "re-importing the export changed the ClinVar keep-list"
    )


# ---------------------------------------------------------------------------
# The upload path, and the report surviving the rerun
# ---------------------------------------------------------------------------


def test_an_upload_migrates_into_the_session_and_parks_its_report():
    """The dict lands, and what the migration cost is kept for the next render.

    Written this way because of what ``st.rerun`` does: it raises immediately and the
    run's output is discarded, so a message drawn just before it is thrown away. Emitting
    the notes inline would have satisfied a reading of the code and shown the user
    nothing — and "unknown keys dropped **and named**" is an acceptance criterion.
    """
    from page_modules import parameter_config

    session = FakeSessionState()
    with patch.object(parameter_config, "st", MagicMock()) as st_mock:
        # The real `st.rerun` raises to end the run; the mock's simply returns, which is
        # what lets this call be inspected at all.
        st_mock.session_state = session
        parameter_config._adopt_uploaded_parameters(_legacy_somatic(), "old_params.json")
        assert st_mock.rerun.called, "the widgets are not rebuilt around the new dict"

    assert session["filter_params"]["sample_type"] == "somatic"
    assert "keep_pathogenic" not in session["filter_params"]

    report = session[parameter_config.MIGRATION_REPORT_KEY]
    assert report["filename"] == "old_params.json"
    # That the file *was* migrated is now said by the adoption notice rather than by the
    # report, because one sentence has to carry the arm clause every replacement route
    # shares (issue #133). The report keeps what only an upload has.
    adoption = session[parameter_config.PARAM_NOTICE_KEY]
    assert "migrated onto" in adoption["confirmation"], adoption["confirmation"]
    assert "old_params.json" in adoption["confirmation"], adoption["confirmation"]
    assert any("notes_for_myself" in note for note in report["notes"]), report["notes"]
    assert "somatic_gene_set" in report["dropped"]


def test_every_dropped_key_reaches_the_screen():
    """Not only the ones a rule had something to say about.

    Some keys are dropped by a rule that needs no sentence — a saved panel name, the other
    arm's gene list where the file already carries ``filter_genes`` — and those are exactly
    the ones that would otherwise be dropped without appearing anywhere. The criterion is
    "unknown keys dropped **and named**", and the report is the only place naming happens.
    """
    from page_modules import parameter_config

    session = FakeSessionState()
    with patch.object(parameter_config, "st", MagicMock()) as st_mock:
        st_mock.session_state = session
        parameter_config._adopt_uploaded_parameters(_legacy_somatic(), "old_params.json")

    dropped = session[parameter_config.MIGRATION_REPORT_KEY]["dropped"]
    assert dropped

    with patch.object(parameter_config, "st", MagicMock()) as st_mock:
        st_mock.session_state = session
        parameter_config.show_migration_report()
        drawn = " ".join(str(call[0][0]) for call in st_mock.info.call_args_list)

    for key in dropped:
        assert key in drawn, f"{key} was dropped without being named on screen"


def test_the_parked_report_is_drawn_once():
    """Drawn on the run after the rerun, and not on every run after that."""
    from page_modules import parameter_config

    session = FakeSessionState({
        parameter_config.MIGRATION_REPORT_KEY: {
            "filename": "old_params.json",
            "notes": ["Dropped 1 key this app has no parameter for: notes_for_myself."],
        }
    })
    with patch.object(parameter_config, "st", MagicMock()) as st_mock:
        st_mock.session_state = session

        parameter_config.show_migration_report()
        assert st_mock.info.call_count == 1
        assert "notes_for_myself" in st_mock.info.call_args[0][0]

        parameter_config.show_migration_report()
        assert st_mock.info.call_count == 1, "the report was drawn a second time"


# ---------------------------------------------------------------------------
# The shipped examples, which are the file format's documentation
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "name,arm", [("example_parameters.json", "somatic"), ("example_parameters.yaml", "germline")]
)
def test_the_shipped_example_files_load_as_written(name, arm):
    """They document the format, so they have to *be* the format.

    Unpinned they are a third copy of the contract, free to drift into showing users a
    shape the app no longer accepts — which is what they were before this ticket: sentinel
    values, ``min_freq``, and catalogues of options nested under keys nothing reads. A
    round trip is the strongest cheap check: it fails if the stamp is wrong, if a key has
    been renamed under them, or if the contract has moved.
    """
    import yaml

    path = STREAMLIT_APP / name
    document = (json.loads if name.endswith(".json") else yaml.safe_load)(path.read_text())

    assert document[SCHEMA_VERSION_KEY] == PARAM_SCHEMA_VERSION

    result = migrate_params(document)
    assert result.migrated is False, f"{name} is not written in the current format"
    assert result.dropped == ()
    assert result.params["sample_type"] == arm

    contract = pipeline_params(arm)
    if arm == "germline":
        contract["vaf_threshold_germline"] = contract["vaf_threshold"]
    assert result.params == contract, (
        f"{name} has drifted from the contract it is supposed to illustrate"
    )


# ---------------------------------------------------------------------------
# The dead code this ticket removes
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "name",
    [
        "map_clinvar_category",
        "CLINVAR_CATEGORY_PRIORITY",
        "CLINVAR_MAIN_CATEGORIES",
        "load_preset_from_file",
        "process_preset_params",
    ],
)
def test_the_superseded_migration_code_is_gone(name):
    """Five pieces, each dead for its own reason.

    The ClinVar mapper rewrote any unrecognised label to ``"other"`` — a *real* CLNSIG
    term, so a value that matched nothing quietly became one that matches something. The
    preset loader was unreachable: every preset declares ``"file": None``. And
    ``process_preset_params`` is what this ticket replaces.
    """
    from page_modules import parameter_config

    assert not hasattr(parameter_config, name)


def test_no_deprecated_parameter_list_remains():
    """It had one member, ``skip_civic``, which is a live parameter on both arms."""
    source = (STREAMLIT_APP / "page_modules" / "parameter_config.py").read_text()
    assert "deprecated_params" not in source
