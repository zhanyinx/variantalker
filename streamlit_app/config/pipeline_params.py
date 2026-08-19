"""The pipeline's own filter defaults, as a named app parameter set.

``PIPELINE_PARAMS`` is what "at parity with the pipeline" *means* in this app,
written down. Its values are ``nextflow.config``'s ``params`` block; its keys are the
long options ``bin/filter_variants.py`` accepts, with the pipeline's ``skip_*``
polarity rather than the app's inverted ``keep_*``. From here on, "the contract" in
this codebase refers to this module.

It is **not** the app's default yet. The app still applies its own re-implementation
of the pipeline's filter, so loading these values would change which rows the app
reports without making them agree with the pipeline — the two changes have to be
separate. ``tests/test_param_contract.py`` asserts the deferral, and the ticket that
routes the filter through the vendored pipeline code is the one that deletes that
assertion.

Why the config block and not argparse
-------------------------------------
``modules/local/annotation/small_variants/main.nf`` passes *every* argument on every
run, so ``filter_variants.py``'s argparse defaults are unreachable in production and
several of them disagree with the config. The config block is the only set of values
a real run ever uses.

Why literals and not parsing
----------------------------
The packaged .dmg/.exe ships neither ``nextflow.config`` nor ``modules/`` nor
``bin/``, so nothing here can be read at runtime — not the values, and not the
config-parameter-to-CLI-option mapping. Both are literals below, and
``tests/test_param_contract.py`` is what keeps them honest: it reads all three
pipeline files and fails when a default moves, when a name changes, or when a new
filter parameter is wired up that the contract does not carry.

Four decisions worth reading before editing (issue #31)
-------------------------------------------------------
* **One VAF key.** The pipeline takes two thresholds and uses one per arm, so the
  contract carries a single ``vaf_threshold`` whose value depends on the arm. The
  app's second key is what let three sources drift to three different values.
* **One gene-list key.** ``filter_genes``, unprefixed. The pipeline needs
  ``_somatic``/``_germline`` names because one command line configures either arm;
  the app knows its arm. Gene *set* selection (the panel dropdown) is UI state that
  resolves to this list, never a parameter of its own.
* **Atomic ClinVar terms.** The pipeline splits ``ClinVar_VCF_CLNSIG`` on
  ``[|/;,]`` and matches the pieces, so composite keep-values such as
  ``Pathogenic/Likely_pathogenic`` cannot match anything. This decision recorded that
  the app's presets carried one anyway, and that stayed true for four tickets. Issue
  #88 acted on it: the composites are out of ``config.vocabularies.CLINVAR_OPTIONS``
  and out of both Broad presets, measured at zero row change, and the guards that keep
  them out run the matcher itself rather than restating the ``[|/;,]`` rule.
* **Per-source spacing.** InterVar and ReNovo values carry spaces; CancerVar and
  ClinVar values carry underscores. That is how the annotations are written, so
  there is no normalisation to apply and applying one would break two sources.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass

#: The two arms, as a tuple. Deliberately not imported from ``constants.SAMPLE_TYPES``
#: — this module is the new canonical statement of the filter's parameters and does not
#: depend on the legacy constants grab-bag — but pinned to it by
#: ``tests/test_param_contract.py`` so the two spellings cannot diverge.
ARMS = ("somatic", "germline")


# ---------------------------------------------------------------------------
# Provenance: where each contract value comes from
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ParamOrigin:
    """One contract key, and the ``nextflow.config`` parameter behind it.

    ``config`` is ``None`` for the app's own additions, which have no pipeline
    counterpart to be checked against, and ``cli`` is then empty too. Where the
    pipeline splits by arm what the app keeps as one key — the VAF thresholds and the
    gene lists — there are two origins for the one ``key``, distinguished by ``arms``.
    """

    key: str
    config: str | None
    cli: str = ""
    #: "list" (comma-separated), "scalar" (number or bool), or "genes" (a file path
    #: in the pipeline, a list of symbols in the app).
    kind: str = "scalar"
    arms: tuple = ARMS

    def split(self, literal) -> list:
        """A config literal as the list of terms it denotes."""
        if literal is None:
            return []
        return [term for term in str(literal).split(",")]

    def expected_from(self, literal):
        """The contract value this origin's config literal should produce.

        Declared exceptions are applied here, so that a departure is only tolerated
        where ``DECLARED_EXCEPTIONS`` names it.
        """
        if self.kind == "genes":
            # Both filter_genes_* default to null: no gene filtering.
            return [] if literal is None else literal
        if self.kind == "list":
            terms = self.split(literal)
            for exception in exceptions_for(self.key):
                terms = [
                    exception.contract_value if t == exception.config_value else t for t in terms
                ]
            return terms
        return literal


@dataclass(frozen=True)
class DeclaredException:
    """A value the contract deliberately spells differently from the config."""

    key: str
    config_value: str
    contract_value: str
    reason: str


#: The contract's *only* departure from the ``nextflow.config`` literal.
#:
#: The config asks ClinVar to keep ``Likely pathogenic``, spaced. ClinVar's own
#: ``CLNSIG`` values are underscored, so after the pipeline splits the field on
#: ``[|/;,]`` the spaced form matches 0 of the 181,566 rows in the GERSOM reference —
#: it is a keep-term that keeps nothing. The contract carries the underscored form so
#: the clause does what the config plainly intends.
#:
#: Declared rather than hidden: ``contract_findings`` consults this list instead of
#: tolerating any ClinVar mismatch, so a *second* silent departure fails the guard,
#: and a stale exception — one the config no longer departs from — fails it too.
DECLARED_EXCEPTIONS = (
    DeclaredException(
        key="filter_clinvar",
        config_value="Likely pathogenic",
        contract_value="Likely_pathogenic",
        reason=(
            "ClinVar CLNSIG values are underscored; the config's spaced form matches "
            "0 of 181,566 reference rows (issue #28)"
        ),
    ),
)


def exceptions_for(key: str) -> tuple:
    """Declared exceptions applying to one contract key.

    The single place this relation is resolved. ``ParamOrigin.expected_from`` needs it
    to build the value it expects; the guard needs it to check that each declared
    exception still departs from something. Both go through here rather than
    re-deriving the join in opposite directions.
    """
    return tuple(e for e in DECLARED_EXCEPTIONS if e.key == key)


def origins_for(key: str) -> tuple:
    """Every declared origin of one contract key (two, where the pipeline splits arms)."""
    return tuple(o for o in ORIGINS if o.key == key)


#: Every contract key, with the config parameter and CLI option behind it. The guard
#: walks this to compare values, and to check each name against ``main.nf`` and
#: against ``filter_variants.py``'s argparse.
ORIGINS = (
    ParamOrigin("min_depth", "filter_min_depth", "min_depth"),
    # Two config parameters, one contract key: the value differs per arm, the meaning
    # does not. `arms` records which arm each mapping applies to.
    ParamOrigin("vaf_threshold", "filter_vaf_threshold", "vaf_threshold", arms=("somatic",)),
    ParamOrigin(
        "vaf_threshold",
        "filter_vaf_threshold_germline",
        "vaf_threshold_germline",
        arms=("germline",),
    ),
    ParamOrigin(
        "filter_variant_classification",
        "filter_var_classification",
        "filter_variant_classification",
        kind="list",
    ),
    ParamOrigin(
        "filter_cancervar", "filter_cancervar", "filter_cancervar", kind="list", arms=("somatic",)
    ),
    ParamOrigin(
        "filter_civic",
        "filter_civic_evidence_level",
        "filter_civic",
        kind="list",
        arms=("somatic",),
    ),
    ParamOrigin("filter_escat", "filter_escat", "filter_escat", kind="list", arms=("somatic",)),
    ParamOrigin(
        "filter_intervar", "filter_intervar", "filter_intervar", kind="list", arms=("germline",)
    ),
    ParamOrigin("filter_renovo", "filter_renovo", "filter_renovo", kind="list", arms=("germline",)),
    ParamOrigin("filter_clinvar", "filter_clinvar", "filter_clinvar", kind="list"),
    ParamOrigin(
        "filter_genes",
        "filter_genes_somatic",
        "filter_genes_somatic",
        kind="genes",
        arms=("somatic",),
    ),
    ParamOrigin(
        "filter_genes",
        "filter_genes_germline",
        "filter_genes_germline",
        kind="genes",
        arms=("germline",),
    ),
    ParamOrigin("skip_civic", "skip_civic", "skip_civic"),
    ParamOrigin("skip_pathogenic", "skip_pathogenic_retention", "skip_pathogenic"),
    # App-only. The pipeline has no population-frequency filter, so 1.0 — keep
    # everything — is the only value at which this layer is neutral.
    ParamOrigin("max_freq_population", None),
)


# ---------------------------------------------------------------------------
# The contract
# ---------------------------------------------------------------------------

#: Keys shared by both arms. Values are nextflow.config lines 88-96 (plus the one
#: declared exception above); the guard is what keeps that true.
_COMMON = {
    "min_depth": 50,                                            # filter_min_depth
    # filter_var_classification. An *exclude* list, which is the pipeline's meaning:
    # these three classifications are dropped. The app's parameter of the same name is
    # an include list today — the divergence that makes stating this worthwhile.
    "filter_variant_classification": ["Silent", "IGR", "RNA"],
    # filter_clinvar, atomic and underscored (see DECLARED_EXCEPTIONS).
    "filter_clinvar": ["Pathogenic", "Likely_pathogenic"],
    # filter_genes_{somatic,germline} are both null: no gene filtering.
    "filter_genes": [],
    # skip_civic. On **both** arms, even though germline_filters() takes no skip_civic
    # argument: the pipeline passes --skip_civic for either arm because it also decides
    # the output column set, and the app's column resolver takes it for the same reason.
    # Somatic is the only arm where it additionally changes which rows pass.
    "skip_civic": False,
    "skip_pathogenic": False,                                   # skip_pathogenic_retention
    # App-only, neutral. Not a pipeline parameter at any value.
    "max_freq_population": 1.0,
}

#: The somatic contract. CancerVar, CIViC and ESCAT are the somatic guideline
#: sources — the three arguments ``somatic_filters()`` takes and
#: ``germline_filters()`` does not.
PIPELINE_SOMATIC_PARAMS = {
    # Arm identity, not a filter setting: it selects which filter runs rather than
    # tuning one. Hence no ParamOrigin — the pipeline takes it per run, from the sample
    # sheet, not from the params block — and hence its place in the guard's
    # NON_CONTRACT_CLI_OPTIONS.
    "sample_type": "somatic",
    **_COMMON,
    "vaf_threshold": 0.01,                                      # filter_vaf_threshold
    "filter_cancervar": ["Tier_II_potential", "Tier_I_strong"],
    "filter_civic": ["A", "B", "C"],                            # filter_civic_evidence_level
    "filter_escat": ["IA", "IB", "IC", "IIA", "IIB"],
}

#: The germline contract. InterVar and ReNovo, spaced as the annotations write them,
#: and **no ESCAT**: ``germline_filters()`` takes no ESCAT argument, and the app's
#: extra germline ESCAT clause is the single largest divergence on real data.
PIPELINE_GERMLINE_PARAMS = {
    "sample_type": "germline",
    **_COMMON,
    "vaf_threshold": 0.2,                                       # filter_vaf_threshold_germline
    "filter_intervar": ["Pathogenic", "Likely pathogenic"],
    "filter_renovo": ["LP Pathogenic", "IP Pathogenic", "HP Pathogenic"],
}

#: The contract, by sample type. Use :func:`pipeline_params` to get a copy.
PIPELINE_PARAMS = {
    "somatic": PIPELINE_SOMATIC_PARAMS,
    "germline": PIPELINE_GERMLINE_PARAMS,
}


def pipeline_params(sample_type: str) -> dict:
    """A deep copy of the contract for ``sample_type``.

    A copy because the app's parameter pages mutate whatever dict they are handed —
    session state is edited in place — and the keep-lists are nested, so a shallow
    copy would leave a widget able to append to the contract's own list and silently
    redefine what parity means for the rest of the session.
    """
    if sample_type not in PIPELINE_PARAMS:
        raise ValueError(
            f"unknown sample_type {sample_type!r}; expected one of {list(PIPELINE_PARAMS)}"
        )
    return copy.deepcopy(PIPELINE_PARAMS[sample_type])
