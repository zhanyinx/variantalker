"""Drift guard for the parameter contract.

``config/pipeline_params.py`` claims to hold the pipeline's *own* filter defaults —
the values in ``nextflow.config``'s ``params`` block — under the pipeline's CLI
long-option names. That claim is only worth anything while it is provably still
true, so this module is the proof.

Three things are pinned, not one:

1. **The values.** Each contract value equals its ``nextflow.config`` literal, with
   exactly one declared departure (``Likely_pathogenic``, see
   ``config.pipeline_params.DECLARED_EXCEPTIONS``). An undeclared departure fails,
   and so does a declared exception that no longer departs from anything.
2. **The names.** Each config parameter really is passed to
   ``bin/filter_variants.py`` under the long option the contract keys itself by.
   The mapping is read out of ``modules/local/annotation/small_variants/main.nf``
   and out of the script's own argparse — the app never infers it at runtime,
   because the packaged .dmg/.exe carries neither file.
3. **The shape.** The contract-specific decisions issue #31 made — one VAF key
   rather than two, atomic ClinVar terms, per-source spacing, a neutral frequency
   layer, ESCAT somatic-only.

``test_a_changed_config_default_fails_the_guard`` is what makes this falsifiable:
it mutates a copy of the config text and asserts the checker reports it.

Both the config and the module skip when absent: this guard needs the pipeline
checkout, and the packaged app deliberately ships neither ``nextflow.config`` nor
this test. Nothing here imports pandas, streamlit, or the filter code — the
comparison is on config text and plain dicts.
"""

import ast
import re
from pathlib import Path

import pytest

from config.pipeline_params import ARMS

STREAMLIT_APP = Path(__file__).resolve().parents[1]
REPO_ROOT = STREAMLIT_APP.parent
NEXTFLOW_CONFIG = REPO_ROOT / "nextflow.config"
FILTER_MODULE = REPO_ROOT / "modules" / "local" / "annotation" / "small_variants" / "main.nf"
FILTER_SCRIPT = REPO_ROOT / "bin" / "filter_variants.py"
ESCAT_TIERING = REPO_ROOT / "resources" / "escat_tiering.csv"

needs_pipeline_config = pytest.mark.skipif(
    not NEXTFLOW_CONFIG.is_file(),
    reason="nextflow.config not present (packaged app or partial checkout)",
)
needs_filter_module = pytest.mark.skipif(
    not FILTER_MODULE.is_file(),
    reason="pipeline modules/ not present (packaged app or partial checkout)",
)
needs_filter_script = pytest.mark.skipif(
    not FILTER_SCRIPT.is_file(),
    reason="pipeline bin/ not present (packaged app or partial checkout)",
)
needs_escat_tiering = pytest.mark.skipif(
    not ESCAT_TIERING.is_file(),
    reason="pipeline resources/ not present (packaged app or partial checkout)",
)

FIX_HINT = (
    "If the pipeline default genuinely changed, update config/pipeline_params.py to "
    "match it — and expect the parity baselines to move with it."
)


# ---------------------------------------------------------------------------
# Reading nextflow.config
# ---------------------------------------------------------------------------


def _strip_comment(value: str) -> str:
    """Drop a trailing ``//`` comment, ignoring ``//`` inside a quoted string.

    Needed because several params carry URLs and DSN examples in their comments
    (``postgresql://…``), so cutting at the first ``//`` unconditionally would be
    wrong for any value that quotes one.
    """
    quote = None
    for i, char in enumerate(value):
        if quote:
            if char == quote:
                quote = None
        elif char in "\"'":
            quote = char
        elif char == "/" and value[i : i + 2] == "//":
            return value[:i]
    return value


def _coerce(literal: str):
    """A Groovy scalar literal as the Python value it denotes."""
    literal = literal.strip()
    if len(literal) >= 2 and literal[0] == literal[-1] and literal[0] in "\"'":
        return literal[1:-1]
    if literal == "true":
        return True
    if literal == "false":
        return False
    if literal == "null":
        return None
    try:
        return int(literal)
    except ValueError:
        pass
    try:
        return float(literal)
    except ValueError:
        return literal


def parse_params_block(text: str) -> dict:
    """``nextflow.config``'s ``params { … }`` block as a dict of Python values.

    A deliberately small reader: it takes ``name = <scalar literal>`` assignments at
    any depth inside the block and ignores everything else. That is enough for every
    parameter the contract names, and it fails loudly rather than quietly if the
    block cannot be found at all.
    """
    lines = text.splitlines()
    start = next(
        (i for i, line in enumerate(lines) if re.match(r"^\s*params\s*\{", line)),
        None,
    )
    if start is None:
        raise AssertionError("no `params {` block found in nextflow.config")

    params = {}
    depth = 0
    for line in lines[start:]:
        depth += line.count("{") - line.count("}")
        match = re.match(r"^\s*([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.+)$", line)
        if match:
            value = _strip_comment(match.group(2)).strip()
            if value:
                params[match.group(1)] = _coerce(value)
        if depth <= 0:
            break
    assert params, "parsed the params block but found no assignments in it"
    return params


@pytest.fixture(scope="module")
def config_params() -> dict:
    """``nextflow.config``'s ``params`` block, parsed.

    The absence gate lives **here**, not on each test that asks for the fixture. Three of
    the five did not carry ``@needs_pipeline_config`` and so raised ``FileNotFoundError``
    on any checkout without the pipeline — the packaged app's shape — turning a legitimate
    skip into three errors. Measured, not guessed: issue #41's ``make test-cov`` runs the
    suite in a tree with no pipeline in it, which is what found them.

    A gate a call site can forget is a gate that will be forgotten. Requesting the fixture
    is the complete statement of the dependency now, and the marker stays only on the two
    tests that read the config file without going through here.
    """
    if not NEXTFLOW_CONFIG.is_file():
        pytest.skip("nextflow.config not present (packaged app or partial checkout)")
    return parse_params_block(NEXTFLOW_CONFIG.read_text())


@pytest.fixture(scope="module")
def origins():
    from config.pipeline_params import ORIGINS

    return ORIGINS


# ---------------------------------------------------------------------------
# The value check, as a reusable function so it can be run against mutated text
# ---------------------------------------------------------------------------


def contract_findings(config_text: str) -> list:
    """Human-readable disagreements between the contract and ``config_text``.

    Empty means the contract still states the pipeline's defaults. Factored out of
    the tests so ``test_a_changed_config_default_fails_the_guard`` can run the very
    same rule against a deliberately edited config.
    """
    from config.pipeline_params import (
        DECLARED_EXCEPTIONS,
        ORIGINS,
        PIPELINE_PARAMS,
        origins_for,
    )

    params = parse_params_block(config_text)
    findings = []

    for origin in ORIGINS:
        if origin.config is None:  # app-only, no counterpart to compare against
            continue
        if origin.config not in params:
            findings.append(
                f"GONE nextflow.config no longer defines `{origin.config}` "
                f"(contract key `{origin.key}`)"
            )
            continue

        want = origin.expected_from(params[origin.config])
        for arm in origin.arms:
            got = PIPELINE_PARAMS[arm].get(origin.key, "<<missing>>")
            if got != want:
                findings.append(
                    f"DRIFT {arm} `{origin.key}` is {got!r}, but nextflow.config's "
                    f"`{origin.config}` says {want!r}"
                )

    for exception in DECLARED_EXCEPTIONS:
        for origin in origins_for(exception.key):
            literal = params.get(origin.config)
            if literal is None:
                continue
            if exception.config_value not in origin.split(literal):
                findings.append(
                    f"STALE declared exception for `{exception.key}`: nextflow.config no "
                    f"longer carries {exception.config_value!r}, so the departure to "
                    f"{exception.contract_value!r} is no longer an exception to anything"
                )

    return findings


# ---------------------------------------------------------------------------
# 1. The values
# ---------------------------------------------------------------------------


@needs_pipeline_config
def test_contract_values_match_the_config_literals():
    """Every contract value is its ``nextflow.config`` literal, exceptions aside."""
    findings = contract_findings(NEXTFLOW_CONFIG.read_text())
    assert not findings, (
        "the parameter contract no longer states the pipeline's defaults:\n  "
        + "\n  ".join(findings)
        + f"\n{FIX_HINT}"
    )


@needs_pipeline_config
def test_a_changed_config_default_fails_the_guard():
    """The guard must actually catch a moved pipeline default.

    Without this, every check above could be passing vacuously — over a params block
    it failed to parse, or a mapping that names nothing. Each mutation below is a
    plausible real edit, one per value shape the contract carries: a number, a
    boolean, and a comma-separated list.
    """
    original = NEXTFLOW_CONFIG.read_text()
    assert not contract_findings(original), "fixture precondition: the contract is green"

    mutations = [
        ("filter_min_depth            = 50", "filter_min_depth            = 99"),
        ("skip_pathogenic_retention   = false", "skip_pathogenic_retention   = true"),
        (
            'filter_escat                = "IA,IB,IC,IIA,IIB"',
            'filter_escat                = "IA,IB"',
        ),
        (
            'filter_intervar             = "Pathogenic,Likely pathogenic"',
            'filter_intervar             = "Pathogenic"',
        ),
        ("filter_vaf_threshold_germline = 0.2", "filter_vaf_threshold_germline = 0.4"),
        # filter_clinvar specifically, because it is the only key whose comparison runs
        # the exception-substitution branch in ParamOrigin.expected_from. Without this,
        # that branch could swallow any ClinVar change and the guard would look green.
        # Both directions: a term added, and a term changed to something else.
        (
            'filter_clinvar              = "Pathogenic,Likely pathogenic"',
            'filter_clinvar              = "Pathogenic,Likely pathogenic,drug_response"',
        ),
        (
            'filter_clinvar              = "Pathogenic,Likely pathogenic"',
            'filter_clinvar              = "Benign,Likely pathogenic"',
        ),
    ]

    for old, new in mutations:
        assert old in original, (
            f"cannot mutate {old!r}: nextflow.config's formatting changed, so this "
            "falsifiability check is no longer testing anything"
        )
        findings = contract_findings(original.replace(old, new))
        assert findings, (
            f"changing nextflow.config from {old.strip()!r} to {new.strip()!r} left the "
            "guard green — a pipeline default can move without the contract noticing"
        )


def test_the_clinvar_underscore_exception_is_the_only_one(config_params):
    """One declared departure from the config literal, named as such.

    ``nextflow.config`` says ``Likely pathogenic``; the contract says
    ``Likely_pathogenic``, because the spaced form matches 0 of 181,566 reference
    rows (issue #28). The point of this test is that the exception is *declared* —
    ``contract_findings`` consults ``DECLARED_EXCEPTIONS`` rather than tolerating
    any mismatch on ClinVar terms, so a second silent departure would fail.
    """
    from config.pipeline_params import DECLARED_EXCEPTIONS, PIPELINE_PARAMS

    assert len(DECLARED_EXCEPTIONS) == 1, (
        "the contract is meant to carry exactly one deliberate departure from the "
        f"config literal; it now declares {len(DECLARED_EXCEPTIONS)}: "
        f"{[e.key for e in DECLARED_EXCEPTIONS]}"
    )
    exception = DECLARED_EXCEPTIONS[0]
    assert exception.key == "filter_clinvar"
    assert exception.config_value == "Likely pathogenic"
    assert exception.contract_value == "Likely_pathogenic"
    assert exception.reason, "a declared exception must say why it exists"

    # And it is live on both arms, rather than declared and unused.
    assert "Likely pathogenic" in config_params["filter_clinvar"]
    for arm in ARMS:
        terms = PIPELINE_PARAMS[arm]["filter_clinvar"]
        assert "Likely_pathogenic" in terms, arm
        assert "Likely pathogenic" not in terms, arm


# ---------------------------------------------------------------------------
# 2. The names — config parameter to CLI long option
# ---------------------------------------------------------------------------


def _module_cli_pairs(text: str) -> set:
    """``(cli_option, config_param)`` pairs the pipeline module actually passes.

    Two shapes, because ``skip_*`` are argparse flags rather than valued options:

        --min_depth ${params.filter_min_depth}
        if [ "${params.skip_civic}" == "true" ]; then  skip_civic="--skip_civic"
    """
    valued = re.findall(r"--([A-Za-z_][A-Za-z0-9_]*)\s+\"?\$\{params\.([A-Za-z0-9_]+)\}", text)
    flags = re.findall(
        r"\$\{params\.([A-Za-z0-9_]+)\}\"\s*==\s*\"true\"\s*\]\s*;\s*then\s*\n"
        r"\s*[A-Za-z_][A-Za-z0-9_]*=\"--([A-Za-z_][A-Za-z0-9_]*)\"",
        text,
    )
    return set(valued) | {(cli, cfg) for cfg, cli in flags}


def _filter_variants_invocation(text: str) -> str:
    """The ``filter_variants.py`` command, from its first line to the block end."""
    start = text.index("filter_variants.py")
    end = text.index('"""', start)
    return text[start:end]


#: CLI options the invocation passes that are deliberately *not* part of the parameter
#: contract. ``--config``, ``--funcotator``, ``--technology``, ``--annovar_version`` and
#: ``--genome_build`` are argparse-required and never read (issue #18); the rest are
#: per-run plumbing (input path, output path, sample arm, project id, skip-annotation
#: MAF) rather than filter settings. Listing them is what lets the test below insist
#: every *other* option is accounted for, so a newly-wired filter parameter cannot
#: reach the pipeline without the contract learning about it.
NON_CONTRACT_CLI_OPTIONS = {
    "maf",
    "output",
    "sample_type",
    "config",
    "funcotator",
    "technology",
    "annovar_version",
    "genome_build",
    "projectid",
    "maf_skip",
}


@needs_filter_module
def test_config_to_cli_mapping_is_the_one_the_pipeline_uses(origins):
    """Each declared origin is really passed under that long option.

    This is the criterion "the config-to-CLI name mapping is pinned by the guard,
    not inferred at runtime": the app cannot read ``main.nf``, because the packaged
    build does not ship it, so the mapping is a literal in ``pipeline_params.py``
    and this test is what keeps the literal honest.
    """
    pairs = _module_cli_pairs(FILTER_MODULE.read_text())
    assert pairs, "could not parse any `--option ${params.x}` pair out of main.nf"

    for origin in origins:
        if origin.config is None:
            continue
        assert (origin.cli, origin.config) in pairs, (
            f"the contract claims nextflow.config's `{origin.config}` reaches "
            f"filter_variants.py as `--{origin.cli}`, but "
            f"{FILTER_MODULE.relative_to(REPO_ROOT)} does not pass that pair. "
            f"Pairs it does pass: {sorted(pairs)}"
        )


@needs_filter_module
def test_no_filter_parameter_is_wired_without_the_contract(origins):
    """Every option in the invocation is either in the contract or declared plumbing."""
    module_text = FILTER_MODULE.read_text()
    invocation = _filter_variants_invocation(module_text)

    passed = {cli for cli, _ in _module_cli_pairs(invocation)}
    # The skip_* flags are set up above the command and interpolated as shell vars, so
    # they are read from the whole process block rather than the command text alone.
    passed |= {
        cli
        for cli, cfg in _module_cli_pairs(module_text)
        if cfg in ("skip_civic", "skip_pathogenic_retention")
    }
    assert passed, "could not parse the filter_variants.py invocation out of main.nf"

    known = {o.cli for o in origins if o.cli} | NON_CONTRACT_CLI_OPTIONS
    unaccounted = sorted(passed - known)
    assert not unaccounted, (
        f"{FILTER_MODULE.relative_to(REPO_ROOT)} passes filter_variants.py option(s) "
        f"{unaccounted} that the parameter contract neither carries nor lists as "
        "non-contract plumbing. If they affect which rows the pipeline reports, the "
        "app is now off parity by exactly them."
    )


@needs_filter_script
def test_every_contract_option_is_a_real_long_option(origins):
    """The long options the contract keys itself by exist in the script's argparse."""
    script = FILTER_SCRIPT.read_text()
    declared = set(re.findall(r'"--([A-Za-z_][A-Za-z0-9_]*)"', script))
    assert declared, "could not find any long option in bin/filter_variants.py"

    for origin in origins:
        if not origin.cli:  # app-only, no CLI counterpart
            continue
        assert origin.cli in declared, (
            f"the contract keys `{origin.key}` off `--{origin.cli}`, which "
            f"bin/filter_variants.py does not define. Options it defines: "
            f"{sorted(declared)}"
        )


# ---------------------------------------------------------------------------
# 3. The shape — issue #31's contract decisions
# ---------------------------------------------------------------------------


def test_keys_are_cli_long_option_names_with_skip_polarity(origins):
    """Keys read as the pipeline's CLI, including its ``skip_*`` polarity."""
    from config.pipeline_params import PIPELINE_PARAMS

    for arm in ARMS:
        params = PIPELINE_PARAMS[arm]
        assert "skip_pathogenic" in params, arm
        assert params["skip_pathogenic"] is False, arm
        assert "keep_pathogenic" not in params, (
            f"{arm}: `keep_pathogenic` is the app's old inverted spelling; the contract "
            "takes the pipeline's `skip_pathogenic` polarity"
        )
        assert "skip_civic" in params, arm
        assert params["skip_civic"] is False, arm

        # Old app-side spellings that the CLI has no counterpart for.
        for stale in ("somatic_genes", "germline_genes", "somatic_gene_set", "germline_gene_set"):
            assert stale not in params, (
                f"{arm}: `{stale}` is app-side vocabulary. The contract carries one "
                "unprefixed `filter_genes`; gene *set* selection is UI state, not a "
                "filter parameter."
            )


def test_one_gene_list_key_not_a_somatic_germline_pair():
    """`filter_genes`, unprefixed and neutral — the pipeline needs two names, not the app."""
    from config.pipeline_params import PIPELINE_PARAMS

    for arm in ARMS:
        params = PIPELINE_PARAMS[arm]
        assert "filter_genes" in params, arm
        assert params["filter_genes"] == [], (
            f"{arm}: nextflow.config leaves both filter_genes_* at null, so the "
            f"contract's gene list must be empty (no gene filtering), not "
            f"{params['filter_genes']!r}"
        )
        for prefixed in ("filter_genes_somatic", "filter_genes_germline"):
            assert prefixed not in params, f"{arm}: {prefixed}"


def test_one_vaf_threshold_key(config_params):
    """The redundant second VAF threshold is collapsed to one key.

    Three sources said three things before this: ``nextflow.config`` and both
    germline presets 0.2, and the germline widget's own fallback 0.3. The contract
    carries a single ``vaf_threshold``, per arm, so it cannot join that argument.
    """
    from config.pipeline_params import PIPELINE_PARAMS

    for arm in ARMS:
        assert "vaf_threshold" in PIPELINE_PARAMS[arm], arm
        assert "vaf_threshold_germline" not in PIPELINE_PARAMS[arm], (
            f"{arm}: the contract must carry one VAF key. A second one is how the "
            "three-way disagreement got in."
        )

    assert PIPELINE_PARAMS["somatic"]["vaf_threshold"] == config_params["filter_vaf_threshold"]
    assert (
        PIPELINE_PARAMS["germline"]["vaf_threshold"]
        == config_params["filter_vaf_threshold_germline"]
    )


PARAMETER_CONFIG = STREAMLIT_APP / "page_modules" / "parameter_config.py"

#: The one parameter name a control may still resolve with a literal of its own, and why.
#:
#: ``sample_type`` is arm identity rather than a filter setting — it selects which filter
#: runs — and *every* reader in this app spells that reading the same way,
#: ``.get("sample_type", "somatic")``: a dict that states no arm is somatic. That is a
#: documented convention about which filter to run, not a filter value being invented, and
#: ``complete_params`` deliberately leaves the key alone for the same reason.
_MAY_CARRY_ITS_OWN_DEFAULT = frozenset({"sample_type"})


def _literal_parameter_defaults(path):
    """Every ``…get("<contract key>", <literal>)`` in one module, as (key, literal, line).

    Read from the AST rather than by grep: black wraps these calls across lines, and the
    regex this replaced could only see them on one. Any receiver counts, because the page
    reads its parameters both as ``st.session_state.filter_params`` and through a local
    ``params`` alias, and a rule that only watched one spelling would be half a rule.
    """
    from config.pipeline_params import PIPELINE_PARAMS

    contract_keys = set().union(*(set(params) for params in PIPELINE_PARAMS.values()))
    found = []
    for node in ast.walk(ast.parse(path.read_text())):
        if not isinstance(node, ast.Call):
            continue
        if not (isinstance(node.func, ast.Attribute) and node.func.attr == "get"):
            continue
        if len(node.args) != 2:
            continue
        key = node.args[0]
        if not (isinstance(key, ast.Constant) and key.value in contract_keys):
            continue
        if key.value in _MAY_CARRY_ITS_OWN_DEFAULT:
            continue
        found.append((key.value, ast.unparse(node.args[1]), node.lineno))
    return found


def test_no_control_resolves_a_parameter_from_a_literal_of_its_own():
    """The wider rule the germline VAF fallback was one instance of (issue #280).

    That guard pinned one widget's literal to the contract, and its own message said to
    delete it if the widget was ever rewritten to read the contract directly. It has been:
    the parameters page now completes its dict against the arm's contract **once, at the
    boundary**, before any tab draws, so a control has nothing left to fall back to.

    Pinning the literals was never going to hold, and the reason is that all three of them
    had already drifted by the time anyone looked. ``max_freq_population`` invented ``0.01``
    where both arms say ``1.0``; the somatic VAF widget invented ``0.05`` where the contract
    says ``0.01``; the germline one invented ``0.2`` under ``vaf_threshold_germline``, a key
    the contract does not spell, so the contract's own value was never read on that arm at
    all. A guard that checks each literal against the contract has to be *written* for each
    literal, and the two nobody wrote were the two that were wrong.

    So this forbids the shape instead of auditing its values. The contract is the only place
    a parameter's default is stated, ``complete_params`` is the only thing that applies it,
    and a control that wants a value indexes for it — which raises loudly if the completion
    is ever bypassed, rather than quietly minting a plausible number.
    """
    offenders = _literal_parameter_defaults(PARAMETER_CONFIG)
    assert not offenders, (
        "these controls resolve a filter parameter from a literal of their own rather than "
        f"from the contract: {offenders}. Delete the fallback and index the key — "
        "`complete_params` at the top of `show_parameter_config_page` guarantees it is "
        "there. A literal here is a second opinion about a default, and every literal this "
        "guard replaced had drifted away from the contract it stood in for (issue #280)."
    )


def test_max_freq_population_is_neutral():
    """1.0 — the app-only frequency layer must not filter anything at parity."""
    from config.pipeline_params import PIPELINE_PARAMS

    for arm in ARMS:
        assert PIPELINE_PARAMS[arm]["max_freq_population"] == 1.0, (
            f"{arm}: the pipeline has no population-frequency filter at all, so any "
            "value below 1.0 puts the app off parity by whatever it removes"
        )


def test_max_freq_population_has_no_config_counterpart(origins, config_params):
    """It is declared app-only, not silently mapped onto some config parameter."""
    origin = next(o for o in origins if o.key == "max_freq_population")
    assert origin.config is None, (
        "max_freq_population is an app-only layer; the pipeline has no counterpart"
    )
    assert not [name for name in config_params if "freq" in name.lower()], (
        "nextflow.config has grown a frequency parameter — the app's frequency layer "
        "may no longer be app-only, and this contract entry needs revisiting"
    )


def test_clinvar_terms_are_atomic():
    """Post-split terms only; the composite ``/``-joined values are gone.

    The pipeline splits ``ClinVar_VCF_CLNSIG`` on ``[|/;,]`` and matches the pieces,
    so a composite keep-value like ``Pathogenic/Likely_pathogenic`` can never match
    anything.

    This guarded the contract alone, and said so: *the app's presets still carry it; the
    contract must not*. That sentence held for four tickets. Issue #88 made it false by
    taking the composites out of ``CLINVAR_OPTIONS`` and out of both Broad presets, so the
    property is now the app's everywhere — and the guards that hold it there run against
    the matcher itself, in ``test_app_defaults.py``. This one stays as it was: the
    contract is a separate statement of the same values and can drift on its own.
    """
    from config.pipeline_params import PIPELINE_PARAMS

    for arm in ARMS:
        terms = PIPELINE_PARAMS[arm]["filter_clinvar"]
        composite = [t for t in terms if re.search(r"[|/;,]", t)]
        assert not composite, (
            f"{arm}: composite ClinVar keep-terms {composite} cannot match a "
            "post-split value and are dead weight in the contract"
        )
        assert terms, f"{arm}: the ClinVar keep-list must not be empty"


def test_spacing_stays_per_source(config_params):
    """Two sources spaced, two underscored — no global normalisation.

    InterVar and ReNovo really do carry spaces in the data; CancerVar and ClinVar
    really do carry underscores. Normalising either way would break the sources it
    was normalised away from, so the contract keeps each source's own spelling and
    this test pins the mixture rather than a single convention.
    """
    from config.pipeline_params import PIPELINE_PARAMS

    germline = PIPELINE_PARAMS["germline"]
    somatic = PIPELINE_PARAMS["somatic"]

    spaced = {
        "filter_intervar": germline["filter_intervar"],
        "filter_renovo": germline["filter_renovo"],
    }
    underscored = {
        "filter_cancervar": somatic["filter_cancervar"],
        "filter_clinvar": somatic["filter_clinvar"],
    }

    for key, values in spaced.items():
        assert any(" " in v for v in values), (
            f"{key} lost its spaces — it is one of the two sources whose values are "
            f"space-separated in the data: {values}"
        )
        assert not any("_" in v for v in values), f"{key} gained underscores: {values}"
        # These two are the config literal exactly, which is the strongest form of
        # "not normalised".
        assert values == config_params[key].split(",")

    for key, values in underscored.items():
        assert any("_" in v for v in values), (
            f"{key} lost its underscores — it is one of the two sources whose values "
            f"are underscored in the data: {values}"
        )
        assert not any(" " in v for v in values), f"{key} gained spaces: {values}"


def test_arms_agree_with_the_apps_sample_types():
    """``ARMS`` is the tuple form of ``vocabularies.SAMPLE_TYPES``, not a second opinion."""
    from config.vocabularies import SAMPLE_TYPES
    from config.pipeline_params import PIPELINE_PARAMS

    assert ARMS == tuple(SAMPLE_TYPES), (
        "config/pipeline_params.py's ARMS and config/vocabularies.py's SAMPLE_TYPES "
        f"disagree: {ARMS} vs {tuple(SAMPLE_TYPES)}"
    )
    assert set(PIPELINE_PARAMS) == set(ARMS), "the contract must cover every arm"


def test_every_contract_key_has_a_declared_origin():
    """No contract key without provenance.

    ``contract_findings`` walks ``ORIGINS`` and checks each against the config, so a
    key added to ``PIPELINE_PARAMS`` with no ``ParamOrigin`` would be a value nothing
    compares to anything — app-side vocabulary wearing the contract's name. This is the
    other direction of the same rule ``test_vendor_drift.py`` applies to the vendored
    symbols with ``test_no_extra_functions_vendored_silently``.
    """
    from config.pipeline_params import ORIGINS, PIPELINE_PARAMS

    declared = {origin.key for origin in ORIGINS}
    # sample_type is arm identity rather than a filter setting: the pipeline takes it
    # per run from the sample sheet, not from the params block, so it has no origin by
    # design. Named here so that is a decision rather than an omission.
    exempt = {"sample_type"}

    for arm in ARMS:
        unaccounted = sorted(set(PIPELINE_PARAMS[arm]) - declared - exempt)
        assert not unaccounted, (
            f"{arm} contract key(s) {unaccounted} have no ParamOrigin, so nothing checks "
            "them against nextflow.config. Add an origin, or take the key out."
        )


#: Maps a keep-list argument of the pipeline's filter functions to the contract key that
#: supplies it. Used to derive arm membership from the real signatures instead of
#: asserting a hardcoded list.
_KEEP_ARG_TO_KEY = {
    "cancervar_keep": "filter_cancervar",
    "civic_keep": "filter_civic",
    "escat_keep": "filter_escat",
    "intervar_keep": "filter_intervar",
    "renovo_keep": "filter_renovo",
    "clinvar_keep": "filter_clinvar",
}

#: The vendored copy, not ``bin/``: it is always present (the packaged app ships it) and
#: ``test_vendor_drift.py`` already proves it is byte-for-byte the pipeline's, so
#: reading it is as authoritative as reading ``bin/`` and works in more trees.
VENDORED_FILTERS = STREAMLIT_APP / "vendor" / "pipeline_filters.py"


def _filter_arg_names(function: str) -> list:
    """Parameter names of a top-level function in the vendored filter module."""
    tree = ast.parse(VENDORED_FILTERS.read_text())
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name == function:
            return [arg.arg for arg in node.args.args]
    raise AssertionError(f"{function} not found in {VENDORED_FILTERS}")


@pytest.mark.skipif(not VENDORED_FILTERS.is_file(), reason="vendor/ not present")
@pytest.mark.parametrize(
    "arm,function",
    [("somatic", "somatic_filters"), ("germline", "germline_filters")],
)
def test_guideline_keep_lists_follow_the_pipeline_filter_signatures(arm, function):
    """Each arm carries exactly the keep-lists its pipeline filter function takes.

    Derived from ``somatic_filters()`` / ``germline_filters()`` rather than asserted
    from a hand-written list, so a signature change in the pipeline surfaces here. This
    is what pins ESCAT to the somatic arm: ``germline_filters()`` has no ESCAT argument,
    and the app's extra germline ESCAT clause is the largest divergence on real data
    (540 rows, issue #28).
    """
    from config.pipeline_params import PIPELINE_PARAMS

    args = _filter_arg_names(function)
    expected = {_KEEP_ARG_TO_KEY[a] for a in args if a in _KEEP_ARG_TO_KEY}
    assert expected, (
        f"recognised no keep-list argument in {function}{tuple(args)} — the signature "
        "changed and this test is no longer reading it"
    )

    present = {key for key in _KEEP_ARG_TO_KEY.values() if key in PIPELINE_PARAMS[arm]}
    assert present == expected, (
        f"the {arm} contract carries keep-lists {sorted(present)} but "
        f"{function}() takes {sorted(expected)}. Missing: "
        f"{sorted(expected - present)}; extra: {sorted(present - expected)}"
    )


@pytest.mark.skipif(not VENDORED_FILTERS.is_file(), reason="vendor/ not present")
@pytest.mark.parametrize(
    "arm,function",
    [("somatic", "somatic_filters"), ("germline", "germline_filters")],
)
def test_the_pages_guideline_sources_follow_the_pipeline_filter_signatures(arm, function):
    """The parameter page's copy of "which sources this arm ORs together" is guarded.

    ``config/vocabularies.GUIDELINE_SOURCES`` is a written-down copy of a
    relation that really lives in the vendored signatures — the page cannot import them
    without pulling pandas onto a parameter screen. Under this project's own rule (issue
    #28: "Where the app must agree with the pipeline about a list, derive it or guard it
    — never copy it"), a copy is only allowed with a guard, and this is it.

    It matters beyond tidiness. The page uses that mapping to decide when *every*
    guideline source has been emptied, which is the state issue #36 requires be warned
    about. A source missing from the copy is a source whose emptiness the warning ignores
    — the report silently collapses to the pathogenic-rescue floor with nothing said.

    Imported, not parsed. This guard used to read the relation out of
    ``page_modules/parameter_config.py`` with ``ast``, because importing a page pulls in
    streamlit and this guard has to run on a bare pytest. The relation now lives in
    ``config/``, which is streamlit-free on purpose, so the guard can read the object
    itself — and can no longer be fooled by a declaration whose *source* parses while the
    value the app actually uses is built some other way.
    """
    from config.vocabularies import GUIDELINE_SOURCES

    sources = {arm_name: set(keys) for arm_name, keys in GUIDELINE_SOURCES.items()}

    args = _filter_arg_names(function)
    expected = {_KEEP_ARG_TO_KEY[a] for a in args if a in _KEEP_ARG_TO_KEY}
    assert expected, (
        f"recognised no keep-list argument in {function}{tuple(args)} — the signature "
        "changed and this test is no longer reading it"
    )

    assert sources.get(arm) == expected, (
        f"config/vocabularies.py says the {arm} arm's guideline sources are "
        f"{sorted(sources.get(arm, []))}, but {function}() ORs together "
        f"{sorted(expected)}. Missing: {sorted(expected - sources.get(arm, set()))}; "
        f"extra: {sorted(sources.get(arm, set()) - expected)}. A missing source is one "
        "whose emptiness the all-empty warning will not notice."
    )


def test_escat_is_somatic_only():
    """ESCAT is absent from the germline contract, as it is from ``germline_filters``.

    The app's germline presets carry a ``filter_escat`` list, and the app ORs an
    ESCAT clause into its germline guideline block that the pipeline does not have —
    540 rows, the single largest divergence on real data (issue #28). Keeping the key
    out of the germline contract is what stops it being reintroduced by accident.
    """
    from config.pipeline_params import PIPELINE_PARAMS

    assert "filter_escat" in PIPELINE_PARAMS["somatic"]
    assert "filter_escat" not in PIPELINE_PARAMS["germline"], (
        "the pipeline's germline_filters() takes no ESCAT argument; a germline ESCAT "
        "keep-list in the contract is the 540-row divergence waiting to happen again"
    )


#: Every ESCAT level the pipeline's tiering table assigns, measured over all 172 of its rows.
#:
#: ``resources/escat_tiering.csv`` is the only source of the ``ESCAT`` column — the pipeline
#: reads it in ``bin/add_guidelines_escat_to_funcotator.py`` and matches gene, mutation and
#: tissue against it — so this is the whole vocabulary a MAF can carry, and therefore the whole
#: vocabulary the control has anything to offer. It is why issue #89 did not add ESCAT's ``IV``
#: and ``X``: the scale defines them, this table never assigns them, and a keep-list option
#: that cannot keep anything is the shape of issue #88.
#:
#: Frozen here as well as derived below, deliberately. ``resources/`` ships with the pipeline
#: and not with the packaged app, so a check that only reads the file asserts nothing wherever
#: the file is absent — the hole issue #77 found in the tests it deleted, where four assertions
#: sat behind a condition and skipped in silence. The literal runs unconditionally; the file,
#: where it is present, is checked against the literal so the literal cannot go stale unseen.
ESCAT_LEVELS_ASSIGNED = ["IA", "IB", "IC", "IIA", "IIB", "IIIA", "IIIB", "V"]


def test_the_escat_control_offers_exactly_the_levels_the_annotation_assigns():
    """No level the annotation cannot carry, and none it can that the control omits.

    Order matters as well as membership: it is ESCAT's own, strongest evidence first, and
    both the help page's list and the parameter tooltip's "from IA downwards" read it off
    ``ESCAT_OPTIONS``.
    """
    from config.vocabularies import ESCAT_OPTIONS

    assert ESCAT_OPTIONS == ESCAT_LEVELS_ASSIGNED, (
        "the ESCAT control's options no longer match the levels resources/escat_tiering.csv "
        "assigns; an added level keeps nothing, a dropped one hides real variants"
    )


@needs_escat_tiering
def test_the_frozen_escat_levels_still_match_the_tiering_table():
    """The literal above, re-derived from the file it was measured from.

    Guards the guard: without this, ``ESCAT_LEVELS_ASSIGNED`` would be a number typed once
    and never checked again, and the tiering table gaining a level — ESCAT's ``IV``, say, or
    the ``IIIC`` that ``bin/report_raw.py`` already excludes without any row carrying it —
    would leave the control silently unable to offer it.
    """
    import csv

    with ESCAT_TIERING.open(newline="") as handle:
        assigned = {row["ESCAT"].strip() for row in csv.DictReader(handle)}

    assert assigned == set(ESCAT_LEVELS_ASSIGNED), (
        f"resources/escat_tiering.csv now assigns {sorted(assigned)}, but "
        f"ESCAT_LEVELS_ASSIGNED froze {sorted(ESCAT_LEVELS_ASSIGNED)} — reconcile the "
        "control's options, the help page's definitions, and this literal together"
    )


def test_arms_carry_only_the_keys_their_filter_takes():
    """Somatic-only and germline-only guideline sources stay on their own arm."""
    from config.pipeline_params import PIPELINE_PARAMS

    somatic, germline = PIPELINE_PARAMS["somatic"], PIPELINE_PARAMS["germline"]

    for key in ("filter_cancervar", "filter_civic", "filter_escat"):
        assert key in somatic, key
        assert key not in germline, f"{key} is a somatic-only source"
    for key in ("filter_intervar", "filter_renovo"):
        assert key in germline, key
        assert key not in somatic, f"{key} is a germline-only source"

    assert somatic["sample_type"] == "somatic"
    assert germline["sample_type"] == "germline"


def test_variant_classification_carries_the_configs_three_exclusions(config_params):
    """The config's three values, with the pipeline's *exclude* meaning.

    ``filter_var_classification`` names what the pipeline drops. The app's parameter
    of the same name is an include list today (divergence #1), which is why this is
    worth stating explicitly: the contract holds three values, not the sixteen an
    include list would need.
    """
    from config.pipeline_params import PIPELINE_PARAMS

    expected = config_params["filter_var_classification"].split(",")
    assert expected == ["Silent", "IGR", "RNA"]
    for arm in ARMS:
        assert PIPELINE_PARAMS[arm]["filter_variant_classification"] == expected, arm


# ---------------------------------------------------------------------------
# The other copy of these values
# ---------------------------------------------------------------------------

#: ``tests/parity/contract.py``'s ``CONTRACT`` states the same pipeline defaults, in the
#: harness's own shape: flat rather than per-arm, both VAF keys rather than one, and
#: ``genes`` rather than ``filter_genes``. It has to keep that shape — it translates into
#: *both* sides, so it needs the two CLI names the pipeline actually takes.
#:
#: What it must not do is disagree. Two statements of the same defaults with no guard
#: between them is precisely the failure this ticket exists to end, so the mapping below
#: pins them together. Only the harness is guarded against the app here; the app is
#: guarded against ``nextflow.config`` directly by the tests above.
_HARNESS_KEY_TO_CONTRACT = {
    "min_depth": [("somatic", "min_depth"), ("germline", "min_depth")],
    "vaf_threshold": [("somatic", "vaf_threshold")],
    "vaf_threshold_germline": [("germline", "vaf_threshold")],
    "filter_variant_classification": [
        ("somatic", "filter_variant_classification"),
        ("germline", "filter_variant_classification"),
    ],
    "filter_cancervar": [("somatic", "filter_cancervar")],
    "filter_civic": [("somatic", "filter_civic")],
    "filter_escat": [("somatic", "filter_escat")],
    "filter_intervar": [("germline", "filter_intervar")],
    "filter_renovo": [("germline", "filter_renovo")],
    "filter_clinvar": [("somatic", "filter_clinvar"), ("germline", "filter_clinvar")],
    "skip_civic": [("somatic", "skip_civic"), ("germline", "skip_civic")],
    "skip_pathogenic": [("somatic", "skip_pathogenic"), ("germline", "skip_pathogenic")],
    "max_freq_population": [
        ("somatic", "max_freq_population"),
        ("germline", "max_freq_population"),
    ],
}


def test_the_parity_harness_contract_agrees_with_this_one():
    """The harness's ``CONTRACT`` and ``PIPELINE_PARAMS`` state the same defaults."""
    from config.pipeline_params import PIPELINE_PARAMS

    from .parity.contract import CONTRACT

    for harness_key, targets in _HARNESS_KEY_TO_CONTRACT.items():
        assert harness_key in CONTRACT, (
            f"tests/parity/contract.py's CONTRACT no longer has `{harness_key}`; this "
            "mapping is stale and the two contracts are no longer being compared"
        )
        for arm, key in targets:
            assert PIPELINE_PARAMS[arm][key] == CONTRACT[harness_key], (
                f"the parity harness and the app contract disagree: CONTRACT"
                f"[{harness_key!r}] is {CONTRACT[harness_key]!r} but "
                f"PIPELINE_PARAMS[{arm!r}][{key!r}] is {PIPELINE_PARAMS[arm][key]!r}. "
                "Two statements of the pipeline's defaults must not drift apart."
            )


def test_the_harness_contract_mapping_is_complete():
    """Every value-bearing harness key is compared, so none can drift unwatched."""
    from .parity.contract import CONTRACT

    # `genes` is the one key with no direct counterpart: the harness carries a fixture
    # *filename* (or None), the contract a list of symbols. Its neutral value is checked
    # by test_one_gene_list_key_not_a_somatic_germline_pair instead.
    exempt = {"genes"}
    unmapped = sorted(set(CONTRACT) - set(_HARNESS_KEY_TO_CONTRACT) - exempt)
    assert not unmapped, (
        f"tests/parity/contract.py's CONTRACT has key(s) {unmapped} that this guard does "
        "not compare against PIPELINE_PARAMS. Add them to _HARNESS_KEY_TO_CONTRACT."
    )


def test_the_harness_genes_default_is_neutral_like_the_contracts():
    """Both spellings of "no gene filter" agree, in their own shapes."""
    from config.pipeline_params import PIPELINE_PARAMS

    from .parity.contract import CONTRACT

    assert CONTRACT["genes"] is None
    for arm in ARMS:
        assert PIPELINE_PARAMS[arm]["filter_genes"] == []


# ---------------------------------------------------------------------------
# What this ticket deliberately does *not* change
# ---------------------------------------------------------------------------


def test_soft_and_clinical_presets_still_exist_as_deviations():
    """The four existing presets survive as opt-in deviations.

    None of them is at parity, and none is meant to be — so this checks the values that
    would have moved if one had been quietly edited into it.

    Issue #36 changed exactly one key in each: ``filter_variant_classification``, from the
    include list the app's own filter used to read to the exclude list the pipeline's
    parameter of that name actually means. That is a correction, not a deviation being
    softened — an unconverted preset would have said "drop every missense call". Its
    shape is asserted in ``tests/test_app_defaults.py``, which is where the polarity
    change is owned; here the point is only that the presets are still *not* the
    contract.
    """
    from config.presets import (
        CLINICAL_GERMLINE_PARAMS,
        CLINICAL_SOMATIC_PARAMS,
        SOFT_GERMLINE_PARAMS,
        SOFT_SOMATIC_PARAMS,
    )
    from config.pipeline_params import PIPELINE_PARAMS

    for name, preset in (
        ("SOFT_SOMATIC_PARAMS", SOFT_SOMATIC_PARAMS),
        ("SOFT_GERMLINE_PARAMS", SOFT_GERMLINE_PARAMS),
        ("CLINICAL_SOMATIC_PARAMS", CLINICAL_SOMATIC_PARAMS),
        ("CLINICAL_GERMLINE_PARAMS", CLINICAL_GERMLINE_PARAMS),
    ):
        assert isinstance(preset, dict) and preset, name
        # Still app-side vocabulary and still not neutral on frequency: both would
        # have had to change for a preset to have been quietly moved to parity.
        assert preset["keep_pathogenic"] is True, name
        assert preset["max_freq_population"] < 1.0, name
        assert "filter_escat" in preset, name

    assert SOFT_SOMATIC_PARAMS["vaf_threshold"] == 0.05
    assert SOFT_GERMLINE_PARAMS["vaf_threshold_germline"] == 0.2
    assert CLINICAL_SOMATIC_PARAMS["vaf_threshold"] == 0.05
    assert CLINICAL_GERMLINE_PARAMS["vaf_threshold_germline"] == 0.2

    # Each germline preset carries the two germline guideline sources, and each somatic one
    # the two somatic sources the germline filter does not take. Inherited from
    # `test_config.test_default_params_somatic_vs_germline`, which issue #77 deleted along
    # with the `DEFAULT_PARAMS` half it was really written around: this was the only
    # assertion that a *preset* names its arm's sources, and a preset that lost them would
    # arrive at its tab absent, become empty, and be written back — the silent widening
    # `test_app_defaults.test_switching_arms_opens_that_arm_at_parity` measured at a
    # criteria path of 13 against 27. That test covers the arm *switch*; nothing covered
    # the presets. Asserted here rather than back in `test_config.py` because this test
    # already owns "the presets are still deviations" and carries no `nextflow.config` gate.
    for name, preset in (
        ("SOFT_GERMLINE_PARAMS", SOFT_GERMLINE_PARAMS),
        ("CLINICAL_GERMLINE_PARAMS", CLINICAL_GERMLINE_PARAMS),
    ):
        for key in ("filter_intervar", "filter_renovo"):
            assert preset[key], (
                f"{name} no longer names {key}, so the germline tab opens without it"
            )
    for name, preset in (
        ("SOFT_SOMATIC_PARAMS", SOFT_SOMATIC_PARAMS),
        ("CLINICAL_SOMATIC_PARAMS", CLINICAL_SOMATIC_PARAMS),
    ):
        for key in ("filter_cancervar", "filter_civic"):
            assert preset[key], f"{name} no longer names {key}, so the somatic tab opens without it"
    # Still not the contract's three exclusions, which is what "deviation" means here.
    contract_exclusions = PIPELINE_PARAMS["somatic"]["filter_variant_classification"]
    for name, preset in (
        ("SOFT_SOMATIC_PARAMS", SOFT_SOMATIC_PARAMS),
        ("CLINICAL_SOMATIC_PARAMS", CLINICAL_SOMATIC_PARAMS),
    ):
        assert preset["filter_variant_classification"] != contract_exclusions, name


# `test_the_contract_is_not_yet_the_apps_default` stood here. It recorded the deferral
# issue #31 made — the contract written down but not yet loaded, because making it the
# default *before* the filter ran the pipeline's code would have changed the app's
# reported rows while it was still using its own re-implementation. Issue #33 routed the
# filter and issue #36 moved the default, so the deferral is spent and the test was
# #36's to delete, by its own instruction.
#
# What replaces it is `tests/test_app_defaults.py`, which asserts the positive: the app's
# default *is* the contract.


def test_the_contract_is_not_mutable_shared_state():
    """Callers get copies; a preset the UI edits must not rewrite the contract."""
    from config.pipeline_params import PIPELINE_PARAMS, pipeline_params

    first = pipeline_params("somatic")
    first["min_depth"] = 12345
    first["filter_clinvar"].append("Benign")

    assert pipeline_params("somatic")["min_depth"] == PIPELINE_PARAMS["somatic"]["min_depth"]
    assert "Benign" not in PIPELINE_PARAMS["somatic"]["filter_clinvar"], (
        "the nested keep-lists are shared with the contract, so a UI edit to a loaded "
        "preset would silently redefine parity"
    )

    with pytest.raises(ValueError):
        pipeline_params("tumour")
