"""One rule for every column the variant panel reads by name (issue #212).

**The rule.** A column name that carries a character illegal in an R identifier can arrive
spelled two ways, because somewhere between the pipeline and the MAF a step passes the header
through ``make.names``. Every such name the variant panel reads must reach the row through
:func:`config.columns.spelled_in`, never by equality.

**Why a rule and not a list.** Issues #189, #190, #191, #198 and #210 between them measured
twenty-one columns one at a time, and the survey kept finding the same single offender. Issue
#212 measured the *shape* of the problem instead, over 167 byte-distinct real MAFs:

* **Mangling is all-or-nothing per file.** 3 files were mangled wholesale — every name with an
  illegal character rewritten — **0** files were mixed, 157 carried raw names only and 7 carried
  neither. So there is no per-column answer to *"is this column mangled?"*; a file either went
  through ``make.names`` or it did not, and then every exposed name in it moved at once.
* **So the exposed set is exactly the character class.** Of the 58 column names the panel reads,
  **3** contain such a character: the two evidence columns and ``GERP++_RS``. The other 55 cannot
  be rewritten — including ``AAChange.refGene``, ``Func.refGene`` and ``ExonicFunc.refGene``,
  because a dot is *legal* in an R identifier and ``make.names`` leaves those alone.

A list needs re-measuring every time a column is wired. The rule does not, which is what makes
these guards worth having: they fail on the day an exposed column is added, not on the day
someone next thinks to walk the corpus.

**Two guards, because a name reaches a row two ways.**

:func:`test_no_exposed_column_name_is_read_by_equality` reads the panel's own AST and fails a
literal sitting in a row-lookup position. That catches the direct read and cannot catch a name
that arrives through a variable — which is precisely how the ``GERP++_RS`` bug #190 found got in,
since ``_PREDICTOR_TABLE_CONFIG`` held the name and the lookup read a loop variable.

:func:`test_a_wholesale_mangled_file_renders_the_same_panel` closes that by asserting the
measured property directly: pass every column name through ``make.names`` and the panel must
render identically. It does not care how a name reaches a row, so a table-driven read is covered
by the same assertion as a literal one.

Every guard here was made to fail before being trusted, per this repo's standing rule — see
``docs/wayfinder/issue-212/README.md`` for what each mutation was and what it caught.
"""

import ast
import os
import sys
from pathlib import Path
from unittest.mock import MagicMock

import pandas as pd
import pytest

STREAMLIT_APP = Path(__file__).resolve().parents[1]
if str(STREAMLIT_APP) not in sys.path:
    sys.path.insert(0, str(STREAMLIT_APP))

from config.columns import is_mangling_exposed, spelled_in  # noqa: E402

CANCERVAR_FIXTURES = Path(__file__).resolve().parent / "fixtures" / "cancervar"

#: The modules the variant panel draws through. The panel itself, the two evidence sections, the
#: score-context tables and the reference scale — the surfaces map #184 owns.
#:
#: ``config/contaminated_columns.py`` is here because the panel's chromosome check reads columns by
#: name too, and it is the read that decides whether a number is drawn as a warning or as a score.
PANEL_MODULES = (
    "components/variant_detail.py",
    "components/acmg_evidence.py",
    "components/cbp_evidence.py",
    "components/cancervar_markers.py",
    "components/clinical_badges.py",
    "components/predictor_context.py",
    "components/alphamissense.py",
    "components/reference_scales.py",
    "config/contaminated_columns.py",
)


def _row_lookup_literals(source: str):
    """Every string literal sitting where the app looks a column up, as ``(lineno, literal)``.

    Three positions, and only three, because these are the ways a name indexes a row or a frame
    directly:

    * ``row.get("X")`` — including ``.get("X", default)``
    * ``row["X"]``
    * ``"X" in row.index`` / ``"X" in frame.columns``

    A literal handed to :func:`~config.columns.spelled_in` is deliberately **not** collected: that
    is the resolved path, and the whole point of the rule is that an exposed name goes there.

    Restricted to these positions rather than to every literal in the file, because the panel's
    modules are full of prose that ``make.names`` would happily rewrite — every ACMG criterion
    sentence, every tier name, every hex colour. A guard over all literals would report 60 lines of
    English and be turned off within a week.
    """
    hits = []

    class Walker(ast.NodeVisitor):
        def visit_Call(self, node):
            func = node.func
            if isinstance(func, ast.Attribute) and func.attr == "get" and node.args:
                first = node.args[0]
                if isinstance(first, ast.Constant) and isinstance(first.value, str):
                    hits.append((node.lineno, first.value))
            self.generic_visit(node)

        def visit_Subscript(self, node):
            index = node.slice
            if isinstance(index, ast.Constant) and isinstance(index.value, str):
                hits.append((node.lineno, index.value))
            self.generic_visit(node)

        def visit_Compare(self, node):
            for op, comparator in zip(node.ops, node.comparators):
                if not isinstance(op, (ast.In, ast.NotIn)):
                    continue
                if not isinstance(node.left, ast.Constant):
                    continue
                if not isinstance(node.left.value, str):
                    continue
                target = ast.unparse(comparator)
                if target.endswith((".index", ".columns")):
                    hits.append((node.lineno, node.left.value))
            self.generic_visit(node)

    Walker().visit(ast.parse(source))
    return hits


# ---------------------------------------------------------------------------
# The rule, read off the panel's own source
# ---------------------------------------------------------------------------


def test_no_exposed_column_name_is_read_by_equality():
    """A name ``make.names`` rewrites may not index a row directly.

    The failure this prevents is silent and plausible: an exact lookup of ``GERP++_RS`` on one of
    the 3 wholesale-mangled files reports the column *absent*, and the panel's absent path is a
    real, reachable state that says "not in this file" — so a present score reads as a missing one
    and nothing raises. Issue #190 met exactly that, and #194's contamination detector inherited it
    for the same reason.
    """
    offenders = []
    for relative in PANEL_MODULES:
        path = STREAMLIT_APP / relative
        for lineno, literal in _row_lookup_literals(path.read_text()):
            if is_mangling_exposed(literal):
                offenders.append(f"{relative}:{lineno} reads {literal!r} by equality")

    assert not offenders, (
        "these reads would miss the column on a file that passed through make.names; look the "
        "name up with config.columns.spelled_in instead:\n  " + "\n  ".join(offenders)
    )


def test_the_source_walk_actually_finds_the_panels_reads():
    """The guard above is only worth its green if it is reading the modules.

    A walk that collected nothing would pass :func:`test_no_exposed_column_name_is_read_by_equality`
    forever — this repo's recurring failure shape, and one #210 found sitting under a test that
    passed on a falsehood. So the walk is asserted to see the reads that are there, by name, on the
    module that has the most of them.
    """
    found = _row_lookup_literals((STREAMLIT_APP / "components/variant_detail.py").read_text())
    literals = {literal for _lineno, literal in found}

    assert len(found) >= 20, f"the walk found only {len(found)} row lookups in variant_detail.py"
    for expected in ("Hugo_Symbol", "Chromosome", "InterVar", "CancerVar", "tumor_f"):
        assert expected in literals, f"the walk no longer sees the {expected!r} read"

    # `am_class` used to be one of the names checked here and left the panel altogether in issue
    # #204, which also moved the clinical row's reads to their own module. The new module is
    # asserted rather than trusted, for the same reason: a module added to `PANEL_MODULES` whose
    # reads the walk cannot see is a walk that collects nothing, one module later.
    badge_literals = {
        literal
        for _lineno, literal in _row_lookup_literals(
            (STREAMLIT_APP / "components/clinical_badges.py").read_text()
        )
    }
    for expected in ("ClinVar_VCF_CLNSIG", "ESCAT"):
        assert expected in badge_literals, (
            f"the walk no longer sees the {expected!r} read in clinical_badges.py"
        )

    # `components/alphamissense.py` is deliberately **not** asserted the same way: it reads its one
    # column through the `AM_COLUMN` constant, so no literal sits in a lookup position and this
    # walk cannot see it at all. That is the walk's known blind spot — the shape of the `GERP++_RS`
    # bug issue #190 found — and the reason `test_a_wholesale_mangled_file_renders_the_same_panel`
    # exists beside it: the column is registered in `declared_read_columns` and covered there.


def test_the_rule_is_about_the_character_class_and_not_about_a_dot():
    """``AAChange.refGene`` is safe and ``GERP++_RS`` is not, and the guard must tell them apart.

    This is the distinction that makes the rule narrow enough to be true. A dot is legal in an R
    identifier, so ``make.names`` leaves a dotted ANNOVAR name alone; treating "contains a dot" as
    the danger would flag three of the panel's reads that cannot move, and treating "contains
    punctuation" as safe would miss the one that does.
    """
    for safe in ("AAChange.refGene", "Func.refGene", "ExonicFunc.refGene", "Polyphen2_HDIV_pred"):
        assert not is_mangling_exposed(safe)
    for exposed in ("GERP++_RS", "M-CAP_score", "CancerVar: CancerVar and Evidence"):
        assert is_mangling_exposed(exposed)


# ---------------------------------------------------------------------------
# The rule, applied to the panel as a whole
# ---------------------------------------------------------------------------


def make_names(name: str) -> str:
    """R's ``make.names`` over one column name, as the corpus shows it was actually applied.

    Written here rather than imported from ``config.columns`` on purpose. The module under test
    owns the *canonical* direction — reduce a name to its normal form — and re-using it would make
    this a comparison of ``spelled_in`` with itself, which is issue #77's shape and issue #100's.
    This is the pipeline's direction: take a header and mangle it, the way the 3 real files were.

    **The strip is not a convenience, it is the observed behaviour**, and writing it the obvious way
    first got it wrong. Mangling the padding as well turns ``' CancerVar: CancerVar and Evidence '``
    into ``'.CancerVar..CancerVar.and.Evidence.'``, and no real MAF spells anything that way: across
    167 files and 727 distinct column names, **0** begin with a dot. The 3 wholesale-mangled files
    carry ``'CancerVar..CancerVar.and.Evidence'``, so the padding was gone before ``make.names``
    ran. :func:`test_the_mangler_reproduces_the_spelling_the_real_files_carry` pins that against the
    two fixtures rather than leaving it as a claim in a docstring.
    """
    import re

    return re.sub(r"[^A-Za-z0-9_.]", ".", name.strip())


def test_the_mangler_reproduces_the_spelling_the_real_files_carry():
    """The guard below is only as good as its model of what a mangled header looks like.

    Both fixtures are real column names from real MAFs — the padded spelling 93 of 96 files use and
    the dot-mangled one 3 of them use — so this asks the mangler to turn the first into the second.
    A mangler that got this wrong would still make the comparison guard fail or pass; it would just
    be failing about a header no file has.
    """
    padded = " CancerVar: CancerVar and Evidence "
    dotted = "CancerVar..CancerVar.and.Evidence"

    frames = {
        name: _fixture_frame(name)
        for name in ("somatic_cancervar_evidence.maf", "somatic_cancervar_dotted_no_column.maf")
    }
    assert padded in list(frames["somatic_cancervar_evidence.maf"].columns)
    assert dotted in list(frames["somatic_cancervar_dotted_no_column.maf"].columns)

    assert make_names(padded) == dotted, (
        f"the mangler produces {make_names(padded)!r}, which no real MAF carries"
    )
    assert make_names("GERP++_RS") == "GERP.._RS"


def _panel_text(monkeypatch, row: pd.Series) -> str:
    """Everything :func:`render_variant_detail_panel` draws, as one string.

    A fake ``st`` rather than ``AppTest``, because the claim under test is *what text was drawn*
    for two headers of the same file — and each of the panel's tables is a single ``st.markdown``
    whose content element counts cannot see inside (issue #188). ``session_state`` is a real dict
    so the panel's three session reads take their documented absent paths instead of receiving a
    ``MagicMock`` that answers every question.
    """
    from components import acmg_evidence, cbp_evidence, predictor_context, reference_scales
    from components import variant_detail

    drawn = []

    class _Ctx:
        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

    fake = MagicMock()
    fake.session_state = {}
    fake.columns.side_effect = lambda n, *a, **k: [
        _Ctx() for _ in range(n if isinstance(n, int) else len(n))
    ]
    fake.expander.side_effect = lambda *a, **k: _Ctx()
    for method in ("markdown", "caption", "info", "warning", "error", "html"):
        getattr(fake, method).side_effect = (
            lambda text, *a, **k: drawn.append(str(text))
        )
    fake.metric.side_effect = lambda label, value, *a, **k: drawn.append(f"{label}: {value}")

    for module in (
        variant_detail,
        cbp_evidence,
        acmg_evidence,
        predictor_context,
        reference_scales,
    ):
        monkeypatch.setattr(module, "st", fake)

    variant_detail.render_variant_detail_panel(row)
    return "\n".join(drawn)


def _fixture_frame(name: str) -> pd.DataFrame:
    return pd.read_csv(CANCERVAR_FIXTURES / name, sep="\t", low_memory=False)


def declared_read_columns():
    """Every MAF column the variant panel declares it reads, from the panel's own tables.

    The tables are named here; the **columns are not**. That is the whole distinction from
    ``test_contaminated_columns.py``'s registry: adding a column to any table below is covered by
    the guards without touching this file, which is what "a rule instead of a list" buys. What is
    *not* covered is a brand-new table in a new module — so each name below is dereferenced rather
    than looked up defensively, and a rename fails here with an ``AttributeError`` instead of
    quietly shrinking the set.
    """
    from components import (
        alphamissense,
        cancervar_markers,
        clinical_badges,
        predictor_context,
        reference_scales,
        variant_detail,
    )
    from config import contaminated_columns

    columns = set()

    # The three marker columns behind CBP1/CBP2/CBP3, and the sample's own tumour type. Issue #212
    # shipped naming this module as a *limit* — a new table of column names is not covered until
    # someone registers it — because #198 was still an open PR when the reads were enumerated. It
    # merged before this branch did, so the limit is discharged here rather than left as a note.
    # ``resolve_markers`` is handed ``row.get`` and calls it with each of these by exact name, which
    # is a read the AST walk cannot see: the lookup is on a variable, in a different module from the
    # table. All four names are safe, so this registers coverage rather than fixing a violation.
    columns |= {column for _code, column, _type, _noun in cancervar_markers.CRITERION_COLUMNS}
    columns.add("tumor_tissue")

    # The ClinGen reference scale: (display, column, direction, pp3, bp4), column may be None.
    columns |= {entry[1] for entry in reference_scales.CLINGEN_SVI_TABLE if entry[1]}
    # The inputs each classifier's own criterion reads, as that section names them.
    columns |= {
        column
        for inputs in reference_scales.CLASSIFIER_INPUTS.values()
        for _name, column in inputs
    }
    # The clinical row's three badges and the population-frequency chart.
    columns |= {entry[0] for entry in clinical_badges.CLINICAL_ROW}
    columns |= {entry[0] for entry in variant_detail._POP_FREQ_COLUMNS}
    # The three qualifier columns that ride on those badges rather than making claims of their
    # own (issue #200's Corollary A): the review status ClinVar's significance carries, RENOVO's
    # score, and the two spellings of ESCAT's tumour type. None of them is in a table — each is
    # read by name inside its own badge builder — so registering them is what puts them inside
    # the mangled-render comparison below.
    columns |= {
        "ClinVar_VCF_CLNREVSTAT",
        "RENOVO_pls",
        "ESCAT_TISSUE",
        "ESCAT_CANCER",
    }
    # AlphaMissense's own section, which reads one column (issue #203) — and deliberately not the
    # dbNSFP ``AlphaMissense_score`` beside it, which is a different annotation.
    columns.add(alphamissense.AM_COLUMN)
    # Both guideline verdict columns and both evidence columns.
    for classifier in variant_detail._CLASSIFIERS:
        columns.add(classifier.column)
        columns.add(classifier.evidence_column)
    # The score-context tables, read through the resolver already.
    empty = pd.Series(dtype=object)
    columns |= {r.column for r in predictor_context._cbp10_readings(empty, frozenset())}
    columns |= {r.column for r in predictor_context._pp3_bp4_readings(empty, frozenset())}
    # The columns the chromosome check measures.
    columns |= set(contaminated_columns.PREDICTOR_COLUMNS)

    return {column for column in columns if column}


#: A value that is legal in every column the panel reads, so one filler serves them all. The
#: comparison below is raw-header against mangled-header on the *same* values, so what the cell
#: says does not matter — only that the column is present in both and that neither render raises.
_FILLER = "0.5"


def _row_carrying_every_declared_column(fixture_row: pd.Series) -> pd.Series:
    """The fixture's row, widened to carry every column :func:`declared_read_columns` names.

    Widened rather than synthesised from scratch, so the evidence strings, the consequence and the
    identity fields stay real — the panel's arm gate and both parsers have to run for the render to
    be worth comparing.
    """
    row = fixture_row.copy()
    for column in sorted(declared_read_columns()):
        if spelled_in(row.index, column) is None:
            row[column] = _FILLER
    return row


def test_a_wholesale_mangled_file_renders_the_same_panel(monkeypatch):
    """The measured property, asserted directly: mangling the header must change nothing.

    This is the guard that does not depend on *how* a column name reaches a row, so it covers the
    table-driven reads the AST walk cannot see — ``CLINGEN_SVI_TABLE``'s columns,
    ``_CLINICAL_BADGE_CONFIG``'s five, ``_POP_FREQ_COLUMNS``' five and the score-context tables'
    nine. Add ``M-CAP_score`` to any of them behind an exact-name lookup and this fails while the
    AST walk stays green.

    **The row is widened to carry every declared column first, and that is load-bearing.** Driven
    off the fixture alone, this guard could only notice a column the fixture happens to have: a new
    exposed column would read as absent under *both* headers, the renders would match, and the guard
    would pass on the mutation it exists to catch. Written the narrow way first, it did exactly that.

    Both renders are asserted non-trivial, and the header is asserted to actually move. Two empty
    strings are equal, and a mangling that changed no name would compare a render with itself — the
    vacuous-guard shape this repo keeps meeting.
    """
    fixture_row = _fixture_frame("somatic_cancervar_evidence.maf").iloc[5]
    row = _row_carrying_every_declared_column(fixture_row)

    moved = [name for name in row.index if make_names(name) != name]
    assert moved, (
        "mangling this header changes no name, so the comparison below is a render against itself"
    )

    raw = _panel_text(monkeypatch, row)
    mangled = _panel_text(monkeypatch, row.rename(index=make_names))

    assert "CancerVar" in raw and len(raw) > 500, (
        "the unmangled render is too thin to be comparing anything: " + raw[:200]
    )
    assert mangled == raw, (
        "the panel reads differently on a wholesale-mangled header. The columns that moved:\n  "
        + "\n  ".join(f"{name!r} -> {make_names(name)!r}" for name in moved)
    )


def test_every_declared_column_is_actually_in_that_comparison():
    """The row the comparison runs over must carry every column the panel declares it reads.

    Asserted over the **row**, not over :func:`declared_read_columns`, and mutation testing is what
    forced that. Checking only the declared set left the widening itself unguarded: the fixture
    already carries both exposed names, so deleting the widening kept the render comparison green
    and passing — and quietly took with it the one thing the widening is for, catching an exposed
    column the fixture does *not* have. The two mutations only fail together, so one of them had to
    be turned into an assertion.
    """
    fixture_row = _fixture_frame("somatic_cancervar_evidence.maf").iloc[5]
    row = _row_carrying_every_declared_column(fixture_row)
    declared = declared_read_columns()

    assert len(declared) >= 25, f"only {len(declared)} declared columns found: {sorted(declared)}"

    missing = sorted(name for name in declared if spelled_in(row.index, name) is None)
    assert not missing, (
        f"the render comparison is driven over a row with no {missing} column, so a mangled "
        "spelling of it would read as absent under both headers and the comparison would pass"
    )

    exposed = sorted(name for name in declared if is_mangling_exposed(name))
    assert exposed == [
        "CancerVar: CancerVar and Evidence",
        "GERP++_RS",
        "InterVar: InterVar and Evidence",
    ], f"the exposed set the panel declares has changed: {exposed}"


# ---------------------------------------------------------------------------
# One mechanism, not two
# ---------------------------------------------------------------------------


def test_the_evidence_column_is_resolved_by_the_one_resolver():
    """Issue #212's second decision, held against the source that has to keep it.

    The rule the guards above state — *an exposed name reaches the row through*
    :func:`~config.columns.spelled_in` — is only checkable if there is one way for such a name to
    get there. The substring match this replaced was a second way: it violated no rule, resolved the
    padded spelling correctly, and would have gone on doing so, which is exactly why nothing else
    here can notice it coming back.
    """
    import inspect
    import textwrap

    from components import variant_detail

    tree = ast.parse(textwrap.dedent(inspect.getsource(variant_detail._evidence_column)))
    called = {
        node.func.id
        for node in ast.walk(tree)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
    }

    assert "spelled_in" in called, (
        "_evidence_column no longer resolves through config.columns.spelled_in, so the panel has "
        "two ways to find a column whose name make.names would rewrite and the rule the other "
        "guards state is no longer the whole of it"
    )
    assert not [
        node for node in ast.walk(tree) if isinstance(node, (ast.For, ast.comprehension))
    ], "_evidence_column scans the header itself again instead of asking the resolver"


def test_the_dot_mangled_fixture_is_read_as_the_evidence_column(monkeypatch):
    """#189's fixture, now resolved by the rule rather than by a substring match.

    ``somatic_cancervar_dotted_no_column.maf`` carries ``CancerVar..CancerVar.and.Evidence`` and no
    ``CancerVar`` column — the shape of the 3 real files — so the tier the panel badges can only
    come from the evidence string, and it can only get there if the mangled name resolves.
    """
    frame = _fixture_frame("somatic_cancervar_dotted_no_column.maf")
    text = _panel_text(monkeypatch, frame.iloc[0])

    assert "Tier" in text, "the mangled evidence column was not read at all"
    assert "AMP/ASCO/CAP" in text


# ---------------------------------------------------------------------------
# The resolver the rule points at
# ---------------------------------------------------------------------------


def test_the_resolver_handles_the_padding_the_evidence_columns_arrive_with():
    """The strip issue #212 added, and the measurement that required it.

    Without it the resolver failed on 93 of the 96 real files carrying CancerVar's evidence column
    and 56 of 56 carrying InterVar's, because a space is a character ``make.names`` rewrites and
    the padding normalised to leading and trailing dots.
    """
    canonical = "CancerVar: CancerVar and Evidence"

    assert spelled_in([" CancerVar: CancerVar and Evidence "], canonical) == (
        " CancerVar: CancerVar and Evidence "
    ), "the padded spelling 93 of 96 real files use does not resolve"
    assert spelled_in(["CancerVar..CancerVar.and.Evidence"], canonical) == (
        "CancerVar..CancerVar.and.Evidence"
    ), "the dot-mangled spelling 3 real files use does not resolve"
    assert spelled_in(["Chromosome", "CIViC_Evidence_Level"], canonical) is None


def test_what_comes_back_is_the_key_the_row_is_indexed_with():
    """The strip is for the comparison only — the padding has to survive into the answer.

    Returning the canonical name would raise ``KeyError`` on 93 of 96 files, which is a louder
    failure than the one this fixes but no more correct.
    """
    padded = " CancerVar: CancerVar and Evidence "
    row = pd.Series({padded: "CancerVar: 8#Tier_II_potential EVS=[1, 0, 2, 1, 0, 0, 0, 1, 2, 1, 1, 1]"})
    resolved = spelled_in(row.index, "CancerVar: CancerVar and Evidence")

    assert resolved == padded
    assert row[resolved]  # the point: the answer is usable as an index key


def test_a_tie_is_still_refused_after_the_strip():
    """The strip widens what matches, so the refusal has to be re-checked against it.

    ``GERP++_NR`` is a real ANNOVAR column, and guessing between two candidates would put a
    different conservation score under the label of the one CBP10 and PP3 read.
    """
    assert spelled_in(["GERP.._RS", "GERP+-_RS"], "GERP++_RS") is None
    assert spelled_in([" GERP.._RS ", "GERP+-_RS"], "GERP++_RS") is None
    assert spelled_in(["GERP.._RS", "GERP.._NR"], "GERP++_RS") == "GERP.._RS"


def test_the_canonical_spelling_wins_when_a_file_carries_both():
    """No real MAF does, and the tie-break is still not allowed to be arbitrary."""
    assert spelled_in(["GERP.._RS", "GERP++_RS"], "GERP++_RS") == "GERP++_RS"
