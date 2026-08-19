"""Reading structure out of the vendored filter source.

Two modules derive a column contract from ``vendor/pipeline_filters.py`` rather than writing
one down — ``numeric_columns`` (which columns are compared as numbers, so an unreadable one is
a refusal) and ``absent_columns`` (which columns are read at all, so an absent one needs a
neutral fill). They ask different questions, but they ask them of the same two structures:
which parameter is the MAF, and which subscripts are columns of it.

Those two live here so there is one answer. They were ``numeric_columns``'s privates first;
the second caller is what turned them into a shared vocabulary rather than a helper. Neither
derivation is here: what counts as a *numeric* read and what counts as a read *on every path*
are the two modules' own subject matter, and merging them would produce a parser that neither
module could explain.

Everything here is deliberately conservative. A shape it does not understand yields nothing,
so a derivation built on it under-collects rather than inventing a column name — and both
callers turn an empty derivation into a startup failure, which is where a parser that has
stopped parsing becomes visible.
"""

from __future__ import annotations

import ast
from typing import Iterable

#: The annotation that marks a parameter as the MAF. Matched on the trailing name so that
#: ``pd.DataFrame``, ``pandas.DataFrame``, a bare ``DataFrame`` and the string form all count —
#: the vendored code writes the first, but which import style ``bin/`` uses is not something
#: these derivations should have an opinion about.
FRAME_ANNOTATION = "DataFrame"


def functions(tree: ast.AST) -> Iterable[ast.FunctionDef | ast.AsyncFunctionDef]:
    """Every function defined anywhere under ``tree``."""
    return (
        node
        for node in ast.walk(tree)
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
    )


def frame_parameters(function: ast.FunctionDef | ast.AsyncFunctionDef) -> set[str]:
    """The parameter names annotated as a DataFrame.

    Annotation-driven rather than name-driven, so nothing here has to know that the pipeline
    happens to call its frame ``maf``. An unannotated frame parameter yields nothing —
    ``compute_keep(args, out)`` is the real instance, and it neither compares a column nor is
    called with a MAF the app could fill, so the conservative reading costs nothing and a guess
    based on the name would be a rule about the pipeline's spelling habits.
    """
    args = function.args
    return {
        arg.arg
        for arg in [*args.posonlyargs, *args.args, *args.kwonlyargs]
        if annotation_name(arg.annotation) == FRAME_ANNOTATION
    }


def annotation_name(node: ast.AST | None) -> str | None:
    """The trailing name of an annotation: ``pd.DataFrame`` -> ``DataFrame``."""
    if isinstance(node, ast.Attribute):
        return node.attr
    if isinstance(node, ast.Name):
        return node.id
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        return node.value.rsplit(".", 1)[-1].strip("'\" ")
    return None


def is_frame_column(node: ast.Subscript, frames: set[str]) -> bool:
    """``maf["literal"]`` — a frame name subscripted by a string constant.

    A computed key (``maf[column]``) is not readable from the source, so it yields nothing. The
    pipeline has no such subscript today; if one appears, a derivation quietly stops covering
    it rather than inventing a name, and both callers' mutation tests pin that limit on the
    record rather than leaving it a surprise.
    """
    return (
        isinstance(node.value, ast.Name)
        and node.value.id in frames
        and isinstance(node.slice, ast.Constant)
        and isinstance(node.slice.value, str)
    )
