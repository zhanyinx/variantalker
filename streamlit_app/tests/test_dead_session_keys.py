"""No branch in the app waits on a session key that nothing ever sets.

A guard rather than a unit module: nothing here calls the app or imports Streamlit. It parses
the source of every app module and asserts a *shape*, so it holds in a checkout with no
dependencies installed at all.

Issue #167 is why. ``page_modules/help.py`` read a ``help_tab_focus`` key and branched on two
values, each drawing an ``st.info`` naming the Help tab the reader should click, then clearing
the request. **Nothing in the app ever wrote it** — that read and its two clears were the only
three mentions anywhere — so both sentences were unreachable from every route onto the page,
and had been since before the button that might once have written one was deleted at issue
#161. Worse than unreachable: #161 read the code, concluded the mechanism existed, and had to
go and measure ``st.tabs`` to find out that Python cannot select a tab at all. Dead code
shaped like a feature costs the next author a measurement, which is the whole reason this file
exists rather than a comment.

**Why a read-with-no-writer guard would not have caught it**, and what this file checks
instead. ``help_tab_focus`` *was* written — twice, by the two ``= None`` clears inside the very
branches its read gates. A rule about reads and writes is satisfied by a closed loop. So the
rule here is about writes that put a **value** in:

* assigning something that is not the key's own off-value;
* a widget owning the key through ``key="..."``, which is Streamlit writing the user's input;
* **mutating the container the key holds** — ``st.session_state.notes[vkey] = ...``,
  ``.update(...)``, ``.append(...)``.

That third clause is load-bearing and is asserted to be, below. ``variant_notes`` and
``custom_annotations`` are initialised to ``{}`` and never assigned anything else; every note
a user writes arrives by subscript into that dict. Without container mutation counting as a
write they would read exactly like the dead key, and a guard whose first run flags two live
features gets an allowlist on day one and means nothing by day two.

The converse clause is what separates the dead flag this file found in passing from the live
latch two files over. ``auto_load_checked`` is initialised to ``True`` and read only by its own
``not in`` guard — a once-only latch whose meaning *is* presence, and alive. ``params_cached``
was initialised to ``False`` under the comment "Initialize cache tracking", and those two lines
were its only mentions in the repository: nothing set it when parameters were cached, nothing
read it to find out. A flag only ever initialised to its off value and never turned on is a
flag nothing sets, and issue #167 deleted it alongside the key it came for.

**Parsed, not grepped.** Deliberate, and not a style preference: both deletions left a comment
naming the key they removed, so a text rule over these files would be satisfied by the
explanation of the thing it is looking for. This repo has already shipped one guard satisfied
by the filename in its own header comment (issue #162). An AST sees only the code that runs.
"""

from __future__ import annotations

import ast
from pathlib import Path

APP_ROOT = Path(__file__).resolve().parent.parent

#: Not the app's own code. ``vendor/`` is a byte-for-byte copy of the pipeline's filter code
#: and is guarded against drift by ``test_vendor_drift.py`` — editing it to satisfy this file
#: would break that one. ``tests/`` drives session state with fakes, where a key written by no
#: app module is the normal case.
SKIP_DIRS = {"vendor", "tests", "__pycache__", "parity"}

#: Mapping methods reached *on* ``st.session_state`` itself, so ``st.session_state.get`` is not
#: mistaken for a read of a key called ``get``. A session key sharing one of these names would
#: be invisible here; naming one ``update`` is its own problem.
STATE_METHODS = {"get", "pop", "setdefault", "update", "keys", "values", "items", "clear"}

#: Methods that, called on the *object a key holds*, put something in it. ``pop``, ``clear``,
#: ``remove`` and ``discard`` are deliberately absent: taking a note out of a dict is no
#: evidence that anything ever puts one in, which is precisely the mistake that would have
#: passed ``help_tab_focus``.
FILLING_METHODS = {"append", "extend", "insert", "update", "setdefault", "add", "sort"}


def _is_session_state(node: ast.AST) -> bool:
    """``<name>.session_state`` — the receiver is not pinned to ``st``, only its attribute."""
    return (
        isinstance(node, ast.Attribute)
        and node.attr == "session_state"
        and isinstance(node.value, ast.Name)
    )


def _string(node: ast.AST) -> str | None:
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        return node.value
    return None


def _direct_key(node: ast.AST) -> str | None:
    """The key ``node`` names on session state, for the two equivalent spellings.

    ``st.session_state.foo`` and ``st.session_state["foo"]``. Both are in use in this app, and
    a rule that recognised one of them would be walked past by writing the other.
    """
    if isinstance(node, ast.Attribute) and _is_session_state(node.value):
        return None if node.attr in STATE_METHODS else node.attr
    if isinstance(node, ast.Subscript) and _is_session_state(node.value):
        return _string(node.slice)
    return None


def _container_key(node: ast.AST) -> str | None:
    """The key whose *contents* ``node`` addresses, unwrapping index and attribute chains.

    ``st.session_state.custom_annotations[col][vkey]`` addresses the contents of
    ``custom_annotations`` two levels down, and nesting is how the annotation store is
    actually written.
    """
    while isinstance(node, (ast.Subscript, ast.Attribute)):
        key = _direct_key(node)
        if key is not None:
            return key
        node = node.value
    return None


def _is_off_value(node: ast.AST) -> bool:
    """Whether an assigned value only puts the key in its own empty or off state.

    ``None``, ``False``, ``0``, ``""``, an empty container literal, or a bare ``dict()`` /
    ``list()`` / ``set()``. Note ``True`` is a value: a latch initialised to ``True`` is
    carrying a fact, which is what keeps ``auto_load_checked`` off the offender list.
    """
    if isinstance(node, ast.Constant):
        return node.value is None or node.value is False or node.value == 0 or node.value == ""
    if isinstance(node, (ast.List, ast.Set)):
        return not node.elts
    if isinstance(node, ast.Dict):
        return not node.keys
    if isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
        return node.func.id in {"dict", "list", "set"} and not node.args and not node.keywords
    return False


class _Survey:
    """Every literal-keyed use of session state across the app, sorted by what it proves."""

    def __init__(self) -> None:
        self.modules: list[str] = []
        self.reads: dict[str, list[str]] = {}
        self.filled: dict[str, list[str]] = {}  # a value goes in
        self.emptied: dict[str, list[str]] = {}  # only ever cleared or initialised off
        self.widgets: dict[str, list[str]] = {}  # Streamlit owns it via key=
        self.mutated: dict[str, list[str]] = {}  # the container it holds gains contents

    def add(self, bucket: dict[str, list[str]], key: str, site: str) -> None:
        bucket.setdefault(key, []).append(site)

    @property
    def written(self) -> set[str]:
        return set(self.filled) | set(self.emptied) | set(self.widgets) | set(self.mutated)

    def is_set_somewhere(self, key: str) -> bool:
        return key in self.filled or key in self.widgets or key in self.mutated


class _Visitor(ast.NodeVisitor):
    def __init__(self, survey: _Survey, module: str) -> None:
        self.survey = survey
        self.module = module

    def site(self, node: ast.AST) -> str:
        return f"{self.module}:{node.lineno}"

    # -- writes ------------------------------------------------------------------

    def visit_Assign(self, node: ast.Assign) -> None:
        for target in node.targets:
            self._target(target, node.value)
        self.generic_visit(node)

    def visit_AugAssign(self, node: ast.AugAssign) -> None:
        self._target(node.target, node.value)
        self.generic_visit(node)

    def visit_AnnAssign(self, node: ast.AnnAssign) -> None:
        if node.value is not None:
            self._target(node.target, node.value)
        self.generic_visit(node)

    def _target(self, target: ast.AST, value: ast.AST) -> None:
        if isinstance(target, (ast.Tuple, ast.List)):
            for element in target.elts:
                self._target(element, value)
            return
        key = _direct_key(target)
        if key is not None:
            bucket = self.survey.emptied if _is_off_value(value) else self.survey.filled
            self.survey.add(bucket, key, self.site(target))
            return
        # `st.session_state.notes[vkey] = ...` — the container gains an entry, whatever the
        # entry is, so the value is not consulted here.
        if isinstance(target, ast.Subscript):
            container = _container_key(target.value)
            if container is not None:
                self.survey.add(self.survey.mutated, container, self.site(target))

    def visit_Delete(self, node: ast.Delete) -> None:
        for target in node.targets:
            key = _direct_key(target)
            if key is not None:
                self.survey.add(self.survey.emptied, key, self.site(target))
        self.generic_visit(node)

    # -- reads, and the calls that are both -------------------------------------

    def visit_Call(self, node: ast.Call) -> None:
        func = node.func
        if isinstance(func, ast.Attribute):
            if _is_session_state(func.value) and func.attr in STATE_METHODS:
                self._state_method(node, func)
            elif func.attr in FILLING_METHODS:
                container = _container_key(func.value)
                if container is not None:
                    self.survey.add(self.survey.mutated, container, self.site(node))
        elif isinstance(func, ast.Name) and func.id in {"getattr", "hasattr"} and node.args:
            if _is_session_state(node.args[0]) and len(node.args) > 1:
                key = _string(node.args[1])
                if key is not None:
                    self.survey.add(self.survey.reads, key, self.site(node))
        for keyword in node.keywords:
            if keyword.arg == "key":
                key = _string(keyword.value)
                if key is not None:
                    self.survey.add(self.survey.widgets, key, self.site(node))
        self.generic_visit(node)

    def _state_method(self, node: ast.Call, func: ast.Attribute) -> None:
        key = _string(node.args[0]) if node.args else None
        if key is None:
            return
        if func.attr in {"get", "pop", "setdefault"}:
            self.survey.add(self.survey.reads, key, self.site(node))
        if func.attr == "pop":
            self.survey.add(self.survey.emptied, key, self.site(node))
        elif func.attr == "setdefault":
            default = node.args[1] if len(node.args) > 1 else None
            bucket = (
                self.survey.emptied
                if default is None or _is_off_value(default)
                else self.survey.filled
            )
            self.survey.add(bucket, key, self.site(node))

    def visit_Attribute(self, node: ast.Attribute) -> None:
        if (
            _is_session_state(node.value)
            and isinstance(node.ctx, ast.Load)
            and node.attr not in STATE_METHODS
        ):
            self.survey.add(self.survey.reads, node.attr, self.site(node))
        self.generic_visit(node)

    def visit_Subscript(self, node: ast.Subscript) -> None:
        if _is_session_state(node.value) and isinstance(node.ctx, ast.Load):
            key = _string(node.slice)
            if key is not None:
                self.survey.add(self.survey.reads, key, self.site(node))
        self.generic_visit(node)

    def visit_Compare(self, node: ast.Compare) -> None:
        # `"foo" in st.session_state` — the init guard almost every key is created behind.
        for operator, comparator in zip(node.ops, node.comparators):
            if isinstance(operator, (ast.In, ast.NotIn)) and _is_session_state(comparator):
                key = _string(node.left)
                if key is not None:
                    self.survey.add(self.survey.reads, key, self.site(node))
        self.generic_visit(node)


def _survey() -> _Survey:
    survey = _Survey()
    for path in sorted(APP_ROOT.rglob("*.py")):
        relative = path.relative_to(APP_ROOT)
        if set(relative.parts) & SKIP_DIRS:
            continue
        module = relative.as_posix()
        survey.modules.append(module)
        _Visitor(survey, module).visit(ast.parse(path.read_text(encoding="utf-8"), filename=module))
    return survey


def test_every_session_key_the_app_reads_is_set_somewhere():
    """A branch waiting on a key nothing fills is a branch that cannot run."""
    survey = _survey()
    offenders = {
        key: sites for key, sites in survey.reads.items() if not survey.is_set_somewhere(key)
    }

    assert not offenders, "\n".join(
        [
            "These session keys are read but nothing ever puts a value in them, so every "
            "branch behind them is unreachable:",
            *(
                f"  {key!r} read at {', '.join(sites)}"
                + (
                    f" — cleared at {', '.join(survey.emptied[key])}, and nowhere else written"
                    if key in survey.emptied
                    else " — never written at all"
                )
                for key, sites in sorted(offenders.items())
            ),
            "",
            "This is issue #167's shape: `help_tab_focus` was read on the Help page to choose "
            "between two hints, and its only writers were the two clears inside the branches "
            "that read it. Either give the key a writer that means something, or delete the "
            "read and the branches with it — leaving them costs the next author the "
            "measurement it cost issue #161.",
        ]
    )


def test_no_session_key_is_written_and_never_read():
    """The same deadness from the other side: a value stored and never consulted.

    Widget-owned keys are excluded, and have to be: ``st.button(key="recap_reapply")`` is read
    by Streamlit, not by this app, and eleven of them are live today.
    """
    survey = _survey()
    offenders = {
        key: sorted(
            set(
                survey.filled.get(key, [])
                + survey.emptied.get(key, [])
                + survey.mutated.get(key, [])
            )
        )
        for key in survey.written
        if key not in survey.reads and key not in survey.widgets
    }

    assert not offenders, "\n".join(
        [
            "These session keys are written and never read back:",
            *(
                f"  {key!r} written at {', '.join(sites)}"
                for key, sites in sorted(offenders.items())
            ),
            "",
            "Nothing downstream asks for them, so the write is doing no work. `params_cached` "
            "was this, initialised to `False` under a comment about tracking the parameter "
            "cache while the real tracking happened on disk in `param_store.py`; issue #167 "
            "deleted it.",
        ]
    )


def test_the_survey_has_something_to_say():
    """Vacuity, asserted rather than assumed — including for each clause that grants life.

    Both rules above are "no offenders" rules, and a parse that found nothing would satisfy
    them in silence. The three clauses that let a key count as set are each pinned to a live
    example here, because the one that matters most is invisible from the rules themselves:
    strike container mutation out of this file and it flags ``variant_notes`` and
    ``custom_annotations``, the two stores every note and annotation the user types goes into.
    """
    survey = _survey()

    assert len(survey.modules) > 20, (
        f"Only {len(survey.modules)} app modules parsed. The rules here are quantified over "
        "what this walk finds, so a broken walk reads exactly like a clean app."
    )
    assert len(survey.reads) > 15, f"Only {len(survey.reads)} session keys read across the app."
    assert survey.filled, "No session key is assigned a value anywhere. The parse is broken."
    assert survey.widgets, "No widget owns a session key through `key=`. The parse is broken."

    kept_alive_only_by_mutation = {
        key
        for key in survey.reads
        if key not in survey.filled and key not in survey.widgets and key in survey.mutated
    }
    assert kept_alive_only_by_mutation, (
        "No session key is kept alive by container mutation alone, so the clause that "
        "recognises it is no longer doing any work — and the next store initialised to `{}` "
        "and filled by subscript will be reported dead. `variant_notes` and "
        "`custom_annotations` were both this at issue #167."
    )
