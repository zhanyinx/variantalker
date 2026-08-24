"""No script that requires bash may be invoked through `sh`.

`update_db` runs inside an Ubuntu 22.04 image where `/bin/sh` is **dash**. Every `.sh` file here
is a bash script, so handing one to `sh` makes it die on a syntax error -- and `utils_update.py`
does not inspect the result. Two failure modes, and the difference decides the blast radius:

* `function name {` is a **parse** error: dash exits 2 having run nothing.
* `[[ ... ]]` is a **runtime** command-not-found: the script *continues*, with every conditional
  false. `update_clinvar_funcotator.sh` is the only live callee with no `function` keyword, so it
  takes this path and exits **0** having done nothing but `mkdir clinvar`. That empty directory,
  copied over the live annotation database, is the defect this suite exists for.

**Why this half is static.** On macOS `/bin/sh` is GNU bash 3.2, which accepts both bashisms
happily -- so this entire bug class is invisible on a development machine and manifests only in
the container. Reproducing it dynamically would mean a 2.24 GB amd64 image under emulation,
which is not a guard anyone will run. Checking it statically runs everywhere, including on the
machines where the bug cannot be reproduced at all, which is strictly better coverage. It is
also why the container's `/bin/sh` is deliberately left as dash: remapping it would erase the
one environment where this is observable.

Run as a report:  python3 update_db/tests/test_shell_interpreters.py
Run as a guard:   python -m pytest update_db/tests/test_shell_interpreters.py
"""

from __future__ import annotations

import ast
import re
from dataclasses import dataclass
from pathlib import Path

import pytest


def _repo_root() -> Path:
    for parent in Path(__file__).resolve().parents:
        if (parent / "update_db" / "scripts" / "utils_update.py").is_file():
            return parent
    raise AssertionError(f"no checkout containing update_db/scripts found above {__file__}")


REPO_ROOT = _repo_root()
UPDATE_DB = REPO_ROOT / "update_db"

# This guard lives INSIDE the tree it audits, so every walk below must exclude its own
# directory. That is not a tidiness point. This file's prose names `sh` repeatedly, and while it
# was being prototyped the textual cross-check flagged four lines of its own docstrings as
# unexplained shell-outs -- it PASSED while the file sat under `docs/` and FAILED once it sat
# where it ships. The audit's verdict depended on where the audit was stored. Any `rglob` guard
# in this repository has to be run from where it will actually live before it is trusted.
EXCLUDED_DIRS = {"tests", "__pycache__", ".pytest_cache"}


def _sources(suffix: str) -> list[Path]:
    """Files under `update_db/` with this suffix, skipping the guard's own directory."""
    return sorted(
        p
        for p in UPDATE_DB.rglob(f"*{suffix}")
        if not EXCLUDED_DIRS & set(p.relative_to(UPDATE_DB).parts)
    )


# ---------------------------------------------------------------------------------------------
# Half one: what do the shell scripts require?
# ---------------------------------------------------------------------------------------------

BASH_SHEBANGS = ("#!/bin/bash", "#!/usr/bin/env bash")


def shell_scripts() -> dict[str, str]:
    """Relative path -> first line, for every `.sh` under `update_db/`."""
    out = {}
    for script in _sources(".sh"):
        first = script.read_text(errors="replace").splitlines()[0].strip()
        out[str(script.relative_to(REPO_ROOT))] = first
    return out


def requires_bash(shebang: str) -> bool:
    return any(shebang.startswith(s) for s in BASH_SHEBANGS)


# ---------------------------------------------------------------------------------------------
# Half two: which call sites pick the interpreter, and which one do they pick?
# ---------------------------------------------------------------------------------------------

POSIX_INTERPRETERS = {"sh", "/bin/sh"}


@dataclass(frozen=True)
class ShellOut:
    file: str
    line: int
    interpreter: str
    target: str  # the script argument as written in the source, best-effort
    form: str  # "argv" (a list) or "shell-string" (shell=True with a command string)

    def __str__(self) -> str:
        return f"{self.file}:{self.line}  [{self.form:<13}] {self.interpreter} {self.target}"


def _literal_head(node: ast.AST):
    """argv[0] of a `subprocess.run([...])` call, when it is a plain string literal."""
    if isinstance(node, (ast.List, ast.Tuple)) and node.elts:
        first = node.elts[0]
        if isinstance(first, ast.Constant) and isinstance(first.value, str):
            return first.value
    return None


def _render(node: ast.AST) -> str:
    """The second argv element as written -- a name, an f-string, or a literal."""
    try:
        return ast.unparse(node)
    except Exception:  # pragma: no cover - ast.unparse handles everything we produce
        return "<unrenderable>"


def _command_string_head(node: ast.AST):
    """`(interpreter, rest)` for a `shell=True` command written as a string or f-string.

    This second form is the whole reason this walk is not a one-liner. Four of the seven `sh`
    invocations in `utils_update.py` are NOT argv lists -- they are
    `subprocess.run(f"sh {script} ...", shell=True)`. A check that understood only argv lists
    found three of seven, and would have gone green the moment those three were fixed with four
    still broken. (Worse: under `shell=True` the command string is already handed to `/bin/sh`,
    so `f"sh {script}"` is dash launching dash.)
    """
    text = None
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        text = node.value
    elif isinstance(node, ast.JoinedStr) and node.values:
        first = node.values[0]
        if isinstance(first, ast.Constant) and isinstance(first.value, str):
            text = first.value
    if not text:
        return None
    parts = text.split(None, 1)
    if not parts:
        return None
    return parts[0], (parts[1] if len(parts) > 1 else "")


def shell_outs(paths: list[Path] | None = None) -> list[ShellOut]:
    """Every shell-out under `update_db/`, in both call forms `utils_update.py` uses."""
    found: list[ShellOut] = []
    for py in paths or _sources(".py"):
        tree = ast.parse(py.read_text(), filename=str(py))
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call) or not node.args:
                continue
            func = node.func
            name = (
                func.attr
                if isinstance(func, ast.Attribute)
                else func.id if isinstance(func, ast.Name) else None
            )
            if name != "run":
                continue
            rel = str(py.relative_to(REPO_ROOT))

            head = _literal_head(node.args[0])
            if head is not None:
                target = _render(node.args[0].elts[1]) if len(node.args[0].elts) > 1 else ""
                found.append(ShellOut(rel, node.lineno, head, target, "argv"))
                continue

            command = _command_string_head(node.args[0])
            if command is not None:
                interpreter, rest = command
                found.append(
                    ShellOut(rel, node.lineno, interpreter, rest.strip(), "shell-string")
                )
    return found


def posix_invocations() -> list[ShellOut]:
    return [s for s in shell_outs() if s.interpreter in POSIX_INTERPRETERS]


# ---------------------------------------------------------------------------------------------
# The properties.
# ---------------------------------------------------------------------------------------------


def test_every_update_db_shell_script_declares_bash():
    """The premise of the guard below. If a script is ever written as real POSIX `sh`, this
    fails and the next test's blanket rule has to become per-script -- which is the honest
    outcome, rather than the rule quietly over-reaching."""
    scripts = shell_scripts()
    assert scripts, "found no .sh files under update_db/, so this guard is checking nothing"
    posix = {p: s for p, s in scripts.items() if not requires_bash(s)}
    assert not posix, (
        "these scripts do not declare bash, so the interpreter rule below needs revisiting: "
        f"{posix}"
    )


def test_no_update_db_script_is_invoked_through_sh():
    """No call site hands a bash script to a POSIX shell. This now passes outright.

    It shipped from #360 as a **strict `xfail`** naming seven such call sites, deliberately ahead
    of the fix: swapping the interpreter on its own was measured to make the ClinVar failure
    *subtler* rather than absent -- it turns "ClinVar deleted" into "ClinVar overwritten with a
    0-byte VCF and a blank-version config", which Funcotator accepts as a valid data source with
    no records. So each site had to change together with its route's validation, and neither route
    could unmark this on its own. #353 took the four on the Funcotator route and #354 the three on
    the Annovar route; with both landed the test XPASSed, which -- the marker being strict -- turned
    this file red on purpose. That was the signal to remove the marker, and #354 removed it.

    Read the git history of this decorator before adding a `# noqa`-shaped exemption to anything
    below: the marker was load-bearing while it lasted, and its removal is the record that the
    defect is closed rather than tolerated.
    """
    violations = posix_invocations()
    assert not violations, (
        f"{len(violations)} call site(s) run a bash script through a POSIX shell:\n"
        + "\n".join(f"  {v}" for v in violations)
    )


def sh_mentions_by_text() -> dict[str, set]:
    """A deliberately crude textual scan for something that looks like invoking `sh`.

    Its job is to over-report: anything the AST walk fails to explain shows up here.
    """
    pattern = re.compile(r"""(^|[\s"'(\[])(/bin/)?sh(["'\s]|$)""")
    out: dict[str, set] = {}
    for py in _sources(".py"):
        rel = str(py.relative_to(REPO_ROOT))
        for number, line in enumerate(py.read_text().splitlines(), start=1):
            code = line.split("#", 1)[0]
            if pattern.search(code):
                out.setdefault(rel, set()).add(number)
    return out


def test_the_interpreter_audit_sees_every_sh_a_crude_scan_can_find():
    """The anti-vacuity control, and the one that caught this guard's first draft.

    A harness reporting zero findings is indistinguishable from a broken harness -- and one
    reporting *some* findings is indistinguishable from one that under-counts. The first version
    of this file walked only argv lists and found 3 of the 7 `sh` invocations in
    `utils_update.py`; the other 4 are `shell=True` command strings. It would have gone green
    with four call sites still broken.

    So: cross-check the AST walk against a crude regex scan. Every line the text scan flags must
    be a line the AST walk also explains. **A cross-check, never a count** -- a count assertion
    is the thing that has gone stale elsewhere in this repository, and it gets looser as the code
    changes, where a cross-check gets tighter.
    """
    ast_lines: dict[str, set] = {}
    for site in shell_outs():
        ast_lines.setdefault(site.file, set()).add(site.line)

    unexplained: dict[str, list] = {}
    for file, lines in sh_mentions_by_text().items():
        # subprocess.run() calls are frequently split across lines, so accept a match within a
        # few lines of a call the AST found.
        known = ast_lines.get(file, set())
        missed = [n for n in sorted(lines) if not any(abs(n - k) <= 3 for k in known)]
        if missed:
            unexplained[file] = missed

    assert not unexplained, (
        "a crude textual scan found `sh` on lines the AST walk does not explain, so this "
        f"audit is under-counting and could go green while call sites stay broken: {unexplained}"
    )


def test_the_interpreter_audit_finds_the_known_call_sites():
    """Belt and braces on the scan above: assert the walk still sees shell-outs at all, in both
    call forms. A refactor that hid every call site would otherwise satisfy the cross-check
    trivially -- nothing for the text scan to flag, nothing for the AST to explain.

    This one keeps working after the interpreters are fixed, because it looks at the call
    *forms* rather than at which interpreter they name.
    """
    # Deliberately with no explicit path, so this also checks that EXCLUDED_DIRS has not grown
    # to exclude the module the whole audit exists for.
    found = [s for s in shell_outs() if s.file.endswith("scripts/utils_update.py")]
    forms = {site.form for site in found}
    assert forms == {"argv", "shell-string"}, (
        f"expected both call forms in utils_update.py, found {forms or 'none'} -- the default "
        "walk has stopped seeing one of them"
    )


if __name__ == "__main__":
    print("shell scripts under update_db/ and what they declare")
    for path, shebang in shell_scripts().items():
        mark = "bash" if requires_bash(shebang) else "POSIX?"
        print(f"  [{mark:>6}] {path}")
        if not requires_bash(shebang):
            print(f"           first line: {shebang!r}")

    print("\nsubprocess.run call sites under update_db/")
    for site in shell_outs():
        flag = "  <-- POSIX shell" if site.interpreter in POSIX_INTERPRETERS else ""
        print(f"  {site}{flag}")

    bad = posix_invocations()
    print(f"\n{len(bad)} of them hand a bash script to a POSIX shell.")
