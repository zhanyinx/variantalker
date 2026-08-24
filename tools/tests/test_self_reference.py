"""Every self-reference this tree ships resolves, kept true by four rules about shape.

Issue #383, on map #376. The map was chartered on one defect —
``https://github.com/<project>/blob/master/CITATIONS.md``, printed by every run, naming
neither this repository nor any other — and the defect turned out to be one instance of a
class that nothing swept: 47 shipped prose references to a path the public export strips, and
175 bare ``#NNN`` cross-references, 55 of which GitHub silently auto-links to unrelated pull
requests in the public repository. Issue #396 and issue #397 fixed both classes. This module
is what stops them coming back, and issue #382 is the ruling it implements.

**All four rules are keyed on the *shape* of a token, never on what the token names.** That is
the property issue #234 could not get, and the reason there is no exemption list anywhere
below: a rule about form never has to decide whether a link to the samtools repository is a
legitimate third-party reference, because such a link simply does not match it.

1. **Well-formedness.** A ``github.com`` URL carrying ``/blob/`` or ``/tree/`` has two path
   segments before it. One segment is malformed *whoever* it names, so there is nothing to
   exempt. This is the chartering defect, and the one shape
   ``streamlit_app/tests/test_public_repo_name.py`` is structurally blind to: its own
   docstring calls a URL naming neither repository "permanently invisible".

2. **The branch-and-anchor policy** from issue #380, as **two independent halves**, because
   neither subsumes the other and a guard shipping one ships the other's defects for ever.
   (a) a URL naming the public repository carries no ref, unless it is a ``/blob/<ref>`` or
   ``/tree/<ref>`` deep path — the declared exception, which cannot be expressed without one —
   where the ref is the default branch. (b) every anchor on a link into a document in this
   tree resolves against that document's headings, by GitHub's slug rules. Issue #380 measured
   both directions: a ``tree/dev#output`` link *passes* an anchor check, because ``#output``
   genuinely resolves here, and is caught only by (a); the schema's ``#input`` named no ref at
   all and is caught only by (b).

3. **No shipped prose names a path the public export strips.** A backticked token
   **containing a ``/``** may not match a deny-list rule. The ``/`` is load-bearing and it is
   why this rule needs no exemptions: the deny-list names itself, so without that clause the
   rule fires on the thirteen sites that name the deny-list as a *mechanism* rather than as a
   path.

4. **No bare ``#NNN`` in a shipped Markdown document.** The number stays and only the sigil
   goes — ``divergence 6``, ``issue 199`` — because the sigil was doing double duty across two
   independent numeric indices, and dropping only the character would have written a tracker
   reference and an in-document row index identically.

**Where this runs, which was half the ticket.** ``.github/workflows/version-contract.yml``
runs ``pytest tools/tests`` — a *directory*, from the repository root — on every push and
every pull request. Nothing in this repository runs the full test suite; all the other
workflows enumerate their test files, so a guard in a file no workflow names is green by never
executing, which is how issue #355 merged red and sat there. Adding a module to a directory
some workflow globs is therefore the whole of the wiring, and it is the mechanism map #376 and
map #371 agreed to share rather than cut a workflow each. That workflow's job is labelled for
the version contract because issue #375 cut it first; the same is already true of
``update-db-contract.yml``, whose single-named job runs six unrelated contracts out of one
directory.

**No literal below violates a rule it is testing.** Every planted defect is composed at
import time or at call time from :data:`PUBLIC_OWNER`, :data:`PROJECT` and the deny-list read
off disk — never spelled out. That is not fastidiousness: the sweep reads every tracked file
and this module is a tracked file, so a spelled-out fixture would be a real violation and
:func:`test_the_sweep_reads_this_module` asserts the module stays inside its own remit. Issue
#375's guard learned this the hard way, passing on spelled-out fixtures only while its file
was still untracked.

**What this deliberately does not catch** — recorded in issue #382 as residue accepted rather
than overlooked, so a later reader does not mistake silence for coverage: the *semantic* class
(a framework wordmark, advice for profiles that exist elsewhere), staleness (a link that
resolves to a page frozen two years ago), everything needing the network (organisation
existence, tags, and the 34 ``docker://`` sites where a wrong name fails a run outright rather
than a click), a prose path with no defined base, and ``#NNN`` outside Markdown, which GitHub
does not auto-link and which therefore misleads nobody.
"""

from __future__ import annotations

import re
import subprocess
from pathlib import Path

import pytest

TESTS_DIR = Path(__file__).resolve().parent
REPO_ROOT = TESTS_DIR.parents[1]

#: This module's own path, as the sweep sees it. Used by
#: :func:`test_the_sweep_reads_this_module`: a guard module is the likeliest place in the tree
#: for the material it forbids to live, so exempting it would put the hole exactly where the
#: fixtures are.
THIS_MODULE = "tools/tests/test_self_reference.py"

#: The public repository, in halves. The host and the two path segments are never written as
#: one string anywhere in this module, so no fixture below can be a live URL.
HOST = "github.com"
RAW_HOST = "raw.githubusercontent.com"
PUBLIC_OWNER = "zhanyinx"
PROJECT = "variantalker"

#: The one ref a shipped deep path may name (issue #380). Not a variable in any meaningful
#: sense — the policy names the default branch specifically, and pinning to a tag or a commit
#: was ruled out on measurement rather than taste: the public repository's existing tags are
#: scheduled for deletion, a manifest-derived tag URL would 404 on every run made before
#: release day, and this repository's commits have no public existence at all, because the
#: export commits a fresh squashed sync rather than the history.
ALLOWED_REF = "main"

#: Path segments that introduce a ref. ``blob`` and ``tree`` are issue #380's declared
#: exception and are held to :data:`ALLOWED_REF`; the rest name a ref in a form the policy has
#: no exception for, so a shipped self-reference may not use them at all. Kept as one table so
#: the two halves of Rule 2(a) read as one policy.
REF_BEARING = ("blob", "tree", "commit", "commits", "raw", "archive")

#: The tracked files the sweep reads as text. A suffix list rather than a binary sniff, so the
#: population is stated rather than discovered — the same list issue #381's inventory swept
#: with, which is what makes this module's counts comparable to its published ones.
TEXT_SUFFIXES = frozenset(
    {
        ".md", ".py", ".sh", ".groovy", ".config", ".json", ".yml", ".yaml", ".Rmd", ".R",
        ".txt", ".cfg", ".toml", ".iss", ".nf", ".spec", ".plist", ".html", ".css", ".js",
        ".ini", ".in", ".bat", ".ps1", ".xml", ".csv", ".tsv", ".sql", ".pl", ".r",
    }
)
TEXT_NAMES = frozenset({"Makefile", "Dockerfile", "LICENSE", "CITATIONS.md", ".publicignore"})

#: Files the sweep must be reading, named rather than counted. A guard whose file list quietly
#: emptied — a broken ``git ls-files``, a deny-list rule that grew a leading ``*`` — would pass
#: over nothing at all and look identical to a clean tree. These four have each carried a
#: self-reference of a different one of the four rules' shapes. Paths, never line numbers:
#: issue #370's guard caught its own README claim going stale minutes after it was written.
SWEPT_FILES_CARRYING_A_SELF_REFERENCE = (
    "nextflow_schema.json",
    "bin/report.Rmd",
    "README.md",
    "streamlit_app/README.md",
)

URL_RE = re.compile(r"https?://[^\s)>\"'\]`,;]+")

#: A backtick span, single or doubled. Both spellings ship: Markdown uses one backtick, the
#: Python docstrings use RST's two.
SPAN_RE = re.compile(r"`+([^`\n]+?)`+")

#: Punctuation a path token picks up from the sentence around it.
STRIP = " \t,.;:)(\"'*<>[]!?"

MD_LINK_RE = re.compile(r"\[(?P<text>[^\]\n]*)\]\((?P<target>[^)\s]+)(?:\s+\"[^\"]*\")?\)")

#: Rule 4, with its three exclusions in the pattern itself.
#:
#: ``(?<![\w&])`` keeps a six-digit colour literal and an HTML entity out. The other two are
#: the exclusions issue #383 found by re-measuring the *merged* tree rather than a branch:
#: ``](`` and ``]: `` precede an in-document anchor link whose heading begins with a digit, so
#: its GitHub slug does too — ``[step 2](#2-download-annotation-databases)`` and its two
#: siblings in the root README. GitHub renders those as section links and never auto-links
#: them. Both exclusions are deliberately *narrow*: keying on a bare preceding ``(`` would
#: have been simpler and would have blinded the guard to ``(#328)``, one of the commonest
#: forms this class takes.
#:
#: Fenced blocks and inline code spans are excluded before the pattern is applied, by
#: :func:`_markdown_prose`, because a fence spans lines and a pattern cannot see across them.
SIGIL_RE = re.compile(r"(?<![\w&])(?<!\]\()(?<!\]: )#(\d{1,4})\b")


# ---------------------------------------------------------------------------
# The shipped tree
# ---------------------------------------------------------------------------


def deny_rules(root: Path) -> list[str] | None:
    """The public export's deny-list, read from ``.publicignore`` **directly**.

    Never by importing the export script, and this is house precedent rather than preference:
    ``streamlit_app/tests/test_makefile_targets_are_runnable.py`` states it in place, because
    that script is itself on the deny-list and "importing it here would make this module fail
    in the tree it most needs to pass in".

    Returns ``None`` where the deny-list is absent, "which is exactly the exported tree" — and
    there the honest move is to skip and say so rather than to pass over nothing, the shape
    issue #348 set. Rule 3 is not merely unchecked in a public clone, it is **unaskable**: the
    list that defines it is one of the things the list strips.
    """
    path = root / ".publicignore"
    if not path.exists():
        return None
    return [
        line.strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.strip().startswith("#")
    ]


def rule_matches(rule: str, relpath: str) -> bool:
    """The deny-list grammar, as the export applies it: a directory rule, or one exact file."""
    target = rule.lstrip("/")
    if rule.endswith("/"):
        return relpath.startswith(target)
    return relpath == target


def tracked_files(root: Path) -> list[str]:
    """Every file git tracks under *root*, as root-relative paths.

    Through git rather than a filesystem walk, for the reason the export's own test module
    gives: ``.claude/worktrees`` holds full checkouts of this same tree, so a walk answers
    questions about a branch this test is not running on, and it would additionally read
    untracked run artifacts — a local database download, a Nextflow work directory — which
    state nothing this repository ships.
    """
    completed = subprocess.run(
        ["git", "ls-files"], cwd=root, capture_output=True, text=True, check=True
    )
    return [line for line in completed.stdout.splitlines() if line]


def swept_files(root: Path) -> list[str]:
    """The shipped tree, reproduced in-process: tracked files minus the deny-list.

    No export rehearsal and no network. This is the same computation the export's own test
    module performs, and it agrees with issue #381's inventory of the rehearsed tree file for
    file — but that module is deny-listed, so a shipping guard cannot import it and has to
    re-derive it here.

    Raises rather than returning nothing, because a sweep over an empty tree passes every rule
    below and is indistinguishable from a clean tree. A shallow or non-git checkout is the
    realistic way to arrive here.
    """
    tracked = tracked_files(root)
    if not tracked:
        raise AssertionError(
            f"git ls-files listed nothing under {root} — this guard would pass over any tree "
            "at all, so the checkout is what fails here, not the rules"
        )
    rules = deny_rules(root)
    if rules is None:
        return tracked
    return [path for path in tracked if not any(rule_matches(r, path) for r in rules)]


def is_text(relpath: str) -> bool:
    path = Path(relpath)
    return path.name in TEXT_NAMES or path.suffix in TEXT_SUFFIXES


def read_sources(root: Path, paths: list[str]) -> dict[str, str]:
    sources: dict[str, str] = {}
    for rel in paths:
        if not is_text(rel):
            continue
        try:
            text = (root / rel).read_text(encoding="utf-8")
        except (UnicodeDecodeError, OSError):
            continue
        sources[rel] = text
    return sources


# ---------------------------------------------------------------------------
# Anchors, by GitHub's slug rules
# ---------------------------------------------------------------------------


def heading_slugs(text: str) -> set[str]:
    """The anchors GitHub generates for a Markdown document's headings.

    Slugs, never headings and never line numbers. Both anchors issue #393 landed sit at
    different line numbers in this tree and in the public one, and a heading's *case* has been
    transcribed wrongly on this map at least once — so a guard comparing heading strings goes
    stale where one comparing slugs does not.

    Lines inside a fence are not headings, and the schema's ``#input`` defect is invisible
    without that exclusion.
    """
    slugs: set[str] = set()
    fenced = False
    for line in text.splitlines():
        if line.strip().startswith("```") or line.strip().startswith("~~~"):
            fenced = not fenced
            continue
        if fenced:
            continue
        match = re.match(r"^(#{1,6})\s+(.*?)\s*#*\s*$", line)
        if not match:
            continue
        title = match.group(2)
        title = re.sub(r"`([^`]*)`", r"\1", title)
        title = re.sub(r"\[([^\]]*)\]\([^)]*\)", r"\1", title)
        title = re.sub(r"[*_]", "", title)
        slug = title.strip().lower()
        slug = re.sub(r":[a-z0-9_+-]+:", "", slug)
        slug = re.sub(r"[^\w\s-]", "", slug, flags=re.UNICODE)
        slug = re.sub(r"\s+", "-", slug.strip())
        if slug:
            slugs.add(slug)
    return slugs


def resolve(rel_from: str, target: str) -> str:
    """A relative link target, resolved against the file that carries it."""
    parts: list[str] = []
    for seg in (Path(rel_from).parent / target).as_posix().split("/"):
        if seg in ("", "."):
            continue
        if seg == "..":
            if parts:
                parts.pop()
            continue
        parts.append(seg)
    return "/".join(parts)


# ---------------------------------------------------------------------------
# The four rules, as pure functions over path-to-text
# ---------------------------------------------------------------------------
#
# Pure so that every one of them can be pointed at planted content and shown to fail. A rule
# that reads the real tree directly can only ever be shown green, and this repository has a
# documented history of guards that were green because they were vacuous.


def _url_path_segments(url: str) -> tuple[str, list[str]]:
    """Split a URL into its host and its path segments, query and fragment discarded."""
    without_scheme = re.sub(r"^https?://", "", url)
    host, _, path = without_scheme.partition("/")
    path = path.split("#", 1)[0].split("?", 1)[0]
    return host, [seg for seg in path.split("/") if seg]


def rule_1_malformed_deep_paths(sources: dict[str, str]) -> list[tuple[str, int, str]]:
    """A ``github.com`` deep path with fewer than two segments before ``/blob/`` or ``/tree/``."""
    found: list[tuple[str, int, str]] = []
    for rel, text in sources.items():
        for lineno, line in enumerate(text.splitlines(), 1):
            for match in URL_RE.finditer(line):
                url = match.group(0).rstrip(".,:;\"'")
                host, segments = _url_path_segments(url)
                if host != HOST:
                    continue
                for index, segment in enumerate(segments):
                    if segment in ("blob", "tree"):
                        if index < 2:
                            found.append((rel, lineno, url))
                        break
    return found


def rule_2a_refs(sources: dict[str, str]) -> list[tuple[str, int, str]]:
    """A self-reference naming a ref other than the default branch.

    Scoped to URLs that name the public repository *already*, which is why issue #234's
    objection does not apply here. Issue #234 declined the rule "every URL naming a project
    repository must name the public one" because it needed exemptions for around twenty-five
    genuine third-party links; that rule was keyed on URLs that *ought to* name the public
    repository. This one is keyed on URLs that already do, so every third-party link is out of
    scope by construction rather than by exemption.
    """
    found: list[tuple[str, int, str]] = []
    for rel, text in sources.items():
        for lineno, line in enumerate(text.splitlines(), 1):
            for match in URL_RE.finditer(line):
                url = match.group(0).rstrip(".,:;\"'")
                host, segments = _url_path_segments(url)
                if host not in (HOST, RAW_HOST):
                    continue
                if segments[:2] != [PUBLIC_OWNER, PROJECT]:
                    continue
                if host == RAW_HOST:
                    # The raw host names its ref positionally, with no ``blob``/``tree``
                    # segment to introduce it. The schema's ``$id`` is this shape and names
                    # the default branch, so it is the declared exception rather than a
                    # violation — worth stating, because it is also the only place in this
                    # tree that ever hardcoded a ref at all.
                    if len(segments) > 2 and segments[2] != ALLOWED_REF:
                        found.append((rel, lineno, url))
                    continue
                rest = segments[2:]
                if not rest or rest[0] not in REF_BEARING:
                    continue
                if rest[0] in ("blob", "tree"):
                    if len(rest) < 2 or rest[1] != ALLOWED_REF:
                        found.append((rel, lineno, url))
                else:
                    found.append((rel, lineno, url))
    return found


def rule_2b_anchors(
    sources: dict[str, str], tracked: set[str]
) -> list[tuple[str, int, str]]:
    """An anchor that does not resolve against the headings of the document it points into.

    Three anchor-bearing shapes, all resolved against **this** tree. That is a proxy for the
    public one, and an exact proxy at export time — the export ships this tree verbatim with
    no content rewriting, so it drifts only between syncs, and Rule 2(a) forbidding every ref
    but the default branch is what keeps the drift bounded to one branch.
    """
    slug_cache: dict[str, set[str]] = {}

    def slugs_of(rel: str) -> set[str] | None:
        if rel not in slug_cache:
            if rel in sources:
                slug_cache[rel] = heading_slugs(sources[rel])
            else:
                return None
        return slug_cache[rel]

    def check(rel: str, lineno: int, ref: str, target_doc: str, anchor: str,
              found: list[tuple[str, int, str]]) -> None:
        if not target_doc.endswith(".md") or target_doc not in tracked:
            return
        slugs = slugs_of(target_doc)
        if slugs is None:
            return
        if anchor.lower() not in slugs:
            found.append((rel, lineno, ref))

    found: list[tuple[str, int, str]] = []
    for rel, text in sources.items():
        for lineno, line in enumerate(text.splitlines(), 1):
            if rel.endswith(".md"):
                for match in MD_LINK_RE.finditer(line):
                    target = match.group("target")
                    if "{" in target:
                        # Three of the fifty Markdown-link matches in this tree are not links
                        # — two brace interpolations, and a subscript followed by a call that
                        # reads as a link and has no brace to skip on. Recorded in issue #382
                        # so this guard does not reintroduce the false positive by widening to
                        # Markdown links generally.
                        continue
                    if re.match(r"^(https?:|mailto:|ftp:)", target):
                        continue
                    path_part, _, anchor = target.partition("#")
                    if not anchor:
                        continue
                    target_doc = resolve(rel, path_part) if path_part else rel
                    check(rel, lineno, target, target_doc, anchor, found)
            for match in URL_RE.finditer(line):
                url = match.group(0).rstrip(".,:;\"'")
                if "#" not in url:
                    continue
                host, segments = _url_path_segments(url)
                if host != HOST or segments[:2] != [PUBLIC_OWNER, PROJECT]:
                    continue
                anchor = url.split("#", 1)[1]
                rest = segments[2:]
                if not rest:
                    target_doc = "README.md"
                elif rest[0] == "blob" and len(rest) > 2:
                    target_doc = "/".join(rest[2:])
                else:
                    continue
                check(rel, lineno, url, target_doc, anchor, found)
    return found


def rule_3_denied_paths(
    sources: dict[str, str], rules: list[str]
) -> list[tuple[str, int, str]]:
    """A backticked token containing a ``/`` that the public export would strip.

    Two token definitions, deliberately: the whole backtick span, and each whitespace-split
    token inside it. Issue #396 measured both and they agree exactly, which is how it settled
    a published disagreement between two counts — so a divergence here would mean one of the
    two readings has drifted, not that the tree has.

    The predicate is *deny-listed*, never *absent*. Twelve shipped path tokens name something
    that does not exist and are correct for saying so: gitignored build stamps, a fixture the
    prose says was deleted, a database directory the prose says a fresh clone does not have.
    Existence is not this guard's business.
    """
    found: list[tuple[str, int, str]] = []
    for rel, text in sources.items():
        for lineno, line in enumerate(text.splitlines(), 1):
            for match in SPAN_RE.finditer(line):
                span = match.group(1).strip()
                candidates = [span] + [t.strip(STRIP) for t in span.split()]
                for token in candidates:
                    if "/" in token and any(rule_matches(r, token) for r in rules):
                        found.append((rel, lineno, token))
                        break
    return found


def _markdown_prose(text: str) -> list[tuple[int, str]]:
    """A Markdown document's lines with its code blanked out, positions preserved.

    Fenced blocks drop out entirely and inline spans are replaced by spaces of the same width,
    so a column offset still means something and the two lookbehinds in :data:`SIGIL_RE` still
    see the characters that really precede a match.
    """
    out: list[tuple[int, str]] = []
    fenced = False
    for lineno, line in enumerate(text.splitlines(), 1):
        stripped = line.strip()
        if stripped.startswith("```") or stripped.startswith("~~~"):
            fenced = not fenced
            continue
        if fenced:
            continue
        out.append((lineno, SPAN_RE.sub(lambda m: " " * len(m.group(0)), line)))
    return out


def rule_4_bare_sigils(sources: dict[str, str]) -> list[tuple[str, int, str]]:
    """A bare ``#NNN`` in a shipped Markdown document.

    Markdown only, and that is a scope rather than an oversight: GitHub auto-links ``#NNN`` in
    rendered Markdown and not in a source blob, so the same token in a Python docstring
    misleads nobody. It is also why this class is a hundred-and-seventy-odd sites rather than
    many hundreds.
    """
    found: list[tuple[str, int, str]] = []
    for rel, text in sources.items():
        if not rel.endswith(".md"):
            continue
        for lineno, line in _markdown_prose(text):
            for match in SIGIL_RE.finditer(line):
                found.append((rel, lineno, match.group(0)))
    return found


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def swept() -> list[str]:
    return swept_files(REPO_ROOT)


@pytest.fixture(scope="module")
def sources(swept: list[str]) -> dict[str, str]:
    return read_sources(REPO_ROOT, swept)


@pytest.fixture(scope="module")
def rules() -> list[str]:
    found = deny_rules(REPO_ROOT)
    if found is None:
        pytest.skip(
            "no .publicignore in this tree, which is exactly the exported tree. Rule 3 is "
            "unaskable here rather than merely unchecked: the list that defines which paths "
            "are stripped is one of the things it strips. Skipping and saying so, rather than "
            "passing over nothing — see streamlit_app/tests/test_suite_organisation.py"
        )
    assert found, (
        ".publicignore parsed to no rules at all, so Rule 3 would read nothing and pass. A "
        "deny-list that parses to nothing is the vacuous shape, not a clean tree"
    )
    return found


# ---------------------------------------------------------------------------
# Vacuity: what the sweep read, before anything it concluded
# ---------------------------------------------------------------------------


def test_the_swept_tree_is_the_shipped_tree(swept: list[str]) -> None:
    """The sweep read a real tree, asserted by name and not by count alone.

    Counts on this map are method-dependent — two published pairs disagree honestly, 41 versus
    47 and 174 versus 175, neither of them drift — so a remembered total is the one thing this
    module must not assert. Files and rules instead.
    """
    assert len(swept) > 100, (
        f"the sweep read {len(swept)} files. Every rule below passes over a tree this small "
        "for the wrong reason"
    )
    missing = [p for p in SWEPT_FILES_CARRYING_A_SELF_REFERENCE if p not in swept]
    assert not missing, (
        f"{missing} are not in the swept tree, so the rules that were written for them read "
        "nothing. Either the deny-list grew a rule that strips them, or they were renamed"
    )
    stripped = [p for p in tracked_files(REPO_ROOT) if p not in swept]
    if deny_rules(REPO_ROOT) is None:
        # A public clone. The tracked tree *is* the shipped tree, so stripping nothing is the
        # correct answer rather than a broken sweep — asserting otherwise here is what made
        # this module fail its own export rehearsal, the very shape the deny-list's own house
        # rule warns about: a guard that fails in the tree it most needs to pass in.
        assert not stripped
    else:
        assert stripped, (
            "the deny-list parsed but stripped nothing at all, so this sweep is reading the "
            "private tree and calling it the shipped one"
        )


def test_the_sweep_reads_this_module(swept: list[str]) -> None:
    """This module is inside its own remit.

    A guard module is the likeliest place in the tree for the material it forbids to live,
    since its fixtures are made of that material. Exempting it would put the hole exactly
    where the hole would be worst — which is why every fixture below is composed rather than
    spelled. Issue #375's guard found this by passing on spelled-out fixtures while its file
    was still untracked, and going red on itself the moment it was committed.
    """
    assert THIS_MODULE in swept, (
        f"{THIS_MODULE} is not swept, so the fixtures below are unchecked. If this module was "
        "renamed, rename the constant with it"
    )


def test_an_empty_tree_fails_rather_than_passing(tmp_path: Path) -> None:
    """A sweep that reads nothing is a failure, not a clean tree.

    The realistic way to arrive here is a shallow clone or a checkout git does not track, and
    the two are indistinguishable from a tree with no defects once the file list is empty.
    """
    subprocess.run(["git", "init", "--quiet"], cwd=tmp_path, check=True)
    with pytest.raises(AssertionError, match="listed nothing"):
        swept_files(tmp_path)


def test_a_missing_deny_list_skips_rather_than_passing(tmp_path: Path) -> None:
    """Rule 3 in a public clone: unaskable, so it says so.

    The other three rules still run there — they need no list to state their predicate, and in
    a public clone the tracked tree *is* the shipped tree.
    """
    subprocess.run(["git", "init", "--quiet"], cwd=tmp_path, check=True)
    (tmp_path / "README.md").write_text("nothing here\n", encoding="utf-8")
    subprocess.run(["git", "add", "README.md"], cwd=tmp_path, check=True)
    assert deny_rules(tmp_path) is None
    assert swept_files(tmp_path) == ["README.md"]


# ---------------------------------------------------------------------------
# The four rules, against the tree
# ---------------------------------------------------------------------------


def _report(sites: list[tuple[str, int, str]]) -> str:
    return "\n".join(f"  {rel}:{lineno} → {token}" for rel, lineno, token in sites)


def test_rule_1_every_deep_path_names_a_repository(sources: dict[str, str]) -> None:
    sites = rule_1_malformed_deep_paths(sources)
    assert not sites, (
        "these deep-path URLs have fewer than two path segments before /blob/ or /tree/, so "
        "they name no repository in any organisation and a reader following one gets a 404:\n"
        f"{_report(sites)}"
    )


def test_rule_2a_every_self_reference_names_the_default_branch(
    sources: dict[str, str],
) -> None:
    sites = rule_2a_refs(sources)
    assert not sites, (
        "these self-references name a ref other than the default branch. A shipped link may "
        f"name no ref at all, or {ALLOWED_REF!r} as a deep path; a branch other than the "
        "default one freezes what a reader sees at whenever that branch last moved, and a tag "
        "or a commit has no stable public existence here at all:\n"
        f"{_report(sites)}"
    )


def test_rule_2b_every_anchor_resolves(sources: dict[str, str], swept: list[str]) -> None:
    sites = rule_2b_anchors(sources, set(swept))
    assert not sites, (
        "these anchors do not match any heading slug in the document they point into, so "
        "GitHub drops the reader at the top of the page with no indication anything was "
        f"missed:\n{_report(sites)}"
    )


def test_rule_3_no_shipped_prose_names_a_stripped_path(
    sources: dict[str, str], rules: list[str]
) -> None:
    sites = rule_3_denied_paths(sources, rules)
    assert not sites, (
        "these backticked tokens name a path the public export strips, so a public reader is "
        "sent to a file that is not in their clone. Where the path is a code literal, point "
        "the prose at the symbol instead of the value — a public reader has the symbol:\n"
        f"{_report(sites)}"
    )


def test_rule_4_no_bare_sigil_in_shipped_markdown(sources: dict[str, str]) -> None:
    sites = rule_4_bare_sigils(sources)
    assert not sites, (
        "these bare sigil-plus-number tokens ship in Markdown, where GitHub auto-links any "
        "number that exists in the public repository — silently, to an unrelated pull request, "
        "and newly so every time one is merged. Keep the number and drop the sigil, naming "
        f"what the number indexes if the sentence does not:\n{_report(sites)}"
    )


# ---------------------------------------------------------------------------
# Every rule, shown red on a planted defect
# ---------------------------------------------------------------------------
#
# A guard that has never been seen failing is not evidence. Two of these four rules have zero
# live violations in this tree and always will if they work, so a planted defect is their only
# possible proof — and the table is what makes adding a fifth rule acquire its proof rather
# than merely its assertion.


def _planted(rules: list[str] | None) -> list[tuple[str, str, dict[str, str], int]]:
    """One planted defect per rule, composed rather than spelled.

    Returns ``(rule name, description, sources, expected site count)``. Every URL is built
    from the host and repository halves at call time, so nothing here is a live self-reference
    in this module's own source — which the sweep reads.

    *rules* is ``None`` in a public clone, where Rule 3 is unaskable. Only Rule 3's own entry
    drops out there: the other three rules need no deny-list to state their predicate, so
    their proofs travel with them rather than skipping alongside a rule they do not depend on.
    """
    doc = "# Title\n\n## Input Format\n\n"

    rule_3_entries: list[tuple[str, str, dict[str, str], int]] = []
    if rules is not None:
        denied_dir = next((r for r in rules if r.endswith("/")), None)
        assert denied_dir, (
            "the deny-list has no directory rule, so Rule 3's planted defect cannot be "
            "composed from it and would have to be spelled out in this file"
        )
        denied_path = denied_dir.lstrip("/") + "issue-1/notes.md"
        rule_3_entries = [
            (
                "rule 3",
                "backticked prose naming a path the public export strips",
                {"README.md": f"The notes are in `{denied_path}` for now.\n"},
                1,
            ),
        ]

    return [
        (
            "rule 1",
            "a deep path with one segment before /blob/ — the chartering defect, naming "
            "neither this repository nor any other",
            {"CITATIONS.md": f"See https://{HOST}/{PROJECT}/blob/master/CITATIONS.md\n"},
            1,
        ),
        (
            "rule 2a",
            "a self-reference naming a branch that is not the default one",
            {"docs.md": f"See https://{HOST}/{PUBLIC_OWNER}/{PROJECT}/tree/dev/README.md\n"},
            1,
        ),
        (
            "rule 2a (tag)",
            "a self-reference naming a tag, a ref shape the policy has no exception for",
            {"docs.md": f"https://{HOST}/{PUBLIC_OWNER}/{PROJECT}/archive/v1.tar.gz\n"},
            1,
        ),
        (
            "rule 2b",
            "an anchor naming a section the target document does not have",
            {
                "README.md": doc,
                "schema.json": f'"help": "https://{HOST}/{PUBLIC_OWNER}/{PROJECT}#input"\n',
            },
            1,
        ),
        (
            "rule 2b (relative)",
            "a relative Markdown link whose anchor names no heading in its target",
            {"README.md": doc, "docs/guide.md": "see [the format](../README.md#inpt)\n"},
            1,
        ),
        *rule_3_entries,
        (
            "rule 4",
            "a bare sigil-plus-number that GitHub would auto-link",
            {"README.md": "Fixed by #999, and silenced by #7 in the table above.\n"},
            2,
        ),
    ]


def _run_rule(
    name: str, sources: dict[str, str], rules: list[str] | None
) -> list[tuple[str, int, str]]:
    tracked = set(sources)
    if name.startswith("rule 1"):
        return rule_1_malformed_deep_paths(sources)
    if name.startswith("rule 2a"):
        return rule_2a_refs(sources)
    if name.startswith("rule 2b"):
        return rule_2b_anchors(sources, tracked)
    if name.startswith("rule 3"):
        assert rules is not None, "rule 3's runner reached a tree with no deny-list"
        return rule_3_denied_paths(sources, rules)
    if name.startswith("rule 4"):
        return rule_4_bare_sigils(sources)
    raise AssertionError(f"no runner for {name!r}, so its planted defect proves nothing")


def test_every_rule_is_red_on_a_planted_defect() -> None:
    """Each rule, run against content built out of the defect it exists for.

    Reads the deny-list directly rather than through the ``rules`` fixture, so that the three
    rules needing no deny-list keep their proof in a public clone instead of skipping
    alongside the one rule that does.
    """
    rules = deny_rules(REPO_ROOT)
    for name, description, planted, expected in _planted(rules):
        sites = _run_rule(name, planted, rules)
        assert len(sites) == expected, (
            f"{name} found {len(sites)} sites in a tree planted with {description}, expected "
            f"{expected}. A rule that cannot see its own defect is the vacuous shape:\n"
            f"{_report(sites)}"
        )


def test_every_rule_is_green_once_the_planted_defect_is_removed() -> None:
    """The other half of the proof: red *because of* the defect, not red regardless.

    A rule matching everything would pass the test above and be useless.
    """
    rules = deny_rules(REPO_ROOT)
    for name, _description, planted, _expected in _planted(rules):
        cleaned = {rel: "" for rel in planted}
        sites = _run_rule(name, cleaned, rules)
        assert not sites, (
            f"{name} still reports sites after its planted content was emptied, so the test "
            f"above proved nothing about the defect:\n{_report(sites)}"
        )


# ---------------------------------------------------------------------------
# Negative controls: the shapes each rule must not fire on
# ---------------------------------------------------------------------------
#
# Every one of these is a real shape this tree ships. A guard that fires on them would be
# deleted, which is the failure mode behind half the residue issue #382 accepted.


def test_rule_1_ignores_a_third_party_deep_path() -> None:
    """The property issue #234 could not get: no exemption, because there is no match.

    A third-party deep path is well-formed, so it does not match a well-formedness rule at
    all, and the rule never has to decide whether the link is legitimate.
    """
    sources = {
        "README.md": f"See https://{HOST}/samtools/samtools/blob/develop/NEWS\n"
        f"and https://{HOST}/broadinstitute/gatk/tree/master/scripts\n"
    }
    assert not rule_1_malformed_deep_paths(sources)
    assert not rule_2a_refs(sources), (
        "a third-party repository's branch is not this project's ref to choose. Rule 2(a) is "
        "scoped to URLs that already name the public repository"
    )


def test_rule_2a_accepts_the_two_compliant_shapes() -> None:
    """A bare repository URL, and the declared deep-path exception."""
    sources = {
        "README.md": f"https://{HOST}/{PUBLIC_OWNER}/{PROJECT}\n"
        f"https://{HOST}/{PUBLIC_OWNER}/{PROJECT}/blob/{ALLOWED_REF}/CITATIONS.md\n"
        f"https://{RAW_HOST}/{PUBLIC_OWNER}/{PROJECT}/{ALLOWED_REF}/nextflow_schema.json\n"
    }
    assert not rule_2a_refs(sources)


def test_rule_3_ignores_the_deny_list_named_as_a_mechanism(rules: list[str]) -> None:
    """The ``/`` clause, which is why Rule 3 needs no exemption list.

    Thirteen shipped sites name the deny-list as a *mechanism* rather than as a path. A rule
    without this clause fires on all thirteen and needs them excused by name — the list this
    repository refuses. With it they drop out mechanically, which is a property of the shape.
    """
    bare = next((r.lstrip("/") for r in rules if "/" not in r.lstrip("/")), None)
    assert bare, "the deny-list has no root-level file rule, so this control tests nothing"
    sources = {"Makefile": f"# every rule in `{bare}` still matches something tracked\n"}
    assert not rule_3_denied_paths(sources, rules), (
        f"Rule 3 fired on `{bare}` named as a mechanism. The token contains no '/', so the "
        "clause that excludes it is missing"
    )


def test_rule_4_ignores_its_three_exclusions() -> None:
    """Fences, inline spans, and an in-document anchor whose heading starts with a digit.

    The third was found by re-measuring the *merged* tree rather than a branch: it came into
    existence only when a sibling pull request landed three such links in the root README, and
    the branch that measured a clean zero was right about its own base and wrong about the
    tree. Digit-numbered headings are a standing convention there, so this shape recurs — it
    is a property of the rule, not an exemption for three sites.
    """
    sources = {
        "README.md": (
            "See [step 2](#2-download-annotation-databases) and\n"
            "[step 3](#3-run-the-pipeline) for the order.\n"
            "\n"
            "[ref]: #4-outputs\n"
            "\n"
            "A literal `#999` inside a span stays, and so does a whole block:\n"
            "\n"
            "```\n"
            "curl -H 'X-Thing: #123' https://example.invalid\n"
            "```\n"
            "\n"
            "A colour literal like #1a2b3c and an entity like &#39; are not this shape.\n"
        )
    }
    assert not rule_4_bare_sigils(sources)


def test_rule_4_still_catches_the_narrow_form_a_wider_exclusion_would_lose() -> None:
    """``(#328)`` is a real site, so exclusion three keys on ``](`` and not on a bare ``(``.

    Stated as a test rather than a comment because it is the one way to get that exclusion
    wrong, and the wrong version is the simpler one.
    """
    sites = rule_4_bare_sigils({"README.md": "Reverted the browser (#328) after review.\n"})
    assert len(sites) == 1, (
        "a parenthesised cross-reference is one of the commonest forms this class takes here. "
        "An exclusion keyed on a bare '(' rather than on '](' loses it"
    )


def test_rule_4_reads_markdown_only() -> None:
    """The scope, asserted so it is not mistaken for coverage.

    GitHub auto-links in rendered Markdown, not in a source blob, so the same token in a
    Python docstring misleads nobody and stays. That is what keeps the rule aimed at the
    surface where the defect is real.
    """
    line = "Fixed by #999 in the previous release.\n"
    assert rule_4_bare_sigils({"README.md": line})
    assert not rule_4_bare_sigils({"tools/thing.py": f'"""{line}"""\n'})


def test_the_anchor_slugger_follows_github_rules() -> None:
    """Slugs, against the heading shapes this tree actually ships.

    Held here rather than trusted, because Rule 2(b) is only as good as this function and the
    two anchors issue #393 landed both depend on punctuation and case being handled the way
    GitHub handles them.
    """
    slugs = heading_slugs(
        "## Input Format\n"
        "### Default Filters\n"
        "## 2. Download annotation databases\n"
        "## :link: Map: every self-reference\n"
        "## `nextflow.config`, and what reads it\n"
        "```\n"
        "## Not A Heading\n"
        "```\n"
    )
    assert "input-format" in slugs
    assert "default-filters" in slugs
    assert "2-download-annotation-databases" in slugs
    assert "map-every-self-reference" in slugs
    assert "nextflowconfig-and-what-reads-it" in slugs
    assert "not-a-heading" not in slugs, (
        "a heading-shaped line inside a fence is not a heading, and the schema's broken anchor "
        "is invisible without that exclusion"
    )
