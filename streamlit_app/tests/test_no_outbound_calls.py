"""MAFigate makes no outbound network call, asserted rather than remembered.

Issue #257.

**What this guard covers, exactly.** The *shipped application*, in all three languages it is
written in:

* every tracked Python file under ``streamlit_app/`` except ``tests/`` — ``build/`` included,
  for the reason under :data:`UNSHIPPED_DIRS`;
* every tracked ``.sh``, ``.bat`` and ``Makefile`` under it bar ``tests/`` and the two named
  build scripts, which is the two installer launchers plus the clone channel's route in —
  derived rather than listed, so a new launcher is covered without editing this file;
* ``MAFigateLauncher.swift``, compiled into the ``.app`` and the executable the user's
  double-click actually runs. It was unscanned in this module's first draft, which made the
  claim true of the app's Python and shell and silent about its entry point.

Nothing else. The scope is stated here because a privacy claim asserted over the wrong set of
files is worse than no assertion at all: it reads as a proof and is a coincidence.

**What it deliberately does not cover, and why that is not a loophole.** Two things in this
repository legitimately reach the network, and neither ships with MAFigate:

* ``update_db/`` fetches annotation databases over FTP and HTTP. It is the *pipeline's*
  tooling; the app never imports it and neither installer copies it. It is also this
  module's live demonstration that the scanner works — see
  :func:`test_the_scanner_finds_the_real_network_client_outside_the_app`.
* ``tests/`` — this suite. A test that needed a socket would be a legitimate test, and if
  ``tests/`` were in scope the property could be satisfied by *moving* a phone-home into a
  test file, which ships nowhere and proves nothing.

The build scripts download a portable Python, and they are out of scope by *language* rather
than by exemption: ``build_dmg.sh`` and ``build_installer.bat`` are shell, run on a
maintainer's machine, and are stripped from both bundles. Nothing in this file exempts them,
so a Python build script arriving there meets the guard, which is the right conversation to
have rather than a hole.

Excluding ``tests/`` and ``update_db/`` narrows the claim honestly. What would make the guard
worthless is the inverse — letting one of those files *answer* for the app — so the
exclusions are named paths here rather than a pattern that could quietly grow.

**The one thing the shipped launchers do reach the network for**, stated out loud because a
guard that hid it would be lying: on first launch each one runs ``pip install`` against
``requirements.txt`` to build the virtual environment. That is the install step, not the
application, and it is documented as such in ``build/BUILD_INSTRUCTIONS.md`` and the
README. :data:`PIP_INSTALL_MARKERS` allows exactly that shape — an install of the pinned
requirements, or pip upgrading itself — and refuses any other package fetch, so an update
check dressed as ``pip install --upgrade mafigate`` is caught rather than waved through.

**Why the property is worth a guard at all — issue #229.** The delivery matrix decided
there will be **no update check, ever**. Not "not yet": ever. The argument was that MAFigate
makes zero outbound calls today, so a version check would be its *first* one, and would
falsify the data-posture claim the README makes to every channel on the very day the
installers shipped. That decision is the reason this file exists. A maintainer who arrives
here wanting to add one is not being told "a test says no" — they are being told which
decision they are reopening and what it costs, and the honest way to reopen it is to change
#229's answer and the README's promise first, in the open, and then this guard.

**Why an import scan, and where it is blind.** The claim is about clients, not URLs. The app
renders ClinVar and gnomAD links into the variant panel on purpose — those are hrefs the
*browser* follows if the user clicks them, which the README already distinguishes from
lookups the app performs — so a guard that refused remote URL literals would refuse the
app's own correct behaviour. What it can refuse is the machinery for following one: an
imported client, a dynamically named one, a remote URL handed to a pandas reader, and a
network command line in a launcher.

An import scan alone is not enough, and three rules exist because of shapes that need no
denied import at all. :data:`NETWORK_CALL_NAMES` catches the call however the module reaching
it was spelled — ``asyncio``'s ``open_connection``, or ``import urllib`` followed by
``urllib.request.urlopen``, which works whenever a dependency has already imported the
submodule. :data:`SHELL_OUT_CALLS` catches handing the job to a program:
``subprocess.run(["curl", url])`` imports nothing this file denies. And the shell scan reads
per *command*, not per line, because a loopback address anywhere on a line was found vouching
for every other command on it.

**Blind spots, recorded rather than papered over.** Third-party dependencies are out of reach:
Streamlit itself reports usage statistics by default, which is why the README carries a caveat
and ``test_delivery_channels_copy.py`` holds it there. A URL held in a *variable* is invisible
to the Python half — ``pd.read_csv(UPDATE_URL)`` reads as a local path — and so is a module
named by a variable rather than a literal. The shell half closes its own version of that gap by
refusing a fetch whose target it cannot see at all, which the Python half cannot do without
refusing every local read in the app.

**Local ports are not outbound.** The deleted launcher probed for a free port by binding one
on ``127.0.0.1``; ``launch.sh`` does the same with ``lsof`` and ``launch.bat`` with
``netstat``, and both still do. Nothing here objects to that, and #257 removed an unshipped
file rather than a capability: what it bought is that the *import* of ``socket`` is now gone
from the whole tree, so this property can be asserted over the app rather than carved around
two lines nobody ran.
"""

from __future__ import annotations

import ast
import re
import subprocess
from pathlib import Path

STREAMLIT_APP = Path(__file__).resolve().parents[1]
REPO_ROOT = STREAMLIT_APP.parent

#: The one directory under ``streamlit_app/`` whose Python reaches no user's machine.
#: ``tests/`` is stripped from both bundles, asserted independently by
#: ``parity/test_parity.py::test_harness_is_excluded_from_packaged_builds``.
#:
#: ``build/`` is deliberately **not** here, and that is the whole of #257's argument rather
#: than an oversight. It holds no Python at all now — the build machinery is shell, the mac
#: launcher is Swift, and the alternative Python launcher that used to sit there was the
#: tree's only ``import socket``. Scanning ``build/`` is what makes the deletion stick:
#: restore that file, or write a new Python launcher beside it, and this guard goes red
#: rather than carving an exception around a file nobody ships.
UNSHIPPED_DIRS = ("tests",)

#: Modules whose whole purpose is to talk to another host. Dotted entries match themselves
#: and anything beneath them; ``urllib`` is deliberately **not** here, because
#: ``urllib.parse`` is string manipulation and refusing it would be refusing arithmetic.
NETWORK_MODULES = (
    "aiohttp",
    "boto3",
    "botocore",
    "ftplib",
    "http",
    "httpx",
    "imaplib",
    "nntplib",
    "paramiko",
    "poplib",
    "pycurl",
    "requests",
    "smtplib",
    "socket",
    "socketserver",
    "ssl",
    "telnetlib",
    "urllib.error",
    "urllib.request",
    "urllib.robotparser",
    "urllib3",
    "websocket",
    "websockets",
    "xmlrpc",
)

#: Call names that open a connection whatever the module they were reached through was
#: called. This is the half of the scan that survives an alias, a re-export, or a bare
#: ``import urllib`` leaning on a submodule some dependency already imported.
NETWORK_CALL_NAMES = (
    "create_connection",
    "getaddrinfo",
    "open_connection",
    "socketpair",
    "urlopen",
    "urlretrieve",
)

#: The pandas readers that accept a URL where a path is expected. A remote literal handed to
#: one of these is an outbound fetch with no network import anywhere near it.
PANDAS_READERS = (
    "read_csv",
    "read_excel",
    "read_feather",
    "read_fwf",
    "read_html",
    "read_json",
    "read_orc",
    "read_parquet",
    "read_pickle",
    "read_sas",
    "read_stata",
    "read_table",
    "read_xml",
)

#: URL schemes that name another host.
REMOTE_SCHEMES = ("http://", "https://", "ftp://", "ftps://", "s3://", "gs://", "sftp://")

#: Command-line tools that move bytes between machines. ``lsof``, ``netstat`` and
#: ``findstr`` are absent on purpose: they read local kernel state, which is how both
#: shipped launchers find a free port, and neither opens a connection.
NETWORK_TOOLS = (
    "bitsadmin",
    "curl",
    "ftp",
    "ncat",
    "netcat",
    "nslookup",
    "scp",
    "sftp",
    "ssh",
    "telnet",
    "wget",
)

#: PowerShell and ``certutil`` reach the network through an argument rather than through the
#: name of the program, so they are matched as phrases.
NETWORK_PHRASES = (
    "invoke-webrequest",
    "invoke-restmethod",
    "downloadfile",
    "downloadstring",
    "start-bitstransfer",
    "urlcache",
)

#: Where the shipped launchers may point a network tool: at this machine.
#:
#: ``0.0.0.0`` is **not** here. It reads as "local" and means every interface — the exact
#: confusion issue #182 settled by keeping the loopback bind after #173's machinery came out,
#: and Streamlit's unset ``server.address`` does the same thing while looking innocent.
LOOPBACK_MARKERS = ("127.0.0.1", "localhost", "::1")

#: The only package fetches the shipped launchers may perform — the first-run dependency
#: install, and pip upgrading itself before it. Anything else is a package this app decided
#: to go and get at launch, which is the shape an update check takes here.
#:
#: Both markers are narrow on purpose. A bare ``"requirements"`` would admit
#: ``pip install --upgrade mafigate -c requirements.txt``, and
#: ``"install --upgrade pip"`` was redundant beside ``"--upgrade pip"`` — a marker that
#: cannot be the deciding one is a marker nobody will maintain.
PIP_INSTALL_MARKERS = ("-r ", "--upgrade pip")

#: Where one shell line stops being one command. Splitting here is what stops a loopback
#: address *anywhere on the line* from vouching for every other command on it — see
#: :func:`outbound_shell_commands`.
COMMAND_SEPARATORS = re.compile(r"&&|\|\||[;|&\n]")

#: A URL, and therefore a host, inside a command.
URL_TOKEN = re.compile(r"[a-z][a-z0-9+.\-]*://[^\s\"'`;|&()]+")

#: Shell that a maintainer runs and no user ever does: the two build scripts, which download
#: a portable Python. Named rather than pattern-matched, and asserted to still exist, so a
#: rename turns this red instead of silently exempting the file that replaced it.
BUILD_ONLY_SHELL = (
    "build/mac/build_dmg.sh",
    "build/windows/build_installer.bat",
)

#: Network APIs in the mac bundle's native launcher. It is compiled into the ``.app`` and is
#: the executable the user actually double-clicks, so it reaches a user's machine as surely as
#: the Python does — and being neither Python nor shell, nothing else here would read it.
SWIFT_NETWORK_APIS = (
    "URLSession",
    "URLRequest",
    "NSURLConnection",
    "NWConnection",
    "NWBrowser",
    "CFStream",
    "getaddrinfo",
    "socket(",
    "Network.framework",
    "import Network",
)


def _tracked(*, under=""):
    """Tracked paths, repo-relative, optionally under one prefix.

    Read through git rather than :meth:`Path.rglob` because ``.claude/worktrees`` holds full
    checkouts of this same tree: a filesystem walk from the repo root answers questions
    about a branch this test is not running on.

    **The limit this puts on the claim, stated because it is not obvious.** An *untracked*
    file is invisible here, and ``build_dmg.sh`` copies directories with ``cp -R``, so a
    scratch file left in ``components/`` would reach a DMG cut on that machine without this
    guard ever seeing it. The alternative — walking the app directory — trades that for
    failing on every developer's working copy, which is how a guard gets deleted by whoever
    meets it. So the scan is over what the repository ships, and the mutation sweep for
    issue #257 demonstrates the boundary rather than hiding it: restoring
    the deleted launcher scores SURVIVED until the file is added to the index, and CAUGHT the
    moment it is.
    """
    completed = subprocess.run(
        ["git", "ls-files"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=True,
    )
    paths = [line for line in completed.stdout.splitlines() if line.startswith(under)]
    if not paths:
        raise AssertionError(
            f"git ls-files listed nothing under {under!r} — this guard would then pass over "
            "any tree at all, which is the one failure mode a privacy assertion cannot have"
        )
    return paths


def shipped_shell_files():
    """Every tracked shell entry point that reaches a user's machine, as app-relative paths.

    **Derived, not listed**, because ``tests/README.md``'s house pattern is *derive it or
    guard it — never copy it*, and a copied inventory of launch routes is precisely the thing
    this repository has already been bitten by: #173 found five publishing launch routes
    against three that were correct, because nobody had a list of the routes.

    So the rule is every tracked ``.sh``, ``.bat`` and ``Makefile`` under the app, minus
    ``tests/``, minus the two named build scripts. A ``build/linux/launch.sh`` arriving
    tomorrow is scanned without anyone editing this file.
    """
    shipped = []
    for relative in _tracked(under="streamlit_app/"):
        path = Path(relative).relative_to("streamlit_app")
        if path.parts[0] == "tests" or path.as_posix() in BUILD_ONLY_SHELL:
            continue
        if path.suffix in {".sh", ".bat"} or path.name == "Makefile":
            shipped.append(path.as_posix())
    return sorted(shipped)


def shipped_python_files():
    """Every tracked Python file that reaches a user's machine, as app-relative paths."""
    shipped = []
    for relative in _tracked(under="streamlit_app/"):
        path = Path(relative).relative_to("streamlit_app")
        if path.suffix != ".py" or path.parts[0] in UNSHIPPED_DIRS:
            continue
        shipped.append(path.as_posix())
    return sorted(shipped)


def _dotted(node):
    """``a.b.c`` out of an attribute chain, or ``None`` if it is not a plain one."""
    parts = []
    while isinstance(node, ast.Attribute):
        parts.append(node.attr)
        node = node.value
    if not isinstance(node, ast.Name):
        return None
    parts.append(node.id)
    return ".".join(reversed(parts))


def _is_network_module(name):
    return any(name == denied or name.startswith(denied + ".") for denied in NETWORK_MODULES)


def _string(node):
    """The literal value of a string node — plain, or the *leading* piece of an f-string.

    The leading piece is what settles a scheme: ``f"https://{host}/version"`` starts with
    ``https://`` whatever ``host`` turns out to be. Only the first piece, because a later one
    says nothing about where the request goes.
    """
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        return node.value
    if isinstance(node, ast.JoinedStr) and node.values:
        return _string(node.values[0])
    return None


def network_reaches(source, *, filename="<source>"):
    """Every way ``source`` reaches the network, as ``(line, what)`` pairs.

    One function for all four Python shapes, and it takes text rather than a path so the
    seeded-violation test below can drive the *same* extractor the real scan uses. A guard
    whose demonstration runs through a second, simpler copy of its logic has demonstrated
    the copy.
    """
    found = []
    tree = ast.parse(source, filename=filename)
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                if _is_network_module(alias.name):
                    found.append((node.lineno, f"import {alias.name}"))
        elif isinstance(node, ast.ImportFrom):
            if node.level or not node.module:  # relative: first-party by definition
                continue
            for alias in node.names:
                for candidate in (node.module, f"{node.module}.{alias.name}"):
                    if _is_network_module(candidate):
                        found.append((node.lineno, f"from {node.module} import {alias.name}"))
                        break
        elif isinstance(node, ast.Call):
            found.extend(_call_reaches(node))
    return sorted(set(found))


#: The ways shipped Python could hand the job to a program instead of doing it itself. No
#: shipped module imports ``subprocess`` today, so this rule is latent — and it is here
#: because latent is exactly what an import scan is blind to: ``subprocess.run(["curl", url])``
#: needs no network module at all.
SHELL_OUT_CALLS = (
    "call",
    "check_call",
    "check_output",
    "popen",
    "Popen",
    "run",
    "system",
)


def _literal_strings(node):
    """Every string literal directly inside a call's arguments, list and tuple items too."""
    literals = []
    for argument in [*node.args, *(keyword.value for keyword in node.keywords)]:
        if isinstance(argument, (ast.List, ast.Tuple)):
            literals.extend(_string(item) for item in argument.elts)
        else:
            literals.append(_string(argument))
    return [literal for literal in literals if literal]


def _call_reaches(node):
    """The call-shaped reaches: a connecting call, a named import, a remote read, a shell-out."""
    name = _dotted(node.func) or ""
    tail = name.rsplit(".", 1)[-1]

    if tail in NETWORK_CALL_NAMES:
        return [(node.lineno, f"call to {name}()")]

    if tail in SHELL_OUT_CALLS:
        for literal in _literal_strings(node):
            for word in re.findall(r"[\w.-]+", literal.lower()):
                if word in NETWORK_TOOLS or word == "pip":
                    return [(node.lineno, f"{name}() running {word!r}")]
            if literal.lower().startswith(REMOTE_SCHEMES):
                return [(node.lineno, f"{name}() reaching {literal!r}")]

    if tail in {"import_module", "__import__"} and node.args:
        named = _string(node.args[0])
        if named and _is_network_module(named):
            return [(node.lineno, f"{tail}({named!r})")]

    if tail in PANDAS_READERS:
        # Keyword arguments as well as the first positional: ``read_csv(url)`` and
        # ``read_csv(filepath_or_buffer=url)`` fetch the same bytes, and a scan that knows
        # only one spelling of an equivalent call is blind to whichever one gets written.
        for argument in [*node.args, *(keyword.value for keyword in node.keywords)]:
            target = _string(argument)
            if target and target.lower().startswith(REMOTE_SCHEMES):
                return [(node.lineno, f"{tail}() reading {target!r}")]

    return []


def _runnable_lines(text):
    """The lines of a shell script or Makefile that *run* something.

    Comments and ``echo`` are stripped first, and this is not tidiness. A well-commented
    file is the adversary of any text scan here: ``launch.sh``'s own comment says *"An
    internet connection is required on first launch"*, and the Makefile prints
    ``pip install pytest-cov`` inside an ``echo`` where no install happens. Issue #162's
    launcher guard was green over a deleted dependency check for exactly this reason — it
    matched the header comment that named the file instead of the line that ran it.
    """
    runnable = []
    for number, raw in enumerate(text.splitlines(), start=1):
        line = raw.strip()
        if not line:
            continue
        if line.startswith(("#", ";", "::", "rem ", "REM ", "@rem", "@REM")):
            continue
        stripped = line.lstrip("@-\t ")
        head = stripped.split(maxsplit=1)[0].lower() if stripped else ""
        if head in {"echo", "echo.", "printf", "@echo"}:
            continue
        runnable.append((number, stripped))
    return runnable


def _hosts(command):
    """The hosts a command names, out of the URLs in it.

    The *host*, not the line. An earlier draft asked whether a loopback marker appeared
    anywhere in the line, which three lines got past in review: an update check appended to
    ``launch.sh``'s existing health-check line, ``curl https://example.org -o
    /tmp/localhost.json``, and ``wget https://example.org --bind-address=127.0.0.1``. In each
    one a loopback token sat on the line while the bytes came from somewhere else.
    """
    hosts = []
    for url in URL_TOKEN.findall(command):
        authority = url.split("://", 1)[1].split("/")[0]
        hosts.append(authority.rsplit("@", 1)[-1])
    return hosts


def outbound_shell_commands(text):
    """Every runnable command that moves bytes to a host that is not this one.

    ``(line, command)`` pairs, judged **per command** rather than per line: a line is split on
    ``;``, ``&&``, ``||`` and ``|`` first, so one command's loopback address cannot vouch for
    the next command on the same line.

    A command invoking a network tool passes only by naming a loopback host *in its own URL*.
    A command with no URL at all is reported too — the target of a fetch has to be visible on
    the command that makes it, or this guard is reading a variable it cannot follow. A command
    invoking pip passes only by installing the pinned requirements, or pip itself.
    """
    offenders = []
    for number, line in _runnable_lines(text):
        for command in COMMAND_SEPARATORS.split(line):
            command = command.strip()
            if not command:
                continue
            lowered = command.lower()
            tools = [
                tool for tool in NETWORK_TOOLS if re.search(rf"(?<![\w.-]){tool}\b", lowered)
            ]
            tools += [phrase for phrase in NETWORK_PHRASES if phrase in lowered]
            if tools:
                hosts = _hosts(lowered)
                remote = [
                    host
                    for host in hosts
                    if not any(marker in host for marker in LOOPBACK_MARKERS)
                ]
                if remote or not hosts:
                    offenders.append((number, command))
            # ``_pip\b`` rather than ``_pip``, or every command mentioning ``_pipeline`` — and
            # this tree mentions the pipeline constantly — reads as a pip invocation. Checked
            # alongside the tools rather than instead of them: an earlier draft skipped this
            # whole branch whenever a tool matched, so `curl http://127.0.0.1/health && pip
            # install --upgrade mafigate` was read as a health check and nothing else.
            if re.search(r"(?<!\w)pip\b|_pip\b", lowered) and "install" in lowered:
                if not any(marker in lowered for marker in PIP_INSTALL_MARKERS):
                    offenders.append((number, command))
    return offenders


def swift_network_apis(source):
    """Every network API named in Swift source, as ``(line, api)`` pairs.

    A plain text scan rather than a parse, and that is a real limitation: a comment naming
    ``URLSession`` would be reported. It is the right trade here — there is no Swift parser in
    a suite whose CI installs ``pytest pandas`` and nothing else, the file is 60 lines of
    process supervision, and the failure mode of this direction is a maintainer having to
    reword a comment rather than a phone-home shipping unseen.
    """
    found = []
    for number, line in enumerate(source.splitlines(), start=1):
        for api in SWIFT_NETWORK_APIS:
            if api in line:
                found.append((number, api))
    return found


def _read(app_relative):
    return (STREAMLIT_APP / app_relative).read_text(encoding="utf-8", errors="replace")


def test_the_scope_is_not_empty_and_covers_what_the_installers_copy():
    """Anti-vacuity, in the direction that matters: a scan of nothing passes everything.

    The scan takes the whole app tree minus :data:`UNSHIPPED_DIRS`, which is deliberately
    *wider* than either installer's copy list — so a new shipped directory is picked up
    without anyone editing this file. What that leaves worth asserting is the narrow claim in
    the other direction: **no installer copies something this module declared unshipped.**
    The day ``tests/`` or ``build/`` appears in a copy list, the exclusion above stops being
    a description of the bundles and becomes a hole in this guard, and that is the line that
    goes red. The parse itself is asserted too, because a regex that matches nothing would
    make the claim vacuously true.
    """
    shipped = shipped_python_files()
    # 52 at the time of writing, across eight roots. The floor catches the enumeration
    # *collapsing* — reading one directory, or none — and deliberately not a count an
    # ordinary afternoon moves. It was 40 in the first draft, which is under half a package
    # away from 52: declaring `config/` unshipped dropped the count to 39, so the floor fired
    # instead of the assertion that was supposed to catch it and neither hole check was ever
    # demonstrated. A floor close enough to be tripped by a normal edit is a floor that
    # answers for its neighbours.
    assert len(shipped) >= 20, (
        f"only {len(shipped)} shipped Python files were found; the app had 52 when this was "
        "written, so the enumeration has broken and every assertion below is now reading "
        "almost nothing"
    )
    # A floor rather than an inventory, and the distinction is the point: the list is derived,
    # so it cannot fall behind a *new* launcher, but BUILD_ONLY_SHELL could be widened until it
    # swallowed an existing one. These two are what the .dmg and the .exe actually run.
    shell = shipped_shell_files()
    for launcher in ("build/mac/MAFigate.app/Contents/MacOS/launch.sh", "build/windows/launch.bat"):
        assert launcher in shell, (
            f"{launcher} is what an installer runs on a user's machine and it is not being "
            f"scanned. The derived list read {shell} — check BUILD_ONLY_SHELL has not grown to "
            "cover it."
        )
    assert "MAFigate.py" in shipped
    assert not [path for path in shipped if path.startswith("tests/")]

    dmg = (STREAMLIT_APP / "build" / "mac" / "build_dmg.sh").read_text(encoding="utf-8")
    copy_list = re.search(r"for item in ([^;]+); do", dmg)
    assert copy_list, (
        "could not find the app-source copy list in build_dmg.sh — the scope check below "
        "cannot tell whether a shipped directory is being scanned, so it is now vacuous"
    )

    iss = (STREAMLIT_APP / "build" / "windows" / "installer.iss").read_text(encoding="utf-8")
    sources = re.findall(r'Source:\s*"\.\.\\\.\.\\([^"*\\]+)', iss)
    assert sources, "could not find any app-source Source: directive in installer.iss"

    copied = {item for item in copy_list.group(1).split() if not item.endswith(".txt")}
    copied |= {source.rstrip("\\") for source in sources if not source.endswith(".txt")}
    assert "components" in copied and "vendor" in copied, (
        f"the copy lists parsed as {sorted(copied)}, which does not include the app's own "
        "packages — the two parses have drifted from the scripts and this check is now "
        "comparing the scan against almost nothing"
    )

    declared_unshipped = sorted(item for item in copied if item in UNSHIPPED_DIRS)
    assert not declared_unshipped, (
        f"an installer now copies {declared_unshipped} onto a user's machine, and this guard "
        "excludes it from the scan as unshipped. The exclusion has become a hole: either the "
        "copy list is wrong, or UNSHIPPED_DIRS is."
    )

    # There was a third assertion here — that no copied directory holding tracked Python was
    # outside the scanned set — and it went out because it could not fail. The scanned set is
    # *everything* tracked bar UNSHIPPED_DIRS, so that difference is a subset of
    # UNSHIPPED_DIRS by construction and the assertion above always fired first. Review
    # proved it: adding `tests` to the DMG copy list failed on the line above, never on that
    # one. An assertion that cannot be the one to fail is not a weak check, it is decoration.
    #
    # Answered off the tracked list rather than the filesystem: ``resources/`` is not tracked
    # in this repository at all — it arrives with the pipeline checkout, which is why both
    # installers copy it conditionally and why it is absent from a fresh worktree. A check
    # that read the disk would say something different depending on which checkout it ran in.
    tracked_py_roots = {
        Path(path).relative_to("streamlit_app").parts[0]
        for path in _tracked(under="streamlit_app/")
        if path.endswith(".py")
    }
    assert "resources" not in tracked_py_roots, (
        "resources/ is now tracked here and both installers copy it, so its Python ships and "
        "has to be scanned — it was excluded only because there was nothing to read"
    )


def test_no_shipped_module_reaches_the_network():
    """The property itself: nothing that ships imports, names, or calls a network client.

    See the module docstring for the scope and for #229, which is the decision a red line
    here is asking you to reopen.
    """
    offenders = {}
    for relative in shipped_python_files():
        reaches = network_reaches(_read(relative), filename=relative)
        if reaches:
            offenders[relative] = reaches

    assert not offenders, "\n".join(
        [
            "The shipped app reaches the network. Issue #229 decided MAFigate makes no "
            "outbound call — no update check, ever — because the first one falsifies the "
            "data-posture promise the README makes to every delivery channel:",
            *(
                f"  streamlit_app/{path}:{line} — {what}"
                for path, reaches in sorted(offenders.items())
                for line, what in reaches
            ),
            "If this is deliberate, #229's answer and the README's promise change first.",
        ]
    )


def test_the_shipped_launchers_talk_only_to_loopback():
    """The other half of what ships: the shell entry points, no remote command lines.

    Python is not the only way to phone home, and the launcher is where an update check
    would most naturally land — it already runs before the app, already talks to pip, and
    already has a log file to hide a failure in.
    """
    assert all((STREAMLIT_APP / path).exists() for path in BUILD_ONLY_SHELL), (
        f"a script named in BUILD_ONLY_SHELL is gone: {BUILD_ONLY_SHELL}. Those two names "
        "exempt maintainer-side build scripts from this scan, so a stale name exempts "
        "whatever file replaced it"
    )

    offenders = {}
    for relative in shipped_shell_files():
        found = outbound_shell_commands(_read(relative))
        if found:
            offenders[relative] = found

    assert not offenders, "\n".join(
        [
            "A shipped launcher runs a command against another host. The first-run "
            "`pip install -r requirements.txt` is the one allowed fetch (see "
            "PIP_INSTALL_MARKERS); a loopback address is the other. Neither of these is:",
            *(
                f"  streamlit_app/{path}:{line} — {command}"
                for path, found in sorted(offenders.items())
                for line, command in found
            ),
        ]
    )


def test_the_native_launcher_reaches_nothing():
    """The third language that ships, and the one nothing else here would read.

    ``MAFigateLauncher.swift`` is compiled into the ``.app`` and is the executable the user
    double-clicks; today it supervises ``launch.sh`` and opens nothing. Review found it
    unscanned by every rule above, which made the "no outbound call" claim true of the app's
    Python and shell and silent about its actual entry point.
    """
    swift = [
        Path(path).relative_to("streamlit_app").as_posix()
        for path in _tracked(under="streamlit_app/")
        if path.endswith(".swift")
    ]
    assert swift, (
        "no Swift file is tracked under streamlit_app/ any more. If the native launcher is "
        "gone this test is checking nothing; if it moved, follow it."
    )
    offenders = {path: swift_network_apis(_read(path)) for path in swift}
    offenders = {path: found for path, found in offenders.items() if found}
    assert not offenders, "\n".join(
        [
            "The native launcher names a network API. #229 decided MAFigate makes no outbound "
            "call, and this is the file the user's double-click actually runs:",
            *(
                f"  streamlit_app/{path}:{line} — {api}"
                for path, found in sorted(offenders.items())
                for line, api in found
            ),
        ]
    )


def _unnoticed(cases, scanner):
    """The labels in ``cases`` whose source the scanner reported nothing about."""
    return sorted(label for label, source in cases.items() if not scanner(source))


def _wrongly_refused(cases, scanner):
    """The labels in ``cases`` the scanner reported something about, and should not have."""
    return sorted(label for label, source in cases.items() if scanner(source))


def test_the_scanner_refuses_every_shape_it_claims_to_catch():
    """Each rule, seeded and watched refusing — through the extractor the real scan uses.

    A guard that has never refused anything has not been demonstrated, and this suite has a
    documented history of tests that passed because they had stopped applying. The manual
    demonstration #257 asked for — an outbound call added to a real app module, watched
    going red, removed — is recorded in that issue's notes; this test is the part
    that keeps being run after that afternoon.
    """
    seeded = {
        "a plain import": "import requests\n",
        "a submodule import": "import urllib.request\n",
        "a from-import of a submodule": "from urllib import request\n",
        "a from-import of a member": "from http.client import HTTPSConnection\n",
        "a dynamically named client": "import importlib\nimportlib.import_module('socket')\n",
        "a connecting call": "def go(loop):\n    return loop.open_connection('example.org', 443)\n",
        "a remote read": "import pandas as pd\npd.read_csv('https://example.org/x.csv')\n",
        "a remote read by keyword": (
            "import pandas as pd\npd.read_csv(filepath_or_buffer='https://example.org/x')\n"
        ),
        # Below this line: the shapes review found surviving. Each needed no network import at
        # all, which is the blind spot an import scan has by construction.
        "shelling out to curl": (
            "import subprocess\nsubprocess.run(['curl', 'https://example.org/version'])\n"
        ),
        "shelling out to pip": (
            "import subprocess\n"
            "subprocess.run(['pip', 'install', '--upgrade', 'mafigate'])\n"
        ),
        "os.system with wget": "import os\nos.system('wget https://example.org/v.json')\n",
        "a remote url handed to a runner": (
            "import subprocess\nsubprocess.check_output(['fetcher', 'https://example.org'])\n"
        ),
    }
    survivors = _unnoticed(seeded, network_reaches)
    assert not survivors, (
        f"the Python scanner did not notice: {survivors}. Every entry above is a way to "
        "reach another host, so a shape it cannot see is a shape that ships unnoticed."
    )

    seeded_shell = {
        "a remote fetch": 'curl -s "https://example.org/latest-version"\n',
        "a download": "wget https://example.org/mafigate.json\n",
        "a package fetch": "pip install --upgrade mafigate\n",
        "a powershell fetch": "powershell -c Invoke-WebRequest https://example.org\n",
        # And here the four the per-line draft let through, all of which put a loopback token
        # on a line whose bytes came from somewhere else.
        "a fetch chained onto the health check": (
            'curl -s "https://example.org/latest" > /dev/null 2>&1; '
            'curl -s "http://127.0.0.1:${PORT}/_stcore/health"\n'
        ),
        "a fetch written to a file called localhost": (
            "curl -s https://example.org/latest -o /tmp/localhost.json\n"
        ),
        "a fetch bound to a loopback interface": (
            "wget https://example.org/v.json --bind-address=127.0.0.1\n"
        ),
        "a package fetch after a health check": (
            'curl http://127.0.0.1/health && "$PIP" install --upgrade mafigate\n'
        ),
        "a fetch whose target is hidden in a variable": 'curl -s "${UPDATE_URL}"\n',
    }
    shell_survivors = _unnoticed(seeded_shell, outbound_shell_commands)
    assert not shell_survivors, (
        f"the shell scanner did not notice: {shell_survivors}. The launcher is the likeliest "
        "home for an update check, so these are the lines it most needs to see."
    )

    seeded_swift = {
        "a URLSession fetch": "let task = URLSession.shared.dataTask(with: url)\n",
        "the Network framework": "import Network\n",
    }
    swift_survivors = _unnoticed(seeded_swift, swift_network_apis)
    assert not swift_survivors, (
        f"the Swift scanner did not notice: {swift_survivors}, and the native launcher is "
        "what the user's double-click runs."
    )


def test_the_scanner_leaves_the_app_own_correct_behaviour_alone():
    """The other direction, which is what stops the guard being deleted by whoever meets it.

    Every line here is something the app or a launcher legitimately does today. A guard that
    refused any of them would be refusing correct code, and would be rewritten around rather
    than fixed.
    """
    allowed = {
        "URL parsing": "from urllib.parse import quote\n",
        "a rendered ClinVar link": (
            "import streamlit as st\n"
            "st.markdown('[ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/?term=BRCA1)')\n"
        ),
        # A relative path, and deliberately not an absolute one under a home directory:
        # #247's sweep scans every tracked file that would be published, and its home-paths
        # rule caught the first draft of this line — twice, since the comment explaining the
        # fix quoted the shape it was avoiding.
        "a local read": "import pandas as pd\npd.read_csv('sample.maf', sep='\\t')\n",
        "a docstring naming a client": '"""Do not import requests here."""\n',
        "a local subprocess call": (
            "import subprocess\nsubprocess.run(['git', 'ls-files'], check=True)\n"
        ),
        "a local pipeline invocation": (
            "import subprocess\nsubprocess.run([python, 'bin/filter_variants.py', '--maf', p])\n"
        ),
    }
    tripped = _wrongly_refused(allowed, network_reaches)
    assert not tripped, (
        f"the Python scanner refused correct code: {tripped}. The app renders links the "
        "browser follows and parses URLs; neither is a call it makes."
    )

    allowed_shell = {
        "the loopback health check": 'curl -s "http://127.0.0.1:${PORT}/_stcore/health"\n',
        "the first-run install": '"${VENV_PIP}" install -q -r "${REQUIREMENTS}"\n',
        "pip upgrading itself": '"${VENV_PIP}" install --upgrade pip\n',
        "a printed hint": 'echo "   requirements.txt lists it: pip install pytest-cov"\n',
        "a commented-out fetch": "# curl https://example.org would be an update check\n",
        "the local port probe": 'netstat -an | findstr ":%PORT% "\n',
        "the health check with a redirect": (
            'if curl -s "http://127.0.0.1:${PORT}/_stcore/health" > /dev/null 2>&1; then\n'
        ),
        "the local port probe on mac": 'if ! lsof -i :${port} > /dev/null 2>&1; then\n',
        "a make recipe installing the requirements": (
            "\t@$(PYTHON) -m pip install -r requirements.txt\n"
        ),
    }
    tripped_shell = _wrongly_refused(allowed_shell, outbound_shell_commands)
    assert not tripped_shell, (
        f"the shell scanner refused correct code: {tripped_shell}. Every one of these is a "
        "line a shipped launcher runs, or a comment about one."
    )

    allowed_swift = {
        "process supervision": (
            "let task = Process()\ntask.executableURL = URL(fileURLWithPath: script)\n"
        ),
        "a local file URL": 'let url = URL(fileURLWithPath: "/tmp/mafigate.log")\n',
    }
    tripped_swift = _wrongly_refused(allowed_swift, swift_network_apis)
    assert not tripped_swift, (
        f"the Swift scanner refused correct code: {tripped_swift}. `URL(fileURLWithPath:)` and "
        "`Process()` are how the native launcher starts a local script."
    )


def test_the_scanner_finds_the_real_network_client_outside_the_app():
    """The scan, demonstrated on real code rather than on a fixture, at the scope boundary.

    ``update_db/`` fetches annotation databases over FTP and HTTP. It is out of scope — it
    is the pipeline's tooling and ships with neither installer — and it is the reason the
    boundary has to be deliberate rather than lucky: point this scanner one directory wider
    and it goes red on code that is doing its job.
    """
    outside = [path for path in _tracked(under="update_db/") if path.endswith(".py")]
    reaching = {
        path: network_reaches((REPO_ROOT / path).read_text(encoding="utf-8", errors="replace"))
        for path in outside
    }
    reaching = {path: found for path, found in reaching.items() if found}

    assert reaching, (
        "no file under update_db/ imports a network client any more, so this module has "
        "lost its live demonstration that the scanner works on real code. Point it at "
        "another out-of-scope file that does, or delete this test and keep the seeded one."
    )
