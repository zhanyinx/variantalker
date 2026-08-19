"""Which build of MAFigate this is: the route it arrived by, and the commit it was cut from.

Issue #263. A version number cannot identify a build. The same ``APP_VERSION`` reaches a
user as a macOS .dmg, as a Windows .exe or as a clone of the source, and with **no update
check ever** — decided in #229, because such a check would be this app's first outbound call
and would falsify its own privacy claim — the About dialog is the only thing that tells a
user, and therefore the maintainer reading their bug report, what they are actually running.

Two facts, one module
---------------------
The channel and the commit are written by whichever build script produced the artifact, into
a **gitignored** ``config/build_stamp.py`` (see ``build/build_stamp.py``, which writes it, and
``config/.gitignore``, which keeps it out of git). This module is the only reader of that
file, and the reason it exists rather than the app importing the stamp directly is the case
where the stamp is **absent**: a plain clone, a ``setup.sh`` install, and every test run and
CI job. That is not an error state to be handled at each call site — it is the ordinary state
of the source checkout channel, and it has a name here rather than a ``try`` at each reader.

So the absence *is* the answer: no stamp means :data:`SOURCE_CHANNEL`, and nothing writes
that value. A build script that stamped ``source-checkout`` would give one channel two
representations and leave the default path — the one every contributor takes — unexercised.

Why the raw tokens reach the screen
-----------------------------------
About prints ``channel windows-installer · build a1b2c3d`` rather than prose like *Windows
installer*. Three reasons, and the first is the load-bearing one: a label table would be a
second spelling of the same fact, and ``tests/test_app_identity.py``'s one-surface sweep can
only look for a channel it can name. Second, this line exists to be quoted into a bug report,
where an unambiguous token beats a sentence. Third, an unfamiliar channel is passed through
as written (see :func:`identify`) — a mis-stamped build says something visibly odd instead of
being quietly rewritten into a claim that is false.

``tests/test_build_identity.py`` holds the two halves together: that the writer's channels
are these channels, that the stamp is never committed, that both build scripts write it, and
that this module answers *source checkout* when it is not there.
"""

from __future__ import annotations

import importlib
from typing import Any, NamedTuple

#: The module a build script writes and this one reads. Named as a string rather than
#: imported at module scope on purpose: it is absent in every source checkout, which is the
#: majority case, and an ``ImportError`` at import time is not a thing a page should catch.
STAMP_MODULE = "config.build_stamp"

#: What the app reports when nothing stamped this tree. Not a placeholder for a real channel —
#: running the working tree *is* a way MAFigate reaches a user (#229's clone route), and it is
#: the way every contributor and every CI run reaches it.
SOURCE_CHANNEL = "source-checkout"

#: What the app reports for a build identifier it was not given. Its own value rather than an
#: empty string, so About always names three facts and a bug report is never missing a field
#: silently; and its own *word*, so the one-surface sweep has something to look for.
#:
#: A build machine with no ``git`` writes an empty commit and this fills in, which is why the
#: sentinel lives here and not in the writer: one spelling, in the module that renders it.
UNRECORDED_BUILD = "unrecorded"

#: The channels a build script may stamp — the two installers, and nothing else. There is no
#: hosted channel: the streamlit.io deployment was removed and #229 ruled it out of scope
#: permanently, so a stamp naming one would describe a route that does not exist.
#:
#: ``build/build_stamp.py`` declares this list too, and refuses any other value. The two are
#: held equal by ``tests/test_build_identity.py`` in both directions rather than shared,
#: because that writer is stdlib-only and must run on a build machine before this app is
#: importable at all — the same rule ``build/version.py`` and ``vendor/_sync.py`` follow.
INSTALLER_CHANNELS = ("macos-dmg", "windows-installer")


class BuildIdentity(NamedTuple):
    """A channel and a build identifier, both always present and always non-empty."""

    channel: str
    build: str

    def __str__(self) -> str:
        """The one line About prints, and the one line a bug report quotes."""
        return f"channel {self.channel} · build {self.build}"


def _stated(value: Any) -> str:
    """*value* if it is a non-blank string, else ``""``.

    A stamp is a generated file, so the failures worth surviving are the generated ones: a
    key the writer stopped writing, a value it wrote empty because ``git`` could not answer,
    or something that is not a string at all. Each of those reads as *not stated* and falls
    back, rather than putting ``None`` or a tuple on a clinical screen.
    """
    return value.strip() if isinstance(value, str) and value.strip() else ""


def identify(channel: Any, build: Any) -> BuildIdentity:
    """What a stamp's two values mean, including their absence. The whole rule, and pure.

    An unfamiliar *channel* is passed through as written rather than checked against
    :data:`INSTALLER_CHANNELS`. The writer refuses anything that is not a real route, so the
    only way here with a strange one is a hand-edited generated file — and a maintainer told
    ``channel hosted`` learns that the build is broken, where one told ``channel
    source-checkout`` by a silent rewrite is told something false about a shipped installer.
    """
    return BuildIdentity(
        channel=_stated(channel) or SOURCE_CHANNEL,
        build=_stated(build) or UNRECORDED_BUILD,
    )


def build_identity() -> BuildIdentity:
    """:func:`identify` applied to this checkout's build stamp, if it has one.

    Reached through :func:`importlib.import_module` rather than a ``from … import`` so that
    the absent case is an ordinary ``ImportError`` on a name — which is what a source
    checkout raises — and so a test can stand a stamp up, or take one away, through
    ``sys.modules`` without writing a generated file into the tree it is testing.
    """
    try:
        stamp = importlib.import_module(STAMP_MODULE)
    except ImportError:
        return identify(None, None)

    return identify(
        getattr(stamp, "BUILD_CHANNEL", None),
        getattr(stamp, "BUILD_COMMIT", None),
    )
