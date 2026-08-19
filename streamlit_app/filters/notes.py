"""What kind of thing the filter is telling the user, carried with the sentence itself.

Every message the filter produces travels to one slot on the data page, and until issue #151
they travelled as bare strings with a single predicate deciding how each was drawn. That gave
two Streamlit levels for messages of three quite different kinds, and the mismatch was
measurable rather than theoretical: on both reference MAFs, on their own arm, under all four
shipped presets, the **only** thing the slot ever drew was a yellow warning box reporting a
population-frequency filter that had worked perfectly — 4 runs out of 4. Set against issue
#136's measurement that the escalated warning fires on the file's own arm for 2 of 173
placeable real MAFs, the slot's ordinary output was a warning about success, so yellow could
not mean anything by the time a user met a real one.

Three levels, keyed to **what the message is about** rather than to how loud it is:

``ERROR``
    Rows are missing from this report and nothing on screen shows their absence. The user
    cannot recover this by re-reading their parameters or the table.
``WARNING``
    This report is not the one the user asked for — a filter input was replaced by a stand-in,
    or a restriction they chose could not be applied.
``INFO``
    The filter ran and did what was asked; this says what it did.

The glyph, and what issue #150 was protecting
--------------------------------------------
This is the one account of it; the builders that changed glyph point here rather than each
retelling it. ``st.error``, ``st.warning`` and ``st.info`` all draw **no icon of their own**
unless one is passed — verified, not assumed — so the emoji in the message text is the only
icon the user sees, and glyph and box can disagree.

Before #151 they did, deliberately. #150 set ``ℹ️`` on two notes meaning *a filter you asked
for did not run* so that the two would agree with **each other**, accepting that both then sat
in a yellow box, because a single predicate was the only level the renderer had. #151 found a
third note meaning the same thing and carrying ``⚠️``, so the agreement #150 bought was already
incomplete. All three are now ``WARNING`` with ``⚠️``: they still agree with each other, and now
with their box too.

So the two warning tiers own a severity glyph — ``❌`` and ``⚠️`` — while the info tier's glyph
marks the **topic** instead (``🌍`` population frequency, ``🧬`` the gene list), which is the one
thing a coloured box cannot say. Guarded as a prohibition rather than a list of permitted
topical glyphs: a list of topics stops applying the moment a filter gains a note.

Why the level is a field and not a predicate
--------------------------------------------
:func:`~filters.absent_columns.is_escalated` reads the phrase ``PATHOGENIC RETENTION
DEGRADED`` out of the message text, which works because that phrase *is* the copy — it is what
the sentence says. Extending the same trick to a second level would have meant testing the
leading emoji, and the emoji is user-facing copy that an editorial pass can legitimately
change. A box that repaints itself because someone improved a glyph is the drift shape issue
#151 asked to avoid. So the module that writes the words also states what kind of note it is,
and the renderer branches on that with no test on copy at all.

Kept free of pandas and of Streamlit, like its neighbours, so the tests that assert which
level each message carries run without either installed.
"""

from typing import NamedTuple

#: Rows are missing from this report and their absence is invisible on screen.
ERROR = "error"
#: This report is not the one the user asked for.
WARNING = "warning"
#: The filter ran as asked; this note says what it did.
INFO = "info"

#: Every level a note may carry, in descending severity. The renderer is checked against this
#: tuple rather than against a list written beside it, so a fourth level cannot be added here
#: and left undrawn — see ``tests/test_filter_notes.py``.
LEVELS = (ERROR, WARNING, INFO)


class Note(NamedTuple):
    """One message for the filter's slot on the data page, and what kind of message it is.

    A :class:`NamedTuple` rather than a dataclass because these are stashed in
    ``st.session_state`` to survive the rerender that follows a silent load, and a plain tuple
    is the least surprising thing to leave there.

    Attributes:
        level: one of :data:`LEVELS`.
        text: the sentence as the user reads it, glyph included. The glyph stays in the text
            because the two warning tiers mark severity (``❌``, ``⚠️``) while the info tier
            marks *topic* (``🌍`` for population frequency, ``🧬`` for the gene list) — a
            distinction the box cannot carry and the level should not have to encode.
    """

    level: str
    text: str

    @classmethod
    def error(cls, text: str) -> "Note":
        """Rows are missing and the user cannot see that they are."""
        return cls(ERROR, text)

    @classmethod
    def warning(cls, text: str) -> "Note":
        """The report is not the one the user asked for."""
        return cls(WARNING, text)

    @classmethod
    def info(cls, text: str) -> "Note":
        """The filter ran as asked; this is what it did."""
        return cls(INFO, text)
