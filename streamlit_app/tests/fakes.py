"""What the unit suite shares: stand-ins for the parts of Streamlit a mock gets wrong,
and the small readers that keep a shared shape from being copied into every caller.

Not a ``test_*.py``, so pytest does not collect it and ``tests/README.md``'s table does not
name it — the same arrangement as ``tests/parity/harness.py`` and ``tests/parity/contract.py``,
which are libraries the parity modules import rather than instruments in their own right.

This module exists because ``FakeSessionState`` had two readers and lived in one of them.
``tests/test_app_identity.py`` reached into ``tests/test_components.py`` to borrow it, which
made a *collected* module double as a library: the coupling is invisible to
``tests/test_suite_organisation.py``, and a module renamed for reasons of its own would break
a file that has nothing to do with UI components. Borrowing beat copying; the target was
wrong.
"""

from __future__ import annotations


class FakeSessionState(dict):
    """A session_state that behaves like Streamlit's: dict access plus attributes.

    A bare `Mock()` will not do here. Its `.get()` returns the same value for every key,
    so the sidebar's status block reads that one value as the file name, the row frame and
    the parameter dict at once — and a `str` standing in for `filter_params` raises rather
    than failing an assertion.

    ``in`` matters as much as ``.get``: a ``Mock`` raises ``TypeError`` on membership tests,
    and the header checks ``"cache_timestamp" in st.session_state`` before drawing its
    cache-restored banner.
    """

    def __getattr__(self, key):
        try:
            return self[key]
        except KeyError as exc:
            raise AttributeError(key) from exc

    def __setattr__(self, key, value):
        self[key] = value


def note_texts(diagnostics) -> tuple[str, ...]:
    """The prose of a :class:`~filters.variant_filters.Diagnostics`'s notes, without the levels.

    Three modules assert what the filter *says* — ``test_absent_columns``,
    ``test_gene_lists``, ``test_filter_app_extras`` — and since issue #151 the notes are
    ``(level, text)`` pairs rather than strings. Each of those modules owns the wording of some
    sentence and none of them owns the levels: which level a note carries is asserted in
    ``tests/test_filter_notes.py``, next to the renderer it decides. Keeping the two apart is
    what stops a re-levelling from reading as a wording change, and a rephrasing from reading
    as a re-levelling.

    Here rather than in each caller because it arrived as three byte-identical copies, which is
    the one duplication this repo has a documented history with.
    """
    return tuple(note.text for note in diagnostics.notes)


def page_config_kwargs() -> dict:
    """What the app hands ``st.set_page_config``, read by driving the real function.

    Two modules ask this dict a question and they are not the same question:
    ``test_app_identity.py`` asks what the About dialog says the app *is*, and
    ``test_public_repo_name.py`` asks where the two clickable menu items *land*. Both need the
    same shape to ask it — drive ``setup_page_config`` with Streamlit mocked out and keep the
    keyword arguments — so the shape lives here rather than in the second module to arrive.

    The call-count assertion belongs here rather than in either caller. A second
    ``set_page_config`` call is a fact about the app, and reading the first call's keywords
    while a second one existed would quietly answer a question about a page nobody sees.
    """
    from unittest.mock import patch

    import MAFigate

    with patch("MAFigate.st") as mock_st:
        MAFigate.setup_page_config()
        assert mock_st.set_page_config.call_count == 1, (
            "setup_page_config called st.set_page_config "
            f"{mock_st.set_page_config.call_count} times, not once — whichever call this "
            "helper returned, it is no longer the whole story."
        )
        return mock_st.set_page_config.call_args.kwargs
