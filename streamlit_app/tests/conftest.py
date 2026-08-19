"""Suite-wide fixtures. Currently one, and it exists to protect the developer's own machine.

``config/file_history.py`` writes to ``~/.mafigate`` — the real one — and the app writes it
from ``_open_the_file_just_read``, the tail all three load doors share. That function is
driven for real by several test modules that have nothing to do with the file history:
``test_file_chooser.py``, ``test_data_page_sections.py``, ``test_arm_mismatch.py`` and
``test_messages_reach_the_user.py`` among them. Measured before this file existed: one full
run left seven fixture names — ``germline_reference.maf``, ``sample_b.maf``,
``next_sample.maf`` — in the developer's real ``~/.mafigate/opened_files.json``.

The fix is here rather than in those modules on purpose. Stubbing the store in each one that
happens to drive the tail today leaves the *next* module that drives it writing to a real home
directory, silently and with nothing failing — which is the same shape as the Clear Cache
button that deletes your parameters when a test forgets to stub it. One autouse fixture covers
every test that exists and every test that will.

``test_file_history.py`` asserts that this fixture is in force, so it cannot quietly stop
applying.
"""

from __future__ import annotations

import pytest


@pytest.fixture(autouse=True)
def keep_the_file_history_out_of_your_home(monkeypatch, tmp_path):
    """Point the recents list at a temporary directory for the duration of every test.

    Both constants, not just the file: :func:`config.file_history.record_open_file` calls
    ``mkdir`` on the directory, so redirecting the file alone would still create
    ``~/.mafigate`` on a machine that had none.
    """
    from config import file_history

    home = tmp_path / ".mafigate"
    monkeypatch.setattr(file_history, "HISTORY_DIR", home)
    monkeypatch.setattr(file_history, "HISTORY_FILE", home / "opened_files.json")
