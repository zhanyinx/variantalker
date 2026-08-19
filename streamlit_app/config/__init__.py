"""
Configuration package for MAFigate application.
"""

from .columns import resolve_visible_columns

# One name, matching what this package exported before: consumers import from
# `config.columns` directly, so anything more here would be re-export with no reader.
#
# `from .constants import *` used to sit above, and it is gone. Nothing read a constant off
# the package — every consumer names the module it wants — and a star import is how a name
# added to `constants.py` becomes visible in three places at once. `config/constants.py` now
# has a stated surface and a test that pins it; a wildcard re-export is the one thing that
# could quietly widen it again.
__all__ = ['resolve_visible_columns']
