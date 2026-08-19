"""What the app calls itself, the pages it has, and how tall its tables are.

This file was 4,021 lines and held five unrelated kinds of thing. Each kind now has its
own home, and the rule for this one is narrow enough to state: **a constant belongs here
only if it is about the application as an application** — its name, its version, its
navigation, its presentation defaults. Anything that is about *variants* belongs to the
domain module for that kind:

* ``config/gene_panels.py`` — the predefined panels and the option list naming them;
* ``config/vocabularies.py`` — the values a filter control may offer, and the normaliser
  that cleans a stored value against them;
* ``config/presets.py`` — the four presets and the settings the app opens with;
* ``config/columns.py`` — which columns a filtered table has, what each one means, and
  which of them a MAF must carry.

The split is by **kind** rather than by consumer, and the reason is not the one a first
census suggested. Counting mentions made every page look like a reader of everything;
counting *uses* showed that ``page_modules/data_loading.py`` imported nineteen of these
constants without reading any of them. On the true count most groups have one or two
readers — so a by-consumer split would have pushed preset dicts and vocabulary lists back
*inside* page modules, which is exactly what made those files oversized. Kind is the axis
that keeps the pages thin.

Seven names were deleted rather than rehoused — ``APP_DESCRIPTION``, ``SUPPORTED_FORMATS``,
``MAX_FILE_SIZE``, ``SIDEBAR_WIDTH``, ``ITEMS_PER_PAGE``, ``COLORS`` and
``EXPORT_FORMATS``. Not one had a reader anywhere in the app or the suite, and two of them
lied about the app while they sat here: ``MAX_FILE_SIZE`` read as an upload limit that has
never been enforced, and ``EXPORT_FORMATS`` named CSV/TSV/Excel as the export formats while
being *imported by three modules and used by none* of them — the app's real export formats
are decided at each download button.

``APP_TAGLINE`` arrived after that sweep and is deliberately *not* a resurrection of
``APP_DESCRIPTION``: the deleted name had no reader, this one has two, and both are the
places the app describes itself out loud (issue #71). It is one sentence in one file
precisely because it was three sentences in three files — a header line, an About blurb
and a sidebar footer, agreeing in meaning and disagreeing in words.

``tests/test_config.py`` pins this module's surface, so the file cannot quietly become a
junk drawer again.
"""

# Application metadata
APP_NAME = "MAFigate"

#: The version, and the **only** place it is written (issue #260). Everything that has to
#: stamp a version on something reads it from here through ``build/version.py``: the
#: Windows installer's version and output filename, the macOS bundle's two ``Info.plist``
#: keys, the DMG's filename, and the ``mafigate-v<version>`` tag a release is cut from.
#: ``tests/test_installer_version.py`` fails if any of them grows a literal again.
#:
#: **It says 1.0.0 and it used to say 2.0.0.** That is not a typo and not a downgrade of
#: anything released: no MAFigate has ever been released. The 2.0.0 was an internal literal
#: that, by this repo's own account of it in ``config/param_migration.py``, "never moved at
#: all" — a number with no lineage behind it. The number that *has* been published is
#: 1.0.0: the public README tells readers the desktop installers are "coming with the first
#: `mafigate-v1.0.0` release", and #229's delivery matrix, this ticket's parent (#240) and
#: the release ticket (#265) all name that tag. Between an unreleased internal literal and
#: a promise already on a public page, the promise wins — a first public release numbered
#: 2.0.0 would also owe its readers a 1.x history that does not exist.
APP_VERSION = "1.0.0"

#: The app's one self-description, written for the clinician who opens it rather than for
#: the developer who built it: it says what the app *does* to a file. Read by
#: ``render_header`` (under the title, on every page) and by the About dialog — those two
#: cannot drift, because neither owns the words.
#:
#: The sentence was taken from the Home page's own expander, which says the same thing at
#: greater length. That is a borrowing, **not** an invariant: ``page_modules/home.py`` holds
#: its own prose and nothing checks the two against each other. Asserting the overlap was
#: considered and rejected — a test comparing two pieces of prose holds the *rendering*,
#: which this suite deliberately does not do, and Home's copy is expected to be edited by
#: later UX work. If they diverge, prefer editing whichever one has stopped being true.
APP_TAGLINE = "Filters annotated MAF files down to the variants you need to look at."

# Page configurations
PAGES = {
    "home": {"title": "🏠 Home", "icon": "🏠", "description": "Welcome to MAFigate"},
    "parameter_config": {
        "title": "⚙️ Parameter Configuration",
        "icon": "⚙️",
        "description": "Configure filtering parameters",
    },
    "data_loading": {
        "title": "📊 Data Loading & Analysis",
        "icon": "📊",
        "description": "Load and analyze MAF files",
    },
    "help": {
        "title": "❓ Help & Documentation",
        "icon": "❓",
        "description": "Column information and user guide",
    },
}

# Presentation defaults
TABLE_HEIGHT = 500  # Increased height for better table display
