"""How a MAF spells its chromosomes, settled once when the file is opened.

A MAF writes ``chr1`` or it writes ``1``, and until issue #149 this app answered that
question six times and got it right at three.

Four sites **added** the prefix — ``components/variant_detail.py``'s UCSC link and its
``Position`` line, and ``components/variant_table.py``'s row selector and *View details*
button — and only the first of those four asked whether the value already had one. Two
more **strip** it, for the two consumers that want the number bare:
``components/charts.py``'s chromosome axis and the gnomAD coordinate link beside the UCSC
one. Those two were right all along, and stay as they are: an unconditional
``.replace("chr", "")`` is correct under either spelling, and correct under this module
too, which only ever makes the prefixed spelling *more* certain.

The answer is not the app's to invent
-------------------------------------
The pipeline already decides this, one step upstream of where MAFigate joins in.
``bin/add_guidelines_escat_to_funcotator.py`` prefixes the whole column while it is
annotating::

    if not any(["chr" in str(x) for x in maf["Chromosome"].values]):
        maf["Chromosome"] = "chr" + maf["Chromosome"].astype(str)

So a MAF that has been all the way through the pipeline carries ``chr1`` **by
construction**, and by the time ``filter_maf`` — the thing this app is the analogue of —
is handed the file, the question has already been settled. This module puts MAFigate at
the same point in the sequence: the reader reads the file exactly as the pipeline's
reader does (``utils.read_maf``, which normalises nothing and says so), and then this
runs, once, before the frame becomes ``maf_data``.

Measured, because the guard looked like dead code and is not
------------------------------------------------------------
Over **187 byte-distinct real MAFs on this machine, 336,170 rows** (the ~290 the corpus
notes record, less the 101 unhydrated OneDrive placeholders):

* **174 files** spell every chromosome ``chr1``.
* **9 files** spell every one bare — 8,056 rows.
* **0 files mix the two**, which is exactly what a whole-column rule upstream predicts.

Seven of the nine bare files are GATK test resources and a TCGA demo. The other two are
the institute's own germline MAFs — 267 real columns with ClinVar and gnomAD annotations,
bare because they stop short of the ESCAT step that would have prefixed them. So a user
can open a bare-spelled file here, and the guard this replaces was written by somebody
who had met one. (Those two are named in the measurement record and not here: their
filenames are case identifiers, which the publication scan gate keeps out of the tree.
Nothing a reader can act on is lost — the files are on one clinical share.)

Three consequences the render sites could not have reached
----------------------------------------------------------
**The note key.** ``components/variant_table._KEY_COLUMNS`` builds a note's identity from
the raw ``Chromosome`` value, so before this the same variant carried *two different note
keys* depending on which spelling its file used — and issue #67 decided a note is about
the variant, for the institute, surfacing across MAFs by design. No render-time helper
could have fixed that; it is not a rendering.

**The export.** ``filtered_data`` descends from ``maf_data``, so a download taken from a
bare-spelled file now spells ``chr1``. That is a change to the user's own data, which is
why ``page_modules/data_loading`` says so on the ~5% of files where this fires rather
than doing it silently.

**The verdict is untouched.** ``vendor/pipeline_filters.py`` never reads ``Chromosome``,
so no keep-or-drop decision moves whichever spelling arrives. This is a rendering and
identity concern only.

Where this deviates from the pipeline's line, and why
-----------------------------------------------------
In one place, deliberately. ``"chr" + column.astype(str)`` turns an empty cell into the
literal string ``chrnan``, which ``config.missing_values.says_nothing`` does not
recognise — so the detail panel would print ``chrnan`` at a clinician where it prints an
em dash today. No row in the corpus is a blank chromosome, so this costs nothing
observable, and the alternative is a new false string on screen. The *decision* is copied
exactly: whole-column, all-or-nothing, fired only when nothing in the column already says
``chr``.
"""

from __future__ import annotations

import pandas as pd

from config.missing_values import says_nothing_over

#: The column this module is about, named for the two reads below and for the tests, which
#: would otherwise spell it as a literal in twenty places.
#:
#: Deliberately **not** offered as the app's one name for this column. Three other modules
#: spell it inline — ``components/variant_table._KEY_COLUMNS``, ``components/charts.py``
#: and ``components/variant_detail.py`` — and importing this constant into them would say
#: they depend on this module, which they do not: they read a column, they do not decide
#: its spelling. A single name for every reader of ``Chromosome`` is a different change
#: from issue #149's, and it is not one this module should make by the back door.
CHROMOSOME = "Chromosome"


def normalise_chromosome_spelling(maf: pd.DataFrame) -> bool:
    """Give this MAF's chromosomes the ``chr`` prefix, if none of them has it already.

    Args:
        maf: a freshly read MAF. Modified in place — it is not shared and nothing caches
            it, and copying a 300-column frame to change one column is not free.

    Returns:
        bool: whether anything was re-spelled. ``page_modules/data_loading`` uses this to
        decide whether to tell the user, since the change reaches their download.

    Does nothing at all when the column is absent (``validate_required_columns`` refuses
    such a file a moment later, and refusing is its job, not this one's), when any value
    already says ``chr``, or when every value is blank.
    """
    if CHROMOSOME not in maf.columns:
        return False

    column = maf[CHROMOSOME]
    rendered = column.astype(str)

    # The pipeline's own test, character for character: substring rather than prefix,
    # case-sensitive, and over the *whole* column rather than per row. Copying it means
    # the app cannot reach a different answer than the file's annotator did.
    if rendered.str.contains("chr", regex=False).any():
        return False

    says_something = ~says_nothing_over(column)
    if not says_something.any():
        return False

    # Built as an object column and assigned wholesale rather than written through
    # ``.loc``: a bare all-numeric chromosome column reads as ``int64``, and writing
    # strings into a subset of one leaves pandas to upcast mid-assignment. Building the
    # replacement keeps every blank cell exactly the missing value it was — assigning
    # ``rendered`` wholesale would have put the string ``"nan"`` in the export.
    replacement = column.astype(object).copy()
    replacement[says_something] = "chr" + rendered[says_something]
    maf[CHROMOSOME] = replacement
    return True
