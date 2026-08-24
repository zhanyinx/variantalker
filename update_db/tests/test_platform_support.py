"""The manual route is Linux-only, and `update_db/README.md` says so with numbers.

Issue #370 reported that macOS `zcat` cannot open a `.gz` at all -- it is the `compress`-era
tool, which appends `.Z` to the name it is given. The ruling on that ticket was **not** to make
the spellings portable: the databases these scripts maintain need several hundred gigabytes, so a
laptop is ruled out by disk whatever its userland does, and a portable `zcat` would buy a route
nobody could use. Instead `README.md` gained a `Platform support` section that declares the manual
route Linux-only and states, precisely, which GNU spellings break and how many call sites each has.

That declaration is the thing this file guards. A README sentence with a number in it is exactly
what rots in this repository -- so if someone changes the scripts, the counts here fail and name
the README as the file to update. Two deliberate choices in how that is asserted:

**Counts and files, never line numbers.** This was vindicated immediately: issue #354's
`set -eo pipefail` header landed on the Annovar route while this file was being written, shifting
every `zcat` line in `update_clinvar_annovar.sh` from 98/148 to 111/161. The counts did not move.
A guard that pinned line numbers would have failed on a change with nothing to do with portability.

**Comment lines are excluded, and that is load-bearing rather than tidiness.** Issue #353 added
header comments to the Funcotator scripts that discuss `zcat` by name, so a naive scan counts
**ten** `zcat` sites where there are eight. Getting that wrong would have put a wrong number in
the README and then guarded it.

The reverse trap -- a guard whose own prose is counted by its own scan -- is the one
`test_shell_interpreters.py` records paying for. It does not bite here because the scan reads only
`*.sh` files while this file is Python, but the walk still excludes `tests/` explicitly, and
`test_the_scan_is_not_hiding_shell_scripts_in_its_own_tree` proves that exclusion hides nothing.
"""

from __future__ import annotations

import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
UPDATE_DB = REPO_ROOT / "update_db"
README = UPDATE_DB / "README.md"

# The distribution `README.md`'s `Platform support` section states in prose. Keyed by the path the
# README names, so a failure message points at the same words a reader would go and fix.
EXPECTED_ZCAT = {
    "annovar4cancervar_intervar/update_clinvar_annovar.sh": 2,
    "funcotator/update_clinvar/update_clinvar_funcotator.sh": 2,
    "funcotator/update_dbsnp/update_dbsnp.sh": 4,
}
EXPECTED_SED_I = {
    "funcotator/update_clinvar/update_clinvar_funcotator.sh": 4,
}

# The floor the README quotes, in GB. Stated as "at least 500 GB" in two files.
EXPECTED_DISK_FLOOR_GB = 500


def _shell_scripts() -> list[Path]:
    """Every shell script the manual route runs, excluding this guard's own directory."""
    return sorted(p for p in UPDATE_DB.rglob("*.sh") if "tests" not in p.relative_to(UPDATE_DB).parts)


def _call_sites(pattern: str) -> dict[str, int]:
    """Count non-comment lines matching `pattern`, per script, relative to `update_db/`."""
    found: dict[str, int] = {}
    regex = re.compile(pattern)
    for script in _shell_scripts():
        hits = sum(
            1
            for line in script.read_text().splitlines()
            # A `#` in column one (after indentation) is a comment. #353's headers discuss
            # `zcat` at length, and counting them would inflate eight sites to ten.
            if not line.lstrip().startswith("#") and regex.search(line)
        )
        if hits:
            found[script.relative_to(UPDATE_DB).as_posix()] = hits
    return found


# ---------------------------------------------------------------------------------------------
# The properties.
# ---------------------------------------------------------------------------------------------


def test_the_zcat_call_sites_are_where_the_readme_says():
    """Eight `zcat` sites, across the two ClinVar scripts and `update_dbsnp.sh`."""
    assert _call_sites(r"\bzcat\b") == EXPECTED_ZCAT, (
        "the `zcat` call sites have moved. update_db/README.md's `Platform support` section "
        "states this distribution in prose -- if the spellings were made portable, or a new one "
        "was added, that section is now wrong and is the thing to fix."
    )


def test_the_sed_in_place_call_sites_are_where_the_readme_says():
    """Four `sed -i` sites, all on the Funcotator ClinVar script.

    These are the dangerous half. BSD `sed` reads the expression as `-i`'s backup-suffix
    argument, exits 1 and leaves the file unchanged -- and two of the four are
    `s/DATE/$date/g` on `clinvar_vcf.config`, the file that *names* the database. That is
    issue #354's stale-data-reported-as-fresh shape, which is why the count is pinned.
    """
    assert _call_sites(r"\bsed\s+-i\b") == EXPECTED_SED_I, (
        "the `sed -i` call sites have moved. update_db/README.md's `Platform support` section "
        "names this file and this count; update it, or -- if these became portable -- say so."
    )


def test_the_readme_declares_the_manual_route_unsupported_on_macos():
    """The declaration itself: the section exists, and macOS is refused rather than warned about."""
    text = README.read_text()

    assert "## Where this runs" in text, "the disk/HPC requirement section has been removed"
    assert "### Platform support" in text, (
        "the `Platform support` section has been removed, but the scripts still use GNU-only "
        "spellings -- a reader on macOS now gets no warning at all"
    )

    macos_row = next((ln for ln in text.splitlines() if ln.startswith("| macOS")), None)
    assert macos_row is not None, "the platform table no longer has a macOS row"
    assert "**Not supported**" in macos_row, (
        f"the macOS row no longer refuses the manual route: {macos_row!r}. Softening this to a "
        "warning is the failure mode -- the route does not work there, it is not merely untested."
    )


def test_both_readmes_quote_the_same_disk_floor():
    """The disk requirement is the reason for the ruling, so both routes' docs must state it."""
    for name in ("README.md", "README_DOCKER.md"):
        text = (UPDATE_DB / name).read_text()
        assert f"{EXPECTED_DISK_FLOOR_GB} GB" in text, (
            f"update_db/{name} no longer states the {EXPECTED_DISK_FLOOR_GB} GB floor. It is the "
            "whole reason the manual route is Linux-only rather than merely buggy on macOS, and "
            "the container route needs it just as much -- the image supplies software, not disk."
        )


def test_the_readmes_account_of_the_failure_being_loud_is_still_true():
    """`Platform support` claims the failure is loud on **both** routes. It says so because of
    `set -eo pipefail`, so the claim is only as true as those headers.

    This property has already earned its keep, and the story is worth keeping. Its first version
    asserted the *asymmetry* that held while it was being written -- `pipefail` on the Funcotator
    route (#353) and not yet on the Annovar one (#354) -- and said in its own failure message that
    failing would be good news. Issue #354 then merged between the branch being pushed and CI
    running it, the property failed on the first CI run, and the message named the README paragraph
    to delete. Which is what happened. A doc claim that was true for about five minutes was caught
    by the guard written to protect it rather than by a reader years later.

    So the claim is now the stronger one -- both routes, no asymmetry -- and what would break it is
    a `pipefail` header being *removed*.
    """
    routes = {
        "Funcotator ClinVar (#353)": "funcotator/update_clinvar/update_clinvar_funcotator.sh",
        "Funcotator dbSNP (#353)": "funcotator/update_dbsnp/update_dbsnp.sh",
        "Annovar ClinVar (#354)": "annovar4cancervar_intervar/update_clinvar_annovar.sh",
    }
    missing = [
        label
        for label, relpath in routes.items()
        if "set -eo pipefail" not in (UPDATE_DB / relpath).read_text()
    ]
    assert not missing, (
        f"these scripts have lost their `set -eo pipefail` header: {missing}. "
        "update_db/README.md's `Platform support` section says the failure is loud on both routes "
        "and that a Mac run is wasted rather than corrupting -- without `pipefail` that is false "
        "again, because a pipeline takes its exit status from its last command and the failed "
        "`zcat` goes unnoticed."
    )


def test_the_scan_is_not_hiding_shell_scripts_in_its_own_tree():
    """The `tests/` exclusion must hide nothing, or the counts above are quietly partial.

    `test_shell_interpreters.py` records what this costs: an `rglob` guard living inside the tree
    it audits passed while it sat under `docs/` and failed once it shipped. Nothing here relies on
    that going well -- it is asserted.
    """
    assert not list((UPDATE_DB / "tests").rglob("*.sh")), (
        "a shell script has appeared under update_db/tests/, which the scan in this file skips. "
        "Either move it or widen the walk -- as written, its call sites are invisible."
    )
    # And the scan must actually be reaching the tree: a walk that found nothing would make
    # every count above vacuously equal to an empty dict.
    #
    # Named scripts rather than a total, and that is the whole point (issue 407). This was
    # `== 9`, which is a statement about one checkout rather than about the repository: the
    # export strips the gnomAD preprocessing directory, so the published tree holds eight
    # and the anchor failed there — red CI on the repository the release is cut from, on a
    # claim that was never about anything a reader cares about. The deny-list cannot be
    # consulted to excuse it either, since `.publicignore` does not travel: only a claim true
    # in **both** trees can hold, so the anchor names the scripts that are always present.
    #
    # Third time this shape has bitten — see issues 383 and 375, both anchors written as
    # totals, both false in the public clone. Assert files and rules, never a total.
    found = {path.name for path in _shell_scripts()}
    always_present = {
        "update_clinvar_annovar.sh",
        "update_cosmic.sh",
        "update_clinvar_funcotator.sh",
        "createSqliteCosmicDb.sh",
        "updateCosmicDataBase.sh",
        "update_dbsnp.sh",
        "getGencode.sh",
        "get_new_hgnc.sh",
    }
    missing = sorted(always_present - found)
    assert not missing, (
        f"the scan is not reaching {missing}. Every script named here ships in both the private "
        "tree and the published one, so this is a walk that has stopped working rather than a "
        "tree that has changed — and without it every count in this module goes vacuously true."
    )
