# Where the per-database Funcotator scripts came from

The scripts in this directory are derived from GATK's own Funcotator data-source scripts,
which live upstream at
[`scripts/funcotator/data_sources/`](https://github.com/broadinstitute/gatk/tree/master/scripts/funcotator/data_sources)
in the GATK repository. `update_gencode/getGencode.sh` and
`update_gencode/fixGencodeOrdering.py` still carry their upstream filenames; the rest were
renamed as they were adapted.

They were taken from an older revision and then modified locally — `getGencode.sh`, for
instance, gained a `--build` argument that upstream's copy does not have, and its
`LATEST_RELEASE` was moved on from upstream's value. They have **not** been kept in step
with upstream since, so the two have diverged in both directions.

## Why that matters

Upstream ships these scripts alongside an empty file named
`ALL_SCRIPTS_ARE_UNSUPPORTED.txt` — the name is the whole message. So when an external
source changes its download format, there is unlikely to be an upstream fix to pull in:
the COSMIC and GENCODE updates disabled in [../README.md](../README.md) would have to be
rewritten here.

Consult upstream when adapting one of these scripts to a source that has changed, but
expect to do the work rather than to find it done.
