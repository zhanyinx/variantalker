# Utils for variantalker

Bash scripts that **build the samplesheet** for a variantalker run, so that you point at a
folder of caller output rather than typing a csv by hand. Nothing here is part of the
pipeline: each script walks a directory tree and writes a csv to `--output` (default
`out.csv`), which you then pass to `nextflow run ... --input`.

The root [README](../README.md) does not mention this directory, so this file is the only
index of it.

| script | builds the samplesheet for | needs |
|---|---|---|
| `create_inputfile_variantalker.sh` | `--analysis annotation` — the [4-column input](../README.md#input-format), which is why it needs a tissue | `-i` a folder of per-sample folders holding `*hard-filtered*vcf.gz` and `*cnv.vcf.gz`, `-t` the sample tissue |
| `create_inputfile_biomarkers.sh` | `--analysis biomarkers` — the [3-column input](../docs/biomarkers/README.md#input) | `-i` a variantalker annotation folder (containing `germline/` and `somatic/`); optionally `-d` a Dragen folder with `*.tmb.metrics.csv` and `*.microsat_output.json`, and `-r` a Dragen RNA folder with `*.sf` |
| `create_inputfile_clonal_tmb_input.sh` | `--clonal_tmb_input` — the [clonal_evolution format](https://github.com/zhanyinx/clonal_evolution#input) | `-i` a folder of per-sample folders holding two crams, one named `*_tumor.cram`; `-a` a variantalker annotation output folder, e.g. `path2/somatic` |

Each script self-documents. Run it with no arguments for a usage line, or `-h` for the
full option list — that is the authority on its flags, not this table:

```bash
utils/create_inputfile_variantalker.sh -h
utils/create_inputfile_biomarkers.sh -h
utils/create_inputfile_clonal_tmb_input.sh -h
```

Both long and short option names work (`--input` / `-i`). Check the csv before you use it:
these scripts infer patient and sample names from directory names, and a tree laid out
differently from what the script expects produces a csv that is wrong rather than empty.
