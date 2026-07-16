# comutplotlib configuration

This folder holds **institution / study-specific vocabulary** as editable data
files, so that adapting comutplotlib to a new cohort does not require changing
Python source. Everything here is a **template**: it ships with the sensible
defaults the tool was originally developed against, and you are expected to
copy and adapt it.

## Files

| File | Consumed by | Purpose |
|---|---|---|
| `gene_aliases.json` | `comutplotlib/gistic.py` | Maps deprecated/aliased gene symbols to current HGNC-approved symbols for GISTIC indices. |
| `meta_data_rows.json` | `comutplotlib/comut_argparse.py` | Default `--meta-data-rows` and `--meta-data-rows-per-sample`. |
| `sif_classification.json` | `comutplotlib/sample_classification.py` (via `sif.py`) | Ordered rules that harmonize free-text SIF fields (cancer type, platform, sample type, material) into short category codes. |

## Overriding the defaults

Two mechanisms, in order of precedence:

1. **A private config directory.** Set the environment variable
   `COMUTPLOTLIB_CONFIG_DIR` to a folder containing files with the **same
   names** as above. Any file present there is loaded instead of the shipped
   template; missing files fall back to the template. This is the recommended
   way to share a config with internal collaborators without forking the code:

   ```bash
   export COMUTPLOTLIB_CONFIG_DIR=/path/to/your/comutplotlib_config
   ```

2. **CLI flags** for the pieces that expose them, e.g. `--meta-data-rows`,
   `--meta-data-rows-per-sample`.

## Rule format (`sif_classification.json`)

Each classifier is an ordered list of rules; the **first match wins**.

| `op` | Meaning |
|---|---|
| `contains` | any string in `values` is a substring of `field` |
| `equals` | `field == value`, or `field in values` |
| `na_or_in` | `field` is a float/NaN, or `field in values` |
| `all` / `any` | combine nested `conditions` |
| `default` | always matches (place last) |

`return` is either a literal code (e.g. `"BRCA"`) or `"@field"` to return that
field's raw value. `transform: "replace_space_dash"` replaces spaces with
dashes in the returned value.

## Provenance notes

- `gene_aliases.json` is derived from the HGNC approved-symbol database
  (<https://www.genenames.org>). **Please record the HGNC release** when you
  update it.
- `sif_classification.json` was extracted verbatim from the previously
  hard-coded logic in `sif.py`; it encodes conventions from the original
  development cohorts (neuro-oncology / breast) and should be reviewed before
  use on other cohorts.

