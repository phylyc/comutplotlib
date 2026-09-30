# ComutPlotLib - Genomic Comutation Plots
Welcome to **ComutPlotLib**, a Python tool for generating **genomic comutation plots**. These plots visualize **co-occurrence patterns** of genomic alterations—**mutations, copy number variations (CNVs), and other genetic events**—across multiple patients or samples.  

## **Example Output**  
![Sample Comut Plot Output](https://raw.githubusercontent.com/phylyc/comutplotlib/main/demo/comut_test.png)  
This comutation plot (central panel) visualizes the mutation landscape across a patient cohort. Rows represent genes, and columns correspond to patients. Each cell indicates a gene’s mutation status in a patient: rectangles denote copy-number variations (CNVs), and ellipses indicate short nucleotide variations (SNVs), with multiple SNVs shown as subdivided wedges. Colors encode mutation types and functional effects.

The top panel displays tumor mutation burden (TMB) per patient, with high TMB (≥10/Mb) highlighted in red. The mutational signature panel shows the relative fraction of exposures to different mutational signatures for each patient or sample. The right panel summarizes mutation recurrence, showing SNV and/or CNV frequencies per gene, reporting the percentage of patients with high-level CNVs or at least low-level CNVs.

The bottom panel presents patient- and sample-level metadata. For patients with multiple samples, metadata cells are subdivided accordingly.

The library also supports case-control plots, featuring fold-change of the mutational frequency between two cohorts:

![Sample Comut Plot Output](https://raw.githubusercontent.com/phylyc/comutplotlib/main/demo/comut_test.control.png)




## **Features**  
✔ **Visualizes SNVs and CNVs in a single plot**  
✔ **Summarizes mutational burden, recurrence, and metadata**  
✔ **Customizable layout and annotation**  
✔ **Compatible with MAF, GISTIC, and SIF files**  
✔ **Integrates with GATK Funcotator, GISTIC 2.0, and other genomic tools**  


---

## **Installation**  
To install **ComutPlotLib**, run:  
```bash
pip install -r requirements.txt
```
Alternatively, use the provided installation script:
```
bash install.sh
```

## **Usage**

Run **ComutPlotLib** with:
```
python comut_argparse.py --output output_plot.png --maf input.maf
```
For a full list of options, use:
```
python comut_argparse.py --help
```
📌 For detailed examples, refer to the ![demo folder](https://github.com/phylyc/comutplotlib/tree/main/demo).


### **Mutational signatures**

Signature exposures are ingested with `--signatures` (and `--control-signatures`):
a table with one row per patient/sample and one column per signature. Columns are
coloured and ordered by the etiology they belong to (`clock-like`, `APOBEC`,
`MMR`, `HRD`, `PolE/D/H`, `Smoking`, `UV`, `Treatment`, `Other`, `Error`,
`Unknown`).

Add `--group-signatures-by-etiology` to **sum all signatures of the same etiology
into a single stacked-bar category**, which is useful when many individual
signatures make the panel hard to read:

```bash
python call.py -o out.png --maf input.maf \
  --signatures exposures.tsv \
  --group-signatures-by-etiology
# SBS1, SBS5 -> clock-like;  SBS2, SBS13 -> APOBEC;  SBS4, SBS29, SBS92 -> Smoking; ...
```

Signatures that are not part of any known etiology keep their own category and are
sorted to the end. The flag applies to the case and the control cohort alike.


### **Grid sub-panels**

Large cohorts can be split into a 2-D **grid of smaller comut sub-panels**,
stratified by sample metadata (grid columns) and/or explicit gene sets (grid
rows). Row (gene) and column (sample) orderings stay globally consistent, and all
shared scales and legends are global for cross-stratum comparability.

```
python comut_argparse.py --output grid.png --maf input.maf --sif input.sif \
  --column-group-by "Sample Type" \
  --gene-groups "SNV genes:EGFR,KRAS,TP53;CNV genes:MYC,CDKN2A"
```

Relevant options:
- `--column-group-by` — comma-separated metadata keys to stratify samples into
  grid columns (auto-added to `--meta-data-rows`).
- `--hide-grouped-meta-data` — hide rows named by `--column-group-by` from the
  metadata table after they have been used to create the grid columns.
- `--gene-groups` — explicit gene groups forming grid rows
  (`GroupA:GeneA,GeneB;GroupB:GeneC`). The **order of the groups determines the
  top-to-bottom row order**. Ungrouped genes go to a trailing `Other` row (unless
  `--drop-ungrouped-genes`).
- `--group-order` — highest-priority explicit ordering of group keys (then by
  group size, then key tuple). See the note below.
- `--column-group-labels` — space-separated display labels for the grid columns,
  applied in grid-column order (see below).
- `--na-group-label` / `--other-gene-group-label` — labels for missing-value and
  ungrouped groups.

**Using `--group-order`.** It is a **flat, comma-separated list of individual
group *values*** (not a dict, *not* the `--column-group-by` column names, and
*not* joined keys). A value is:

- for **column groups**: a metadata *value* of a stratifying column
  (e.g. `neg`, `pos`). Missing / `nan` / `unknown` values fall under the
  `--na-group-label` (default `NA`);
- for **gene groups**: the group *name* from `--gene-groups`
  (e.g. `SNV genes`, or the `--other-gene-group-label`, default `Other`).

Grid columns are **nested in `--column-group-by` order**: the first key is the
outermost split, the second nests inside it, and so on. `--group-order` is
applied at *every* level; values you list come first (in that order), then any
remaining values by group size (descending) then name. So the order of
`--column-group-by` determines the sort *hierarchy*, while `--group-order`
determines the value order within each level.

Examples — a single status column with values `pos`, `neg`, and `unknown`
(→ `NA`), sorted `neg` first:

```
--column-group-by "HER2 Status inferred" \
--group-order "neg,pos"          # => columns: neg, pos, NA
```

Two status columns — `HR Status inferred` is the outer split (listed first),
`HER2 Status inferred` nests inside it; `neg` before `pos` at both levels:

```
--column-group-by "HR Status inferred,HER2 Status inferred" \
--group-order "neg,pos"
# => neg|neg, neg|pos, neg|NA, pos|neg, pos|pos, pos|NA, NA|...
```

Swap the two `--column-group-by` entries to nest by HER2 first instead. You only
list the individual values once (`neg,pos`); they apply at every level.

The grid can also draw **group titles**: add `column group label` and/or
`gene group label` to `--panels-to-plot` to render the column-group headers
(above the top marginals) and gene-group labels (in the far-left gutter). They are
off by default.

**Naming the grid columns.** By default a column-group header shows the raw
metadata value(s) (multiple `--column-group-by` keys are joined by ` | `). Use
`--column-group-labels` to give them readable titles. It takes **space-separated
labels** (quote labels containing spaces) applied **positionally in grid-column
order**, i.e. the order produced by `--group-order`:

```
--column-group-by "Sample Type" \
--group-order "BM,EM" \
--column-group-labels "Brain metastases" "Extracranial metastases"
```

Passing `--column-group-labels` implies adding the `column group label` panel to
`--panels-to-plot`. Fewer labels than grid columns leaves the trailing columns
with their auto-generated titles; pass `""` to keep an individual one. With a
control cohort, both cohorts label the same group identically, even if the
control resolves the groups in a different order.

Omitting these options reproduces the standard single-panel figure exactly.


### **Input Files**:
ComutPlotLib requires at least one of the following input files:
1. **Mutation Annotation Format (MAF)**: 
   - Output of GATK Funcotator
   - Contains mutation calls and annotations
2. **GISTIC output**: 
   - From GISTIC 2.0
   - Provides copy number alteration calls (file: all_thresholded.by_gene.txt)

Sample information can be provided via
- **Sample Information File (SIF)**:
  - Tab-separated metadata file
  - Contains sample attributes (e.g., tumor purity, platform, histology)
  - See ![sample_annotation.py](https://raw.githubusercontent.com/phylyc/comutplotlib/main/comutplotlib/sample_annotation.py) for required columns


## **Demo & Examples**:

🔬 Try the demo:

- Generate your own synthetic data: 
```
python make_data.py
```
- Generate a plot
```
bash call_comut.sh
```
This will generate an example comutation plot using the synthetic test datasets.


## **Dependencies**: 
Required dependencies are listed in requirements.txt. Install them via:
```
pip install -r requirements.txt
```


## **Author**:
Developed by Philipp Hähnel.


## **License**  
This project is licensed under the **MIT License**. See the [LICENSE](./LICENSE) file for details.  
