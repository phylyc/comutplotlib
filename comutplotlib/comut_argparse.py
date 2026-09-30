import argparse
from collections import defaultdict

from comutplotlib.config import load_config
from comutplotlib.comut_panels import DEFAULT_PANELS, ComutPanels, control

# Default metadata rows are shipped as editable config (see comutplotlib/config/).
_META_ROWS = load_config("meta_data_rows.json", default={"rows": [], "rows_per_sample": []})


def validate_args(args):
    """Ensure required arguments are correctly specified."""
    def remove(panel):
        if panel in args.panels_to_plot:
            args.panels_to_plot.remove(panel)

    if args.maf is None and args.gistic is None:
        raise ValueError("Either --maf or --gistic must be specified.")

    if args.control_maf is None and args.control_gistic is None:
        remove(ComutPanels.recurrence_fold_change)
        for panel in args.panels_to_plot:
            if panel.endswith("control"):
                remove(panel)
    if ComutPanels.recurrence not in args.panels_to_plot:
        remove(ComutPanels.total_recurrence_overall)
    if control(ComutPanels.recurrence) not in args.panels_to_plot:
        remove(control(ComutPanels.total_recurrence_overall))
    if ComutPanels.recurrence_fold_change not in args.panels_to_plot:
        remove(ComutPanels.total_recurrence_fold_change)

    if args.cohort_label is None:
        remove(ComutPanels.cohort_label)
    if args.control_cohort_label is None:
        remove(control(ComutPanels.cohort_label))

    if args.snv_interesting_genes is None and args.cnv_interesting_genes is None:
        remove(ComutPanels.model_annotation)
    if ComutPanels.model_annotation not in args.panels_to_plot:
        remove(ComutPanels.model_annotation_legend)

    if args.signatures is None:
        remove(ComutPanels.mutational_signatures)
    if ComutPanels.mutational_signatures not in args.panels_to_plot:
        remove(ComutPanels.mutational_signatures_legend)

    if args.sif is None:
        remove(ComutPanels.meta_data)
    if ComutPanels.meta_data not in args.panels_to_plot:
        remove(ComutPanels.meta_data_legend)

    if args.gene_meta_data is None:
        remove(ComutPanels.gene_meta_data)

    # Grid stratification keys must be present as metadata rows so they are loaded.
    if getattr(args, "column_group_by", None):
        meta_rows = list(args.meta_data_rows) if args.meta_data_rows is not None else []
        for key in args.column_group_by:
            if key not in meta_rows:
                meta_rows.append(key)
        args.meta_data_rows = meta_rows

    # Explicit column-group labels are only visible when the corresponding panel
    # is drawn, so requesting labels implies requesting the panel.
    if getattr(args, "column_group_labels", None):
        if not getattr(args, "column_group_by", None):
            raise ValueError("--column-group-labels requires --column-group-by.")
        if ComutPanels.column_group_label not in args.panels_to_plot:
            args.panels_to_plot.append(ComutPanels.column_group_label)


def print_args(args):
    if args.verbose:
        print("Calling ComutPlotLib")
        print("Arguments:")
        for key, value in vars(args).items():
            print(f"  {key}: {value}")
        print()


def parse_comma_separated(value):
    """Convert a comma-separated string into a list."""
    return value.split(",") if value else None


def parse_palette(value):
    """ Parse a color palette definition from a semicolon-separated string. Format:
        value='key1:r,g,b;key2:r,g,b|metaColumn1>key3:r,g,b;key4:r,g,b|metaColumn2>key5:r,g,b'
        returns {'global': {key1: (r, g, b), key2: {r, g, b)}, 'local': {metaColumn1: {key3: (r, g, b), key4: (r, g, b)}, metaColumn2: {key5: (r, g, b)}}}
    """
    if not value:
        return None
    palette = {}
    meta_palette = {}
    palettes = value.split("|")
    if len(palettes[0]) > 0:
        for entry in palettes[0].split(";"):
            key, rgb = entry.split(":")
            palette[key] = tuple(map(float, rgb.split(",")))
    if len(palettes) > 1:
        for meta_entries in palettes[1:]:
            column, meta_hash = meta_entries.split(">")
            meta_palette[column] = {}
            for entry in meta_hash.split(";"):
                key, rgb = entry.split(":")
                meta_palette[column][key] = tuple(map(float, rgb.split(",")))
    return {"global": palette, "local": meta_palette}


def parse_gene_categories(value):
    """ Parse a color palette definition from a semicolon-separated string. Format:
        value='x,y,z|gene1:x,y,z;gene2:x,y,z'
        returns {'global': [x, y,  z], 'local': {gene1: [x, y, z], gene2: [x, y, z]}
    """
    if not value:
        return None
    global_categories = []
    gene_categories = defaultdict(list)
    cats = value.split("|")
    if len(cats[0]) > 0:
        for cat in cats[0].split(","):
            global_categories.append(cat)
    if len(cats) > 1:
        for gene_entry in cats[1].split(";"):
            gene, gene_cats = gene_entry.split(":")
            for cat in gene_cats.split(","):
                gene_categories[gene].append(cat)
    return {"global": global_categories, "local": gene_categories}


def parse_ground_truth_genes(value):
    """Parse ground truth genes into a dictionary format."""
    if not value:
        return None
    return {k: v.split(",") for k, v in (x.split(":") for x in value.split(";"))}


def parse_gene_groups(value):
    """Parse gene groups for grid rows: 'GroupA:GeneA,GeneB;GroupB:GeneC'."""
    if not value:
        return None
    return {k: v.split(",") for k, v in (x.split(":") for x in value.split(";"))}



def parse_args():
    """Parse command-line arguments for ComutPlot."""
    parser = argparse.ArgumentParser(
        prog="ComutPlot",
        description="Generate genomic comutation plots.",
        epilog="For more details, refer to the documentation.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "-o", "--output", type=str, required=True,
        help="Path to the output file."
    )

    parser.add_argument(
        "--maf", type=str, action='append', default=None,
        help="Path to a MAF file (output from GATK Funcotate). Can be specified multiple times."
    )
    parser.add_argument(
        "--sif", type=str, action='append', default=None,
        help="Path to a SIF input file. See documentation for format details."
    )
    parser.add_argument(
        "--gistic", type=str, action='append', default=None,
        help="Path to a GISTIC output file (e.g., 'all_thresholded.by_gene.txt')."
    )
    parser.add_argument(
        "--signatures", type=str, action='append', default=None,
        help="Path to a file containing mutational signature exposures (index: patient, columns: signatures)."
    )
    parser.add_argument(
        "--group-signatures-by-etiology", default=False, action=argparse.BooleanOptionalAction,
        help="Aggregate mutational signature exposures by etiology (e.g. 'clock-like', "
             "'APOBEC', 'MMR', ...), summing all signatures that belong to the same "
             "etiology into a single stacked-bar category. Signatures that are not part "
             "of any known etiology keep their own category. Applies to both --signatures "
             "and --control-signatures."
    )
    parser.add_argument(
        "--cohort-label", type=str, default=None,
        help="Name of the cohort for panel title."
    )
    parser.add_argument(
        "--cohort-tag", type=str, default=None,
        help="Short tag for the cohort, used for in-plot annotations such as the "
             "fold-change panel arrows. Defaults to 'case'."
    )

    parser.add_argument(
        "--control-maf", type=str, action='append', default=None,
        help="Path to a MAF file (output from GATK Funcotate). Can be specified multiple times."
    )
    parser.add_argument(
        "--control-sif", type=str, action='append', default=None,
        help="Path to a SIF input file. See documentation for format details."
    )
    parser.add_argument(
        "--control-gistic", type=str, action='append', default=None,
        help="Path to a GISTIC output file (e.g., 'all_thresholded.by_gene.txt')."
    )
    parser.add_argument(
        "--control-signatures", type=str, action='append', default=None,
        help="Path to a file containing mutational signature exposures (index: patient, columns: signatures)."
    )
    parser.add_argument(
        "--control-cohort-label", type=str, default=None,
        help="Name of the control cohort for panel title."
    )
    parser.add_argument(
        "--control-cohort-tag", type=str, default=None,
        help="Short tag for the control cohort, used for in-plot annotations such "
             "as the fold-change panel arrows. Defaults to 'control'."
    )
    parser.add_argument(
        "--control-position", type=str, choices=["left", "right"], default="right",
        help="Side on which to place the control cohort (only relevant when a "
             "control cohort is supplied). 'right' (default) keeps the current "
             "layout; 'left' swaps the case and control cohorts, including the "
             "case/control fold-change arrows."
    )

    parser.add_argument(
        "--gene-meta-data", type=str, action='append', default=None,
        help="Path to a file containing a table with binary gene annotations, one gene per row."
    )
    parser.add_argument(
        "--gene-meta-columns", type=parse_comma_separated, default=None,
        help="Comma-separated names of columns to plot for gene meta data."
    )

    parser.add_argument(
        "--model-names", type=parse_comma_separated, default=["CNV burden model", "SNV burden model"],
        help="Comma-separated names of models corresponding to SNV and CNV interesting genes."
    )
    parser.add_argument(
        "--model-significances", type=str, action="append", default=None,
        help="Path to files containing gene significance values for each model."
    )

    parser.add_argument(
        "--by", type=str, choices=["Sample", "Patient"], default="Patient",
        help="Aggregation level for mutations (Sample or Patient)."
    )
    parser.add_argument(
        "--label-columns", action="store_true",
        help="Display labels at the bottom of the plot."
    )
    parser.add_argument(
        "--column-order", type=parse_comma_separated, default=None,
        help="Comma-separated order of columns in the plot."
    )
    parser.add_argument(
        "--control-column-order", type=parse_comma_separated, default=None,
        help="Comma-separated order of columns in the plot."
    )
    parser.add_argument(
        "--index-order", type=parse_comma_separated, default=None,
        help="Comma-separated order of gene indices in the plot."
    )
    parser.add_argument(
        "--column-sort-by", type=parse_comma_separated, default=["COMUT"],
        help="Comma-separated order of 'COMUT', 'TMB', and meta data rows to sort by, in order."
    )
    parser.add_argument(
        "--index-sort-by", type=str, default="COMUT",
        choices=["COMUT", "fold-change", "high-fold-change", "low-fold-change"],
        help="How to sort the gene (row) index. 'COMUT' is the standard mutation-based "
             "ordering. The 'fold-change' modes require a control cohort and sort genes by "
             "the case-vs-control fold change: 'fold-change' uses the total fold change across "
             "all categories; 'high-fold-change' uses the selected high fold change (high/mid "
             "CNV or SNV); 'low-fold-change' uses the selected low fold change (low CNV). "
             "For the 'high'/'low' modes only the categories selected by --recurrence-categories "
             "for each gene are considered."
    )
    parser.add_argument(
        "--index-sort-by-direction", type=str, default="ascending",
        choices=["ascending", "descending"],
        help="Direction of the --index-sort-by ordering (default: descending)."
    )
    parser.add_argument(
        "--drop-empty-columns", action="store_true",
        help="Exclude samples/patients with no mutation data."
    )

    parser.add_argument(
        "--meta-data-rows", type=parse_comma_separated,
        default=list(_META_ROWS.get("rows", [])),
        help="Comma-separated list of metadata fields to display.")
    parser.add_argument(
        "--meta-data-rows-per-sample", type=parse_comma_separated,
        default=list(_META_ROWS.get("rows_per_sample", [])),
        help="Metadata fields to display per sample."
    )
    parser.add_argument(
        "--meta-data-position", type=str, choices=["top", "bottom"], default="bottom",
        help="Place the sample metadata panel below the comutation plot, or above it below mutational signatures."
    )

    parser.add_argument(
        "--interesting-gene", type=str, default=None,
        help="Find genes co-mutated at least a given percentage (set by --interesting-gene-comut-percent-threshold)."
    )
    parser.add_argument(
        "--interesting-gene-comut-percent-threshold", type=float, default=None,
        help="Minimum co-mutation percentage threshold for interesting genes."
    )
    parser.add_argument(
        "--interesting-genes", type=parse_comma_separated, default=None,
        help="Comma-separated list of interesting genes."
    )
    parser.add_argument(
        "--snv-interesting-genes", type=parse_comma_separated, default=None,
        help="Comma-separated list of SNV interesting genes."
    )
    parser.add_argument(
        "--cnv-interesting-genes", type=parse_comma_separated, default=None,
        help="Comma-separated list of CNV interesting genes."
    )
    parser.add_argument(
        "--scale-recurrence", action="store_true",
        help=" "
    )
    parser.add_argument(
        "--recurrence-categories", type=parse_gene_categories, default={"global": ["snv", "amp", "del"], "local": {}},
        help="List of mutation categories (snv/amp/del) to plot recurrence for; "
             "optional gene level categories can be added via the format x,y,z|gene1:x,y,z;gene2:x,y,z"
    )
    parser.add_argument(
        "--total-recurrence-threshold", type=float, default=None,
        help="Minimum percentage of patients with a mutation in a gene for it to be plotted."
    )
    parser.add_argument(
        "--snv-recurrence-threshold", type=int, default=5,
        help="Minimum number of patients with the same protein change for annotation."
    )
    parser.add_argument(
        "--min-fold-change", type=float, default=None,
        help="Minimum total fold-change of fraction of patients carrying a mutation in a gene for it to be plotted."
    )
    parser.add_argument(
        "--max-fold-change", type=float, default=None,
        help="Maximum total fold-change of fraction of patients carrying a mutation in a gene for it to be plotted."
    )
    parser.add_argument(
        "--ground-truth-genes", type=parse_ground_truth_genes, default=None,
        help="Highlight genes with specified colors. Format: 'color1:gene1,gene2,...;color2:gene3,gene4,...'."
    )
    parser.add_argument(
        "--collapse-cytobands", default=False, action="store_true",
        help="Collapse genes in the same cytoband to one gene label. Keeps genes in --ground-truth-genes uncollapsed."
    )
    parser.add_argument(
        "--collapse-cytobands-range", type=int, default=0,
        help="Range for collapsing also neighboring cytobands."
    )

    parser.add_argument(
        "--low-amp-threshold", type=float, default=1,
        help="Threshold for low-level amplification."
    )
    parser.add_argument(
        "--mid-amp-threshold", type=float, default=1.5,
        help="Threshold for mid-level amplification."
    )
    parser.add_argument(
        "--high-amp-threshold", type=float, default=2,
        help="Threshold for high-level amplification."
    )
    parser.add_argument(
        "--baseline", type=float, default=0,
        help="Copy-number value representing the neutral baseline."
    )
    parser.add_argument(
        "--low-del-threshold", type=float, default=-1,
        help="Threshold for low-level deletion."
    )
    parser.add_argument(
        "--mid-del-threshold", type=float, default=-1.5,
        help="Threshold for mid-level deletion."
    )
    parser.add_argument(
        "--high-del-threshold", type=float, default=-2,
        help="Threshold for high-level deletion."
    )
    parser.add_argument(
        "--show-low-level-cnvs", default=True, action=argparse.BooleanOptionalAction,
        help="Include low-level CNVs in the plot (use --no-show-low-level-cnvs to hide them)."
    )

    parser.add_argument(
        "--panels-to-plot", type=parse_comma_separated,
        default=list(DEFAULT_PANELS),
        help="Comma-separated list of plot panels to include.")
    parser.add_argument(
        "--palette", type=parse_palette, default=None,
        help="Define additional colors for categories. Global palette in the front, "
             "followed by meta data-specific palettes. Format: "
             "'key1:r,g,b;key2:r,g,b|metaColumn1>key3:r,g,b;key4:r,g,b|metaColumn2>key5:r,g,b'."
    )
    parser.add_argument(
        "--max-xfigsize", type=int, default=None,
        help="Maximum x-axis figure size; the comutation plot scales to fit."
    )
    parser.add_argument(
        "--max-xfigsize-scale", type=float, default=1,
        help="Parameter by which xfigsize scales with the number of columns."
    )

    # --- Grid sub-panels (see GRID_SUBPANELS_PLAN.md) ---
    parser.add_argument(
        "--column-group-by", type=parse_comma_separated, default=None,
        help="Comma-separated metadata keys to stratify samples into grid columns "
             "(e.g. 'Histology,Platform'). Keys are auto-added to --meta-data-rows."
    )
    parser.add_argument(
        "--hide-grouped-meta-data", action="store_true",
        help="Hide metadata rows used by --column-group-by from the metadata table."
    )
    parser.add_argument(
        "--gene-groups", type=parse_gene_groups, default=None,
        help="Explicit gene groups forming grid rows. Format: "
             "'GroupA:GeneA,GeneB;GroupB:GeneC'. The order of the groups here "
             "determines the top-to-bottom row order. Ungrouped genes go to a "
             "trailing 'Other' row."
    )
    parser.add_argument(
        "--group-order", type=parse_comma_separated, default=None,
        help="Comma-separated ordering of individual GROUP VALUES (not the "
             "--column-group-by column names, and not joined keys). For column "
             "groups these are the metadata VALUES (e.g. 'neg,pos'; missing/"
             "'unknown' values fall under --na-group-label); for gene groups they "
             "are the group NAMES from --gene-groups. The list is applied at every "
             "hierarchy level. Grid columns are nested in --column-group-by order "
             "(first key = outermost split); within each level, listed values come "
             "first (in this order), then the rest by group size (desc) then name. "
             "Example: --column-group-by 'HR Status,HER2 Status' "
             "--group-order 'neg,pos'."
    )
    parser.add_argument(
        "--drop-ungrouped-genes", default=False, action=argparse.BooleanOptionalAction,
        help="Drop genes not covered by --gene-groups instead of collecting them into 'Other'."
    )
    parser.add_argument(
        "--column-group-labels", type=str, nargs="+", default=None,
        help="Display labels for the grid columns created by --column-group-by, "
             "applied positionally in grid-column order (i.e. the order resulting "
             "from --group-order). Quote labels containing spaces, e.g. "
             "--column-group-labels 'HR+ / HER2-' 'Triple negative'. Fewer labels "
             "than grid columns leaves the trailing columns with their "
             "auto-generated labels; pass '' to keep an individual one. Implies "
             f"adding the '{ComutPanels.column_group_label}' panel to --panels-to-plot."
    )
    parser.add_argument(
        "--na-group-label", type=str, default="NA",
        help="Group label for samples with a missing stratification value."
    )
    parser.add_argument(
        "--other-gene-group-label", type=str, default="Other",
        help="Label for the trailing gene group of ungrouped genes."
    )

    args = parser.parse_args()
    validate_args(args)
    return args
