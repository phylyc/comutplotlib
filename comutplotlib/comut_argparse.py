import argparse
from collections import defaultdict


def validate_args(args):
    """Ensure required arguments are correctly specified."""
    def remove(panel):
        if panel in args.panels_to_plot:
            args.panels_to_plot.remove(panel)

    if args.maf is None and args.gistic is None:
        raise ValueError("Either --maf or --gistic must be specified.")

    if args.control_maf is None and args.control_gistic is None:
        remove("recurrence fold change")
        for panel in args.panels_to_plot:
            if panel.endswith("control"):
                remove(panel)
    if "recurrence" not in args.panels_to_plot:
        remove("total recurrence overall")
    if "recurrence control" not in args.panels_to_plot:
        remove("total recurrence overall control")
    if "recurrence fold change" not in args.panels_to_plot:
        remove("total recurrence fold change")

    if args.cohort_label is None:
        remove("cohort label")
    if args.control_cohort_label is None:
        remove("cohort label control")

    if args.snv_interesting_genes is None and args.cnv_interesting_genes is None:
        remove("model annotation")
    if "model annotation" not in args.panels_to_plot:
        remove("model annotation legend")

    if args.mutsig is None:
        remove("mutational signatures")
    if "mutational signatures" not in args.panels_to_plot:
        remove("mutational signatures legend")

    if args.sif is None:
        remove("meta data")
    if "meta data" not in args.panels_to_plot:
        remove("meta data legend")

    if args.gene_meta_data is None:
        remove("gene meta data")


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
        "--mutsig", type=str, action='append', default=None,
        help="Path to a file containing mutational signature exposures (index: patient, columns: signatures)."
    )
    parser.add_argument(
        "--cohort-label", type=str, default=None,
        help="Name of the cohort for panel title."
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
        "--control-mutsig", type=str, action='append', default=None,
        help="Path to a file containing mutational signature exposures (index: patient, columns: signatures)."
    )
    parser.add_argument(
        "--control-cohort-label", type=str, default=None,
        help="Name of the control cohort for panel title."
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
        "--drop-empty-columns", action="store_true",
        help="Exclude samples/patients with no mutation data."
    )

    parser.add_argument(
        "--meta-data-rows", type=parse_comma_separated,
        default=[
            "WGD",
            "Ploidy",
            "Tumor Purity",
            "Subclonal Fraction",
            "Contamination",
            "Material",
            "Platform",
            "has matched N",
            "HR Status",
            "HER2 Status",
            "Chemo Tx",
            "XRT",
            "Targeted Tx",
            "Hormone Tx",
            "Immuno Tx ICI",
            "ADC",
            "Age at P Dx",
            "Sex",
            "Sample Type",
            "Histology"
        ],
        help="Comma-separated list of metadata fields to display.")
    parser.add_argument(
        "--meta-data-rows-per-sample", type=parse_comma_separated,
        default=[
            "WGD",
            "Ploidy",
            "Tumor Purity",
            "Subclonal Fraction",
            "Contamination",
            "Material",
            "Platform",
            "HR Status",
            "HER2 Status",
            "Sample Type",
        ],
        help="Metadata fields to display per sample."
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
        "--low-amp-threshold", type=int, default=1,
        help="Threshold for low-level amplification."
    )
    parser.add_argument(
        "--high-amp-threshold", type=int, default=2,
        help="Threshold for high-level amplification."
    )
    parser.add_argument(
        "--low-del-threshold", type=int, default=-1,
        help="Threshold for low-level deletion."
    )
    parser.add_argument(
        "--high-del-threshold", type=int, default=-2,
        help="Threshold for high-level deletion."
    )
    parser.add_argument(
        "--show-low-level-cnvs", default=True, action="store_true",
        help="Include low-level CNVs in the plot."
    )

    parser.add_argument(
        "--panels-to-plot", type=parse_comma_separated,
        default=[
            "comutation",
            "comutation control",
            "cohort label",
            "cohort label control",
            "tmb",
            "tmb control",
            "tmb legend",
            "mutational signatures",
            "mutational signatures control",
            "mutational signatures legend",
            "recurrence",
            "recurrence control",
            "recurrence fold change",
            "total recurrence fold change",
            # "total recurrence",
            # "total recurrence control",
            "total recurrence overall",
            "total recurrence overall control",
            "gene meta data",
            "cytoband",
            "gene names",
            "snv legend",
            "cnv legend",
            "model annotation",
            "model annotation legend",
            "meta data",
            "meta data control",
            "meta data legend"
        ],
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
        "--max-xfigsize-scale", type=int, default=1,
        help="Parameter by which xfigsize scales with the number of columns."
    )

    args = parser.parse_args()
    validate_args(args)
    return args
