import argparse

parser = argparse.ArgumentParser(
    prog="ComutPlot",
    description="Plotting",
    epilog="Text at the bottom of help")

parser.usage = 'use to...'

parser.add_argument("-o", "--output", type=str, required=True, help="Path to the output file.")
parser.add_argument("--maf", type=str, required=False, default=None, help="Path to the MAF input file.")
parser.add_argument("--sif", type=str, required=False, default=None, help="Path to the SIF input file.")
parser.add_argument("--gistic", type=str, required=True, default=None, help="Path to the MAF input file.")
parser.add_argument("--meta_data_rows_per_sample", type=str, required=False, default=None, help="list of meta data rows separated by sample")
parser.add_argument("--by", type=str, choices=["Sample", "Patient"], required=False, default="Patient")
parser.add_argument("--label_columns", type=str, required=False, default=False, help="labels at the bottom of plot")
parser.add_argument("--column_order", type=str, required=False, default=None, help="order of the columns on plot")
parser.add_argument("--index_order", type=str, required=False, default=None, help="order of the index on plot")
parser.add_argument("--meta", type=str, required=False, default=None, help="Path to the meta input file.")
parser.add_argument("--genes", type=str, required=False, default=None, help="Comma separated List of genes.")
parser.add_argument("--snv_interesting_genes", type=str, required=False, default=None,
                    help="Comma separated list of snv genes.")
parser.add_argument("--cnv_interesting_genes", type=str, required=False, default=None,
                    help="Comma separated list of cnv genes")
parser.add_argument("--interesting_gene_comut_threshold_percent", type=float, required=False, default=None,
                    help="threshold of significance for interesting genes")
parser.add_argument("--ground_truth_genes", type=str, required=False, default=None,
                    help="dictionary of genes to compare to")
parser.add_argument("--low_amp_threshold", type=int, required=False, default=1, help="threshold for low amplification")
parser.add_argument("--high_amp_threshold", type=int, required=False, default=2, help="threshold for high amplification")
parser.add_argument("--low_del_threshold", type=int, required=False, default=-1, help="threshold for low deletion")
parser.add_argument("--high_del_threshold", type=int, required=False, default=-2, help="threshold for high deletion")
parser.add_argument("--total_prevalence_threshold", type=float, required=False, default=None,
                    help="threshold for total prevalence")
parser.add_argument("--model_names", type=str, required=False, default=None, help="names of models.")
parser.add_argument("--parts_to_plot", type=str, required=False, default=None, help="parts to be plotted")
parser.add_argument("--max_xfigsize", type=int, required=False, defualt=None, help="maximum x figure size")
parser.add_argument("--yfigsize", type=int, required=False, default=None, help ="y figure size")
parser.add_argument("--sub_folders", type=str, required=False, default=None, help="sub folders")
parser.add_argument("--file_name_prefix", type=str, required=False, default="", help="file name prefix")
parser.add_argument("--file_name_suffix", type=str, required=False, default="", help="file name suffix")
parser.add_argument("--regions", type=str, required=False, default=None, help="Path to the covered regions input file.")
parser.add_argument("--mutsig", type=str, required=False, default=None, help="Path to mutational signatures input file.")
parser.add_argument("--cytoband", type=str, required=False, default=None, help="Path to the MAF input file.")
parser.add_argument("--tmb", type=str, required=False, default=None, help="Path to the TMB input file.")
parser.add_argument("--prev", type=str, required=False, default=None, help="Path to the prevalence input file.")

args = parser.parse_args()
args.genes = args.genes.split(",")
args.snv.genes = args.snv_interesting_genes.split(",")
args.cnv.genes = args.cnv_interesting_genes.split(",")

if args.meta_data_rows is None:
        args_meta_data_rows = [
    "Sample Type",
    "Material",
    "Contamination",
    "Tumor Purity",
    "Platform",
    "has_matched_N",
    "Sex",
    "Histology"
]
else:
    args.meta_data_rows = args.meta_data_rows.split(",")

if args.meta_data_rows_per_sample is None:
        args_meta_data_rows_per_sample = [
            "Sample Type",
            "Material",
            "Platform",
            "Contamination",
            "Tumor Purity"
        ]
else:
    args.meta_data_rows_per_sample = args.meta_data_rows_per_sample.split(",")

if args.model_names is None:
        args_model_names = [
            "SNV burden model",
            "CNV burden model",
        ]
else:
    args.model_names = args.model_names.split(",")

if args.parts_to_plot is None:
    args_parts_to_plot = [
        "mutation burden",
        "mutation burden legend",
        "recurrence",
        "prevalence",
        "total prevalence",
        "cytoband",
        "comutations",
        "comutations legend",
        "meta data",
        "meta data legend",
    ]
else:
    args.parts_to_plot = args.parts_to_plot.split(",")

print(args.output)
print(args.genes)
print(args.snv.genes)
print(args.cnv.genes)

# {"MYC": [1,1,1]} RGB for palette argument
