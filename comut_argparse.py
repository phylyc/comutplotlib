import argparse

parser = argparse.ArgumentParser(
    prog="ComutPlot",
    description="Plotting",
    epilog="Text at the bottom of help",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.usage = 'use to...'

parser.add_argument("-o", "--output", type=str, required=True, help="Path to the output file.")
parser.add_argument("--maf", type=str, required=False, default=None, help="Path to the MAF input file.")
parser.add_argument("--sif", type=str, required=False, default=None, help="Path to the SIF input file.")
parser.add_argument("--gistic", type=str, required=False, default=None, help="Path to the MAF input file.")
parser.add_argument("--meta_data_rows_per_sample", type=str, required=False, default="Sample Type,Material,Contamination,Tumor Purity,Platform", help="Comma separated list of SIF columns to plot per sample.")
parser.add_argument("--by", type=str, choices=["Sample", "Patient"], required=False, default="Patient")
parser.add_argument("--label_columns", type=str, required=False, default=False, help="labels at the bottom of plot")
parser.add_argument("--column_order", type=str, required=False, default=None, help="order of the columns on plot")
parser.add_argument("--index_order", type=str, required=False, default=None, help="order of the index on plot")
parser.add_argument("--meta", type=str, required=False, default=None, help="Path to the meta input file.")
parser.add_argument("--meta_data_rows", type=str, required=False, default="Sample Type,Material,Contamination,Tumor Purity,Platform,has_matched_N,Sex,Histology", help="Comma separated list of SIF columns to plot.")
parser.add_argument("--genes", type=str, required=False, default=None, help="Comma separated list of genes.")
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
parser.add_argument("--model_names", type=str, required=False, default="SNV burden model,CNV burden model", help="names of models.")
parser.add_argument("--parts_to_plot", type=str, required=False, default="mutation burden,mutation burden legend,recurrence,prevalence,total prevalence,cytoband,comutations,comutations legend,meta data,meta data legend", help="parts to be plotted")
parser.add_argument("--max_xfigsize", type=int, required=False, default=None, help="maximum x figure size")
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
if args.genes is not None:
    args.genes = args.genes.split(",")
if args.snv_interesting_genes is not None:
    args.snv_interesting_genes = args.snv_interesting_genes.split(",")
if args.cnv_interesting_genes is not None:
    args.cnv_interesting_genes = args.cnv_interesting_genes.split(",")

args.meta_data_rows = args.meta_data_rows.split(",")
args.meta_data_rows_per_sample = args.meta_data_rows_per_sample.split(",")
args.model_names = args.model_names.split(",")
args.parts_to_plot = args.parts_to_plot.split(",")

print(args.output)
print(args.genes)
print(args.parts_to_plot)
print(args.cnv_interesting_genes)

# {"MYC": [1,1,1]} RGB for palette argument
