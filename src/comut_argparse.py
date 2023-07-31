import argparse


def parse_args():

    parser = argparse.ArgumentParser(
        prog="ComutPlot",
        description="Plotting",
        epilog="Text at the bottom of help",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.usage = 'use to...'

    def validate_flags(args):
        if args.maf is None and args.gistic is None:
            raise ValueError("Either --maf or --gistic must be specified.")

    parser.add_argument("-o", "--output", type=str, required=True, help="Path to the output file.")
    parser.add_argument("--maf", type=str, required=False, default=None, action='append', help="Path to the MAF input file.")
    parser.add_argument("--sif", type=str, required=False, default=None, action='append', help="Path to the SIF input file.")
    parser.add_argument("--gistic", type=str, required=False, default=None, action='append', help="Path to the MAF input file.")
    parser.add_argument("--mutsig", type=str, required=False, default=None, action='append', help="Path to mutational signatures input file.")

    parser.add_argument("--model_names", type=str, required=False, default="SNV burden model,CNV burden model", help="names of models.")
    parser.add_argument("--model_significances", type=str, required=False, default=None, action="append", help="Path to the model significances input file.")

    parser.add_argument("--by", type=str, choices=["Sample", "Patient"], required=False, default="Patient")
    parser.add_argument("--label_columns", type=bool, required=False, default=True, help="labels at the bottom of plot")
    parser.add_argument("--column_order", type=str, required=False, default=None, help="order of the columns on plot")
    parser.add_argument("--index_order", type=str, required=False, default=None, help="order of the index on plot")

    parser.add_argument("--meta", type=str, required=False, default=None, help="Path to the meta input file.")
    parser.add_argument("--meta_data_rows", type=str, required=False, default="Sample Type,Material,Contamination,Tumor Purity,Platform,has_matched_N,Sex,Histology", help="Comma separated list of SIF columns to plot.")
    parser.add_argument("--meta_data_rows_per_sample", type=str, required=False,
                        default="Sample Type,Material,Contamination,Tumor Purity,Platform",
                        help="Comma separated list of SIF columns to plot per sample.")

    parser.add_argument("--interesting_gene", type=str, required=False, default=None, help="Interesting gene.")
    parser.add_argument("--interesting_gene_comut_percent_threshold", type=float, required=False, default=None,
                        help="threshold of significance for interesting genes")
    parser.add_argument("--interesting_genes", type=str, required=False, default=None, help="Comma separated list of genes.")
    parser.add_argument("--snv_interesting_genes", type=str, required=False, default=None,
                        help="Comma separated list of snv genes.")
    parser.add_argument("--cnv_interesting_genes", type=str, required=False, default=None,
                        help="Comma separated list of cnv genes")
    parser.add_argument("--total_prevalence_threshold", type=float, required=False, default=None,
                        help="Minimum percentage of patients to have a mutation in this gene to be plotted")
    parser.add_argument("--ground_truth_genes", type=str, required=False, default=None, action="append",
                        help="dictionary of color in palette as keys with list of genes to be colored as values")

    parser.add_argument("--low_amp_threshold", type=int, required=False, default=1, help="threshold for low amplification")
    parser.add_argument("--high_amp_threshold", type=int, required=False, default=2, help="threshold for high amplification")
    parser.add_argument("--low_del_threshold", type=int, required=False, default=-1, help="threshold for low deletion")
    parser.add_argument("--high_del_threshold", type=int, required=False, default=-2, help="threshold for high deletion")

    parser.add_argument("--panels_to_plot", type=str, required=False, default="tmb,mutational signatures,recurrence,prevalence,total prevalence,total prevalence overall,cytoband,gene names,model annotation,comutation,mutsig legend,snv legend,cnv legend,model annotation legend,meta data,meta data legend", help="parts to be plotted")
    parser.add_argument("--max_xfigsize", type=int, required=False, default=None, help="maximum x figure size; central comutation plot will be scaled to fit this size")

    args = parser.parse_args()
    validate_flags(args)

    if args.interesting_genes is not None:
        args.interesting_genes = args.interesting_genes.split(",")
    if args.snv_interesting_genes is not None:
        args.snv_interesting_genes = args.snv_interesting_genes.split(",")
    if args.cnv_interesting_genes is not None:
        args.cnv_interesting_genes = args.cnv_interesting_genes.split(",")
    if args.ground_truth_genes is not None:
        args.ground_truth_genes = {
            k: v.split(",") for k, v in [x.split(":") for x in args.ground_truth_genes]
        }
    args.meta_data_rows = args.meta_data_rows.split(",")
    args.meta_data_rows_per_sample = args.meta_data_rows_per_sample.split(",")
    args.model_names = args.model_names.split(",")
    args.panels_to_plot = args.panels_to_plot.split(",")

    return args
