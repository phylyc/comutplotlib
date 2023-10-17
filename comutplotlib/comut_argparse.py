import argparse


def parse_args():

    parser = argparse.ArgumentParser(
        prog="ComutPlot",
        description="Plotting",
        epilog="Text at the bottom of help",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.usage = '''
    ComutPlot is a powerful tool for generating interactive genomic comutation plots.

    Usage:
      python your_script.py -o OUTPUT_FILE [--maf MAF_FILE] [--sif SIF_FILE] [--gistic GISTIC_FILE]
                             [--by {Sample, Patient}] [--label_columns {True, False}]
                             [--column_order COLUMN_ORDER] [--index_order INDEX_ORDER]
                             [--meta META_FILE] [--meta_data_rows META_DATA_ROWS]
                             [--meta_data_rows_per_sample META_DATA_ROWS_PER_SAMPLE]
                             [--genes GENES] [--snv_interesting_genes SNV_INTERESTING_GENES]
                             [--cnv_interesting_genes CNV_INTERESTING_GENES]
                             [--interesting_gene_comut_threshold_percent GENE_COMUT_THRESHOLD_PERCENT]
                             [--ground_truth_genes GROUND_TRUTH_GENES]
                             [--low_amp_threshold LOW_AMP_THRESHOLD]
                             [--high_amp_threshold HIGH_AMP_THRESHOLD]
                             [--low_del_threshold LOW_DEL_THRESHOLD]
                             [--high_del_threshold HIGH_DEL_THRESHOLD]
                             [--total_prevalence_threshold TOTAL_PREVALENCE_THRESHOLD]
                             [--model_names MODEL_NAMES]
                             [--panels_to_plot PANELS_TO_PLOT]
                             [--max_xfigsize MAX_XFIGSIZE] [--yfigsize YFIGSIZE]
                             [--sub_folders SUB_FOLDERS] [--file_name_prefix FILE_NAME_PREFIX]
                             [--file_name_suffix FILE_NAME_SUFFIX]
                             [--regions REGIONS] [--mutsig MUTSIG] [--cytoband CYTOBAND]
                             [--tmb TMB] [--prev PREV]

    Required Arguments:
      -o, --output OUTPUT_FILE
          Path to the output file.
    
    Semi-Optional Arguments (at least one of the two is required):
      --maf MAF_FILE
          Path to the MAF input file (direct output from Funcotator GATK).

      --gistic GISTIC_FILE
          Path to the GISTIC input file.

    Optional Arguments:
      --sif SIF_FILE
          Path to the SIF input file.

      --by {Sample, Patient}
          Group the data by Sample or Patient. (Default: Patient)

      --label_columns {True, False}
          Include labels at the bottom of the plot. (Default: True)

      --column_order COLUMN_ORDER
          Order of the columns on the plot.

      --index_order INDEX_ORDER
          Order of the index on the plot.

      --meta META_FILE
          Path to the meta input file.

      --meta_data_rows META_DATA_ROWS
          Comma separated list of SIF columns to plot.

      --meta_data_rows_per_sample META_DATA_ROWS_PER_SAMPLE
          Comma separated list of SIF columns to plot per sample.

      --genes GENES
          Comma separated list of genes.

      --snv_interesting_genes SNV_INTERESTING_GENES
          Comma separated list of SNV genes.

      --cnv_interesting_genes CNV_INTERESTING_GENES
          Comma separated list of CNV genes.

      --interesting_gene_comut_threshold_percent GENE_COMUT_THRESHOLD_PERCENT
          Threshold of significance for interesting genes.

      --ground_truth_genes GROUND_TRUTH_GENES
          Dictionary of genes to compare to. Format: key1:value1,value2,...;key2:value3,value4,...

      --low_amp_threshold LOW_AMP_THRESHOLD
          Threshold for low amplification.

      --high_amp_threshold HIGH_AMP_THRESHOLD
          Threshold for high amplification.

      --low_del_threshold LOW_DEL_THRESHOLD
          Threshold for low deletion.

      --high_del_threshold HIGH_DEL_THRESHOLD
          Threshold for high deletion.

      --total_prevalence_threshold TOTAL_PREVALENCE_THRESHOLD
          Threshold for total prevalence.

      --model_names MODEL_NAMES
          Comma separated list of names of models.

      --panels_to_plot PANELS_TO_PLOT
          Comma separated list of parts to be plotted.

      --max_xfigsize MAX_XFIGSIZE
          Maximum x figure size.

      --yfigsize YFIGSIZE
          y figure size.

      --sub_folders SUB_FOLDERS
          Sub folders.

      --file_name_prefix FILE_NAME_PREFIX
          File name prefix.

      --file_name_suffix FILE_NAME_SUFFIX
          File name suffix.

      --regions REGIONS
          Path to the covered regions input file.

      --mutsig MUTSIG
          Path to mutational signatures input file.

      --cytoband CYTOBAND
          Path to the cytoband input file.

      --tmb TMB
          Path to the TMB input file.

      --prev PREV
          Path to the prevalence input file.
    '''
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
    if args.column_order is not None:
        args.column_order = args.column_order.split(",")
    if args.index_order is not None:
        args.index_order = args.index_order.split(",")
    args.meta_data_rows = args.meta_data_rows.split(",")
    args.meta_data_rows_per_sample = args.meta_data_rows_per_sample.split(",")
    args.model_names = args.model_names.split(",")
    args.panels_to_plot = args.panels_to_plot.split(",")

    return args
