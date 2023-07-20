# comut_plot

**DESCRIPTION**:
ComutPlot is a tool for plotting and visualizing genomic data. It offers
a flexible and customizable way to generate informative plots for various
genomic analyses. This Python script uses the argparse library to parse
command-line arguments and allows users to specify different input files, 
options, and configurations.

**INSTALLATION**: 
see install.sh and setup.cfg

**USAGE GUIDE**: 
To use ComutPlot, simply execute the script with the desired options 
and arguments. Below are the available command-line options:

General Options:
-o, --output: Path to the output file. (Required)
Input Files:
--maf: Path to the MAF input file. This option can be specified multiple times to provide multiple MAF files.
--sif: Path to the SIF input file. This option can be specified multiple times to provide multiple SIF files.
--gistic: Path to the GISTIC input file. This option can be specified multiple times to provide multiple GISTIC files.
--meta: Path to the meta input file.
--regions: Path to the covered regions input file.
--mutsig: Path to mutational signatures input file.
--cytoband: Path to the cytoband input file.
--tmb: Path to the TMB (tumor mutational burden) input file.
--prev: Path to the prevalence input file.
Other Options:
--by: Choose between "Sample" or "Patient" to define how data should be grouped (default: "Patient").
--label_columns: Labels at the bottom of the plot.
--column_order: Order of the columns on the plot.
--index_order: Order of the index on the plot.
--meta_data_rows: Comma-separated list of SIF columns to plot.
--meta_data_rows_per_sample: Comma-separated list of SIF columns to plot per sample.
--genes: Comma-separated list of genes.
--snv_interesting_genes: Comma-separated list of SNV genes.
--cnv_interesting_genes: Comma-separated list of CNV genes.
--interesting_gene_comut_threshold_percent: Threshold of significance for interesting genes.
--ground_truth_genes: Dictionary of genes to compare to.
--low_amp_threshold: Threshold for low amplification.
--high_amp_threshold: Threshold for high amplification.
--low_del_threshold: Threshold for low deletion.
--high_del_threshold: Threshold for high deletion.
--total_prevalence_threshold: Threshold for total prevalence.
--model_names: Names of models.
--parts_to_plot: Parts to be plotted.
--max_xfigsize: Maximum x figure size.
--yfigsize: Y figure size.
--sub_folders: Sub folders.
--file_name_prefix: File name prefix.
--file_name_suffix: File name suffix.

**EXAMPLES**:

1. Generate a plot with a single MAF file and specific output path:
    python script.py --maf file1.maf --output output_plot.png

2. Use multiple input files (maf, sif, gistic) and set additional options:
    python script.py --maf file1.maf --sif file2.sif --gistic file3.gistic --by Sample --label_columns Labels --output output_plot.png

3. Provide multiple interesting genes and customize the plot layout:
    python script.py --genes geneA,geneB,geneC --snv_interesting_genes geneX,geneY --parts_to_plot mutation burden,comutations,cytoband --output output_plot.png

**AUTHOR**:
ComutPlot is authored by [Philipp Hähnel].
For any questions or support related to ComutPlot, please contact [Philipp_Hahnel@DFCI.HARVARD.EDU].

Contribution Guidelines:

Project Structure:

Dependencies:

Known Issues or Limitations: