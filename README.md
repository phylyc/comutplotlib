# comut_plot

**DESCRIPTION**:
ComutPlot is a tool for plotting and visualizing genomic data. It offers
a flexible and customizable way to generate informative plots for various
genomic analyses. This Python script uses the argparse library to parse
command-line arguments and allows users to specify different input files, 
options, and configurations.

**INSTALLATION**: 
see install.sh, setup.cfg, and required.txt

**USAGE GUIDE**: 
To use ComutPlot, simply execute the script with the desired options and arguments. 
For a comprehensive list of available options and their descriptions, run the following command:
`python comut_argparse.py --help`

**EXAMPLES**:

1. Generate a plot with a single MAF file and specific output path: 
`python comut_argparse.py --maf file1.maf --output output_plot.png`

2. Use multiple input files (maf, sif, gistic) and set additional options: 
`python comut_argparse.py --maf file1.maf --sif file2.sif --gistic file3.gistic --by Sample --label_columns Labels --output output_plot.png`

3. Provide multiple interesting genes and customize the plot layout: 
`python comut_argparse.py --genes geneA,geneB,geneC --snv_interesting_genes geneX,geneY --parts_to_plot mutation burden,comutations,cytoband --output output_plot.png`

**AUTHOR**:
ComutPlot is authored by [Philipp Hähnel].
For any questions or support related to ComutPlot, please contact [Philipp_Hahnel@DFCI.HARVARD.EDU].

Contribution Guidelines:

Project Structure:

Dependencies:

Known Issues or Limitations: