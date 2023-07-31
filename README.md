# ComutPlot - Interactive Genomic Comutation Plots
Welcome to ComutPlot, a powerful tool for generating interactive genomic comutation plots. 
Comutation plots are a widely used visualization method to explore co-occurrence patterns 
of genomic alterations, such as mutations, copy number variations (CNVs), and other genetic 
events, across multiple samples or patients.

![Sample Comut Plot Output](https://github.com/phylyc/comut_plot/demo/comut_test.pdf)

### DESCRIPTION:
ComutPlot is a tool for plotting and visualizing genomic data. It offers
a flexible and customizable way to generate informative plots for various
genomic analyses. This Python script uses the argparse library to parse
command-line arguments and allows users to specify different input files, 
options, and configurations.

### INSTALLATION: 
To install comut_plot, please refer to the following files:
- `install.sh`: Shell script for installation
- `setup.cfg`: Configuration file
- `setup.py`: Setup script
- `required.txt`: List of required dependencies

### USAGE GUIDE: 
To use ComutPlot, simply execute the script with the desired options and arguments.  
For a comprehensive list of available options and their descriptions, run the following command:
`python comut_argparse.py --help`

#### INPUT FILES:
There are three main inputs to provide. Most of them are direct outputs of other tools to minimize special formatting needs:
1. **MAF file**: This should be the output of Funcotator GATK. Refer to `mutation_annotation.py` for related fields.
2. **GISTIC file**: The output of GISTIC 2.0, which stands for Genomic Identification of Significant Targets in Cancer.
3. **SIF file**: A tab-separated file. Please see `sample_annotation.py` for column names that need to be specified.

*Note: At least one of the MAF or GISTIC files must be provided as input for the program to run.*

#### EXAMPLES:

- Generate a plot with a single MAF file and specific output path:   
`python comut_argparse.py --maf file1.maf --output output_plot.png`  

- Use multiple input files (MAF, SIF, GISTIC) and set additional options:   
`python comut_argparse.py --maf file1.maf --sif file2.sif --gistic file3.gistic --by Sample --label_columns True --output output_plot.png`  

- Provide multiple interesting genes and customize the plot layout:   
`python comut_argparse.py --genes geneA,geneB,geneC --snv_interesting_genes geneX,geneY --parts_to_plot mutation burden,comutations,cytoband --output output_plot.png`  

- For more complicated demo figures, please run the code in bash/call_comut.sh:  
`bash call_comut.sh`

Contribution Guidelines (don't need): 
if want to suggest new features, create a new issue; 
if want to contribute, create a branch & make a pull request;

### DEPENDENCIES: 
See `requirements.txt` for the required dependencies

### AUTHOR:
ComutPlot is authored by **[Philipp Hähnel]**.  
For any questions or support related to ComutPlot, please contact **[Philipp_Hahnel@DFCI.HARVARD.EDU](mailto:Philipp_Hahnel@DFCI.HARVARD.EDU)**.

We hope you find ComutPlot to be a valuable tool for your genomic data analysis and visualization needs. Happy plotting!