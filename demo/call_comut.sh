#!/bin/bash

python make_data.py
python ../call.py \
  -o "comut_test.png" \
  --maf "test.maf.tsv" \
  --gistic "test.all_thresholded.by_genes.txt" \
  --sif "test.sif.tsv" \
  --signatures "test.mutsig.tsv" \
  --snv-interesting-genes "Gene 1,Gene 2,Gene 3,Gene 4,Gene 5,Gene 6,Gene 7,Gene 8,Gene 9,Gene 10,Gene 11,Gene 12" \
  --cnv-interesting-genes "Gene 9,Gene 10,Gene 11,Gene 12,Gene 13,Gene 14,Gene 15,Gene 16,Gene 17,Gene 18,Gene 19,Gene 20" \
  --ground-truth-genes "control:Gene 7,Gene 9,Gene 11,Gene 13,Gene 15;case:Gene 10,Gene 8" \
  --meta-data-rows "WGD,Ploidy,Tumor Purity,Subclonal Fraction,Contamination,Material,Platform,has matched N,HR Status,HER2 Status,Chemo Tx,XRT,Targeted Tx,Hormone Tx,Immuno Tx ICI,ADC,Sex,Sample Type,Histology" \
  --palette "|Material>FFPE:0,0,0" \
  --recurrence-categories "max" \
  --column-sort-by "COMUT,TMB" \
  --max-xfigsize 8

python ../call.py \
  -o "comut_test.tiny.png" \
  --maf "test.maf.tsv" \
  --gistic "test.all_thresholded.by_genes.txt" \
  --sif "test.sif.tsv" \
  --signatures "test.mutsig.tsv" \
  --snv-interesting-genes "Gene 1,Gene 2,Gene 3,Gene 4,Gene 5,Gene 6,Gene 7,Gene 8,Gene 9,Gene 10,Gene 11,Gene 12" \
  --cnv-interesting-genes "Gene 9,Gene 10,Gene 11,Gene 12,Gene 13,Gene 14,Gene 15,Gene 16,Gene 17,Gene 18,Gene 19,Gene 20" \
  --ground-truth-genes "control:Gene 7,Gene 9,Gene 11,Gene 13,Gene 15;case:Gene 10,Gene 8" \
  --meta-data-rows "WGD,Ploidy,Tumor Purity,Subclonal Fraction,Contamination,Material,Platform,has matched N,HR Status,HER2 Status,Chemo Tx,XRT,Targeted Tx,Hormone Tx,Immuno Tx ICI,ADC,Sex,Sample Type,Histology" \
  --palette "|Material>FFPE:0,0,0" \
  --recurrence-categories "max" \
  --column-sort-by "COMUT,TMB" \
  --panels-to-plot "tmb,comutation" \
  --max-xfigsize 8

python ../call.py \
  -o "comut_test.control.png" \
  --maf "test.maf.tsv" \
  --gistic "test.all_thresholded.by_genes.txt" \
  --sif "test.sif.tsv" \
  --signatures "test.mutsig.tsv" \
  --cohort-label "CASE" \
  --control-cohort-label "CONTROL" \
  --control-maf "control.maf.tsv" \
  --control-gistic "control.all_thresholded.by_genes.txt" \
  --control-sif "control.sif.tsv" \
  --control-signatures "control.mutsig.tsv" \
  --snv-interesting-genes "Gene 1,Gene 2,Gene 3,Gene 4,Gene 5,Gene 6,Gene 7,Gene 8,Gene 9,Gene 10,Gene 11,Gene 12" \
  --cnv-interesting-genes "Gene 9,Gene 10,Gene 11,Gene 12,Gene 13,Gene 14,Gene 15,Gene 16,Gene 17,Gene 18,Gene 19,Gene 20" \
  --ground-truth-genes "control:Gene 7,Gene 9,Gene 11,Gene 13,Gene 15;case:Gene 10,Gene 8" \
  --meta-data-rows "WGD,Ploidy,Tumor Purity,Subclonal Fraction,Contamination,Material,Platform,has matched N,HR Status,HER2 Status,Chemo Tx,XRT,Targeted Tx,Hormone Tx,Immuno Tx ICI,ADC,Sex,Sample Type,Histology" \
  --palette "|Material>FFPE:0,0,0" \
  --recurrence-categories "snv,amp,del|Gene 10:snv,amp;Gene 11:amp;Gene 12:amp;Gene 9:amp;Gene 15:del;Gene 18:del;Gene 8:snv;Gene 13:max;Gene 19:del;Gene 16:del;Gene 17:del;Gene 14:amp;Gene 20:del;Gene 6:max;Gene 1:max;Gene 2:max" \
  --column-sort-by "COMUT" \
  --max-xfigsize 6

# Same control comparison, but with the control cohort placed on the LEFT
# (case and control swapped, including the fold-change "case"/"control" arrows).
python ../call.py \
  -o "comut_test.control.left.png" \
  --maf "test.maf.tsv" \
  --gistic "test.all_thresholded.by_genes.txt" \
  --sif "test.sif.tsv" \
  --signatures "test.mutsig.tsv" \
  --cohort-label "CASE" \
  --control-cohort-label "CONTROL" \
  --control-maf "control.maf.tsv" \
  --control-gistic "control.all_thresholded.by_genes.txt" \
  --control-sif "control.sif.tsv" \
  --control-signatures "control.mutsig.tsv" \
  --control-position left \
  --snv-interesting-genes "Gene 1,Gene 2,Gene 3,Gene 4,Gene 5,Gene 6,Gene 7,Gene 8,Gene 9,Gene 10,Gene 11,Gene 12" \
  --cnv-interesting-genes "Gene 9,Gene 10,Gene 11,Gene 12,Gene 13,Gene 14,Gene 15,Gene 16,Gene 17,Gene 18,Gene 19,Gene 20" \
  --ground-truth-genes "control:Gene 7,Gene 9,Gene 11,Gene 13,Gene 15;case:Gene 10,Gene 8" \
  --meta-data-rows "WGD,Ploidy,Tumor Purity,Subclonal Fraction,Contamination,Material,Platform,has matched N,HR Status,HER2 Status,Chemo Tx,XRT,Targeted Tx,Hormone Tx,Immuno Tx ICI,ADC,Sex,Sample Type,Histology" \
  --palette "|Material>FFPE:0,0,0" \
  --recurrence-categories "snv,amp,del|Gene 10:snv,amp;Gene 11:amp;Gene 12:amp;Gene 9:amp;Gene 15:del;Gene 18:del;Gene 8:snv;Gene 13:max;Gene 19:del;Gene 16:del;Gene 17:del;Gene 14:amp;Gene 20:del;Gene 6:max;Gene 1:max;Gene 2:max" \
  --column-sort-by "COMUT" \
  --max-xfigsize 6

# Grid sub-panels: split the figure into a 2-D grid of comut sub-panels,
# stratified by a sample metadata key (grid columns) and explicit gene sets
# (grid rows). Ungrouped genes are collected into a trailing "Other" row.
python ../call.py \
  -o "comut_test.grid.png" \
  --maf "test.maf.tsv" \
  --gistic "test.all_thresholded.by_genes.txt" \
  --sif "test.sif.tsv" \
  --signatures "test.mutsig.tsv" \
  --snv-interesting-genes "Gene 1,Gene 2,Gene 3,Gene 4,Gene 5,Gene 6,Gene 7,Gene 8,Gene 9,Gene 10,Gene 11,Gene 12" \
  --cnv-interesting-genes "Gene 9,Gene 10,Gene 11,Gene 12,Gene 13,Gene 14,Gene 15,Gene 16,Gene 17,Gene 18,Gene 19,Gene 20" \
  --meta-data-rows "WGD,Ploidy,Tumor Purity,Material,Platform,Sex,Sample Type,Histology" \
  --palette "|Material>FFPE:0,0,0" \
  --recurrence-categories "max" \
  --column-sort-by "COMUT,TMB" \
  --column-group-by "Sample Type" \
  --group-order "BM,EM" \
  --column-group-labels "Brain metastases" "Extracranial metastases" \
  --gene-groups "SNV genes:Gene 1,Gene 2,Gene 3,Gene 4,Gene 5,Gene 6;CNV genes:Gene 13,Gene 14,Gene 15,Gene 16,Gene 17,Gene 18" \
  --max-xfigsize 12

