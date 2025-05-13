#!/bin/bash

python make_data.py
python ../call.py \
  -o "comut_test.png" \
  --maf "test.maf.tsv" \
  --gistic "test.all_thresholded.by_genes.txt" \
  --sif "test.sif.tsv" \
  --mutsig "test.mutsig.tsv" \
  --snv-interesting-genes "Gene 1,Gene 2,Gene 3,Gene 4,Gene 5,Gene 6,Gene 7,Gene 8,Gene 9,Gene 10,Gene 11,Gene 12" \
  --cnv-interesting-genes "Gene 9,Gene 10,Gene 11,Gene 12,Gene 13,Gene 14,Gene 15,Gene 16,Gene 17,Gene 18,Gene 19,Gene 20" \
  --ground-truth-genes "control:Gene 7,Gene 9,Gene 11,Gene 13,Gene 15;case:Gene 10,Gene 8" \
  --meta-data-rows "WGD,Ploidy,Tumor Purity,Subclonal Fraction,Contamination,Material,Platform,has matched N,HR Status,HER2 Status,Chemo Tx,XRT,Targeted Tx,Hormone Tx,Immuno Tx ICI,ADC,Sex,Sample Type,Histology" \
  --palette "|Material>FFPE:0,0,0" \
  --column-sort-by "COMUT,TMB" \
  --max-xfigsize 8
python ../call.py \
  -o "comut_test.control.png" \
  --maf "test.maf.tsv" \
  --gistic "test.all_thresholded.by_genes.txt" \
  --sif "test.sif.tsv" \
  --mutsig "test.mutsig.tsv" \
  --control-maf "control.maf.tsv" \
  --control-gistic "control.all_thresholded.by_genes.txt" \
  --control-sif "control.sif.tsv" \
  --control-mutsig "control.mutsig.tsv" \
  --snv-interesting-genes "Gene 1,Gene 2,Gene 3,Gene 4,Gene 5,Gene 6,Gene 7,Gene 8,Gene 9,Gene 10,Gene 11,Gene 12" \
  --cnv-interesting-genes "Gene 9,Gene 10,Gene 11,Gene 12,Gene 13,Gene 14,Gene 15,Gene 16,Gene 17,Gene 18,Gene 19,Gene 20" \
  --ground-truth-genes "control:Gene 7,Gene 9,Gene 11,Gene 13,Gene 15;case:Gene 10,Gene 8" \
  --meta-data-rows "WGD,Ploidy,Tumor Purity,Subclonal Fraction,Contamination,Material,Platform,has matched N,HR Status,HER2 Status,Chemo Tx,XRT,Targeted Tx,Hormone Tx,Immuno Tx ICI,ADC,Sex,Sample Type,Histology" \
  --palette "|Material>FFPE:0,0,0" \
  --column-sort-by "COMUT" \
  --max-xfigsize 8

