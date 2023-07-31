#!/bin/bash

python make_data.py
python ../call.py \
  -o "comut_test.pdf" \
  --maf "test.maf.tsv" \
  --gistic "test.all_thresholded.by_genes.txt" \
  --sif "test.sif.tsv" \
  --mutsig "test.mutsig.tsv" \
  --snv_interesting_genes "Gene 1,Gene 2,Gene 3,Gene 4,Gene 5,Gene 6,Gene 7,Gene 8,Gene 9,Gene 10,Gene 11,Gene 12" \
  --cnv_interesting_genes "Gene 9,Gene 10,Gene 11,Gene 12,Gene 13,Gene 14,Gene 15,Gene 16,Gene 17,Gene 18,Gene 19,Gene 20" \
  --ground_truth_genes "control:Gene 7,Gene 9,Gene 11,Gene 13,Gene 15" \
  --ground_truth_genes "case:Gene 10,Gene 8" \
  --max_xfigsize 8

