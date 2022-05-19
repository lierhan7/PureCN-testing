#!/bin/bash
set -e

PATH=/home/lierhan/Program/miniconda3/bin:$PATH
export PURECN="/home/lierhan/Program/miniconda3/lib/R/library/PureCN/extdata"

output_dir=""
coverage_files_list_file=""
normal_panel_file="normals.merged.min5.vcf.gz"

Rscript ${PURECN}/NormalDB.R \
  --out-dir ${output_dir} \
  --coverage-files ${coverage_files_list_file} \
  --normal-panel ${normal_panel_file} \
  --genome hg19 \
  --assay foundation
