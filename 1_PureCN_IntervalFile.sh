#!/bin/bash
set -e

PATH=/home/lierhan/Program/miniconda3/bin:$PATH
export PURECN="/home/lierhan/Program/miniconda3/lib/R/library/PureCN/extdata"

target_bed="target.bed"
reference_genome="hg19.fa"
target_intervals_file="target_intervals.txt"
target_intervals_bed="target_intervals.bed"

Rscript ${PURECN}/IntervalFile.R \
  --in-file ${target_bed} \
  --fasta ${reference_genome} \
  --out-file target_intervals.txt \
  --off-target \
  --genome hg19 \
  --export ${target_intervals_bed} \
  --mappability wgEncodeCrgMapabilityAlign100mer.bigWig \
  --reptiming wgEncodeUwRepliSeqK562WaveSignalRep1.bigWig