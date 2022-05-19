#!/bin/bash
set -e

PATH=/home/lierhan/Program/miniconda3/bin:$PATH
export PURECN="/home/lierhan/Program/miniconda3/lib/R/library/PureCN/extdata"

output_dir=""
normaldb_file="normalDB_foundation_hg19.rds"
mapping_bias_file="mapping_bias_f1cdx_hg19.rds"
target_intervals_file="target_intervals.txt"
snp_blacklist_file="hg19_simpleRepeats.bed"

for tumor_sample in SRR948953 SRR948954 SRR948956 SRR948955 SRR948957 SRR948958
do
	# Production pipeline run
	Rscript $PURECN/PureCN.R --out ${output_dir} \
		--tumor ${tumor_sample}.mkdup_coverage_loess.txt.gz \
		--sampleid ${tumor_sample} \
		--vcf ${tumor_sample}.mutect.vcf \
		--stats-file ${tumor_sample}.mutect.stats.txt \
		--fun-segmentation PSCBS \
		--normaldb ${normaldb_file} \
		--mapping-bias-file ${mapping_bias_file} \
		--intervals ${target_intervals_file} \
		--snp-blacklist ${snp_blacklist_file} \
		--genome hg19 \
		--model betabin \
		--force --post-optimize --seed 123
done