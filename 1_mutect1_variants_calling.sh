#!/bin/bash
set -e

java_7=""	# mutect-1.1.7.jar need java 1.7
mutect="mutect-1.1.7.jar"
reference_genome="hg19.fa"
dbsnp_file="dbsnp_138.hg19.vcf"
cosmic_file="hg19_CosmicCombinedMuts_v77.vcf"

for tumor_sample in SRR948953 SRR948954 SRR948956 SRR948955 SRR948957 SRR948958
do
	tumor_bam="${tumor_sample}.mkdup.bam"
	mutect_stats_file="{tumor_sample}.mutect.stats.txt"
	mutect_vcf_file="{tumor_sample}.mutect.vcf"
	${java} -Xmx4g -jar ${mutect} --analysis_type MuTect --reference_sequence ${reference_genome} --dbsnp ${dbsnp_file} --cosmic ${cosmic_file} --input_file:tumor ${tumor_bam} --out ${mutect_stats_file} --vcf ${mutect_vcf_file}
done