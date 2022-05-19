#!/bin/bash
set -e

java_7=""	# mutect-1.1.7.jar need java 1.7
java_8=""	# GenomeAnalysisTK-3.8 need java 1.8
mutect="mutect-1.1.7.jar"
gatk="GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar"
reference_genome="hg19.fa"
dbsnp_file="dbsnp_138.hg19.vcf"
cosmic_file="hg19_CosmicCombinedMuts_v77.vcf"

# 1. Calling variants for all normal samples by artifact_detection_mode and select variants
for normal_sample in SRR948862 SRR948864 SRR948863 SRR948865 SRR948867 SRR948866 SRR948868 SRR948869 SRR948872 SRR948871 SRR948870 SRR948873 SRR948874 SRR948875 SRR948876 SRR948877 SRR948878 SRR948879 SRR948880 SRR948881
do
	normal_bam="${normal_sample}.mkdup.bam"
	mutect_stats_file="{normal_sample}.mutect.stats.txt"
	mutect_vcf_file="{normal_sample}.mutect.vcf"
	mutect_pon_vcf_file="${normal_sample}.mutect_pon.vcf"
	# Calling variants for all normal samples by artifact_detection_mode
	${java_7} -Xmx4g -jar ${mutect} --analysis_type MuTect --reference_sequence ${reference_genome} --artifact_detection_mode --dbsnp ${dbsnp_file} --cosmic ${cosmic_file} -dt None --input_file:tumor ${normal_bam} --out ${mutect_stats_file} --vcf ${mutect_vcf_file}
	# select variants
	${java_8} -jar -Xmx4g ${gatk} --analysis_type SelectVariants -R ${reference_genome} -V ${mutect_vcf_file} -o ${mutect_pon_vcf_file}
done

# 2. Combine variants
${java_8} -Xmx24g -jar ${gatk} \
  -T CombineVariants \
  -nt 4 --minimumN 5 --genotypemergeoption UNSORTED \
  -R ${reference_genome} \
  -V SRR948862.mutect_pon.vcf \
  -V SRR948864.mutect_pon.vcf \
  -V SRR948863.mutect_pon.vcf \
  -V SRR948865.mutect_pon.vcf \
  -V SRR948867.mutect_pon.vcf \
  -V SRR948866.mutect_pon.vcf \
  -V SRR948868.mutect_pon.vcf \
  -V SRR948869.mutect_pon.vcf \
  -V SRR948872.mutect_pon.vcf \
  -V SRR948871.mutect_pon.vcf \
  -V SRR948870.mutect_pon.vcf \
  -V SRR948873.mutect_pon.vcf \
  -V SRR948874.mutect_pon.vcf \
  -V SRR948875.mutect_pon.vcf \
  -V SRR948876.mutect_pon.vcf \
  -V SRR948877.mutect_pon.vcf \
  -V SRR948878.mutect_pon.vcf \
  -V SRR948879.mutect_pon.vcf \
  -V SRR948880.mutect_pon.vcf \
  -V SRR948881.mutect_pon.vcf \
  -o normals.merged.min5.vcf

bgzip normals.merged.min5.vcf
tabix -p vcf normals.merged.min5.vcf.gz