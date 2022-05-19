#!/bin/bash
set -e

PATH=/home/lierhan/Program/miniconda3/bin:$PATH
export PURECN="/home/lierhan/Program/miniconda3/lib/R/library/PureCN/extdata"
target_intervals_file="target_intervals.txt"

for normal_sample in SRR948862 SRR948864 SRR948863 SRR948865 SRR948867 SRR948866 SRR948868 SRR948869 SRR948872 SRR948871 SRR948870 SRR948873 SRR948874 SRR948875 SRR948876 SRR948877 SRR948878 SRR948879 SRR948880 SRR948881
do
        normal_bam="${normal_sample}.mkdup.bam"
        nohup Rscript ${PURECN}/Coverage.R --out-dir ${normal_sample} --bam ${normal_bam} --intervals ${target_intervals_file} >${normal_sample}.PureCN.coverage.log 2>&1 &
done

for tumor_sample in 
do
        tumor_bam="${tumor_sample}.mkdup.bam"
        nohup Rscript ${PURECN}/Coverage.R --out-dir ${tumor_sample} --bam ${tumor_bam} --intervals ${target_intervals_file} >${tumor_sample}.PureCN.coverage.log 2>&1 &
done