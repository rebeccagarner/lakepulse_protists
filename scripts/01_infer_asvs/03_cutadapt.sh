#!/bin/bash
cutadapt --version  # 3.1
for sample in $(cat samples)
do
	echo "On sample: $sample"
	cutadapt -a GGCTTAATTTGACTCAACRCG...ATAACAGGTCTGTGATGCCC -A GGGCATCACAGACCTGTTAT...CGYGTTGAGTCAAATTAAGCC -m 160 -M 230 --discard-untrimmed -o ${sample}_R1_trimmed.fastq.gz -p ${sample}_R2_trimmed.fastq.gz ${sample}_R1.fastq.gz ${sample}_R2.fastq.gz >> lp2017-2019watercolumn18s_cutadapt_stats.txt 2>&1
done
