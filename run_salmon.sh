#!/bin/bash
for letter in "A" "B" "C" "D"; do
    for i in 1 2 3 4; do
        ./src/salmon quant -i gencode.v36.index/ -l ISF -1 ../../Chisanga_data/fqs/SEQC2012-ILM-AGR-${letter}-${i}_R1.fastq.gz -2 ../../Chisanga_data/fqs/SEQC2012-ILM-AGR-${letter}-${i}_R2.fastq.gz -o out/SEQC2012-ILM-AGR-${letter}-${i}.out -g ../../needle/gencode.v36.annotation.gtf
    done
done
