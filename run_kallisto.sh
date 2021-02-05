#!/bin/bash
for letter in "A" "B" "C" "D"; do
    for i in 1 2 3 4; do
        ./kallisto quant -i gencode.v36.index -o SEQC2012-ILM-AGR-${letter}-${i}/ ../Chisanga_data/fqs/SEQC2012-ILM-AGR-${letter}-${i}_R1.fastq.gz ../Chisanga_data/fqs/SEQC2012-ILM-AGR-${letter}-${i}_R2.fastq.gz -g ../data/gencode.v36.annotation.gtf
    done
done
