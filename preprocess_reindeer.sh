#!/bin/bash  
# Based on https://github.com/kamimrcht/REINDEER/blob/master/reproduce_manuscript_results/bcalm_2585.sh

bcalm="PATH_TO_BCALM"

mkdir -p bcalm_results
mkdir -p tmp
input_file=$1
#get fastq.gz and launch bcalm on each file
while read -r filename threshold; do
    #filename = $(cut -d'\t' -f1 <<<< $line)
    #threshold = $(cut -d'\t' -f2 <<<< $line)
    $bcalm -in /project/archive-index-data/rnaseq_benchmark_real_data/$filename -kmer-size 21 -abundance-min "$threshold" -out-dir bcalm2 -out-tmp bcalm2 -nb-cores 32
    mv $filename.unitigs.fa bcalm_results
    ls $PWD/$filename.unitigs.fa >> fof_reindeer.lst
done < $input_file
