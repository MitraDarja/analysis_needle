#!/bin/bash

mkdir w_19
mkdir w_23
mkdir w_39
mkdir kallisto
mkdir salmon
mkdir reindeer
mkdir reindeer/bcalm2

kallisto="PATH to executable"
salmon="PATH to executable"
bcalm="PATH to executable"
reindeer="PATH to executable"
needle="PATH to executable"
SEQC_DIR="PATH to SEQC DIR"

# Run to create the index of the transcriptome for the simulated data
$kallisto index -i kallisto/gencode.v36.index -k 19  data/gencode.v36.pc_transcripts.fa.gz
$salmon index -i salmon/gencode.v36.index -k 19 -t data/gencode.v36.pc_transcripts.fa.gz

for letter in "A" "B" "C" "D"; do
    for number in 1 2 3 4; do
        $kallisto quant -i kallisto/gencode.v36.index -o kallisto/SEQC2012-ILM-AGR-${letter}-${number}/ $SEQC_DIR/fqs/SEQC2012-ILM-AGR-${letter}-${number}_R1.fastq.gz $SEQC_DIR/fqs/SEQC2012-ILM-AGR-${letter}-${number}_R2.fastq.gz
        $salmon quant -i salmon/gencode.v36.index/ -l ISF -1 $SEQC_DIR/fqs/SEQC2012-ILM-AGR-${letter}-${number}_R1.fastq.gz -2 $SEQC_DIR/fqs/SEQC2012-ILM-AGR-${letter}-${number}_R2.fastq.gz -o salmon/SEQC2012-ILM-AGR-${letter}-${number}.out

        # Preprocess Reindeer
        $bcalm -out reindeer/bcalm2/SEQC_${letter}_${number}_ -in $SEQC_DIR/fqs/SEQC2012-ILM-AGR-${letter}-${number}_R1.fastq.gz -in $SEQC_DIR/fqs/SEQC2012-ILM-AGR-${letter}-${number}_R2.fastq.gz -kmer-size 19 -abundance-min 13 -out-dir reindeer/bcalm2 -out-tmp reindeer/bcalm2 -nb-cores 32
    done
done

# Preprocess Needle
$needle minimiser -o w_19/ -k 19 -w 19 --include data/gencode.v36.pc_transcripts.fa.gz $SEQC_DIR/fqs/SEQC2012-ILM-AGR-*.fastq.gz --paired
$needle minimiser -o w_23/ -k 19 -w 23 --include data/gencode.v36.pc_transcripts.fa.gz $SEQC_DIR/fqs/SEQC2012-ILM-AGR-*.fastq.gz --paired
$needle minimiser -o w_39/ -k 19 -w 39 --include data/gencode.v36.pc_transcripts.fa.gz $SEQC_DIR/fqs/SEQC2012-ILM-AGR-*.fastq.gz --paired

# Create Reindeer index
realpath $(ls -v reindeer/bcalm2/SEQC*unitigs.fa) > reindeer/fof_seqc.lst
$reindeer --index -f reindeer/fof_seqc.lst -k 19 -o reindeer/out_seqc

# Query Reindeer index
$reindeer --query -l reindeer/out_seqc -q data/gencode.v36.pc_transcripts.fa.gz -o reindeer/seqc_query
gunzip data/gencode.v36.pc_transcripts_k_23.fa.gz
python3 reindeer_estimate.py reindeer/seqc_query/query_results/out_query_Reindeer_P40_gencode_0.out data/gencode.v36.pc_transcripts_k_23.fa reindeer/expressions_seqc.out

# Create Needle index
$needle needle ibfmin -o w_19/SEQC_data/gencode.v36.pc_transcripts.fa.gz -l 15 -f 0.05
$needle needle ibfmin -o w_23/SEQC_ $(ls -v w_23/SEQC2012-ILM-AGR-*minimiser) -l 15 -f 0.05
$needle needle ibfmin -o w_39/SEQC_ $(ls -v w_39/SEQC2012-ILM-AGR-*minimiser) -l 15 -f 0.05

# Query Needle index
$needle estimate -i w_19/SEQC_ data/gencode.v36.pc_transcripts.fa.gz  -o  w_19/expressions_SEQC.out
$needle estimate -i w_23/SEQC_ data/gencode.v36.pc_transcripts.fa.gz  -o  w_23/expressions_SEQC.out
$needle estimate -i w_39/SEQC_ data/gencode.v36.pc_transcripts.fa.gz  -o  w_39/expressions_SEQC.out

# Query Needle with norm
$needle estimate -i w_19/SEQC_ data/gencode.v36.pc_transcripts.fa.gz  -o  w_19/expressions_SEQC_norm.out -m
$needle estimate -i w_23/SEQC_ data/gencode.v36.pc_transcripts.fa.gz  -o  w_23/expressions_SEQC_norm.out -m
$needle estimate -i w_39/SEQC_ data/gencode.v36.pc_transcripts.fa.gz  -o  w_39/expressions_SEQC_norm.out -m

# Evaluation
python3 summary_seqcdata.py $SEQC_DIR

# Create Needle count
$needle genome -k 19 -w 19 -o w_19/gencode data/gencode.v36.pc_transcripts.fa.gz
$needle genome -k 19 -w 23 -o w_23/gencode data/gencode.v36.pc_transcripts.fa.gz
$needle genome -k 19 -w 39 -o w_39/gencode data/gencode.v36.pc_transcripts.fa.gz

$needle count --genome w_19/gencode.genome -k 19 -w 19 --include data/gencode.v36.pc_transcripts.fa.gz $SEQC_DIR/fqs/SEQC2012-ILM-AGR-*.fastq.gz -o w_19/ --paired
$needle count --genome w_23/gencode.genome -k 19 -w 23 --include data/gencode.v36.pc_transcripts.fa.gz $SEQC_DIR/fqs/SEQC2012-ILM-AGR-*.fastq.gz -o w_23/ --paired
$needle count --genome w_39/gencode.genome -k 19 -w 39 --include data/gencode.v36.pc_transcripts.fa.gz $SEQC_DIR/fqs/SEQC2012-ILM-AGR-*.fastq.gz -o w_39/ --paired

python3 summary_seqcdata_count.py $SEQC_DIR
