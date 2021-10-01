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

# Run to create the index of the transcriptome for the simulated data
$kallisto index -i kallisto/100.index -k 19  data/100.fa
$salmon index -i salmon/100.index -k 19 -t data/100.fa

for i in {1..256}; do
    $kallisto quant -i kallisto/100.index -o kallisto/Test-${i}/ data/Test_${i}/sample_01_1.fasta.gz data/Test_${i}/sample_01_2.fasta.gz
	$kallisto quant -i kallisto/100.index -o kallisto/Test-${i}_2/ data/Test_${i}/sample_02_1.fasta.gz data/Test_${i}/sample_02_2.fasta.gz
    $salmon quant -i salmon/100.index/ -l ISF -1 data/Test_${i}/sample_01_1.fasta.gz -2 data/Test_${i}/sample_01_2.fasta.gz  -o salmon/Test_${i}.out
    $salmon quant -i salmon/100.index/ -l ISF -1 data/Test_${i}/sample_02_1.fasta.gz -2 data/Test_${i}/sample_02_2.fasta.gz  -o salmon/Test_${i}_2.out

    # Preprocess Reindeer and Needle
    $bcalm -out reindeer/bcalm2/Test_${i}_sample_01_ -in data/Test_${i}/sample_01_1.fasta.gz -in data/Test_${i}/sample_01_2.fasta.gz  -kmer-size 19 -out-dir reindeer/bcalm2 -out-tmp reindeer/bcalm2 -nb-cores 32
    $bcalm -out reindeer/bcalm2/Test_${i}_sample_02_ -in data/Test_${i}/sample_02_1.fasta.gz -in data/Test_${i}/sample_02_2.fasta.gz  -kmer-size  19 -out-dir reindeer/bcalm2 -out-tmp reindeer/bcalm2 -nb-cores 32
    $needle minimiser -k 19 -w 19  data/Test_${i}/sample*.fasta.gz -o w_19/Test_${i}_ --paired
    $needle minimiser -k 19 -w 23  data/Test_${i}/sample*.fasta.gz -o w_23/Test_${i}_ --paired
    $needle minimiser -k 19 -w 39  data/Test_${i}/sample*.fasta.gz -o w_39/Test_${i}_ --paired
done

# Create Reindeer index
realpath $(ls -v reindeer/bcalm2/Test*unitigs.fa) > reindeer/fof_simulated.lst
$reindeer --index -f reindeer/fof_simulated.lst -k 19 -o reindeer/out_simulated

# Query Reindeer index
$reindeer --query -l reindeer/out_simulated -q data/100.fa -o reindeer/simulated_query
python3 reindeer_estimate.py reindeer/simulated_query/query_results/out_query_Reindeer_P40_100_0.out data/100.fa reindeer/expressions_simulated.out

# Create Needle index
$needle ibfmin -o w_19/Simulated_ $(ls -v w_19/Test_*)  -e 5 -e 10 -e 15 -e 20 -e 30 -e 40 -e 60 -e 80 -e 120 -e 160 -e 240 -e 320  -f 0.05
$needle ibfmin -o w_23/Simulated_ $(ls -v w_23/Test_*)  -e 5 -e 10 -e 15 -e 20 -e 30 -e 40 -e 60 -e 80 -e 120 -e 160 -e 240 -e 320  -f 0.05
$needle ibfmin -o w_39/Simulated_ $(ls -v w_39/Test_*)  -e 5 -e 10 -e 15 -e 20 -e 30 -e 40 -e 60 -e 80 -e 120 -e 160 -e 240 -e 320  -f 0.05

# Query Needle index
$needle estimate -i w_19/Simulated_ data/100.fa  -o  w_19/expressions_simulated.out
$needle estimate -i w_23/Simulated_ data/100.fa  -o  w_23/expressions_simulated.out
$needle estimate -i w_39/Simulated_ data/100.fa  -o  w_39/expressions_simulated.out

# Evaluation

python3 compare_simulated.py 1 data/ 512 $(ls -v kallisto/Test-*/abundance.tsv)
python3 compare_simulated.py 2 data/ 512 $(ls -v salmon/Test_*.out/quant.sf)
python3 compare_simulated.py 3 data/ 512 reindeer/expressions_simulated.out
mv Needle_Reindeer_DE.npy Reindeer.npy
python3 compare_simulated.py 3 data/ 512 w_19/expressions_simulated.out
mv Needle_Reindeer_DE.npy Needle_19.npy
python3 compare_simulated.py 3 data/ 512 w_23/expressions_simulated.out
mv Needle_Reindeer_DE.npy Needle_23.npy
python3 compare_simulated.py 3 data/ 512 w_39/expressions_simulated.out
mv Needle_Reindeer_DE.npy Needle_39.npy

python3 compare_simulated_cov.py 1 data/ 512 $(ls -v kallisto/Test-*/abundance.tsv)
python3 compare_simulated_cov.py 2 data/ 512 $(ls -v salmon/Test_*.out/quant.sf)
python3 compare_simulated_cov.py 3 data/ 512 reindeer/expressions_simulated.out
mv Needle_Reindeer_Cov.npy Reindeer_Cov.npy
python3 compare_simulated_cov.py 3 data/ 512 w_19/expressions_simulated.out
mv Needle_Reindeer_Cov.npy Needle_19_Cov.npy
python3 compare_simulated_cov.py 3 data/ 512 w_23/expressions_simulated.out
mv Needle_Reindeer_Cov.npy Needle_23_Cov.npy
python3 compare_simulated_cov.py 3 data/ 512 w_39/expressions_simulated.out
mv Needle_Reindeer_Cov.npy Needle_39_Cov.npy

python3 boxplot.py
