This is the analysis of the tool [Needle](https://github.com/seqan/needle).

# Accuracy Analysis

The purpose of this analysis is to determine the accuracy of needle. In this analysis Needle is compared with the
state-of-the-art tools kallisto and salmon as well as Reindeer.

## Data
Two data sets were used, one simulated and one real data set. The simulated data set can be created by using the provided Rscript `simulate.R`.
Run:
```
Rscript simulate.R data/100.fa
```

The real data set is from the SEQC studies, more precisely, the data provided by Chinsanga et al. (https://github.com/ShiLab-Bioinformatics/GeneAnnotation) was analysed. In order to repeat the analysis, please download this data.

## Analysis Preparations

Please, download and install the tools: [kallisto](https://github.com/pachterlab/kallisto) (v0.46.2), [salmon](https://github.com/COMBINE-lab/salmon)(v1.5.1), [REINDEER](https://github.com/kamimrcht/REINDEER)(v1.0.2) and [Needle](https://github.com/seqan/needle)(v1.0.0).
Then run the following commands.

### kallisto

```
# Run to create the index of the transcriptome for the simulated data
kallisto index -i 100.index -k 19  /data/100.fa

# Run for every input experiment of the simulated data set, i goes from 1 to 256
kallisto quant -i 100.index -o Test-${i}/ /data/Test_${i}/sample_01_1.fasta.gz /data/Test_${i}/sample_01_2.fasta.gz
kallisto quant -i 100.index -o Test-${i}_2/ /data/Test_${i}/sample_02_1.fasta.gz /data/Test_${i}/sample_02_2.fasta.gz

# Run to create the index of the transcriptome
kallisto index -i gencode.v36.index -k 19  /data/gencode.v36.pc_transcripts.fa.gz

# Run for every input experiment of the real data set, letter of ["A", "B", "C", "D"] and number of [1, 2, 3, 4]
kallisto quant -i gencode.v36.index -o SEQC2012-ILM-AGR-${letter}-${number}/ SEQC2012-ILM-AGR-${letter}-${i}_R1.fastq.gz SEQC2012-ILM-AGR${letter}-${number}_R2.fastq.gz
```

### salmon

```
# Run to create the index of the transcriptome for the simulated data
salmon index -t /data/100.fa -i 100.index -k 19

# Run for every input experiment of the simulated data set, i goes from 1 to 256
salmon quant -i 100/ -l ISF -1 /data/Test_${i}/sample_01_1.fasta.gz -2 /data/Test_${i}/sample_01_2.fasta.gz  -o Test_${i}.out
salmon quant -i 100/ -l ISF -1 /data/Test_${i}/sample_02_1.fasta.gz -2 /data/Test_${i}/sample_02_2.fasta.gz  -o Test_${i}_2.out

# Run to create the index of the transcriptome
salmon index -t /data/gencode.v36.pc_transcripts.fa.gz -i gencode.v36.index -k 19

# Run for every input experiment of the real data set, letter of ["A", "B", "C", "D"] and number of [1, 2, 3, 4]
salmon quant -i gencode.v36.index/ -l ISF -1 SEQC2012-ILM-AGR-${letter}-${number}_R1.fastq.gz -2 SEQC2012-ILM-AGR-${letter}-${number}_R2.fastq.gz -o out/SEQC2012-ILM-AGR-${letter}-${number}.out
```

### REINDEER

```
# Run bcalm as preprocess,  i goes from 1 to 256
bcalm -out Test_${i}/sample_01_ -in /data/Test_${i}/sample_01*.fasta.gz -in /data/Test_${i}/sample_01*.fasta.gz  -kmer-size 19 -out-dir bcalm2 -out-tmp bcalm2 -nb-cores 32
bcalm -out Test_${i}/sample_01_ -in /data/Test_${i}/sample_02*.fasta.gz -in /data/Test_${i}/sample_02*.fasta.gz  -kmer-size  19 -out-dir bcalm2 -out-tmp bcalm2 -nb-cores 32

# Construct index, fof_simulated.lst is a file containing all files created by bcalm
reindeer --index -f reindeer/fof_simulated.lst -k 19 -o reindeer/out_simulated

# Query
./Reindeer --query -l reindeer/out_simulated -q data/100.fa -o reindeer/simulated_query


# Run bcalm as preprocess for every input experiment of the real data set, letter of ["A", "B", "C", "D"] and number of [1, 2, 3, 4]
bcalm -in SEQC2012-ILM-AGR-${letter}-${i}_R1.fastq.gz -in SEQC2012-ILM-AGR-${letter}-${i}_R2.fastq.gz -kmer-size 19 -abundance-min 13 -out-dir bcalm2 -out-tmp bcalm2 -nb-cores 32

# Construct index, fof_seqc.lst is a file containing all files created by bcalm
reindeer --index -f reindeer/fof_seqc.lst -k 19 -o reindeer/out_seqc

# Query
./Reindeer --query -l reindeer/out_seqc -q data/gencode.v36.pc_transcripts.fa.gz -o reindeer/seqc_query
```

### Needle

```
# Prerprocess, i goes from 1 to 256
needle minimiser -k 19 -w 19  data/Test_${i}/sample*.fasta.gz -o w_19_Test_${i}_ --paired

# Build index for simulated data
needle ibfmin -o w_19_Simulated_ $(ls -v w_19/Test_*)  -e 5 -e 10 -e 15 -e 20 -e 30 -e 40 -e 60 -e 80 -e 120 -e 160 -e 240 -e 320  -f 0.05

# Query
needle estimate -i w_19_Simulated_ data/100.fa  -o  w_19/expressions_simulated.out

# Preprocess for every input experiment of the real data set
needle minimiser -o w_19_ -k 19 -w 19 -g data/gencode.v36.pc_transcripts.fa.gz SEQC2012-ILM-AGR-*.fastq.gz --paired

# Build index for SEQC data
needle ibfmin -o w_19_SEQC_ $(ls -v w_19_SEQC2012-ILM-AGR-*minimiser) -l 15 -f 0.05

# Query
needle estimate -i w_19_Level_SEQC_ data/gencode.v36.pc_transcripts.fa.gz  -o  w_19_expressions_SEQC.out
```

If the window size of 23 or 39 is used, `-w 23` or `-w 39` has to be used.

# Analysis

## Simulated Data
The difference between the ground truth and the predicted expressions for the differential expression and the coverages on the simulated data set can be obtained by running the following commands:

```
python3 compare_simulated.py [X] data/ 512 [DATA]
python3 compare_simulated_cov.py [X] data/ 512 [DATA]
```

## SEQC Data
[X] presents the method to analysis. 1 stands for kallisto, 2 for salmon and 3 for needle estimate or
REINDEER. [DATA] stands for the output data of the commands in section Analysis Preparations (please input the data in order, you can use the bash command `ls -v`).

The spearman correlation with the gene expression according to Taqman QT-PCR and according to the microarray expression
values can be obtained by the following command:

```
# Taqman
python3 compare_seqc_microarray.py [X] 0 data/Taqman-raw.txt 16 [DATA]

# Microarray
python3 compare_seqc_microarray.py [X] 1 data/BeadChip-Log2.table data/entrez_id2gene_id.txt 16 [DATA]
```

[X] presents the method to analysis. 0 stands for needle count, 1 for kallisto, 2 for salmon and 3 for needle estimate or
REINDEER. [DATA] stands for the output data of the commands in section Analysis Preparations. For needle count those should be
all `*.out` files, for kallisto all abundance.tsv files, for salmon all quant.genes.sf files and for needle estimate one expression.out file. For REINDEER the following command has to be used before the analysis can be run with [DATA] being reindeer.out.

```
python3 reindeer_estimate.py [REINDEER_OUTPUT] data/gencode.v36.pc_transcripts_k_23.fa reindeer.out
```

Note: Taqman-raw.txt and BeadChip-Log2.table were obtained from the provided data from Chinsanga et al.

The mean squared error value in regard of the titration monotoncity can be obtained by the following commands.

```
python3 check_titration_monotonicity.py [X] 16 [DATA]
```

[X] [DATA] are equivalent to the explanation above.

# Space And Speed Analysis

The purpose of this analysis is to determine the space and speed consumption of Needle compared to its main competitor
REINDEER.
Time was measured by /usr/bin/time -v. (Tests were run with either 1 thread, 4 or 32 threads, add `-t | --threads` to the
commands to do the same.)

## Data

For the analysis, a well known data set containing breast, brain and blood RNA sequence files was used. The complete list
of experiments used by us can be found in `data\large_dataset.lst`. Please download these files as fastq.gz files.

## Preprocess

The data was then preprocessed and all k-mers with occurrences below a certain threshold were disregarded. (Thresholds
were based on the file size and follow the recommendation of [Salmon et al.](10.1038/nbt.3442)) The preprocess can be done by running the
scripts `preprocess.sh` and `preprocess_reindeer.sh`.  (Please adjust the paths beforehand.)

```
bash preprocess.sh
preprocess_reindeer.sh

```

## Build & Query

Then adjust the paths in `run_large_dataset.sh` and run:

```
bash run_large_dataset.sh
```

In the folders w_21, w_25, w_41 and reindeer the indexes and their measurement of `/usr/bin/time -v` for the construction and the queries can be found.
