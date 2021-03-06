This is the analysis of the tool Needle (https://github.com/MitraDarja/needle).

# Accuracy Analysis

The purpose of this analysis is to determine the accuracy of needle. In this analysis Needle is compared with the
state-of-the-art tools kallisto and salmon. Moreover, `needle count` is used for a comparison, to check the validity
of determining the expression value of a gene by taking the median of the count values of a transcript's minimizers
(if a gene consits of multiple transcripts the mean of all expression values is taken) without having to account for
errors arising from the probabilistic and approximative structure of Needle (namely the probabilistic data structure of
the interleaved bloom filters and the approximation following the expression levels).
Moreover, REINDEER as a direct competitor of needle was used for the comparison.

## Data
Used was data from the SEQC studies, more precisely, the data provided by Chinsanga et al. (https://github.com/ShiLab-Bioinformatics/GeneAnnotation) was analysed. In order to repeat the analysis, please download this data.

## Analysis Preparations

Please, download and install the tools: kallisto(https://github.com/pachterlab/kallisto), salmon(https://github.com/COMBINE-lab/salmon), REINDEER(https://github.com/kamimrcht/REINDEER) and Needle (https://github.com/MitraDarja/needle).
Then run the following commands.
Time was measured by /usr/bin/time -v. (Tests were run with either 1 thread or 32 threads, add `-t | --threads` to the
commands to do the same.)

### kallisto

```
# Run to create the index of the transcriptome
kallisto index -i gencode.v36.index -k 19  analysis_needle/data/gencode.v36.pc_transcripts.fa.gz

# Run for every input experiment, letter of ["A", "B", "C", "D"] and number of [1, 2, 3, 4]
kallisto quant -i gencode.v36.index -o SEQC2012-ILM-AGR-${letter}-${number}/ SEQC2012-ILM-AGR-${letter}-${i}_R1.fastq.gz SEQC2012-ILM-AGR${letter}-${number}_R2.fastq.gz
```

### salmon

```
# Run to create the index of the transcriptome
salmon index -t analysis_needle/data/gencode.v36.pc_transcripts.fa.gz -i gencode.v36.index -k 19

# Run for every input experiment, letter of ["A", "B", "C", "D"] and number of [1, 2, 3, 4]
salmon quant -i gencode.v36.index/ -l ISF -1 SEQC2012-ILM-AGR-${letter}-${number}_R1.fastq.gz -2 SEQC2012-ILM-AGR-${letter}-${number}_R2.fastq.gz -o out/SEQC2012-ILM-AGR-${letter}-${number}.out
```

### REINDEER

### Needle count

```
needle count -k 19 -w 19 --paired -g analysis_needle/data/gencode.v36.pc_transcripts.fa.gz  $(ls -v SEQC2012-ILM-AGR-*fastq.gz)
```

### Needle estimate

```
# Build the index method 1
needle ibf -k 19 -w 19 -c --paired -g analysis_needle/data/gencode.v36.pc_transcripts.fa.gz  $(ls -v SEQC2012-ILM-AGR-*fastq.gz)
-e 2 -e 4 -e 8 -e 16 -e 32 -e 64 -e 128 -e 256 -e 512 -e 1024 -e 2048 -e 4096 -e 8192 -e 16384 -e 3276 -u 0 -b 100000000 -l 15

# Estimate with method 1
needle estimate -i w_19/ -c -k 19 -w 19 -o w_19/expressions_15.out -e 2 -e 4 -e 8 -e 16 -e 32 -e 64 -e 128 -e 256 -e 512 -e 1024 -e 2048 -e 4096 -e 8192 -e 16384 -e 32768 analysis_needle/data/gencode.v36.pc_transcripts.fa.gz

# Build the index method 2
needle ibf -k 19 -w 19 -c --paired -g analysis_needle/data/gencode.v36.pc_transcripts.fa.gz  $(ls -v SEQC2012-ILM-AGR-*fastq.gz) -u 0 -b 100000000 -l 15

# Estimate with method 2
needle estimate -i w_19/ -c -k 19 -w 19 -o w_19/expressions_levels_15.out -e 0 -e 1 -e 2 -e 3 -e 4 -e 5 -e 6 -e 7 -e 8 -e 9 -e 10 -e 11 -e 12 -e 13 -e 14  analysis_needle/data/gencode.v36.pc_transcripts.fa.gz -d w_19/IBF_Levels.levels
```

If the window size of 23 is used, `-w 23` has to be used.

# Analysis

The spearman correlation with the gene expression according to Taqman QT-PCR and according to the microarray expression
values can be obtained by the following command:

```
# Taqman
python3 compare_seqc_microarray.py [X] 0 Taqman-raw.txt 16 [DATA]

# Microarray
python3 compare_seqc_microarray.py [X] 1 BeadChip-Log2.table data/entrez_id2gene_id.txt 16 [DATA]
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
Time was measured by /usr/bin/time -v. (Tests were run with either 1 thread or 32 threads, add `-t | --threads` to the
commands to do the same.)

## Data

For the analysis, a well known data set containing breast, brain and blood RNA sequence files was used. The complete list
of experiments used by us can be found in `data\samples.in`. Please download these files as fastq.gz files.

## Preprocess

The data was then preprocessed and all k-mers with occurrences below a certain threshold were disregarded. (Thresholds
were based on the file size and follow the recommendation of Salmon et al.)

```
# REINDEER

# Needle

```

## Build

## Query
