This is the analysis of the tool [Needle](https://github.com/seqan/needle).


# Download Data & Tools

Please, download and install the tools: [kallisto](https://github.com/pachterlab/kallisto) (v0.46.2), [salmon](https://github.com/COMBINE-lab/salmon)(v1.5.1), [REINDEER](https://github.com/kamimrcht/REINDEER)(v1.0.2) and [Needle](https://github.com/seqan/needle)(v1.0.1).
Then run the following commands.

For the accuracy analysis, you need the following two data sets: one data set is a simulated data set, the other one a real data set.
The simulated data set can be created by using the provided Rscript `simulate.R`.
Run:
```
Rscript simulate.R data/100.fa
```

The real data set is from the SEQC studies, more precisely, the data provided by Chinsanga et al. (https://github.com/ShiLab-Bioinformatics/GeneAnnotation) was analysed. In order to repeat the analysis, please download this data.

For the speed and space analysis, a well known data set containing breast, brain and blood RNA sequence files was used. The complete list
of experiments used by us can be found in `data\large_dataset.lst`. Please download these files as fastq.gz files.

# Accuracy Analysis

The purpose of this analysis is to determine the accuracy of needle. In this analysis Needle is compared with the
state-of-the-art tools kallisto and salmon as well as Reindeer.

## Simulated Data
The difference between the ground truth and the predicted expressions for the differential expression and the coverages
on the simulated data set can be obtained by running the following command. Please adaptt the paths in `run_simulated.sh` beforehand.

```
bash run_simulated.sh
```
This creates two Boxplots: `Boxplot.png` and `Boxplot_Cov.png` showing the estimated fold changes of the differential expression and the different coverages.
(It also creates `Boxplot_03.png` and `Boxplot_Cov_03.png` for a FPR of 0.3.)

## SEQC Data
The correlation values with the expression values of the RT-PCR and the microarray and the mean square error of the fold changes for all tools can be obtained
by running the following command:
```
bash run_seqc.sh
```

The results are written in `Results_SEQC_Data_set.txt` (for a FPR of 0.3 in `Results_SEQC_Data_set.txt`).

# Space And Speed Analysis

The purpose of this analysis is to determine the space and speed consumption of Needle compared to its main competitor
REINDEER.
Time was measured by /usr/bin/time -v. (Tests were run with either 1 thread, 4 or 32 threads, add `-t | --threads` to the
commands to do the same.)


## Preprocess

The data was then preprocessed and all k-mers with occurrences below a certain threshold were disregarded. (Thresholds
were based on the file size and follow the recommendation of [Salmon et al.](10.1038/nbt.3442)) The preprocess can be done by running the
scripts `preprocess.sh` and `preprocess_reindeer.sh`.  (Please adjust the paths beforehand.)

```
bash preprocess.sh 
preprocess_reindeer.sh data/samples.in

```

## Build & Query

Then adjust the paths in `run_large_dataset.sh` and run:

```
bash run_large_dataset.sh
```

In the folders w_21, w_25, w_41 and reindeer the indexes and their measurement of `/usr/bin/time -v` for the construction and the queries can be found.
A summary of the results can be found in the file `Results_Large_Data_set.txt` containing a summary of the times, ram usages and sizes of the build and the query.
If you want the same analysis, run the scripts with same names and the `_thread4` extension.

# Using quantification for differential gene expression analysis

For this analysis, we need to know the names of the RNA-seq experiments and need to know, which tissue they are from. This information is provided in `data/sras_1742.lst`, `data/brain_1742.lst`, `data/breast_1742.lst` and `data/blood_1742.lst`.
Moreover, we need to estimate the expression of all proteincoding genes on a needle index, we picked here one we created under the space and speed analysis.

```
./needle estimate -i w_41/SRR_ data/gencode.v36.pc_transcripts.fa.gz -o gencode.out
```

Then you can run the following command to get all over expressed genes between one type of tissues compared to the other types according to t-test:
```
python3 ttest.py gencode.out data/sras_1742.lst data/breast_1742.lst data/blood_1742.lst data/brain_1742.lst
```

You can find the genes then in `brain.lst`, `breast.lst` and `blood.lst`, which you can then analyse with ShinyGo or any other gene ontology tool of your choice.

# Extracting experiments based on quantification

For this analysis, we copied the names of the breast cancer experiments from [REINDEER](https://github.com/kamimrcht/REINDEER/blob/master/reproduce_manuscript_results/data/cancer_dataset) to our repo (`data/cancer_dataset.lst`) and copied the list of experiments and their tissue of origin form [here](https://www.cs.cmu.edu/%7Eckingsf/software/bloomtree/srr-list.txt) to our repo (`data/sras_types.lst`). With this information, we can check how well a distinction based on the quantification of the oncogenes `CCND1`, `ERBB2`,`FOXM1` and `MYC` works by running the following command:

```
python3 oncogene_analysis.py data/cancer_dataset.lst gencode.out data/sras_1742.lst data/sras_types.lst
```
