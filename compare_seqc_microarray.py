# Use:
# Taqman
# python3 compare_seqc_microarray.py [X] 0 data/Taqman-raw.txt 16 [DATA]

# Microarray
# python3 compare_seqc_microarray.py [X] 1 data/BeadChip-Log2.table data/entrez_id2gene_id.txt 16 [DATA]

# [X] presents the method to analysis. 1 stands for kallisto, 2 for salmon and 3 for needle estimate or REINDEER.
# [DATA] stands for the output data of the commands in section Analysis Preparations
# (please input the data in order, you can use the bash command ls -v).]

import numpy as np
import sys
from scipy import stats

# TPM values are used for kallisto and salmon in order to have a length correction.
def get_exp_value(line, method):
    if (method == 0):
        return [int(x) for x in line.split()[1:]]
    elif (method == 1):
        return  [float(line.split()[4])]
    elif (method == 2):
        return [float(line.split()[3])]

# Needle does not estimate where a sequence is mostly likely coming from, so if transcripts share some sequence,
# reads overlaping this shared sequence are considered for multiple transcripts. That is why, the mean should be taken.
# Kallisto and salmmon on the other hand, approximate the most likely transcript a read is coming from and therefore
# do not count reads more than once, therefore sum is the appropiate method.
def get_transcript_exp(expressions, method, it):
    if (method == 0):
        return np.mean(expressions, axis = 0)
    elif (method != 3):
        return np.sum(expressions, axis = 0)
    else:
        return np.mean(expressions, axis = 0)[it]

def read_file(file, values, method):
    with open(file, 'r') as f:
        for line in f:
            if (line[0] != "t") & (line[0] != "N"):
                gene = line.split()[0].split('|')[5]
                exp_list = get_exp_value(line, method)
                if gene in values:
                    values[gene].append(exp_list[0])
                else:
                    values.update({gene:exp_list})

def read_needle_estimate(estimate_file, estimate):
    # Read estimate file
    with open(estimate_file, "r") as f:
        for line in f:
            gene = line.split()[0].split('|')[5]
            exp_list = [float(x) for x in line.split()[1:]]
            if gene in estimate:
                estimate[gene].append(exp_list)
            else:
                estimate.update({gene:[exp_list]})

def get_exp_mean(method, microrray, seqc_file, seqc_file2, files, num_files):
    seqc_values = {}
    iterator = 0
    Letter = "A"
    Letters = ["B", "C", "D"]
    pearson = []
    spearman = []
    it = 0
    miss = 0

    # Get Gene ids for microarray data
    gene_ids = {}
    if microrray:
        with open(seqc_file2, 'r') as f:
            for line in f:
                if line[0] != "q":
                    line_table = line.split()
                    gene = line_table[1]
                    if gene not in gene_ids:
                        gene_ids.update({line_table[0]:gene})

    with open(seqc_file, 'r') as f:
        for line in f:
            if line[0] != "E":
                line_table = line.split()
                if microrray:
                     if line_table[0] in gene_ids:
                         gene =  gene_ids[line_table[0]]
                         exp_list_A = [float(x) for x in [line_table[5], line_table[9], line_table[13], line_table[17]]]
                         exp_list_B = [float(x) for x in [line_table[6], line_table[10], line_table[14], line_table[18]]]
                         exp_list_C = [float(x) for x in [line_table[7], line_table[11], line_table[15], line_table[19]]]
                         exp_list_D = [float(x) for x in [line_table[8], line_table[12], line_table[16], line_table[20]]]
                     else:
                        continue
                else:
                    gene = line_table[1]
                    exp_list_A = [float(x) for x in [line_table[2], line_table[4], line_table[6], line_table[8]]]
                    exp_list_B = [float(x) for x in [line_table[10], line_table[12], line_table[14], line_table[16]]]
                    exp_list_C = [float(x) for x in [line_table[18], line_table[20], line_table[22], line_table[24]]]
                    exp_list_D = [float(x) for x in [line_table[26], line_table[28], line_table[30], line_table[32]]]

                if gene not in seqc_values:
                    dict = {"A": [exp_list_A], "B": [exp_list_B], "C": [exp_list_C],"D": [exp_list_D] }
                    seqc_values.update({gene: dict})
                else:
                    seqc_values[gene]["A"].append(exp_list_A)
                    seqc_values[gene]["B"].append(exp_list_B)
                    seqc_values[gene]["C"].append(exp_list_C)
                    seqc_values[gene]["D"].append(exp_list_D)

    values = {}
    if (method == 3):
        read_needle_estimate(files[0], values)

    counti = 0
    for f in range(num_files):
        if (method != 3):
            values = {}
            read_file(files[f], values, method)

        seqc = []
        expressions = []
        for gene in seqc_values:
            if gene in values:
                exps = np.array(values[gene])
                exps_value = get_transcript_exp(exps, method, f)
                exps2 = np.array(seqc_values[gene][Letter])
                if ((np.mean(exps2, axis = 0)[it]) >= 0.01):
                    expressions.append(exps_value)
                    seqc.append((np.mean(exps2, axis = 0)[it]))
                    counti = counti +1
                else:
                    miss +=1
            else:
                miss += 1

        if (stats.pearsonr(seqc, expressions)[1] <= 0.003):
            pearson.append(stats.pearsonr(seqc, expressions)[0] )
        else:
            pearson.append(0)
        if (stats.spearmanr(seqc, expressions)[1] <= 0.003):
            spearman.append(stats.spearmanr(seqc, expressions)[0] )
        else:
            spearman.append(0)

        miss = 0
        iterator += 1
        it += 1
        if (it == 4) & (iterator < num_files):
            Letter = Letters[0]
            Letters = Letters[1:]
            it = 0
    print(counti, miss)
    print(np.mean(pearson), np.var(pearson))
    return (round(np.mean(spearman),3), round(np.var(spearman),3))

if __name__ == "__main__":
    method = int(sys.argv[1]) # 0: needle count 1: kallisto 2: salmon, 3: needle estimate or reindeer after running reinder_estimate.py
    microrray = int(sys.argv[2]) # 0: seqc 1: microarrayy
    j = 3
    files = []
    seqc_file2 = ""
    if microrray:
        seqc_file = sys.argv[j]
        seqc_file2 = sys.argv[j+1]
        num_files = int(sys.argv[j+2])
        j += 3
    else:
        seqc_file = sys.argv[j]
        num_files = int(sys.argv[j+1])
        j += 2

    if (method != 3):
        for i in range(j, j+num_files):
            files.append(sys.argv[i])
    else:
        files.append(sys.argv[j])

    results = get_exp_mean(method, microrray, seqc_file, seqc_file2, files, num_files)
    print(results)
