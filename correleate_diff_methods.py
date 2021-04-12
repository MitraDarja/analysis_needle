# Calculates the pearson and spearman correlation between two experiments.
# Usage: python3 correlation_diff_methods.py method_1 method_2 num_files files_method_1 files_method_2

import numpy as np
import sys
from scipy import stats

# No normalization of sequencing data is necessary, because one looks only at one data set at a time
def get_exp_value(line, method):
    if (method == 0):
        return [int(x) for x in line.split()[1:]]
    elif (method == 1):
        return  [float(line.split()[3])]
    elif (method == 2):
        return [float(line.split()[4])]

# Needle does not estimate where a sequence is mostly likely coming from, so if transcripts share some sequence,
# reads overlaping this shared sequence are considered for multiple transcripts. That is why, the mean should be taken.
# Kallisto and salmmon on the other hand, approximate the most likely transcript a read is coming from and therefore
# do not count reads more than once, therefore sum is the appropiate method.
def get_transcript_exp(expressions, method, it):
    if (method == 0):
        return np.mean(expressions, axis = 0)
    elif (method != 4):
        return np.sum(expressions, axis = 0)
    else:
        return np.mean(expressions, axis = 0)[it]


def read_file(file, values, method):
    with open(file, 'r') as f:
        for line in f:
            if (line[0] != "t") & (line[0] != "N"):
                gene = line.split()[0].split('|')[1]
                exp_list = get_exp_value(line, method)
                if gene in values:
                    values[gene].append(exp_list[0])
                else:
                    values.update({gene:exp_list})

def read_needle_estimate(estimate_file, estimate):
    # Read estimate file
    with open(estimate_file, "r") as f:
        for line in f:
            gene = line.split()[0].split('|')[1]
            exp_list = [int(x) for x in line.split()[1:]]
            if gene in estimate:
                estimate[gene].append(exp_list)
            else:
                estimate.update({gene:[exp_list]})

method_1 = int(sys.argv[1]) # 0: needle count 1: kallisto 2: salmon, 4: needle estimate
method_2 = int(sys.argv[2]) # 0: needle count 1: kallisto 2: salmon, 4: needle estimate
method1_values = {}
method2_values = {}
j = 3
num_files = int(sys.argv[j])
files_1 = []
if (method_1 != 4):
    for i in range(j+1, j+1+num_files):
        files_1.append(sys.argv[i])
    j += 1+num_files
else:
    read_needle_estimate(sys.argv[j+1], method1_values)
    j += 2

files_2 = []
if (method_2 != 4):
    for i in range(j, j+num_files):
        files_2.append(sys.argv[i])
else:
    read_needle_estimate(sys.argv[j], method2_values)

pearson = []
spearman = []
miss = 0
for i in range(num_files):
    if (method_1 != 4):
        method1_values = {}
        read_file(files_1[i], method1_values, method_1)
    if (method_2 != 4):
        method2_values = {}
        read_file(files_2[i], method2_values, method_2)

    method_1_exps = []
    method_2_exps = []
    for gene in method1_values:
        if gene in method2_values:
            exps = np.array(method1_values[gene])
            method_1_exps.append(get_transcript_exp(exps, method_1, i))
            exps2 = np.array(method2_values[gene])
            method_2_exps.append(get_transcript_exp(exps2, method_2, i))
        else:
            miss += 1

    if (stats.pearsonr(method_1_exps, method_2_exps)[1] <= 0.003):
        pearson.append(stats.pearsonr(method_1_exps, method_2_exps)[0] )
    else:
        pearson.append(0)
    if (stats.spearmanr(method_1_exps, method_2_exps)[1] <= 0.003):
        spearman.append(stats.spearmanr(method_1_exps, method_2_exps)[0] )
    else:
        spearman.append(0)

    print (miss, len(method_1_exps), len(method_2_exps))
    miss = 0

print(pearson)
print(np.mean(pearson), np.var(pearson))
print(spearman)
print(np.mean(spearman), np.var(spearman))
