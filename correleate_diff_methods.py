
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

method_1 = int(sys.argv[1]) # 0: needle count 1: kallisto 2: salmon
method_2 = int(sys.argv[2]) # 0: needle count 1: kallisto 2: salmon
j = 3
num_files = int(sys.argv[j])
files_1 = []
for i in range(j+1, j+1+num_files):
    files_1.append(sys.argv[i])

j += 1+num_files
files_2 = []
for i in range(j, j+num_files):
    files_2.append(sys.argv[i])

pearson = []
spearman = []
miss = 0

gene_lengths = {}
for i in range(num_files):
    method1_values = {}
    method2_values = {}
    with open(files_1[i], 'r') as f:
        for line in f:
            if (line[0] != "t") & (line[0] != "N"):
                gene = line.split()[0].split('|')[0]
                exp_list = get_exp_value(line, method_1)
                length = int(line.split()[1])
                if gene in method1_values:
                    method1_values[gene].append(exp_list[0])
                    gene_lengths[gene].append(length)
                else:
                    method1_values.update({gene:exp_list})
                    gene_lengths.update({gene:[length]})

    with open(files_2[i], 'r') as f:
        for line in f:
            if (line[0] != "t") & (line[0] != "N"):
                gene = line.split()[0].split('|')[0]
                exp_list = get_exp_value(line, method_2)
                if gene in method2_values:
                    method2_values[gene].append(exp_list[0])
                else:
                    method2_values.update({gene:exp_list})

    method_1_exps = []
    method_2_exps = []
    for gene in method1_values:
        if gene in method2_values:
            exps = np.array(method1_values[gene])
            method_1_exps.append(np.mean(exps, axis = 0))
            exps2 = np.array(method2_values[gene])
            method_2_exps.append(np.mean(exps2, axis = 0))
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
