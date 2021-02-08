import matplotlib.pyplot as plt
import numpy as np
import sys
from scipy import stats


def get_exp_value(line, method):
    if (method == 0):
        return [int(x) for x in line.split()[1:]]
    elif (method == 1):
        return  [float(line.split()[4])]
    elif (method == 2):
        return [float(line.split()[3])]

method = int(sys.argv[1]) # 0: needle count 1: kallisto 2: salmon
seqc_file = sys.argv[2]
num_files = int(sys.argv[3])
files = []
for i in range(4, 4+num_files):
    files.append(sys.argv[i])

seqc_values = {}
iterator = 0
Letter = "A"
Letters = ["B", "C", "D"]
pearson = []
spearman = []
it = 0
miss = 0
with open(seqc_file, 'r') as f:
    for line in f:
        if line[0] != "E":
            line_table = line.split()
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

gene_lengths = {}
for file in files:
    values = {}
    with open(file, 'r') as f:
        for line in f:
            if line[0] != "t":
                gene = line.split()[0].split('|')[5]
                exp_list = get_exp_value(line, method)
                length = int(line.split()[1])
                if gene in values:
                    values[gene].append(exp_list[0])
                    gene_lengths[gene].append(length)
                else:
                    values.update({gene:exp_list})
                    gene_lengths.update({gene:[length]})

    seqc = []
    expressions = []
    for gene in seqc_values:
        if gene in values:
            exps = np.array(values[gene])
            expressions.append(np.mean(exps, axis = 0))
            exps2 = np.array(seqc_values[gene][Letter])
            seqc.append((np.mean(exps2, axis = 0)[it]))
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

    print (miss, len(seqc), len(expressions))
    miss = 0
    iterator += 1
    it += 1
    if (it == 4) & (iterator < num_files):
        Letter = Letters[0]
        Letters = Letters[1:]
        it = 0



print(pearson)
print(np.mean(pearson), np.var(pearson))
print(spearman)
print(np.mean(spearman), np.var(spearman))

spearman_means = []
spearman = np.array(spearman)
before = 0
for i in range(4,20,4):
    spearman_means.append(np.mean(spearman[before:i]))
    before = i
print(spearman_means)
print(np.mean(spearman_means), np.var(spearman_means))
