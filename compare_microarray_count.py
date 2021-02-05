import numpy as np
import sys
from scipy import stats


seqc_file = sys.argv[1]
seqc_file2 = sys.argv[2]
num_needle_files = int(sys.argv[3])
files = []
for i in range(4, 4+num_needle_files):
    files.append(sys.argv[i])


seqc_values = {}
iterator = 0
Letter = "A"
Letters = ["B", "C", "D"]
pearson = []
spearman = []
it = 0
miss = 0
gene_ids = {}
with open(seqc_file2, 'r') as f:
    for line in f:
        #print(line)
        if line[0] != "q":
            line_table = line.split(',')
            gene = line_table[3]
            if gene not in gene_ids:
                gene_ids.update({line_table[0]:gene})

with open(seqc_file, 'r') as f:
    for line in f:
        if line[0] != "E":
            line_table = line.split()
            if line_table[0] in gene_ids:
                gene = gene_ids[line_table[0]]
                exp_list_A = [float(x) for x in [line_table[5], line_table[9], line_table[13], line_table[17]]]
                exp_list_B = [float(x) for x in [line_table[6], line_table[10], line_table[14], line_table[18]]]
                exp_list_C = [float(x) for x in [line_table[7], line_table[11], line_table[15], line_table[19]]]
                exp_list_D = [float(x) for x in [line_table[8], line_table[12], line_table[16], line_table[20]]]
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
    needle_values = {}
    with open(file, 'r') as f:
        for line in f:
            gene = line.split()[0].split('|')[5]
            exp_list = [int(x) for x in line.split()[1:]]
            length = int(line.split('|')[6])
            if gene in needle_values:
                needle_values[gene].append(exp_list[0])
                gene_lengths[gene].append(length)
            else:
                needle_values.update({gene:exp_list})
                gene_lengths.update({gene:[length]})

    seqc = []
    needle = []
    for gene in seqc_values:
        if gene in needle_values:
            exps = np.array(needle_values[gene])
            needle.append(np.mean(exps, axis = 0) ) # Chisana uses + 0.5
            exps2 = np.array(seqc_values[gene][Letter])
            # How Chisanga did it
            #means = [np.mean(exps2, axis = 1)]
            #pos_max = means.index(max(means))
            #seqc.append(np.log2(exps2[pos_max][it]+1)) # /np.mean(gene_lengths[gene], axis = 0)
            #seqc.append(np.log2(np.mean(exps2, axis = 0)[it]/np.mean(gene_lengths[gene], axis = 0) +1))
            seqc.append(np.mean(exps2, axis = 0)[it])
        else:
            miss += 1

    if (stats.pearsonr(seqc, needle)[1] <= 0.003):
        pearson.append(stats.pearsonr(seqc, needle)[0] )
    else:
        pearson.append(0)
    if (stats.spearmanr(seqc, needle)[1] <= 0.003):
        spearman.append(stats.spearmanr(seqc, needle)[0] )
    else:
        spearman.append(0)

    print (miss, len(seqc), len(needle))
    miss = 0
    iterator += 1
    it += 1
    if (it == 4) & (iterator < num_needle_files):
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
