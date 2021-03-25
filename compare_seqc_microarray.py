# Use:
# seqc: python3 compare_seqc_microarray.py method(0, 1, 2 for needle count, kallisto or salmon) 0 secq_expression num_files data
# needle: python3 compare_seqc_microarray.py 0 0 ../Chisanga_data/Taqman-raw.txt 16 ../needle/build/Chisanga_data2/SEQC2012-ILM-AGR-*.out
# kallisto: python3 compare_seqc_microarray.py 1 0 ../Chisanga_data/Taqman-raw.txt 16 ../kallisto/SEQC2012-ILM-AGR-*/abundance.tsv
# salmon: python3 compare_seqc_microarray.py 2 0 ../Chisanga_data/Taqman-raw.txt 16 ../salmon-1.4.0/build/out/SEQC2012-ILM-AGR-*-*.out/quant.sf

# microrray: python3 compare_seqc_microarray.py method(0, 1, 2 for needle count, kallisto or salmon) 1 secq_expression file with entred2gene_id_info num_files data
# kallisto: python3 compare_seqc_microarray.py 1 1 ../Chisanga_data/BeadChip-Log2.table data/entrez_id2gene_id.txt 16 ../kallisto/SEQC2012-ILM-AGR-*/abundance.tsv
# salmon: python3 compare_seqc_microarray.py 2 1 ../Chisanga_data/BeadChip-Log2.table data/entrez_id2gene_id.txt 16 ../salmon-1.4.0/build/out/SEQC2012-ILM-AGR-*-*.out/quant.sf

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

method = int(sys.argv[1]) # 0: needle count 1: kallisto 2: salmon
microrray = int(sys.argv[2]) # 0: seqc 1: microarrayy
j = 3
if microrray:
    seqc_file = sys.argv[j]
    seqc_file2 = sys.argv[j+1]
    num_files = int(sys.argv[j+2])
    files = []
    for i in range(j+3, j+3+num_files):
        files.append(sys.argv[i])
else:
    seqc_file = sys.argv[j]
    num_files = int(sys.argv[j+1])
    files = []
    for i in range(j+2, j+2+num_files):
        files.append(sys.argv[i])

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
            #print(line)
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

gene_lengths = {}
for file in files:
    values = {}
    with open(file, 'r') as f:
        for line in f:
            if (line[0] != "t") & (line[0] != "N"):
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
