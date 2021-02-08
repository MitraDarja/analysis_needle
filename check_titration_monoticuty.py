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
num_files = int(sys.argv[2])
files = []
for i in range(3, 3+num_files):
    files.append(sys.argv[i])

iterator = 0
Letter = "A"
Letters = ["B", "C", "D"]
it = 0
values = {}
miss = 0
gene_count = 0

for i in range(4):
    values.update({i:{}})

for file in files:
    gene_lengths = {}
    with open(file, 'r') as f:
        for line in f:
            if (line[0] == "N") | (line[0] == "t"):
                continue
            gene =  line.split()[0].split('|')[5]
            exp_list = get_exp_value(line, method)
            length = int(line.split('|')[6])
            if gene in values[it]:
                if Letter in values[it][gene]:
                    values[it][gene][Letter].append(exp_list[0])
                    gene_lengths[gene].append(length)
                else:
                    values[it][gene].update({Letter:exp_list})
                    gene_lengths.update({gene:[length]})
            else:
                values[it].update({gene:{Letter:exp_list}})
                gene_lengths.update({gene:[length]})

        iterator += 1
        it += 1
        if (it == 4) & (iterator < num_files):
            Letter = Letters[0]
            Letters = Letters[1:]
            it = 0

norm_all = {}
for i in range(4):
    norm_all.update({i:{"A":[0,0,0,0], "B":[0,0,0,0], "C":[0,0,0,0], "D":[0,0,0,0]}})
for it in range(4):
    for gene in values[it]:
        gene_count +=1
        gene_expressions = []
        for l in "ACDB": # Not ABCD, because than it is not in the right order
            exps = np.array(values[it][gene][l])
            norm_all[it][l] += np.mean(exps, axis = 0)/np.mean(gene_lengths[gene])

expected_fold_change_c_d  = []
fold_change_c_d = []
for it in range(4):
    for gene in values[it]:
        gene_count +=1
        gene_expressions = []
        for l in "ACDB": # Not ABCD, because than it is not in the right order
            exps = np.array(values[it][gene][l])
            if method == 0:
                gene_expressions.append(np.mean(exps, axis = 0)/np.mean(norm_all[it][l],axis=0))
            else:
                gene_expressions.append(np.mean(exps, axis = 0))
        if not (all(gene_expressions[j] <= gene_expressions[j + 1] for j in range(len(gene_expressions)-1))  | all(gene_expressions[j] >= gene_expressions[j + 1] for j in range(len(gene_expressions)-1))):
            miss += 1
        if all(gene_expressions[j] > 0 for j in range(len(gene_expressions))):
            x = float(gene_expressions[3])/gene_expressions[0]
            #print(gene_expressions, x)
            expected_fold_change_c_d.append(np.log2(((3*x)+1)/(x + 3))) # Following equation 1 from Chisanga et al.
            fold_change_c_d.append(float(gene_expressions[2])/gene_expressions[1]) # Following equation 1 from Chisanga et al.



print(miss)
print(gene_count, gene_count/4)
plt.plot(expected_fold_change_c_d, fold_change_c_d)
plt.savefig("check_titration_"+str(method)+".png")
