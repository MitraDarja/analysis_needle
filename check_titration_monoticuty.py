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
fold_change_a_b = []
fold_change_c_d = []
error = []
for it in range(4):
    for gene in values[it]:
        gene_expressions = []
        if len(values[it][gene][l]) > 0:
            gene_count +=1
            for l in "ACDB": # Not ABCD, because than it is not in the right order
                exps = np.array(values[it][gene][l])
                if method == 0:
                    gene_expressions.append(np.mean(exps, axis = 0)/np.mean(norm_all[it][l],axis=0))
                else:
                    gene_expressions.append(np.mean(exps, axis = 0))
            if not (all(gene_expressions[j] <= gene_expressions[j + 1] for j in range(len(gene_expressions)-1))  | all(gene_expressions[j] >= gene_expressions[j + 1] for j in range(len(gene_expressions)-1))):
                miss += 1
            if all(gene_expressions[j] > 0 for j in range(len(gene_expressions))):
                #x = float(gene_expressions[3])/gene_expressions[0]
                #print(np.power(2,x))
                #x = np.power(2,x)
                #print(gene_expressions, x)
                #expected_fold_change_c_d.append(np.log2((3*x)+1) - np.log2(x + 3)) # Following equation 1 from Chisanga et al.
                #fold_change_c_d.append(np.log2(float(gene_expressions[2])/gene_expressions[1])) # Following equation 1 from Chisanga et al.
                #error.append((expected_fold_change_c_d[-1]-fold_change_c_d[-1])* (expected_fold_change_c_d[-1]-fold_change_c_d[-1]))
                expected_fold_change_c_d.append(np.log2(float(gene_expressions[0] + (3*gene_expressions[3]))/((gene_expressions[3] + (3*gene_expressions[0]))))) # A+3B/3A+B = B/A
                fold_change_a_b.append(np.log2(float(gene_expressions[3])/gene_expressions[0]))
                fold_change_c_d.append(np.log2(float(gene_expressions[2])/gene_expressions[1]))
                error.append((expected_fold_change_c_d[-1]-fold_change_c_d[-1])* (expected_fold_change_c_d[-1]-fold_change_c_d[-1]))



print(miss)
print(gene_count)
print(np.mean(error), len(error))
#plt.plot(expected_fold_change_c_d, color="red")
plt.plot(fold_change_a_b, fold_change_c_d, 'o', color = "black")
plt.savefig("check_titration_"+str(method)+".png")
