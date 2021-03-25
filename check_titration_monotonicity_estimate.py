import numpy as np
import sys
from scipy import stats


expression_file = sys.argv[1]

iterator = 0
Letter = "A"
Letters = ["B", "C", "D"]
it = 0
values = {}
gene_count = 0

for i in range(4):
    values.update({i:{}})

values = {}
with open(expression_file, 'r') as f:
    for line in f:
        gene = line.split()[0].split('|')[5]
        exp_list = [int(x) for x in line.split()[1:]]

        if gene not in values:
            dict = {"A": [exp_list[0:4]], "B": [exp_list[4:8]], "C": [exp_list[8:12]],"D": [exp_list[12:16]] }
            values.update({gene: dict})

        else:
            values[gene]["A"].append(exp_list[0:4])
            values[gene]["B"].append(exp_list[4:8])
            values[gene]["C"].append(exp_list[8:12])
            values[gene]["D"].append(exp_list[12:16])

# Normalization, a normalization should be performed because different experiments are compared to each other
norm_all = {}
norm_all.update({"A":[0,0,0,0], "B":[0,0,0,0], "C":[0,0,0,0], "D":[0,0,0,0]})
for it in range(4):
    for gene in values:
        for l in "ABCD":
            exps = np.array(np.mean(values[gene][l], axis = 0)[it])
            norm_all[l][it] += exps

for i in range(4):
    for letter in "ABCD":
        norm_all[letter][i] = norm_all[letter][i]/1000000.0


# Calculate MSE
error = []
for it in range(4):
    for gene in values:
        gene_expressions = []
        gene_count +=1
        for letter in "ABCD":
            exps = np.array(values[gene][letter])
            gene_expressions.append(np.mean(exps, axis = 0)[it]/norm_all[letter][it])
        # + 1 for dealing with zeros
        expected_fold_change_c_d = np.log2(1+float(gene_expressions[0] + (3*gene_expressions[3]))/((gene_expressions[3] + 1 + (3*gene_expressions[0])))) # A+3B/3A+B = B/A
        fold_change_a_b = np.log2(1+float(gene_expressions[3])/(gene_expressions[0]+1))
        fold_change_c_d = np.log2(1+float(gene_expressions[2])/(gene_expressions[1]+1))
        error.append((expected_fold_change_c_d-fold_change_c_d)* (expected_fold_change_c_d-fold_change_c_d))

print(gene_count)
print(np.mean(error), np.var(error), len(error))
