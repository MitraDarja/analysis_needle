import numpy as np
import sys
from scipy import stats


# TPM values of kallisto and salmon are used, because different experiments are compared to each other
def get_exp_value(line, method):
    if (method == 0):
        return [int(x) for x in line.split()[1:]]
    elif (method == 1):
        return  [float(line.split()[4])]
    elif (method == 2):
        return [float(line.split()[3])]

method = int(sys.argv[1]) # 0: needle count 1: kallisto 2: salmon 3: needle estimate or REINDEER
num_files = int(sys.argv[2])

iterator = 0
Letter = "A"
Letters = ["B", "C", "D"]
it = 0
values = {}
gene_count = 0

for i in range(4):
    values.update({i:{}})

if (method != 3):
    files = []
    for i in range(3, 3+num_files):
        files.append(sys.argv[i])

    # Read in values for every file
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
else:
    values = {}
    with open(sys.argv[3], 'r') as f:
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
if (method == 0) | (method == 3):
    norm_all.update({"A":[0,0,0,0], "B":[0,0,0,0], "C":[0,0,0,0], "D":[0,0,0,0]})
    for it in range(4):
        if (method == 0):
            for gene in values[it]:
                for l in "ABCD":
                    exps = np.array(values[it][gene][l])
                    norm_all[l][it] += np.mean(exps, axis = 0)
        else:
            for gene in values:
                for l in "ABCD":
                    exps = np.array(np.mean(values[gene][l], axis = 0)[it])
                    norm_all[l][it] += exps
    for i in range(4):
        for letter in "ABCD":
            norm_all[letter][i] = norm_all[letter][i]/1000000.0
else:
    norm_all.update({"A":[1,1,1,1], "B":[1,1,1,1], "C":[1,1,1,1], "D":[1,1,1,1]})

# Calculate MSE
error = []
for it in range(4):
    if (method != 3):
        for gene in values[it]:
            gene_expressions = []
            gene_count +=1
            for letter in "ABCD":
                exps = np.array(values[it][gene][letter])
                gene_expressions.append(np.mean(exps, axis = 0)/norm_all[letter][it])
            # + 1 for dealing with zeros
            expected_fold_change_c_d = np.log2(1+float(gene_expressions[0] + (3*gene_expressions[3]))/((gene_expressions[3] + 1 + (3*gene_expressions[0])))) # A+3B/3A+B = B/A
            fold_change_a_b = np.log2(1+float(gene_expressions[3])/(gene_expressions[0]+1))
            fold_change_c_d = np.log2(1+float(gene_expressions[2])/(gene_expressions[1]+1))
            error.append((expected_fold_change_c_d-fold_change_c_d)* (expected_fold_change_c_d-fold_change_c_d))
    else:
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
