import numpy as np
import sys
from scipy import stats


def get_exp_value(line, method):
    if (method == 0):
        return [int(x) for x in line.split()[1:]]
    else if (method == 1):
        return  [float(line.split()[4])]
    else if (method == 2):
        return [float(line.split()[3])]

num_files = int(sys.argv[1])
method = int(sys.argv[2]) # 0: needle count 1: kallisto 2: salmon
files = []
for i in range(3, 3+num_files):
    files.append(sys.argv[i])

iterator = 0
Letter = "A"
Letters = ["B", "C", "D"]
it = 0
values = {}
miss = 0

for i in range(1,5):
    values.update({i:{})



for file in files:
    with open(file, 'r') as f:
        for line in f:
            gene = line.split()[0].split('|')[5]
            exp_list = get_exp_value(line, method)
            length = int(line.split('|')[6])
            if gene in values[it]:
                values[it][gene][letter].append(exp_list[0])
            else:
                values[it].update({gene:{letter:exp_list}})
        iterator += 1
        it += 1
        if (it == 4) & (iterator < num_needle_files):
            Letter = Letters[0]
            Letters = Letters[1:]
            it = 0


for it in range(1,5):
    for gene in values[it][l]:
        gene_expressions = []
        for l in "ACDB": # Not ABCD, because than it is not in the right order
            exps = np.array(values[it][gene][l])
            gene_expressions.append(np.mean(exps, axis = 0))
        if !(all(gene_expressions[j] <= gene_expressions[j + 1] for j in range(len(gene_expressions)-1)) || all(gene_expressions[j] >= gene_expressions[j + 1] for j in range(len(gene_expressions)-1))):
            miss += 1

print(miss)
