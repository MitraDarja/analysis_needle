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

for i in range(4):
    values.update({i:{}})

for file in files:
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
                else:
                    values[it][gene].update({Letter:exp_list})
            else:
                values[it].update({gene:{Letter:exp_list}})

        iterator += 1
        it += 1
        if (it == 4) & (iterator < num_files):
            Letter = Letters[0]
            Letters = Letters[1:]
            it = 0

for it in range(4):
    for gene in values[it]:
        gene_expressions = []
        for l in "ACDB": # Not ABCD, because than it is not in the right order
            exps = np.array(values[it][gene][l])
            gene_expressions.append(np.mean(exps, axis = 0))
        if not (all(gene_expressions[j] <= gene_expressions[j + 1] for j in range(len(gene_expressions)-1))  | all(gene_expressions[j] >= gene_expressions[j + 1] for j in range(len(gene_expressions)-1))):
            miss += 1

print(miss)
