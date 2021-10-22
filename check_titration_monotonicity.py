# Use:
# python3 check_titration_monotonicity.py [X] 16 [DATA]
# [X] presents the method to analysis. 1 for kallisto, 2 for salmon and 3 for needle estimate
# or REINDEER.
# [DATA] stands for the output data of the commands in section Analysis Preparations. For kallisto those should be all abundance.tsv files,
# for salmon all quant.genes.sf files and for needle estimate one expression.out file.
# For REINDEER the script `reindeer_estimate.py` has to be used before the analysis can be run with [DATA] being reindeer.out.

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

def get_mse(method, num_files, files):
    iterator = 0
    Letter = "A"
    Letters = ["B", "C", "D"]
    it = 0
    values = {}
    gene_count = 0

    for i in range(4):
        values.update({i:{}})

    if (method != 3):
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
        with open(files[0], 'r') as f:
            for line in f:
                gene = line.split()[0].split('|')[5]
                exp_list = [float(x) for x in line.split()[1:]]

                if gene not in values:
                    dict = {"A": [exp_list[0:4]], "B": [exp_list[4:8]], "C": [exp_list[8:12]],"D": [exp_list[12:16]] }
                    values.update({gene: dict})

                else:
                    values[gene]["A"].append(exp_list[0:4])
                    values[gene]["B"].append(exp_list[4:8])
                    values[gene]["C"].append(exp_list[8:12])
                    values[gene]["D"].append(exp_list[12:16])

    # Calculate MSE
    error = []
    for it in range(4):
        if (method != 3):
            for gene in values[it]:
                gene_expressions = []
                gene_count +=1
                for letter in "ABCD":
                    exps = np.array(values[it][gene][letter])
                    gene_expressions.append(np.mean(exps, axis = 0))
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
                    gene_expressions.append(np.mean(exps, axis = 0)[it])
                # + 1 for dealing with zeros
                expected_fold_change_c_d = np.log2(1+float(gene_expressions[0] + (3*gene_expressions[3]))/((gene_expressions[3] + 1 + (3*gene_expressions[0])))) # A+3B/3A+B = B/A
                fold_change_a_b = np.log2(1+float(gene_expressions[3])/(gene_expressions[0]+1))
                fold_change_c_d = np.log2(1+float(gene_expressions[2])/(gene_expressions[1]+1))

                error.append((expected_fold_change_c_d-fold_change_c_d)* (expected_fold_change_c_d-fold_change_c_d))

    return(gene_count, round(np.mean(error), 1), round(np.var(error), 2))

if __name__ == "__main__":
    method = int(sys.argv[1]) # 0: needle count 1: kallisto 2: salmon 3: needle estimate or REINDEER
    num_files = int(sys.argv[2])
    files = []
    if (method != 3):
        files = []
        for i in range(3, 3+num_files):
            files.append(sys.argv[i])
    else:
        files.append(sys.argv[3])


    results = get_mse(method, num_files, files)
    print(results[0])
    print(results[1:])
