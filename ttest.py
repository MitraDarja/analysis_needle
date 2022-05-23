import csv
import numpy as np
import sys

from scipy import stats

estimate_file = sys.argv[1]
experiment_names_file = sys.argv[2]
data_set_file_breast = sys.argv[3]
data_set_file_blood = sys.argv[4]
data_set_file_brain = sys.argv[5]

# Read all experiments
experiment_names = []
with open(experiment_names_file, "r") as f:
    for line in f:
        experiment_names.append(line.strip())

# Read experiments, which belong to the breast tissue
data_set_breast = []
with open(data_set_file_breast, "r") as f:
    for line in f:
        if (line.split("/")[-1].strip().split(".")[0] in experiment_names):
            data_set_breast.append(experiment_names.index(line.strip()))

# Read experiments, which belong to the blood tissue
data_set_blood = []
with open(data_set_file_blood, "r") as f:
    for line in f:
        if (line.split("/")[-1].strip().split(".")[0] in experiment_names):
            data_set_blood.append(experiment_names.index(line.strip()))

# Read experiments, which belong to the brain tissue
data_set_brain = []
with open(data_set_file_brain, "r") as f:
    for line in f:
        if (line.split("/")[-1].strip().split(".")[0] in experiment_names):
            data_set_brain.append(experiment_names.index(line.strip()))

exps = {}
with open(estimate_file, "r") as f:
    for line in f:
        genename = line.split()[0].split("|")[5]
        exp = [int(a) for a in line.split()[1:]]
        if (genename in exps):
            exps[genename].append(exp)
        else:
            exps.update({genename:[exp]})

corrected_p = 0.05/len(exps)
print(corrected_p)
print(len(data_set_breast), len(data_set_blood), len(data_set_brain))
print(len(exps))
diff_genes = []
diff_genes_blood = []
diff_genes_brain = []
genes_exp = {}
for gene in exps:
    exps2 = np.array(exps[gene])
    exps2 = np.mean(exps2, axis = 0)
    genes_exp.update({gene:exps2})
    in_data_breast = []
    not_data_breast = []
    in_data_blood = []
    not_data_blood = []
    in_data_brain = []
    not_data_brain = []
    for index_exp in range(len(exps2)):
        if (index_exp in data_set_breast):
            in_data_breast.append(exps2[index_exp])
        else:
            not_data_breast.append(exps2[index_exp])

        if (index_exp in data_set_blood):
            in_data_blood.append(exps2[index_exp])
        else:
            not_data_blood.append(exps2[index_exp])
        if (index_exp in data_set_brain):
            in_data_brain.append(exps2[index_exp])
        else:
            not_data_brain.append(exps2[index_exp])

    t_stat, p_val = stats.ttest_ind(in_data_breast, not_data_breast)
    if ((p_val <= corrected_p) & (t_stat > 0) ):
        diff_genes.append(gene)
    t_stat, p_val = stats.ttest_ind(in_data_blood, not_data_blood)
    if ((p_val <= corrected_p) & (t_stat > 0)):
        diff_genes_blood.append(gene)
    t_stat, p_val = stats.ttest_ind(in_data_brain, not_data_brain)
    if ((p_val <= corrected_p) & (t_stat > 0)):
        diff_genes_brain.append(gene)

csv_file = "Gencode.csv"
csv_columns = genes_exp.keys()
with open(csv_file, 'w') as csvfile:
    for genename in genes_exp:
        csvfile.write(genename + "," + ','.join([str(elem) for elem in genes_exp[genename]]) + "\n")


print(len(diff_genes), len(diff_genes_blood), len(diff_genes_brain))

with open("breast.lst", "w") as o:
    for gene in diff_genes:
        o.write(gene+"\n")

with open("blood.lst", "w") as o:
    for gene in diff_genes_blood:
        o.write(gene+"\n")

with open("brain.lst", "w") as o:
    for gene in diff_genes_brain:
        o.write(gene+"\n")
