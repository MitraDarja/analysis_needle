import numpy as np
import sys

data_set_file = sys.argv[1]
estimate_file = sys.argv[2]
experiment_names_file = sys.argv[3]
type_file = sys.argv[4]

experiment_names = []
with open(experiment_names_file, "r") as f:
    for line in f:
        experiment_names.append(line.strip())

data_set = []
with open(data_set_file, "r") as f:
    for line in f:
        if (line.strip() in experiment_names):
            data_set.append(experiment_names.index(line.strip()))

type_names = []
with open(type_file, "r") as f:
    for line in f:
        tissue = line.split()[0]
        sra_id = line.split()[1]
        if (tissue == "breast"):
            if (sra_id.strip() in experiment_names):
                type_names.append(experiment_names.index(sra_id.strip()))

genes_exp = {}
oncogenes = ["CCND1", "ERBB2","FOXM1",  "MYC"]

with open(estimate_file, "r") as f:
    for line in f:
        exp = [int(a) for a in line.split()[1:]]
        if (line.split("|")[5] in genes_exp):
            genes_exp[line.split("|")[5]].append(np.array(exp))
        else:
            genes_exp.update({line.split("|")[5]:[np.array(exp)]})

zero_exp_values = []
all_exps = []
all_genes = []
for gene in genes_exp:
    cancer_exp = []
    notcancer_exp= []
    exps2 = np.array(genes_exp[gene])
    exps2 = np.mean(exps2, axis = 0)
    all_exps.append(exps2)
    all_genes.append(gene)

all_exps = np.array(all_exps)
mean_exp = np.mean(all_exps, axis = 1)



tps=[]
tns=[]
fps=[]
fns=[]
names = ["True Positives", "False Positives", "True Negatives", "False Negatives"]

for gene in oncogenes:
    tp = 0
    tn = 0
    fp = 0
    fn = 0
    for index in range(1742):
        if (gene in oncogenes):
            if (all_exps[all_genes.index(gene)][index] > mean_exp[index]):
                if index in data_set:
                    tp+=1
                else:
                    fp+=1
            else:
                if index in data_set:
                    fn+=1
                else:
                    tn+=1

    if (gene in oncogenes):
        tps.append(tp)
        tns.append(tn)
        fps.append(fp)
        fns.append(fn)

tps=np.array(tps)
tns=np.array(tns)
fps=np.array(fps)
fns=np.array(fns)


print(tps, fps, tns, fns)
fpr = [fps[i]*1.0/(fps[i]+tns[i]) for i in range(4)]
fnr = [fns[i]*1.0/(fns[i]+tps[i]) for i in range(4)]
print("False Positive Rate:", fpr)
print("False Negative Rate:",fnr)
