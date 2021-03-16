import numpy as np
import sys
from scipy import stats


num_files = int(sys.argv[1])
files = []
j = 2
for i in range(j, j+num_files):
    files.append(sys.argv[i])

estimate_file = sys.argv[j+num_files]

pearson = []
spearman = []
estimate = {}
count = 0

# Read estimate file
with open(estimate_file, "r") as f:
    for line in f:
        count +=1
        transcript = line.split()[0].split('|')[0]
        exp_list = [int(x) for x in line.split()[1:]]
        estimate.update({transcript: exp_list})

for i in range(num_files):
    count = []
    estimate_values = []
    with open(files[i], "r") as f:
        for line in f:
            transcript = line.split()[0].split('|')[0]
            exp = int(line.split()[1:][0])
            if estimate[transcript][i] > 0:
                count.append(exp)
                estimate_values.append(estimate[transcript][i])

    if (stats.pearsonr(count, estimate_values)[1] <= 0.005):
        pearson.append(stats.pearsonr(count, estimate_values)[0] )
    else:
        pearson.append(0)
    if (stats.spearmanr(count, estimate_values)[1] <= 0.005):
        spearman.append(stats.spearmanr(count,estimate_values)[0] )
    else:
        spearman.append(0)

print(pearson)
print(np.mean(pearson), np.var(pearson))
print(spearman)
print(np.mean(spearman), np.var(spearman))
