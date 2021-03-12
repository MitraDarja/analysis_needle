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
estimate = []
count = 0

# Read estimate file
with open(estimate_file, "r") as f:
    for line in f:
        count +=1
        exp_list = [int(x) for x in line.split()[1:]]
        estimate.append(exp_list)

estimate = np.array(estimate)

for i in range(num_files):
    count = []
    with open(files[i], "r") as f:
        for line in f:
            exp = int(line.split()[1:][0])
            count.append(exp)

    if (stats.pearsonr(count, estimate[:,i])[1] <= 0.005):
        pearson.append(stats.pearsonr(count, estimate[:,i])[0] )
    else:
        pearson.append(0)
    if (stats.spearmanr(count, estimate[:,i])[1] <= 0.005):
        spearman.append(stats.spearmanr(count,estimate[:,i])[0] )
    else:
        spearman.append(0)

print(pearson)
print(np.mean(pearson), np.var(pearson))
print(spearman)
print(np.mean(spearman), np.var(spearman))
