# Use:
# seqc: python3 compare_simulated.py method(0, 1, 2 for needle count, kallisto or salmon) secq_expression num_files data

import numpy as np
import sys
from scipy import stats

# Here it should have an impact that there is a normalization in kallisto and salmon, when using TPM, while needle does not correct for amount, but it is known that coverage is the same in the simulated files
def get_exp_value(line, method):
    if (method == 0):
        return [int(x) for x in line.split()[1:]][0]
    elif (method == 1):
        return  float(line.split()[4])
    elif (method == 2):
        return float(line.split()[3])

method = int(sys.argv[1]) # 0: needle count 1: kallisto 2: salmon
j = 2
dir = sys.argv[j]
num_files = int(sys.argv[j+1])
files = []
for i in range(j+2, j+2+num_files):
    files.append(sys.argv[i])

mse = []
fpr = []
fnr = []
for i in range(0, len(files), 2):
    values_1 = {}
    values_2 = {}
    expected_values = {}
    errors = []
    tp = 0
    tn = 0
    fp = 0
    fn = 0
    with open(dir + "Test_"+str(i+1)+".tsv", 'r') as f:
        for line in f:
            if line[0] != "t":
                transcript = line.split()[0].split('|')[0]
                fold_change = float(line.split()[1])
                expected_values.update({transcript:fold_change})

    with open(files[i], 'r') as f:
        for line in f:
            if (line[0] != "t") & (line[0] != "N"):
                transcript = line.split()[0].split('|')[0]
                exp_list = get_exp_value(line, method)
                values_1.update({transcript:exp_list})


    with open(files[i+1], 'r') as f:
        for line in f:
            if (line[0] != "t") & (line[0] != "N"):
                transcript = line.split()[0].split('|')[0]
                exp_list = get_exp_value(line, method)
                values_2.update({transcript:exp_list})

    for transcript in expected_values:
        if (transcript in values_1) & (transcript in values_2):
            fold_change = (values_1[transcript] + 1)/(values_2[transcript]+ 1) # Log2 drastically improves results of kallisto and salmon, but why?
            errors.append((fold_change-expected_values[transcript]) * (fold_change-expected_values[transcript]))
    mean_square_error = np.mean(errors)
    mse.append(mean_square_error)

    for transcript in values_1:
        if max((values_1[transcript] + 1)/(values_2[transcript]+ 1), (values_2[transcript] + 1)/(values_1[transcript]+ 1)) > 1.15:
            if transcript in expected_values:
                tp += 1
            else:
                fp += 1
        else:
            if transcript in expected_values:
                fn += 1
            else:
                tn += 1
    fpr.append(float(fp)/(fp+tn))
    fnr.append(float(fn)/(fn+tp))

print("Mean Squared error:\n", mse)
print(np.mean(mse), np.var(mse))
print("False positive rate:\n", fpr)
print(np.mean(fpr), np.var(fpr))
print("False negatiive rate:\n", fnr)
print(np.mean(fnr), np.var(fnr))
