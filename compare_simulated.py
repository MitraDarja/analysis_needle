# Use:
# seqc: python3 compare_simulated.py method(0, 1, 2 for needle count, kallisto or salmon) secq_expression num_files data
# python3 compare_simulated.py 0 data/ 20 $(ls -v analysis_needle_simulation/Test_*)
# python3 compare_simulated.py 1 data/ 20 $(ls -v ../kallisto/Test-*/abundance.tsv)
# python3 compare_simulated.py 2 data/ 20 $(ls -v ../salmon-1.4.0/build/out/Test_*.out/quant.sf)

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
for i in range(0, len(files), 2):
    values_1 = {}
    values_2 = {}
    expected_values = {}
    errors = []
    per_million_1 = 0
    per_million_2 = 0
    count = 0
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
                per_million_1 += exp_list


    with open(files[i+1], 'r') as f:
        for line in f:
            if (line[0] != "t") & (line[0] != "N"):
                transcript = line.split()[0].split('|')[0]
                exp_list = get_exp_value(line, method)
                values_2.update({transcript:exp_list})
                per_million_2 += exp_list

    if (method == 0):
        per_million_1 = per_million_1/1000000.0
        per_million_2 = per_million_2/1000000.0
    else:
        per_million_1 = 1
        per_million_2 = 1
    for transcript in expected_values:
        if (transcript in values_1) & (transcript in values_2):
            count +=1
            values_1[transcript] = values_1[transcript]/per_million_1
            values_2[transcript] = values_2[transcript]/per_million_2
            fold_change = (values_1[transcript] + 1)/(values_2[transcript]+ 1) # Log2 drastically improves results of kallisto and salmon, but why?
            errors.append((fold_change-expected_values[transcript]) * (fold_change-expected_values[transcript]))
        else:
            print(transcript, (transcript in values_1), (transcript in values_2))
    mean_square_error = np.mean(errors)
    mse.append(mean_square_error)
    print(count)


print("Mean Squared error:\n", mse)
print(np.mean(mse), np.var(mse))
