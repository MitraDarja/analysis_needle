# Use:
# seqc: python3 compare_simulated.py method(0, 1, 2 for needle count, kallisto or salmon) secq_expression num_files data
# python3 compare_simulated.py 0 data/ 40 $(ls -v analysis_needle_simulation/Test_*.out)
# python3 compare_simulated.py 1 data/ 40 $(ls -v ../kallisto/Test-*/abundance.tsv)
# python3 compare_simulated.py 2 data/ 40 $(ls -v ../salmon-1.4.0/build/out/Test_*.out/quant.sf)

import numpy as np
import sys
from scipy import stats

# Here two experiments are compared, but the coverage difference is known, so unnormalized data can be used.
def get_exp_value(line, method):
    if (method == 0):
        return [int(x) for x in line.split()[1:]][0]
    elif (method == 1):
        return  float(line.split()[3])
    elif (method == 2):
        return float(line.split()[4])

def read_needle_estimate(estimate_file, estimate):
    # Read estimate file
    with open(estimate_file, "r") as f:
        for line in f:
            if line[0] == ">":
                line = line[1:]
            gene =  line.split()[0].split('|')[0]
            exp_list = [float(x) for x in line.split()[1:]]
            estimate.update({gene:exp_list})

method = int(sys.argv[1]) # 0: needle count 1: kallisto 2: salmon 3: needle estimate/reindeer
j = 2
dir = sys.argv[j]
num_files = int(sys.argv[j+1])
files = []
if (method != 3):
    for i in range(j+2, j+2+num_files):
        files.append(sys.argv[i])

values = {}
if (method == 3):
    read_needle_estimate(sys.argv[j+2], values)

mse = []
files_no = 1
count = 0
for i in range(0, num_files, 2):
    values_1 = {}
    values_2 = {}
    expected_values = {}
    errors = []
    max = 0
    max_transcript = ""
    with open(dir + "Test_"+str(files_no)+"/sim_tx_info.txt", 'r') as f:
        for line in f:
            if line[0] != "t":
                transcript = line.split()[0].split('|')[0]
                fold_change = float(line.split()[1])
                expected_values.update({transcript:fold_change})
    files_no +=1

    if (method != 3):
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
                values_1[transcript] = values_1[transcript]
                values_2[transcript] = values_2[transcript]
                # The distinction here is necessary because kallisto and salmon have some cases, where
                # values_2[transcript]==0 and values_1[transcript] is a bigger number, which then leads to crass outliers,
                # which skew the mean. The comparison is still fair, because the number of elements in error is similar
                # for all three experiments.
                if (method == 0):
                    if (values_2[transcript] >= 0):
                        # + 1, in case of zero values
                        fold_change = (values_1[transcript] + 1)/(values_2[transcript] + 1)
                        errors.append((fold_change-expected_values[transcript]) * (fold_change-expected_values[transcript]))
                        count +=1
                else:
                    if (values_2[transcript] > 0):
                        # + 1, in case of values below 1
                        fold_change = (values_1[transcript] + 1)/(values_2[transcript] + 1)
                        errors.append((fold_change-expected_values[transcript]) * (fold_change-expected_values[transcript]))
                        count +=1
    else:
        for transcript in expected_values:
            if (transcript in values):
                if (values[transcript][i+1] >= 1):
                    # + 1, in case of zero values
                    fold_change = (values[transcript][i])/(values[transcript][i+1])
                    errors.append((fold_change-expected_values[transcript]) * (fold_change-expected_values[transcript]))
                    count +=1
                elif ((values[transcript][i+1] == 0) & (values[transcript][i] == 0)):
                    fold_change =1
                    errors.append((fold_change-expected_values[transcript]) * (fold_change-expected_values[transcript]))
                    count +=1
    mean_square_error = np.mean(errors)
    mse.append(mean_square_error)

print("Mean Squared error of fold change:")
print(np.mean(mse), np.var(mse), count)

mse = []
count = 0
expected_value = 1 # normalization should have "fixed" different coverages
for i in range(1, 32, 2):
    values_1 = {}
    values_2 = {}
    values_3 = {}
    values_4 = {}
    errors = []

    if (method != 3):
        with open(files[i], 'r') as f:
            for line in f:
                if (line[0] != "t") & (line[0] != "N"):
                    transcript = line.split()[0].split('|')[0]
                    exp_list = get_exp_value(line, method)
                    values_1.update({transcript:exp_list})

        with open(files[i+32], 'r') as f:
            for line in f:
                if (line[0] != "t") & (line[0] != "N"):
                    transcript = line.split()[0].split('|')[0]
                    exp_list = get_exp_value(line, method)
                    values_2.update({transcript:exp_list})

        with open(files[i+64], 'r') as f:
            for line in f:
                if (line[0] != "t") & (line[0] != "N"):
                    transcript = line.split()[0].split('|')[0]
                    exp_list = get_exp_value(line, method)
                    values_3.update({transcript:exp_list})

        with open(files[i+96], 'r') as f:
            for line in f:
                if (line[0] != "t") & (line[0] != "N"):
                    transcript = line.split()[0].split('|')[0]
                    exp_list = get_exp_value(line, method)
                    values_4.update({transcript:exp_list})

        for transcript in expected_values:
            if (transcript in values_1) & (transcript in values_2):
                values_1[transcript] = values_1[transcript]
                values_2[transcript] = values_2[transcript]

                expected_value = 2
                if ((values_1[transcript] >= 1) & (values_2[transcript] >= 1)):
                        fold_change = values_2[transcript]/values_1[transcript]
                        errors.append((fold_change-expected_value) * (fold_change-expected_value))
                        count +=1
                expected_value = 3
                if ((values_1[transcript] >= 1) & (values_3[transcript] >= 1)):
                        fold_change = values_3[transcript]/values_1[transcript]
                        errors.append((fold_change-expected_value) * (fold_change-expected_value))
                        count +=1
                expected_value = 4
                if ((values_1[transcript] >= 1) & (values_4[transcript] >= 1)):
                        fold_change = values_4[transcript]/values_1[transcript]
                        errors.append((fold_change-expected_value) * (fold_change-expected_value))
                        count +=1

                expected_value = 2
                if ((values_2[transcript] >= 1) & (values_3[transcript] >= 1)):
                        fold_change = values_3[transcript]/values_2[transcript]
                        errors.append((fold_change-expected_value) * (fold_change-expected_value))
                        count +=1
                expected_value = 3
                if ((values_2[transcript] >= 1) & (values_4[transcript] >= 1)):
                        fold_change = values_4[transcript]/values_2[transcript]
                        errors.append((fold_change-expected_value) * (fold_change-expected_value))
                        count +=1
                expected_value = 4
                if ((values_3[transcript] >= 1) & (values_4[transcript] >= 1)):
                        fold_change = values_4[transcript]/values_3[transcript]
                        errors.append((fold_change-expected_value) * (fold_change-expected_value))
                        count +=1
    else:
        for transcript in expected_values:
            if (transcript in values):
                for j in range(3):
                    for l in range(j+1,4):
                        expected_value = l + 1
                        if ((values[transcript][i + (32*j)] >= 1) & (values[transcript][i+(32*l)] >= 0)):
                                fold_change = values[transcript][i + (32*l)]/values[transcript][i+(32*j)]
                                errors.append((fold_change-expected_value) * (fold_change-expected_value))
                                count +=1
                        elif ((values[transcript][i + (32*j)] == 0) & (values[transcript][i+(32*l)] == 0)):
                                fold_change = 1
                                errors.append((fold_change-expected_value) * (fold_change-expected_value))
                                count +=1

    mean_square_error = np.mean(errors)
    mse.append(mean_square_error)

print("Mean Squared error of coverage:\n", mse)
print(np.mean(mse), np.var(mse), count)
