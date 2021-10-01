# Use:
# python3 compare_simulated.py [X] data/ 512 [DATA]
# [X] presents the method to analysis. 1 stands for kallisto, 2 for salmon and 3 for needle estimate or REINDEER.
# [DATA] stands for the output data of the commands in section Analysis Preparations (please input the data in order,
# you can use the bash command ls -v).

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
steps = num_files/4
files = []
if (method != 3):
    for i in range(j+2, j+2+num_files):
        files.append(sys.argv[i])

values = {}
if (method == 3):
    read_needle_estimate(sys.argv[j+2], values)

mse = [[], [], [], []]
files_no = 1
count = 0
mse2 = [[[], [], [], [], []], [[], [], [], [], []], [[], [], [], [], []], [[], [], [], [], []]]
max = 0
max_transcript = ""
max_fc = 0
max_1 = 0
max_2= 0
max_step= 0
for i in range(0, num_files, 2):
    values_1 = {}
    values_2 = {}
    expected_values = {}

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
                if (expected_values[transcript] == 0.5):
                    fc = 1
                elif (expected_values[transcript] < 0.5):
                    fc = 0
                elif (expected_values[transcript] == 2):
                    fc = 2
                elif (expected_values[transcript] == 3):
                    fc = 3
                elif (expected_values[transcript] == 4):
                    fc = 3
                else:
                    fc= 4
                # The distinction here is necessary because kallisto and salmon have some cases, where
                # values_2[transcript]==0 and values_1[transcript] is a bigger number, which then leads to crass outliers,
                # which skew the mean. The comparison is still fair, because the number of elements in error is similar
                # for all three experiments.
                if (method == 0):
                    if (values_2[transcript] >= 0):
                        # + 1, in case of zero values
                        fold_change = (values_1[transcript] + 1)/(values_2[transcript] + 1)
                        mse2[int(i/steps)][fc].append(fold_change)
                        #mse2[int(i/steps)][fc].append((fold_change-expected_values[transcript]) * (fold_change-expected_values[transcript]))
                        count +=1
                else:
                    if (values_2[transcript] > 0):
                        # + 1, in case of values below 1
                        fold_change = (values_1[transcript] + 1)/(values_2[transcript] + 1)
                        mse2[int(i/steps)][fc].append(fold_change)
                        #mse2[int(i/steps)][fc].append((fold_change-expected_values[transcript]) * (fold_change-expected_values[transcript]))
                        count +=1

                        if ((fold_change-expected_values[transcript]) * (fold_change-expected_values[transcript]) > max):
                            max = (fold_change-expected_values[transcript]) * (fold_change-expected_values[transcript])
                            max_transcript = transcript
                            max_fc = expected_values[transcript]
                            max_1 = values_1[transcript]
                            max_2 = values_2[transcript]

    else:
        for transcript in expected_values:
            if (transcript in values):
                if (expected_values[transcript] == 0.5):
                    fc = 1
                elif (expected_values[transcript] < 0.5):
                    fc = 0
                elif (expected_values[transcript] == 2):
                    fc = 2
                elif (expected_values[transcript] == 3):
                    fc = 3
                elif (expected_values[transcript] == 4):
                    fc = 3
                else:
                    fc= 4
                if (values[transcript][i+1] >= 1):
                    # + 1, in case of zero values
                    fold_change = (values[transcript][i])/(values[transcript][i+1])
                    mse2[int(i/steps)][fc].append(fold_change)
                    #mse2[int(i/steps)][fc].append((fold_change-expected_values[transcript]) * (fold_change-expected_values[transcript]))
                    count +=1
                    if (fold_change > max):
                        max = (fold_change-expected_values[transcript]) * (fold_change-expected_values[transcript])
                        max_transcript = transcript
                        max_fc = expected_values[transcript]
                        max_1 = values[transcript][i]
                        max_2 = values[transcript][i+1]
                        max_step = int(i/steps)
                elif ((values[transcript][i+1] == 0) & (values[transcript][i] == 0)):
                    fold_change =1
                    mse2[int(i/steps)][fc].append(fold_change)
                    #mse2[int(i/steps)][fc].append((fold_change-expected_values[transcript]) * (fold_change-expected_values[transcript]))
                    count +=1

print("MAx: ", max, max_transcript, max_fc, max_1, max_2, max_step)

for m in range(4):
    for m2 in range(5):
        print(m, m2, np.mean(mse2[m][m2]), np.var(mse2[m][m2]), len(mse2[m][m2]))

if (method == 1):
    out = "Kallisto_DE"
elif (method == 2):
    out = "Salmon_DE"
elif (method == 3):
    out = "Needle_Reindeer_DE"
np.save(out, mse2)
