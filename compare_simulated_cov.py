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
steps = int(num_files/4)
files = []
if (method != 3):
    for i in range(j+2, j+2+num_files):
        files.append(sys.argv[i])

values = {}
if (method == 3):
    read_needle_estimate(sys.argv[j+2], values)

expected_values = {}
errors = []
files_no = 1
with open(dir + "Test_"+str(files_no)+"/sim_tx_info.txt", 'r') as f:
    for line in f:
        if line[0] != "t":
            transcript = line.split()[0].split('|')[0]
            fold_change = float(line.split()[1])
            expected_values.update({transcript:fold_change})

results_20_40 = []
results_20_60 = []
results_20_80 = []
results_40_60 = []
results_40_80 = []
results_60_80 = []
results = [results_20_40, results_20_60, results_20_80, results_40_60, results_40_80, results_60_80]
expecteds = [2, 3, 4, 1.5, 2, 1.33]
count = 0
for i in range(1, 128, 2):
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

        with open(files[i+steps], 'r') as f:
            for line in f:
                if (line[0] != "t") & (line[0] != "N"):
                    transcript = line.split()[0].split('|')[0]
                    exp_list = get_exp_value(line, method)
                    values_2.update({transcript:exp_list})

        with open(files[i+(steps*2)], 'r') as f:
            for line in f:
                if (line[0] != "t") & (line[0] != "N"):
                    transcript = line.split()[0].split('|')[0]
                    exp_list = get_exp_value(line, method)
                    values_3.update({transcript:exp_list})

        with open(files[i+(steps*3)], 'r') as f:
            for line in f:
                if (line[0] != "t") & (line[0] != "N"):
                    transcript = line.split()[0].split('|')[0]
                    exp_list = get_exp_value(line, method)
                    values_4.update({transcript:exp_list})

        for transcript in expected_values:
            if (transcript in values_1) & (transcript in values_2):
                if ((values_1[transcript] >= 1) & (values_2[transcript] >= 1)):
                        fold_change = values_2[transcript]/values_1[transcript]
                        #results[0].append((fold_change-expecteds[0]) * (fold_change-expecteds[0]))
                        results[0].append(fold_change)
                        count +=1
                if ((values_1[transcript] >= 1) & (values_3[transcript] >= 1)):
                        fold_change = values_3[transcript]/values_1[transcript]
                        #results[1].append((fold_change-expecteds[1]) * (fold_change-expecteds[1]))
                        results[1].append(fold_change)
                        count +=1
                if ((values_1[transcript] >= 1) & (values_4[transcript] >= 1)):
                        fold_change = values_4[transcript]/values_1[transcript]
                        #results[2].append((fold_change-expecteds[2]) * (fold_change-expecteds[2]))
                        results[2].append(fold_change)
                        count +=1
                if ((values_2[transcript] >= 1) & (values_3[transcript] >= 1)):
                        fold_change = values_3[transcript]/values_2[transcript]
                        #results[3].append((fold_change-expecteds[3]) * (fold_change-expecteds[3]))
                        results[3].append(fold_change)
                        count +=1
                if ((values_2[transcript] >= 1) & (values_4[transcript] >= 1)):
                        fold_change = values_4[transcript]/values_2[transcript]
                        #results[4].append((fold_change-expecteds[4]) * (fold_change-expecteds[4]))
                        results[4].append(fold_change)
                        count +=1
                if ((values_3[transcript] >= 1) & (values_4[transcript] >= 1)):
                        fold_change = values_4[transcript]/values_3[transcript]
                        #results[5].append((fold_change-expecteds[5]) * (fold_change-expecteds[5]))
                        results[5].append(fold_change)
                        count +=1
    else:
        for transcript in expected_values:
            if (transcript in values):
                it = 0
                for j in range(3):
                    for l in range(j+1,4):
                        expected_value = l + 1
                        if ((values[transcript][i + int(steps*j)] >= 1) & (values[transcript][i+int(steps*l)] >= 0)):
                                fold_change = values[transcript][i + int(steps*l)]/values[transcript][i+int(steps*j)]
                                results[it].append(fold_change)
                                #results[it].append((fold_change-expecteds[it]) * (fold_change-expecteds[it]))
                                count +=1
                        elif ((values[transcript][i + int(steps*j)] == 0) & (values[transcript][i+int(steps*l)] == 0)):
                                fold_change = 1
                                #results[it].append((fold_change-expecteds[it]) * (fold_change-expecteds[it]))
                                results[it].append(fold_change)
                                count +=1
                        it+=1

for r in range(len(results)):
    print(np.mean(results[r]), np.var(results[r]))

if (method == 1):
    out = "Kallisto_Cov"
elif (method == 2):
    out = "Salmon_Cov"
elif (method == 3):
    out = "Needle_Reindeer_Cov"
np.save(out, results)
