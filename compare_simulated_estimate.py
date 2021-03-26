import numpy as np
import sys
from scipy import stats

expression_file = sys.argv[1]
dir = sys.argv[2]
num_files = int(sys.argv[3])

expressions = {}
with open(expression_file, 'r') as f:
    for line in f:
        transcript = line.split()[0].split('|')[0]
        exp_list = [int(x) for x in line.split()[1:]]
        expressions.update({transcript:exp_list})

mse = []
files_no = 1
count = 0
for i in range(0, num_files, 2):
    values_1 = {}
    values_2 = {}
    expected_values = {}
    errors = []
    per_million_1 = 0
    per_million_2 = 0
    max = 0
    max_transcript = ""
    with open(dir + "Test_"+str(files_no)+"/sim_tx_info.txt", 'r') as f:
        for line in f:
            if line[0] != "t":
                transcript = line.split()[0].split('|')[0]
                fold_change = float(line.split()[1])
                expected_values.update({transcript:fold_change})
    files_no +=1
    for transcript in expressions:
        values_1.update({transcript: expressions[transcript][i]})
        values_2.update({transcript: expressions[transcript][i+1]})
        per_million_1 += expressions[transcript][i]
        per_million_2 += expressions[transcript][i+1]

    # TPM factor for both experiments
    per_million_1 = per_million_1/1000000.0
    per_million_2 = per_million_2/1000000.0
    for transcript in expected_values:
        if (transcript in values_1) & (transcript in values_2):
            values_1[transcript] = values_1[transcript]/per_million_1
            values_2[transcript] = values_2[transcript]/per_million_2
            # The distinction here is necessary because kallisto and salmon have some cases, where
            # values_2[transcript]==0 and values_1[transcript] is a bigger number, which then leads to crass outliers,
            # which skew the mean. The comparison is still fair, because the number of elements in error is similar
            # for all three experiments.
            if (values_2[transcript] >= 0):
                # + 1, in case of zero values
                fold_change = (values_1[transcript] + 1)/(values_2[transcript] + 1)
                errors.append((fold_change-expected_values[transcript]) * (fold_change-expected_values[transcript]))
                count +=1

    mean_square_error = np.mean(errors)
    mse.append(mean_square_error)


print("Mean Squared error:\n", mse)
print(np.mean(mse), np.var(mse), count)
