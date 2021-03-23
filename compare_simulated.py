# Use:
# seqc: python3 compare_simulated.py method(0, 1, 2 for needle count, kallisto or salmon) secq_expression num_files data
# python3 compare_simulated.py 0 data/ 40 $(ls -v analysis_needle_simulation/Test_*.out)
# python3 compare_simulated.py 1 data/ 40 $(ls -v ../kallisto/Test-*/abundance.tsv)
# python3 compare_simulated.py 2 data/ 40 $(ls -v ../salmon-1.4.0/build/out/Test_*.out/quant.sf)

import numpy as np
import sys
from scipy import stats

# Here two experiments are compared, so a normalization like TPM should be used. kallisto and salmon provide a TPM
# normalization and return the TPM value for every transcript, but the transcript file provided contained more
# transcript than we will check during one comparison, therefore for kallisto and salmon not the TPM value is considered
# but the estimated counts. A TPM-normalization is done for all three methods later in this script.
def get_exp_value(line, method):
    if (method == 0):
        return [int(x) for x in line.split()[1:]][0]
    elif (method == 1):
        return  float(line.split()[3])
    elif (method == 2):
        return float(line.split()[4])

method = int(sys.argv[1]) # 0: needle count 1: kallisto 2: salmon
j = 2
dir = sys.argv[j]
num_files = int(sys.argv[j+1])
files = []
for i in range(j+2, j+2+num_files):
    files.append(sys.argv[i])

mse = []
files_no = 1
for i in range(0, len(files), 2):
    values_1 = {}
    values_2 = {}
    expected_values = {}
    errors = []
    per_million_1 = 0
    per_million_2 = 0
    max = 0
    max_transcript = ""
    count = 0
    with open(dir + "Test_"+str(files_no)+".tsv", 'r') as f:
        for line in f:
            if line[0] != "t":
                transcript = line.split()[0].split('|')[0]
                fold_change = float(line.split()[1])
                expected_values.update({transcript:fold_change})
    files_no +=1

    with open(files[i], 'r') as f:
        for line in f:
            if (line[0] != "t") & (line[0] != "N"):
                transcript = line.split()[0].split('|')[0]
                exp_list = get_exp_value(line, method)
                values_1.update({transcript:exp_list})
                if (method == 0):
                    per_million_1 += exp_list
                # Kallisto and salmon need transcript length correction, for needle done by taking median
                else:
                    per_million_1 += exp_list/int(line.split('|')[6])


    with open(files[i+1], 'r') as f:
        for line in f:
            if (line[0] != "t") & (line[0] != "N"):
                transcript = line.split()[0].split('|')[0]
                exp_list = get_exp_value(line, method)
                values_2.update({transcript:exp_list})
                if (method == 0):
                    per_million_2 += exp_list
                # Kallisto and salmon need transcript length correction, for needle done by taking median
                else:
                    per_million_2 += exp_list/int(line.split('|')[6])

    # TPM factor for both files
    per_million_1 = per_million_1/1000000.0
    per_million_2 = per_million_2/1000000.0
    for transcript in expected_values:
        if (transcript in values_1) & (transcript in values_2):
            # Rounding of values for kallisto, which has sometimes really small values, which lead to
            # nonsensical fold changes, because there are unreasonably large
            values_1[transcript] = round(values_1[transcript], 2)/per_million_1
            values_2[transcript] = round(values_2[transcript], 2)/per_million_2
            # + 1, in case of zero values
            fold_change = (values_1[transcript] + 1)/(values_2[transcript] + 1)
            errors.append((fold_change-expected_values[transcript]) * (fold_change-expected_values[transcript]))
            count +=1
            if ((fold_change-expected_values[transcript]) * (fold_change-expected_values[transcript])) > max:
                max = (fold_change-expected_values[transcript]) * (fold_change-expected_values[transcript])
                max_transcript = transcript
        else:
            if (transcript in values_1):
                fold_change = values_1[transcript] + 1
                errors.append((fold_change-expected_values[transcript]) * (fold_change-expected_values[transcript]))
            elif (transcript in values_2):
                fold_change = values_2[transcript]+ 1
                errors.append((fold_change-expected_values[transcript]) * (fold_change-expected_values[transcript]))
            else:
                # If transcript is not found in both experiments, add expected value as error
                # This is necessary because salmon does not give an answer for some transcripts for some reason
                errors.append((expected_values[transcript]) * (expected_values[transcript]))
    mean_square_error = np.mean(errors)
    mse.append(mean_square_error)
    #print(max, max_transcript, count)


print("Mean Squared error:\n", mse)
print(np.mean(mse), np.var(mse))
