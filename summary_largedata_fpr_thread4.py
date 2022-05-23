# Can be called after run_large_dataset.sh was run.

import os

def read_time_file(input_file):
    with open(input_file, 'r') as f:
        for line in f:
            if "Elapsed (wall clock)" in line:
                time = line.split()[-1]
            elif ("Maximum resident set size" in line):
                ram_usage = round(int(line.split()[-1])/1000000, 1)
    return (time, ram_usage)

def get_index_size_Needle(input_dir):
    size = 0
    size += os.path.getsize(input_dir + "IBF_Data")
    size += os.path.getsize(input_dir + "IBF_Levels.levels")
    size += os.path.getsize(input_dir + "IBF_FPRs.fprs")
    for i in range(15):
        size += os.path.getsize(input_dir + "IBF_Level_"+str(i))
    size = round(size/1000000000, 1)
    return size

def get_index_size_Reindeer(input_dir):
    size = 0
    size += os.path.getsize(input_dir + "reindeer_index.gz")
    size += os.path.getsize(input_dir + "reindeer_matrix_eqc.gz")
    size += os.path.getsize(input_dir + "reindeer_matrix_eqc_info")
    size += os.path.getsize(input_dir + "reindeer_matrix_eqc_position.gz")
    size = round(size/1000000000, 1)
    return size

times = []
ram_usages = []
sizes = []
times_query = []
ram_usages_query = []

# Get times and maximal RAM usage
time, ram_usage = read_time_file("reindeer/Time.txt")
times.append(time)
ram_usages.append(ram_usage)
time, ram_usage = read_time_file("reindeer/Time_Thread4_Log.txt")
times.append(time)
ram_usages.append(ram_usage)
time, ram_usage = read_time_file("w_21/SRR_Thread4_FPR03_Time.txt")
times.append(time)
ram_usages.append(ram_usage)
time, ram_usage = read_time_file("w_21/SRR_Thread4_FPR03_Compressed_Time.txt")
times.append(time)
ram_usages.append(ram_usage)
time, ram_usage = read_time_file("w_25/SRR_Thread4_FPR03_Time.txt")
times.append(time)
ram_usages.append(ram_usage)
time, ram_usage = read_time_file("w_25/SRR_Thread4_FPR03_Compressed_Time.txt")
times.append(time)
ram_usages.append(ram_usage)
time, ram_usage = read_time_file("w_41/SRR_Thread4_FPR03_Time.txt")
times.append(time)
ram_usages.append(ram_usage)
time, ram_usage = read_time_file("w_41/SRR_Thread4_FPR03_Compressed_Time.txt")
times.append(time)
ram_usages.append(ram_usage)


# Get sizes
sizes.append(get_index_size_Reindeer("reindeer/out_large_thread4/"))
sizes.append(get_index_size_Reindeer("reindeer/out_large_log/"))
sizes.append(get_index_size_Needle("w_21/SRR_Thread4_FPR03_"))
sizes.append(get_index_size_Needle("w_21/SRR_Thread4_FPR03_Compressed"))
sizes.append(get_index_size_Needle("w_25/SRR_Thread4_FPR03_"))
sizes.append(get_index_size_Needle("w_25/SRR_Thread4_FPR03_Compressed"))
sizes.append(get_index_size_Needle("w_41/SRR_Thread4_FPR03_"))
sizes.append(get_index_size_Needle("w_41/SRR_Thread4_FPR03_Compressed"))

# Query information

# 1
time, ram_usage = read_time_file("reindeer/Time_Thread4_Query_1.txt")
times_query.append(time)
ram_usages_query.append(ram_usage)
time, ram_usage = read_time_file("w_21/Time_Thread4_FPR03_Query_1.txt")
times_query.append(time)
ram_usages_query.append(ram_usage)
time, ram_usage = read_time_file("w_21/Time_Thread4_Compressed_Query_1.txt")
times_query.append(time)
ram_usages_query.append(ram_usage)
time, ram_usage = read_time_file("w_25/Time_Thread4_FPR03_Query_1.txt")
times_query.append(time)
ram_usages_query.append(ram_usage)
time, ram_usage = read_time_file("w_25/Time_Thread4_Compressed_Query_1.txt")
times_query.append(time)
ram_usages_query.append(ram_usage)
time, ram_usage = read_time_file("w_41/Time_Thread4_FPR03_Query_1.txt")
times_query.append(time)
ram_usages_query.append(ram_usage)
time, ram_usage = read_time_file("w_41/Time_Thread4_Compressed_Query_1.txt")
times_query.append(time)
ram_usages_query.append(ram_usage)


# 100
time, ram_usage = read_time_file("reindeer/Time_Thread4_Query_100.txt")
times_query.append(time)
ram_usages_query.append(ram_usage)
time, ram_usage = read_time_file("w_21/Time_Thread4_FPR03_Query_100.txt")
times_query.append(time)
ram_usages_query.append(ram_usage)
time, ram_usage = read_time_file("w_21/Time_Thread4_Compressed_Query_100.txt")
times_query.append(time)
ram_usages_query.append(ram_usage)
time, ram_usage = read_time_file("w_25/Time_Thread4_FPR03_Query_100.txt")
times_query.append(time)
ram_usages_query.append(ram_usage)
time, ram_usage = read_time_file("w_25/Time_Thread4_Compressed_Query_100.txt")
times_query.append(time)
ram_usages_query.append(ram_usage)
time, ram_usage = read_time_file("w_41/Time_Thread4_FPR03_Query_100.txt")
times_query.append(time)
ram_usages_query.append(ram_usage)
time, ram_usage = read_time_file("w_41/Time_Thread4_Compressed_Query_100.txt")
times_query.append(time)
ram_usages_query.append(ram_usage)

# 1000
time, ram_usage = read_time_file("reindeer/Time_Thread4_Query_1000.txt")
times_query.append(time)
ram_usages_query.append(ram_usage)
time, ram_usage = read_time_file("w_21/Time_Thread4_FPR03_Query_1000.txt")
times_query.append(time)
ram_usages_query.append(ram_usage)
time, ram_usage = read_time_file("w_21/Time_Thread4_Compressed_Query_1000.txt")
times_query.append(time)
ram_usages_query.append(ram_usage)
time, ram_usage = read_time_file("w_25/Time_Thread4_FPR03_Query_1000.txt")
times_query.append(time)
ram_usages_query.append(ram_usage)
time, ram_usage = read_time_file("w_25/Time_Thread4_Compressed_Query_1000.txt")
times_query.append(time)
ram_usages_query.append(ram_usage)
time, ram_usage = read_time_file("w_41/Time_Thread4_FPR03_Query_1000.txt")
times_query.append(time)
ram_usages_query.append(ram_usage)
time, ram_usage = read_time_file("w_41/Time_Thread4_Compressed_Query_1000.txt")
times_query.append(time)
ram_usages_query.append(ram_usage)

print(times)
print(ram_usages)
print(sizes)
print(times_query)
print(ram_usages_query)
methods = ["REINDEER", "REINDEER log", "Meedle", "Needle Compressed","Needle (w=25)", "Needle (w=25) Compressed", "Needle (w=41)", "Needle (w=41) Compressed"]
with open("Results_Large_Data_set_FPR03.txt", 'w') as o:
    o.write("Build    \t" + "TIME" +"\t" + "RAM" + "\t" + "Index Size" + "\n")
    for i in range(len(times)):
        o.write(methods[i] + "\t" + str(times[i]) +"\t" + str(ram_usages[i]) + "\t" + str(sizes[i]) + "\n")

    o.write("\n")
    o.write("Query\tREINDEER\tNeedle\tNeedle Compressed\tNeedle (w=25)\tNeedle (w=25)Compressed\tNeedle (w=41)\tNeedle (w=41) Compressed\n")
    j = 0
    o.write("1 Time\t" + times_query[j] + "\t" + times_query[j+1] + "\t" + times_query[j+2] +"\t" + times_query[j+3]  +"\t" + times_query[j+4] +"\t" + times_query[j+5] +"\t" + times_query[j+6]+"\n")
    o.write("1 RAM\t" + str(ram_usages_query[j]) + "\t" + str(ram_usages_query[j+1]) + "\t" + str(ram_usages_query[j+2]) +"\t" + str(ram_usages_query[j+3])  + "\t" +str(ram_usages_query[j+4])  + "\t"+ str(ram_usages_query[j+5]) + "\t"  + str(ram_usages_query[j+6])  +"\n")
    j = 7
    o.write("100 Time\t" + times_query[j] + "\t" + times_query[j+1] + "\t" + times_query[j+2] +"\t" + times_query[j+3]  + times_query[j+4] +"\t" + times_query[j+5] +"\t" + times_query[j+6]+"\n")
    o.write("100 RAM\t" + str(ram_usages_query[j]) + "\t" + str(ram_usages_query[j+1]) + "\t" + str(ram_usages_query[j+2]) +"\t" + str(ram_usages_query[j+3])  + "\t" +str(ram_usages_query[j+4])  + "\t"+ str(ram_usages_query[j+5]) + "\t"  + str(ram_usages_query[j+6])+"\n")
    j = 14
    o.write("1000 Time\t" + times_query[j] + "\t" + times_query[j+1] + "\t" + times_query[j+2] +"\t" + times_query[j+3]  + times_query[j+4] +"\t" + times_query[j+5] +"\t" + times_query[j+6]+"\n")
    o.write("1000 RAM\t" + str(ram_usages_query[j]) + "\t" + str(ram_usages_query[j+1]) + "\t" + str(ram_usages_query[j+2]) +"\t" + str(ram_usages_query[j+3]) + "\t" +str(ram_usages_query[j+4])  + "\t"+ str(ram_usages_query[j+5]) + "\t"  + str(ram_usages_query[j+6]) +"\n")
