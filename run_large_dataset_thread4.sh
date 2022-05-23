#!/bin/bash

# Index building 4 Threads
needle="Path to Needle executable"

/usr/bin/time -v $needle ibfmin -l 15 -f 0.05 w_41/SRR*minimiser -o w_41/SRR_Thread4_ -t 4 &> w_41/SRR_Thread4_Time.txt
/usr/bin/time -v $needle ibfmin -l 15 -f 0.05 w_41/SRR*minimiser -o w_41/SRR_Thread4_Compressed -t 4 &> w_41/SRR_Thread4_Compressed_Time.txt
/usr/bin/time -v $needle ibfmin -l 15 -f 0.3 w_41/SRR*minimiser -o w_41/SRR_Thread4_FPR03_ -t 4&> w_41/SRR_Thread4_FPR03_Time.txt
/usr/bin/time -v $needle ibfmin -l 15 -f 0.3 w_41/SRR*minimiser -o w_41/SRR_Thread4_FPR03_Compressed -c -t 4 &> w_41/SRR_Thread4_FPR03_Compressed_Time.txt

/usr/bin/time -v $needle ibfmin -l 15 -f 0.05 w_25/SRR*minimiser -o w_25/SRR_Thread4_ -t 4 &> w_25/SRR_Thread4_Time.txt
/usr/bin/time -v $needle ibfmin -l 15 -f 0.05 w_25/SRR*minimiser -o w_25/SRR_Thread4_Compressed -t 4 &> w_25/SRR_Thread4_Compressed_Time.txt
/usr/bin/time -v $needle ibfmin -l 15 -f 0.3 w_25/SRR*minimiser -o w_25/SRR_Thread4_FPR03_ -t 4 &> w_25/SRR_Thread4_FPR03_Time.txt
/usr/bin/time -v $needle ibfmin -l 15 -f 0.3 w_25/SRR*minimiser -o w_25/SRR_Thread4_FPR03_Compressed -c -t 4 &> w_25/SRR_Thread4_FPR03_Compressed_Time.txt

/usr/bin/time -v $needle ibfmin -l 15 -f 0.05 w_21/SRR*minimiser -o w_21/SRR_Thread4_ -t 4 &> w_21/SRR_Thread4_Time.txt
/usr/bin/time -v $needle ibfmin -l 15 -f 0.05 w_21/SRR*minimiser -o w_21/SRR_Thread4_Compressed -t 4 &> w_21/SRR_Thread4_Compressed_Time.txt
/usr/bin/time -v $needle ibfmin -l 15 -f 0.3 w_21/SRR*minimiser -o w_21/SRR_Thread4_FPR03_ -t 4 &> w_21/SRR_Thread4_FPR03_Time.txt
/usr/bin/time -v $needle ibfmin -l 15 -f 0.3 w_21/SRR*minimiser -o w_21/SRR_Thread4_FPR03_Compressed -c -t 4 &> w_21/SRR_Thread4_FPR03_Compressed_Time.txt

reindeer="Path to reindeer executable"

/usr/bin/time -v $reindeer --index -f reindeer/fof_large_data.lst -k 21 -o reindeer/out_large_thread4 -t 4 &> reindeer/Time.txt
/usr/bin/time -v $reindeer --index -f reindeer/fof_large_data.lst -k 21 -o reindeer/out_large_log --log-count -t 4 &> reindeer/Time_Thread4_Log.txt

# Querying 4 Threads

/usr/bin/time -v $needle estimate -i w_41/SRR_Thread4_ data/query_1.fa -o  w_41/expressions_SRR_Thread4_1.out  -t 4  &> w_41/Time_Thread4_Query_1.txt
/usr/bin/time -v $needle estimate -i w_41/SRR_Thread4_ data/query_100.fa -o  w_41/expressions_SRR_Thread4_100.out -t 4  &> w_41/Time_Thread4_Query_100.txt
/usr/bin/time -v $needle estimate -i w_41/SRR_Thread4_ data/query_1000.fa -o  w_41/expressions_SRR_Thread4_1000.out -t 4  &> w_41/Time_Thread4_Query_1000.txt
/usr/bin/time -v $needle estimate -i w_25/SRR_Thread4_ data/query_1.fa -o  w_25/expressions_SRR_Thread4_1.out -t 4  &> w_25/Time_Thread4_Query_1.txt
/usr/bin/time -v $needle estimate -i w_25/SRR_Thread4_ data/query_100.fa -o  w_25/expressions_SRR_Thread4_100.out -t 4  &> w_25/Time_Thread4_Query_100.txt
/usr/bin/time -v $needle estimate -i w_25/SRR_Thread4_ data/query_1000.fa -o  w_25/expressions_SRR_Thread4_1000.out -t 4  &> w_25/Time_Thread4_Query_1000.txt
/usr/bin/time -v $needle estimate -i w_21/SRR_Thread4_ data/query_1.fa -o  w_21/expressions_SRR_Thread4_1.out -t 4 &> w_21/Time_Thread4_Query_1.txt
/usr/bin/time -v $needle estimate -i w_21/SRR_Thread4_ data/query_100.fa -o  w_21/expressions_SRR_Thread4_100.out -t 4 &> w_21/Time_Thread4_Query_100.txt
/usr/bin/time -v $needle estimate -i w_21/SRR_Thread4_ data/query_1000.fa -o  w_21/expressions_SRR_Thread4_1000.out -t 4 &> w_21/Time_Thread4_Query_1000.txt

/usr/bin/time -v $needle estimate -i w_41/SRR_Thread4_Compressed data/query_1.fa -o  w_41/expressions_SRR_Thread4_1.out -t 4  &> w_41/Time_Compressed_Query_1.txt
/usr/bin/time -v $needle estimate -i w_41/SRR_Thread4_Compressed data/query_100.fa -o  w_41/expressions_SRR_Thread4_100.out -t 4 &> w_41/Time_Compressed_Query_100.txt
/usr/bin/time -v $needle estimate -i w_41/SRR_Thread4_Compressed data/query_1000.fa -o  w_41/expressions_SRR_Thread4_1000.out -t 4  &> w_41/Time_Compressed_Query_1000.txt
/usr/bin/time -v $needle estimate -i w_25/SRR_Thread4_Compressed data/query_1.fa -o  w_25/expressions_SRR_Thread4_1.out -t 4  &> w_25/Time_Compressed_Query_1.txt
/usr/bin/time -v $needle estimate -i w_25/SRR_Thread4_Compressed data/query_100.fa -o  w_25/expressions_SRR_Thread4_100.out -t 4 &> w_25/Time_Compressed_Query_100.txt
/usr/bin/time -v $needle estimate -i w_25/SRR_Thread4_Compressed data/query_1000.fa -o  w_25/expressions_SRR_Thread4_1000.out -t 4  &> w_25/Time_Compressed_Query_1000.txt
/usr/bin/time -v $needle estimate -i w_21/SRR_Thread4_Compressed data/query_1.fa -o  w_21/expressions_SRR_Thread4_1.out -t 4  &> w_21/Time_Compressed_Query_1.txt
/usr/bin/time -v $needle estimate -i w_21/SRR_Thread4_Compressed data/query_100.fa -o  w_21/expressions_SRR_Thread4_100.out -t 4 &> w_21/Time_Compressed_Query_100.txt
/usr/bin/time -v $needle estimate -i w_21/SRR_Thread4_Compressed data/query_1000.fa -o  w_21/expressions_SRR_Thread4_1000.out -t 4  &> w_21/Time_Compressed_Query_1000.txt

/usr/bin/time -v $needle estimate -i w_41/SRR_Thread4_FPR03_ data/query_1.fa -o  w_41/expressions_SRR_Thread4_FPR03_1.out -t 4 &> w_41/Time_Thread4_FPR03_Query_1.txt
/usr/bin/time -v $needle estimate -i w_41/SRR_Thread4_FPR03_ data/query_100.fa -o  w_41/expressions_SRR_Thread4_FPR03_100.out -t 4  &> w_41/Time_Thread4_FPR03_Query_100.txt
/usr/bin/time -v $needle estimate -i w_41/SRR_Thread4_FPR03_ data/query_1000.fa -o  w_41/expressions_SRR_Thread4_FPR03_1000.out -t 4  &> w_41/Time_Thread4_FPR03_Query_1000.txt
/usr/bin/time -v $needle estimate -i w_25/SRR_Thread4_FPR03_ data/query_1.fa -o  w_25/expressions_SRR_Thread4_FPR03_1.out -t 4 &> w_25/Time_Thread4_FPR03_Query_1.txt
/usr/bin/time -v $needle estimate -i w_25/SRR_Thread4_FPR03_ data/query_100.fa -o  w_25/expressions_SRR_Thread4_FPR03_100.out -t 4  &> w_25/Time_Thread4_FPR03_Query_100.txt
/usr/bin/time -v $needle estimate -i w_25/SRR_Thread4_FPR03_ data/query_1000.fa -o  w_25/expressions_SRR_Thread4_FPR03_1000.out -t 4  &> w_25/Time_Thread4_FPR03_Query_1000.txt
/usr/bin/time -v $needle estimate -i w_21/SRR_Thread4_FPR03_ data/query_1.fa -o  w_21/expressions_SRR_Thread4_FPR03_1.out -t 4 &> w_21/Time_Thread4_FPR03_Query_1.txt
/usr/bin/time -v $needle estimate -i w_21/SRR_Thread4_FPR03_ data/query_100.fa -o  w_21/expressions_SRR_Thread4_FPR03_100.out -t 4 &> w_21/Time_Thread4_FPR03_Query_100.txt
/usr/bin/time -v $needle estimate -i w_21/SRR_Thread4_FPR03_ data/query_1000.fa -o  w_21/expressions_SRR_Thread4_FPR03_1000.out -t 4 &> w_21/Time_Thread4_FPR03_Query_1000.txt

/usr/bin/time -v $needle estimate -i w_41/SRR_Thread4_FPR03_Compressed data/query_1.fa -o  w_41/expressions_SRR_Thread4_FPR03_Compressed_1.out -t 4  &> w_41/Time_Thread4_FPR03_Compressed_Query_1.txt
/usr/bin/time -v $needle estimate -i w_41/SRR_Thread4_FPR03_Compressed data/query_100.fa -o  w_41/expressions_SRR_Thread4_FPR03_Compressed_100.out -t 4 &> w_41/Time_Thread4_FPR03_Compressed_Query_100.txt
/usr/bin/time -v $needle estimate -i w_41/SRR_Thread4_FPR03_Compressed data/query_1000.fa -o  w_41/expressions_SRR_Thread4_FPR03_Compressed_1000.out -t 4  &> w_41/Time_Thread4_FPR03_Compressed_Query_1000.txt
/usr/bin/time -v $needle estimate -i w_25/SRR_Thread4_FPR03_Compressed data/query_1.fa -o  w_25/expressions_SRR_Thread4_FPR03_Compressed_1.out -t 4 &> w_25/Time_Thread4_FPR03_Compressed_Query_1.txt
/usr/bin/time -v $needle estimate -i w_25/SRR_Thread4_FPR03_Compressed data/query_100.fa -o  w_25/expressions_SRR_Thread4_FPR03_Compressed_100.out -t 4 &> w_25/Time_Thread4_FPR03_Compressed_Query_100.txt
/usr/bin/time -v $needle estimate -i w_25/SRR_Thread4_FPR03_Compressed data/query_1000.fa -o  w_25/expressions_SRR_Thread4_FPR03_Compressed_1000.out -t 4  &> w_25/Time_Thread4_FPR03_Compressed_Query_1000.txt
/usr/bin/time -v $needle estimate -i w_21/SRR_Thread4_FPR03_Compressed data/query_1.fa -o  w_21/expressions_SRR_Thread4_FPR03_Compressed_1.out -t 4  &> w_21/Time_Thread4_FPR03_Compressed_Query_1.txt
/usr/bin/time -v $needle estimate -i w_21/SRR_Thread4_FPR03_Compressed data/query_100.fa -o  w_21/expressions_SRR_Thread4_FPR03_Compressed_100.out -t 4 &> w_21/Time_Thread4_FPR03_Compressed_Query_100.txt
/usr/bin/time -v $needle estimate -i w_21/SRR_Thread4_FPR03_Compressed data/query_1000.fa -o  w_21/expressions_SRR_Thread4_FPR03_Compressed_1000.out -t 4 &> w_21/Time_Thread4_FPR03_Compressed_Query_1000.txt


/usr/bin/time -v $reindeer --query -l reindeer/out_large_thread4 -o reindeer/Thread4_Query_1_ -q data/query_1.fa -t 4 &> reindeer/Time_Thread4_Query_1.txt
/usr/bin/time -v $reindeer --query -l reindeer/out_large_thread4 -o reindeer/Thread4_Query_100_ -q data/query_100.fa -t 4 &> reindeer/Time_Thread4_Query_100.txt
/usr/bin/time -v $reindeer --query -l reindeer/out_large_thread4 -o reindeer/Thread4_Query_1000_ -q data/query_1000.fa -t 4 &> reindeer/Time_Thread4_Query_1000.txt


# Get summary in Results_* files for FPR 0.05 and FPR 0.3
python3 summary_largedata_thread4.py
python3 summary_largedata_fpr_thread4.py
