#!/bin/bash

# Index building
needle="Path to Needle executable"

/usr/bin/time -v $needle ibfmin -l 15 -f 0.05 w_41/SRR*minimiser -o w_41/SRR_ &> w_41/SRR_Time.txt
/usr/bin/time -v $needle ibfmin -l 15 -f 0.05 w_41/SRR*minimiser -o w_41/SRR_Compressed -c &> w_41/SRR_Compressed_Time.txt
/usr/bin/time -v $needle ibfmin -l 15 -f 0.3 w_41/SRR*minimiser -o w_41/SRR_FPR03_ &> w_41/SRR_FPR03_Time.txt
/usr/bin/time -v $needle ibfmin -l 15 -f 0.3 w_41/SRR*minimiser -o w_41/SRR_FPR03_Compressed -c &> w_41/SRR_FPR03_Compressed_Time.txt

/usr/bin/time -v $needle ibfmin -l 15 -f 0.05 w_25/SRR*minimiser -o w_25/SRR_ &> w_25/SRR_Time.txt
/usr/bin/time -v $needle ibfmin -l 15 -f 0.05 w_25/SRR*minimiser -o w_25/SRR_Compressed -c &> w_25/SRR_Compressed_Time.txt
/usr/bin/time -v $needle ibfmin -l 15 -f 0.3 w_25/SRR*minimiser -o w_25/SRR_FPR03_ &> w_25/SRR_FPR03_Time.txt
/usr/bin/time -v $needle ibfmin -l 15 -f 0.3 w_25/SRR*minimiser -o w_25/SRR_FPR03_Compressed -c  &> w_25/SRR_FPR03_Compressed_Time.txt

/usr/bin/time -v $needle ibfmin -l 15 -f 0.05 w_21/SRR*minimiser -o w_21/SRR_ &> w_21/SRR_Time.txt
/usr/bin/time -v $needle ibfmin -l 15 -f 0.05 w_21/SRR*minimiser -o w_21/SRR_Compressed -c &> w_21/SRR_Compressed_Time.txt
/usr/bin/time -v $needle ibfmin -l 15 -f 0.05 w_21/SRR*minimiser -o w_21/SRR_Thread4_ -t 4 &> w_21/SRR_Thread4_Time.txt
/usr/bin/time -v $needle ibfmin -l 15 -f 0.05 w_21/SRR*minimiser -o w_21/SRR_Thread4_Compressed -t 4 &> w_21/SRR_Thread4_Compressed_Time.txt
/usr/bin/time -v $needle ibfmin -l 15 -f 0.3 w_21/SRR*minimiser -o w_21/SRR_FPR03_ &> w_21/SRR_FPR03_Time.txt
/usr/bin/time -v $needle ibfmin -l 15 -f 0.3 w_21/SRR*minimiser -o w_21/SRR_FPR03_Compressed -c &> w_21/SRR_FPR03_Compressed_Time.txt

reindeer="Path to reindeer executable"

/usr/bin/time -v $reindeer --index -f reindeer/fof_large_data.lst -k 21 -o reindeer/out_large_thread4 -t 4 &> reindeer/Time.txt
/usr/bin/time -v $reindeer --index -f reindeer/fof_large_data.lst -k 21 -o reindeer/out_large_log --log-count &> reindeer/Time_Log.txt


# Querying

/usr/bin/time -v $needle estimate -i w_41/SRR_ data/query_1.fa -o  w_41/expressions_SRR_1.out  &> w_41/Time_Query_1.txt
/usr/bin/time -v $needle estimate -i w_41/SRR_ data/query_100.fa -o  w_41/expressions_SRR_100.out &> w_41/Time_Query_100.txt
/usr/bin/time -v $needle estimate -i w_41/SRR_ data/query_1000.fa -o  w_41/expressions_SRR_1000.out &> w_41/Time_Query_1000.txt
/usr/bin/time -v $needle estimate -i w_25/SRR_ data/query_1.fa -o  w_25/expressions_SRR_1.out &> w_25/Time_Query_1.txt
/usr/bin/time -v $needle estimate -i w_25/SRR_ data/query_100.fa -o  w_25/expressions_SRR_100.out &> w_25/Time_Query_100.txt
/usr/bin/time -v $needle estimate -i w_25/SRR_ data/query_1000.fa -o  w_25/expressions_SRR_1000.out &> w_25/Time_Query_1000.txt
/usr/bin/time -v $needle estimate -i w_21/SRR_ data/query_1.fa -o  w_21/expressions_SRR_1.out &> w_21/Time_Query_1.txt
/usr/bin/time -v $needle estimate -i w_21/SRR_ data/query_100.fa -o  w_21/expressions_SRR_100.out &> w_21/Time_Query_100.txt
/usr/bin/time -v $needle estimate -i w_21/SRR_ data/query_1000.fa -o  w_21/expressions_SRR_1000.out &> w_21/Time_Query_1000.txt


/usr/bin/time -v $needle estimate -i w_41/SRR_Compressed data/query_1.fa -o  w_41/expressions_SRR_1.out  &> w_41/Time_Compressed_Query_1.txt
/usr/bin/time -v $needle estimate -i w_41/SRR_Compressed data/query_100.fa -o  w_41/expressions_SRR_100.out &> w_41/Time_Compressed_Query_100.txt
/usr/bin/time -v $needle estimate -i w_41/SRR_Compressed data/query_1000.fa -o  w_41/expressions_SRR_1000.out &> w_41/Time_Compressed_Query_1000.txt
/usr/bin/time -v $needle estimate -i w_25/SRR_Compressed data/query_1.fa -o  w_25/expressions_SRR_1.out &> w_25/Time_Compressed_Query_1.txt
/usr/bin/time -v $needle estimate -i w_25/SRR_Compressed data/query_100.fa -o  w_25/expressions_SRR_100.out &> w_25/Time_Compressed_Query_100.txt
/usr/bin/time -v $needle estimate -i w_25/SRR_Compressed data/query_1000.fa -o  w_25/expressions_SRR_1000.out &> w_25/Time_Compressed_Query_1000.txt
/usr/bin/time -v $needle estimate -i w_21/SRR_Compressed data/query_1.fa -o  w_21/expressions_SRR_1.out &> w_21/Time_Compressed_Query_1.txt
/usr/bin/time -v $needle estimate -i w_21/SRR_Compressed data/query_100.fa -o  w_21/expressions_SRR_100.out &> w_21/Time_Compressed_Query_100.txt
/usr/bin/time -v $needle estimate -i w_21/SRR_Compressed data/query_1000.fa -o  w_21/expressions_SRR_1000.out &> w_21/Time_Compressed_Query_1000.txt

/usr/bin/time -v $needle estimate -i w_41/SRR_FPR03_ data/query_1.fa -o  w_41/expressions_SRR_FPR03_1.out  &> w_41/Time_FPR03_Query_1.txt
/usr/bin/time -v $needle estimate -i w_41/SRR_FPR03_ data/query_100.fa -o  w_41/expressions_SRR_FPR03_100.out &> w_41/Time_FPR03_Query_100.txt
/usr/bin/time -v $needle estimate -i w_41/SRR_FPR03_ data/query_1000.fa -o  w_41/expressions_SRR_FPR03_1000.out &> w_41/Time_FPR03_Query_1000.txt
/usr/bin/time -v $needle estimate -i w_25/SRR_FPR03_ data/query_1.fa -o  w_25/expressions_SRR_FPR03_1.out &> w_25/Time_FPR03_Query_1.txt
/usr/bin/time -v $needle estimate -i w_25/SRR_FPR03_ data/query_100.fa -o  w_25/expressions_SRR_FPR03_100.out &> w_25/Time_FPR03_Query_100.txt
/usr/bin/time -v $needle estimate -i w_25/SRR_FPR03_ data/query_1000.fa -o  w_25/expressions_SRR_FPR03_1000.out &> w_25/Time_FPR03_Query_1000.txt
/usr/bin/time -v $needle estimate -i w_21/SRR_FPR03_ data/query_1.fa -o  w_21/expressions_SRR_FPR03_1.out &> w_21/Time_FPR03_Query_1.txt
/usr/bin/time -v $needle estimate -i w_21/SRR_FPR03_ data/query_100.fa -o  w_21/expressions_SRR_FPR03_100.out &> w_21/Time_FPR03_Query_100.txt
/usr/bin/time -v $needle estimate -i w_21/SRR_FPR03_ data/query_1000.fa -o  w_21/expressions_SRR_FPR03_1000.out &> w_21/Time_FPR03_Query_1000.txt


/usr/bin/time -v $needle estimate -i w_41/SRR_FPR03_Compressed data/query_1.fa -o  w_41/expressions_SRR_FPR03_Compressed_1.out  &> w_41/Time_FPR03_Compressed_Query_1.txt
/usr/bin/time -v $needle estimate -i w_41/SRR_FPR03_Compressed data/query_100.fa -o  w_41/expressions_SRR_FPR03_Compressed_100.out &> w_41/Time_FPR03_Compressed_Query_100.txt
/usr/bin/time -v $needle estimate -i w_41/SRR_FPR03_Compressed data/query_1000.fa -o  w_41/expressions_SRR_FPR03_Compressed_1000.out &> w_41/Time_FPR03_Compressed_Query_1000.txt
/usr/bin/time -v $needle estimate -i w_25/SRR_FPR03_Compressed data/query_1.fa -o  w_25/expressions_SRR_FPR03_Compressed_1.out &> w_25/Time_FPR03_Compressed_Query_1.txt
/usr/bin/time -v $needle estimate -i w_25/SRR_FPR03_Compressed data/query_100.fa -o  w_25/expressions_SRR_FPR03_Compressed_100.out &> w_25/Time_FPR03_Compressed_Query_100.txt
/usr/bin/time -v $needle estimate -i w_25/SRR_FPR03_Compressed data/query_1000.fa -o  w_25/expressions_SRR_FPR03_Compressed_1000.out &> w_25/Time_FPR03_Compressed_Query_1000.txt
/usr/bin/time -v $needle estimate -i w_21/SRR_FPR03_Compressed data/query_1.fa -o  w_21/expressions_SRR_FPR03_Compressed_1.out &> w_21/Time_FPR03_Compressed_Query_1.txt
/usr/bin/time -v $needle estimate -i w_21/SRR_FPR03_Compressed data/query_100.fa -o  w_21/expressions_SRR_FPR03_Compressed_100.out &> w_21/Time_FPR03_Compressed_Query_100.txt
/usr/bin/time -v $needle estimate -i w_21/SRR_FPR03_Compressed data/query_1000.fa -o  w_21/expressions_SRR_FPR03_Compressed_1000.out &> w_21/Time_FPR03_Compressed_Query_1000.txt


/usr/bin/time -v ./Reindeer --query -l reindeer/out_large_thread4 -o reindeer/Query_1_ -q data/query_1.fa &> reindeer/Time_Query_1.txt
/usr/bin/time -v ./Reindeer --query -l reindeer/out_large_thread4 -o reindeer/Query_100_ -q data/query_100.fa &> reindeer/Time_Query_100.txt
/usr/bin/time -v ./Reindeer --query -l reindeer/out_large_thread4 -o reindeer/Query_1000_ -q data/query_1000.fa &> reindeer/Time_Query_1000.txt
