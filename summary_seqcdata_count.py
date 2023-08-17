import check_titration_monotonicity as ctm
import compare_seqc_microarray as csm
import sys

input_dir = sys.argv[1]
out_file = "Results_SEQC_Data_set_count.txt"

kallisto_files = []
salmon_files = []
Needle_file = []
Needle_23_file = []
Needle_39_file = []

for letter in "ABCD":
    for i in range(1, 5):
        kallisto_files.append("kallisto/SEQC2012-ILM-AGR-" + letter +"-"+ str(i) +"/abundance.tsv")
        salmon_files.append("salmon/SEQC2012-ILM-AGR-" + letter +"-" + str(i) +".out"+ "/quant.sf")
        Needle_file.append("w_19/Count_SEQC2012-ILM-AGR-" + letter +"-"+str(i) +"_R2.fastq.count.out")
        Needle_23_file.append("w_23/Count_SEQC2012-ILM-AGR-" + letter +"-"+str(i) +"_R2.fastq.count.out")
        Needle_39_file.append("w_39/Count_SEQC2012-ILM-AGR-" + letter +"-"+str(i) +"_R2.fastq.count.out")

kallisto_seqc = csm.get_exp_mean(2, 0, input_dir+"/Taqman-raw.txt", "", kallisto_files, 16)
kallisto_microrray = csm.get_exp_mean(2, 1, input_dir+"/BeadChip-Log2.table", "data/entrez_id2gene_id.txt", kallisto_files, 16)
kallisto_mse = ctm.get_mse(2, 16, kallisto_files)

salmon_seqc = csm.get_exp_mean(1, 0, input_dir+"/Taqman-raw.txt", "", salmon_files, 16)
salmon_microrray = csm.get_exp_mean(1, 1, input_dir+"/BeadChip-Log2.table", "data/entrez_id2gene_id.txt", salmon_files, 16)
salmon_mse = ctm.get_mse(1, 16, salmon_files)

needle_seqc = csm.get_exp_mean(0, 0, input_dir+"/Taqman-raw.txt", "", Needle_file, 16)
needle_microrray = csm.get_exp_mean(0, 1, input_dir+"/BeadChip-Log2.table", "data/entrez_id2gene_id.txt", Needle_file, 16)
needle_mse = ctm.get_mse(0, 16, Needle_file)

needle23_seqc = csm.get_exp_mean(0, 0, input_dir+"/Taqman-raw.txt", "", Needle_23_file, 16)
needle23_microrray = csm.get_exp_mean(0, 1, input_dir+"/BeadChip-Log2.table", "data/entrez_id2gene_id.txt", Needle_23_file, 16)
needle23_mse = ctm.get_mse(0, 16, Needle_23_file)

needle39_seqc = csm.get_exp_mean(0, 0, input_dir+"/Taqman-raw.txt", "", Needle_39_file, 16)
needle39_microrray = csm.get_exp_mean(0, 1, input_dir+"/BeadChip-Log2.table", "data/entrez_id2gene_id.txt", Needle_39_file, 16)
needle39_mse = ctm.get_mse(0, 16, Needle_39_file)

with open(out_file, 'w') as o:
    o.write("\tSEQC\tMicroarray\tMSE\n")
    o.write("kallisto\t"+ str(kallisto_seqc[0]) + "\t" + str(kallisto_microrray[0]) + "\t" + str(kallisto_mse[1]) + "\n")
    o.write("salmon\t"+ str(salmon_seqc[0]) + "\t" + str(salmon_microrray[0]) + "\t" + str(salmon_mse[1]) + "\n")
    o.write("Needle\t"+ str(needle_seqc[0]) + "\t" + str(needle_microrray[0]) + "\t" + str(needle_mse[1]) + "\n")
    o.write("Needle (w=23)\t"+ str(needle23_seqc[0]) + "\t" + str(needle23_microrray[0]) + "\t" + str(needle23_mse[1]) + "\n")
    o.write("Needle (w=39)\t"+ str(needle39_seqc[0]) + "\t" + str(needle39_microrray[0]) + "\t" + str(needle39_mse[1]) + "\n")
