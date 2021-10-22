import check_titration_monotonicity as ctm
import compare_seqc_microarray as csm
import sys

input_dir = sys.argv[1]
out_file = "Results_SEQC_Data_set.txt"
out_file_03 = "Results_SEQC_Data_set_FPR03.txt"

kallisto_files = []
salmon_files = []
reindeer_file = ["reindeer/expressions_seqc.out"]
Needle_file = ["w_19/expressions_SEQC.out"]
Needle_23_file = ["w_23/expressions_SEQC.out"]
Needle_39_file = ["w_39/expressions_SEQC.out"]
Needle_file_norm = ["w_19/expressions_SEQC_norm.out"]
Needle_23_file_norm = ["w_23/expressions_SEQC_norm.out"]
Needle_39_file_norm = ["w_39/expressions_SEQC_norm.out"]

for letter in "ABCD":
    for i in range(1, 5):
        kallisto_files.append("kallisto/SEQC2012-ILM-AGR-" + letter +"-"+ str(i) +"/abundance.tsv")
        salmon_files.append("salmon/SEQC2012-ILM-AGR-" + letter +"-" + str(i) +".out"+ "/quant.sf")

kallisto_seqc = csm.get_exp_mean(1, 0, input_dir+"/Taqman-raw.txt", "", kallisto_files, 16)
kallisto_microrray = csm.get_exp_mean(1, 1, input_dir+"/BeadChip-Log2.table", "data/entrez_id2gene_id.txt", kallisto_files, 16)
kallisto_mse = ctm.get_mse(1, 16, kallisto_files)

salmon_seqc = csm.get_exp_mean(2, 0, input_dir+"/Taqman-raw.txt", "", salmon_files, 16)
salmon_microrray = csm.get_exp_mean(2, 1, input_dir+"/BeadChip-Log2.table", "data/entrez_id2gene_id.txt", salmon_files, 16)
salmon_mse = ctm.get_mse(2, 16, salmon_files)

reindeer_seqc = csm.get_exp_mean(3, 0, input_dir+"/Taqman-raw.txt", "", reindeer_file, 16)
reindeer_microrray = csm.get_exp_mean(3, 1, input_dir+"/BeadChip-Log2.table", "data/entrez_id2gene_id.txt", reindeer_file, 16)
reindeer_mse = ctm.get_mse(3, 16, reindeer_file)

needle_seqc = csm.get_exp_mean(3, 0, input_dir+"/Taqman-raw.txt", "", Needle_file, 16)
needle_microrray = csm.get_exp_mean(3, 1, input_dir+"/BeadChip-Log2.table", "data/entrez_id2gene_id.txt", Needle_file, 16)
needle_mse = ctm.get_mse(3, 16, Needle_file_norm)

needle23_seqc = csm.get_exp_mean(3, 0, input_dir+"/Taqman-raw.txt", "", Needle_23_file, 16)
needle23_microrray = csm.get_exp_mean(3, 1, input_dir+"/BeadChip-Log2.table", "data/entrez_id2gene_id.txt", Needle_23_file, 16)
needle23_mse = ctm.get_mse(3, 16, Needle_23_file_norm)

needle39_seqc = csm.get_exp_mean(3, 0, input_dir+"/Taqman-raw.txt", "", Needle_39_file, 16)
needle39_microrray = csm.get_exp_mean(3, 1, input_dir+"/BeadChip-Log2.table", "data/entrez_id2gene_id.txt", Needle_39_file, 16)
needle39_mse = ctm.get_mse(3, 16, Needle_39_file_norm)

with open(out_file, 'w') as o:
    o.write("\tSEQC\tMicroarray\tMSE\n")
    o.write("kallisto\t"+ str(kallisto_seqc[0]) + "\t" + str(kallisto_microrray[0]) + "\t" + str(kallisto_mse[1]) + "\n")
    o.write("salmon\t"+ str(salmon_seqc[0]) + "\t" + str(salmon_microrray[0]) + "\t" + str(salmon_mse[1]) + "\n")
    o.write("Reindeer\t"+ str(reindeer_seqc[0]) + "\t" + str(reindeer_microrray[0]) + "\t" + str(reindeer_mse[1]) + "\n")
    o.write("Needle\t"+ str(needle_seqc[0]) + "\t" + str(needle_microrray[0]) + "\t" + str(needle_mse[1]) + "\n")
    o.write("Needle (w=23)\t"+ str(needle23_seqc[0]) + "\t" + str(needle23_microrray[0]) + "\t" + str(needle23_mse[1]) + "\n")
    o.write("Needle (w=39)\t"+ str(needle39_seqc[0]) + "\t" + str(needle39_microrray[0]) + "\t" + str(needle39_mse[1]) + "\n")


Needle_file = ["w_19/expressions_SEQC_fpr03.out"]
Needle_23_file = ["w_23/expressions_SEQC_fpr03.out"]
Needle_39_file = ["w_39/expressions_SEQC_fpr03.out"]
Needle_file_norm = ["w_19/expressions_SEQC_norm_fpr03.out"]
Needle_23_file_norm = ["w_23/expressions_SEQC_norm_fpr03.out"]
Needle_39_file_norm = ["w_39/expressions_SEQC_norm_fpr03.out"]

needle_seqc = csm.get_exp_mean(3, 0, input_dir+"/Taqman-raw.txt", "", Needle_file, 16)
needle_microrray = csm.get_exp_mean(3, 1, input_dir+"/BeadChip-Log2.table", "data/entrez_id2gene_id.txt", Needle_file, 16)
needle_mse = ctm.get_mse(3, 16, Needle_file_norm)

needle23_seqc = csm.get_exp_mean(3, 0, input_dir+"/Taqman-raw.txt", "", Needle_23_file, 16)
needle23_microrray = csm.get_exp_mean(3, 1, input_dir+"/BeadChip-Log2.table", "data/entrez_id2gene_id.txt", Needle_23_file, 16)
needle23_mse = ctm.get_mse(3, 16, Needle_23_file_norm)

needle39_seqc = csm.get_exp_mean(3, 0, input_dir+"/Taqman-raw.txt", "", Needle_39_file, 16)
needle39_microrray = csm.get_exp_mean(3, 1, input_dir+"/BeadChip-Log2.table", "data/entrez_id2gene_id.txt", Needle_39_file, 16)
needle39_mse = ctm.get_mse(3, 16, Needle_39_file_norm)

with open(out_file_03, 'w') as o:
    o.write("\tSEQC\tMicroarray\tMSE\n")
    o.write("Needle\t"+ str(needle_seqc[0]) + "\t" + str(needle_microrray[0]) + "\t" + str(needle_mse[1]) + "\n")
    o.write("Needle (w=23)\t"+ str(needle23_seqc[0]) + "\t" + str(needle23_microrray[0]) + "\t" + str(needle23_mse[1]) + "\n")
    o.write("Needle (w=39)\t"+ str(needle39_seqc[0]) + "\t" + str(needle39_microrray[0]) + "\t" + str(needle39_mse[1]) + "\n")
