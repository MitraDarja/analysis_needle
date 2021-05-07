import statistics
import sys

reindeer_file = sys.argv[1]
fasta_file = sys.argv[2]
out_file = sys.argv[3]

genes = {}
# Necessary, so the complete fasta id can be stored
with open(fasta_file, 'r') as f:
    for line in f:
        if (line[0] == ">"):
            genes.update({line.split('|')[0] : line[:-1]})

with open(reindeer_file, 'r') as f:
    with open(out_file, 'w') as o:
        # One line gives the result for every dataset
        for line in f:
            estimations = []
            transcript_name = line.split('|')[0]
            if transcript_name in genes:
                # Each dataset result is divided by a ' ', so results contains all results for all datasets.
                # First entry in a line is the gene, so can be ignored.
                results = line.split()[1:]
                for res in results:
                    kmers_occ = [] # store the occurrences of the k-mers
                    # * means not found
                    if (res != '*'):
                        # Different counts are divided by ','
                         different_counts = res.split(',')
                         for counts in different_counts:
                             occurrences = []
                             # Entry looks like X-Y:Z is, where X and are the position of the k-mers and Z the count value
                             # kmers = Y-X, so number of k-mers with that count value
                             kmers = int(counts.split(':')[0].split('-')[1]) - int(counts.split(':')[0].split('-')[0])
                             occ = 0
                             # occ_s is the occurence (Z above), can be * if not found
                             occ_s = counts.split(':')[1]
                             if (occ_s != '*'):
                                 occ = int(occ_s)
                            # Add all occurences to kmers_occ
                             occurrences = [occ] * kmers
                             kmers_occ += occurrences
                    if (kmers_occ == []):
                        kmers_occ = [0]
                    # Take median of all k-mer counts as espression value for one dataset
                    estimations.append(int(statistics.median(kmers_occ)))
                # Write out the expression values
                o.write(genes[transcript_name])
                o.write('\t')
                for e in estimations:
                    o.write(str(int(e)))
                    o.write('\t')
                o.write('\n')
