import statistics
import sys

reindeer_file = sys.argv[1]
fasta_file = sys.argv[2]
out_file = sys.argv[3]

genes = {}
with open(fasta_file, 'r') as f:
    for line in f:
        if (line[0] == ">"):
            genes.update({line.split('|')[0] : line[:-1]})

with open(reindeer_file, 'r') as f:
    with open(out_file, 'w') as o:
        for line in f:
            estimations = []
            transcript_name = line.split('|')[0]
            if transcript_name in genes:
                results = line.split()[1:]
                for res in results:
                    kmers_occ = []
                    if (res != '*'):
                         different_counts = res.split(',')
                         for counts in different_counts:
                             occurrences = []
                             kmers = int(counts.split(':')[0].split('-')[1]) - int(counts.split(':')[0].split('-')[0])
                             occ = 0
                             occ_s = counts.split(':')[1]
                             if (occ_s != '*'):
                                 occ = int(occ_s)
                             occurrences = [occ] * kmers
                             kmers_occ += occurrences
                    if (kmers_occ == []):
                        kmers_occ = [0]
                    estimations.append(int(statistics.median(kmers_occ)))
                o.write(genes[transcript_name])
                o.write('\t')
                for e in estimations:
                    o.write(str(int(e)))
                    o.write('\t')
                o.write('\n')
