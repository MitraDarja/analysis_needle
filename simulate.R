library(Biostrings)
library("edgeR")
library("polyester")

# Script to create the simulated data set

# Data are the proteincoding transcript sequences from https://www.gencodegenes.org/human/
fastapath = "data/gencode.v36.pc_transcripts.fa.gz"

fastafile = system.file('extdata', fastapath, package='polyester')
fasta = readDNAStringSet(fastapath)

trials = 20
coverage = 40
bp_length = 75
howManyGenes = 0.1 * length(fasta) # 10 % of transcripts are differently expressed
for (j in (1:trials)) {
  out = c()
  # Pick randomly transcripts, which should be differently expressed and then add randomly a fold change value to them
  whichGenes = sample(1:length(fasta),howManyGenes)
  amount = sample(c(0.33,0.5,2,3),howManyGenes,replace = TRUE) # Fold change randomly picked
  fold_changes = matrix(c(rep(1,length(fasta)),rep(1,length(fasta))), nrow=length(fasta))

  for(i in (1:howManyGenes)){
      out = rbind(out, c(names(fasta[whichGenes[i]]),amount[i]))
      fold_changes[whichGenes[i],1] = amount[i]
  }


  write.table(out, file=paste('data/Test_',as.character(j),'.tsv', sep=""),col.names = FALSE, row.names = FALSE, quote=FALSE, sep='\t')

  # Calculate coverage
  readspertx = round(coverage * width(fasta) / bp_length)
  d = paste('data/Test_',as.character(j),sep = "")
  dir.create(d)
  simulate_experiment(fastapath, reads_per_transcript=readspertx, num_reps = c(1,1), outdir = d,
                      fold_changes=fold_changes, seed=142, readlen = bp_length, paired = TRUE, gzip = TRUE)

  fasta = readDNAStringSet(fastapath)
}
