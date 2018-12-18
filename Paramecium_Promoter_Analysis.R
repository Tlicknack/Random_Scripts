#This script will look at the upstream region of genes and determine their nt composition

library(seqinr)

setwd("All_Aurelias_Upstream/")
all_files = list.files(".")
all_files = all_files[-1]
df_gc = data.frame(matrix(ncol=2))
colnames(df_gc) = c("Gene", "GC")
counter = 0

for(file in all_files){
  seqs = read.alignment(all_files[2], forceToLower=F, format = "fasta")
  
  for(i in 1:length(seqs$seq)){
    name = seqs$nam[i]
    gc_content = GC(s2c(seqs$seq[i]))
    df_gc[(counter),] = c(name, gc_content)
    counter=counter+1
  }
}

