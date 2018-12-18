#DONE

#This script will create fasta files for the upstream region of 50 genes 
  #This was originally written for bins of genes based on expression level
  #ONLY Ptet

library(seqinr)

ptet_fa = read.fasta("/N/u/tlicknac/Carbonate/Paramecium_FASTA/ptetraurelia_mac_51.fa", as.string=T, forceDNAtolower=F)
ptet_gff = read.table("/N/u/tlicknac/Carbonate/Paramecium_GFF/ptet-gene.tab", as.is = T, sep="\t")

high = read.table("ptet-top50_expressed_genes.txt", header=F, sep=",")
high = high$V1
middle = read.table("ptet-mid50_expressed_genes.txt", header=F, sep=",")
middle = middle$V1
low = read.table("ptet-low50_expressed_genes.txt", header=F, sep=",")
low = low$V1

#-----
vlabels = c()
vseqs = c()

for(i in 1:length(high)){
  rowi = ptet_gff[which(ptet_gff$V5 == high[i]),]
  scaf_seq = as.character(ptet_fa[which(getName(ptet_fa) == rowi$V1)])
  
  if(rowi$V4 == "+"){
    seq = tolower(substr(scaf_seq, rowi$V2-200, rowi$V2))
  }
  if(rowi$V4 == "-"){
    seq = substr(scaf_seq, rowi$V3, rowi$V3+200)
    seq = tolower(c2s(rev(comp(s2c(seq)))))
  }
  
  label = paste(rowi$V5, "|", rowi$V1, "|", rowi$V4, sep="")
  vlabels = append(vlabels, label)
  vseqs = append(vseqs, seq)
}

write.fasta(sequences = as.list(vseqs), names = vlabels, file.out = "high_expression_ptet.fasta")

#-----
vlabels = c()
vseqs = c()

for(i in 1:length(middle)){
  rowi = ptet_gff[which(ptet_gff$V5 == middle[i]),]
  scaf_seq = as.character(ptet_fa[which(getName(ptet_fa) == rowi$V1)])
  
  if(rowi$V4 == "+"){
    seq = tolower(substr(scaf_seq, rowi$V2-200, rowi$V2))
  }
  if(rowi$V4 == "-"){
    seq = substr(scaf_seq, rowi$V3, rowi$V3+200)
    seq = tolower(c2s(rev(comp(s2c(seq)))))
  }
  
  label = paste(rowi$V5, "|", rowi$V1, "|", rowi$V4, sep="")
  vlabels = append(vlabels, label)
  vseqs = append(vseqs, seq)
}

write.fasta(sequences = as.list(vseqs), names = vlabels, file.out = "mid_expression_ptet.fasta")

#-----
vlabels = c()
vseqs = c()

for(i in 1:length(low)){
  rowi = ptet_gff[which(ptet_gff$V5 == low[i]),]
  scaf_seq = as.character(ptet_fa[which(getName(ptet_fa) == rowi$V1)])
  
  if(rowi$V4 == "+"){
    seq = tolower(substr(scaf_seq, rowi$V2-200, rowi$V2))
  }
  if(rowi$V4 == "-"){
    seq = substr(scaf_seq, rowi$V3, rowi$V3+200)
    seq = tolower(c2s(rev(comp(s2c(seq)))))
  }
  
  label = paste(rowi$V5, "|", rowi$V1, "|", rowi$V4, sep="")
  vlabels = append(vlabels, label)
  vseqs = append(vseqs, seq)
}

write.fasta(sequences = as.list(vseqs), names = vlabels, file.out = "low_expression_ptet.fasta")

