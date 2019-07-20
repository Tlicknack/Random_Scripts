# DONE
# Get upstream sequence for a list of genes of interest


library(seqinr)
setwd("/media/tlicknac/Seagate Expansion Drive/Candidate-Promoters/Sharp-Broad-Pdec_Promotes/")

gff = read.table("/media/tlicknac/Seagate Expansion Drive/Paramecium_Annotations/ptet-gene.tab", skip = 2, as.is = T)
fasta = read.fasta("/media/tlicknac/Seagate Expansion Drive/Paramecium_Assemblies/ptetraurelia_mac_51.fa", as.string = T, forceDNAtolower = T)
genes = c("PTET.51.1.P0240290", "PTET.51.1.P0010170", "PTET.51.1.P0250075", "PTET.51.1.P1060020", "PTET.51.1.P1170148")


for(geneID in genes){
  gffRow = gff[which(gff$V5 == geneID),]
  scafSeq = as.character(fasta[which(getName(fasta) == gffRow$V1)])
  
  if(gffRow$V4 == "+"){                                       
    upstream = substr(scafSeq, gffRow$V2-200, gffRow$V2) 
  } 
  if(gffRow$V4 == "-"){
    upstream = substr(scafSeq, gffRow$V3, gffRow$V3+200) 
    upstream = c2s(rev(comp(s2c(upstream))))                      
  }
  write.fasta(sequences = upstream, names = paste(geneID, "_", gffRow$V1, "-", gffRow$V2, ":", gffRow$V3, sep = ""),file.out = paste(geneID, "_Upstream.fasta", sep = "") ) 
}