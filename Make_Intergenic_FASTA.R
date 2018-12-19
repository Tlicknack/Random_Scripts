#DONE

#This script creates a fasta file of only intergenic sequences, gives names with scaffold and genes surrounding the sequence

#INPUT: 
  #genic GFF
  #FASTA file
#OUTPUT
  #intergenic FASTA


#source("http://bioconductor.org/biocLite.R")
#biocLite("seqinr")
library("seqinr")

#species_fasta = read.fasta("/N/u/tlicknac/Carbonate/Paramecium_FASTA/biaurelia_V1-4_assembly_v1.fa", as.string=TRUE)
species_fasta = read.fasta("/N/u/tlicknac/Carbonate/Paramecium_FASTA/ptetraurelia_mac_51.fa", as.string=TRUE)
#species_fasta = read.fasta("/N/u/tlicknac/Carbonate/Paramecium_FASTA/sexaurelia_AZ8-4_assembly_v1.fasta", as.string=TRUE)
#species_fasta = read.fasta("/N/u/tlicknac/Carbonate/Paramecium_FASTA/pdecaurelia_mac_223_v1.0.fa", as.string=TRUE)
species_fasta = read.fasta("/N/u/tlicknac/Carbonate/Paramecium_FASTA/caudatum_43c3d_assembly_v1.fa", as.string=TRUE)


#species_intergenic_gff = read.table("/N/u/tlicknac/Carbonate/Paramecium_GFF/pbi-intergenic.tab", skip=2)   #gff's ended up with first row of NAs
species_intergenic_gff = read.table("/N/u/tlicknac/Carbonate/Paramecium_GFF/ptet-intergenic.tab", skip=2)   
#species_intergenic_gff = read.table("/N/u/tlicknac/Carbonate/Paramecium_GFF/psex-intergenic.tab", skip=2)   
#species_intergenic_gff = read.table("/N/u/tlicknac/Carbonate/Paramecium_GFF/pdec-intergenic.tab", skip=2)   
species_intergenic_gff = read.table("/N/u/tlicknac/Carbonate/Paramecium_GFF/pcaud-intergenic.tab", skip=2)   

len_intergenic_gff = nrow(species_intergenic_gff)

for(k in 1:len_intergenic_gff){                                                        #iterate through intergenic gff                                      
  #get what i need from gff
  scaf = as.character(species_intergenic_gff[[1]][k])                                                #scaf is first column                            
  start_seq = species_intergenic_gff[[2]][k]                                          
  end_seq = species_intergenic_gff[[3]][k]
  annotation = paste(">", scaf, "|5'+_", species_intergenic_gff[[4]][k], "_5'-_", species_intergenic_gff[[5]][k], "_3'+_", species_intergenic_gff[[6]][k], "_3'-_", species_intergenic_gff[[7]][k], sep="")
  #get what i need from fasta
  scaf_seq = species_fasta[[scaf]][1]                                                                                                                       
  intergenic_seq = as.character(substr(scaf_seq, start_seq, end_seq))                                                                                                     
  cat(annotation, "\n", intergenic_seq, "\n", file="/N/u/tlicknac/Carbonate/Paramecium_FASTA/Intergenic_FASTA/ptet-intergenic.fasta", append=TRUE)  #CHANGE EACH TIME
  cat("Row completed: ", k, "\n")
}
