#This script will return basic genome composition for intergenic regions of Aurelia genome
#INPUT: Paramecium-intergenic.fasta
#OUTPUT: table of genicID with GC content and length

#source("http://bioconductor.org/biocLite.R")
#biocLite("seqinr")
library("seqinr")

species_intergenic_fasta = read.fasta("Paramecium_FASTA/ptet-intergenic.fasta", as.string=T, forceDNAtolower=T)
scafs = getName(species_intergenic_fasta)
final_table = data.frame(matrix(nrow=length(scafs), ncol=3), row.names=NULL)
colnames(final_table) = c("Scaffold+ID", "LengthOfIntergenic", "GCcontent")

for(i in 1:length(scafs)){
  seq = getSequence(species_intergenic_fasta[[i]])
  len_seq = length(seq)
  gc_seq = round(GC(seq), digits=3)
  info = c(scafs[i], len_seq, gc_seq)
  final_table[i,] = info
  cat("Row completed: ", i, "\n")
}
final_table = final_table[complete.cases(final_table),]

#Now with the table...
summary_GC = summary(as.numeric(final_table$GCcontent))
summary_Len = summary(as.numeric(final_table$LengthOfIntergenic))
hist(as.numeric(final_table$LengthOfIntergenic), xlim = c(0,1000), density=5, breaks=1000, main = "Frequency of Intergenic Lengths in Ptet", xlab="Length (nts)", ylab="Number of Sequences")
hist(as.numeric(final_table$GCcontent), density=5, breaks=100, main = "Frequency of GC Content of Scaffold in Ptet", xlab = "GC%", ylab = "Number of Sequences")
scatter.smooth(x=as.numeric(final_table$LengthOfIntergenic), y=as.numeric(final_table$GCcontent))
cor(as.numeric(final_table$GCcontent), as.numeric(final_table$LengthOfIntergenic))  #tet=.006, sex=.41

col_lengths = as.numeric(final_table$LengthOfIntergenic)

seqs1 = final_table[which(col_lengths>0 & col_lengths<100),]
seqs2 = final_table[which(col_lengths>100 & col_lengths<200),]
seqs3 = final_table[which(col_lengths>200 & col_lengths<300),]
seqs4 = final_table[which(col_lengths>300 & col_lengths<400),]
seqs5 = final_table[which(col_lengths>400 & col_lengths<500),]
seqs6 = final_table[which(col_lengths>500 & col_lengths<600),]
seqs7 = final_table[which(col_lengths>600 & col_lengths<700),]
seqs8 = final_table[which(col_lengths>700 & col_lengths<800),]
seqs9 = final_table[which(col_lengths>800 & col_lengths<900),]
seqs10 = final_table[which(col_lengths>900 & col_lengths<1000),]
seqs11 = final_table[which(col_lengths>1000),]

#seqs = grep("^seqs*", ls(), value=T)   #chop down with seqs = seqs[-1] if seqs itself gets captured
#lseqs = list()
#for(seqi in seqs){
#  counter=1
#  df = data.frame(as.numeric(seqi["GCcontent"]))
#  rep = rep(x=counter, each=nrow(seqi))
#  df = cbind(df, rep)
#  colnames(df) = c("GCcontent", "Group")
#  lseqs[[seqi]] = df
#  counter=counter+1  
#}

df1 = data.frame(as.numeric(seqs1$GCcontent))
rep1 = rep(x=1, each=nrow(df1))
df1 = cbind(df1, rep1)
colnames(df1) = c("GCcontent", "Group")

df2 = data.frame(as.numeric(seqs2$GCcontent))
rep2 = rep(x=2, each=nrow(df2))
df2 = cbind(df2, rep2)
colnames(df2) = c("GCcontent", "Group")

df3 = data.frame(as.numeric(seqs3$GCcontent))
rep3 = rep(x=3, each=nrow(df3))
df3 = cbind(df3, rep3)
colnames(df3) = c("GCcontent", "Group")

df4 = data.frame(as.numeric(seqs4$GCcontent))
rep4 = rep(x=4, each=nrow(df4))
df4 = cbind(df4, rep4)
colnames(df4) = c("GCcontent", "Group")

df5 = data.frame(as.numeric(seqs5$GCcontent))
rep5 = rep(x=5, each=nrow(df5))
df5 = cbind(df5, rep5)
colnames(df5) = c("GCcontent", "Group")

df6 = data.frame(as.numeric(seqs6$GCcontent))
rep6 = rep(x=6, each=nrow(df6))
df6 = cbind(df6, rep6)
colnames(df6) = c("GCcontent", "Group")

df7 = data.frame(as.numeric(seqs7$GCcontent))
rep7 = rep(x=7, each=nrow(df7))
df7 = cbind(df7, rep7)
colnames(df7) = c("GCcontent", "Group")

df8 = data.frame(as.numeric(seqs8$GCcontent))
rep8 = rep(x=8, each=nrow(df8))
df8 = cbind(df8, rep8)
colnames(df8) = c("GCcontent", "Group")

df9 = data.frame(as.numeric(seqs9$GCcontent))
rep9 = rep(x=9, each=nrow(df9))
df9 = cbind(df9, rep9)
colnames(df9) = c("GCcontent", "Group")

df10 = data.frame(as.numeric(seqs10$GCcontent))
rep10 = rep(x=10, each=nrow(df10))
df10 = cbind(df10, rep10)
colnames(df10) = c("GCcontent", "Group")

df11 = data.frame(as.numeric(seqs11$GCcontent))
rep11 = rep(x=11, each=nrow(df11))
df11 = cbind(df11, rep11)
colnames(df11) = c("GCcontent", "Group")


mean1 = mean(as.numeric(seqs1$GCcontent))  #bi=13.66%
mean2 = mean(as.numeric(seqs2$GCcontent))  #bi=16.08%
mean3 = mean(as.numeric(seqs3$GCcontent))  #bi=17.15%
mean4 = mean(as.numeric(seqs4$GCcontent))  #bi=17.79%
mean5 = mean(as.numeric(seqs5$GCcontent))  #bi=18.46%
mean6 = mean(as.numeric(seqs6$GCcontent))  #bi=18.95%
mean7 = mean(as.numeric(seqs7$GCcontent))  #bi=19.4%
mean8 = mean(as.numeric(seqs8$GCcontent))  #bi=19.52%
mean9 = mean(as.numeric(seqs9$GCcontent))  #bi=20.09%
mean10 = mean(as.numeric(seqs10$GCcontent))  #bi=20.04%
mean11 = mean(as.numeric(seqs11$GCcontent))  #bi=22.37%

ttest = t.test(x=, y=)  #ttest has elements: statistics, parameter, p.value, conf.int, estimate, null.value, alternative, method, data.name

write.csv(final_table, file="pbi-intergenic-info.csv")
