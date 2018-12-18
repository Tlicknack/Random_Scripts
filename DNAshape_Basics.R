source("https://bioconductor.org/biocLite.R")
biocLite("GenomicRanges")
library(GenomicRanges)
biocLite("BSgenome")
library(BSgenome)
biocLite("DNAshapeR")
library(DNAshapeR)
biocLite("BSgenome.Scerevisiae.UCSC.sacCer3")


BSgenome::getBSgenome(available.genomes()[78])  #BSGenome
#Use BSgenome to get a genome of interest... in this case, 78 is S. cerivisiae

gr <- GRanges(seqnames = c("chrI") , strand = c("+", "-", "+"), ranges = IRanges(start = c(100, 200, 300), width = 100))   #GenomicRanges
#Create a GRanges object...
  #Seqnames is a vector of chromosomes you want
  #Strans is which strand you want
  #ranges is where on that scaffold you want to extract from... c(100, 200) means go to positions 100 and 200 
  #width is how big the DNA range should be
  
  #So taken together, the above code says: go to chromosome 1 positions 100, 200, and 300, extract 100 nts to the right on the +, -, + strandes respectively. Return that object

getFasta(gr, BSgenome = Scerevisiae, width = 100, filename = "tmp.fa")  #DNAshapeR

fn = "tmp.fa"  #DNAshapeR
pred = getShape(fn)  #DNAshapeR


#Example with Weibo's ribosomal motif
library(seqinr)
ribo = read.table("/home/tlicknac/Desktop/Paramecium_Genome_Data/Data/ribo_motif.txt", header = F, sep = "|", as.is = T)
genes = ribo$V1
seqs = strsplit(ribo$V2, split = " ")
seqs = sapply(seqs, "[[", 9)

write.fasta(sequences = as.list(seqs), names = genes, file.out = "ribo_motif_region.fa")
ribo_fasta = "ribo_motif_region.fa"
ribo_shape = getShape(ribo_fasta)
ribo = read.fasta(ribo_fasta)
ribo_seqs = paste(getSequence(ribo))
ribo_genes = getName(ribo)
motif_shapes = as.data.frame(ribo_shape$MGW[,12:19])

firstSeq = TRUE
for(i in 1:length(ribo_genes)){
  name = ribo_genes[i]
  score = motif_shapes[i,]
  sequence = substr(paste(ribo_seqs[[i]], collapse = ""), 12, 19)
  line = c(name, sequence, score)
  if(firstSeq == T){
    df_ribo_shape = data.frame(matrix(line, ncol=10))
    names(df_ribo_shape) = c("Gene", "Motif", "Pos1", "Pos2", "Pos3", "Pos4", "Pos5", "Pos6", "Pos7","Pos8")
    firstSeq=F
  } else{
      df_ribo = data.frame((matrix(line, ncol=10)))
      names(df_ribo) = c("Gene", "Motif", "Pos1", "Pos2", "Pos3", "Pos4", "Pos5", "Pos6", "Pos7","Pos8")
      df_ribo_shape = rbind(df_ribo_shape, df_ribo)
  }
}
#Have to fix this df if shape info isnt returned for each motif... this round, 11 motifs lacked info and slid the entire df down 11 rows
df = apply(df_ribo_shape, 2, as.character)
write.csv(df, file="ribo.csv", row.names=F)
#Make fixes
df_ribo_shape = read.csv("ribo.csv", header=T)
length(unique(df_ribo_shape$Motif))  #39 unique motifs




