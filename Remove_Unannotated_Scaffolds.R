# Modify fasta files to remove scaffolds that arent included in the gff annotation

fasta = read.fasta("~/Paramecium_FASTA/pdecaurelia_mac_223_v1.0.fa", as.string=T, forceDNAtolower=T)  
gff = read.table("~/Paramecium_GFF/pdec-full.gff", header=F, sep="\t")

fasta_mod = fasta[intersect(getName(fasta), gff$V1)]

write.fasta(names=as.list(getName(fasta_mod)), sequences=getSequence(fasta_mod), file.out="~/Paramecium_FASTA/No_Annotation_FASTA/pdec_mod.fa")