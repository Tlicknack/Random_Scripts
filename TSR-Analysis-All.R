source("/N/dc2/projects/ParameciumPromoters/Para-promoters-analysis/scripts/buildTSSDb.R")
source("/N/dc2/projects/ParameciumPromoters/Para-promoters-analysis/scripts/tssToList.R")
sampleDb <- buildTSSDb(speciesNames=c("Pdec","Poct","Pnov", "Pjen"), "/N/dc2/scratch/tlicknac/Data-TSR/TSR-Pdec.txt", "/N/dc2/scratch/tlicknac/Data-TSR/TSR-Poct.txt", "/N/dc2/scratch/tlicknac/Data-TSR/TSR-Pnov.txt", "/N/dc2/scratch/tlicknac/Data-TSR/TSR-Pjen.txt")
  # First layer of sampleDb is each species name
  # Within each species name is each gene for which we have a TSR
  # For each gene, we have TSR info on the strand, nTSSs, nTAGs, tsrPeak, tsrWdth, tsrSI, tsrMSI

source("/N/dc2/projects/ParameciumPromoters/Para-promoters-analysis/scripts/loadOrthoParalogs.R")
orthoPara <- loadOrthoParalogs("/N/dc2/projects/ParameciumPromoters/Para-promoters-analysis/paralog_datasets/all_aurelias-cut-mod.poff")
  # First layer is each WGD family
  # Each WGD family has all 14 aurelia spp names
  # Each species has 1 potential paralogs

source("/N/dc2/projects/ParameciumPromoters/Para-promoters-analysis/scripts/buildGeneDb.R")
geneDb <- buildGeneDb(c("Pdec", "Pjen", "Pnov", "Poct"), "~/Paramecium_GFF/pdec-full.gff", "~/Paramecium_GFF/pjen-full.gff", "~/Paramecium_GFF/pnov-full.gff", "~/Paramecium_GFF/poct-full.gff")

source("/N/dc2/projects/ParameciumPromoters/Para-promoters-analysis/scripts/loadWGD3paralogs.R")
wgd3paralogs = loadWGD3paralogs(c("Pdec", "Pjen", "Poct", "Pnov"), "~/Paramecium_WGD3_Trees/pdecaurelia_223_annotation_v1.0.WGD.tree", "~/Paramecium_WGD3_Trees/pjenningsi_M_annotation_v1.0.WGD.tree", "~/Paramecium_WGD3_Trees/poctaurelia_K8_annotation_v1.0.WGD.tree", "~/Paramecium_WGD3_Trees/pnovaurelia_TE_annotation_v1.0.WGD.tree")





#############################################################################################

## buildTSSDb.R"
## Builds a TSS database from selected paramecium species and returns an .RData object that can be accesed programatically. 

buildTSSDb <- function(speciesNames, ...) {
            if (is.character(speciesNames)==FALSE) {
               stop("speciesNames must be of class 'character'")
               }
            sp.len <- length(speciesNames)
            x.1 <- list(...)
            print(x.1)
            x <- length(x.1)
            my.seq <- 1:x
            my.list <- vector(mode="list", length=x)
            my.slots <- speciesNames
            #print(my.slots) #for debugging
            names(my.list) <- my.slots
            #print(names(my.list)) #for debugging
                 for (i in 1:x) {
                      this.table  <- read.table(file=x.1[[i]][1], header=TRUE)
                      genes.table <- this.table[complete.cases(this.table[,12]),]
                      for (j in 1:nrow(genes.table)) {
                          my.gene <- genes.table[j,12]
                          my.list[[i]][j] <- tssToList(genes.table, my.gene)
                          print(j)
                          }
                       names(my.list[[i]]) <- genes.table[,12]
              }

           return(my.list)
}

## tssToList.R
## Creates a list from the contents of a standard TSRsetMerged.txt file from TSRchitect

tssToList <- function(x, geneID) {
           library(GenomicRanges)
           this.table  <- x
           this.list <- vector(mode="list", length=1)
           this.list.2 <- vector(mode="list", length=9)
           names(this.list) <- geneID
           names(this.list.2) <- c("seq", "coords", "strand", "nTSSs", "nTAGs", "tsrPeak", "tsrWdth", "tsrSI", "tsrMSI")
           this.gene <- this.table[this.table$featureID==geneID,]
                        if (nrow(this.gene) == 1) { 
                        this.list.2$seq <- this.gene$seq
                        this.list.2$coords <- IRanges(start=this.gene$start, end=this.gene$end)
                        this.list.2$strand <- this.gene$strand
                        this.list.2$nTSSs <- this.gene$nTSSs
                        this.list.2$nTAGs <- this.gene$nTAGs
                        this.list.2$tsrPeak <- this.gene$tsrPeak
                        this.list.2$tsrWdth <- this.gene$tsrWdth
                        this.list.2$tsrSI <- this.gene$tsrSI
                        this.list.2$tsrMSI <- this.gene$tsrMSI
                        this.list[[1]] <- this.list.2
                        }
           return(this.list)
}

## loadOrthoParalogs.R
## Parses the poff file and deposits the results in a ParaOrtho object

loadOrthoParalogs <- function(OPtable) {
  if (is.character(OPtable)==FALSE) {
    stop("Argument OPtable must be of class 'character'")
  }
  message("Importing ortho-paralog table.")
  myOP.table <- read.table(file=OPtable, header=TRUE)
  sppVec <- c("WGD_ID", "nbGenes", "Pbi", "Pdec", "Pdodec", "Pjenn", "Pnov", "Poct", "Ppent", "Pprim", "Pquad", "Psept", "Psex", "Pson", "Ptet", "Ptre", "Pcaud", "Pmult")
  colnames(myOP.table) <- sppVec
  spp.table <- myOP.table[,-1:-2]
  spp.table <- spp.table[,-15:-16]
  sppIds <- c("Pbi", "Pdec", "Pdodec", "Pjenn", "Pnov", "Poct", "Ppent", "Pprim", "Pquad", "Psept", "Psex", "Pson", "Ptet", "Ptre")
  colnames(myOP.table) <- sppVec          
  opList <- vector(mode="list", length=nrow(myOP.table))
  names(opList) <- myOP.table[,1]
  for (i in 1:nrow(myOP.table)) {
    this.id <- myOP.table[i,1]
    this.vec <- spp.table[i,]
    #print(head(this.vec))
    this.list <- strsplit(as.character(unlist(this.vec)), split=",")
    this.vector <- unlist(this.list)
    small.list <- vector(mode="list", length=14)
    names(small.list) <- sppIds
    this.ma <- matrix(this.vector, nrow=14, ncol=2, byrow=TRUE)
    for (j in 1:length(sppIds)) {
      this.spp <- sppIds[j]
      small.vec <- vector(mode="character", length=2)
      names(small.vec) <- c("Para1","Para2")
      small.list[[j]] <- small.vec
      #print(j)
      small.list[[j]] <- this.ma[j,]
      names(small.list[[j]]) <- c("Para1", "Para2")
      #print(small.vec)
      #print(this.ma[j,])
    }
    opList[[i]] <- small.list
  }
  
  return(opList)
}

## buildGeneDb.R
## Builds a database of genes selected paramecium species and returns a list object

library(rtracklayer)

buildGeneDb <- function(speciesNames, ...) {
  if (is.character(speciesNames)==FALSE) {
    stop("speciesNames must be of class 'character'")
  }
  
  message("Loading gene annotations to object.")
  
  sp.len <- length(speciesNames)
  x.1 <- list(...)
  print(x.1)
  x <- length(x.1)
  
  if (x != sp.len) {
    stop("speciesNames and gff files must be of the same length")
  }
  
  my.seq <- 1:x
  my.list <- vector(mode="list", length=x)
  my.slots <- speciesNames
  #print(my.slots) #for debugging
  names(my.list) <- my.slots
  #print(names(my.list)) #for debugging
  for (i in 1:x) {
    this.gff  <- readGFF(filepath=x.1[[i]])
    this.df <- this.gff[this.gff$type=="gene",]
    genes.table <- this.df$ID
    my.list[[i]] <- genes.table
  }
  
  return(my.list)
}

# loadWGD3paralogs.R (me) (DONE)
#  Builds a database of genes from the 3rd most recent WGD

loadWGD3paralogs <- function(speciesNames, ...){
  if (is.character(speciesNames)==FALSE) {
    stop("speciesNames must be of class 'character'")
  }
  x.1 = list(...)                                                             #make input a list

  if(length(speciesNames) != length(x.1)){
    stop("speciesName and wgd files must be of the same length")
  }
  message("Importing WGD3 paralog table.")

  wgdDb = vector(mode="list", length = length(speciesNames))                  # FINAL LIST
  names(wgdDb) = speciesNames
  
  for(tree in x.1){                                # Iterate through each species given
    thisWGD = read.table(file=tree, header=T)                  # read in the WGD3 data frame... should have many thousand rows containing 1-8 genes per row
    thisWGD$NB = NULL                                              # remove the first, useless row
    tmpDb = vector(mode="list", length = nrow(thisWGD))            # make a temporary list that will be attached to each species's object
    
    wgdNames = replicate(n= nrow(thisWGD), expr= "WGD_x")          # make names for each WGD gene family for each i'th species: wgdDb$Pdec$WGD1
    for(q in 1:length(wgdNames)){
      wgdNames[q] = gsub("_x", as.character(q), wgdNames[q])
    }
    names(tmpDb) = wgdNames                                        # assign the vector of names to the temporary list
    
    for(j in 1:nrow(thisWGD)){                                       # Iterate through each row in the WGD3 data frame we imported
      print(j)
      thisWgdRow = thisWGD[j,]
      wgdGenes = thisWgdRow[which(nchar(as.vector(unlist(thisWgdRow))) > 1)]  # remove any item in each row that has less than 1 character (i.e. is just a ".")
      print(wgdGenes)
      tmpDb[[j]] = wgdGenes                                                   # Assign the remaining gene names to the j'th oject in the list (i.e. j=1, tmpDb$WGD1)
      
      tmpNames = replicate(n=length(tmpDb[[j]]), expr="ID_X")     # make new names for ID's in each WGD3 gene family, as sometimes only ID7 + ID8 remain, and that number is meaningless alone
      for(r in 1:length(tmpDb[[j]])){
        tmpNames[r] = gsub("_X", as.character(r), tmpNames[r])
      }
      names(tmpDb[[j]]) = tmpNames                                # assign this vector of ID names to each WGD3 gene family object
    }
    wgdDb[[i]] = tmpDb
  }
  return(wgdDb)
}

########################
# loadWGD3paralogs.R (me) (DONE)
# Slightly different than the above fxn. this one makes each species have objects as gene names, and each gene name contains vector of WGD3 paralogs
loadWGD3paralogs <- function(speciesNames, ...){
  if (is.character(speciesNames)==FALSE) {
    stop("speciesNames must be of class 'character'")
  }
  x.1 = list(...)                                                             #make input a list
  
  if(length(speciesNames) != length(x.1)){
    stop("speciesName and wgd files must be of the same length")
  }
  message("Importing WGD3 paralog table.")
  
  wgdDb = vector(mode="list", length = length(speciesNames))                  # FINAL LIST
  names(wgdDb) = speciesNames
  
  for(tree in x.1){                                # Iterate through each species given
    thisWGD = read.table(file=tree, header=T)                  # read in the WGD3 data frame... should have many thousand rows containing 1-8 genes per row
    thisWGD$NB = NULL                                              # remove the first, useless row
    tmpDb = vector(mode="list")            # make a temporary list that will be attached to each species's object
    
    for(j in 1:nrow(thisWGD)){                                       # Iterate through each row in the WGD3 data frame we imported
      thisWgdRow = thisWGD[j,]
      wgdGenes = thisWgdRow[which(nchar(as.vector(unlist(thisWgdRow))) > 1)]  # remove any item in each row that has less than 1 character (i.e. is just a ".")
      
      for(wgdGene in wgdGenes){
        currentGene = as.character(unlist(wgdGene))
        otherGenes = as.character(unlist(wgdGenes[which(wgdGenes != currentGene)]))
        tmpDb[[currentGene]] = otherGenes
      }
    }
    wgdDb[[i]] = tmpDb
  }
  return(wgdDb)
}