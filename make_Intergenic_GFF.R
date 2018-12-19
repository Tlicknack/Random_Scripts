#DONE

#This script takes an annotation (GFF) and makes a modified GFF with intergenic coordinates

#INPUT: GFF
#OUTPUT: intergenic GFF

library("seqinr")

input_file = "/N/u/tlicknac/Carbonate/Paramecium_GFF/pcaud-gene.tab"
gff = read.table(input_file, sep="\t", as.is=T)

vscaf = unique(gff$V1)                        #vector of all scaffolds                                           
dffinal = data.frame(matrix(ncol=7))                                                  
colnames(dffinal) = c("scaffold", "start_position", "end_position", "5'+strand", "5'-strand", "3'+strand", "3'-strand")
dffinal_overlap = data.frame(matrix(ncol=7)) 
colnames(dffinal_overlap) = c("scaffold", "start_position", "end_position", "5'+strand", "5'-strand", "3'+strand", "3'-strand")

for(scaf in vscaf){
  	cat("Now working on scaffold: ", scaf, "\n")
  	ts = gff[which(gff$V1==scaf), ]       			#data frame of gff with only scaf of interest                                                                                          
	
  	#First gene in each scaffold
  	five_minus_1 = NA                                                      
  	five_plus_1 = NA
  	three_minus_1 = NA
  	three_plus_1 = NA
  
  	row1 = ts[1,]                        #curent row                                   
  	start_row1 = as.numeric(row1[2])     #start position                                            
  	#end_row1 = as.numeric(row1[3])       #end position                                                
  	row1_orient = as.character(row1[4])  #gene orientation                                                   
  	row1_id = as.character(row1[5])      #geneID
	
	if(row1_orient == "+"){               #if gene is going left to right.... this sequence is the 5' upstream
		five_plus_1 = row1_id	
	} else{                               #if gene goes right to left.... this sequence is 3' downstream
		three_minus_1 = row1_id	
	}

	dfrow1 = data.frame(matrix(c(ts[1,1], start_row1-200, start_row1, five_plus_1, five_minus_1, three_plus_1, three_minus_1), byrow=F, ncol=7))   #Arbitrarily make 200nts the cutoff
	colnames(dfrow1) = c("scaffold", "start_position", "end_position", "5'+strand", "5'-strand", "3'+strand", "3'-strand")
	dffinal = rbind(dffinal, dfrow1)      #bind final df to first row

	#Middle Genes
	i=2
	while(i < nrow(ts)-1){
	  cat("Intergenic region number: ", i, "\n")
	  five_minus_i = NA                                                      
  	five_plus_i = NA
  	three_minus_i = NA
  	three_plus_i = NA
	
	  rowi = ts[i-1,]                           #first row
	  rowi2 = ts[i,]                            #second row
	  end_rowi = as.numeric(rowi[3])
	  start_rowi2 = as.numeric(rowi2[2])
	  rowi_orient = as.character(rowi[4])
	  rowi2_orient = as.character(rowi2[4])
	  rowi_id = as.character(rowi[5])
	  rowi2_id = as.character(rowi2[5])
		
	  if(rowi_orient == "-"){                 #sequence between 1st and 2nd row.... if first row gene is on minus strand, then sequence is 5' upstream
      five_minus_i = rowi_id
    }
    if(rowi_orient == "+"){                 #if first row gene is on plus strand, then sequence is 3' downstream
      three_plus_i = rowi_id 
    }
    if(rowi2_orient == "-"){
      three_minus_i = rowi2_id
    }
    if(rowi2_orient == "+"){
      five_plus_i = rowi2_id
    }	
		
	  if(start_rowi2 > end_rowi){
	    df_nonoverlapping = data.frame()  #initiate data frame to add to first row
		  df_nonoverlapping = data.frame(matrix(c(ts$V1[i], end_rowi, start_rowi2,  five_plus_i, five_minus_i, three_plus_i,  three_minus_i), byrow=F,  ncol=7 )) #df: scaf, start/end of intergenic, genes
		  colnames(df_nonoverlapping) = c("scaffold", "start_position", "end_position", "5'+strand", "5'-strand", "3'+strand", "3'-strand")
		  dffinal = rbind(dffinal, df_nonoverlapping)
		  
	  } else{
	    df_overlapping = data.frame()
	    df_overlapping = data.frame(matrix(c(ts$V1[i], start_rowi2, end_rowi,  five_plus_i, five_minus_i, three_plus_i,  three_minus_i), byrow=F,  ncol=7 ))
	    colnames(df_overlapping) = c("scaffold", "start_position", "end_position", "5'+strand", "5'-strand", "3'+strand", "3'-strand")
	    dffinal_overlap = rbind(dffinal_overlap, df_overlapping)
	  }
	  i=i+1
	}
	
	#Last Gene
	five_minus_n = NA                                                      
	five_plus_n = NA
	three_minus_n = NA
	three_plus_n = NA
	
	rown = ts[nrow(ts),]                                                  
	#start_rown = as.numeric(rown[2])                                               
	end_rown = as.numeric(rown[3])                                                   
	rown_orient = as.character(rown[4])                                                    
	rown_id = as.character(rown[5])     
	
	if(rown_orient == "+"){
	  three_plus_n = rown_id
	} else{
	  five_minus_n = rown_id
	}
	dfrown = data.frame(matrix(c(ts[nrow(ts),1], end_rown, end_rown+200, five_plus_n, five_minus_n, three_plus_n, three_minus_n), byrow=F, ncol=7))   #Arbitrarily make 200nts the cutoff
	colnames(dfrown) = c("scaffold", "start_position", "end_position", "5'+strand", "5'-strand", "3'+strand", "3'-strand")
	dffinal = rbind(dffinal, dfrown)      #bind final df to first row
}	

dffinal = dffinal[-1,]
dffinal_overlap = dffinal_overlap[-1,]

write.table(dffinal, file="pcaud-intergenic.tab", sep="\t", row.names=F)                  #15157 intergenic regions  
write.table(dffinal_overlap, file="pcaud-genic_overlap.tab", sep="\t", row.names=F)       #3781 regions of genic overlap ....   19.9%	
