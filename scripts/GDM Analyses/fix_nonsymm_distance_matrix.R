files = list.files(path="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/BY_SPECIES/",
                   pattern="DISTANCES.txt.converted.txt",
                   full.names = T,recursive = T)

#files="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/BY_SPECIES/CONVERTED/NITENS/NITENS_Idw_abundance_clipped.ascAGGFACT-1.MORPH-AND-GENE-DISTANCES.txt.converted.txt"

for (string in files) {
  print(string)
  
  #string = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/BY_SPECIES/CONVERTED/BELLII/BELLII_Idw_abundance_clipped.ascAGGFACT-1.MORPH-AND-GENE-DISTANCES.txt.converted.txt"
  symmdf = read.csv(string)
  if (ncol(symmdf) == 1) {
    symmdf = read.csv(string,sep="\t")
  }
  rownames(symmdf) = make.unique(as.character(symmdf[,1]))
  colnames(symmdf) = make.unique(colnames(symmdf))
  symmdf = symmdf[,-1]
  fixRowCols = function(symmdf) {
    fixdf = symmdf
    
    print("ORIGINAL DIM:")
    print(dim(fixdf))
    finish = F
    while (finish == F) {
      ## get set of cols
      colset = unique(sort(colnames(fixdf)))
      ## get set of rows
      rowset = unique(sort(rownames(fixdf)))
      ## which cols not in rows
      colnotrow = setdiff(colset,rowset)
      ## which rows not in cols
      rownotcol = setdiff(rowset,colset)
      
      if(length(colnotrow) == 0 && length(rownotcol) == 0) {
        finish = T
      }
      
      
      if(length(rownotcol) != 0) {
        ## for rows not in cols, make col with row name
        for (row in rownotcol) {
          #print(row)
          fixdf = cbind(fixdf,NA)
          colnames(fixdf)[ncol(fixdf)] = row
        }
        print("NEW DIM 1:")
        print(dim(fixdf))
      }
      
      if(length(rownotcol) != 0) {
        ## for cols not in rows, make row with col name
        for (col in colnotrow) {
          #print(col)
          fixdf = t(fixdf)
          fixdf = cbind(fixdf,NA)
          colnames(fixdf)[ncol(fixdf)] = col
          fixdf = t(fixdf)
        }
        print("NEW DIM 2:")
        print(dim(fixdf))
      }
      
      ## verify that the row and col nums are identical
      if (dim(fixdf)[1] == dim(fixdf)[2]) {
        ## verify that the row and col sets are identical
        if ( length(setdiff(rownames(fixdf),colnames(fixdf)))==0  && 
             length(setdiff(colnames(fixdf),rownames(fixdf)))==0 ) {
          finish = T
        } else {
          print("RESTART1")
          fixdf = t(fixdf)
        }
      } else {
        print("RESTART2")
        fixdf=t(fixdf)
      }
    }
    ## sort the df by rows and cols so that they have the same order
    fixdf = fixdf[order(rownames(fixdf)) , order(colnames(fixdf))]
    ## identify the rows that are empty
    emptyrows = which(rowSums(is.na(fixdf))==ncol(fixdf))
    ## identify the cols that are empty
    emptycols = which(colSums(is.na(fixdf))==nrow(fixdf))
    ## df [i,j] = df[j,i] unless i=j then df[i,j] = 0
    for(i in emptyrows){
      #print(i)
      for(j in 1:ncol(fixdf)) {
        #print(j)
        if (i == j) {
          fixdf[i,j] = 0
        } else {
          fixdf[i,j] = fixdf[j,i]
        }
      }
    }
    for(j in emptycols){
      #print(j)
      for(i in 1:nrow(fixdf)) {
        #print(j)
        if (j == i) {
          fixdf[i,j] = 0
        } else {
          fixdf[i,j] = fixdf[j,i]
        }
      }
    }
    ## check whether top and bottom are the same
    for (i in 1:nrow(fixdf)) {
      for(j in 1:ncol(fixdf)) {
        
        ijna=F
        jina=F
        
        if(is.na(fixdf[i,j])) {
          #print("IS NA [i,j]")
          ijna = T
          
        }
        if(is.na(fixdf[j,i])) {
          #print("IS NA [j,i]")
          jina = T
        }
        
        if(ijna==T && jina==F) {
          fixdf[i,j] = fixdf[j,i]
        }
        else if(ijna==F && jina==T) {
          fixdf[j,i] = fixdf[i,j]
        }
        else if (ijna==T && jina==T) {
          #print(paste("NA [",i,",",j,"]",sep=""))
          next
        } else {
          
          if(fixdf[i,j] != fixdf[j,i]) {
            
            #print(paste("BAD [",i,",",j,"]",sep=""))
            next
            
          }  
          
        }
        
      }
    }
    return(fixdf)
  }
  
  fixdf = fixRowCols(symmdf)
  write.csv(fixdf,string)
  
}
