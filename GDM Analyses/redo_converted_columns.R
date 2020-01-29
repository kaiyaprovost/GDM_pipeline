overwrite=T
files = list.files(path="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/BY_SPECIES/OTHERTEXT/",
                   #path="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/FULL_ENV/DISTANCES/",
                   pattern=".converted$",recursive = T,
                   full.names = T)
#files = c("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/BY_SPECIES/OTHERTEXT/CRISSALE/")

## done: bel bil bru cri 
#x <- file.info(files)
#files = files[match(1:length(files),rank(x$size))]

for (f in 1:length(files)) {
#for (f in 50:1) {

  print(paste(f,"of",length(files)))
  temp = files[f]
  
  print(temp)
  
  if (overwrite==F && (file.exists(paste(temp,".txt",sep="")))) {
    print("SKIP")
  } else {
    
    #temp = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/BY_SPECIES/BIL/TEXTFILES/test_converted.txt"
    
    df = read.csv(temp,sep=",",header=F)
    tabs=F
    if (ncol(df) == 1) {
      tabs=T
      df = read.csv(temp,sep="\t",header=F)
      
    }
    
    ## go through the columns -- check and see if the column name contains "\t"
    ## row 1, col 1
    endcol = ncol(df)
    toremove = c()
    print("COLUMNS")
    for (i in 2:endcol) {
      inds = as.character(df[1,i])
      column = df[,i]
      #print(inds)
      
      if (tabs==F) {
        split = strsplit(inds,",")[[1]]
      } else {
        split = strsplit(inds,"\t")[[1]]
      }
      
      numinds = length(split)
      
      
      
      if (numinds > 1) {
        new = column
        
        toremove = c(toremove,i)
        
        for (j in 2:numinds) {
          if (j == 2) {
            new = cbind(as.character(new),as.character(column))
            
          } else {
            new = cbind(new,as.character(column))
            
          }
          
          #print(new)
          
        }
        
        new[1,] = split
        
        df = cbind(df,new)
        
      }
      
      
      
    }
    
    #print(toremove)
    
    print("ROWS")
    df = t(df)
    endcol = ncol(df)
    for (i in 2:endcol) {
      inds = as.character(df[1,i])
      column = df[,i]
      #print(inds)
      
      if (tabs==F) {
        split = strsplit(inds,",")[[1]]
      } else {
        split = strsplit(inds,"\t")[[1]]
      }
      
      numinds = length(split)
      
      
      
      if (numinds > 1) {
        new = column
        
        toremove = c(toremove,i)
        
        for (j in 2:numinds) {
          if (j == 2) {
            new = cbind(as.character(new),as.character(column))
            
          } else {
            new = cbind(new,as.character(column))
            
          }
          
          #print(new)
          
        }
        
        new[1,] = split
        
        df = cbind(df,new)
        
      }
      
      
      
    }
    
    df = t(df)
    
    toremove = unique(toremove)
    
    # for (k in toremove) {
    #   #print(paste("Row",k,"of",endcol))
    #   
    #   inds = as.character(df[k,1])
    #   row = df[k,]
    #   #print(inds)
    #   
    #   if (tabs==F) {
    #     split = strsplit(inds,",")[[1]]
    #     
    #   } else {
    #     split = strsplit(inds,"\t")[[1]]
    #     
    #   }
    #   
    #   numinds = length(split)
    #   
    #   new = row
    #   
    #   for (j in 2:numinds) {
    #     new = rbind(new,row)
    #   }
    #   
    #   #print(new)
    #   
    #   
    #   
    #   new[,1] = split
    #   
    #   df = rbind(df,new)
    #   
    #   #print(new)
    #   
    #   
    #   
    # }
    
    print("SUBSETTING")
    
    if(!(is.null(toremove))) {
    
    df = subset(df, select = -toremove )
    df = t(df)
    df = subset(df, select = -toremove )
    df = t(df)
    }
    
    rownames(df) = df[,1]
    colnames(df) = df[1,]
    
    df = subset(df, select = -1 )
    df = t(df)
    df = subset(df, select = -1 )
    df = t(df)
    
    df=as.data.frame(df)
    
    print("WRITING")
    write.csv(df,file=paste(temp,".txt",sep=""))
    
  }
}
