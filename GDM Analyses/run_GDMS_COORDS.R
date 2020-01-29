rm(list=ls())

library(vegan)
library(sp)
library(sgdm)
library(rgdal)
library(raster)
library(parallel)
library(iterators)
library(gdm)
library(foreach)
library(ecodist)
library(doParallel)

print(version)

loadMorphData = function(spp,doGENE=F) {
  if (doGENE == F) {
    morph = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/measurements/aggregated datasets/KLP_Master_Spreadsheet_Morphology_AGGREGATED_19_July_2019.csv"
    morphdf = read.csv(morph)
    morphdf = morphdf[,2:40]
    #morphdata = morphdf[,c(3:9)]
  } else {
    morph = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/CATALOG_BY_SPECIES.txt"
    morphdf = read.csv(morph,sep = "\t")
    colnames(morphdf)[1] = ("CATALOG.NUMBER")
  }
  morphdf = unique(morphdf)
  morphdf = morphdf[morphdf$CATALOG.NUMBER != "",]
  morphlocs = morphdf[,c("LAT","LONG")]
  morphspp = morphdf$SPP
  return(list(morphdf,morphlocs,morphspp))
}

loadEnvData = function(morphlocs,morphspp,morphdf,doGENE=F,spp) {
  ## pure environments-- only worldclim i think
  envfile = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/enm_layers/ENMS_multilayer_onlyworldclim.tif"
  env = stack(envfile)
  envdata = as.data.frame(extract(env,morphlocs[,2:1]))
  envdata$SPP = morphspp
  envdata$CATALOG.NUMBER = morphdf$CATALOG.NUMBER
  #return(list(env,envdata))
  inputdata = merge(morphdf,envdata)
  data = unique(inputdata)
  data = data[data$CATALOG.NUMBER != "",]
  rownames(data) = data$CATALOG.NUMBER
  if (doGENE == F) {
    bioclim = data[,c(40:58)] 
  } else {
    bioclim = data[,c(5:23)] 
  }
  colnames(bioclim) = c("BIO1","BIO2","BIO3","BIO4","BIO5","BIO6","BIO7","BIO8","BIO9","BIO10","BIO11","BIO12","BIO13","BIO14","BIO15","BIO16","BIO17","BIO18","BIO19")
  bioclim = bioclim[complete.cases(bioclim),]
  if (doGENE == F) {
    m = cor(bioclim,use = "pairwise.complete.obs")
    xy <- t(combn(colnames(m),2))
    z = data.frame(xy,dist = m[xy])
    z[abs(z$dist) >= 0.7,]
    ## to remove because over 0.7: "BIO2","BIO3","BIO4","BIO8","BIO9","BIO10","BIO15","BIO16","BIO17"
    bioclim = bioclim[,c("BIO1","BIO5","BIO6","BIO7","BIO11","BIO12","BIO13","BIO14","BIO18","BIO19")]
    temp_pca <-
      prcomp(bioclim[,c("BIO1","BIO5","BIO6","BIO7","BIO11")],center = TRUE,scale. = TRUE)
    temp_scores <- predict(temp_pca)[,1:3]
    colnames(temp_scores) <- c("PC1T","PC2T","PC3T")
    prec_pca <-
      prcomp(bioclim[,c("BIO12","BIO13","BIO14","BIO18","BIO19")],center = TRUE,scale. = TRUE)
    summary(prec_pca)
    prec_scores <- predict(prec_pca)[,1:3]
    colnames(prec_scores) <- c("PC1P","PC2P","PC3P")
    toprint_temp = rbind(temp_pca$rotation,summary(temp_pca)$importance)
    toprint_prec = rbind(prec_pca$rotation,summary(prec_pca)$importance)
    write.csv(toprint_temp,paste("temperature_pca.csv",sep = ""))
    write.csv(toprint_prec,paste("precipitation_pca.csv",sep = ""))
  } else {
    temp_pca = read.csv(paste("temperature_pca.csv",sep = ""),row.names = 1)
    temp_pca = temp_pca[1:5,]
    prec_pca = read.csv(paste("precipitation_pca.csv",sep = ""),row.names = 1)
    prec_pca = prec_pca[1:5,]
    temp_scores = bioclim[,c(1:5)]
    prec_scores = bioclim[,c(6:10)]
    for (i in 1:nrow(prec_scores)) {
      row = t(prec_scores[i,])
      newrow = colSums(row * prec_pca)
      prec_scores[i,] = newrow
      row = t(temp_scores[i,])
      newrow = colSums(row * temp_pca)
      temp_scores[i,] = newrow
    }
    temp_scores <- temp_scores[,1:3]
    colnames(temp_scores) <- c("PC1T","PC2T","PC3T")
    prec_scores = prec_scores[,1:3]
    colnames(prec_scores) <- c("PC1P","PC2P","PC3P")
  }
  bioclim_final <-
    cbind(bioclim[,1:2],temp_scores,prec_scores)
  bioclim_final = bioclim_final[,3:ncol(bioclim_final)]
  for (x in 2:ncol(bioclim_final)) {
    bioclim_final[,x] = scales::rescale(bioclim_final[,x])
  }
  #return(bioclim_final)
  data = data[data$SPP == spp,]
  gene = data[,c("CATALOG.NUMBER","LAT","LONG")]
  colnames(gene) = c("Sample","Latitude","Longitude")
  rownames(gene) <- gene$Sample
  samples = gene$Sample
  return(list(bioclim_final,data,gene,samples))
}

generateMatrixCsvString = function(matVariable,doPCA=F,doGENE=F,doChroms=F,spp,whichpca,whichchrom) {
  if(matVariable=="STR") {
    csvstring=paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/BY_SPECIES/CSVS/",spp,"/",spp,"_distancematrix_STR.csv",sep = "")
  } else if (matVariable=="PRES") {
    csvstring=list.files(path = paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/BY_SPECIES/OTHERTEXT/",spp,sep = ""),
                         pattern = glob2rx("*worldclim.asc*AGGFACT-1.*MORPH-AND-GENE-DISTANCES.txt"),full.names = T)
  } else if (matVariable=="MID") {
    csvstring=list.files(path = paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/BY_SPECIES/OTHERTEXT/",spp,sep = ""),
                         pattern = glob2rx("*MID.asc*AGGFACT-1.*MORPH-AND-GENE-DISTANCES.txt"),
                         full.names = T)
  } else if (matVariable=="LGM") {
    csvstring=list.files(path = paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/BY_SPECIES/OTHERTEXT/",spp,sep = ""),
                         pattern = glob2rx("*LGM.asc*AGGFACT-1.*MORPH-AND-GENE-DISTANCES.txt"),full.names = T)
    # /Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/BY_SPECIES/OTHERTEXT/BELLII/AddedPredThin_BestModel_Vireo bellii_addedLayers_LGM.asc_-5.38855934143066_-61.0209884643555.rescaled.ascAGGFACT-1.MORPH-AND-GENE-DISTANCES.txt
  } else if (matVariable=="IBD") {
    csvstring=paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/BY_SPECIES/CSVS/",spp,"/",spp,"_distancematrix_geog.csv",sep = "")
  } else if (matVariable=="ENV") {
    csvstring=paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/BY_SPECIES/CSVS/",spp,"/",spp,"_distancematrix_env.csv",sep = "")
  } else if (matVariable=="ABUN") {
    csvstring=paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/BY_SPECIES/OTHERTEXT/",
                    spp,"/",spp,"_Idw_abundance_clipped.ascAGGFACT-1.MORPH-AND-GENE-DISTANCES.txt",sep = "")
    # /Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/BY_SPECIES/OTHERTEXT/BEL/BELLII_Idw_abundance_clipped.ascAGGFACT-1.MORPH-AND-GENE-DISTANCES.txt
  } else if (matVariable=="GENE") {
    if (doPCA == F) {
      if (doGENE == F) {
        csvstring = paste( "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/BY_SPECIES/CSVS/",
                           spp,"/",spp,"_distancematrix_morph.csv",sep = "")
      } else {
        if (doChroms==F) {
          csvstring = paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/BY_SPECIES/CSVS/",
                            spp,"/",spp,"_distancematrix_FULLGENOME.csv.converted",sep = "")
        } else {
          csvstring = list.files(path=paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/BY_SPECIES/CSVS/",spp,"/",sep=""),
                                 pattern=paste("Tgut_",whichchrom,".csv.converted",sep=""),full.names = T)
        }
      }
    } else {
      csvstring = paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/BY_SPECIES/CSVS/",
                        spp,"/",spp,"_distancematrix_pc",whichpca,"morph.csv",sep = "")
      
    }
  } else {
    csvstring=NULL
  }
  print(paste("LOADED FILE:",csvstring))
  return(csvstring)
}

readDataChangeHeaderGDM = function(csvstring,makeUnique=F) {
  df = read.csv(csvstring,header = F,sep = "\t",stringsAsFactors = F)
  #df = read.csv(paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/BY_SPECIES/CONVERTED/",spp,"/",spp,"_Idw_abundance_clipped.ascAGGFACT-1.MORPH-AND-GENE-DISTANCES.txt.converted.txt",sep = ""),header = F,sep = "\t",stringsAsFactors=F)
  if (ncol(df) == 1) {
    #df = read.csv(paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/BY_SPECIES/CONVERTED/",spp,"/",spp,"_Idw_abundance_clipped.ascAGGFACT-1.MORPH-AND-GENE-DISTANCES.txt.converted.txt",sep = ""),header = T,sep = ",")
    df = read.csv(csvstring,header = F,sep = ",",stringsAsFactors = F)
  }
  ## make sure only using the first one
  colnames(df) = NULL
  colnames(df) = as.character(df[1,])
  splits = strsplit(as.character(df[1,]),"/")
  firstsplits = sapply(splits,FUN=function(x){return(x[1])})
  firstsplits[is.na(firstsplits)] = ""
  colnames(df) = firstsplits
  df[1,] = firstsplits
  colnames(df)[1] = "Sample"
  df[,1] = factor(df[,1])
  df = df[-1,]
  
  colnames(df) = make.unique(colnames(df))
  rownames(df) = colnames(df)[2:ncol(df)]
  
  if(makeUnique==T) {
    df = as.data.frame(unique(df))
    df = as.data.frame(t(unique(df)))
    df = as.data.frame(t(unique(df)))
    df = as.data.frame(unique(df))
    df = df[which(rownames(df) %in% colnames(df)),]
    
  }

  #rownames(df) = make.unique(stringr::str_replace_all(rownames(df),"[^[:alnum:]]",""))
  #colnames(df) = make.unique(stringr::str_replace_all(colnames(df),"[^[:alnum:]]",""))
  #df$Sample = stringr::str_replace_all(df$Sample,"[^[:alnum:]]","")
  
  
  #df = unique(df)
  
  #rownames(df) = stringr::str_replace_all(rownames(df),"[A-IN-Z]","")
  #colnames(df)[2:ncol(df)] = stringr::str_replace_all(colnames(df)[2:ncol(df)],"[A-IN-Z]","")
  #df$Sample = stringr::str_replace_all(df$Sample,"[A-IN-Z]","")
  
  return(df)
}

generateFormatted = function(univariate=F,columnnumber=NULL,bioclim_final,distPredList,GENE,dogeo=F,IBD) {
  if (univariate==F) {
    formatted <- formatsitepair(bioData = GENE,bioFormat = 3,siteColumn = "Sample",XColumn = "Longitude",YColumn = "Latitude",
                                predData = bioclim_final,distPreds = distPredList)
    formatted$distance = scales::rescale(formatted$distance)
    selected_col = "FULL"
  } else {
    i=columnnumber
    if (i <= 6) {
      selected_col = colnames(bioclim_final)[i+1]
      print(paste("Selected Column:",selected_col))
      bioclim_subset = bioclim_final[,c("Sample",selected_col,"Latitude","Longitude")]
      if (dogeo==F) {
        formatted <- formatsitepair(bioData = GENE,bioFormat = 3,siteColumn = "Sample",XColumn = "Longitude",YColumn = "Latitude",
                                    predData = bioclim_subset)
      } else {
        formatted <- formatsitepair(bioData = GENE,bioFormat = 3,siteColumn = "Sample",XColumn = "Longitude",YColumn = "Latitude",
                                    predData = bioclim_subset,distPreds = list(as.matrix(IBD)))
      }
      formatted$distance = scales::rescale(formatted$distance)
    } else {
      #todolist = list(as.matrix(IBD),as.matrix(LGM),as.matrix(PRES),as.matrix(ENV),as.matrix(MID),as.matrix(STR),as.matrix(ABUN))
      #names(todolist) = c("IBD","LGM","PRES","ENV","MID","STR","ABUN")
      todolist = distPredList
      selected_col = names(todolist)[i-6]
      if(dogeo==T) {
        ibdnum = which(names(todolist)=="IBD")
        todolist = todolist[unique(c(ibdnum,i-6))]
      } else {
        todolist = todolist[unique(c(i-6))]
      }
      bioclim_subset = bioclim_final[,c("Sample","Latitude","Longitude")]
      formatted <- formatsitepair(bioData = GENE,bioFormat = 3,siteColumn = "Sample",XColumn = "Longitude",YColumn = "Latitude",
                                  predData = bioclim_subset,
                                  distPreds = todolist)
      formatted$distance = scales::rescale(formatted$distance)
    }
  }
  return(list(formatted,selected_col))
}

createModels = function(formatted,univariate,dogeo) {
  numvars=((ncol(formatted)-6)/2)
  print(numvars)
  if (numvars > 2) {
    geotag = F; fulltag="full"
  } else {
    geotag = dogeo;fulltag="selected column"
  }
  model <- gdm::gdm(formatted,geo = geotag,splines = NULL,knots = NULL)
  if(is.null(model)) {
    model <- gdm(formatted,geo = F,splines = rep(10,numvars),knots = NULL) 
    counter=-1
  } else {
    counter=-2 
  }
  while(is.null(model)) {
    counter = counter + 1
    model <- gdm(formatted,geo = F,knots = rep(counter,3*numvars))
    if(counter >= 100) {
      break 
    }
  }
  print(paste("Final Counter (-2 = default, -1 = 10 splines):",counter))
  print(paste("Explained by:",fulltag,model$explained)) ## need to save this
  print(paste("Null deviance: ",model$nulldeviance))
  print(paste("GDM deviance: ",model$gdmdeviance)) 
  print(paste("Model splines:",model$splines))
  #residuals = model$predicted-model$observed
  return(model)
}

doImportanceTesting = function(doImp,model,formatted,distPredList,doPCA,doGENE,spp) {
  start_time <- Sys.time()
  print(start_time)
  if(model$splines[1] == 3) {
    modTest <- gdm.varImp(formatted,geo = F,fullModelOnly = F,splines=model$splines)
  } else {
    modTest <- gdm.varImp(formatted,geo = F,fullModelOnly = T,splines=model$splines)
  }
  end_time <- Sys.time()
  print(end_time - start_time)
  rownames(modTest[[2]])[(length(model$predictors) - 6):length(model$predictors)] = distPredList
  rownames(modTest[[3]])[(length(model$predictors) - 6):length(model$predictors)] = distPredList
  rownames(modTest[[4]])[(length(model$predictors) - 6):length(model$predictors)] = distPredList
  if (doPCA == F) {
    if (doGENE == F) {
      suffix="_allmorph" 
    } else {
      
      if(doChroms==F) {
      suffix="_gene" } else {suffix=paste("_chr",whichchrom,sep="") }
    }
  } else {
    suffix=paste("_pc",whichpca,sep = "") 
  }
  #Deviance
  modTest[[1]]
  write.csv(modTest[[1]],paste(spp,"_variable_deviance",suffix,".csv",sep = ""))
  #Importance
  modTest[[2]]
  write.csv(modTest[[2]],paste(spp,"_variable_importance",suffix,".csv",sep = ""))
  #Significance
  modTest[[3]]
  write.csv(modTest[[3]],paste(spp,"_variable_significance",suffix,".csv",sep = ""))
  #Permuts Calculated
  modTest[[4]]
  write.csv(modTest[[4]],paste(spp,"_variable_permutations",suffix,".csv",sep = ""))
  png(paste(spp,"_importance",suffix,".png",sep = ""),height = 500,width = 1200)
  if(model$splines[1] == 3) {
    par(mfrow = c(3,4)) 
  }
  barplot(modTest[[2]][,1] / colSums(modTest[[2]],na.rm = T)[1],main = "full Model",las = 2)
  if(model$splines[1] == 3) {
    lapply(2:ncol(modTest[[2]]),FUN = function(x) {barplot(modTest[[2]][,x] / colSums(modTest[[2]],na.rm = T)[x], main = paste("full Model -",x - 1), las = 2 )})
  }
  dev.off()
}

runGDMs = function(doimp = F, doPCA = F,doGENE = F,whichpca = 1,univariate = T,dogeo = T,doChroms = F,whichchrom="Z",
                   matList=c("STR","PRES","MID","LGM","IBD","ENV","ABUN","GENE"),
                   spplist=c("BELLII","BILINEATA","BRUNNEICAPILLUS","CRISSALE","CURVIROSTRE","FLAVICEPS","FUSCA","MELANURA","NITENS","SINUATUS")) {
  print(paste("Do Importance?",doimp,"Do Gene?",doGENE,"Do Pca?",doPCA,"Which PCA?",whichpca))
  print(paste("Do univariate?",univariate,"Do IBD?",dogeo,"Do chromosomes?",doChroms,"Which chroms?",whichchrom))
  for (spp in spplist) {
    print(spp)
    print("loading morph")
    morphlist = loadMorphData(spp,doGENE)
    
    morphdf = morphlist[[1]]
    morphlocs = morphlist[[2]]
    morphspp = morphlist[[3]]
    setwd("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/")
    print("loading environment")
    envlist = loadEnvData(morphlocs,morphspp,morphdf,doGENE,spp)
    bioclim_final = envlist[[1]]
    data = envlist[[2]]
    gene = envlist[[3]]
    samples = envlist[[4]]
    
    gsplit = strsplit(as.character(gene$Sample),"/")
    gene$Sample = sapply(gsplit,FUN=function(x){return(x[1])})
    gene=unique(gene)
    
    print("relabeling")
    #colnames(bioclim_final) = (stringr::str_replace_all(colnames(bioclim_final),"[^[:alnum:]]",""))
    #bioclim_final$Sample = (stringr::str_replace_all(rownames(bioclim_final),"[^[:alnum:]]",""))
    #rownames(bioclim_final) = (stringr::str_replace_all(rownames(bioclim_final),"[^[:alnum:]]",""))
    
    bsplit = strsplit(as.character(bioclim_final$Sample),"/")
    #bioclim_final$Sample = sapply(bsplit,FUN=function(x){return(x[1])})
    #bioclim_final$Sample = stringr::str_replace_all(bioclim_final$Sample,"[^[:alnum:]]","")
    
    #bioclim_final=unique(bioclim_final)
    
    #rownames(gene) = stringr::str_replace_all(rownames(gene),"[^[:alnum:]]","")
    #colnames(gene) = stringr::str_replace_all(colnames(gene),"[^[:alnum:]]","")
    #gene$Sample = stringr::str_replace_all(gene$Sample,"[^[:alnum:]]","")
    bioclim_final = merge(bioclim_final,gene)
    
    rownames(bioclim_final) = make.unique(bioclim_final$Sample)
    #bioclim_final$Sample = stringr::str_replace_all(bioclim_final$Sample,"[^[:alnum:]]","")
    
    #gene$Sample = stringr::str_replace_all(gene$Sample,"[A-IN-Z]","")
    #rownames(gene) = stringr::str_replace_all(rownames(gene),"[A-IN-Z]","")
    #bioclim_final$Sample = stringr::str_replace_all(bioclim_final$Sample,"[A-IN-Z]","")
    #rownames(bioclim_final) = stringr::str_replace_all(rownames(bioclim_final),"[A-IN-Z]","")

    ## load matList 
    print("loading matlist")
    toexit=F
    for (matVariable in matList) {
      print(matVariable)
      csvstring = generateMatrixCsvString(matVariable,doPCA,doGENE,doChroms,spp,whichpca,whichchrom)
      if(length(csvstring) > 1) {
        df = lapply(csvstring,FUN=readDataChangeHeaderGDM)
        #numgenes = length(csvstring)
        assign(matVariable,df)
      } else if (length(csvstring) == 1) {
        
        if(matVariable=="STR") {
          df = readDataChangeHeaderGDM(csvstring,makeUnique = F)
        } else {
          df = readDataChangeHeaderGDM(csvstring,makeUnique = T)
        }
        
        
        assign(matVariable,df)
      } else {
        print("FILE DOES NOT EXIST, SKIPPING")
        toexit = T
      }
      
    }
    
    if(toexit==T) {
      next
    }
    
    print("cross-referencing matlists")
    # generate the ones to keep
    tokeep = c()
    for (matVariable in matList) {
      print(matVariable)
      print("Samples here:")
      print(sort(get(matVariable)$Sample))
      if (is.null(tokeep)) { tokeep = intersect(get(matVariable)$Sample,get(matVariable)$Sample) } else { tokeep = intersect(tokeep,get(matVariable)$Sample) }
      print(paste("temp keep:",length(tokeep)))
      print(sort(tokeep))
    }
    tokeep = intersect(tokeep,bioclim_final$Sample)
    tokeep = c(tokeep,"Sample")
    print(paste("Keeping:",length(tokeep)))
    ## remove the ones that don't match 
    for (matVariable in matList) {
      #print(matVariable)
      toset = get(matVariable)[which(get(matVariable)$Sample %in% tokeep),
                                 which(colnames(get(matVariable)) %in% tokeep)]
      assign(matVariable,toset)
    }
    
    print("done removing matches")
        bioclim_final = bioclim_final[which(bioclim_final$Sample %in% tokeep),]
    distPredList = lapply(matList,FUN = function(x) {df = get(matVariable); df = as.matrix(df); return(df) })
    names(distPredList) = matList
    distPredList = distPredList[names(distPredList)!="GENE"]
    if(univariate==F) {
      formattedlist = generateFormatted(univariate,columnnumber=NULL,bioclim_final,distPredList,GENE,dogeo,IBD)
      formatted = formattedlist[[1]]
      selected_col = formattedlist[[2]]
      ## run models 
      model = createModels(formatted,univariate,dogeo)
      model$predictors[7:length(model$predictors)] = names(distPredList)
      
      if (doPCA == F) {
        if (doGENE == F) {
          suffix="_allmorph" 
        } else {
          
          if(doChroms==F) {
            suffix="_gene" } else {suffix=paste("_chr",whichchrom,sep="") }
        }
      } else {
        suffix=paste("_pc",whichpca,sep = "") 
      }
      
      splinestring = paste(spp,"_splines",suffix,".png",sep = "") 
      png(splinestring,height = 900,width = 1500) 
      plot(model,plot.layout = c(2,6))
      dev.off()
      if(doimp==T) {
        doImportanceTesting(doImp,model,formatted,distPredList,doPCA,doGENE,spp)
      }
    } else {
      for (columnnumber in 1:(6+length(distPredList))) {
        #names(bioclim_final)[columnnumber]
        formattedlist = generateFormatted(univariate,columnnumber,bioclim_final,distPredList,GENE,dogeo,IBD)
        formatted = formattedlist[[1]]
        selected_col = formattedlist[[2]]
        model = createModels(formatted,univariate,dogeo)
        if (!(is.null(model))) {
          length(model$predictors)
          # model$predictors = selected_col
          if(dogeo==T) {geotext="IBD"} else {geotext=""}
          
          pngtext = paste(substr(spp,1,3),"_",selected_col,"_splines",geotext,suffix,".png",sep = "")
          png(pngtext,height = 450,width = 600) 
          plot(model,plot.layout = c(1,4))
          mtext(text=paste(spp,selected_col,"Expl:",(model$explained)))
          dev.off()
        }
      }
    }
  }
}


chromlist=c("1","1A","1B","2","3","4","4A","5","6","7","8","9","10",
            "11","12","13","14","15","16","17","18","19","20","21",
            "22","23","24","25","26","27","28","LG2","LGE22","mtDNA","Z")
spplist=c("BELLII","BILINEATA","BRUNNEICAPILLUS","CRISSALE","CURVIROSTRE","FLAVICEPS","FUSCA","MELANURA","NITENS","SINUATUS")
matList=c("GENE","STR","IBD","ENV","ABUN","PRES","MID","LGM")

print("full genes without importance")
runGDMs(doimp = F,doPCA = F,doGENE = T,whichpca = 1,univariate = F,dogeo = F,doChroms = F,whichchrom="Z",matList,spplist="BELLII")
  
print("chroms without importance")
for (chrom in chromlist) {
  print(chrom)
  runGDMs(doimp = F,doPCA = F,doGENE = T,whichpca = 1,univariate = F,dogeo = F,doChroms = T,whichchrom=chrom,matList,spplist)
}
print("chroms with importance")
for (chrom in chromlist) {
  print(chrom)
  runGDMs(doimp = T,doPCA = F,doGENE = T,whichpca = 1,univariate = F,dogeo = F,doChroms = T,whichchrom=chrom,matList,spplist)
}


# 
# doimp = F
# doPCA = F
# doGENE = F
# whichpca = 1
# univariate = T
# dogeo = T
# doChroms = F
# options = c("FFF1","FTF1","FTF2","FTF3",
#             "FFT1"#,
#             #"TFF1","TTF1","TTF2","TTF3",
#             #"TFT1"
# )
# for (opt in options) 
# {
#   print(opt)
#   doimp = as.logical(substr(opt,1,1))
#   doPCA = as.logical(substr(opt,2,2))
#   doGENE = as.logical(substr(opt,3,3))
#   whichpca = as.numeric(substr(opt,4,4))
#   print(paste("Do Importance?",doimp,"Do Gene?",doGENE,"Do Pca?",doPCA,"Which PCA?",whichpca))
#   for (spp in (c("BELLII",## weird gene
#                  "BRUNNEICAPILLUS",
#                  "CRISSALE",## weird pca1
#                  "CURVIROSTRE",####
#                  "FLAVICEPS",####### -- problem with GENE,weird data?
#                  "FUSCA",
#                  "MELANURA",## weirc pca1
#                  "NITENS",## weird gene
#                  "SINUATUS",## weird with pca 1
#                  "BILINEATA" ##### -- problem with GENE,no distance matrix
#   ))) 
#   {
#     print(spp)
#     if (doGENE == F) 
#     {
#       morph = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/morphology/measurements/aggregated datasets/KLP_Master_Spreadsheet_Morphology_AGGREGATED_19_July_2019.csv"
#       morphdf = read.csv(morph)
#       morphdf = morphdf[,2:40]
#       morphdf = unique(morphdf)
#       morphdf = morphdf[morphdf$CATALOG.NUMBER != "",]
#       morphlocs = morphdf[,c("LAT","LONG")]
#       morphdata = morphdf[,c(3:9)]
#       morphspp = morphdf$SPP
#       
#     }
#     else 
#     {
#       morph = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/CATALOG_BY_SPECIES.txt"
#       morphdf = read.csv(morph,sep = "\t")
#       morphdf = unique(morphdf)
#       colnames(morphdf)[1] = ("CATALOG.NUMBER")
#       morphdf = morphdf[morphdf$CATALOG.NUMBER != "",]
#       morphlocs = morphdf[,c("LAT","LONG")]
#       morphspp = morphdf$SPP
#       
#     }
#     setwd("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/")
#     ## pure environments-- only worldclim i think
#     envfile = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/enm_layers/ENMS_multilayer_onlyworldclim.tif"
#     env = stack(envfile)
#     envdata = as.data.frame(extract(env,morphlocs[,2:1]))
#     envdata$SPP = morphspp
#     envdata$CATALOG.NUMBER = morphdf$CATALOG.NUMBER
#     inputdata = merge(morphdf,envdata)
#     data = unique(inputdata)
#     data = data[data$CATALOG.NUMBER != "",]
#     rownames(data) = data$CATALOG.NUMBER
#     if (doGENE == F) 
#     {
#       bioclim = data[,c(40:58)] 
#     }
#     else 
#     {
#       bioclim = data[,c(5:23)] 
#     }
#     colnames(bioclim) = c("BIO1","BIO2","BIO3","BIO4","BIO5","BIO6","BIO7","BIO8","BIO9","BIO10","BIO11","BIO12","BIO13","BIO14","BIO15","BIO16","BIO17","BIO18","BIO19")
#     bioclim = bioclim[complete.cases(bioclim),]
#     if (doGENE == F) 
#     {
#       m = cor(bioclim,use = "pairwise.complete.obs")
#       xy <- t(combn(colnames(m),2))
#       z = data.frame(xy,dist = m[xy])
#       z[abs(z$dist) >= 0.7,]
#       ## to remove because over 0.7: "BIO2","BIO3","BIO4","BIO8","BIO9","BIO10","BIO15","BIO16","BIO17"
#       bioclim = bioclim[,c("BIO1","BIO5","BIO6","BIO7","BIO11","BIO12","BIO13","BIO14","BIO18","BIO19")]
#       temp_pca <-
#         prcomp(bioclim[,c("BIO1","BIO5","BIO6","BIO7","BIO11")],center = TRUE,scale. = TRUE)
#       temp_scores <- predict(temp_pca)[,1:3]
#       colnames(temp_scores) <- c("PC1T","PC2T","PC3T")
#       prec_pca <-
#         prcomp(bioclim[,c("BIO12","BIO13","BIO14","BIO18","BIO19")],center = TRUE,scale. = TRUE)
#       summary(prec_pca)
#       prec_scores <- predict(prec_pca)[,1:3]
#       colnames(prec_scores) <- c("PC1P","PC2P","PC3P")
#       toprint_temp = rbind(temp_pca$rotation,summary(temp_pca)$importance)
#       toprint_prec = rbind(prec_pca$rotation,summary(prec_pca)$importance)
#       write.csv(toprint_temp,paste("temperature_pca.csv",sep = ""))
#       write.csv(toprint_prec,paste("precipitation_pca.csv",sep = ""))
#       
#     }
#     else 
#     {
#       temp_pca = read.csv(paste("temperature_pca.csv",sep = ""),row.names = 1)
#       temp_pca = temp_pca[1:5,]
#       prec_pca = read.csv(paste("precipitation_pca.csv",sep = ""),row.names = 1)
#       prec_pca = prec_pca[1:5,]
#       temp_scores = bioclim[,c(1:5)]
#       prec_scores = bioclim[,c(6:10)]
#       for (i in 1:nrow(prec_scores)) 
#       {
#         row = t(prec_scores[i,])
#         newrow = colSums(row * prec_pca)
#         prec_scores[i,] = newrow
#         row = t(temp_scores[i,])
#         newrow = colSums(row * temp_pca)
#         temp_scores[i,] = newrow
#         
#       }
#       temp_scores <- temp_scores[,1:3]
#       colnames(temp_scores) <- c("PC1T","PC2T","PC3T")
#       prec_scores = prec_scores[,1:3]
#       colnames(prec_scores) <- c("PC1P","PC2P","PC3P")
#       
#     }
#     bioclim_final <-
#       cbind(bioclim[,1:2],temp_scores,prec_scores)
#     bioclim_final = bioclim_final[,3:ncol(bioclim_final)]
#     for (x in 2:ncol(bioclim_final)) 
#     {
#       bioclim_final[,x] = scales::rescale(bioclim_final[,x])
#       
#     }
#     data = data[data$SPP == spp,]
#     gene = data[,c("CATALOG.NUMBER","LAT","LONG")]
#     colnames(gene) = c("Sample","Latitude","Longitude")
#     rownames(gene) <- gene$Sample
#     samples = gene$Sample
#     if (doPCA == F) 
#     {
#       if (doGENE == F) 
#       {
#         GENE = read.csv(paste( "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/BY_SPECIES/CSVS/",spp,"/",spp,"_distancematrix_morph.csv",sep = ""),header = T)
#         
#       }
#       else 
#       {
#         GENE = read.csv(paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/BY_SPECIES/CSVS/",spp,"/",spp,"_distancematrix_FULLGENOME.csv.converted",sep = ""),header = T)
#         
#       }
#       
#     }
#     else 
#     {
#       if(doChroms==F) 
#       {
#         GENE = read.csv(paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/BY_SPECIES/CSVS/",spp,"/",spp,"_distancematrix_pc",whichpca,"morph.csv",sep = ""),header = T)
#         
#       }
#       else 
#       {
#         print("FAIL")
#         
#       }
#       
#     }
#     
#     
#     # Predictor dissimilarity matrices
#     STR <- read.csv(paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/BY_SPECIES/CSVS/",spp,"/",spp,"_distancematrix_STR.csv",sep = ""))#[-remove,-remove]
#     PRES = read.csv(list.files(path = paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/BY_SPECIES/CONVERTED/",spp,sep = ""),pattern = glob2rx("*worldclim.asc*AGGFACT*MORPH-AND-GENE-DISTANCES.txt.converted.txt"),full.names = T),header = T)
#     MID = read.csv(list.files(path = paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/BY_SPECIES/CONVERTED/",spp,sep = ""),
#                               pattern = glob2rx("*MID.asc*AGGFACT*MORPH-AND-GENE-DISTANCES.txt.converted.txt"),
#                               full.names = T),
#                    header = T,
#                    row.names = NULL)
#     LGM = read.csv(list.files(path = paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/BY_SPECIES/CONVERTED/",spp,sep = ""),pattern = glob2rx("*LGM.asc*AGGFACT*MORPH-AND-GENE-DISTANCES.txt.converted.txt"),full.names = T),header = T)
#     IBD = read.csv(paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/BY_SPECIES/CSVS/",spp,"/",spp,"_distancematrix_geog.csv",sep = ""),header = T)
#     ENV <- read.csv(paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/BY_SPECIES/CSVS/",spp,"/",spp,"_distancematrix_env.csv",sep = ""))#[-remove,-remove]
#     ABUN = read.csv(paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/BY_SPECIES/CONVERTED/",spp,"/",spp,"_Idw_abundance_clipped.ascAGGFACT-1.MORPH-AND-GENE-DISTANCES.txt.converted.txt",sep = ""),header = T,sep = "\t")
#     if (ncol(ABUN) == 1) 
#     {
#       ABUN = read.csv(paste("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/BY_SPECIES/CONVERTED/",spp,"/",spp,"_Idw_abundance_clipped.ascAGGFACT-1.MORPH-AND-GENE-DISTANCES.txt.converted.txt",sep = ""),header = T,sep = ",")
#       
#     }
#     colnames(STR)[1] = "Sample"
#     colnames(PRES)[1] = "Sample"
#     colnames(MID)[1] = "Sample"
#     colnames(LGM)[1] = "Sample"
#     colnames(IBD)[1] = "Sample"
#     colnames(ENV)[1] = "Sample"
#     colnames(ABUN)[1] = "Sample"
#     colnames(GENE)[1] = "Sample"
#     # GDM
#     rownames(STR) = stringr::str_replace_all(rownames(STR),"[^[:alnum:]]","")
#     rownames(PRES) = stringr::str_replace_all(rownames(PRES),"[^[:alnum:]]","")
#     rownames(MID) = stringr::str_replace_all(rownames(MID),"[^[:alnum:]]","")
#     rownames(LGM) = stringr::str_replace_all(rownames(LGM),"[^[:alnum:]]","")
#     rownames(IBD) = stringr::str_replace_all(rownames(IBD),"[^[:alnum:]]","")
#     rownames(GENE) = stringr::str_replace_all(rownames(GENE),"[^[:alnum:]]","")
#     rownames(gene) = stringr::str_replace_all(rownames(gene),"[^[:alnum:]]","")
#     rownames(ENV) = stringr::str_replace_all(rownames(ENV),"[^[:alnum:]]","")
#     rownames(bioclim_final) = stringr::str_replace_all(rownames(bioclim_final),"[^[:alnum:]]","")
#     rownames(ABUN) = stringr::str_replace_all(rownames(ABUN),"[^[:alnum:]]","")
#     colnames(STR) = stringr::str_replace_all(colnames(STR),"[^[:alnum:]]","")
#     colnames(PRES) = stringr::str_replace_all(colnames(PRES),"[^[:alnum:]]","")
#     colnames(MID) = stringr::str_replace_all(colnames(MID),"[^[:alnum:]]","")
#     colnames(LGM) = stringr::str_replace_all(colnames(LGM),"[^[:alnum:]]","")
#     colnames(IBD) = stringr::str_replace_all(colnames(IBD),"[^[:alnum:]]","")
#     colnames(GENE) = stringr::str_replace_all(colnames(GENE),"[^[:alnum:]]","")
#     colnames(gene) = stringr::str_replace_all(colnames(gene),"[^[:alnum:]]","")
#     colnames(ENV) = stringr::str_replace_all(colnames(ENV),"[^[:alnum:]]","")
#     colnames(bioclim_final) = stringr::str_replace_all(colnames(bioclim_final),"[^[:alnum:]]","")
#     colnames(ABUN) = stringr::str_replace_all(colnames(ABUN),"[^[:alnum:]]","")
#     bioclim_final$Sample = rownames(bioclim_final)
#     STR$Sample = stringr::str_replace_all(STR$Sample,"[^[:alnum:]]","")
#     PRES$Sample = stringr::str_replace_all(PRES$Sample,"[^[:alnum:]]","")
#     MID$Sample = stringr::str_replace_all(MID$Sample,"[^[:alnum:]]","")
#     LGM$Sample = stringr::str_replace_all(LGM$Sample,"[^[:alnum:]]","")
#     IBD$Sample = stringr::str_replace_all(IBD$Sample,"[^[:alnum:]]","")
#     GENE$Sample = stringr::str_replace_all(GENE$Sample,"[^[:alnum:]]","")
#     gene$Sample = stringr::str_replace_all(gene$Sample,"[^[:alnum:]]","")
#     ENV$Sample = stringr::str_replace_all(ENV$Sample,"[^[:alnum:]]","")
#     bioclim_final$Sample = stringr::str_replace_all(bioclim_final$Sample,"[^[:alnum:]]","")
#     ABUN$Sample = stringr::str_replace_all(ABUN$Sample,"[^[:alnum:]]","")
#     bioclim_final = merge(bioclim_final,gene)
#     rownames(bioclim_final) = bioclim_final$Sample
#     bioclim_final$Sample = stringr::str_replace_all(bioclim_final$Sample,"[^[:alnum:]]","")
#     ## the ones to keep
#     tokeep = intersect(GENE$Sample,GENE$Sample)
#     tokeep = intersect(tokeep,STR$Sample) ## good
#     tokeep = intersect(tokeep,PRES$Sample) ## good
#     tokeep = intersect(tokeep,MID$Sample) ## not good
#     tokeep = intersect(tokeep,LGM$Sample)
#     tokeep = intersect(tokeep,IBD$Sample) 
#     tokeep = intersect(tokeep,ENV$Sample)
#     tokeep = intersect(tokeep,bioclim_final$Sample)
#     tokeep = intersect(tokeep,ABUN$Sample)
#     tokeep = c(tokeep,"Sample")
#     print(paste("Keeping:",length(tokeep)))
#     ## then we need to make sure that we trim the larger matrices
#     STR = STR[which(STR$Sample %in% tokeep),which(colnames(STR) %in% tokeep)]
#     PRES = PRES[which(PRES$Sample %in% tokeep),which(colnames(PRES) %in% tokeep)]
#     MID = MID[which(MID$Sample %in% tokeep),which(colnames(MID) %in% tokeep)]
#     LGM = LGM[which(LGM$Sample %in% tokeep),which(colnames(LGM) %in% tokeep)]
#     IBD = IBD[which(IBD$Sample %in% tokeep),which(colnames(IBD) %in% tokeep)]
#     GENE = GENE[which(GENE$Sample %in% tokeep),which(colnames(GENE) %in% tokeep)]
#     ENV = ENV[which(ENV$Sample %in% tokeep),which(colnames(ENV) %in% tokeep)]
#     bioclim_final = bioclim_final[which(bioclim_final$Sample %in% tokeep),]
#     ABUN = ABUN[which(ABUN$Sample %in% tokeep),which(colnames(ABUN) %in% tokeep)]
#     
#     
#     
#     if (univariate==F) 
#     {
#       
#       formatted <- formatsitepair(bioData = GENE,bioFormat = 3,siteColumn = "Sample",XColumn = "Longitude",YColumn = "Latitude",
#                                   predData = bioclim_final,distPreds = list(as.matrix(IBD),as.matrix(LGM),as.matrix(PRES),as.matrix(ENV),as.matrix(MID),as.matrix(STR),as.matrix(ABUN)))
#       formatted$distance = scales::rescale(formatted$distance)
#       ## create the GDM model
#       model <- gdm::gdm(formatted,geo = F,splines = NULL,knots = NULL)
#       model$nulldeviance
#       model$gdmdeviance
#       if(is.null(model)) 
#       {
#         
#         model <- gdm(formatted,geo = F,splines = rep(10,13),knots = NULL) 
#         counter=-1
#         
#       }
#       else 
#       {
#         counter=-2 
#       }
#       while(is.null(model)) 
#       {
#         counter = counter + 1
#         model <- gdm(formatted,geo = F,knots = rep(counter,3*13))
#         if(counter >= 100) 
#         {
#           break 
#         }
#         
#       }
#       print(paste("Final Counter (-2 = default, -1 = 10 splines):",counter))
#       print(paste("Explained by full:",model$explained)) ## need to save this
#       print(paste("Null deviance:",model$nulldeviance)) ## need to save this
#       print(paste("Full deviance (lower is better):",model$gdmdeviance)) ## need to save this
#       length(model$predictors)
#       model$predictors[(length(model$predictors) - 6):length(model$predictors)] = c("IBD","LGM","PRES","ENV","MID","STR","ABUN")
#       if (doPCA == F) 
#       {
#         if (doGENE == F) 
#         {
#           png(paste(spp,"_splines.png",sep = ""),height = 900,width = 1500) 
#         }
#         
#         else 
#         {
#           png(paste(spp,"_splines_gene.png",sep = ""),height = 900,width = 1500) 
#         }
#         
#       }
#       
#       else 
#       {
#         png(paste(spp,"_splines_pc",whichpca,".png",sep = ""),height = 900,width = 1500) 
#       }
#       plot(model,plot.layout = c(2,6))
#       dev.off()
#       
#       residuals = model$predicted-model$observed
#       
#       ## model importance tests need to be done for GENE with only the first model,not all the models.
#       if (doimp == T) 
#       {
#         
#         start_time <- Sys.time()
#         print(start_time)
#         if (counter == -2) 
#         {
#           modTest <- gdm.varImp(formatted,geo = F,fullModelOnly = F)
#           
#         }
#         else if (counter == -1) 
#         {
#           modTest <- gdm.varImp(formatted,geo = F,fullModelOnly = T,splines = rep(10,13))
#           
#         }
#         else 
#         {
#           modTest <- gdm.varImp(formatted,geo = F,fullModelOnly = T,knots = rep(counter,3*13))
#           
#         }
#         end_time <- Sys.time()
#         print(end_time - start_time)
#         
#         rownames(modTest[[2]])[(length(model$predictors) - 6):length(model$predictors)] = c("IBD","LGM","PRES","ENV","MID","STR","ABUN")
#         rownames(modTest[[3]])[(length(model$predictors) - 6):length(model$predictors)] = c("IBD","LGM","PRES","ENV","MID","STR","ABUN")
#         rownames(modTest[[4]])[(length(model$predictors) - 6):length(model$predictors)] = c("IBD","LGM","PRES","ENV","MID","STR","ABUN")
#         
#         if (doPCA == F) 
#         {
#           if (doGENE == F) 
#           {
#             #Deviance
#             modTest[[1]]
#             write.csv(modTest[[1]],paste(spp,"_variable_deviance_test.csv",sep = ""))
#             #Importance
#             modTest[[2]]
#             write.csv(modTest[[2]],paste(spp,"_variable_importance_test.csv",sep = ""))
#             #Significance
#             modTest[[3]]
#             write.csv(modTest[[3]],paste(spp,"_variable_significance_test.csv",sep = ""))
#             #Permuts Calculated
#             modTest[[4]]
#             write.csv(modTest[[4]],paste(spp,"_variable_permutations_test.csv",sep = ""))
#             
#           }
#           
#           else 
#           {
#             #Deviance
#             modTest[[1]]
#             write.csv(modTest[[1]],paste(spp,"_variable_deviance_gene.csv",sep = ""))
#             #Importance
#             modTest[[2]]
#             write.csv(modTest[[2]],paste(spp,"_variable_importance_gene.csv",sep = ""))
#             #Significance
#             modTest[[3]]
#             write.csv(modTest[[3]],paste(spp,"_variable_significance_gene.csv",sep = ""))
#             #Permuts Calculated
#             modTest[[4]]
#             write.csv(modTest[[4]],paste(spp,"_variable_permutations_gene.csv",sep = ""))
#             
#           }
#           
#         }
#         
#         else 
#         {
#           #Deviance
#           modTest[[1]]
#           write.csv(modTest[[1]],paste(spp,"_variable_deviance_pc",whichpca,".csv",sep = ""))
#           #Importance
#           modTest[[2]]
#           write.csv(modTest[[2]],paste(spp,"_variable_importance_pc",whichpca,".csv",sep = ""))
#           #Significance
#           modTest[[3]]
#           write.csv(modTest[[3]],paste(spp,"_variable_significance_pc",whichpca,".csv",sep = ""))
#           #Permuts Calculated
#           modTest[[4]]
#           write.csv(modTest[[4]],paste(spp,"_variable_permutations_pc",whichpca,".csv",sep = ""))
#           
#         }
#         if (doPCA == F) 
#         {
#           if (doGENE == F) 
#           {
#             png(paste(spp,"_importance.png",sep = ""),height = 500,width = 1200)
#           }
#           
#           else 
#           {
#             png(paste(spp,"_importance_gene.png",sep = ""),height = 500,width = 1200 )
#           }
#           
#         }
#         else 
#         {
#           png(paste(spp,"_importance_pc",whichpca,".png",sep = ""),height = 500,width = 1200) 
#         }
#         if(counter == -2) 
#         {
#           par(mfrow = c(3,4)) 
#         }
#         barplot(modTest[[2]][,1] / colSums(modTest[[2]],na.rm = T)[1],main = "full Model",las = 2)
#         if(counter == -2) 
#         {
#           lapply(2:ncol(modTest[[2]]),FUN = function(x) 
#           {
#             barplot(modTest[[2]][,x] / colSums(modTest[[2]],na.rm = T)[x], main = paste("full Model -",x - 1), las = 2 )
#           }
#           )
#         }
#         dev.off()
#         
#       }
#       
#     }
#     
#     else if (univariate == T) 
#     {
#       
#       for(i in 1:13) 
#       {
#         if (i <= 6) 
#         {
#           selected_col = colnames(bioclim_final)[i+1]
#           bioclim_subset = bioclim_final[,c("Sample",selected_col,"Latitude","Longitude")]
#           
#           if (dogeo==F) 
#           {
#             formatted <- formatsitepair(bioData = GENE,bioFormat = 3,siteColumn = "Sample",XColumn = "Longitude",YColumn = "Latitude",
#                                         predData = bioclim_subset)
#             
#           }
#           else 
#           {
#             formatted <- formatsitepair(bioData = GENE,bioFormat = 3,siteColumn = "Sample",XColumn = "Longitude",YColumn = "Latitude",
#                                         predData = bioclim_subset,distPreds = list(as.matrix(IBD)))
#             
#           }
#           formatted$distance = scales::rescale(formatted$distance)
#           
#           
#         }
#         else 
#         {
#           
#           todolist = list(as.matrix(IBD),as.matrix(LGM),as.matrix(PRES),as.matrix(ENV),as.matrix(MID),as.matrix(STR),as.matrix(ABUN))
#           names(todolist) = c("IBD","LGM","PRES","ENV","MID","STR","ABUN")
#           
#           selected_col = names(todolist)[i-6]
#           
#           if(dogeo==T) 
#           {
#             
#             todolist = todolist[unique(c(1,i-6))]
#             
#             
#           }
#           else 
#           {
#             
#             todolist = todolist[unique(c(i-6))]
#             
#             
#           }
#           
#           bioclim_subset = bioclim_final[,c("Sample","Latitude","Longitude")]
#           formatted <- formatsitepair(bioData = GENE,bioFormat = 3,siteColumn = "Sample",XColumn = "Longitude",YColumn = "Latitude",
#                                       predData = bioclim_subset,
#                                       distPreds = todolist)
#           
#           formatted$distance = scales::rescale(formatted$distance)
#           
#         }
#         
#         if(dogeo==F) 
#         {
#           numvars = 1
#         }
#         else 
#         {
#           if(i!= 7) 
#           {
#             numvars=2
#           }
#           else 
#           {
#             numvars=1
#           }
#           
#           model <- gdm(formatted,geo = F,splines = NULL,knots = NULL)
#           if(is.null(model)) 
#           {
#             
#             model <- gdm(formatted,geo = F,splines = rep(10,numvars),knots = NULL) 
#             counter=-1
#             
#           }
#           else 
#           {
#             counter=-2 
#           }
#           while(is.null(model)) 
#           {
#             counter = counter + 1
#             model <- gdm(formatted,geo = F,knots = rep(counter,3*numvars))
#             if(counter >= 100) 
#             {
#               break 
#             }
#             
#           }
#           print(paste("Final Counter (-2 = default, -1 = 10 splines):",counter))
#           print(paste("Explained by",selected_col,":",(model$explained))) ## need to save this
#           
#           if(dogeo==T) 
#           {
#             geotext="IBD"
#           }
#           else 
#           {
#             geotext=""
#           }
#           
#           if (!(is.null(model))) 
#           {
#             
#             length(model$predictors)
#             model$predictors = selected_col
#             if (doPCA == F) 
#             {
#               if (doGENE == F) 
#               {
#                 png(paste(substr(spp,1,3),"_",selected_col,"_splines",geotext,".png",sep = ""),height = 450,width = 600) 
#               }
#               
#               else 
#               {
#                 png(paste(substr(spp,1,3),"_",selected_col,"_splines",geotext,"_gene.png",sep = ""),height = 450,width = 600) 
#               }
#               
#             }
#             
#             else 
#             {
#               png(paste(substr(spp,1,3),"_",selected_col,"_splines",geotext,"_pc",whichpca,".png",sep = ""),height = 450,width = 600) 
#             }
#             plot(model,plot.layout = c(1,4))
#             mtext(text=paste(spp,selected_col,"Expl:",(model$explained)))
#             dev.off()
#             
#           }
#           else 
#           {
#             print("IS NULL CANNOT RUN")
#           }
#           
#         }
#         
#       }
#       
#     }
#     
#   }
# }
