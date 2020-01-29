#install.packages("gdm")
library(gdm)
library(parallel)
library(iterators)
library(foParallel)
library(foreach)
library(sp)
library(raster)

## just gonna try with bellii for now
## using lucas' names to make sure i understand how works -- will need to update

#data <- read.csv("Icterus_ingroup_environmentaldata.csv",header=T)
data = unique(inputdata)

#rownames(data) <- data$Sample
rownames(data) = data$CATALOG.NUMBER
data = data[data$SPP=="BELLII",]

# Genetic data (response variable)

#gene <- data[,c("Sample","Latitude","Longitude")]
gene = data[,c("CATALOG.NUMBER","LAT","LONG")]
colnames(gene) = c("Sample","Latitude","Longitude")

rownames(gene) <- gene$Sample
#samples = c("A1","A10","A3","A4","A5","A6","A8","A9","B10","B4","B5","B7","B9","C1","C10","C2","C3","C4","C5","C6","C7","C9","D1","D10","D3","D4","D6","D7","D8","D9","E1","E10","E2","E3","E4","E5","E6","E7","E8","E9","F1","F2","F3","F5","F6","F7","F8","F9","G1","G2","G3","G4","G6","G7","G8","G9","H1","H2","H3","H4","H5","H6","H8","H9")
samples = gene$Sample
#gene <- gene[samples,]

#remove = c(1,3,7,14,16,17,23,25,31,33,34,41,43,45,49,50,52,53,57,58,62)
#remove = which(is.na(morphdf$BILL.HEIGHT))
#gene <- gene[-remove,]



#GENE <- read.csv("Icterus_ingroup_IBS_dissimilarity.csv",header=T)[-remove,-remove]
GENE = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/BELLII_distancematrix_morph.csv",
                header=T)#[-remove,-remove]
colnames(GENE)[1] = "Sample"
#GENE <- cbind(gene$Sample,GENE)
#colnames(GENE) <- c("Sample","A10","A4","A5","A6","A9","B10","B4","B5","B7","B9","C10","C4","C5","C6","C7","C9","D10","D4","D6","D7","D8","D9","E10","E4","E5","E6","E7","E8","E9","F2","F5","F7","F8","F9","G3","G7","G8","G9","H3","H4","H5","H8","H9")

# Environmental data

#bioclim <- data[samples,c(1,3,4,9:28)][-remove,]
bioclim = data[,c(40:58)]
colnames(bioclim) = c("BIO1","BIO2","BIO3","BIO4","BIO5",
                      "BIO6","BIO7","BIO8","BIO9","BIO10",
                      "BIO11","BIO12","BIO13","BIO14","BIO15",
                      "BIO16","BIO17","BIO18","BIO19")

# bioclim <- cbind(data[,c(2,6,7)][-remove,],rep(1,55))
# colnames(bioclim) <- c("siteCol","Latitude","Longitude","dummy")

temp_pca <- prcomp(bioclim[,c("BIO1","BIO2","BIO3","BIO4","BIO5","BIO6","BIO7","BIO8","BIO9","BIO10","BIO11")],center=TRUE,scale.=TRUE)
summary(temp_pca)
biplot(temp_pca,cex=0.5)
temp_scores <- predict(temp_pca)[,1:4]
colnames(temp_scores) <- c("PC1T","PC2T","PC3T","PC4T")

prec_pca <- prcomp(bioclim[,c("BIO12","BIO13","BIO14","BIO15","BIO16","BIO17","BIO18","BIO19")],center=TRUE,scale.=TRUE)
summary(prec_pca)
biplot(prec_pca,cex=0.5)
prec_scores <- predict(prec_pca)[,1:4]
colnames(prec_scores) <- c("PC1P","PC2P","PC3P","PC4P")

bioclim_final <- cbind(bioclim[,1:3],temp_scores,prec_scores)

# Predictor dissimilarity matrices

# ENV <- read.csv("environmental_distance.csv")[-remove,-remove]
# ENV <- cbind(pheno$siteCol,ENV)
# colnames(ENV) <- c("siteCol",82717,82718,82719,82720,82721,82722,106308,144750,156482,254754,328897,328900,388899,398723,398724,398730,398733,398737,398738,398748,398750,398753,398754,398757,398758,398760,398761,398763,398764,398765,399356,399359,521922,776163,776164,776165,776166,778430,778431,787600,813830,26638,26641,26653,26656,26666,26667,26673,26675,26677,26687,26697,26699,26705,1822)

# PREC <- read.csv("precipitation.dist.csv")[-remove,-remove]
# PREC <- cbind(pheno$siteCol,PREC)
# colnames(PREC) <- c("siteCol",82717,82718,82719,82720,82721,82722,106308,144750,156482,254754,328897,328900,388899,398723,398724,398730,398733,398737,398738,398748,398750,398753,398754,398757,398758,398760,398761,398763,398764,398765,399356,399359,521922,776163,776164,776165,776166,778430,778431,787600,813830,26638,26641,26653,26656,26666,26667,26673,26675,26677,26687,26697,26699,26705,1822)
# 
# TEMP <- read.csv("temperature.dist.csv")[-remove,-remove]
# TEMP <- cbind(pheno$siteCol,TEMP)
# colnames(TEMP) <- c("siteCol",82717,82718,82719,82720,82721,82722,106308,144750,156482,254754,328897,328900,388899,398723,398724,398730,398733,398737,398738,398748,398750,398753,398754,398757,398758,398760,398761,398763,398764,398765,399356,399359,521922,776163,776164,776165,776166,778430,778431,787600,813830,26638,26641,26653,26656,26666,26667,26673,26675,26677,26687,26697,26699,26705,1822)

#LGM <- read.csv("LGM_resistance.csv",header=T)[-remove,-remove]
LGM = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/BELLII_distancematrix_str.csv",
               header=T)
colnames(LGM)[1] = "Sample"
#LGM <- cbind(gene$Sample,LGM)
#colnames(LGM) <- c("Sample","A10","A4","A5","A6","A9","B10","B4","B5","B7","B9","C10","C4","C5","C6","C7","C9","D10","D4","D6","D7","D8","D9","E10","E4","E5","E6","E7","E8","E9","F2","F5","F7","F8","F9","G3","G7","G8","G9","H3","H4","H5","H8","H9")


#IBD <- read.csv("Present_plain_resistance.csv",header=T)[-remove,-remove]
IBD = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/BELLII_distancematrix_geog.csv",
               header=T)
colnames(IBD)[1] = "Sample"
#IBD <- cbind(gene$Sample,IBD)
#colnames(IBD) <- c("Sample","A10","A4","A5","A6","A9","B10","B4","B5","B7","B9","C10","C4","C5","C6","C7","C9","D10","D4","D6","D7","D8","D9","E10","E4","E5","E6","E7","E8","E9","F2","F5","F7","F8","F9","G3","G7","G8","G9","H3","H4","H5","H8","H9")

PRES = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/BELLII_distancematrix_env.csv",
                header=T)
colnames(PRES)[1] = "Sample"
#PRES <- read.csv("Present_resistance.csv",header=T)[-remove,-remove]
#PRES <- cbind(gene$Sample,PRES)
#colnames(PRES) <- c("Sample","A10","A4","A5","A6","A9","B10","B4","B5","B7","B9","C10","C4","C5","C6","C7","C9","D10","D4","D6","D7","D8","D9","E10","E4","E5","E6","E7","E8","E9","F2","F5","F7","F8","F9","G3","G7","G8","G9","H3","H4","H5","H8","H9")


# GDM

#GENE = merge(GENE,gene) ## not sure if here is where it goes

bioclim_final$Sample = rownames(bioclim_final)
bioclim_final = merge(bioclim_final,gene)
rownames(bioclim_final) = bioclim_final$Sample
#bioclim_final = bioclim_final[,c(-1)]

rownames(bioclim_final) = stringr::str_replace_all(rownames(bioclim_final), "[^[:alnum:]]", "")
colnames(bioclim_final) = stringr::str_replace_all(colnames(bioclim_final), "[^[:alnum:]]", "")
rownames(GENE) = stringr::str_replace_all(rownames(GENE), "[^[:alnum:]]", "")
colnames(GENE) = stringr::str_replace_all(colnames(GENE), "[^[:alnum:]]", "")
rownames(IBD) = stringr::str_replace_all(rownames(IBD), "[^[:alnum:]]", "")
colnames(IBD) = stringr::str_replace_all(colnames(IBD), "[^[:alnum:]]", "")
rownames(LGM) = stringr::str_replace_all(rownames(LGM), "[^[:alnum:]]", "")
colnames(LGM) = stringr::str_replace_all(colnames(LGM), "[^[:alnum:]]", "")
rownames(PRES) = stringr::str_replace_all(rownames(PRES), "[^[:alnum:]]", "")
colnames(PRES) = stringr::str_replace_all(colnames(PRES), "[^[:alnum:]]", "")
bioclim_final$Sample = stringr::str_replace_all(bioclim_final$Sample, "[^[:alnum:]]", "")
GENE$Sample = stringr::str_replace_all(GENE$Sample, "[^[:alnum:]]", "")
LGM$Sample = stringr::str_replace_all(LGM$Sample, "[^[:alnum:]]", "")
IBD$Sample = stringr::str_replace_all(IBD$Sample, "[^[:alnum:]]", "")
PRES$Sample = stringr::str_replace_all(PRES$Sample, "[^[:alnum:]]", "")

formatted <- formatsitepair(GENE, 3, siteColumn="Sample", 
                           XColumn="Longitude", YColumn="Latitude", 
                           predData = bioclim_final, 
                           distPreds=list(as.matrix(IBD),as.matrix(LGM),as.matrix(PRES)))
names(formatted)
formatted$distance = scales::rescale(formatted$distance)

model <- gdm(formatted, geo=FALSE, splines=NULL, knots=NULL)
summary(model)
str(model)

model$explained
model$predictors
model$coefficients

length(model$predictors)
plot(model)
plot(model, plot.layout=c(4,3))

modTest <- gdm.varImp(formatted, geo=F, nPerm=100, parallel=T, cores=4)

#Deviance
modTest[[1]]
write.csv(modTest[[1]],"variable_importance_test.csv")

#Importance

#Significance
modTest[[3]]
write.csv(modTest[[3]],"variable_significance_test.csv")


modTest_b <- gdm.varImp(formatted,fullModelOnly = FALSE, geo=F, nPerm=100, parallel=T, cores=4)

barplot(modTest[[2]][,1],main = "full Model")


barplot(modTest[[2]][,2],main = "full Model - 1")
barplot(modTest[[2]][,3],main = "full Model - 2")
barplot(modTest[[2]][,4],main = "full Model - 3")
barplot(modTest[[2]][,5],main = "full Model - 4")
barplot(modTest[[2]][,6],main = "full Model - 5")
barplot(modTest[[2]][,7],main = "full Model - 6")
barplot(modTest[[2]][,8],main = "full Model - 7")
barplot(modTest[[2]][,9],main = "full Model - 8")
barplot(modTest[[2]][,10],main = "full Model - 9")

save.image("PTO.RData")
