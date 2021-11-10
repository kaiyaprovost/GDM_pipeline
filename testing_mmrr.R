library(raster)
library(ecodist)
library(PopGenReport)

gen="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/BY_SPECIES/CSVS/BELLII/BELLII_distancematrix_FULLGENOME.csv.converted"
gen_mat = read.table(gen,header=T,sep=",",row.names = 1)

env="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/BY_SPECIES/CSVS/BELLII/BELLII_distancematrix_env.csv.converted"
ibe = read.table(env,header=T,sep=",",row.names = 1)

ibd_string = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/BY_SPECIES/CSVS/BELLII/BELLII_distancematrix_geog.csv.converted"
ibd = read.table(ibd_string,sep=",",header=T,row.names = 1)

str_string="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/Distances/BY_SPECIES/CSVS/BELLII/BELLII_distancematrix_str.csv.converted"
str = read.table(str_string,sep=",",header=T,row.names = 1)


ibe2=ibe[which(rownames(ibe) %in% rownames(gen_mat)),which(colnames(ibe) %in% colnames(gen_mat))]
ibd2=ibe[which(rownames(ibd) %in% rownames(gen_mat)),which(colnames(ibd) %in% colnames(gen_mat))]
str2=ibe[which(rownames(str) %in% rownames(gen_mat)),which(colnames(str) %in% colnames(gen_mat))]

ibe3 = as.dist(scales::rescale(as.matrix(ibe2)),upper=T,diag=T)
ibd3 = as.dist(scales::rescale(as.matrix(ibd2)),upper=T,diag=T)
str3 = as.dist(scales::rescale(as.matrix(str2)),upper=T,diag=T)
gen3 = as.dist(scales::rescale(as.matrix(gen_mat)),upper=T,diag=T)

## this can solve either ibd or ibe on linear, but not both. 
## it seems like which one of these is NA depends on the order
MRM(formula = gen3~ibe3+ibd3,method="linear") ## no 
MRM(formula = gen3~ibd3,method="linear")
MRM(formula = gen3~ibe3,method="linear")
MRM(formula = gen3~str3,method="linear")
MRM(formula = gen3~str3+ibe3,method="linear") ## no
MRM(formula = gen3~str3+ibd3,method="linear") ## no
MRM(formula = gen3~str3+ibd3+ibe3,method="linear") ## no

MMRR<-function(Y,X,nperm=999){
  #compute regression coefficients and test statistics
  nrowsY<-nrow(Y)
  y<-unfold(Y)
  if(is.null(names(X))) {names(X)<-paste("X",1:length(X),sep="")}
  Xmats<-sapply(X,unfold)
  fit<-lm(y~Xmats)
  coeffs<-fit$coefficients
  summ<-summary(fit)
  r.squared<-summ$r.squared
  tstat<-summ$coefficients[,"t value"]
  
  
  if(length(tstat)<length(names(X))+1){
    tstat[seq(length(tstat)+1,length(names(X))+1)] = NA
  }
  
  Fstat<-summ$fstatistic[1]
  tprob<-rep(1,length(tstat))
  Fprob<-1
  
  #perform permutations
  for(i in 1:nperm){
    rand<-sample(1:nrowsY)
    Yperm<-Y[rand,rand]
    yperm<-unfold(Yperm)
    fit<-lm(yperm~Xmats)
    summ<-summary(fit)
    Fprob<-Fprob+as.numeric(summ$fstatistic[1]>=Fstat)
    tprob<-tprob+as.numeric(abs(summ$coefficients[,"t value"])>=abs(tstat))
  }
  
  #return values
  tp<-tprob/(nperm+1)
  Fp<-Fprob/(nperm+1)
  names(r.squared)<-"r.squared"
  names(coeffs)<-c("Intercept",names(X))
  
  
  names(tstat)<-paste(c("Intercept",names(X)),"(t)",sep="")
  names(tp)<-paste(c("Intercept",names(X)),"(p)",sep="")
  names(Fstat)<-"F-statistic"
  names(Fp)<-"F p-value"
  return(list(r.squared=r.squared,
              coefficients=coeffs,
              tstatistic=tstat,
              tpvalue=tp,
              Fstatistic=Fstat,
              Fpvalue=Fp))
}
# unfold converts the lower diagonal elements of a matrix into a vector
# unfold is called by MMRR
unfold<-function(X){
  x<-vector()
  for(i in 2:nrow(X)) x<-c(x,X[i,1:i-1])
  return(x)
}


## this way will run, but only because I forced the program to insert a NA
MMRR(Y=as.matrix(gen3),X=list("IBD"=as.matrix(ibd3),"IBE"=as.matrix(ibe3)))
MMRR(Y=as.matrix(gen3),X=list("IBE"=as.matrix(ibe3)))
MMRR(Y=as.matrix(gen3),X=list("IBD"=as.matrix(ibd3)))

     