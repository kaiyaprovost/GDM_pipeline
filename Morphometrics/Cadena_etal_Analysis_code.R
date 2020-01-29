
################################################################################################################
################################################################################################################
################################################################################################################
## This script runs analyses for species delimitaion of Darwin's Finches in the genus Geospiza using
## morphological data obtained by H. S. Swarth and used by McKay and Zink (2015. Sisyphean evolution in Darwin's
## finches. Biological Reviews 90: 689-698). The script uses Gaussian mixture models to detect morphological
## groups. Following cKay and Zink (2015), the focus is on six traits: wing length, tail length, bill length,
## bill depth, bill width, and tarsus length.
##
## The main sections of the script are: 
##
## 1) Preliminaries (loading packages and reading, transforming and rotating [with PCA] data).
## 2) Data analysis using Gaussian mixture models on logarithmically transformed data.
################################################################################################################
################################################################################################################
################################################################################################################


################################################################################################################
################################################################################################################
# 1) Preliminaries (loading packages and reading, transforming and rotating [with PCA] data)
################################################################################################################
################################################################################################################

################################################################################################################
# 1.1) Load packages

#load packages
library(mclust)
library(mvtnorm)
library(ellipse)
library(clustvarsel)


################################################################################################################
# 1.2) Read data

#read the morphological data and examine the resulting data frame
path = "/Users/kprovost/Documents/Dissertation/CHAPTER3_TRAITS/morphology/"
fileprefix1 = "Master_Spreadsheet_Morphology_round2_"
fileprefix2 = "KLP_Master_Spreadsheet_Morphology_"
date1 = "21may2019"
date2 = "21_May_2019"
agg <- read.table(paste(path,fileprefix2,"AGGREGATED_",date2,".csv",sep=""), header=T, sep=",")
names(agg)[c(10,12)] = c("TAIL","LONG")
dim(agg)
head(agg)
summary(agg)
summary(agg$SEX)

# morpho.data <- read.table("data.csv", header=T, sep=",")
# dim(morpho.data)
# head(morpho.data)
# summary(morpho.data)
# summary(morpho.data$Sex)

################################################################################################################
# 1.3) Transform data

#perform logarithmic transformation of the six traits of interest:
#Wing, Tail, Blength, Bdepth, Bwidth, Tarsus
# morpho.data.ln <- log(morpho.data[,c(9:14)])
#change column names and examine resulting data frame
# colnames(morpho.data.ln)
# colnames(morpho.data.ln) <- c("Ln (Wing length)", "Ln (Tail length)", "Ln (Bill length)",
# "Ln (Bill depth)", "Ln (Bill width)", "Ln (Tarsus)")
# dim(morpho.data.ln)
# summary(morpho.data.ln)

agg.ln <- log(agg[,c(4:10,41:47)])
agg.ln$SPP = agg$SPP
agg.ln = agg.ln[complete.cases(agg.ln),]
#change column names and examine resulting data frame
colnames(agg.ln)
colnames(agg.ln) = c("LOG_BILL.HEIGHT",
                     "LOG_BILL.LENGTH",
                     "LOG_BILL.WIDTH",
                     "LOG_TARSUS.LENGTH",
                     "LOG_WING.LENGTH.TO.PRIMARIES",
                     "LOG_WING.LENGTH.TO.SECONDARIES",
                     "LOG_TAIL",
                     "LOG_KIPPSINDEX_CALC",
                     "LOG_BEAKBASEAREA_CALC",
                     "LOG_BEAKVOL_CALC",
                     "LOG_BEAKLATERALSURFACE_CALC",
                     "LOG_BEAKTOTALSURFACE_CALC",
                     "LOG_BEAKAREAVOLUME_CALC",
                     "LOG_BEAKTOTAREAVOLUME_CALC",
                     "SPP")
dim(agg.ln)
summary(agg.ln)


################################################################################################################
# 1.4) Principal component analysis (i.e., rotate data)

# 1.5.1) perform PCA using the covariance matrix of the logarithmic
#transformation of the six traits of interest (Wing, Tail, Blength,
#Bdepth, Bwidth, Tarsus) and examine the results
# morpho.data.ln.pca <- prcomp(morpho.data.ln, center = T, scale. = F) #PCA using the covariance matrix
agg.ln.pca <- prcomp(agg.ln, center = T, scale. = F) #PCA using the covariance matrix


#examine PCA results
#get the list of attributes of the R object containing the PCA results (which
#is an R object of class "prcomp")
# attributes(morpho.data.ln.pca)
# morpho.data.ln.pca$scale
# morpho.data.ln.pca$center

attributes(agg.ln.pca)
agg.ln.pca$scale
agg.ln.pca$center

#variance explained by each principal component: 
# summary(morpho.data.ln.pca)
summary(agg.ln.pca)

#summary of the principal components: 
# summary(morpho.data.ln.pca$x)
summary(agg.ln.pca$x)


#examine the "rotation" element of the R object containing the PCA results,
#it shows the coefficients (or "loadings") of each trait on each principal component
# morpho.data.ln.pca$rotation
agg.ln.pca$rotation


################################################################################################################
################################################################################################################
# 2) Data analysis using Gaussian mixture models on morphological axes defined by PCA on a covariance matrix,
# and using variable selection.
################################################################################################################
################################################################################################################

################################################################################################################
# 2.1) Data analysis

#backward variable selection using the PCA of the logarithmic transformation of the six traits of interest
#(Wing, Tail, Blength, Bdepth, Bwidth and Tarsus), and examine results:
mclust.options() #check Mclust options
OptMc <- mclust.options() #save default
mclust.options(hcUse="VARS") #change default as needed
# morpho.data.ln.pca.varsel.back <- clustvarsel(morpho.data.ln.pca$x, G=1:30, search=c("greedy"), direction = c("backward"))
# attributes(morpho.data.ln.pca.varsel.back)
# summary(morpho.data.ln.pca.varsel.back)
# morpho.data.ln.pca.varsel.back$subset
# #PC1 PC2 PC3 PC4 
# #  1   2   3   4 

agg.ln.pca.varsel.back <- clustvarsel(agg.ln.pca$x, G=1:30, search=c("greedy"), direction = c("backward"))
attributes(agg.ln.pca.varsel.back)
summary(agg.ln.pca.varsel.back)
agg.ln.pca.varsel.back$subset
"PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 
1   2   3   4   5   6   7   8   9 "

# morpho.data.ln.pca.varsel.back$steps.info
# #  Variable proposed      BIC BIC difference Type of step Decision
# #1               PC6 5772.130      -64.22597       Remove Accepted
# #2               PC5 3877.614      -29.43604       Remove Accepted
# #3               PC4 3877.614       55.08995       Remove Rejected
# morpho.data.ln.pca.varsel.back$search
# morpho.data.ln.pca.varsel.back$direction

agg.ln.pca.varsel.back$steps.info
"   Variable proposed Type of step BICclust Model  G      BICdiff Decision
1               PC14       Remove 2292.818   EII 19 -12223.59543 Accepted
2               PC13       Remove 1769.704   EII 20 -12011.35821 Accepted
3               PC12       Remove 1373.135   EII 22 -12092.02268 Accepted
4               PC12          Add 1373.135   EII 22 -12488.59240 Rejected
5               PC11       Remove 4124.194   EEV  7 -15194.90594 Accepted
6               PC11          Add 1373.135   EII 22 -15194.90594 Rejected
7               PC10       Remove 2792.715   VVV  5    -37.57227 Accepted
8               PC10          Add 4124.194   EEV  7    -37.57227 Rejected
9                PC5       Remove 2551.142   EEV  8     35.72071 Rejected
10              PC10          Add 4124.194   EEV  7    -37.57227 Rejected"

agg.ln.pca.varsel.back$search
agg.ln.pca.varsel.back$direction

#forward variable selection using logarithmic transformation of the six traits of interest
#(Wing, Tail, Blength, Bdepth, Bwidth and Tarsus), and examine results:
mclust.options() #check Mclust options
#OptMc <- mclust.options() #save default
mclust.options(hcUse="VARS") #change default as needed
# morpho.data.ln.pca.varsel.for <- clustvarsel(morpho.data.ln.pca$x, G=1:30, search=c("greedy"), direction = c("forward"))
# attributes(morpho.data.ln.pca.varsel.for)
# summary(morpho.data.ln.pca.varsel.for)
# morpho.data.ln.pca.varsel.for$subset
#PC2 PC1 PC3 PC4 
#  2   1   3   4 
agg.ln.pca.varsel.for <- clustvarsel(agg.ln.pca$x, G=1:30, search=c("greedy"), direction = c("forward"))
attributes(agg.ln.pca.varsel.for)
summary(agg.ln.pca.varsel.for)
agg.ln.pca.varsel.for$subset
"PC10  PC9  PC2  PC1  PC4  PC7 
  10    9    2    1    4    7 "

#morpho.data.ln.pca.varsel.for$steps.info
#  Variable proposed       BIC BIC difference Type of step Decision
#1               PC2  932.7921      226.01345          Add Accepted
#2               PC1  657.5991      415.82129          Add Accepted
#3               PC3 2079.5264      126.18147          Add Accepted
#4               PC3 2079.5264      126.18318       Remove Rejected
#5               PC4 3877.6139       55.08995          Add Accepted
#6               PC4 3877.6139       62.14693       Remove Rejected
#7               PC5 3877.6139      -29.43604          Add Rejected
#8               PC4 3877.6139       62.14693       Remove Rejected
# morpho.data.ln.pca.varsel.for$search
# morpho.data.ln.pca.varsel.for$direction

agg.ln.pca.varsel.for$steps.info
"   Variable proposed Type of step  BICclust Model  G    BICdiff Decision
1               PC10          Add 1535.9808     V  3  166.92918 Accepted
2                PC9          Add 2830.5509   VVI  4  137.29499 Accepted
3                PC3          Add 2868.6326   EEV  6   72.67624 Accepted
4                PC3       Remove 2830.5509   VVI  4   72.67624 Rejected
5                PC2          Add 2595.1038   EEV  7   56.32131 Accepted
6                PC2       Remove 2868.6326   EEV  6   56.32131 Rejected
7                PC1          Add 1985.6026   VEV  8  200.23962 Accepted
8                PC9       Remove  731.5358   EEV 10   96.79175 Rejected
9                PC4          Add 2052.7799   VEV  8   56.05747 Accepted
10               PC3       Remove 2218.4803   VEE 12 -131.10587 Accepted
11               PC7          Add 2669.4146   VEE 13   69.73445 Accepted
12              PC10       Remove 1273.8577   VEE 10   26.50525 Rejected
13               PC8          Add 3537.3245   VEI 13  -34.90462 Rejected
14              PC10       Remove 1273.8577   VEE 10   26.50525 Rejected"
agg.ln.pca.varsel.for$search
agg.ln.pca.varsel.for$direction

#based on the results above, carry out the Mclust analysis using four characters:
#PC2, PC1, PC3 and PC4:
mclust.options() #check Mclust options and make sure hcUSE="VARS", otherwise change it.
mclust.options(hcUse="VARS")
#Run mclust analysis
# Mcluster.morpho.data.ln.pca.subset <- Mclust(morpho.data.ln.pca$x[,1:4], G=1:30)
# summary(Mcluster.morpho.data.ln.pca.subset)
# attributes(Mcluster.morpho.data.ln.pca.subset)
# #help(mclustModelNames) #in this help page there is information about model names
# #plot(Mcluster.morpho.data.ln.pca.subset$BIC)
# summary(Mcluster.morpho.data.ln.pca.subset$data)
# dim(Mcluster.morpho.data.ln.pca.subset$data)

Mcluster.agg.ln.pca.subset <- Mclust(agg.ln.pca$x[,c(1,2,4,7,9)], G=1:30)
summary(Mcluster.agg.ln.pca.subset)
attributes(Mcluster.agg.ln.pca.subset)
#help(mclustModelNames) #in this help page there is information about model names
#plot(Mcluster.agg.ln.pca.subset$BIC)
summary(Mcluster.agg.ln.pca.subset$data)
dim(Mcluster.agg.ln.pca.subset$data)


#extract BIC values for the best model conditional on the number of groups
# BIC.Best.Model.Per.G <- apply(Mcluster.morpho.data.ln.pca.subset$BIC, 1, max, na.rm=T)
#Mcluster.morpho.data.ln.pca.subset$BIC
BIC.Best.Model.Per.G <- apply(Mcluster.agg.ln.pca.subset$BIC, 1, max, na.rm=T)
#Mcluster.agg.ln.pca.subset$BIC
min(BIC.Best.Model.Per.G) ## G1 -- 10 or 11, should be 11 but curious for 10


# #Obtain second best Mclust model, with 7 morphogroups 
# Mcluster.morpho.data.ln.pca.subset.G7 <- Mclust(morpho.data.ln.pca$x[,1:4], G=7)
# summary(Mcluster.morpho.data.ln.pca.subset.G7)
# attributes(Mcluster.morpho.data.ln.pca.subset.G7)
# #help(mclustModelNames) #in this help page there is information about model names
# Mcluster.morpho.data.ln.pca.subset.G7$BIC
# summary(Mcluster.morpho.data.ln.pca.subset.G7$data)
# dim(Mcluster.morpho.data.ln.pca.subset.G7$data)

#Obtain second best Mclust model, with 11 morphogroups 
Mcluster.agg.ln.pca.subset.G11 <- Mclust(agg.ln.pca$x[,c(1,2,4,7,9)], G=11)
summary(Mcluster.agg.ln.pca.subset.G11)
attributes(Mcluster.agg.ln.pca.subset.G11)
#help(mclustModelNames) #in this help page there is information about model names
Mcluster.agg.ln.pca.subset.G11$BIC
summary(Mcluster.agg.ln.pca.subset.G11$data)
dim(Mcluster.agg.ln.pca.subset.G11$data)

# #Calculate empirical support for the hypothesis of species limits by Lack (1947) 
# H.Lack <- MclustDA(Mcluster.morpho.data.ln.pca.subset$data, class=as.vector(morpho.data$Taxon), G = 1, modelType = "MclustDA")
# #attributes(H.Lack)
# summary(H.Lack)
# #H.Lack$models
# H.Lack$bic

#Calculate empirical support for the hypothesis of species limits by Lack (1947) 
H.Lack <- MclustDA(Mcluster.agg.ln.pca.subset$data, class=as.vector(agg.ln$SPP), G = 1, modelType = "MclustDA")
#attributes(H.Lack)
summary(H.Lack)
#H.Lack$models
H.Lack$bic

#Calculate empirical support for the hypothesis of species limits proposed by the currently accepted taxonomy 
H.Current.taxonomy <- MclustDA(Mcluster.morpho.data.ln.pca.subset$data, class=as.vector(morpho.data$New_Taxonomy), G = 1, modelType = "MclustDA")
#attributes(H.Current.taxonomy)
summary(H.Current.taxonomy)
#H.Current.taxonomy$models
H.Current.taxonomy$bic

################################################################################################################
# 2.2) Graph empirical support for different hypotheses

max.BIC <- max(BIC.Best.Model.Per.G)

par(mar=c(5, 4.5, 2, 2) + 0.1)
plot(1:30, max.BIC-BIC.Best.Model.Per.G[1:30], type="n", bty="n", xlim=c(1,30), ylim=c(900,0), yaxt="n", xaxt="n",
	xlab="Number of morphological groups", ylab=expression(paste("Empirical support (",Delta, "BIC )", sep="")), main="",
	cex.axis=1.5, cex.lab=1.5, cex.main=1.5)
points(2:30, max.BIC-BIC.Best.Model.Per.G[2:30], cex=2, pch=21, col="black", lwd=2)
#show the best Mclust models
points(Mcluster.morpho.data.ln.pca.subset$G, max.BIC-max(BIC.Best.Model.Per.G), pch=19, cex=2)
points(7, max.BIC-sort(BIC.Best.Model.Per.G, decreasing=T)[2], pch=19, cex=2)
#show the hypothesis by McKay and Zink (2015)
points(1, max.BIC-BIC.Best.Model.Per.G[1], pch=24, bg="black", cex=2)
#show the hypothesis by Lack (1947)
points(length(unique(morpho.data$Taxon)), max.BIC-H.Lack$bic, pch=25, bg="black", cex=2)
#show the hypothesis according to current taxonomy
points(length(unique(morpho.data$New_Taxonomy)), max.BIC-H.Current.taxonomy$bic, pch=15, bg="black", cex=2)
#add abscissa
axis(1, at=seq(1,30,1), labels=F, tcl=-0.4)
axis(1, at=c(1,seq(5,30,5)), labels=T, tcl=-0.7, cex.axis=1.5)
axis(2, at=seq(1000,0,-50), labels=F, tcl=-0.4)
axis(2, at=seq(1000,0,-200), tcl=-0.7, cex.axis=1.5)
#add legend
legend(10.3,400, c("Sisyphean evolution", "(McKay & Zink 2015)", "Previous taxonomy", "(Lack 1947, Rising et al. 2011)", "Current taxonomy", "(Remsen et al. 2017)", expression(paste("Best Mclust models (", Delta,BIC<=1.26, ")", sep="")), expression(paste("Other Mclust models (", Delta, BIC>20, ")", sep=""))),
pch=c(24,NA,25,NA,15,NA,19,21), pt.bg=c("black", NA, "black", NA,"black",NA, "black", NA), pt.cex=2, cex=1.3, pt.lwd=c(1,0,1,0,1,0,1,2))
#add arrows
arrows(7, 1000, 
	7, max.BIC-sort(BIC.Best.Model.Per.G, decreasing=T)[2], length=0, lty=3, lwd=1.5)
arrows(Mcluster.morpho.data.ln.pca.subset$G, 1000, 
	Mcluster.morpho.data.ln.pca.subset$G, max.BIC-max(BIC.Best.Model.Per.G), length=0, lty=3, lwd=1.5)
arrows(length(unique(morpho.data$Taxon)), 1000, 
	length(unique(morpho.data$Taxon)),max.BIC-H.Lack$bic, length=0, lty=3, lwd=1.5)
arrows(length(unique(morpho.data$New_Taxonomy)), 1000, 
	length(unique(morpho.data$New_Taxonomy)),max.BIC-H.Current.taxonomy$bic, length=0, lty=3, lwd=1.5)
arrows(1, 1000, 
	1,max.BIC-BIC.Best.Model.Per.G[1], length=0, lty=3, lwd=1.5)
#add horizontal line
abline(h=max(BIC.Best.Model.Per.G)-10, lty=2)


################################################################################################################
# 2.3) Plot phenotypic traits and specimens to show morphological groups according
# to the best Mclust model with eight morphological groups

#define morphogroup colors
morpho.groups.colors <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

col.BestModel <- rep(morpho.groups.colors[1], length(morpho.data$New_Taxonomy)) 
col.BestModel[Mcluster.morpho.data.ln.pca.subset$classification==2] <- morpho.groups.colors[2]
col.BestModel[Mcluster.morpho.data.ln.pca.subset$classification==3] <- morpho.groups.colors[3]
col.BestModel[Mcluster.morpho.data.ln.pca.subset$classification==4] <- morpho.groups.colors[4]
col.BestModel[Mcluster.morpho.data.ln.pca.subset$classification==5] <- morpho.groups.colors[5]
col.BestModel[Mcluster.morpho.data.ln.pca.subset$classification==6] <- morpho.groups.colors[6]
col.BestModel[Mcluster.morpho.data.ln.pca.subset$classification==7] <- morpho.groups.colors[7]
col.BestModel[Mcluster.morpho.data.ln.pca.subset$classification==8] <- morpho.groups.colors[8]


#gather for all morphological groups in the best Mclust model estimates of
#i) the mean vector and ii) the variance-covariance matrix
m.all <- Mcluster.morpho.data.ln.pca.subset$parameters$mean #the mean vector
Z.all <- Mcluster.morpho.data.ln.pca.subset$parameters$variance$sigma #the variance-covariance matrices

#Make scatterplot. Repeat this code for the six combinations of traits changin below trait.x and trait.y combining 1-4
trait.x <- 1
trait.y <- 2
plot(Mcluster.morpho.data.ln.pca.subset$data[,trait.x], Mcluster.morpho.data.ln.pca.subset$data[,trait.y],
xlab=colnames(Mcluster.morpho.data.ln.pca.subset$data)[trait.x], ylab=colnames(Mcluster.morpho.data.ln.pca.subset$data)[trait.y],
cex.axis=1.5, cex.lab=1.5, bty="n", col=col.BestModel, pch=21,
cex=1.5, lwd=2)
for(i in 1:Mcluster.morpho.data.ln.pca.subset$G)
{
	points(ellipse(Z.all[c(trait.x,trait.y),c(trait.x,trait.y),i], centre = m.all[c(trait.x,trait.y),i], level = 0.95, npoints = 10000), type="l")
}

#plot loadings PC1 and PC2
#arrows
plot(0,0, xlim=c(-1,1), ylim=c(-1,1), asp=1, type="n", bty="n",
xlab="Loadings PC1", ylab="Loadings PC2",
cex.axis=1.5, cex.lab=1.5)
for (i in 1:6){
arrows(0, 0, x1 = morpho.data.ln.pca$rotation[,1][i], y1 = morpho.data.ln.pca$rotation[,2][i], length = 0.1, angle = 20, code = 2, col = "black", lwd=1)
}

#plot circle for equal contribution of traits.
x <- seq(-(2/6)^0.5, (2/6)^0.5, 1e-4) 
equal_contribution<-((2/6)-x^2)^0.5
points(x, equal_contribution, type="l", lty=1, col="gray70")
points(x, -equal_contribution, type="l", lty=1, col="gray70")

#Name arrows by trait
rownames(morpho.data.ln.pca$rotation)[c(3,4,5)]
place_out.x <- rep(1.5,6)
place_out.x[3] <- 1
place_out.x[4] <- 1.1
place_out.x[5] <- 1.35
place_out.y <- rep(1.3,6)
place_out.y[3] <- 1.15
place_out.y[4] <- 1.25
place_out.y[5] <- 0.3
na <- c("Wing length", "Tail length", "Bill length", "Bill depth", "Bill", "Tarsus length")
for (i in c(3,4,5)){
text(place_out.x[i]*morpho.data.ln.pca$rotation[i,1], place_out.y[i]*morpho.data.ln.pca$rotation[i,2], labels = na[i], cex=1)
}
place_out.x[5] <- 1.35
place_out.y[5] <- 0.8
na <- c("Wing length", "Tail length", "Bill length", "Bill depth", "width", "Tarsus length")
for (i in c(5)){
text(place_out.x[i]*morpho.data.ln.pca$rotation[i,1], place_out.y[i]*morpho.data.ln.pca$rotation[i,2], labels = na[i], cex=1)
}

#plot loadings PC3 and PC4
#arrows
plot(0,0, xlim=c(-1,1), ylim=c(-1,1), asp=1, type="n", bty="n",
xlab="Loadings PC3", ylab="Loadings PC4",
cex.axis=1.5, cex.lab=1.5)
for (i in 1:6){
arrows(0, 0, x1 = morpho.data.ln.pca$rotation[,3][i], y1 = morpho.data.ln.pca$rotation[,4][i], length = 0.1, angle = 20, code = 2, col = "black", lwd=1)
}

#plot circle for equal contribution of traits.
x <- seq(-(2/6)^0.5, (2/6)^0.5, 1e-4)
equal_contribution <- ((2/6)-x^2)^0.5
points(x, equal_contribution, type="l", lty=1, col="gray70")
points(x, -equal_contribution, type="l", lty=1, col="gray70")

#Name arrows by trait
rownames(morpho.data.ln.pca$rotation)[c(3,4,5)]
place_out.x <- rep(1,6)
place_out.y <- rep(1.1,6)
place_out.y[2] <- 1.15
na <- c("Wing length", "Tail length", "Bill length", "Bill depth", "Bill width", "Tarsus length")
for (i in c(2,6)){
text(place_out.x[i]*morpho.data.ln.pca$rotation[i,3], place_out.y[i]*morpho.data.ln.pca$rotation[i,4], labels = na[i], cex=1)
}

#plot only legend
par(mar=c(0, 0, 0, 0))
plot(c(0,1), c(0,1), type="n", axes=F, bty="n", 
    xlab="", ylab="")
legend(0.1, 1,
	paste("Morphological group", 1:8),
	col=morpho.groups.colors[1:8], 	pch=21, pt.lwd=2, pt.cex=1.5, cex=1.45, bty="o")


################################################################################################################
# 2.4) Plot phenotypic traits and specimens to show morphological groups according
# to the best Mclust model with seven morphological groups

#define morphogroup colors
morpho.groups.colors <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

col.BestModel <- rep(morpho.groups.colors[1], length(morpho.data.raw$New_Taxonomy))
col.BestModel[Mcluster.morpho.data.ln.pca.subset.G7$classification==2] <- morpho.groups.colors[2]
col.BestModel[Mcluster.morpho.data.ln.pca.subset.G7$classification==3] <- morpho.groups.colors[3]
col.BestModel[Mcluster.morpho.data.ln.pca.subset.G7$classification==4] <- morpho.groups.colors[4]
col.BestModel[Mcluster.morpho.data.ln.pca.subset.G7$classification==5] <- morpho.groups.colors[5]
col.BestModel[Mcluster.morpho.data.ln.pca.subset.G7$classification==6] <- morpho.groups.colors[6]
col.BestModel[Mcluster.morpho.data.ln.pca.subset.G7$classification==7] <- morpho.groups.colors[7]


#gather for all morphological groups in the best Mclust model estimates of
#i) the mean vector and ii) the variance-covariance matrix
m.all <- Mcluster.morpho.data.ln.pca.subset.G7$parameters$mean #the mean vector
Z.all <- Mcluster.morpho.data.ln.pca.subset.G7$parameters$variance$sigma #the variance-covariance matrices

#Make scatterplot. Repeat this code for the six combinations of traits changin below trait.x and trait.y combining 1-4
trait.x <- 1
trait.y <- 2
plot(Mcluster.morpho.data.ln.pca.subset.G7$data[,trait.x], Mcluster.morpho.data.ln.pca.subset.G7$data[,trait.y],
xlab=colnames(Mcluster.morpho.data.ln.pca.subset.G7$data)[trait.x], ylab=colnames(Mcluster.morpho.data.ln.pca.subset.G7$data)[trait.y],
cex.axis=1.5, cex.lab=1.5, bty="n", col=col.BestModel, pch=21,
cex=1.5, lwd=2)
for(i in 1:Mcluster.morpho.data.ln.pca.subset.G7$G)
{
	points(ellipse(Z.all[c(trait.x,trait.y),c(trait.x,trait.y),i], centre = m.all[c(trait.x,trait.y),i], level = 0.95, npoints = 10000), type="l")
}

#plot loadings PC1 and PC2
#arrows
plot(0,0, xlim=c(-1,1), ylim=c(-1,1), asp=1, type="n", bty="n",
xlab="Loadings PC1", ylab="Loadings PC2",
cex.axis=1.5, cex.lab=1.5)
for (i in 1:6){
arrows(0, 0, x1 = morpho.data.ln.pca$rotation[,1][i], y1 = morpho.data.ln.pca$rotation[,2][i], length = 0.1, angle = 20, code = 2, col = "black", lwd=1)
}

#plot circle for equal contribution of traits.
x <- seq(-(2/6)^0.5, (2/6)^0.5, 1e-4) 
equal_contribution<-((2/6)-x^2)^0.5
points(x, equal_contribution, type="l", lty=1, col="gray70")
points(x, -equal_contribution, type="l", lty=1, col="gray70")

#name arrows by trait
rownames(morpho.data.ln.pca$rotation)[c(3,4,5)]
place_out.x <- rep(1.5,6)
place_out.x[3] <- 1
place_out.x[4] <- 1.1
place_out.x[5] <- 1.35
place_out.y <- rep(1.3,6)
place_out.y[3] <- 1.15
place_out.y[4] <- 1.25
place_out.y[5] <- 0.3
na <- c("Wing length", "Tail length", "Bill length", "Bill depth", "Bill", "Tarsus length")
for (i in c(3,4,5)){
text(place_out.x[i]*morpho.data.ln.pca$rotation[i,1], place_out.y[i]*morpho.data.ln.pca$rotation[i,2], labels = na[i], cex=1)
}
place_out.x[5] <- 1.35
place_out.y[5] <- 0.8
na <- c("Wing length", "Tail length", "Bill length", "Bill depth", "width", "Tarsus length")
for (i in c(5)){
text(place_out.x[i]*morpho.data.ln.pca$rotation[i,1], place_out.y[i]*morpho.data.ln.pca$rotation[i,2], labels = na[i], cex=1)
}

#plot loadings PC3 and PC4
#arrows
plot(0,0, xlim=c(-1,1), ylim=c(-1,1), asp=1, type="n", bty="n",
xlab="Loadings PC3", ylab="Loadings PC4",
cex.axis=1.5, cex.lab=1.5)
for (i in 1:6){
arrows(0, 0, x1 = morpho.data.ln.pca$rotation[,3][i], y1 = morpho.data.ln.pca$rotation[,4][i], length = 0.1, angle = 20, code = 2, col = "black", lwd=1)
}

#plot circle for equal contribution of traits.
x <- seq(-(2/6)^0.5, (2/6)^0.5, 1e-4)
equal_contribution <- ((2/6)-x^2)^0.5
points(x, equal_contribution, type="l", lty=1, col="gray70")
points(x, -equal_contribution, type="l", lty=1, col="gray70")

#name arrows by trait
rownames(morpho.data.ln.pca$rotation)[c(3,4,5)]
place_out.x <- rep(1,6)
place_out.y <- rep(1.1,6)
place_out.y[2] <- 1.15
na <- c("Wing length", "Tail length", "Bill length", "Bill depth", "Bill width", "Tarsus length")
for (i in c(2,6)){
text(place_out.x[i]*morpho.data.ln.pca$rotation[i,3], place_out.y[i]*morpho.data.ln.pca$rotation[i,4], labels = na[i], cex=1)
}

#plot legend
par(mar=c(0, 0, 0, 0))
plot(c(0,1), c(0,1), type="n", axes=F, bty="n", 
    xlab="", ylab="")
legend(0.1, 1,
	paste("Morphological group", 1:7),
	col=morpho.groups.colors[1:7], 	pch=21, pt.lwd=2, pt.cex=1.5, cex=1.45, bty="o")

################################################################################################################
# 2.5)Compare assignment of specimens to groups between the best Mclust model 
# with eight groups and the hypothesis according to the current taxonomy
# and with the hypothesis by Lack (1947)

#define morphogroup colors
morpho.groups.colors.hist <- c("#999999", "black", "#E69F00", "black", "#56B4E9", "black", "#009E73", "black", "#F0E442", "black", "#0072B2", "black", "#D55E00", "black", "#CC79A7")

#set up the layout of the plot
par(mfrow=c(5,3))
par(oma=c(4,4,0,0))
par(mar=c(2,2,2,1))

hist(Mcluster.morpho.data.ln.pca.subset$classification[morpho.data$New_Taxonomy=="G. acutirostris"], 
breaks=seq(0.75, Mcluster.morpho.data.ln.pca.subset$G+0.25, 0.5),
main=expression(italic("G. acutirostris")), xlab="", ylab="",
cex.lab=1.5, cex.axis=1.5, col=morpho.groups.colors.hist, xaxt="n") 
axis(1, at=1:Mcluster.morpho.data.ln.pca.subset$G, labels=F)
axis(1, at=seq(1,Mcluster.morpho.data.ln.pca.subset$G, 3), labels=T, cex.axis=1.5)

hist(Mcluster.morpho.data.ln.pca.subset$classification[morpho.data$New_Taxonomy=="G. conirostris"], 
breaks=seq(0.75, Mcluster.morpho.data.ln.pca.subset$G+0.25, 0.5),
main=expression(italic("G. conirostris")), xlab="", ylab="",
cex.lab=1.5, cex.axis=1.5, col=morpho.groups.colors.hist, xaxt="n") 
axis(1, at=1:Mcluster.morpho.data.ln.pca.subset$G, labels=F)
axis(1, at=seq(1,Mcluster.morpho.data.ln.pca.subset$G, 3), labels=T, cex.axis=1.5)

hist(Mcluster.morpho.data.ln.pca.subset$classification[morpho.data$New_Taxonomy=="G. difficilis"], 
breaks=seq(0.75, Mcluster.morpho.data.ln.pca.subset$G+0.25, 0.5),
main=expression(italic("G. difficilis")), xlab="", ylab="",
cex.lab=1.5, cex.axis=1.5, col=morpho.groups.colors.hist, xaxt="n") 
axis(1, at=1:Mcluster.morpho.data.ln.pca.subset$G, labels=F)
axis(1, at=seq(1,Mcluster.morpho.data.ln.pca.subset$G, 3), labels=T, cex.axis=1.5)

hist(Mcluster.morpho.data.ln.pca.subset$classification[morpho.data$New_Taxonomy=="G. fortis"], 
breaks=seq(0.75, Mcluster.morpho.data.ln.pca.subset$G+0.25, 0.5),
main=expression(italic("G. fortis")), xlab="", ylab="",
cex.lab=1.5, cex.axis=1.5, col=morpho.groups.colors.hist, xaxt="n") 
axis(1, at=1:Mcluster.morpho.data.ln.pca.subset$G, labels=F)
axis(1, at=seq(1,Mcluster.morpho.data.ln.pca.subset$G, 3), labels=T, cex.axis=1.5)

hist(Mcluster.morpho.data.ln.pca.subset$classification[morpho.data$New_Taxonomy=="G. fuliginosa"], 
breaks=seq(0.75, Mcluster.morpho.data.ln.pca.subset$G+0.25, 0.5),
main=expression(italic("G. fuliginosa")), xlab="", ylab="",
cex.lab=1.5, cex.axis=1.5, col=morpho.groups.colors.hist, xaxt="n") 
axis(1, at=1:Mcluster.morpho.data.ln.pca.subset$G, labels=F)
axis(1, at=seq(1,Mcluster.morpho.data.ln.pca.subset$G, 3), labels=T, cex.axis=1.5)

hist(Mcluster.morpho.data.ln.pca.subset$classification[morpho.data$New_Taxonomy=="G. magnirostris"], 
breaks=seq(0.75, Mcluster.morpho.data.ln.pca.subset$G+0.25, 0.5),
main=expression(italic("G. magnirostris")), xlab="", ylab="",
cex.lab=1.5, cex.axis=1.5, col=morpho.groups.colors.hist, xaxt="n") 
axis(1, at=1:Mcluster.morpho.data.ln.pca.subset$G, labels=F)
axis(1, at=seq(1,Mcluster.morpho.data.ln.pca.subset$G, 3), labels=T, cex.axis=1.5)

hist(Mcluster.morpho.data.ln.pca.subset$classification[morpho.data$New_Taxonomy=="G. propinqua"], 
breaks=seq(0.75, Mcluster.morpho.data.ln.pca.subset$G+0.25, 0.5),
main=expression(italic("G. propinqua")), xlab="", ylab="",
cex.lab=1.5, cex.axis=1.5, col=morpho.groups.colors.hist, xaxt="n") 
axis(1, at=1:Mcluster.morpho.data.ln.pca.subset$G, labels=F)
axis(1, at=seq(1,Mcluster.morpho.data.ln.pca.subset$G, 3), labels=T, cex.axis=1.5)

hist(Mcluster.morpho.data.ln.pca.subset$classification[morpho.data$New_Taxonomy=="G. scandens"], 
breaks=seq(0.75, Mcluster.morpho.data.ln.pca.subset$G+0.25, 0.5),
main=expression(italic("G. scandens")), xlab="", ylab="",
cex.lab=1.5, cex.axis=1.5, col=morpho.groups.colors.hist, xaxt="n") 
axis(1, at=1:Mcluster.morpho.data.ln.pca.subset$G, labels=F)
axis(1, at=seq(1,Mcluster.morpho.data.ln.pca.subset$G, 3), labels=T, cex.axis=1.5)

hist(Mcluster.morpho.data.ln.pca.subset$classification[morpho.data$New_Taxonomy=="G. septentrionalis"], 
breaks=seq(0.75, Mcluster.morpho.data.ln.pca.subset$G+0.25, 0.5),
main=expression(italic("G. septentrionalis")), xlab="", ylab="",
cex.lab=1.5, cex.axis=1.5, col=morpho.groups.colors.hist, xaxt="n") 
axis(1, at=1:Mcluster.morpho.data.ln.pca.subset$G, labels=F)
axis(1, at=seq(1,Mcluster.morpho.data.ln.pca.subset$G, 3), labels=T, cex.axis=1.5)

# Histograms of assignment of specimens to groups between the best Mclust model
# with eigth morphological groups and the hypothesis by Lack (1947)

#define morphogroup colors
morpho.groups.colors.hist <- c("#999999", "black", "#E69F00", "black", "#56B4E9", "black", "#009E73", "black", "#F0E442", "black", "#0072B2", "black", "#D55E00", "black", "#CC79A7")


hist(Mcluster.morpho.data.ln.pca.subset$classification[morpho.data$Taxon=="G. conirostris"], 
breaks=seq(0.75, Mcluster.morpho.data.ln.pca.subset$G+0.25, 0.5),
main=expression(italic("G. conirostris")), xlab="", ylab="",
cex.lab=1.5, cex.axis=1.5, col=morpho.groups.colors.hist, xaxt="n") 
axis(1, at=1:Mcluster.morpho.data.ln.pca.subset$G, labels=F)
axis(1, at=seq(1,Mcluster.morpho.data.ln.pca.subset$G, 3), labels=T, cex.axis=1.5)

hist(Mcluster.morpho.data.ln.pca.subset$classification[morpho.data$Taxon=="G. difficilis"], 
breaks=seq(0.75, Mcluster.morpho.data.ln.pca.subset$G+0.25, 0.5),
main=expression(italic("G. difficilis")), xlab="", ylab="",
cex.lab=1.5, cex.axis=1.5, col=morpho.groups.colors.hist, xaxt="n") 
axis(1, at=1:Mcluster.morpho.data.ln.pca.subset$G, labels=F)
axis(1, at=seq(1,Mcluster.morpho.data.ln.pca.subset$G, 3), labels=T, cex.axis=1.5)

hist(Mcluster.morpho.data.ln.pca.subset$classification[morpho.data$Taxon=="G. fortis"], 
breaks=seq(0.75, Mcluster.morpho.data.ln.pca.subset$G+0.25, 0.5),
main=expression(italic("G. fortis")), xlab="", ylab="",
cex.lab=1.5, cex.axis=1.5, col=morpho.groups.colors.hist, xaxt="n") 
axis(1, at=1:Mcluster.morpho.data.ln.pca.subset$G, labels=F)
axis(1, at=seq(1,Mcluster.morpho.data.ln.pca.subset$G, 3), labels=T, cex.axis=1.5)

hist(Mcluster.morpho.data.ln.pca.subset$classification[morpho.data$Taxon=="G. fuliginosa"], 
breaks=seq(0.75, Mcluster.morpho.data.ln.pca.subset$G+0.25, 0.5),
main=expression(italic("G. fuliginosa")), xlab="", ylab="",
cex.lab=1.5, cex.axis=1.5, col=morpho.groups.colors.hist, xaxt="n") 
axis(1, at=1:Mcluster.morpho.data.ln.pca.subset$G, labels=F)
axis(1, at=seq(1,Mcluster.morpho.data.ln.pca.subset$G, 3), labels=T, cex.axis=1.5)

hist(Mcluster.morpho.data.ln.pca.subset$classification[morpho.data$Taxon=="G. magnirostris"], 
breaks=seq(0.75, Mcluster.morpho.data.ln.pca.subset$G+0.25, 0.5),
main=expression(italic("G. magnirostris")), xlab="", ylab="",
cex.lab=1.5, cex.axis=1.5, col=morpho.groups.colors.hist, xaxt="n") 
axis(1, at=1:Mcluster.morpho.data.ln.pca.subset$G, labels=F)
axis(1, at=seq(1,Mcluster.morpho.data.ln.pca.subset$G, 3), labels=T, cex.axis=1.5)

hist(Mcluster.morpho.data.ln.pca.subset$classification[morpho.data$Taxon=="G. scandens"], 
breaks=seq(0.75, Mcluster.morpho.data.ln.pca.subset$G+0.25, 0.5),
main=expression(italic("G. scandens")), xlab="", ylab="",
cex.lab=1.5, cex.axis=1.5, col=morpho.groups.colors.hist, xaxt="n") 
axis(1, at=1:Mcluster.morpho.data.ln.pca.subset$G, labels=F)
axis(1, at=seq(1,Mcluster.morpho.data.ln.pca.subset$G, 3), labels=T, cex.axis=1.5)

mtext('Morphological group', side = 1, outer = TRUE, line = 2, cex = 1.3)
mtext('Specimens', side = 2, outer = TRUE, line = 2, cex = 1.3)

################################################################################################################
# 2.6)Compare assignment of specimens to groups between the best Mclust model 
# with seven groups and the hypothesis according to the current taxonomy 
# and with the hypothesis by Lack (1947)

#define morphogroup colors
morpho.groups.colors.hist <- c("#999999", "black", "#E69F00", "black", "#56B4E9", "black", "#009E73", "black", "#F0E442", "black", "#0072B2", "black", "#D55E00")

par(mfrow=c(5,3))
par(oma=c(4,4,0,0))
par(mar=c(2,2,2,1))

hist(Mcluster.morpho.data.ln.pca.subset.G7$classification[morpho.data$New_Taxonomy=="G. acutirostris"], 
breaks=seq(0.75, Mcluster.morpho.data.ln.pca.subset.G7$G+0.25, 0.5),
main=expression(italic("G. acutirostris")), xlab="", ylab="",
cex.lab=1.5, cex.axis=1.5, col=morpho.groups.colors.hist, xaxt="n") 
axis(1, at=1:Mcluster.morpho.data.ln.pca.subset.G7$G, labels=F)
axis(1, at=seq(1,Mcluster.morpho.data.ln.pca.subset.G7$G, 3), labels=T, cex.axis=1.5)

hist(Mcluster.morpho.data.ln.pca.subset.G7$classification[morpho.data$New_Taxonomy=="G. conirostris"], 
breaks=seq(0.75, Mcluster.morpho.data.ln.pca.subset.G7$G+0.25, 0.5),
main=expression(italic("G. conirostris")), xlab="", ylab="",
cex.lab=1.5, cex.axis=1.5, col=morpho.groups.colors.hist, xaxt="n") 
axis(1, at=1:Mcluster.morpho.data.ln.pca.subset.G7$G, labels=F)
axis(1, at=seq(1,Mcluster.morpho.data.ln.pca.subset.G7$G, 3), labels=T, cex.axis=1.5)

hist(Mcluster.morpho.data.ln.pca.subset.G7$classification[morpho.data$New_Taxonomy=="G. difficilis"], 
breaks=seq(0.75, Mcluster.morpho.data.ln.pca.subset.G7$G+0.25, 0.5),
main=expression(italic("G. difficilis")), xlab="", ylab="",
cex.lab=1.5, cex.axis=1.5, col=morpho.groups.colors.hist, xaxt="n") 
axis(1, at=1:Mcluster.morpho.data.ln.pca.subset.G7$G, labels=F)
axis(1, at=seq(1,Mcluster.morpho.data.ln.pca.subset.G7$G, 3), labels=T, cex.axis=1.5)

hist(Mcluster.morpho.data.ln.pca.subset.G7$classification[morpho.data$New_Taxonomy=="G. fortis"], 
breaks=seq(0.75, Mcluster.morpho.data.ln.pca.subset.G7$G+0.25, 0.5),
main=expression(italic("G. fortis")), xlab="", ylab="",
cex.lab=1.5, cex.axis=1.5, col=morpho.groups.colors.hist, xaxt="n") 
axis(1, at=1:Mcluster.morpho.data.ln.pca.subset.G7$G, labels=F)
axis(1, at=seq(1,Mcluster.morpho.data.ln.pca.subset.G7$G, 3), labels=T, cex.axis=1.5)

hist(Mcluster.morpho.data.ln.pca.subset.G7$classification[morpho.data$New_Taxonomy=="G. fuliginosa"], 
breaks=seq(0.75, Mcluster.morpho.data.ln.pca.subset.G7$G+0.25, 0.5),
main=expression(italic("G. fuliginosa")), xlab="", ylab="",
cex.lab=1.5, cex.axis=1.5, col=morpho.groups.colors.hist, xaxt="n") 
axis(1, at=1:Mcluster.morpho.data.ln.pca.subset.G7$G, labels=F)
axis(1, at=seq(1,Mcluster.morpho.data.ln.pca.subset.G7$G, 3), labels=T, cex.axis=1.5)

hist(Mcluster.morpho.data.ln.pca.subset.G7$classification[morpho.data$New_Taxonomy=="G. magnirostris"], 
breaks=seq(0.75, Mcluster.morpho.data.ln.pca.subset.G7$G+0.25, 0.5),
main=expression(italic("G. magnirostris")), xlab="", ylab="",
cex.lab=1.5, cex.axis=1.5, col=morpho.groups.colors.hist, xaxt="n") 
axis(1, at=1:Mcluster.morpho.data.ln.pca.subset.G7$G, labels=F)
axis(1, at=seq(1,Mcluster.morpho.data.ln.pca.subset.G7$G, 3), labels=T, cex.axis=1.5)

hist(Mcluster.morpho.data.ln.pca.subset.G7$classification[morpho.data$New_Taxonomy=="G. propinqua"], 
breaks=seq(0.75, Mcluster.morpho.data.ln.pca.subset.G7$G+0.25, 0.5),
main=expression(italic("G. propinqua")), xlab="", ylab="",
cex.lab=1.5, cex.axis=1.5, col=morpho.groups.colors.hist, xaxt="n") 
axis(1, at=1:Mcluster.morpho.data.ln.pca.subset.G7$G, labels=F)
axis(1, at=seq(1,Mcluster.morpho.data.ln.pca.subset.G7$G, 3), labels=T, cex.axis=1.5)

hist(Mcluster.morpho.data.ln.pca.subset.G7$classification[morpho.data$New_Taxonomy=="G. scandens"], 
breaks=seq(0.75, Mcluster.morpho.data.ln.pca.subset.G7$G+0.25, 0.5),
main=expression(italic("G. scandens")), xlab="", ylab="",
cex.lab=1.5, cex.axis=1.5, col=morpho.groups.colors.hist, xaxt="n") 
axis(1, at=1:Mcluster.morpho.data.ln.pca.subset.G7$G, labels=F)
axis(1, at=seq(1,Mcluster.morpho.data.ln.pca.subset.G7$G, 3), labels=T, cex.axis=1.5)

hist(Mcluster.morpho.data.ln.pca.subset.G7$classification[morpho.data$New_Taxonomy=="G. septentrionalis"], 
breaks=seq(0.75, Mcluster.morpho.data.ln.pca.subset.G7$G+0.25, 0.5),
main=expression(italic("G. septentrionalis")), xlab="", ylab="",
cex.lab=1.5, cex.axis=1.5, col=morpho.groups.colors.hist, xaxt="n") 
axis(1, at=1:Mcluster.morpho.data.ln.pca.subset.G7$G, labels=F)
axis(1, at=seq(1,Mcluster.morpho.data.ln.pca.subset.G7$G, 3), labels=T, cex.axis=1.5)

# Histograms of assignment of specimens to groups between the best Mclust model
# with eigth morphological groups and the hypothesis by Lack (1947)

#define morphogroup colors
morpho.groups.colors.hist <- c("#999999", "black", "#E69F00", "black", "#56B4E9", "black", "#009E73", "black", "#F0E442", "black", "#0072B2", "black", "#D55E00")

hist(Mcluster.morpho.data.ln.pca.subset.G7$classification[morpho.data$Taxon=="G. conirostris"], 
breaks=seq(0.75, Mcluster.morpho.data.ln.pca.subset.G7$G+0.25, 0.5),
main=expression(italic("G. conirostris")), xlab="", ylab="",
cex.lab=1.5, cex.axis=1.5, col=morpho.groups.colors.hist, xaxt="n") 
axis(1, at=1:Mcluster.morpho.data.ln.pca.subset.G7$G, labels=F)
axis(1, at=seq(1,Mcluster.morpho.data.ln.pca.subset.G7$G, 3), labels=T, cex.axis=1.5)

hist(Mcluster.morpho.data.ln.pca.subset.G7$classification[morpho.data$Taxon=="G. difficilis"], 
breaks=seq(0.75, Mcluster.morpho.data.ln.pca.subset.G7$G+0.25, 0.5),
main=expression(italic("G. difficilis")), xlab="", ylab="",
cex.lab=1.5, cex.axis=1.5, col=morpho.groups.colors.hist, xaxt="n") 
axis(1, at=1:Mcluster.morpho.data.ln.pca.subset.G7$G, labels=F)
axis(1, at=seq(1,Mcluster.morpho.data.ln.pca.subset.G7$G, 3), labels=T, cex.axis=1.5)

hist(Mcluster.morpho.data.ln.pca.subset.G7$classification[morpho.data$Taxon=="G. fortis"], 
breaks=seq(0.75, Mcluster.morpho.data.ln.pca.subset.G7$G+0.25, 0.5),
main=expression(italic("G. fortis")), xlab="", ylab="",
cex.lab=1.5, cex.axis=1.5, col=morpho.groups.colors.hist, xaxt="n") 
axis(1, at=1:Mcluster.morpho.data.ln.pca.subset.G7$G, labels=F)
axis(1, at=seq(1,Mcluster.morpho.data.ln.pca.subset.G7$G, 3), labels=T, cex.axis=1.5)

hist(Mcluster.morpho.data.ln.pca.subset.G7$classification[morpho.data$Taxon=="G. fuliginosa"], 
breaks=seq(0.75, Mcluster.morpho.data.ln.pca.subset.G7$G+0.25, 0.5),
main=expression(italic("G. fuliginosa")), xlab="", ylab="",
cex.lab=1.5, cex.axis=1.5, col=morpho.groups.colors.hist, xaxt="n") 
axis(1, at=1:Mcluster.morpho.data.ln.pca.subset.G7$G, labels=F)
axis(1, at=seq(1,Mcluster.morpho.data.ln.pca.subset.G7$G, 3), labels=T, cex.axis=1.5)

hist(Mcluster.morpho.data.ln.pca.subset$classification[morpho.data$Taxon=="G. magnirostris"], 
breaks=seq(0.75, Mcluster.morpho.data.ln.pca.subset.G7$G+0.25, 0.5),
main=expression(italic("G. magnirostris")), xlab="", ylab="",
cex.lab=1.5, cex.axis=1.5, col=morpho.groups.colors.hist, xaxt="n") 
axis(1, at=1:Mcluster.morpho.data.ln.pca.subset.G7$G, labels=F)
axis(1, at=seq(1,Mcluster.morpho.data.ln.pca.subset.G7$G, 3), labels=T, cex.axis=1.5)

hist(Mcluster.morpho.data.ln.pca.subset.G7$classification[morpho.data$Taxon=="G. scandens"], 
breaks=seq(0.75, Mcluster.morpho.data.ln.pca.subset.G7$G+0.25, 0.5),
main=expression(italic("G. scandens")), xlab="", ylab="",
cex.lab=1.5, cex.axis=1.5, col=morpho.groups.colors.hist, xaxt="n") 
axis(1, at=1:Mcluster.morpho.data.ln.pca.subset.G7$G, labels=F)
axis(1, at=seq(1,Mcluster.morpho.data.ln.pca.subset.G7$G, 3), labels=T, cex.axis=1.5)

mtext('Morphological group', side = 1, outer = TRUE, line = 2, cex = 1.3)
mtext('Specimens', side = 2, outer = TRUE, line = 2, cex = 1.3)


################################################################################################################
# 2.7)Examine co-occurrence of morphological groups (as defined by the best Mclust model 
# with seven groups) across islands.

BestModel.G7.cooccurrence.matrix <- matrix(NA, ncol=Mcluster.morpho.data.ln.pca.subset.G7$G, nrow=Mcluster.morpho.data.ln.pca.subset.G7$G)

for(i in 1:Mcluster.morpho.data.ln.pca.subset.G7$G)
{
	BestModel.G7.Islands.FocalMorphogroup <- unique(morpho.data$Island[Mcluster.morpho.data.ln.pca.subset.G7$classification==i])
	BestModel.G7.Islands.FocalMorphogroup <- as.vector(BestModel.G7.Islands.FocalMorphogroup)

	index.cooccurrence.with.FocalMorphogroup <- !is.na(match(morpho.data$Island, BestModel.G7.Islands.FocalMorphogroup))
	morphogroups.cooccurring.with.FocalMorphogroup <- Mcluster.morpho.data.ln.pca.subset.G7$classification[index.cooccurrence.with.FocalMorphogroup]
	unique.morphogroups.cooccurring.with.FocalMorphogroup <- sort(unique(morphogroups.cooccurring.with.FocalMorphogroup))
	Islands.of.cooccurrence.with.FocalMorphogroup <- morpho.data$Island[index.cooccurrence.with.FocalMorphogroup]
	Islands.of.cooccurrence.with.FocalMorphogroup <- as.vector(Islands.of.cooccurrence.with.FocalMorphogroup)

	BestModel.G7.cooccurrence.matrix[unique.morphogroups.cooccurring.with.FocalMorphogroup,i] <- 
	rowSums(table(morphogroups.cooccurring.with.FocalMorphogroup, Islands.of.cooccurrence.with.FocalMorphogroup)>0)
}

BestModel.G7.cooccurrence.matrix[is.na(BestModel.G7.cooccurrence.matrix)] <- 0

#examine results
BestModel.G7.cooccurrence.matrix


################################################################################################################
# 2.8)Examine co-occurrence of morphological groups (as defined by the best Mclust model 
# with eigth groups) across islands.

BestModel.G8.cooccurrence.matrix <- matrix(NA, ncol=Mcluster.morpho.data.ln.pca.subset$G, nrow=Mcluster.morpho.data.ln.pca.subset$G)

for(i in 1:Mcluster.morpho.data.ln.pca.subset$G)
{
	BestModel.G8.Islands.FocalMorphogroup <- unique(morpho.data$Island[Mcluster.morpho.data.ln.pca.subset$classification==i])
	BestModel.G8.Islands.FocalMorphogroup <- as.vector(BestModel.G8.Islands.FocalMorphogroup)

	index.cooccurrence.with.FocalMorphogroup <- !is.na(match(morpho.data$Island, BestModel.G8.Islands.FocalMorphogroup))
	morphogroups.cooccurring.with.FocalMorphogroup <- Mcluster.morpho.data.ln.pca.subset$classification[index.cooccurrence.with.FocalMorphogroup]
	unique.morphogroups.cooccurring.with.FocalMorphogroup <- sort(unique(morphogroups.cooccurring.with.FocalMorphogroup))
	Islands.of.cooccurrence.with.FocalMorphogroup <- morpho.data$Island[index.cooccurrence.with.FocalMorphogroup]
	Islands.of.cooccurrence.with.FocalMorphogroup <- as.vector(Islands.of.cooccurrence.with.FocalMorphogroup)

	BestModel.G8.cooccurrence.matrix[unique.morphogroups.cooccurring.with.FocalMorphogroup,i] <- 
	rowSums(table(morphogroups.cooccurring.with.FocalMorphogroup, Islands.of.cooccurrence.with.FocalMorphogroup)>0)
}


BestModel.G8.cooccurrence.matrix[is.na(BestModel.G8.cooccurrence.matrix)] <- 0

#examine results
BestModel.G8.cooccurrence.matrix