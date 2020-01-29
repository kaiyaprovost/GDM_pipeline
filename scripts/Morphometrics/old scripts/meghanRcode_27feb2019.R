# Exploratory Data Analysis
#2/13/2019, Meghan Forcellati, AMNH Ornithology 
#modified 2/20/2019

#PCA BY THE END OF TODAY
"
NO PCA'S FOR EACH DESERT
PCA ON EVERYTHING (ALL SPECIES CONGLOMERATED INTO 2, LETS YOU COMPARE
SPECIES AND DESERT)

MOST INTEREST: WITHIN SPECIES DESERT COMPARISON 

SPLIT OUT EVERY SPECIES, COMPARE WITHIN DESERT. 
EACH DESERT WILL BE COMPARABLE. ON SAME PC AXES. 
MOST AMOUNT OF ABILITY TO TELL APART THE DESERTS IF YOU DO EACH 

TESTING SIGNIFICANCE - TTEST ON THEM
LINEAR MODEL - WHICH DESERT? DOES IT PREDICT PC SCORE? 
TTEST - IS PC1 SIGNIFICANTLY DIFFERENT? 
  COLOR IN WHICH INDIVIDUAL IS WHICH DESERT 

LARGER NUMBER = MORE IMPORTANT TO PC ANALYSIS
  PC1 IS USUALLY SIZE
  PC2 IS USUALLY SHAPE
"

#library imports
library("corrplot")
library("ggplot2")
library("doBy")
library("plyr")
library("dplyr")
library("tibble")
#Set working directory

#Import the file
prefix="Master_Spreadsheet_Morphology_round2_"
date="27feb2019-2"
suffix=".csv"
bird_data=paste(prefix, date, suffix, sep="")
rawdata<-read.csv(bird_data)
View(rawdata)

#Separate out AMNH specimens
prefixes=substr(rawdata$CATALOG.NUMBER, 0,  4)
wip_data <- tibble::add_column(rawdata, prefixes)
refined_data <- dplyr::filter(rawdata, prefixes!="DMNH", prefixes != "AMNH" )
refined_data$prefixes <- NULL
View(refined_data)


#get number averages per individual
dplyr::tbl_df(refined_data)
View(refined_data)
summary(as.numeric(refined_data$TAIL))
agg1 <- refined_data %>% group_by(CATALOG.NUMBER) %>% summarise(BILL.LENGTH=mean(BILL.LENGTH), BILL.HEIGHT=mean(BILL.LENGTH),
                        BILL.WIDTH=mean(BILL.LENGTH), TARSUS.LENGTH=mean(TARSUS.LENGTH), WING.LENGTH.TO.PRIMARIES=mean(WING.LENGTH.TO.PRIMARIES), WING.LENGTH.TO.SECONDARIES=mean(WING.LENGTH.TO.SECONDARIES), TAIL=mean(as.numeric(TAIL)))
View(agg1)

#make string columns stuck together by data frame
agg2 <- refined_data %>% group_by(CATALOG.NUMBER, GENUS, SPP, WHICH.SIDE.OF.CFB)
  agg2<- summarise(agg2)
View(agg2)

#combine.
good_data<-merge.data.frame(agg1, agg2)
View(good_data)

good_data <- dplyr::filter(good_data, WHICH.SIDE.OF.CFB!="UNCLEAR")
View(good_data)
#Various graphs which may be useful
ggplot(good_data, aes(x=TARSUS.LENGTH, color=SPP))+geom_histogram()
ggplot(good_data, aes(x=SPP, y=BILL.LENGTH, color=SPP))+geom_boxplot()
ggplot(good_data, aes(x=SPP, y=TARSUS.LENGTH, color=SPP)) +geom_boxplot()       
ggplot(good_data, aes(x=WHICH.SIDE.OF.CFB, y=TARSUS.LENGTH, color=SPP)) +geom_boxplot()

#Separating specific species/genus for the data set
fusca<-good_data[good_data$SPP=="FUSCA",]
ggplot(fusca, aes(x=WHICH.SIDE.OF.CFB, y=WING.LENGTH.TO., color=WHICH.SIDE.OF.CFB)) +geom_boxplot()

bellii<-good_data[good_data$SPP=="BELLII",]
ggplot(bellii, aes(x=WHICH.SIDE.OF.CFB, y=TARSUS.LENGTH, color=WHICH.SIDE.OF.CFB)) +geom_boxplot()
ggplot(good_data, aes(x=WHICH.SIDE.OF.CFB, color=WHICH.SIDE.OF.CFB))+geom_bar()
ggplot(good_data, aes(fill=WHICH.SIDE.OF.CFB, x=SPP, y=WHICH.SIDE.OF.CFB)) + 
  geom_bar(position="dodge", stat='identity')



#finding the sample size per desert per bird
#SORT DATA by counts
df <- good_data %>%
  group_by(WHICH.SIDE.OF.CFB, SPP) %>%
  summarise(counts = n()) 
view(df)


# plot everything
ggplot(df, aes(SPP, counts)) +   
  geom_bar(aes(fill = WHICH.SIDE.OF.CFB), position = "dodge", stat="identity")

#DO PCA
#SEPARATE EACH SPECIES
print(good_data$SPP)
#EVERYTHING TO NUMERIC
good_data$BILL.HEIGHT<- as.numeric(good_data$BILL.HEIGHT)
good_data$BILL.WIDTH<- as.numeric(good_data$BILL.WIDTH)
good_data$BILL.LENGTH<- as.numeric(good_data$BILL.LENGTH)
good_data$TARSUS.LENGTH<- as.numeric(good_data$BILL.HEIGHT)
good_data$WING.LENGTH.TO.PRIMARIES<- as.numeric(good_data$WING.LENGTH.TO.PRIMARIES)
good_data$WING.LENGTH.TO.SECONDARIES<- as.numeric(good_data$WING.LENGTH.TO.SECONDARIES)
good_data$WING.TAIL<- as.numeric(good_data$TAIL)
View(good_data)

#bellii
#data
bellii<-good_data[good_data$SPP=="BELLII", good_data$WHICH.SIDE.OF.CFB!="UNCLEAR"]
bellii.WHICH.SIDE.OF.CFB<-bellii[,11]
View(bellii)

#pca
scaled_bellii.pca<-prcomp(bellii[2:8],
                                  center=TRUE,
                                  scale.=TRUE)
View(scaled_bellii.pca)

print(scaled_bellii.pca)
summary(scaled_bellii.pca)
plot(scaled_bellii.pca, type="lines")

bellii.results=scaled_bellii.pca$x
View(bellii.results)
colnames(bellii.results)
bellii.PC1 = bellii.results[,1]
bellii.PC2 = bellii.results[,2]
bellii.PC3 = bellii.results[,3]
plot(bellii.PC1, bellii.PC2)
ggplot(bellii, aes(bellii.PC1, bellii.PC2, color=WHICH.SIDE.OF.CFB)) +geom_point()
ggplot(bellii, aes(bellii.PC1, bellii.PC2, color=WHICH.SIDE.OF.CFB)) +geom_boxplot()

View(bellii.results)

#ttest
bellii_ttest <- cbind(bellii, bellii.results)
View(bellii_ttest)

t.test(bellii_ttest$PC1 ~ bellii_ttest$WHICH.SIDE.OF.CFB)
sd(bellii_ttest$PC1[bellii_ttest$WHICH.SIDE.OF.CFB=="SONORAN"])
sd(bellii_ttest$PC1[bellii_ttest$WHICH.SIDE.OF.CFB=="CHIHUAHUAN"])
t.test(bellii_ttest$PC2~ bellii_ttest$WHICH.SIDE.OF.CFB)
sd(bellii_ttest$PC2[bellii_ttest$WHICH.SIDE.OF.CFB=="SONORAN"])
sd(bellii_ttest$PC2[bellii_ttest$WHICH.SIDE.OF.CFB=="CHIHUAHUAN"])
t.test(bellii_ttest$PC3 ~ bellii_ttest$WHICH.SIDE.OF.CFB)
sd(bellii_ttest$PC3[bellii_ttest$WHICH.SIDE.OF.CFB=="SONORAN"])
sd(bellii_ttest$PC3[bellii_ttest$WHICH.SIDE.OF.CFB=="CHIHUAHUAN"])

#BILINEATA
bilineata<-good_data[good_data$SPP=="BILINEATA",]

bilineata.WHICH.SIDE.OF.CFB<-bilineata[,11]
View(bilineata)

#pca
scaled_bilineata.pca<-prcomp(bilineata[2:8],
                          center=TRUE,
                          scale.=TRUE)

print(scaled_bilineata.pca)
summary(scaled_bilineata.pca)
plot(scaled_bilineata.pca, type="lines")

bilineata.results=scaled_bilineata.pca$x
colnames(bilineata.results)
bilineata.PC1 = bilineata.results[,1]
bilineata.PC2 = bilineata.results[,2]
bilineata.PC3 = bilineata.results[,3]
plot(bilineata.PC1, bilineata.PC2)
ggplot(bilineata, aes(bilineata.PC1, bilineata.PC2, color=WHICH.SIDE.OF.CFB)) +geom_point()
ggplot(bilineata, aes(bilineata.PC1, bilineata.PC2, color=WHICH.SIDE.OF.CFB)) +geom_boxplot()

#ttest
bilineata_ttest <- cbind(bilineata, bilineata.results)
View(bilineata_ttest)

t.test(bilineata_ttest$PC1 ~ bilineata_ttest$WHICH.SIDE.OF.CFB)
sd(bilineata_ttest$PC1[bilineata_ttest$WHICH.SIDE.OF.CFB=="SONORAN"])
sd(bilineata_ttest$PC1[bilineata_ttest$WHICH.SIDE.OF.CFB=="CHIHUAHUAN"])
t.test(bilineata_ttest$PC2~ bilineata_ttest$WHICH.SIDE.OF.CFB)
sd(bilineata_ttest$PC2[bilineata_ttest$WHICH.SIDE.OF.CFB=="SONORAN"])
sd(bilineata_ttest$PC2[bilineata_ttest$WHICH.SIDE.OF.CFB=="CHIHUAHUAN"])
t.test(bilineata_ttest$PC3 ~ bilineata_ttest$WHICH.SIDE.OF.CFB)
sd(bilineata_ttest$PC3[bilineata_ttest$WHICH.SIDE.OF.CFB=="SONORAN"])
sd(bilineata_ttest$PC3[bilineata_ttest$WHICH.SIDE.OF.CFB=="CHIHUAHUAN"])

"
#BRUNNEICAPILLUS
brunneicapillus<-good_data[good_data$SPP=='BRUNNEICAPILLUS',]
brunneicapillus.WHICH.SIDE.OF.CFB<-bilineata[,11]
View(brunneicapillus)

#pca
scaled_brunneicapillus.pca<-prcomp(brunneicapillus[2:8],
                             center=TRUE,
                             scale.=TRUE)

print(scaled_brunneicapillus.pca)
summary(scaled_brunneicapillus.pca)
plot(scaled_brunneicapillus.pca, type='lines')

brunneicapillus.results=scaled_brunneicapillus.pca$x
colnames(brunneicapillus.results)
brunneicapillus.PC1 = brunneicapillus.results[,1]
brunneicapillus.PC2 = brunneicapillus.results[,2]
brunneicapillus.PC3 = brunneicapillus.results[,3]
plot(brunneicapillus.PC1, brunneicapillus.PC2)
ggplot(brunneicapillus, aes(brunneicapillus.PC1, brunneicapillus.PC2, color=WHICH.SIDE.OF.CFB)) +geom_point()
ggplot(brunneicapillus, aes(brunneicapillus.PC1, brunneicapillus.PC2, color=WHICH.SIDE.OF.CFB)) +geom_boxplot()

WARNING: NO BRUNEICAPILLUS IN DATASET 
"

#SINUATUS
sinuatus<-good_data[good_data$SPP=="SINUATUS",]
sinuatus.WHICH.SIDE.OF.CFB<-sinuatus[,11]
View(sinuatus)

#pca
scaled_sinuatus.pca<-prcomp(sinuatus[2:8],
                             center=TRUE,
                             scale.=TRUE)

print(scaled_sinuatus.pca)
summary(scaled_sinuatus.pca)
plot(scaled_sinuatus.pca, type="lines")

sinuatus.results=scaled_sinuatus.pca$x
colnames(sinuatus.results)
sinuatus.PC1 = sinuatus.results[,1]
sinuatus.PC2 = sinuatus.results[,2]
sinuatus.PC3 = sinuatus.results[,3]
plot(sinuatus.PC1, sinuatus.PC2)
ggplot(sinuatus, aes(sinuatus.PC1, sinuatus.PC2, color=WHICH.SIDE.OF.CFB)) +geom_point()
ggplot(sinuatus, aes(sinuatus.PC1, sinuatus.PC2, color=WHICH.SIDE.OF.CFB)) +geom_boxplot()

#ttest
sinuatus_ttest <- cbind(sinuatus, sinuatus.results)
View(sinuatus_ttest)

t.test(sinuatus_ttest$PC1 ~ sinuatus_ttest$WHICH.SIDE.OF.CFB)
sd(sinuatus_ttest$PC1[sinuatus_ttest$WHICH.SIDE.OF.CFB=="SONORAN"])
sd(sinuatus_ttest$PC1[sinuatus_ttest$WHICH.SIDE.OF.CFB=="CHIHUAHUAN"])
t.test(sinuatus_ttest$PC2~ sinuatus_ttest$WHICH.SIDE.OF.CFB)
sd(sinuatus_ttest$PC2[sinuatus_ttest$WHICH.SIDE.OF.CFB=="SONORAN"])
sd(sinuatus_ttest$PC2[sinuatus_ttest$WHICH.SIDE.OF.CFB=="CHIHUAHUAN"])
t.test(sinuatus_ttest$PC3 ~ sinuatus_ttest$WHICH.SIDE.OF.CFB)
sd(sinuatus_ttest$PC3[sinuatus_ttest$WHICH.SIDE.OF.CFB=="SONORAN"])
sd(sinuatus_ttest$PC3[sinuatus_ttest$WHICH.SIDE.OF.CFB=="CHIHUAHUAN"])

#CRISSALE
crissale<-good_data[good_data$SPP=="CRISSALE",]
crissale.WHICH.SIDE.OF.CFB<-crissale[,11]
View(crissale)

#pca
scaled_crissale.pca<-prcomp(crissale[2:8],
                            center=TRUE,
                            scale.=TRUE)

print(scaled_crissale.pca)
summary(scaled_crissale.pca)
plot(scaled_crissale.pca, type="lines")

crissale.results=scaled_crissale.pca$x
colnames(crissale.results)
crissale.PC1 = crissale.results[,1]
crissale.PC2 = crissale.results[,2]
crissale.PC3 = crissale.results[,3]
plot(crissale.PC1, crissale.PC2)
ggplot(crissale, aes(crissale.PC1, crissale.PC2, color=WHICH.SIDE.OF.CFB)) +geom_point()
ggplot(crissale, aes(crissale.PC1, crissale.PC2, color=WHICH.SIDE.OF.CFB)) +geom_boxplot()

#ttest
crissale_ttest <- cbind(crissale, crissale.results)
View(crissale_ttest)

t.test(crissale_ttest$PC1 ~ crissale_ttest$WHICH.SIDE.OF.CFB)
sd(crissale_ttest$PC1[crissale_ttest$WHICH.SIDE.OF.CFB=="SONORAN"])
sd(crissale_ttest$PC1[crissale_ttest$WHICH.SIDE.OF.CFB=="CHIHUAHUAN"])
t.test(crissale_ttest$PC2~ crissale_ttest$WHICH.SIDE.OF.CFB)
sd(crissale_ttest$PC2[crissale_ttest$WHICH.SIDE.OF.CFB=="SONORAN"])
sd(crissale_ttest$PC2[crissale_ttest$WHICH.SIDE.OF.CFB=="CHIHUAHUAN"])
t.test(crissale_ttest$PC3 ~ crissale_ttest$WHICH.SIDE.OF.CFB)
sd(crissale_ttest$PC3[crissale_ttest$WHICH.SIDE.OF.CFB=="SONORAN"])
sd(crissale_ttest$PC3[crissale_ttest$WHICH.SIDE.OF.CFB=="CHIHUAHUAN"])


#CURVIROSTRE
curvirostre<-good_data[good_data$SPP=="CURVIROSTRE",]
curvirostre.WHICH.SIDE.OF.CFB<-curvirostre[,11]
View(curvirostre)

#pca
scaled_curvirostre.pca<-prcomp(curvirostre[2:8],
                            center=TRUE,
                            scale.=TRUE)

print(scaled_curvirostre.pca)
summary(scaled_curvirostre.pca)
plot(scaled_curvirostre.pca, type="lines")

curvirostre.results=scaled_curvirostre.pca$x
colnames(curvirostre.results)
curvirostre.PC1 = curvirostre.results[,1]
curvirostre.PC2 = curvirostre.results[,2]
curvirostre.PC3 = curvirostre.results[,3]
plot(curvirostre.PC1, curvirostre.PC2)
ggplot(curvirostre, aes(curvirostre.PC1, curvirostre.PC2, color=WHICH.SIDE.OF.CFB)) +geom_point()
ggplot(curvirostre, aes(curvirostre.PC1, curvirostre.PC2, color=WHICH.SIDE.OF.CFB)) +geom_boxplot()

#ttest
curvirostre_ttest <- cbind(curvirostre, curvirostre.results)
View(curvirostre_ttest)

t.test(curvirostre_ttest$PC1 ~ curvirostre_ttest$WHICH.SIDE.OF.CFB)
sd(curvirostre_ttest$PC1[curvirostre_ttest$WHICH.SIDE.OF.CFB=="SONORAN"])
sd(curvirostre_ttest$PC1[curvirostre_ttest$WHICH.SIDE.OF.CFB=="CHIHUAHUAN"])
t.test(curvirostre_ttest$PC2~ curvirostre_ttest$WHICH.SIDE.OF.CFB)
sd(curvirostre_ttest$PC2[curvirostre_ttest$WHICH.SIDE.OF.CFB=="SONORAN"])
sd(curvirostre_ttest$PC2[curvirostre_ttest$WHICH.SIDE.OF.CFB=="CHIHUAHUAN"])
t.test(curvirostre_ttest$PC3 ~ curvirostre_ttest$WHICH.SIDE.OF.CFB)
sd(curvirostre_ttest$PC3[curvirostre_ttest$WHICH.SIDE.OF.CFB=="SONORAN"])
sd(curvirostre_ttest$PC3[curvirostre_ttest$WHICH.SIDE.OF.CFB=="CHIHUAHUAN"])


#FLAVICEPS
flaviceps<-good_data[good_data$SPP=="FLAVICEPS",]
flaviceps.WHICH.SIDE.OF.CFB<-flaviceps[,11]
View(flaviceps)

#pca
scaled_flaviceps.pca<-prcomp(flaviceps[2:8],
                               center=TRUE,
                               scale.=TRUE)

print(scaled_flaviceps.pca)
summary(scaled_flaviceps.pca)
plot(scaled_flaviceps.pca, type="lines")

flaviceps.results=scaled_flaviceps.pca$x
colnames(flaviceps.results)
flaviceps.PC1 = flaviceps.results[,1]
flaviceps.PC2 = flaviceps.results[,2]
flaviceps.PC3 = flaviceps.results[,3]
plot(flaviceps.PC1, flaviceps.PC2)
ggplot(flaviceps, aes(flaviceps.PC1, flaviceps.PC2, color=WHICH.SIDE.OF.CFB)) +geom_point()
ggplot(flaviceps, aes(flaviceps.PC1, flaviceps.PC2, color=WHICH.SIDE.OF.CFB)) +geom_boxplot()

#ttest
flaviceps_ttest <- cbind(flaviceps, flaviceps.results)
View(flaviceps_ttest)

t.test(flaviceps_ttest$PC1 ~ flaviceps_ttest$WHICH.SIDE.OF.CFB)
sd(flaviceps_ttest$PC1[flaviceps_ttest$WHICH.SIDE.OF.CFB=="SONORAN"])
sd(flaviceps_ttest$PC1[flaviceps_ttest$WHICH.SIDE.OF.CFB=="CHIHUAHUAN"])
t.test(flaviceps_ttest$PC2~ flaviceps_ttest$WHICH.SIDE.OF.CFB)
sd(flaviceps_ttest$PC2[flaviceps_ttest$WHICH.SIDE.OF.CFB=="SONORAN"])
sd(flaviceps_ttest$PC2[flaviceps_ttest$WHICH.SIDE.OF.CFB=="CHIHUAHUAN"])
t.test(flaviceps_ttest$PC3 ~ flaviceps_ttest$WHICH.SIDE.OF.CFB)
sd(flaviceps_ttest$PC3[flaviceps_ttest$WHICH.SIDE.OF.CFB=="SONORAN"])
sd(flaviceps_ttest$PC3[flaviceps_ttest$WHICH.SIDE.OF.CFB=="CHIHUAHUAN"])


#FUSCA
fusca<-good_data[good_data$SPP=="FUSCA",]
fusca.WHICH.SIDE.OF.CFB<-fusca[,11]
View(fusca)

#pca
scaled_fusca.pca<-prcomp(fusca[2:8],
                             center=TRUE,
                             scale.=TRUE)

print(scaled_fusca.pca)
summary(scaled_fusca.pca)
plot(scaled_fusca.pca, type="lines")

fusca.results=scaled_fusca.pca$x
colnames(fusca.results)
fusca.PC1 = fusca.results[,1]
fusca.PC2 = fusca.results[,2]
fusca.PC3 = fusca.results[,3]
plot(fusca.PC1, fusca.PC2)
ggplot(fusca, aes(fusca.PC1, fusca.PC2, color=WHICH.SIDE.OF.CFB)) +geom_point()
ggplot(fusca, aes(fusca.PC1, fusca.PC2, color=WHICH.SIDE.OF.CFB)) +geom_boxplot()

#ttest
fusca_ttest <- cbind(fusca, fusca.results)
View(fusca_ttest)

t.test(fusca_ttest$PC1 ~ fusca_ttest$WHICH.SIDE.OF.CFB)
sd(fusca_ttest$PC1[fusca_ttest$WHICH.SIDE.OF.CFB=="SONORAN"])
sd(fusca_ttest$PC1[fusca_ttest$WHICH.SIDE.OF.CFB=="CHIHUAHUAN"])
t.test(fusca_ttest$PC2~ fusca_ttest$WHICH.SIDE.OF.CFB)
sd(fusca_ttest$PC2[fusca_ttest$WHICH.SIDE.OF.CFB=="SONORAN"])
sd(fusca_ttest$PC2[fusca_ttest$WHICH.SIDE.OF.CFB=="CHIHUAHUAN"])
t.test(fusca_ttest$PC3 ~ fusca_ttest$WHICH.SIDE.OF.CFB)
sd(fusca_ttest$PC3[fusca_ttest$WHICH.SIDE.OF.CFB=="SONORAN"])
sd(fusca_ttest$PC3[fusca_ttest$WHICH.SIDE.OF.CFB=="CHIHUAHUAN"])


#MELANURA
melanura<-good_data[good_data$SPP=="MELANURA",]
melanura.WHICH.SIDE.OF.CFB<-melanura[,11]
View(melanura)

#pca
scaled_melanura.pca<-prcomp(melanura[2:8],
                             center=TRUE,
                             scale.=TRUE)

print(scaled_melanura.pca)
summary(scaled_melanura.pca)
plot(scaled_melanura.pca, type="lines")

melanura.results=scaled_melanura.pca$x
colnames(melanura.results)
melanura.PC1 = melanura.results[,1]
melanura.PC2 = melanura.results[,2]
melanura.PC3 = melanura.results[,3]
plot(melanura.PC1, melanura.PC2)
ggplot(melanura, aes(melanura.PC1, melanura.PC2, color=WHICH.SIDE.OF.CFB)) +geom_point()
ggplot(melanura, aes(melanura.PC1, melanura.PC2, color=WHICH.SIDE.OF.CFB)) +geom_boxplot()

#ttest
melanura_ttest <- cbind(melanura, melanura.results)
View(melanura_ttest)

t.test(melanura_ttest$PC1 ~ melanura_ttest$WHICH.SIDE.OF.CFB)
sd(melanura_ttest$PC1[melanura_ttest$WHICH.SIDE.OF.CFB=="SONORAN"])
sd(melanura_ttest$PC1[melanura_ttest$WHICH.SIDE.OF.CFB=="CHIHUAHUAN"])
t.test(melanura_ttest$PC2~ melanura_ttest$WHICH.SIDE.OF.CFB)
sd(melanura_ttest$PC2[melanura_ttest$WHICH.SIDE.OF.CFB=="SONORAN"])
sd(melanura_ttest$PC2[melanura_ttest$WHICH.SIDE.OF.CFB=="CHIHUAHUAN"])
t.test(melanura_ttest$PC3 ~ melanura_ttest$WHICH.SIDE.OF.CFB)
sd(melanura_ttest$PC3[melanura_ttest$WHICH.SIDE.OF.CFB=="SONORAN"])
sd(melanura_ttest$PC3[melanura_ttest$WHICH.SIDE.OF.CFB=="CHIHUAHUAN"])


#NITENS
nitens<-good_data[good_data$SPP=="NITENS",]
nitens.WHICH.SIDE.OF.CFB<-nitens[,11]
View(nitens)

#pca
scaled_nitens.pca<-prcomp(nitens[2:8],
                             center=TRUE,
                             scale.=TRUE)

print(scaled_nitens.pca)
summary(scaled_nitens.pca)
plot(scaled_nitens.pca, type="lines")

nitens.results=scaled_nitens.pca$x
colnames(nitens.results)
nitens.PC1 = nitens.results[,1]
nitens.PC2 = nitens.results[,2]
nitens.PC3 = nitens.results[,3]
plot(nitens.PC1, nitens.PC2)
ggplot(nitens, aes(nitens.PC1, nitens.PC2, color=WHICH.SIDE.OF.CFB)) +geom_point()
ggplot(nitens, aes(nitens.PC1, nitens.PC2, color=WHICH.SIDE.OF.CFB)) +geom_boxplot()

#ttest
nitens_ttest <- cbind(nitens, nitens.results)
View(nitens_ttest)

t.test(nitens_ttest$PC1 ~ nitens_ttest$WHICH.SIDE.OF.CFB)
sd(nitens_ttest$PC1[nitens_ttest$WHICH.SIDE.OF.CFB=="SONORAN"])
sd(nitens_ttest$PC1[nitens_ttest$WHICH.SIDE.OF.CFB=="CHIHUAHUAN"])
t.test(nitens_ttest$PC2~ nitens_ttest$WHICH.SIDE.OF.CFB)
sd(nitens_ttest$PC2[nitens_ttest$WHICH.SIDE.OF.CFB=="SONORAN"])
sd(nitens_ttest$PC2[nitens_ttest$WHICH.SIDE.OF.CFB=="CHIHUAHUAN"])
t.test(nitens_ttest$PC3 ~ nitens_ttest$WHICH.SIDE.OF.CFB)
sd(nitens_ttest$PC3[nitens_ttest$WHICH.SIDE.OF.CFB=="SONORAN"])
sd(nitens_ttest$PC3[nitens_ttest$WHICH.SIDE.OF.CFB=="CHIHUAHUAN"])


#ALL ttests & data tables



