# Exploratory Data Analysis
#2/13/2019, Meghan Forcellati, AMNH Ornithology 
#modified 2/20/2019

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
date="20feb2019-2"
suffix=".csv"
bird_data=paste(prefix, date, suffix, sep="")
rawdata<-read.csv(bird_data)
View(rawdata)

#Separate out AMNH specimens
prefixes=substr(rawdata$CATALOG.NUMBER, 0,  4)
wip_data <- tibble::add_column(rawdata, prefixes)
refined_data <- dplyr::filter(rawdata, prefixes!="DMNH" )
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

originalpca = prcomp(na.omit(good_data[,2:8]), scale=TRUE, center=TRUE)
plot(originalpca)
print(originalpca)
#rotation = correlation. bill length negative correlated w pc1, pc2, pc3. everything flipped by x-axis.
#width: normal x, flipped y, normal z
#rotated but not saved
#center and scale: if measurements on different scales, ie. apartment sizes in ny state. number floors vs sqftge
#num floors is 1-6, sqftge = 800-8000
#800-8000 would always be biggest axis no matter what. prcomp can convert things into same center/scale, however
print(scaledpca)
#z-scored it, basically. slightly same numbers. 
pca$center
#avg values for centering
pca$scale
#scaling factor
pca$rotation
pca$sdev
pca$x

originalpca
#obviously, length is longest for pc1
scaledpca
#now it lines up width/height more because there's more variation there. 
results = scaledpca$x
summary(scaledpca)
#Proportion of Variance tells you how successful your pca is on those variables 
#if variables dont explain data, might not be worth analyzing
#pc1 is over 65% data. pc2 32%. 
#cumulative proportion = % explained over time
#plot(scaledpca) plots the variance
plot(scaledpca, type="l")
#high pc1, moderate pc2. want higher dropoff between 1 and 2. pc1 = 99.1 is awesome because it collapses
#everything down. 

#now that I've taken my variables, rotated them, extracted the max out of them
#use plot function for pc1 vs pc2 and see what that tells them

colnames(results) #if $ doesn't work
PC1 = results[,1]
PC2 = results[,2]
PC3 = results[,3]
plot(PC1, PC2) #probably has to do with species. 

par(mfrow=c(1, 2))
plot(bill_height, bill_length) #PC1 describes bill height variation. PC2 bill length variation. Scales are different
#graphs look a little different but 

plot(PC1, PC3)
#garbage bc pc1 doesnt describe anything

  
