#Meghan Forcellati
#2/27/2019
#Temporary ttest for tail measurements to ID significance where needed (remeasurement process)

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

print(unique(good_data$SPP))

#FLAVICEPS
flaviceps<-good_data[good_data$SPP=="FLAVICEPS",]
ggplot(flaviceps, aes(flaviceps$TAIL, flaviceps$WHICH.SIDE.OF.CFB, color=WHICH.SIDE.OF.CFB)) +geom_point()
t.test(flaviceps$TAIL ~ flaviceps$WHICH.SIDE.OF.CFB)

#CRISSALE
crissale<-good_data[good_data$SPP=="CRISSALE",]
ggplot(crissale, aes(crissale$TAIL, crissale$WHICH.SIDE.OF.CFB, color=WHICH.SIDE.OF.CFB)) +geom_point()
t.test(crissale$TAIL ~ crissale$WHICH.SIDE.OF.CFB)

#MELANURA
melanura<-good_data[good_data$SPP=="MELANURA",]
ggplot(melanura, aes(melanura$TAIL, melanura$WHICH.SIDE.OF.CFB, color=WHICH.SIDE.OF.CFB)) +geom_point()
t.test(melanura$TAIL ~ melanura$WHICH.SIDE.OF.CFB)


#NITENS
nitens<-good_data[good_data$SPP=="NITENS",]
ggplot(nitens, aes(nitens$TAIL, nitens$WHICH.SIDE.OF.CFB, color=WHICH.SIDE.OF.CFB)) +geom_point()
t.test(nitens$TAIL ~ nitens$WHICH.SIDE.OF.CFB)

#BELLII
bellii<-good_data[good_data$SPP=="BELLII",]
ggplot(bellii, aes(bellii$TAIL, bellii$WHICH.SIDE.OF.CFB, color=WHICH.SIDE.OF.CFB)) +geom_point()
t.test(bellii$TAIL ~ bellii$WHICH.SIDE.OF.CFB)

#FUSCA
fusca<-good_data[good_data$SPP=="FUSCA",]
ggplot(fusca, aes(fusca$TAIL, fusca$WHICH.SIDE.OF.CFB, color=WHICH.SIDE.OF.CFB)) +geom_point()
t.test(fusca$TAIL ~ fusca$WHICH.SIDE.OF.CFB)

#SINUATUS
sinuatus<-good_data[good_data$SPP=="SINUATUS",]
ggplot(sinuatus, aes(sinuatus$TAIL, sinuatus$WHICH.SIDE.OF.CFB, color=WHICH.SIDE.OF.CFB)) +geom_point()
t.test(sinuatus$TAIL ~ sinuatus$WHICH.SIDE.OF.CFB)

#BILINEATA
bilineata<-good_data[good_data$SPP=="BILINEATA",]
ggplot(bilineata, aes(bilineata$TAIL, bilineata$WHICH.SIDE.OF.CFB, color=WHICH.SIDE.OF.CFB)) +geom_point()
t.test(bilineata$TAIL ~ bilineata$WHICH.SIDE.OF.CFB)

#CURVIROSTRE
curvirostre<-good_data[good_data$SPP=="CURVIROSTRE",]  
ggplot(curvirostre, aes(curvirostre$TAIL, curvirostre$WHICH.SIDE.OF.CFB, color=WHICH.SIDE.OF.CFB)) +geom_point()
t.test(curvirostre$TAIL ~ curvirostre$WHICH.SIDE.OF.CFB)