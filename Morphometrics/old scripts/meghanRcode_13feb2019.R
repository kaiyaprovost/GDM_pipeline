# Exploratory Data Analysis
#2/13/2019, Meghan Forcellati, AMNH Ornithology 

library("corrplot")
library("ggplot2")
library("doBy")
library("dplyr")
library("tibble")
#set working directory
prefix="Master_Spreadsheet_Morphology_round2_"
date="13feb2019"
suffix=".csv"
bird_data=paste(prefix, date, suffix, sep="")

rawdata<-read.csv(bird_data)

#test
View(rawdata)

"#Looking at just bill length, bill width by catalogue number. 
agg = aggregate(cbind(BILL.LENGTH, BILL.WIDTH) ~ CATALOG.NUMBER, 
                data=rawdata,
                FUN=mean, 
                na.action=na.pass) #make columns stuck together by data frame

View(agg)"

#Let's avoid the AMNH specimens for now because some of them are relatively sketch.
"
Problem: hard to
1. Create column with all prefixes
2. Filter our any rows with value of AMNH
3. remove the extra column"
prefixes=substr(rawdata$CATALOG.NUMBER, 0,  4)
wip_data <- tibble::add_column(rawdata, prefixes)
refined_data <- dplyr::filter(rawdata, prefixes!="AMNH", prefixes!="DMNH", )
refined_data$prefixes <- NULL
View(refined_data)

#take avg's for numeric values, any values you want
agg=refined_data %>%
  group_by(~CATALOG.NUMBER, 
           ~BILL.LENGTH, ~BILL.WIDTH, ~BILL.HEIGHT, ~TARSUS.LENGTH, 
           ~WING.LENGTH.TO.PRIMARIES, ~WING.LENGTH.TO.SECONDARIES, 
           ~TAIL)
summarise_()
View(agg)

#make string columns stuck together by data frame
agg2 <-refined_data %>% 
  group_by_(~CATALOG.NUMBER, ~GENUS, ~SPP, ~WHICH.SIDE.OF.CFB) %>%
  summarize_()
View(agg2)

#combine agg and agg2
good_data<-merge.data.frame(agg, agg2)
View(good_data)


#graph auriparus flaviceps bill length vs height by which side of barrier
ggplot(good_data, aes(x=TARSUS.LENGTH, color=SPP))+geom_histogram()
ggplot(good_data, aes(x=SPP, y=BILL.LENGTH, color=SPP))+geom_boxplot()
ggplot(good_data, aes(x=SPP, y=TARSUS.LENGTH, color=SPP)) +geom_boxplot()       
ggplot(good_data, aes(x=WHICH.SIDE.OF.CFB, y=TARSUS.LENGTH, color=SPP)) +geom_boxplot()
fusca<-good_data[good_data$SPP=="FUSCA",]
ggplot(fusca, aes(x=WHICH.SIDE.OF.CFB, y=TARSUS.LENGTH, color=WHICH.SIDE.OF.CFB)) +geom_boxplot()

bellii<-good_data[good_data$SPP=="BELLII",]
ggplot(bellii, aes(x=WHICH.SIDE.OF.CFB, y=TARSUS.LENGTH, color=WHICH.SIDE.OF.CFB)) +geom_boxplot()


bilineata<-good_data[good_data$SPP=="BILINEATA",]
ggplot(bilineata, aes(x=WHICH.SIDE.OF.CFB, y=TARSUS.LENGTH, color=WHICH.SIDE.OF.CFB)) +geom_boxplot()