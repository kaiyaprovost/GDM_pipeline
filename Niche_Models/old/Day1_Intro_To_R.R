library(data.table)
library(tidyverse)
library(mapdata)
library(maps)
library(ggmap)
library(magrittr)
library(devtools)
library(ggplot2)
library(ggmap)
devtools::install_github("dkahle/ggmap")
devtools::install_github("hadley/ggplot2")

##fread is faster read, gives all the ram up front, not chunk by chunk
ebd = fread("~/Documents/Classes/Spatial Bioinformatics/ebd_trim.txt")
#ebd = read.delim("~/Documents/Classes/Spatial Bioinformatics/ebd_trim.txt",
#                 sep="\t")

## note
## readcsv2 is faster but doesn't put in dataframe

## toy matrix to play with
y = cbind(seq(1:5), 
          seq(1:5),
          seq(1:5),
          seq(1:5),
          seq(1:5))
y = as.data.frame(y)

li1 = apply(y,1,log) ##takes the log of items by row
li2 = apply(y,2,log) ## takes the log of items by col -- transposes? thinks transposes because of the log function?

li1 = apply(y,1,sum) ##takes the sum of items by row
li2 = apply(y,2,sum) ## takes the sum of items by col

li3 = lapply(y[,1],log) ##returns list, takes the log of items in col 1
li4 = sapply(y[,1],log) ##returns vector, takes the log of items in col 1

## replicate sapply 10 times
rep = replicate(10,log(y[,1]),simplify="array")

## time to learn ggplot 
virus_data = read.delim("/Users/kprovost/Documents/Classes/Spatial Bioinformatics/spatial_bioinformatics-master/data/R_workshop_data.txt")
head(virus_data,n=10)
virus_data$frequency=virus_data$count/virus_data$coverage
## OR
virus_data$frequency=with(virus_data,count/coverage)

## make scatterplot, basic
## tell it data, tell it where to plot points
ggplot(data=virus_data) +
  geom_point(mapping=aes(x=position,y=frequency))
## add labels
+ xlab("Genome Position") + ylab("SNP Frequency")

## to make them transparent need to add alpha to aes
ggplot(data=virus_data) +
  geom_point(mapping=aes(x=position,y=frequency),alpha=0.5) + 
  xlab("Genome Position") + ylab("SNP Frequency")

## to make them sort by different things? make multiple panels
## this is called faceting -- adding filter to it 
## in this case by host type? 
## first want only to show things if they are hosts, not vectors
ggplot(data = filter(virus_data, host_type == 'host')) +
  ## then map them scatterplot as usual, with labels
  geom_point(mapping = aes(x = position, y = frequency)) +
  xlab("Genome Position") + ylab("SNP Frequency") +
  ## then tell it to sort by the host type, in this case to different panels
  facet_wrap(~ host)


## now make a map!
## first thing is you want to get a lat-long of a point
loc = as.data.frame(cbind(-73.180088,44.532371))
loc2 = cbind(loc,c(-108.4382689, 31.7000670))
colnames(loc) = c("lon","lat")
## now you can plot it on a google maps image
bkmap = get_map(location=loc,maptype="satellite",source="google",zoom=25) #14
## now plot the map with a point where the location is 
ggmap(bkmap) + 
  ## add the point
  geom_point(data=loc,color="red",size=4)

## now do it with terrain 
bkmap3 = get_map(location=loc,maptype="terrain",source="google",zoom=17) ## 12 is zoomed out more than 14
ggmap(bkmap3) + 
  ## add the point
  geom_point(data=loc,color="red",size=4)

bkbox = make_bbox(lon=loc$lon,lat=loc$lat,f=0.1)


