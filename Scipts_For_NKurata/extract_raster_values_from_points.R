## load packages
library(dismo)
library(rgdal)
library(raster)

## read points
points=read.csv("~/cell_locations_raster_enms.csv")
points=round(points,5)
points=unique(points)

## make coordinates file from lat and long
coor = coordinates(points[c("LONG","LAT")])

## read raster want to extract from
ras = raster("~/amphispiza_maxent.asc")

## extract values from the raster via the lat and long
## Note: extract requires LONG LAT, not LAT LONG
extracted = data.frame(coor,extract(ras,coor,small=T)); names(extracted) = c("LONG","LAT","EXTRACTED")

extracted = extracted[order(extracted$LONG),]

write.table(extracted,"~/suitability.csv")
