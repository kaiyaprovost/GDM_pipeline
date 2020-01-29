#code for converting rater to ascii
library(raster)
f <- "path/to/downloaded/file.tif"
r <- raster(f)
ra <- aggregate(r, fact=2)  ## By default aggregates using mean, but see fun=
writeRaster(ra, "path/to/outfile.asc", format="ascii")

#my code
f <- "C:/users/xphil/Documents/RGGS GIS/NLCD2011/"
r <- raster(f)
ra <- aggregate(r, fact=2)
writeRaster(ra, "C:/users/xphil/Documents/RGGS GIS/Global_wc/testBio2R.asc", format="ascii")

#batch convert from Peter
#set wd first
bio <- stack(list.files(pattern = '.tif', full.names = T))
plot(bio[[1]])  #plot a graph to check out how it looks
e<- c(-90, -45, 25, 65) #specify extent of cropped region
bio3<- crop(bio2, e)  #crop to e (bio2 was Peter's first attemp
#at cropping, but it was still too big, so he cropped again and 
#saved to bio3)
plot(bio3[[]])   #plot that new crop
setwd('E:/name') #where to save output file
for(i in 1:nlayers(bio3)) #for each interation of the 19 bioclim layers
{writeRaster(bio3[[i]], filename = names(bio3[[i]]), format='ascii')
}