library(rgdal)
library(rgeos)

path = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/OUTPUTS/2.processed/data/"
setwd(path)

listf = list.files(pattern="_thin1.csv")

for (file in listf) {

  print(file)
#file = "Toxostoma_curvirostre_ALL_thin1_USA.csv"
shapefile = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/ECOLOGY/enm_layers/mexstates/mexstates.shp"

tester=readOGR(shapefile)
points=read.table(file,row.names = NULL,header=T,sep=",")
rawpt = points
coordinates(points) <- ~longitude+latitude
proj4string(points) <- CRS(proj4string(tester))
ovdf = over(points, tester)
notover = which((is.na(ovdf[,1])))

notoverdf = rawpt[notover,]

notoverdf = notoverdf[notoverdf$latitude>=25.83,]

write.table(notoverdf,paste(file,"_USA.csv",sep=""),row.names = F,sep = ",")

}

