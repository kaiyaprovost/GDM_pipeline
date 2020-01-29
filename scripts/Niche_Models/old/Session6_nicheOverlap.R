#install.packages(c("ENMeval","phyloclim"))
require(dismo)
require(ENMeval)
require(phyloclim)
require(sp)
require(rgdal)
require(rgeos)
require(devtools)
#install_github("danlwarren/ENMTools")
require(ENMTools)

#################################################################################################################
#### Niche comparisons using Schoener's D. 
# Load two rasters to compare
# Model maxent rasters
path = "/Users/kprovost/Documents/Classes/GIS/Session6/Lemurdata/"
r1<-raster(paste(path,'E_m_flavifrons.asc',sep="")) # can be any raster type (i.e. .asci, .tif, .bil)
r2<-raster(paste(path,'E_m_macaco.asc',sep=""))

# Calculate niche overlap between the rasters. Change stat to 'I' for Moran's I test, and 
# to 'D' for Schoener's D. Remember to use ?nicheOverlap to see the arguments for this function. 
Overlap <- nicheOverlap(r1, r2, stat='I', mask=TRUE, checkNegatives=TRUE)
# This function is used by ENMeval to calculate the niche overlap between candidate models during the
# model tuning process. Set overlap to TRUE to get pairwise Schoener's D values for all candidate models. 

#################################################################################################################
#### Maxent model tuning using ENMeval
#### You don't have to do this if you have already tuned them
# Load environmental rasters as stack. Change pattern to .asc for ascii.
Env<-stack(list.files(path = "Pathway/to/environmental/data", pattern = '\\.tif$', full.names = T))
# Load locality information with 3 columns of "Species", "Longitude", "Latitude"
locs<-read.csv('Pathway/to/locality/data.csv')

# Here, method refers to data partitioning method. Categoricals refers to the names of any environmental data
# that are categorical. fc refers to feature classes to use. bg.coords refers to the background points used every time. You
# can set this as an object (i.e. a .csv file). RMvalues refers to the regularization multipliers to use. 
# Remember to check out ?ENMevaluate
res <- ENMevaluate(occ = locs, env = Env, method='block', categoricals=NULL, 
                   fc = c("L", 'LQ', "H", "LQH"), 
                   bg.coords=NULL, RMvalues=seq(1, 5, 0.5), overlap = T)

###################################################################################################################
#### Identity Test (same as Warren 2008, ENMtools)
# Load environmental variables
env<-stack(list.files(path="/Users/kprovost/Documents/Classes/GIS/Session6/Lemurdata/", pattern = '\\.asc$', full.names=T))
# Load occurrence records for both species
setwd("/Users/kprovost/Documents/Classes/GIS/Session6/Lemurdata/") # These should be csv files of records where columns are: "Species, X, Y".
flav <- read.csv("e_m_flav.csv")
maca <- read.csv("e_m_maca.csv")
# Change the species columns to just the species' names
flav[1] <- as.factor('flavifrons')
maca[1] <- as.factor('macaco')
# row bind them so all occurrences are in 3 rows of Species, X, Y
sites<-rbind(maca, flav)
species <- c('flavifrons','macaco')
# Change the column names of sites
colnames(sites)<-c("species","longitude","latitude")
samples <- sites[grep(paste(species, collapse = "|"), sites$species), ]
# Tell R where maxent is (the copy that in with dismo).
maxent.exe <- paste(system.file(package="dismo"),"/java/maxent.jar", sep = "")
### ?niche.equivalency.test
setwd('/Users/kprovost/Documents/Classes/GIS/Session6/')
nicheEquivalency<-niche.equivalency.test(p = samples, env = env, app=maxent.exe, dir = 'NicheEquivalence')
# Note that you can also perform the background test using bg.similarity.test. For more see ?bg.similarity.test
## This will take like 20 minutes, a bg test will take like 40
print(nicheEquivalency)

setwd("/Users/kprovost/Documents/Classes/GIS/Session6/")
png("nicheEquivalency.png",width=600,height=480,units="px")
plot(nicheEquivalency)
dev.off()

#################################################################################################################
#### Get background data from minimum convex polygon of occurrence data.
# Load occurrence data
setwd("/Users/kprovost/Documents/Classes/GIS/Session6/Lemurdata/")
flav <- read.csv("e_m_flav.csv")
maca <- read.csv("e_m_maca.csv")
# Load environmental data
env <- stack(list.files(path="/Users/kprovost/Documents/Classes/GIS/Session6/Lemurdata/", pattern = '\\.asc$', full.names=T))
# Convert a raster to a polygon shape
envPoly <- rasterToPolygons(env[[1]], fun=NULL, na.rm=T)
# To create a minimum convex polygon, this function need to be sourced. Highlight these lines and run them.
mcp <- function (xy) {
  xy <- as.data.frame(coordinates(xy))
  coords.t <- chull(xy[, 1], xy[, 2])
  xy.bord <- xy[coords.t, ]
  xy.bord <- rbind(xy.bord[nrow(xy.bord), ], xy.bord)
  return(SpatialPolygons(list(Polygons(list(Polygon(as.matrix(xy.bord))), 1))))
}
# Now run this function on the occurrence records
MCP.locs.flav <- mcp(flav[,2:3])
MCP.locs.maca <- mcp(maca[,2:3])
# Clip the polygon of the entire study region by the MCP polygon
background.area.flav <- gIntersection(envPoly, MCP.locs.flav)
background.area.maca <- gIntersection(envPoly, MCP.locs.maca)
# Convert this back to a raster with the same attributes as the environmental data
MCP.raster.flav <- rasterize(background.area.flav, env[[1]])
MCP.raster.maca <- rasterize(background.area.maca, env[[1]])
# Generate background points from this MCP raster. remember you want 100 times more points than occurrence localities.
bg.points.flav <- randomPoints(mask = MCP.raster.flav, n = (100*nrow(flav)))
bg.points.maca <- randomPoints(mask = MCP.raster.maca, n = (100*nrow(maca)))

###################################################################################################################
#### Background test  Remember to use: ?bg.similarity.test
setwd('pathway/to/occurrence/records')
flav <- read.csv("e_m_flav.csv")
maca <- read.csv("e_m_maca.csv")
# Change the species columns to just the species' names
flav[1] <- as.factor('flavifrons') # Change this to your species 1's name
maca[1] <- as.factor('macaco') # Change this to your species 2's name
# row bind them so all occurrences are in 3 rows of Species, X, Y
sites<-rbind(maca, flav)
species <- c('flavifrons','macaco')
# Change the column names of sites
colnames(sites)<-c("species","longitude","latitude")
samples <- sites[grep(paste(species, collapse = "|"), sites$species), ]
# Load environmental data
env <- stack(list.files(path='Pathway/to/session/6/data/Session6_data/Lemur_layers', pattern = '\\.asc$', full.names=T))
env.flav<-mask(env, MCP.raster.flav)
env.maca<-mask(env, MCP.locs.maca)
# Tell R where maxent is (the copy that in with dismo).
maxent.exe <- paste(system.file(package="dismo"),"/java/maxent.jar", sep = "")
# set the working directory
setwd('Pathway/to/save/test/Session6_data')
# Perform background similarity test. Remember to use ?bg.similarity.test
bg.test <- bg.similarity.test(p = samples, env = env, app = maxent.exe, dir = 'background', n=2)
