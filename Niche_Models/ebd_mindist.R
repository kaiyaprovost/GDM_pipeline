library(raster)
library(devtools)
install_github("rsh249/rasterExtras")
library(rasterExtras)
library(data.table)
nclus = 24;

## this script will measure the density of points and 
## output an image of that density 


#path = "/Users/kprovost/Documents/Classes/Finished Classes/Spatial Bioinformatics/ebd_trim.csv"
path = "/Users/kprovost/Documents/Dissertation/CHAPTER3_TRAITS/ECOLOGY/OUTPUTS/2.processed/data/Amphispiza_bilineata_thin1.csv"

ebd = data.table::fread(path, sep = ',')

Env = raster::stack(
  list.files(
    path = '/Users/kprovost/Documents/Classes/Spatial Bioinformatics/spatial_bioinformatics-master/ENM/wc2-5/',
    pattern = "\\.bil$",
    full.names = T
  )
)
ext = raster::extent(c(-125,-60, 10, 50)) ## make sure this will play nice with your points
Env = raster::crop(Env, ext)
bg = Env[[1]] ## just for plotting
clim = Env

#ebd = fread('ebd_trim2.csv', sep = '\t')
#clim = stack("../climgrids/bio.gri")
#clim = clim[[1]];
#ext = extent(c(-125, -100, 25, 45));
#clim = crop(clim,ext);

## before was 4:3 and 27:26
stackit = rasterExtras::dist2point(clim, ebd[,2:3], parallel=TRUE, nclus=nclus, maxram = 256)

merged = stackit;

save.image('ebd.mindist.RData')

plot(merged, col = c('black', viridis::magma(9999)));
usa1 = getData("GADM", cou='US', level=1)
plot(usa1, add=TRUE, border='white')
dev.off()



png('ebird_dist.png', height = 8, width = 8, res = 600, units = 'in')
plot(merged, col = c('black', viridis::magma(9999)));
usa1 = getData("GADM", cou='US', level=1)
plot(usa1, add=TRUE, border='white')
dev.off()


merged.scaled = 1-(merged/maxValue(merged)); ##Scaled 0 = high distance, 1 = low distance... Proportional to probability?
png('ebird_distscaled.png', height = 8, width = 8, res = 600, units = 'in')
plot(merged.scaled, col = c('black', viridis::magma(9999)));
usa1 = getData("GADM", cou='US', level=1)
plot(usa1, add=TRUE, border='white')
dev.off()

save.image('ebd.mindist.RData');

q('no');
