
pca     = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/lostruct_2021_done/not_for_gdms/Melozone.fusca.fixed.sorted.nospace.sorted.pca.csv"
coords  = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/lostruct_2021_done/not_for_gdms/Melozone.fusca.fixed.sorted.nospace.sorted.coords.csv"
regions = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/lostruct_2021_done/not_for_gdms/Melozone.fusca.fixed.sorted.nospace.sorted.regions.csv"

pca_df = read.table(pca,header=T,sep=",")
coords_df = read.table(coords,header=T,sep=",")
regions_df = read.table(regions,header=T,sep=",")

head(pca_df)
head(coords_df)
head(regions_df)

pca_regions = cbind(regions_df,pca_df)
pca_regions_coords = cbind(coords_df,pca_regions[complete.cases(pca_regions),])

black = pca_regions_coords[pca_regions_coords$ccols=="black",]
green = pca_regions_coords[pca_regions_coords$ccols=="#1B9E77",]
purple = pca_regions_coords[pca_regions_coords$ccols=="#7570B3",]
orange = pca_regions_coords[pca_regions_coords$ccols=="#D95F02",]

summary(black)
png("test.png")
plot(pca_regions_coords$MDS1,pca_regions_coords$MDS2,col=pca_regions_coords$ccols)
dev.off()

fus_deserts = c("green","green","blue","blue","cyan","green","cyan","green","cyan","green","cyan","green","cyan","cyan","cyan","cyan","green","cyan","cyan","cyan","cyan","blue","blue","green")

par(mfrow=c(2,3))
## max min mean? 

## sort and then plot? 

black = black[order(black$MDS1, black$MDS2), ]
plot(as.numeric(black[1,15:(14+ncol(black[,-c(1:14)])/2)]),
     as.numeric(black[1,(15+ncol(black[,-c(1:14)])/2):ncol(black)]),pch=letters,col=fus_deserts)
plot(as.numeric(black[floor(nrow(black)/2),15:(14+ncol(black[,-c(1:14)])/2)]),
     as.numeric(black[floor(nrow(black)/2),(15+ncol(black[,-c(1:14)])/2):ncol(black)]),pch=letters,col=fus_deserts)
plot(as.numeric(black[nrow(black),15:(14+ncol(black[,-c(1:14)])/2)]),
     as.numeric(black[nrow(black),(15+ncol(black[,-c(1:14)])/2):ncol(black)]),pch=letters,col=fus_deserts)
black = black[order(black$MDS2, black$MDS1), ]
plot(as.numeric(black[1,15:(14+ncol(black[,-c(1:14)])/2)]),
     as.numeric(black[1,(15+ncol(black[,-c(1:14)])/2):ncol(black)]),pch=letters,col=fus_deserts)
plot(as.numeric(black[floor(nrow(black)/2),15:(14+ncol(black[,-c(1:14)])/2)]),
     as.numeric(black[floor(nrow(black)/2),(15+ncol(black[,-c(1:14)])/2):ncol(black)]),pch=letters,col=fus_deserts)
plot(as.numeric(black[nrow(black),15:(14+ncol(black[,-c(1:14)])/2)]),
     as.numeric(black[nrow(black),(15+ncol(black[,-c(1:14)])/2):ncol(black)]),pch=letters,col=fus_deserts)

green = green[order(green$MDS1, green$MDS2), ]
plot(as.numeric(green[1,15:(14+ncol(green[,-c(1:14)])/2)]),
     as.numeric(green[1,(15+ncol(green[,-c(1:14)])/2):ncol(green)]),pch=letters,col=fus_deserts)
plot(as.numeric(green[floor(nrow(green)/2),15:(14+ncol(green[,-c(1:14)])/2)]),
     as.numeric(green[floor(nrow(green)/2),(15+ncol(green[,-c(1:14)])/2):ncol(green)]),pch=letters,col=fus_deserts)
plot(as.numeric(green[nrow(green),15:(14+ncol(green[,-c(1:14)])/2)]),
     as.numeric(green[nrow(green),(15+ncol(green[,-c(1:14)])/2):ncol(green)]),pch=letters,col=fus_deserts)
green = green[order(green$MDS2, green$MDS1), ]
plot(as.numeric(green[1,15:(14+ncol(green[,-c(1:14)])/2)]),
     as.numeric(green[1,(15+ncol(green[,-c(1:14)])/2):ncol(green)]),pch=letters,col=fus_deserts)
plot(as.numeric(green[floor(nrow(green)/2),15:(14+ncol(green[,-c(1:14)])/2)]),
     as.numeric(green[floor(nrow(green)/2),(15+ncol(green[,-c(1:14)])/2):ncol(green)]),pch=letters,col=fus_deserts)
plot(as.numeric(green[nrow(green),15:(14+ncol(green[,-c(1:14)])/2)]),
     as.numeric(green[nrow(green),(15+ncol(green[,-c(1:14)])/2):ncol(green)]),pch=letters,col=fus_deserts)

purple = purple[order(purple$MDS1, purple$MDS2), ]
plot(as.numeric(purple[1,15:(14+ncol(purple[,-c(1:14)])/2)]),
     as.numeric(purple[1,(15+ncol(purple[,-c(1:14)])/2):ncol(purple)]),pch=letters,col=fus_deserts)
plot(as.numeric(purple[floor(nrow(purple)/2),15:(14+ncol(purple[,-c(1:14)])/2)]),
     as.numeric(purple[floor(nrow(purple)/2),(15+ncol(purple[,-c(1:14)])/2):ncol(purple)]),pch=letters,col=fus_deserts)
plot(as.numeric(purple[nrow(purple),15:(14+ncol(purple[,-c(1:14)])/2)]),
     as.numeric(purple[nrow(purple),(15+ncol(purple[,-c(1:14)])/2):ncol(purple)]),pch=letters,col=fus_deserts)
purple = purple[order(purple$MDS2, purple$MDS1), ]
plot(as.numeric(purple[1,15:(14+ncol(purple[,-c(1:14)])/2)]),
     as.numeric(purple[1,(15+ncol(purple[,-c(1:14)])/2):ncol(purple)]),pch=letters,col=fus_deserts)
plot(as.numeric(purple[floor(nrow(purple)/2),15:(14+ncol(purple[,-c(1:14)])/2)]),
     as.numeric(purple[floor(nrow(purple)/2),(15+ncol(purple[,-c(1:14)])/2):ncol(purple)]),pch=letters,col=fus_deserts)
plot(as.numeric(purple[nrow(purple),15:(14+ncol(purple[,-c(1:14)])/2)]),
     as.numeric(purple[nrow(purple),(15+ncol(purple[,-c(1:14)])/2):ncol(purple)]),pch=letters,col=fus_deserts)

orange = orange[order(orange$MDS1, orange$MDS2), ]
plot(as.numeric(orange[1,15:(14+ncol(orange[,-c(1:14)])/2)]),
     as.numeric(orange[1,(15+ncol(orange[,-c(1:14)])/2):ncol(orange)]),pch=letters,col=fus_deserts)
plot(as.numeric(orange[floor(nrow(orange)/2),15:(14+ncol(orange[,-c(1:14)])/2)]),
     as.numeric(orange[floor(nrow(orange)/2),(15+ncol(orange[,-c(1:14)])/2):ncol(orange)]),pch=letters,col=fus_deserts)
plot(as.numeric(orange[nrow(orange),15:(14+ncol(orange[,-c(1:14)])/2)]),
     as.numeric(orange[nrow(orange),(15+ncol(orange[,-c(1:14)])/2):ncol(orange)]),pch=letters,col=fus_deserts)
orange = orange[order(orange$MDS2, orange$MDS1), ]
plot(as.numeric(orange[1,15:(14+ncol(orange[,-c(1:14)])/2)]),
     as.numeric(orange[1,(15+ncol(orange[,-c(1:14)])/2):ncol(orange)]),pch=letters,col=fus_deserts)
plot(as.numeric(orange[floor(nrow(orange)/2),15:(14+ncol(orange[,-c(1:14)])/2)]),
     as.numeric(orange[floor(nrow(orange)/2),(15+ncol(orange[,-c(1:14)])/2):ncol(orange)]),pch=letters,col=fus_deserts)
plot(as.numeric(orange[nrow(orange),15:(14+ncol(orange[,-c(1:14)])/2)]),
     as.numeric(orange[nrow(orange),(15+ncol(orange[,-c(1:14)])/2):ncol(orange)]),pch=letters,col=fus_deserts)

