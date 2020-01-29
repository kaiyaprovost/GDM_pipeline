library(rworldmap)

points <- 
  read.csv("~/Documents/Dissertation/VertNet /allSpecies_pivot_11July2017.txt")


newmap = getMap(resolution="low")

png(filename="~/Documents/Dissertation/VertNet /test.png")
plot(newmap,xlim=c(-125,-80), ylim=c(10,45))
points(points$Long,points$lat,cex=0.5,col="red")
text(points$Long,points$lat,cex=0.5,labels=points$count,
     col="red")
dev.off()

spp = unique(points$species)
for(sp in spp){
  #print(sp)
  x = points[points$species==sp,]
  #write.csv(x,file=paste("~/Documents/Dissertation/VertNet /",sp,"_testTable.csv"))
  minlat = max(min(x$lat))
  minlong = max(min(x$Long),-125)
  maxlat = min(max(x$lat),45)
  maxlong = min(max(x$Long))
  ramp = colorRampPalette(c("blue","red"))
  cols = ramp(max(x$count))
  pdf(file=paste("~/Documents/Dissertation/VertNet /",sp,"_test.pdf"))
  #png(filename=paste("~/Documents/Dissertation/VertNet /",sp,"_test.png"),
  #    h=400,w=600)
  plot(newmap,xlim=c(minlong-1,maxlong+1),ylim=c(minlat-1,maxlat+1),
       main=sp)
  points(x$Long,x$lat,col="black",cex=1,pch=22,
         bg=cols[x$count])
  text(x$Long,x$lat,col="yellow",cex=0.4,labels=x$count)
  dev.off()
}

