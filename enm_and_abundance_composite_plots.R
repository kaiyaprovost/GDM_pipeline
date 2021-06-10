wc_folder = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/worldclim_thresh/"
md_folder = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/mid_thresh/"
lg_folder = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/lgm_thresh/"

wc_files = list.files(path=wc_folder,pattern="asc$",full.names = T)
md_files = list.files(path=md_folder,pattern="asc$",full.names = T)
lg_files = list.files(path=lg_folder,pattern="asc$",full.names = T)


library(raster)
library(RColorBrewer)
yellow=colorRampPalette(brewer.pal(9,"YlOrRd")[3:9])
red=colorRampPalette(brewer.pal(12,"RdYlBu"))
col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#4393C3", "#2166AC", "#053061"))
colors=col2(11)

wc_stack = stack(wc_files)
wc_mean = calc(wc_stack,fun=mean,na.rm=T)
wc_sum = calc(wc_stack,fun=sum,na.rm=T)
wc_sum[is.na(wc_stack[[1]])] = NA

md_stack = stack(md_files)
#md_mean = calc(md_stack,fun=mean,na.rm=T)
md_sum = calc(md_stack,fun=sum,na.rm=T)
md_sum[is.na(md_stack[[1]])] = NA

lg_stack = stack(lg_files)
#lg_mean = calc(lg_stack,fun=mean,na.rm=T)
lg_sum = calc(lg_stack,fun=sum,na.rm=T)
lg_sum[is.na(lg_stack[[1]])] = NA

big_stack = stack(c(wc_files,lg_files,md_files))
#lg_mean = calc(lg_stack,fun=mean,na.rm=T)
big_sum = calc(big_stack,fun=sum,na.rm=T)
big_sum[is.na(big_stack[[1]])] = NA

pdf("stacks_of_stability_all_species.pdf",height=10,width=10)
par(mfrow=c(2,2))
plot(wc_sum,col=(yellow(11)),main="present")
plot(md_sum,col=(yellow(11)),main="mid")
plot(lg_sum,col=(yellow(11)),main="lgm")
plot(big_sum/3,col=(yellow(11)),main="average")
dev.off()

png("big_sum_stability_all_species.png")
#pdf("big_sum_stability_all_species.pdf")
plot(big_sum/3,col=(red(11)),main="Sum of Stability Across Species")
dev.off()



abun_files = list.files(path="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/BBS_data/abundance_asciis/clipped/",
                              pattern="Idw_abundance_clipped.asc$",full.names = T)

wcnt_folder = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER3_TRAITS/usa_nothresh/"
wcnt_files = list.files(path=wcnt_folder,pattern="asc$",full.names = T)

## need wc files that are not thresh

names(abun_files)=c("bel","bil","bru","cri","cur","fla","fus","mel","nit","sin")
names(wcnt_files) = c("bil","fla","bru","sin","fus","nit","mel","cri","cur","bel")

wcnt_files = wcnt_files[c(10,1,3,8,9,2,5,7,6,4)]

final_outputs = list()
logmods = list()

for(i in 1:10) {
  abd = raster(abun_files[i])
  enm = raster(wcnt_files[i])
  
  print(abd)
  print(enm)
  
  ext_a = extent(abd)
  ext_e = extent(enm)
  
  abd_c = crop(abd,ext_e)
  
  enm_c = aggregate(enm,20)
  #abd_cr = disaggregate(abd_c,20)
  
  abd_val = values(abd_c)
  enm_val = values(enm_c)
  
  #plot(enm)
  
  val = data.frame(abd_val,enm_val)
  val = unique(val)
  val = val[complete.cases(val),]
  names(val) = c("abd","enm")
  
  #mod = lm(val$enm~sqrt(sqrt(val$abd)))
  mod=lm(enm~sqrt(sqrt(abd)),data=val)
  
  #newx = (val$abd)
  newx=c(seq(min(val$abd,na.rm=T),max(val$abd,na.rm=T),length.out=500),
         seq(max(val$abd,na.rm=T),100,length.out=100))
  
  newy = predict(mod,data.frame(abd=newx))
  newval=cbind(newx,newy)
  final_outputs[[i]] = newval
  names(final_outputs)[[i]] = names(abun_files)[i]
  
  #png(paste(names(abun_files)[i],"_log_abundance.png"))
  logabd=log10(val$abd)
  
  mod2 = lm(val$enm~logabd)
  
  plot(log10(val[,1]),val[,2],
       col=rgb(0,0,0,0.1),
       main=names(abun_files)[i],xlab="Log Abundance",
       ylab="Habitat Suitability")
  abline(mod2,col="red")
  
  logmods[[i]] = mod2
  names(logmods)[[i]] = names(abun_files)[i]
  
  print(names(abun_files)[i])
  print(summary(mod2))
  
  #dev.off()
  
}

png("all_abundance_plots_curves_legend.png")
colors=c("black","red","orange","goldenrod","green",
         "blue","cyan","magenta","brown","grey")
plot(0,type="n",ylim=c(0,1),xlim=c(0,100),ylab="Habitat Suitability",
     xlab="Abundance")
for (j in 1:10) {
  
  ## first 500 are empirical, last 100 are projected
  
  this=final_outputs[[j]]
  this=this[order(this[,1]),]
  points((this[1:500,1]),(this[1:500,2]),col=adjustcolor(colors[j],alpha.f = 1),pch=16,cex=1)
  points((this[500:600,1]),(this[500:600,2]),col=adjustcolor(colors[j],alpha.f = 1),pch=3,cex=0.1)

}
legend("bottomright",col=colors,lty=1,lwd=3,pch=16,cex=1,legend=names(abun_files))

dev.off()


colors=c("black","black","black",
         "red","red","red",
         "blue","blue","blue",
         "goldenrod")

plot(0,type="n",ylim=c(0.4,0.9),xlim=c(-1,2),xaxt="n",
     ylab=("Habitat Suitability"),
     xlab=("Log Abundance"))
axis(1,labels=c(0.1,0.3,1,3,10,33,100,300),at=log10(c(0.1,0.3,1,3,10,33,100,300)),las=1)

for (j in 1:10) {
  #this=final_outputs[[j]]
  #this=this[order(this[,1]),]
  #lines(log10(this[,1]),(this[,2]),col=j,pch=j,lwd=1+j%%2)
  abline(logmods[[j]],col=colors[j],lwd=1+j%%2,lty=1+(j+2)%%3)
  legend("topleft",legend=names(logmods),
         col=colors,lty=1+(c(1:10)+2)%%3,ncol = 2,lwd=1+c(1:10)%%2,
         bty="n")
}
