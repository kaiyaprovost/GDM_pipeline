library(parallel)

## needs to be long running process
do = seq(1, 10000000)
p = proc.time()
l_works = sapply(do, sqrt) ## applies sqrt to every element in do
x = proc.time() - p
x

## now try do to in para
## set up number nodes
nclus = 4
## make a cluster with that many cores
cl = makeCluster(nclus, type ='FORK');
#"Forking" is the most basic type of paralellization using the system "FORK" command
## just makes 4 new R sessions dependent on this one, locking up this one
p = proc.time()
## split the cluster up by the parallel stuff
## one long vector of numbers and puts in equal parts depending on cluster
splits = clusterSplit(cl, do)
## do the parallel s apply on the cluster
## if you want to do single values, sapply. multiple values, lapply
p_works2 = parSapply(cl, splits, sqrt)
y = proc.time() - p
y
## and turn off the cluster if you're done
## always close when done -- safer not to reuse cluster
stopCluster(cl)


## remember though that you can easily make stuff worse if you do it wrong
## there is overhead to start the daughter processes and get them back
## for instance this one is not split
nclus = 4
cl = makeCluster(nclus, type ='FORK'); #"Forking" is the most basic type of paralellization using the system "FORK" command
p = proc.time()
p_works2 = parSapply(cl, do, sqrt)
w = proc.time() - p; #No faster than non-parallel
w
stopCluster(cl)

## try SOCK cluster -- more efficient
nclus = 4
cl = makeCluster(nclus, type ='SOCK'); 
p = proc.time()
splits = clusterSplit(cl, do)
p_works2 = parSapply(cl, splits, sqrt)
z = proc.time() - p; #About 25% quicker?
z
stopCluster(cl)

## also know that parSupply returns a list of lists, one for each node
length(p_works2)
length(unlist(p_works2))

## you also instead of splitting long job into parts can do other things

## now parallelize getting data from gbif, this is a bad function for this

#we need a simplified function that accepts a single argument to use the parSapply method
wrap_max = function(x){
  require(spocc)
  require(raster)
  occ=data.frame(occ2df(occ(x, limit=3000, from='gbif')));
  max = dismo::maxent(Env, occ[,2:3])
  return(max)
}

library(raster)

## /Users/kprovost/Documents/Classes/Spatial Bioinformatics/spatial_bioinformatics-master/ENM/wc2-5/
Env = raster::stack(list.files(
  path='/Users/kprovost/Documents/Classes/Spatial Bioinformatics/spatial_bioinformatics-master/ENM/wc2-5/',
  pattern="\\.bil$",
  full.names=T))
ext = raster::extent(c(-100, -40, 25, 45))
Env = raster::crop(Env, ext)

splist = c("Quercus virginiana", "Quercus geminata")

p=proc.time()
max_lin = lapply(splist, wrap_max)
## for two species, gets data, builds models
q = proc.time() - p;
q
max_lin

## now to compare do this with parallel
library(parallel)
nclus = 2
cl = makeCluster(nclus, type ='SOCK'); 
clusterExport(cl, list("Env")) ## export the object Env
## our function needs to be able to see it
p = proc.time()
max_par = parLapply(cl, splist, wrap_max); 
#We want parLapply this time because the data frames don't simplify well
r = proc.time() - p; 
r
stopCluster(cl)

print(max_par)
pr1 = predict(Env,max_par[[1]])
pr2 = predict(Env,max_par[[2]])
plot(pr1)
plot(pr2)

## could adjust this function to do instead of species, subsets of points
## and then do species models 