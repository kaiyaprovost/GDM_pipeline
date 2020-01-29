## landsim
library(landsim)
library(raster)
library(plotly)

## make a population
pop <- make_population(
  habitat = system.file("extdata/test_raster.gri",package="landsim"), ## what does the habitat look like?
  inaccessible.value = "NA", ## what are inaccessible values
  uninhabitable.value = 0.0, ## what are uninhabitable values
  genotypes = c( "aa", "aA", "AA" ) ## what are the possible genotypes 
)

## make a migration object
migr <- migration(
  kern = "gaussian", ## distance and relative weights
  sigma = 100, ## scaling for kern
  radius = 1000, ## maximum distance to migrate, here meters
  normalize = 1.0 ## migrants per unit of input 
)

## set up migration matrix with your population
migr <- setup_migration( migr, population=pop )

## create a "vital rate"
## these can be functions, numbers, etc
## "This class allows assigning and accessing 
## things directly from the environment of the 
## function, without having them clutter up the 
## global environment."

f <- vital( function (N,...) { N+x }, x=3 )
f(10)
f$x
f$x <- 12
f(10)

germ.vital <- vital(
  function (N, ...) {
    ## in this case is a function that calculates the probability of germination
    out <- r0 / ( 1 + migrate(rowSums(N),competition)/carrying.capacity )
    return( cbind( aa=out, aA=s*out, AA=s^2*out ) )
  },
  r0 = 0.01,  # one in ten seeds will germinate at low densities
  s = 1.5,    # multiplicative selective benefit of the A allele
  carrying.capacity = 10,
  competition = migration(
    kern="gaussian",
    sigma=100,
    radius=300,
    normalize=1
  )
)

## set up the migration for germ.vital
germ.vital <- setup_vital(germ.vital,pop)

## set up demography
demog <- demography(
  prob.seed = 0.2, ## probability of plant seeding, per ind, per year
  fecundity = 100, ## mean number of seeds produced per individual per year
  prob.germination = germ.vital, ## vital function for calculating the probability of germination
  prob.survival = 0.9,  ## probability that living individual survives 
  seed.migration = migr, ## migration of the seeds 
  pollen.migration = migration( ## migration of the pollen 
    kern="gaussian",
    sigma=300,
    radius=900,
    normalize=10
  ),
  genotypes = c("aa", "aA", "AA")
)

demog <- setup_demography(demog,pop)









## Now lets do density dependence
## make a population
pop <- make_population( 
  habitat = random_habitat(),
  inaccessible.value = NA,
  uninhabitable.value = NA,
  genotypes = c("aa","aA","AA"),
  N = 0 ## values to initialize the matrix of genotype counts with, (habitable cells) x (genotypes) giving the number of each genotype in each habitable cell
)
## make a raster? based on genotype?
pop$N[,"aa"] <- rpois(nrow(pop$N),values(pop$habitat)[pop$accessible]/4)
pop$N[,"aA"] <- rpois(nrow(pop$N),values(pop$habitat)[pop$accessible]/2)
pop$N[,"AA"] <- rpois(nrow(pop$N),values(pop$habitat)[pop$accessible]/4)
plot(pop$habitat)

## generate the basic demography
basic.migr <- migration(
  kern = "gaussian",
  sigma = 300,
  radius = 1000,
  normalize = 1
)
basic.demog <- demography(
  prob.seed = 0.05,
  fecundity = 200,
  prob.germination = 0.4,
  prob.survival = 0.6,
  pollen.migration = basic.migr,
  seed.migration = basic.migr,
  genotypes = c("aa","aA","AA")
)

## apply the demography to the population
## then simulate
demog <- setup_demography( basic.demog, pop )
sim <- simulate_pop( pop, demog, times=seq(0,10,length.out=11),
                     summaries=list( totals=function(N){colSums(N)} )
)
matplot(sim$summaries[["totals"]],type='l',lty=1, log='y',ylab='number of individuals')
legend("bottomright",lty=1,col=1:3,legend=pop$genotypes)

## now invoke resource constraints on space
demog <- basic.demog
demog$prob.germination <- vital(
  function (N,...) {
    out <- r0 / ( 1 + migrate(competition,x=rowSums(N))/K )
    cbind( aa=out, aA=out, AA=out )
  },
  r0 = 0.4,
  K = values(pop$habitat)[pop$habitable]/5,
  competition = migration(
    kern="gaussian",
    sigma=200,
    radius=400,
    normalize=1
  )
)

## this invokes a carrying capacity
demog <- setup_demography( demog, pop )
sim <- simulate_pop( pop, demog, times=seq(0,100,length.out=101),
                     summaries=list( totals=function(N){colSums(N)} )
)
matplot(sim$summaries[["totals"]],type='l',lty=1, log='y', ylab="number of individuals")
legend("bottomright",lty=1,col=1:3,legend=pop$genotypes)

## and looking at it through time
plot(sim,pop,pause=FALSE)

## now try carrying capacity with death constraints
## death depends on local density
demog <- basic.demog
demog$prob.seed <- 0.01
demog$prob.survival <- vital(
  function (N,...) {
    out <- s0 / ( 1 + migrate(competition,x=rowSums(N))/K )
    cbind( aa=out, aA=out, AA=out )
  },
  s0 = 0.6,
  K = values(pop$habitat)[pop$habitable]/3,
  competition = migration(
    kern="gaussian",
    sigma=200,
    radius=400,
    normalize=1
  )
)

demog <- setup_demography( demog, pop )
sim <- simulate_pop( pop, demog, times=seq(0,100,length.out=101),
                     summaries=list( totals=function(N){colSums(N)} ) )
matplot(sim$summaries[["totals"]],type='l',lty=1, log='y')
legend("bottomright",lty=1,col=1:3,legend=pop$genotypes)
#plot(sim,pop,pause=FALSE)












## moving on to types of selection
## generate a habitat
pop <- make_population(
  habitat = random_habitat(),
  inaccessible.value = NA,
  uninhabitable.value = NA,
  genotypes = c("aa","aA","AA"),
  N = 0
)
pop$N[,"aa"] <- rpois(nrow(pop$N),values(pop$habitat)[pop$habitable])
pop$N[,"aA"] <- rpois(nrow(pop$N),values(pop$habitat)[pop$habitable]*100/sum(values(pop$habitat),na.rm=TRUE))
pop$N[,"AA"] <- 0
plot(pop$habitat)

## create basic demography
## with competition
basic.migr <- migration(
  kern = "gaussian",
  sigma = 300,
  radius = 1000,
  normalize = 1
)
basic.demog <- demography(
  prob.seed = 0.05,
  fecundity = 200,
  prob.germination = vital(
    function (N,...) {
      out <- r0 / ( 1 + migrate(competition,x=rowSums(N))/K )
      cbind( aa=out, aA=out, AA=out )
    },
    r0 = 0.4,
    K = values(pop$habitat)[pop$habitable]/5,
    competition = migration(
      kern="gaussian",
      sigma=200,
      radius=400,
      normalize=1
    )
  ),
  prob.survival = 0.6,
  pollen.migration = basic.migr,
  seed.migration = basic.migr,
  genotypes = c("aa","aA","AA")
)

## create directional selection based on genotype
## prob of germination depends on genotype
demog <- basic.demog
demog$prob.germination <- vital(
  function (N,...) {
    out <- r0 / ( 1 + migrate(competition,x=rowSums(N))/K )
    cbind( aa=s[1]*out, aA=s[2]*out, AA=s[3]*out )
  },
  s = c(aa=1,aA=1.4,AA=1.4^2),
  r0 = 0.4,
  K = values(pop$habitat)[pop$habitable]/5,
  competition = migration(
    kern="gaussian",
    sigma=200,
    radius=400,
    normalize=1
  )
)

## the AA allele has an increasing carrying capacity
demog <- setup_demography( demog, pop )
sim <- simulate_pop( pop, demog, times=seq(0,100,length.out=101),
                     summaries=list( totals=function(N){colSums(N)} ) )
matplot(sim$summaries[["totals"]],type='l',lty=1, log='y', ylab="number of individuals")
## this may generate a warning
legend("bottomright",lty=1,col=1:3,legend=pop$genotypes)
#plot(sim,pop,pause=FALSE)

## now create hybrid disadvantage
demog$prob.germination$s <- c( aa=1, aA=0.5, AA=1 )
pop$N[,"aa"] <- rpois(nrow(pop$N),values(pop$habitat)[pop$habitable]/3)
pop$N[,"aA"] <- rpois(nrow(pop$N),values(pop$habitat)[pop$habitable]/3)
pop$N[,"AA"] <- rpois(nrow(pop$N),values(pop$habitat)[pop$habitable]/3)
demog <- setup_demography( demog, pop )
sim <- simulate_pop( pop, demog, times=seq(0,100,length.out=101),
                     summaries=list( totals=function(N){colSums(N)} ) )
matplot(sim$summaries[["totals"]],type='l',lty=1, log='y', ylab="number of individuals")
legend("bottomright",lty=1,col=1:3,legend=pop$genotypes)
#plot(sim,pop,pause=FALSE)

## now local adaptation
## set up a gradient
## one homozygote good in red, one good in green, heterozygote okay in yellow 
grad <- pop$habitat
values(grad) <- colFromCell(grad,1:ncell(grad))/ncol(grad) - 0.5
plot(grad,main="gradient")

## create the demography
demog <- basic.demog
demog$prob.seed <- 0.01
demog$prob.survival <- vital(
  function (N,...) {
    out <- s0 / ( 1 + migrate(competition,x=rowSums(N))/K )
    cbind( aa=(1-s)*out, aA=out, AA=(1+s)*out )
  },
  s0 = 0.6,
  s = values(grad)[pop$habitable],
  K = values(pop$habitat)[pop$habitable]/3,
  competition = migration(
    kern="gaussian",
    sigma=200,
    radius=400,
    normalize=1
  )
)

demog <- setup_demography( demog, pop )
sim <- simulate_pop( pop, demog, times=seq(0,100,length.out=101),
                     summaries=list( totals=function(N){colSums(N)} ) )
matplot(sim$summaries[["totals"]],type='l',lty=1, log='y', ylab="number of individuals")
legend("bottomright",lty=1,col=1:3,legend=pop$genotypes)
#plot(sim,pop,pause=FALSE)

## soft selection
## before was changing the carrying capacity of genotypes based on factors
## now keep carrying capacity the same overall across genotype
demog <- basic.demog
demog$prob.germination <- vital(
  function (N,...) {
    P <- (N[,2]/2+N[,3])
    P[P>0] <- P[P>0]/(rowSums(N)[P>0])
    out <- r0 / ( 1 + migrate(competition,x=rowSums(N))/K )
    out <- cbind( aa=(1-P*s)*out, aA=out, AA=(1+(1-P)*s)*out )
    if (any(out<0)||any(out>1)||any(!is.finite(out))) { browser () }
    out
  },
  s = 0.4,
  r0 = 0.4,
  K = values(pop$habitat)[pop$habitable]/5,
  competition = migration(
    kern="gaussian",
    sigma=200,
    radius=400,
    normalize=1
  )
)

## again AA is best 
demog <- setup_demography( demog, pop )
sim <- simulate_pop( pop, demog, times=seq(0,100,length.out=101),
                     summaries=list( totals=function(N){colSums(N)} ) )
matplot(sim$summaries[["totals"]],type='l',lty=1, log='y', ylab="number of individuals")
legend("bottomright",lty=1,col=1:3,legend=pop$genotypes)
#plot(sim,pop,pause=FALSE)










## spatial sweeps
## generate a habitat
diam <- 1e4
habitat <- raster(xmn=-diam, xmx=diam, ymn=-diam, ymx=diam, 
                  resolution=100,
                  crs="+proj=utm +zone=11 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
values(habitat) <- pmin(20,(2+rcauchy(ncell(habitat))))
habitat <- 20*migrate_raster( habitat, kern="gaussian", sigma=300, radius=1500 )
values(habitat)[values(habitat)<0] <- NA
plot(habitat)

## add strong selection for A allele, additive
## Aa 10% more likely to germinate than aa, AA 20% more likely
demog <- demography(
  prob.seed = 0.2,
  fecundity = 100,
  prob.germination = germination_fun <- vital(
    function (N, carrying.capacity, ...) {
      out <- r0 / ( 1 + migrate(rowSums(N),competition)/carrying.capacity )
      return( cbind( aa=out, aA=s*out, AA=s^2*out ) )
    },
    r0 = 0.01,  # one in ten seeds will germinate at low densities
    s = 1.5,    # multiplicative selective benefit of the A allele
    competition = migration(
      kern="gaussian",
      sigma=100,
      radius=300,
      normalize=1
    )
  ),
  prob.survival = 0.9,
  pollen.migration = migration(
    kern = function (x) { exp(-sqrt(x)) },
    sigma = 30,
    radius = 500,
    normalize = NULL
  ),
  seed.migration = migration(
    kern = "gaussian",
    sigma = 50,
    radius = 400,
    normalize = 1
  ),
  genotypes = c("aa","aA","AA")
)

## start population with very few A alleles
total.habitat <- sum(values(habitat),na.rm=TRUE)
pop <- population( 
  habitat = habitat,
  genotypes = c("aa","aA","AA"),
  N = cbind( aa=rpois_raster(habitat,only.values=TRUE),
             aA=rpois_raster(habitat*40/total.habitat,only.values=TRUE),
             AA=0 )
)
demog <- setup_demography( demog, pop )
pl(pop)

## cannot find function pl
## think is the function in hte next section

## calculate the number of offspring per individual
## at low density?
base.r <- intrinsic_growth( pop, demog, 
                            carrying.capacity=values(pop$habitat)[pop$habitable] )
pl(pop$habitat, v=base.r-1)

## at higher density?
base.r <- intrinsic_growth( pop, demog, density=500,
                            carrying.capacity=values(pop$habitat)[pop$habitable] )
pl(pop$habitat, v=base.r-1)

## do not understand what this is doing
## run simulation, passing in carrying capacity
plot.times <- seq(0,700,length.out=71)
sim <- simulate_pop( pop, demog, times=plot.times, ## simulate_pop was originally simulate only
                 carrying.capacity=values(pop$habitat)[pop$habitable],
                 summaries=list( totals=function(N){colSums(N)} )
)

## plot this?
for (k in seq_along(plot.times)) {
  pl(pop$habitat,v=sim$N[,,k],main=c("",sprintf("t=%d",plot.times[k])),zlim=range(sim$N,finite=TRUE))
}

## more plotting
matplot( sim$summaries[["totals"]], type='l', xlab='time', ylab='numbers', lty=1 )
legend("topleft",lty=1,col=1:3,legend=colnames(sim$summaries[["totals"]]))








## simple selective sweep on a landscape
## generate the landscape
habitat <- raster(xmn=-1000, xmx=1000, ymn=-1000, ymx=1000, 
                  resolution=100,
                  crs="+proj=utm +zone=11 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
values(habitat) <- exp(2+rnorm(ncell(habitat)))
habitat <- migrate_raster( habitat, kern="gaussian", sigma=100, radius=1500 )
plot(habitat)

## use a plot function
## this is the pl function above, probably
# plotting function for typical layers below
pl <- function (x,v,nr=1,zlim, ...) {
  if (!missing(v)) {
    if (nlayers(x)<NCOL(v)) {
      x <- do.call(stack,list(x)[rep(1,NCOL(v))])
    }
    values(x)[!is.na(values(x))] <- as.numeric(v)
  }
  if (missing(zlim)) {
    if (inherits(x,"Raster")) {
      zlim  <- range(0,values(x),finite=TRUE)
    } else if (inherits(x,"population")) {
      zlim  <- range(0,as.numeric(x$N),finite=TRUE)
    }
  }
  plot(x,nr=nr,zlim=zlim,...)
}

## have strong additive selection for A alleles
## giving you genotype dependent recruitment
## calculate a germination function?
germination.args <- list( 
  r0 = 0.01,  # one in ten seeds will germinate at low densities
  s = 1.5,    # multiplicative selective benefit of the A allele
  competition = migration(
    kern="gaussian",
    sigma=100,
    radius=300,
    normalize=1,
    #do.M=TRUE, ## not working for some reason?
    population=habitat
  )
)
germination_fun <- function (N, carrying.capacity, ...) {
  out <- germination.args$r0 / ( 1 + migrate(rowSums(N),germination.args$competition)/carrying.capacity )
  return( cbind( aa=out, aA=germination.args$s*out, AA=germination.args$s^2*out ) )
}

this.demography <- demography(
  prob.seed = 0.2,
  fecundity = 100,
  prob.germination = germination_fun,
  prob.survival = 0.9,
  pollen.migration = migration(
    kern = function (x) { exp(-sqrt(x)) },
    sigma = 100,
    radius = 1000,
    normalize = NULL,
    #do.M=TRUE, ## not working?
    population=habitat
  ),
  seed.migration = migration(
    kern = "gaussian",
    sigma = 20,
    radius = 400,
    normalize = 1,
    #do.M=TRUE, ## not working?
    population=habitat
  ),
  genotypes = c("aa","aA","AA"),
  mating = mating_tensor( c("aa","aA","AA") )
)

## don't know why do.M not working 

## give pop a few advantageous alleles
total.habitat <- sum(values(habitat),na.rm=TRUE)
pop <- population( 
  habitat = habitat,
  genotypes = c("aa","aA","AA"),
  N = cbind( aa=rpois_raster(habitat,only.values=TRUE),
             aA=rpois_raster(habitat*100/total.habitat,only.values=TRUE),
             AA=0 )
)

pl(pop)

## cehck to see if this population is stable 
exp.pop <- pop
exp.pop$N[] <- 1
demog <- setup_demography( this.demography, pop )

vitals <- generation(exp.pop,demog,
                     carrying.capacity=values(pop$habitat)[pop$habitable],
                     expected=TRUE,return.everything=TRUE)
for (k in seq_along(vitals)) {
  pl(exp.pop$habitat,v=vitals[[k]],main=c("",paste(names(vitals)[k],"density 1"),""))
}