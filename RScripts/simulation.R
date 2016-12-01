###
### Source Functions
###

source("~/Desktop/GitHub/rKeystone/RScripts/sim-fxns.R")

###
### SIMULATION
###


t.start <- Sys.time()                                               # begin timing of complete simulation run

S <- 500                                                            # total number of species in the pool 
nSp <- 50                                                           # number of species in each equilibrium comm
growth <- runif(S, .01, 1)                                          # growth rate for each spp in the pool
mats <- get.adjacency(erdos.renyi.game(S, .2, "gnp", directed = F), sparse = F)   # create interaction matrix of all spp in pool
# note: altered to directed net, .15 = C
nINT <- sum(mats)                                                   # total number of interactions 

# using stein data
SteinInt <- read.csv("~/Desktop/GitHub/microbial-dyn/Data/ecomod-ints.csv", row.names = 1)
INTs <- c(SteinInt[upper.tri(SteinInt)],SteinInt[lower.tri(SteinInt)])
INTstr <- rnorm(nINT, mean(INTs), sd(INTs))
#INTstr <- runif(nINT, -1, 1)                                        # sample interaction strengths 

mats[mats != 0] <- INTstr                                            # fill in interaction strengths 
mats[mats != 0] <- mats[mats != 0] * rbinom(nINT, size = 1, prob = .8)
#diag(mats) <- -1                                                    # self limitation set to -1
diag(mats) <- -rbeta(nrow(mats), 1.1, 5)*5
## 
## Begin simulation of subsampled community dynamics

## get initial abundances
set.seed(10)
iA <- t(sapply(1:1000, function(x) runif(nSp, .1, 1)))
# initialize lists of species and dynamics
spp <- list()
r2 <- list()
# simulation
for(i in 1:1000){
  spp[[i]] <- sample(1:S, nSp)              # sample species from the pool (two samples ensure some overlap)
  isp <- spp[[i]]                                                   # local species community
  parms <- list(alpha = growth[isp], m = mats[isp,isp])             # named parameter list (growth rates and int mat)
  
  # numerical integration of ODE, simulates dynamics of local community
  r2[[i]] <- ode(iA[i,], 1:1000, parms = parms, func = lvmod, events = list(func = ext1, time =  1:1000))
  matplot(r2[[i]][,-1], typ = "l")                                  # plot community dynamics
  print(i)                                                          # which iteration are we on again? 
}

# which runs did not get messed up
use <- sapply(r2, nrow) == 1000 & sapply(r2, function(x) sum(tail(x, 1)[-1] > 0) > 0)

# which species are present in equilibrial communities
eqcomm <- sapply(1:sum(use), function(x) spp[use][[x]][which(r2[use][[x]][1000,-1] > 0)])

t.simend <- Sys.time()                                              # note time initial comm sim ends

###
### CHECK FOR KEYSTONE SPECIES
###

t.key <- Sys.time()                                                 # note time keystone simm starts



# quick test of function
#system.time(
#ks2 <- keystone(2, dyn = r2[use], eqcomm, mats)
#)


# Simulated removal for each species in each equilibrium community
ks1 <- list()
ks2 <- list()
ks3 <- list()
for(i in 1:sum(use)){
  KS <- keystone(i, dyn = r2[use], eqcomm, mats, growth)            # keystone species simulation
  ks1[[i]] <- KS[[1]]                                               # biomass variability and persistence
  ks2[[i]] <- KS[[2]]                                               # change in spp biomass with removal
  ks3[[i]] <- KS[[3]]                                               # who went extinct
  
  cat(paste("\n ------------------|   ", i, "   |------------------ \n"))
}


t.end <- Sys.time()                                                 # what time does the simulation end
t.end - t.start                                                     # total time spent simulating from start to finish
