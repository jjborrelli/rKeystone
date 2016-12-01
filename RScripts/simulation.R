###
### SIMULATION
###


t.start <- Sys.time()                                               # begin timing of complete simulation run

S <- 200                                                            # total number of species in the pool 
growth <- runif(S, .01, 1)                                          # growth rate for each spp in the pool
K <- quantile(1:100, rbeta(S, 1, 2))                                # carrying capacities (not used)
mats <- get.adjacency(erdos.renyi.game(S, .2, "gnp", directed = F), sparse = F)   # create interaction matrix of all spp in pool
# note: altered to directed net, .15 = C
nINT <- sum(mats)                                                   # total number of interactions 

# using stein data
SteinInt <- read.csv("~/Desktop/GitHub/microbial-dyn/Data/ecomod-ints.csv", row.names = 1)
INTs <- c(SteinInt[upper.tri(SteinInt)],SteinInt[lower.tri(SteinInt)])
INTstr <- rnorm(nINT, mean(INTs), sd(INTs))
#INTstr <- runif(nINT, -1, 1)                                        # sample interaction strengths 

mats[mats != 0] <- INTstr                                           # fill in interaction strengths 
#diag(mats) <- -1                                                    # self limitation set to -1
diag(mats) <- -rbeta(nrow(mats), 1.1, 5)*5
## 
## Begin simulation of subsampled community dynamics

# initialize lists of species and dynamics
spp <- list()
r2 <- list()
# simulation
for(i in 1:1000){
  spp[[i]] <- c(sample(71:200, 20), sample(1:70, 30))               # sample species from the pool (two samples ensure some overlap)
  isp <- spp[[i]]                                                   # local species community
  parms <- list(alpha = growth[isp], m = mats[isp,isp], k = K[isp]) # named parameter list (growth rates and int mat; K not used in sim)
  
  # numerical integration of ODE, simulates dynamics of local community
  r2[[i]] <- ode(runif(length(isp), .1, 1), 1:1000, parms = parms, func = lvmod, events = list(func = ext1, time =  1:1000))
  matplot(r2[[i]][,-1], typ = "l")                                  # plot community dynamics
  print(i)                                                          # which iteration are we on again? 
}

# which runs did not get messed up
use <- sapply(r2, nrow) == 1000 & sapply(r2, function(x) sum(tail(x, 1)[-1] > 0) > 0)
sum(use)
median(sapply(r2[use], function(x) sum(x[1000,-1] > 0))   )           # how many spp coexisting in each local comm
# which species are present in equilibrial communities
eqcomm <- sapply(1:sum(use), function(x) spp[use][[x]][which(r2[use][[x]][1000,-1] > 0)])

# how equilibrial are the communities?
dyn <- lapply(r2[use], function(x){x[x < 0] <- 0; x})               # get new object of dynamics list for runs that worked
cv.eq <- sapply(dyn, function(x) apply(x[990:1000,-1][,x[1000,-1] > 0], 2, sd)/(colMeans(x[990:1000, -1][,x[1000,-1] > 0])))
#cv.eq[is.nan(cv.eq)] <- 0
#range(colMeans(cv.eq))
hist(unlist(cv.eq))

# matrix of species found in each local equilibrium community
# can be used to determine compositional similarity of communities
eqmat <- matrix(0, nrow = sum(use), ncol = S)                       # initialize eqmat
for(i in 1:sum(use)){
  eqmat[i,eqcomm[[i]]] <- 1                                         # if the species is present in local comm i it gets a one, 0 otherwise
}

# initial interaction matrices
inmatuse <- lapply(1:sum(use), function(x) mats[spp[use][[x]],spp[use][[x]]])
# equilibrium interaction matrices for all iterations that didn't mess up
matuse <- lapply(1:sum(use), function(i) mats[eqcomm[[i]], eqcomm[[i]]])

# histogram of resulting matrix connectances
plot(sapply(matuse, function(x) sum(x != 0)/(nrow(x)*(nrow(x)-1))),sapply(inmatuse, function(x) sum(x != 0)/(nrow(x)*(nrow(x)-1))))

# compute frequency of interaction types in each equilibrium matrix
ity <- t(sapply(matuse, itypes))
plot(ity[,1]/ity[,2])                                               # plot ratio of competition:mutualism
iity <- t(sapply(inmatuse, itypes))

ity2 <- t(apply(ity, 1, function(x) x/sum(x)))
iity2 <- t(apply(iity, 1, function(x) x/sum(x)))

plot(iity2[,1], ity2[,1])

spi <- lapply(matuse, spints)
spi2 <- do.call(rbind, spi)
es <- spi2[,5] + spi2[,4]
eo <- spi2[,10] + spi2[,9]

ityinsp <- lapply(inmatuse, itypes.sp)
ityinsp2 <- do.call(rbind, ityinsp)
summary(glm(unlist(lapply(r2[use], function(x) x[1000, -1] > 0))~ityinsp2, family = "quasibinomial"))

# get frequency of interaction types for each species in each equilibrial comm
itySP <- lapply(matuse, itypes.sp)
# get mean strength of each interaction type species participate in
istrSP <- lapply(matuse, istr.sp)

# quick test to see if interaction participation influences equilibrium abundance
summary(lm(unlist(lapply(1:sum(use), function(x) r2[use][[x]][1000,-1][r2[use][[x]][1000,-1] > 0]))~do.call(rbind, itySP)))


# how each spp removal affects eigenval of comm
eigkey <- lapply(1:sum(use), function(x) eigenkey(mat = mats, growth = growth, isp = eqcomm[[x]], dyna = dyn[[x]]))        
summary(lm(unlist(eigkey)~do.call(rbind, itySP)))                   # relationship btwn eig and interaction types
summary(lm(unlist(eigkey)~do.call(rbind, istrSP)))                  # relationship btwn eig and interaction strengths


allg <- lapply(matuse, getgraph)                                    # get the network for each local eq comm
plot(unlist(eigkey)~unlist(sapply(allg, degree)))                   # look at relationship between degree and eig

betw <- lapply(allg, betweenness)                                   # get betweenness of each node
clocent <- lapply(allg, closeness)                                  # get closeness centrality
# get neighborhood of each spp going out 2 links
g.neighbors2 <- lapply(1:length(allg), function(x){sapply(graph.neighborhood(allg[[x]], 2), function(y) length(V(y)))})
ecent <- lapply(allg, function(x) eigen_centrality(x)$vector)       # get eigenvector centrality
hscore <- lapply(allg, function(x) hub_score(x)$vector)             # get hub score
p.rank <- lapply(allg, function(x) page_rank(x)$vector)             # get page rank algo



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
