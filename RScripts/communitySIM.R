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
use[is.na(use)] <- FALSE
# which species are present in equilibrial communities
eqcomm <- sapply(1:sum(use), function(x) spp[use][[x]][which(r2[use][[x]][1000,-1] > 0)])
eqab <- sapply(1:sum(use), function(x) r2[use][[x]][1000,which(r2[use][[x]][1000,] > 0)][-1])
t.simend <- Sys.time()                                              # note time initial comm sim ends

m.par <- lapply(1:length(spp[use]), function(x) mats[eqcomm[[x]], eqcomm[[x]]])
gr.l <- lapply(1:length(spp[use]), function(x) growth[eqcomm[[x]]])
jacs <- lapply(1:sum(use), function(x) jacobian.full(eqab[[x]], func = lvmod, parms = list(alpha = gr.l[[x]], m = m.par[[x]])))
eigs <- sapply(jacs, function(x) max(Re(eigen(x)$values)))

jacs2 <- lapply(1:5, function(x) jacobian.full(ge.mult2$eqst[[x]], func = lvmodK2, parms = list(alpha = ge.mult2$eqgr[[x]], m = ge.mult2$eqm[[x]], K = ge.mult2$eqkv[[x]])))
eigs2 <- sapply(jacs2, function(x) max(Re(eigen(x)$values)))
hist(abs(diag(jacs2[[1]])))
hist(abs(diag(jacs[[1]])))

get_bestfit2 <- function(avec){
  fits1 <- lapply(avec, fitsad, sad = "lnorm")
  fits2 <- lapply(avec, fitsad, sad = "power")
  #fits3 <- lapply(avec, fitsad, sad = "powbend")
  fits4 <- lapply(avec, fitsad, sad = "mzsm")
  fits5 <- lapply(avec, fitsad, sad = "poilog")
  fits6 <- lapply(avec, fitsad, sad = "bs")
  fits7 <- lapply(avec, fitsad, sad = "ls")
  fits8 <- lapply(avec, fitsad, sad = "weibull")
  fit1 <- sapply(fits1, AIC)
  fit2 <- sapply(fits2, AIC)
  #fit3 <- sapply(fits3, AIC)
  fit4 <- sapply(fits4, AIC)
  fit5 <- sapply(fits5, AIC)
  fit6 <- sapply(fits6, AIC)
  fit7 <- sapply(fits7, AIC)
  fit8 <- sapply(fits8, AIC)
  #sapply(fits5, function(x) c(coef(x), AICvol = AIC(x)))
  t.wins <- table(apply(cbind(fit1, fit2, fit4, fit5, fit6, fit7, fit8), 1, 
                        function(x){c("lnorm", "power", "mzsm", "poilog", "bs", "ls", "weibull")[which.min(x)]}))
  #tAIC <- cbind(fit1, fit2, fit3, fit4, fit5, fit6, fit7, fit8)
  #colnames(tAIC) <- c("lnorm", "power", "powbend", "mzsm", "poilog", "bs", "ls", "weibull")
  #return(tAIC)
  return(t.wins)
}

get_bestfit3 <- function(avec){
  fits1 <- lapply(avec, fitsad, sad = "lnorm")
  fits2 <- lapply(avec, fitsad, sad = "power")
  #fits3 <- lapply(avec, fitsad, sad = "powbend")
  fits4 <- lapply(avec, fitsad, sad = "mzsm")
  fits5 <- lapply(avec, fitsad, sad = "poilog")
  fits6 <- lapply(avec, fitsad, sad = "bs")
  fits7 <- lapply(avec, fitsad, sad = "ls")
  fits8 <- lapply(avec, fitsad, sad = "weibull")
  fit1 <- sapply(fits1, AIC)
  fit2 <- sapply(fits2, AIC)
  #fit3 <- sapply(fits3, AIC)
  fit4 <- sapply(fits4, AIC)
  fit5 <- sapply(fits5, AIC)
  fit6 <- sapply(fits6, AIC)
  fit7 <- sapply(fits7, AIC)
  fit8 <- sapply(fits8, AIC)
  #sapply(fits5, function(x) c(coef(x), AICvol = AIC(x)))
  #t.wins <- table(apply(cbind(fit1, fit2, fit4, fit5, fit6, fit7, fit8), 1, 
  #                      function(x){c("lnorm", "power", "mzsm", "poilog", "bs", "ls", "weibull")[which.min(x)]}))
  tAIC <- cbind(fit1, fit2, fit4, fit5, fit6, fit7, fit8)
  colnames(tAIC) <- c("lnorm", "power", "mzsm", "poilog", "bs", "ls", "weibull")
  #return(tAIC)
  return(tAIC)
}

gbf2a <- get_bestfit2(lapply(eqab, get_abundvec))
gbf2a
gbf3a <- get_bestfit3(lapply(eqab, get_abundvec))

boxplot(t(apply(gbf3a, 1, function(x) x - min(x))))


