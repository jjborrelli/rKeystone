
###
### LIBRARIES
###


library(igraph)
library(deSolve)
library(ggplot2)
library(reshape2)

###
### FUNCTIONS
###


# Basic Lotka-Volterra model
lvmod <- function(times, state, parms){
  with(as.list(c(state, parms)), {
    dB <- state * parms$alpha + state * parms$m %*% state 
    
    list(dB)
  })
}

# Function to detect extinction (prevents negative abundances)
ext1 <- function (times, states, parms){
  with(as.list(states), {
    states[states < 10^-5] <- 0 
    #if(sum(states >= 100) >= 1){states<-rep(0, length(states))} 
    return(c(states))
  })
}

make_community <- function(S, C, mu.ints, sd.ints){
  growth <- runif(S, .01, 1)     
  mats <- get.adjacency(erdos.renyi.game(S, C, "gnp", directed = F), sparse = F)  
  # note: altered to directed net, .15 = C
  nINT <- sum(mats)                                                  
  
  INTstr <- rnorm(nINT, mu.ints, sd.ints)

  mats[mats != 0] <- INTstr    
  mats[mats != 0] <- mats[mats != 0] * rbinom(nINT, size = 1, prob = .8)
  
  diag(mats) <- -rbeta(nrow(mats), 1.1, 5)*5
  
  return(mats)
}

sim_community <- function(times, state, parms, eq = lvmod, ex = ext1){
  out <- ode(state, times, parms = parms, func = eq, events = list(func = ex, time =  times))
  return(out[,-1])
}

################################


SteinInt <- read.csv("~/Desktop/GitHub/microbial-dyn/Data/ecomod-ints.csv", row.names = 1)
INTs <- c(SteinInt[upper.tri(SteinInt)],SteinInt[lower.tri(SteinInt)])

c1 <- make_community(50, .3, mean(INTs), sd(INTs))
iA <- runif(nrow(c1))
iP <- list(alpha = runif(nrow(c1), .1, 1), m = c1)

sc1 <- sim_community(times = seq(1,500,1), state = iA, parms = iP)
matplot(sc1, typ = "l")

scl <- list()
for(i in 1:length(which(sc1[500,] > 0))){
  iA2 <- sc1[500,]
  iA2[which(sc1[500,] > 0)[i]] <- 0
  
  scl[[i]] <- sim_community(times = seq(1,500,1), state = iA2, parms = iP)
}

spp <- which(sc1[500,] > 0)
bio.t <- sc1[500,spp]
g1 <- graph.adjacency(abs(c1), weighted = T, diag = F)

spaths <- matrix(nrow = length(spp), ncol = length(spp))
ipaths.med <- matrix(0, nrow = length(spp), ncol = length(spp))
ipaths.sum <- matrix(0, nrow = length(spp), ncol = length(spp))
ipaths.mult <- matrix(0, nrow = length(spp), ncol = length(spp))
for(i in 1:length(spp)){
  shp <- shortest_paths(g1, spp[i], spp, mode = "out", output = "epath")$epath
  spaths[i,] <- sapply(shp, length)
  shpv <- shortest_paths(g1, spp[i], spp[-i], mode = "out", output = "vpath")$vpath
  ipa <- lapply(shpv, function(x){
    e1 <- cbind(x[1:(length(x)-1)],x[2:(length(x))])
    ipa <- c()
    for(i in 1:nrow(e1)){
      ipa[i] <- c1[e1[i,1], e1[i,2]]
    }
    return(ipa)
  })
  ipaths.med[i, spp %in% spp[-i]] <- sapply(ipa, median)
  ipaths.sum[i, spp %in% spp[-i]] <- sapply(ipa, sum)
  ipaths.mult[i, spp %in% spp[-i]] <- sapply((abs(sapply(ipa, prod))), function(x) x^(length(x)))
}
spaths
ipaths
((scl[[1]][500,] -sc1[500,])/sc1[500,])[spp]

prs <- expand.grid(1:length(spp),1:length(spp))
res <- matrix(nrow = nrow(prs), ncol = 9)
for(i in 1:nrow(prs)){
  rspp <- which(spp %in% as.vector(spp)[prs[i,1]])
  if(nrow(scl[[rspp]]) < 500){
    res[i,] <- c(NA,NA,NA,NA,NA,NA,NA,NA,NA)
    next
  }
  res[i,] <- c(rem = as.vector(spp)[prs[i,1]], trgt = as.vector(spp)[prs[i,2]],
    bio.rem = as.vector(sc1[500,as.vector(spp)[prs[i,1]]]), bio.trgt = as.vector(sc1[500,as.vector(spp)[prs[i,2]]]),
    dbio.trgt = as.vector(scl[[rspp]][500,as.vector(spp)[prs[i,2]]]),
    pl = spaths[prs[i,1], prs[i,2]], imed = ipaths.med[prs[i,1], prs[i,2]], isum = ipaths.sum[prs[i,1], prs[i,2]], iprod = ipaths.mult[prs[i,1], prs[i,2]])
}
colnames(res) <- c("rem", "trgt", "bio.rem", "bio.trgt", "dbio.trgt", "pl", "imed", "isum", "iprod")
head(res)
