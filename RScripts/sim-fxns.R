### Notes to self




###
### LIBRARIES
###


library(igraph)
library(NetIndices)
library(deSolve)
library(ggplot2)
library(MuMIn)
library(rootSolve)
library(DAAG)

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
    states[states < 10^-30] <- 0 
    if(sum(states >= 100) >= 1){states<-rep(0, length(states))} 
    return(c(states))
  })
}

# computes correlation between interaction strength pairs
icor <- function(x){
  i1 <- x[upper.tri(x)]
  i2 <- t(x)[upper.tri(x)] 
  
  ct <- cor.test(i1, i2)
  return(c(ct$statistic, ct$p.value))
}

# compute frequency of different interaction types for whole matrix
itypes <- function(x){
  i1 <- x[upper.tri(x)]
  i2 <- t(x)[upper.tri(x)] 
  
  comp <- sum(i1 < 0 & i2 < 0)
  mut <- sum(i1 > 0 & i2 > 0)
  pred <- sum(i1 > 0 & i2 < 0 | i1 < 0 & i2 > 0)
  amens <- sum(i1 < 0 & i2  == 0 | i1 == 0 & i2 < 0)
  comm <- sum(i1 > 0 & i2  == 0 | i1 == 0 & i2 > 0)
  
  return(c(comp = comp, mut = mut, pred = pred, amens = amens, comm = comm))
}

# compute frequency of different interaction types each spp participates in
itypes.sp <- function(x){
  mm1 <- matrix(nrow = nrow(x), ncol = 5)
  for(i in 1:nrow(x)){
    i1 <- x[i, -i]
    i2 <- x[-i, i]
    
    comp <- sum(i1 < 0 & i2 < 0)
    mut <- sum(i1 > 0 & i2 > 0)
    pred <- sum(i1 > 0 & i2 < 0 | i1 < 0 & i2 > 0)
    amens <- sum(i1 < 0 & i2  == 0 | i1 == 0 & i2 < 0)
    comm <- sum(i1 > 0 & i2  == 0 | i1 == 0 & i2 > 0)
    
    mm1[i,] <- c(comp = comp, mut = mut, pred = pred, amens = amens, comm = comm)
  }
  return(mm1)
}

# compute mean strength of each type interaction a spp participates in
istr.sp <- function(x){
  mm1 <- matrix(nrow = nrow(x), ncol = 3)
  for(i in 1:nrow(x)){
    i1 <- x[i, -i]
    i2 <- x[-i, i]
    
    comp <- which(i1 < 0 & i2 < 0)
    mut <- which(i1 > 0 & i2 > 0)
    pred <- which(i1 > 0 & i2 < 0 | i1 < 0 & i2 > 0)
    amens <- which(i1 < 0 & i2  == 0 | i1 == 0 & i2 < 0)
    comm <- which(i1 > 0 & i2  == 0 | i1 == 0 & i2 > 0)
    
    mm1[i,] <- c(comp = mean(c(i1[comp],i2[comp])), mut = mean(c(i1[mut],i2[mut])), pred = mean(c(i1[pred],i2[pred])))
    mm1[i,][is.nan(mm1[i,])] <- 0
  }
  return(mm1)
}

spints <- function(x){
  mm1 <- matrix(nrow = nrow(x), ncol = 10)
  for(i in 1:nrow(x)){
    i1 <- x[i, -i] # effect on self
    i2 <- x[-i, i] # effect on others
    
    pos1 <- sum(i1 > 0)
    neg1 <- sum(i1 < 0)
    non1 <- sum(i1 == 0)
    
    spos1 <- mean(i1[i1 > 0])
    sneg1 <- mean(i1[i1 < 0])
    
    
    pos2 <- sum(i2 > 0)
    neg2 <- sum(i2 < 0)
    non2 <- sum(i2 == 0)
    
    spos2 <- mean(i2[i2 > 0])
    sneg2 <- mean(i2[i2 < 0])
    
    
    mm1[i,] <- c(pos1, neg1, non1, spos1, sneg1, pos2, neg2, non2, spos2, sneg2)
  }
  mm1[is.nan(mm1)] <- 0
  return(mm1)
}

# function describes how species removal affects the local stability (magnitude of Max(Re(Lambda))) of the equilibrium comm
eigenkey <- function(mat, growth, isp, dyna){
  eq.biom <- dyna[1000,-1][dyna[1000,-1] > 0]
  j1 <- jacobian.full(eq.biom, lvmod, parms = list(alpha = growth[isp], m = mat[isp,isp]))
  
  ev.init <- max(Re(eigen(j1)$values))                              # initial eigenvalue for the community
  ev.key <- c()
  for(i in 1:length(isp)){
    ispR <- isp[-i]                                                 # species removed from comm
    
    j2 <- jacobian.full(eq.biom[-i], lvmod, parms = list(alpha = growth[ispR], m = mat[ispR,ispR]))
    
    ev.key[i] <- max(Re(eigen(j2)$values))                          # eigenvalue of perturbed community
  }
  return(ev.key)                                                    # return new eigenvalue (should this be the difference?)
}

# function to get network properties (for now just the network)
getgraph <- function(mat){
  mat[mat != 0] <- 1
  diag(mat) <- 0
  g <- graph.adjacency(mat)
  return(g)
}

# function to simulate impact of independent removals of all species in equilibrium comm
keystone <- function(x, dyn, eqcomm, mats, growth){
  rem <- list()                                                     # list for ode results
  pers <- c()                                                       # vector for persistence (could be done outside the loop)
  
  dyna <- dyn[[x]]                                                  # which initial community dynamics are we using
  spp1 <- eqcomm[[x]]                                               # what species from initial community are present at equilibrium
  initial1 <- mats[spp1, spp1]                                      # interaction matrix for equilibrium community
  for(i in 1:nrow(initial1)){
    
    parms <- list(alpha = growth[spp1[-i]], m = initial1[-i,-i])    # parameters of community with species i removed
    states <- dyna[1000,-1][dyna[1000,-1] > 0]                      # initial (eq) abundances from prior dynamics
    states <- states[-i]                                            # remove species i from initial (eq) abundances
    
    # simulation of new perturbed community
    rem[[i]] <- ode(states, 1:1000, parms = parms, func = lvmod, events = list(func = ext1, time =  1:1000))
    matplot(rem[[i]][,-1], typ = "l")                               # plot dynamics
    
    if(nrow(rem[[i]]) == 1000){                                     # if statement to determine if the run worked or crapped out
      pers[i] <- sum(rem[[i]][1000,-1] > 0)                         # what fraction of species have positive abundance
    }else{
      pers[i] <- NA                                                 # if run crashed gives NA
    }
    
    
    print(i)
  }
  
  # initial abundances for the removal sim
  init.biom <- dyna[1000,-1][dyna[1000,-1] > 0]
  
  # change in equilibrium abundance following removal for each species
  delta.eq <- sapply(1:length(rem), function(x) if(nrow(rem[[x]]) == 1000){rem[[x]][1000,-1] - init.biom[-x]}else{rep(NA, length(init.biom[-1]))})
  # change in mean abundance following each removal
  #delta.biom <- sapply(rem, function(x) if(nrow(x) == 1000){mean(x[1000,-1] - mean(init.biom))}else{NA})
  delta.biom <- colMeans(delta.eq)
  
  # get coefficient of variation for each species following each removal (last 200 timesteps)
  vary <- lapply(rem, function(x) if(nrow(x) == 1000){apply(x[800:1000,-1], 2, function(y) sd(y)/mean(y))}else{NA})
  # mean CV following each removal
  mean.vary <- sapply(vary, function(x){x[is.nan(x)] <- 0; mean(x)})
  # get coefficient of variation for each species following each removal (first 50 timesteps)
  init.vary <- lapply(rem, function(x) if(nrow(x) == 1000){apply(x[1:50,-1], 2, function(y) sd(y)/mean(y))}else{NA})
  # mean CV immediately following each removal
  m.init.vary <- sapply(init.vary, mean)
  
  # does each species have positive abundance at equilibrium following each removal
  # is.eq <- t(sapply(rem, function(x) if(nrow(x) == 1000){(x[1000,-1] > 0)*1}else{NA}))
  is.eq <- matrix(NA, length(spp1), length(spp1))
  for(i in 1:length(spp1)){if(nrow(rem[[i]] == 1000)){is.eq[i, -i] <- rem[[i]][1000, -1] > 0}else{is.eq[i,-i] <- rep(NA, length(spp1)[-1])}}
  
  
  # get data matrix for abundance change, variability, and persistence
  dat <- cbind(delta.biom, mean.vary, m.init.vary, pers)  
  
  return(list(dat, t(delta.eq), is.eq))
}


plotCI <- function(fit, confidence = 0.95){
  ciA <- confint(fit, level = confidence)
  modat <- data.frame(ciA[complete.cases(ciA),], coeff <- coef(summary(fit))[,1], rownames(ciA)[complete.cases(ciA)], pval = coef(summary(fit))[,4])
  colnames(modat) <- c("lower", "upper", "coeff", "met", "pval")
  modat$cols <- factor("0.1 < p", levels = c("p <= 0.05", "0.05 < p <= 0.1", "0.1 < p"))
  modat$cols[modat$pval <= 0.05] <- factor("p <= 0.05", levels = c("p <= 0.05", "0.05 < p <= 0.1", "0.1 < p"))
  modat$cols[modat$pval > 0.05 & modat$pval <= 0.1] <- factor("0.05 < p <= 0.1", levels = c("p <= 0.05", "0.05 < p <= 0.1", "0.1 < p"))
  modat$cols[modat$pval > 0.1] <- factor("0.1 < p", levels = c("p <= 0.05", "0.05 < p <= 0.1", "0.1 < p"))
  
  
  ggplot(modat) + geom_segment(aes(x = lower, y = met, xend = upper, yend = met, col = cols)) + geom_vline(aes(xintercept = 0)) + 
    geom_point(aes(x = coeff, y = met, col = cols)) + xlab("Coefficient") + ylab("Variable") + 
    scale_color_manual(name = "Significance", limits = levels(cols), values = c("blue", "darkgreen", "grey"), drop = F) + 
    theme_bw()
}