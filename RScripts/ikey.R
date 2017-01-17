
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

get_eqcomm <- function(S, C, INTs, t1, plot = FALSE){
  cond <- FALSE
  while(!cond){
    c1 <- make_community(S, C, mean(INTs), sd(INTs))
    iA <- runif(nrow(c1))
    iP <- list(alpha = runif(nrow(c1), .1, 1), m = c1)
    
    sc1 <- sim_community(times = t1, state = iA, parms = iP)
    if(nrow(sc1) == max(t1)){cond <- TRUE}
  }
  if(plot){matplot(sc1, typ = "l")}
  return(list(comm.mat = c1, comm.dyn = sc1, init.parms = iP))
}

spp_remove <- function(sc1, iP, t1, track.progress = FALSE){
  scl <- list()
  eqS <- which(tail(sc1, 1) > 0)
  for(i in 1:length(eqS)){
    iA2 <- as.vector(tail(sc1, 1))
    iA2[eqS[i]] <- 0
    
    scl[[i]] <- sim_community(times = t1, state = iA2, parms = iP)
    if(track.progress){cat(i, "\t")}
  }
  return(scl)
}

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

path_ints <- function(c1, sc1, scl, iP){
  spp <- which(sc1[nrow(sc1),] > 0)
  g1 <- graph.adjacency(abs(c1), weighted = T, diag = F)
  odegs <- degree(g1, mode = "out")
  idegs <- degree(g1, mode = "in")
  betw <- betweenness(g1, weights = E(g1)$weight)
  gr <- iP$alpha
  itysp <- itypes.sp(c1)
  
  spaths <- matrix(nrow = length(spp), ncol = length(spp))
  ipaths.med <- matrix(0, nrow = length(spp), ncol = length(spp))
  ipaths.sum <- matrix(0, nrow = length(spp), ncol = length(spp))
  ipaths.mult <- matrix(0, nrow = length(spp), ncol = length(spp))
  ipaths.sign <- matrix(0, nrow = length(spp), ncol = length(spp))
  for(i in 1:length(spp)){
    shp <- shortest_paths(g1, spp[i], spp, mode = "out", output = "epath")$epath
    spaths[i,] <- sapply(shp, length)
    shpv <- shortest_paths(g1, spp[i], spp[-i], mode = "out", output = "vpath")$vpath
    ipa <- lapply(shpv, function(x){
      if(length(x) == 0){
        return(0)
      }else{
        e1 <- cbind(x[1:(length(x)-1)],x[2:(length(x))])
        ipa <- c()
        for(i in 1:nrow(e1)){
          ipa[i] <- c1[e1[i,1], e1[i,2]]
        }
        return(ipa)
      }
    })
    
    ipaths.med[i, spp %in% spp[-i]] <- sapply(ipa, median)
    ipaths.sum[i, spp %in% spp[-i]] <- sapply(ipa, sum)
    ipaths.mult[i, spp %in% spp[-i]] <- sapply(ipa, function(x) abs(prod(x))^(1/length(x))*sign(prod(x)))
    ipaths.sign[i, spp %in% spp[-i]] <- sapply(ipa, function(x) paste0(sign(x), collapse=","))
  
  }
  
  prs <- expand.grid(1:length(spp),1:length(spp))
  res <- matrix(nrow = nrow(prs), ncol = 20)
  ipathl <- list()
  i.ty <- list()
  for(i in 1:nrow(prs)){
    rspp <- which(spp %in% as.vector(spp)[prs[i,1]])
    if(nrow(scl[[rspp]]) < nrow(sc1)){
      res[i,] <- c(rspp, rem = as.vector(spp)[prs[i,1]], trgt = as.vector(spp)[prs[i,2]],NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
      ipathl[[i]] <- NA
      i.ty[[i]] <- rep(NA, 10)
      next
    }
    res[i,] <- c(index = rspp, rem = as.vector(spp)[prs[i,1]], trgt = as.vector(spp)[prs[i,2]],
                 bio.rem = as.vector(sc1[500,as.vector(spp)[prs[i,1]]]), bio.trgt = as.vector(sc1[500,as.vector(spp)[prs[i,2]]]),
                 dbio.trgt = as.vector(scl[[rspp]][500,as.vector(spp)[prs[i,2]]]),
                 pl = spaths[prs[i,1], prs[i,2]], imed = ipaths.med[prs[i,1], prs[i,2]], isum = ipaths.sum[prs[i,1], prs[i,2]], iprod = ipaths.mult[prs[i,1], prs[i,2]],
                 odeg.rem = odegs[as.vector(spp)[prs[i,1]]], odeg.trgt = odegs[as.vector(spp)[prs[i,2]]],
                 ideg.rem = idegs[as.vector(spp)[prs[i,1]]], ideg.trgt = idegs[as.vector(spp)[prs[i,2]]],
                 bet.rem = betw[as.vector(spp)[prs[i,1]]], bet.trgt = betw[as.vector(spp)[prs[i,2]]], 
                 r.rem = gr[as.vector(spp)[prs[i,1]]], r.trgt = gr[as.vector(spp)[prs[i,2]]], 
                 self.rem = diag(c1)[as.vector(spp)[prs[i,1]]], self.trgt = diag(c1)[as.vector(spp)[prs[i,2]]])
    ipathl[[i]] <- ipaths.sign[prs[i,1], prs[i,2]]
    i.ty[[i]] <- c(itysp[as.vector(spp)[prs[i,1]],], itysp[as.vector(spp)[prs[i,2]],])
  }
  colnames(res) <- c("index", "rem", "trgt", "bio.rem", "bio.trgt", "dbio.trgt", "pl", "imed", "isum", "iprod", "odeg.rem", "odeg.trgt", "ideg.rem", "ideg.trgt",
                     "bet.rem", "bet.trgt", "r.rem", "r.trgt", "self.rem", "self.trgt")
  resdf <- as.data.frame(res)
  resdf$isign <- unlist(ipathl)
  intstyp <- do.call("rbind", i.ty)
  colnames(intstyp) <- c("comp.rem","mut.rem","pred.rem","amens.rem","comm.rem","comp.trgt","mut.trgt","pred.trgt","amens.trgt","comm.trgt")
  resdf <- data.frame(resdf, intstyp)
  return(resdf)
}

change_in_diversity <- function(eq.init, eq.removed, meas = "shannon"){
  if(all(is.na(eq.removed))){return(NA)}
  dr <- vegan::diversity(eq.removed, index = meas)
  di <- vegan::diversity(eq.init, index = meas)
  
  deltaD <- ((dr - di)/di)*100
  return(deltaD)
}

cdiv <- function(srcoms, eqi, meas = "shannon"){
  eqsr <- lapply(srcoms, function(x){if(nrow(x) == 2000){x[2000,which(eqi != 0)]}else{NA}})
  D <- sapply(eqsr, change_in_diversity, eq.init = eqi[which(eqi != 0)], meas = meas)
  return(D)
}


##############################################################################
##############################################################################
##############################################################################


SteinInt <- read.csv("~/Desktop/GitHub/microbial-dyn/Data/ecomod-ints.csv", row.names = 1)
#SteinInt <- read.csv("C:/Users/jjborrelli/Desktop/GitHub/microbial-dyn/Data/ecomod-ints.csv", row.names = 1)
INTs <- c(SteinInt[upper.tri(SteinInt)],SteinInt[lower.tri(SteinInt)])


e1 <- get_eqcomm(S = 50, C = .05, INTs = INTs, t1 = 1:2000, plot = F)
sr1 <- spp_remove(sc1 = e1$comm.dyn, iP = e1$init.parms, t1 = 1:2000)

pi1 <- path_ints(c1 = e1$comm.mat, sc1 = e1$comm.dyn, scl = sr1, iP = e1$init.parms)
head(pi1)

summary(lm((cdiv(sr1, e1$comm.dyn[2000,], meas = "shannon")/100)~itypes.sp(e1$comm.mat)[unique(pi1$rem),]))

t0 <- Sys.time()
C.list <- rep(c(.05,.15,.25), each = 20)
D.all <- list()
i.all <- list()
pi1 <- list()
for(i in 1:60){
  e1 <- get_eqcomm(S = 50, C = C.list[i], INTs = INTs, t1 = 1:2000, plot = F)
  sr1 <- spp_remove(sc1 = e1$comm.dyn, iP = e1$init.parms, t1 = 1:2000)
  
  pi1[[i]] <- path_ints(c1 = e1$comm.mat, sc1 = e1$comm.dyn, scl = sr1, iP = e1$init.parms)
  
  pi1[[i]]$N <- ifelse(nrow(e1$comm.dyn) == 2000, sum(e1$comm.dyn[2000,] != 0), NA)
  pi1[[i]]$Cin <- ifelse(nrow(e1$comm.dyn) == 2000, C.list[i], NA)
  
  if(nrow(e1$comm.dyn) == 2000){
    c.alt <- sum(sign(abs(e1$comm.dyn[which(e1$comm.dyn[2000,] != 0), which(e1$comm.dyn[2000,] != 0)])))/(sum(e1$comm.dyn[2000,] != 0)*sum(e1$comm.dyn[2000,] != 0))
  }else{
    c.alt <- NA
  }
  pi1[[i]]$Cadj <- c.alt
  
  D.all[[i]] <- cdiv(sr1, e1$comm.dyn[2000,], meas = "shannon")
  i.all[[i]] <- itypes.sp(e1$comm.mat)[unique(pi1$rem),]
  
  print(paste(Sys.time(), "-----------------", i, "-----------------"))
}
t1 <- Sys.time()
diff1 <- ((pi1$dbio.trgt - pi1$bio.trgt)/pi1$bio.trgt)*100
diff2 <- ((pi1$dbio.trgt)/pi1$bio.trgt)
hist(diff2[pi1$pl != 0])
aggregate(diff1, list(pi1$isign), median)
boxplot((log(abs(diff1)) *sign(diff1))~pi1$isign, las = 2)

pi1$diff <- ((pi1$dbio.trgt - pi1$bio.trgt))

pi2 <- do.call("rbind", pi1)
pi2$diff <- ((pi2$dbio.trgt - pi2$bio.trgt))/pi2$bio.trgt
form1 <- diff~bio.rem+bio.trgt+pl+imed+isum+iprod+odeg.rem+odeg.trgt+ideg.rem+ideg.trgt+bet.rem+bet.trgt+r.rem+r.trgt+self.rem+self.trgt+comp.rem+mut.rem+pred.rem+amens.rem+comm.rem+comp.trgt+mut.trgt+pred.trgt+amens.trgt+comm.trgt

form2 <- diff~imed+r.rem

cart1 <- rpart::rpart(form1, data = pi2[complete.cases(pi2),][which(pi2$bio.trgt <= 101),], method = "anova")
plot(cart1)
text(cart1, use.n=TRUE, all=TRUE, cex=.8)

pi3 <- pi2[which(pi2$rem != pi2$trgt),] 
pi3$ext <-pi3$dbio.trgt == 0

form1 <- ext~bio.rem+bio.trgt+pl+imed+isum+iprod+odeg.rem+odeg.trgt+ideg.rem+ideg.trgt+bet.rem+bet.trgt+r.rem+r.trgt+self.rem+self.trgt+comp.rem+mut.rem+pred.rem+amens.rem+comm.rem+comp.trgt+mut.trgt+pred.trgt+amens.trgt+comm.trgt

cart1 <- rpart::rpart(form1, data = pi3[complete.cases(pi3),][which(pi3$bio.trgt <= 101),], method = "anova")
plot(cart1)
text(cart1, use.n=TRUE, all=TRUE, cex=.8)
