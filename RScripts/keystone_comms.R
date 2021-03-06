#save.image("~/Desktop/keystoneDATA.rData")
#load("~/Desktop/keystoneDATA.rData")


library(igraph)
library(NetIndices)
library(deSolve)
library(ggplot2)
library(rootSolve)
library(reshape2)
library(boot)
library(data.table)

#################################################################################################
#################################################################################################
#################################################################################################

# Basic Lotka-Volterra model
lvmodK <- function(times, state, parms){
  with(as.list(c(state, parms)), {
    dB <- state * parms$alpha * (1 - state/(parms$K/sum(state > 0))) + state * parms$m %*% state 
    #dB <- state * parms$alpha * (1 - state/(parms$K)) + state * parms$m %*% state 
    
    list(dB)
  })
}

# Function to detect extinction (prevents negative abundances)
ext1 <- function (times, states, parms){
  with(as.list(states), {
    states[states < 10^-5] <- 0 
    
    return(c(states))
  })
}

# Function to get equilibrium communities
get_eq <- function(mats, times, INTs, Rmax = 1, Kval = 20){
  dyn <- list()
  mtmats <- list()
  grs <- list()
  for(i in 1:length(mats)){
    t1 <- mats[[i]]
    diag(t1) <- 0  #-rbeta(length(diag(t1)), 1.1, 5)*5
    t1[t1 == 1] <- abs(rnorm(sum(t1 == 1), 0, .5))
      #abs(rnorm(sum(t1 == 1), mean(INTs), sd(INTs))) #runif(sum(t1 == 1), 0, 1) 
    t1[t1 == -1] <- -abs(rnorm(sum(t1 == -1), 0, .5))
      #-abs(rnorm(sum(t1 == -1), mean(INTs), sd(INTs))) # runif(sum(t1 == -1), -1, 0) 
    
    gr <- runif(nrow(t1), .1, Rmax)
    parms <- list(alpha = gr, m = t1, K = Kval)
    
    # numerical integration of ODE, simulates dynamics of local community
    test <- ode(runif(nrow(t1), .1, .5), 1:times, parms = parms, func = lvmodK, events = list(func = ext1, time =  1:times))
    
    if(nrow(test) == times){
      dyn[[i]] <- test[,-1]
      mtmats[[i]] <- t1
      grs[[i]] <- gr
    }else{
      dyn[[i]] <- NA
      mtmats[[i]] <- NA
      grs[[i]] <- NA
    }
    
    matplot(test[,-1], type = "l", main = i)
    print(length(dyn[!is.na(dyn)]))
  }
  
  ncomm <- sum(!is.na(dyn))
  ## initial interaction matrices for communities that worked
  mtmats <- mtmats[!is.na(dyn)]#[sapply(dyn[!is.na(dyn)], function(x) sum(is.na(x)) == 0) & !sapply(dyn, function(x) sum(x < 0) > 0)]
  ## growth rates of spp for communities that worked
  grs <- grs[!is.na(dyn)]#[sapply(dyn[!is.na(dyn)], function(x) sum(is.na(x)) == 0)]
  ## dynamics of communities that worked
  dyn <- dyn[!is.na(dyn)]#[sapply(dyn[!is.na(dyn)], function(x) sum(is.na(x)) == 0)]
  ## species with positive biomass
  spp1 <- lapply(dyn, function(x) (x[times,] != 0))
  ## equilibrium matrices
  eqmat <- lapply(1:ncomm, function(x) mtmats[[x]][spp1[[x]], spp1[[x]]])
  eqgrs <- lapply(1:ncomm, function(x) grs[[x]][spp1[[x]]])
  eqst <- lapply(1:ncomm, function(x) dyn[[x]][times, spp1[[x]]])
  
  
  return(list(spp = spp1, eqm = eqmat, eqgr = eqgrs, eqst = eqst))
}


remove_all <- function(eqcomm, tmax){
  eqgrs <- eqcomm$eqgr
  eqmat <- eqcomm$eqm
  eqst <- eqcomm$eqst
  
  spr <- list()
  for(comm in 1:length(eqmat)){
    par <- list(alpha = eqgrs[[comm]], m = eqmat[[comm]], K = 20)
    st1 <- eqst[[comm]]
    
    dflist <- list()
    for(i in 1:nrow(par$m)){
      dflist[[i]] <- remove.sp(i, parms = par, states = st1, times = tmax)
      print(i)
    }
    
    dflist <- dflist[!is.na(dflist)]
    spr[[comm]] <- data.frame(comm = comm, do.call(rbind, dflist))
    
    cat("---------------------|| ", comm, " ||---------------------", "\n")
  }
  
  sprl <- rbindlist(spr)
  return(sprl)
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

coefvar <- function(dyna){
  cv10 <- apply(dyna[1:5,], 2, function(x) sd(x)/mean(x))
  cv100 <- apply(dyna[1:100,], 2, function(x) sd(x)/mean(x))
  cv1k <- apply(dyna[1:1000,], 2, function(x) sd(x)/mean(x))
  cvfi <- apply(tail(dyna, 100), 2, function(x) sd(x)/mean(x))
  
  return(data.frame(cv10 = cv10, cv100 = cv100, cv1k = cv1k, cvfi100 = cvfi))
}

remove.sp <- function(sp, parms, states, times = 2000){
  rbio <- states[sp]
  rbio.rel <- states[sp]/sum(states)
  states[sp] <- 0
  test <- ode(states, 1:times, parms = parms, func = lvmodK, events = list(func = ext1, time =  1:times))
  print(matplot(test[,-1], typ = "l"))
  
  if(nrow(test) == times){
    bi <- betweenness(graph.adjacency(abs(sign(parms$m))))
    br <- rep(NA, length(bi))
    br[-sp] <- betweenness(graph.adjacency(abs(sign(parms$m[-sp,-sp]))))
    bf <- rep(NA, length(bi))
    bf[test[times,-1]!=0] <- betweenness(graph.adjacency(abs(sign(parms$m[test[times,-1]!=0,test[times,-1]!=0]))))
    
    di <- degree(graph.adjacency(abs(sign(parms$m))))
    dr <- rep(NA, length(di))
    dr[-sp] <- degree(graph.adjacency(abs(sign(parms$m[-sp,-sp]))))
    df <- rep(NA, length(di))
    df[test[times,-1]!=0] <- degree(graph.adjacency(abs(sign(parms$m[test[times,-1]!=0,test[times,-1]!=0]))))
    
    
    sppdf <- data.frame(spR = sp, spT = 1:nrow(parms$m), coefvar(test[,-1]), sprBIO = rbio, sprREL = rbio.rel, iBio = test[1,-1],
                        t10Bio = test[10,-1], t100Bio = test[100,-1], fBio = test[times,-1], iBr = test[1,-1]/sum(test[1,-1]),
                        fBr = test[times,-1]/sum(test[times,-1]), iBet = bi, rBet = br, fBet = bf, iDeg = di, rDeg = dr, fDeg = df,
                        gr = parms$alpha, self = diag(parms$m), 
                        exts = (test[times,-1] != 0))
  }else{
    sppdf <- NA
  }
  
  return(sppdf)
}


impacts <- function(ra, ge){
  persist <- aggregate(ra$exts, list(spr = ra$spR, comm = ra$comm), function(x) sum(x)/(length(x)-1))
  dbio <- aggregate(ra$fBio - ra$iBio, list(spr = ra$spR, comm = ra$comm), function(x) c(median(x[x < 0], na.rm = T),median(x[x > 0], na.rm = T)))
  
  eigs <- lapply(unique(persist$comm), function(x){
    gr <- ge$eqgr[[x]]
    m <- ge$eqm[[x]]
    st1 <- ge$eqst[[x]]
    mre <- c(NA)
    mim <- c(NA)
    for(i in persist$spr[persist$comm == x]){
      jf <- jacobian.full(st1[-i], func = lvmodK, parms = list(alpha = gr[-i], m = m[-i,-i], K = 20))
      mre[i] <- max(Re(eigen(jf)$values))
      mim[i] <- Im(eigen(jf)$values)[which.max(Re(eigen(jf)$values))]
    }
    return(data.frame(comm = x, sp = persist$spr[persist$comm == x], mre = mre[!is.na(mre)], mim = mim[!is.na(mre)]))
  })
  
  degs <- aggregate(ra$iDeg, list(spr = ra$spR, comm = ra$comm), function(x) x)
  degs <- apply(degs, 1, function(x) x$x[x$spr[1]])
  
  rbio <- aggregate(ra$sprBIO, list(spr = ra$spR, comm = ra$comm), function(x) median(x))
  relbio <- aggregate(ra$sprREL, list(spr = ra$spR, comm = ra$comm), function(x) median(x))
  
  c10 <- aggregate(ra$cv10, list(spr = ra$spR, comm = ra$comm), function(x) median(x, na.rm = T))$x
  c100 <- aggregate(ra$cv100, list(spr = ra$spR, comm = ra$comm), function(x) median(x, na.rm = T))$x
  c1k <- aggregate(ra$cv1k, list(spr = ra$spR, comm = ra$comm), function(x) median(x, na.rm = T))$x
  cvfi <- aggregate(ra$cvfi100, list(spr = ra$spR, comm = ra$comm), function(x) median(x, na.rm = T))$x
  
  ra$fBio[ra$fBio < 10^-5] <- 0
  d1 <- aggregate(ra$iBio, list(spr = ra$spR, comm = ra$comm), function(x) vegan::diversity(x))$x
  d2 <- aggregate(ra$fBio, list(spr = ra$spR, comm = ra$comm), function(x) vegan::diversity(x))$x
  
  d1a <- aggregate(ra$iBio, list(spr = ra$spR, comm = ra$comm), function(x) vegan::diversity(x, "simpson"))$x
  d2a <- aggregate(ra$fBio, list(spr = ra$spR, comm = ra$comm), function(x) vegan::diversity(x, "simpson"))$x
  
  d1b <- aggregate(ra$iBio, list(spr = ra$spR, comm = ra$comm), function(x) vegan::diversity(x, "inv"))$x
  d2b <- aggregate(ra$fBio, list(spr = ra$spR, comm = ra$comm), function(x) vegan::diversity(x, "inv"))$x
  
  eigs <- do.call(rbind, eigs)
  
  dfall <- data.frame(persist, dbioN = dbio$x[,1], dbioP = dbio$x[,2], c10, c100, c1k, cvfi,
                      div = (d2-d1)/d1, simp = (d2a-d1a)/d1a, invs = (d2b-d1b)/d1b, 
                      degs = degs, eigre = eigs$mre, eigim = eigs$mim, rbio = rbio$x, rrel = relbio$x)
  
  return(dfall)
}

#################################################################################################
#################################################################################################
#################################################################################################
# SteinInt <- read.csv("C:/Users/jjborrelli/Desktop/GitHub/microbial-dyn/Data/ecomod-ints.csv", row.names = 1)
#SteinInt <- read.csv("~/Desktop/GitHub/microbial-dyn/Data/ecomod-ints.csv", row.names = 1)
#INTs <- c(SteinInt[upper.tri(SteinInt)],SteinInt[lower.tri(SteinInt)])

S = 100
multityp <- lapply(1:5, function(x){
  p1 <- runif(1,0,1)
  p2 <- runif(1, p1, 1)
  c1 <- runif(1, .1, .3)
  mats <- get.adjacency(erdos.renyi.game(S, c1, "gnp", directed = F), sparse = F)
  #tat <- tatoosh*sample(c(1,-1), length(tatoosh), replace = T, prob = c(p1,1-p1))
  tat <- mats*sample(c(-1,1,0), length(mats), replace = T, prob = c(p1,p2-p1,1-(p2)))
  return((tat))
})

ge.mult <- get_eq(multityp, times = 4000, INTs = INTs, Rmax = 1)
ge.mult1 <- get_eq(multityp, times = 4000, INTs = INTs, Rmax = 2)
ge.mult2 <- get_eq(multityp, times = 4000, INTs = INTs, Rmax = 3)

a <- sapply(ge.mult$eqst, length)
b <- sapply(ge.mult1$eqst, length)
d <- sapply(ge.mult2$eqst, length)
anova(lm(c(a,b,d)~c(rep("a",length(a)),rep("b",length(b)),rep("d",length(d)))))
boxplot(c(a,b,d)~c(rep("a",length(a)),rep("b",length(b)),rep("d",length(d))))

sapply(ge.mult$eqst, function(x) {(20/length(x))})
sapply(ge.mult1$eqst, function(x) {(20/length(x))})
sapply(ge.mult2$eqst, function(x) {(20/length(x))})


ra.mult <- remove_all(ge.mult, tmax = 2000)

im.mult <- impacts(ra.mult, ge.mult)

pcM <- princomp(im.mult[-c(1,2,12,15,16,17)])
loadings(pcM)

dim(im.mult)

e2 <- list()
eA <-c()
for(x in 1:length(ge.comp$eqm)){

  tm <- ge.comp$eqm[[x]]
  tst <- ge.comp$eqst[[x]]
  tgr <- ge.comp$eqgr[[x]]
  
  jf1 <- jacobian.full(y = tst, func = lvmodK, parms = list(alpha = tgr, m = tm, K = 20))
  jf2 <- sign(jf1)
  nL <- sum(jf2 != 0)
  eA[x] <- max(Re(eigen(jf1)$values))
  e1 <- c()
  for(i in 1:1000){
    nI <- abs(rnorm(nL))
    jf3 <- jf2
    jf3[jf2 != 0] <- rnorm(sum(jf2 != 0), jf2[jf2 != 0], .1)
    #jf3 <- jf3 * jf2
    
    e1[i] <- max(Re(eigen(jf3)$values))
  }
  
  e2[[x]] <- (e1)
}


#################################################################################################
#################################################################################################
#################################################################################################

S = 100
comptyp <- lapply(1:5, function(x){
  c1 <- runif(1, .1, .3)
  mats <- get.adjacency(erdos.renyi.game(S, c1, "gnp", directed = F), sparse = F)
  tat <- mats*-1
  return((tat))
})

ge.comp <- get_eq(comptyp, times = 4000, INTs = INTs)
ra.comp <- remove_all(ge.comp)
im.comp <- impacts(ra.comp, ge.comp)

pcC <- princomp(im.comp[-c(1,2,12,15,16,17)][complete.cases(im.comp),])
loadings(pcC)

predtyp <- lapply(1:5, function(x){
  c1 <- runif(1, .1, .3)
  mats <- get.adjacency(erdos.renyi.game(S, c1, "gnp", directed = T), sparse = F)
  for(i in 1:nrow(mats)){for(j in 1:ncol(mats)){if(mats[i,j] == 1){mats[j,i] <- -1}else{next}}}
  return(mats)
})

ge.pred <- get_eq(predtyp, times = 4000, INTs = INTs)
ra.pred <- remove_all(ge.pred, tmax = 1000)
im.pred <- impacts(ra.pred, ge.pred)

pcP <- princomp(im.pred[-c(1,2,12,15,16,17)])
loadings(pcP)

mututyp <- lapply(1:5, function(x){
  c1 <- .2#runif(1, .1, .3)
  mats <- get.adjacency(erdos.renyi.game(S, c1, "gnp", directed = F), sparse = F)
  return((mats))
})

ge.mutu <- get_eq(mututyp, times = 4000, INTs = INTs*.2)
ra.mutu <- remove_all(ge.mutu)
im.mutu <- impacts(ra.mutu, ge.mutu)


tats <- lapply(1:100, function(x){
  p1 <- runif(1,0,.1)
  p2 <- runif(1, p1, 1)
  mats <- get.adjacency(erdos.renyi.game(S, .2, "gnp", directed = F), sparse = F)
  #tat <- tatoosh*sample(c(1,-1), length(tatoosh), replace = T, prob = c(p1,1-p1))
  tat <- mats*sample(c(0,1,-1), length(mats), replace = T, prob = c(p1,p2-p1,1-(p2)))
  return((tat))
})

ge.tat <- get_eq(tats, times = 2000, INTs = INTs)
ra.tat <- remove_all(ge.tat, tmax = 1000)
im.tat <- impacts(ra.tat, ge.tat)
#################################################################################################
#################################################################################################
#################################################################################################

testA <- rbind(cbind(im.mult, typ = "mult"), cbind(im.comp, typ = "comp"), cbind(im.pred, typ = "pred"))
pcA <- princomp(testA[-c(1,2,10,13,14)])
summary(pcA)
loadings(pcA)
plot(pcA$scores[,1]~testA$rbio)
plot(pcA$scores[,1]~testA$rrel)
cor.test(pcA$scores[,1],testA$rrel)
plot(pcA$scores[,1:2])

ggplot(data.frame(pc1 = pcA$scores[,1], rela = testA$rrel, typ = testA$typ), aes(x = rela, y = pc1)) + geom_point(size = 3, aes(col = typ)) + facet_grid(~typ)

ggplot(testA, aes(x = dbioN, y = dbioP)) + geom_point() + facet_grid(~typ)
ggplot(testA, aes(x = dbioN, y = dbioP)) + geom_point(aes(col = typ)) 
