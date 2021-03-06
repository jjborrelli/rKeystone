library(igraph)
library(deSolve)
library(ggplot2)
library(rootSolve)


###########################################
###########################################
###########################################
###########################################


# Lotka-Volterra model 
lvm <- function(times, state, parms){
  with(as.list(c(state, parms)), {
    dB <- (state * parms$alpha) + (state * parms$m %*% state) 
    list(dB)
  })
}

lvm2 <- function(times, state, parms){
  with(as.list(c(state, parms)), {
    dB <- (state * parms$alpha) - (parms$aii * state) + (state * parms$m %*% state) 
    list(dB)
  })
}

# Extinction function 1: Lower threshold
ext1 <- function (times, states, parms){
  with(as.list(states), {
    states[states < 10^-10] <- 0 
    return(c(states))
  })
}

# Extinction function 1: Lower threshold + stochasticity
ext2 <- function(times, states, parms){
  with(as.list(states), {
    states[states < 10^-5] <- 0 
    states[states > 0] <- abs(rnorm(1:sum(states > 0), states[states > 0], .01))
    return(c(states))
  })
}

# Fill interaction matrices with strengths
fill_mat <- function(mat, dis, p1 = .5, p2 = 1){
  t1 <- mat
  diag(t1) <- 0  
  
  if(dis == "beta"){
    t1[t1 != 0] <- rbeta(sum(t1 != 0), p1, p2) * sign(mat[mat != 0])
  }else if(dis == "unif"){
    t1[t1 == 1] <- runif(sum(t1 == 1), 0, p1)
    t1[t1 == -1] <- runif(sum(t1 == -1), -p2, 0) 
  }else if(dis == "norm"){
    t1[t1 == 1] <- abs(rnorm(sum(t1 == 1), 0, p1))
    t1[t1 == -1] <- -abs(rnorm(sum(t1 == -1), 0, p2))
  }else{
    t1[t1 == 1] <- runif(sum(t1 == 1), 0, 1)
    t1[t1 == -1] <- runif(sum(t1 == -1), -1, 0)
  }
   
  return(t1)
}

# Simulate dynamics
## returns connected community and deterministic and stochastic dynamics
isim <- function(S, tf, efun = ext1, idis = "beta", dp1 = 1, dp2 = 4, Rmax = 1, self = 1, plot = FALSE){
  cond <- FALSE
  while(!cond){
    
    p1 <- runif(1,0,1)
    p2 <- runif(1, p1, 1)
    c1 <- runif(1, .1, .8)
    
    mats <- get.adjacency(erdos.renyi.game(S, c1, "gnp", directed = F), sparse = F)
    
    multityp <- mats*sample(c(-1,1,0), length(mats), replace = T, prob = c(p1,p2-p1,1-(p2)))
    
    multityp.fill <- fill_mat(multityp, dis = idis, p1 = dp1, p2 = dp2)
    #diag(multityp.fill) <- 0
    #self <- runif(length(diag(multityp.fill)), 0, 1)
    diag(multityp.fill) <- runif(length(diag(multityp.fill)), -self, 0)
    a.i <- runif(nrow(multityp.fill), .1, .5)
    
    par1 <- list(alpha = runif(nrow(multityp.fill), 0, Rmax), m = multityp.fill)
    dyn <-(ode(a.i, times = 1:tf, func = lvm, parms = par1, events = list(func = efun, time =  1:tf)))
    
    if(any(is.na(dyn))){cond <- FALSE;next}
    
    if(nrow(dyn) == tf){
      spp <- dyn[tf, -1] > 10^-10
      conn <- is.connected(graph.adjacency(abs(sign(multityp.fill[spp,spp]))))
      cond <- ifelse(sum(spp) >= 25, conn, FALSE)
    }else{
      cond <- FALSE
    }
    
  }
  
  if(plot){
    matplot(dyn[,-1], typ = "l", main = sum(dyn[tf,-1] > 10^-10))
  }
  
  return(list(m = multityp.fill, dyn1 = dyn, alpha = par1$alpha))
}

# Simulate dynamics for tatoosh interaction web
## returns connected community and deterministic and stochastic dynamics
isim_tatoosh <- function(tf, efun = ext1, idis = "beta", dp1 = 1, dp2 = 4, Rmax = 1, self = 1, plot = FALSE, mats = tatoosh){
  cond <- FALSE
  while(!cond){
    
    multityp <- mats
    
    multityp.fill <- fill_mat(multityp, dis = idis, p1 = dp1, p2 = dp2)
    #diag(multityp.fill) <- 0
    #self <- runif(length(diag(multityp.fill)), 0, 1)
    diag(multityp.fill) <- runif(length(diag(multityp.fill)), -self, 0)
    a.i <- runif(nrow(multityp.fill), .1, .5)
    
    par1 <- list(alpha = runif(nrow(multityp.fill), 0, Rmax), m = multityp.fill)
    dyn <-evalWithTimeout(ode(a.i, times = 1:tf, func = lvm, parms = par1, events = list(func = efun, time =  1:tf)), timeout = 150, onTimeout = "silent")
    
    if(is.null(dyn)){cond <- FALSE;next}
    if(any(is.na(dyn))){cond <- FALSE;next}
    
    if(nrow(dyn) == tf){
      spp <- dyn[tf, -1] > 10^-10
      conn <- is.connected(graph.adjacency(abs(sign(multityp.fill[spp,spp]))))
      cond <- ifelse(sum(spp) >= 25, conn, FALSE)
      
    }else{
      cond <- FALSE
    }
    
  }
  
  if(plot){
    matplot(dyn[,-1], typ = "l", main = sum(dyn[tf,-1] > 10^-10))
  }
  
  return(list(m = multityp.fill, dyn1 = dyn, alpha = par1$alpha))
}


isim_competition <- function(S, tf, efun = ext1, idis = "beta", dp1 = 1, dp2 = 4, Rmax = 1, self = 1, plot = FALSE){
  cond <- FALSE
  while(!cond){
    c1 <- runif(1, .1, .8)
    
    mats <- get.adjacency(erdos.renyi.game(S, c1, "gnp", directed = F), sparse = F)
    
    multityp <- mats*-1
    
    multityp.fill <- fill_mat(multityp, dis = idis, p1 = dp1, p2 = dp2)
    #diag(multityp.fill) <- 0
    #self <- runif(length(diag(multityp.fill)), 0, 1)
    diag(multityp.fill) <- runif(length(diag(multityp.fill)), -self, 0)
    a.i <- runif(nrow(multityp.fill), .1, .5)
    
    par1 <- list(alpha = runif(nrow(multityp.fill), 0, Rmax), m = multityp.fill)
    dyn <-(ode(a.i, times = 1:tf, func = lvm, parms = par1, events = list(func = efun, time =  1:tf)))
    
    if(any(is.na(dyn))){cond <- FALSE;next}
    
    if(nrow(dyn) == tf){
      spp <- dyn[tf, -1] > 10^-10 
      spp[spp] <- spp[spp] & colSums(abs(sign(par1$m[spp, spp])))-1 != 0 & rowSums(abs(sign(par1$m[spp, spp])))-1 != 0
      
      if(sum(spp) == 0){next}
      conn <- TRUE#is.connected(graph.adjacency(abs(sign(multityp.fill[spp,spp]))))
      
      cond <- ifelse(sum(spp) >= 25, conn, FALSE)
      
    }else{
      cond <- FALSE
    }
    
  }
  
  if(plot){
    matplot(dyn[,-1], typ = "l", main = sum(dyn[tf,-1] > 10^-10))
  }
  
  return(list(m = multityp.fill[spp,spp], dyn1 = dyn[,c(TRUE, spp)], alpha = par1$alpha[spp]))
}


isim_predation <- function(S, tf, efun = ext1, idis = "beta", dp1 = 1, dp2 = 4, Rmax = 1, self = 1, plot = FALSE){
  cond <- FALSE
  while(!cond){
    c1 <- runif(1, .1, .8)
    
    mats <- get.adjacency(erdos.renyi.game(S, c1, "gnp", directed = F), sparse = F)
    multityp <- mats
    for(i in 1:nrow(multityp)){
      s1 <- sample(c(1, -1), length(i:ncol(multityp)), replace = T)
      multityp[i, i:ncol(multityp)] <- s1*multityp[i, i:ncol(multityp)]
      multityp[i:ncol(multityp), i] <- s1*-1*multityp[i:ncol(multityp),i]
    }
    
    multityp.fill <- fill_mat(multityp, dis = idis, p1 = dp1, p2 = dp2)
    #diag(multityp.fill) <- 0
    #self <- runif(length(diag(multityp.fill)), 0, 1)
    diag(multityp.fill) <- runif(length(diag(multityp.fill)), -self, 0)
    a.i <- runif(nrow(multityp.fill), .1, .5)
    
    par1 <- list(alpha = runif(nrow(multityp.fill), 0, Rmax), m = multityp.fill)
    dyn <-(ode(a.i, times = 1:tf, func = lvm, parms = par1, events = list(func = efun, time =  1:tf)))
    
    if(any(is.na(dyn))){cond <- FALSE;next}
    
    if(nrow(dyn) == tf){
      spp <- dyn[tf, -1] > 10^-10 
      conn <- is.connected(graph.adjacency(abs(sign(multityp.fill[spp,spp]))))
      cond <- ifelse(sum(spp) >= 25, conn, FALSE)
      
    }else{
      cond <- FALSE
    }
    
  }
  
  if(plot){
    matplot(dyn[,-1], typ = "l", main = sum(dyn[tf,-1] > 10^-10))
  }
  
  return(list(m = multityp.fill, dyn1 = dyn, alpha = par1$alpha))
}
###########################################
###########################################
###########################################
###########################################


persist <- function(dyn, tf, alp, mat){
  p <- sum(dyn[1,-1] > 10^-10) - sum(dyn[tf,-1] > 10^-10)
  e_tr <- apply(dyn[,-1], 1, function(x){
    p1 <- list(alpha = alp[x>10^-10], m = mat[x>10^-10,x>10^-10])
    return(max(Re(eigen(jacobian.full(x[x>10^-10], func = lvm, parms = p1))$values)))
  })
  et <- sum(e_tr > 0)
  ev <- e_tr[tf]
  return(c(ext = p, tteq = et, ls = ev))
}



biodiff <- function(dyn, tf){
  ch <- dyn[tf,-1] - dyn[1,-1]
  totalbd <- (sum(dyn[tf,-1]) - sum(dyn[1,-1]))/sum(dyn[1,-1])*100
  npos <- sum(ch > 0)
  mpos <- mean(ch[ch > 0])
  mp <- max(ch[ch > 0])
  nneg <- sum(ch < 0)
  mneg <- mean(ch[ch < 0])
  mn<- min(ch[ch < 0])
  return(c(tot = totalbd, npos = npos, mpos = mpos, maxp = mp, nneg = nneg, mneg = mneg, maxn = mn))
}

cvar <- function(dyn, ntimes = 100){
  cvi <- apply(dyn[1:ntimes,-1], 2, function(x) sd(x)/mean(x))
  cvf <- apply(dyn[(nrow(dyn) - ntimes):nrow(dyn),-1], 2, function(x) sd(x)/mean(x))
  diff <- cvf - cvi
  return(cbind(cvi, cvf, diff))
}

###########################################
###########################################
###########################################
###########################################

key_effect <- function(init, mod = NULL,  plots = TRUE){ 
  nt <- nrow(init$dyn1)
  spp <- init$dyn1[nt,-1] > 10^-10
  mat <- init$m[spp,spp]
  alp <- init$alpha[spp]
  ia <- init$dyn1[nt,-1][init$dyn1[nt,-1] > 10^-10]
  pars <- list(alpha = alp, m = mat) 
  
  div <- c()
  mmch <- list()
  per <- matrix(nrow = nrow(mat), ncol = 3)
  colnames(per) <- c("ext", "tteq", "ls")
  bdiff <- matrix(nrow = nrow(mat), ncol = 7)
  colnames(bdiff) <- c("tot", "npos", "mpos", "maxp", "nneg", "mneg", "maxn")
  cvx <- list()
  
  for(i in 1:nrow(mat)){
    ia <- init$dyn1[nt,-1][spp]
    ia[i] <- 0
    
    dyn <-(ode(ia, times = 1:nt, func = lvm, parms = pars, events = list(func = ext1, time =  1:nt)))
    
    if(nrow(dyn) != nt){
      cvx[[i]] <- rbind(rep(NA, 3), rep(NA, 3))
      next
    }
    
    if(plots){
      matplot(dyn[,-1], typ = "l", main = i)
    }
    
    per[i,] <- persist(dyn, nt, alp, mat)
    bdiff[i,] <- biodiff(dyn, nt)
    cvx[[i]] <- cvar(dyn, 100)
    div[i] <- vegan::diversity(dyn[nt,-1][dyn[nt,-1] > 10^-10])
    
    if(!is.null(mod)){
      ch <- dyn[nt,-1][dyn[1,-1] > 0] - dyn[1,-1][dyn[1,-1] > 0]
      mmch[[i]] <- cbind(aggregate(ch, list(mod[-i]), mean), rmod = mod[i], node = i)
    }
  }
  cova <- t(sapply(cvx, function(x) apply(x, 2, median, na.rm = T)))
  impact <- cbind(per, bdiff, cova, div)
  
  if(!is.null(mod)){
    mod.imp <- do.call(rbind, mmch)
    return(list(impact, mod.imp))
  }
  
  return(impact)
}

###########################################
###########################################
###########################################
###########################################

# info on species level interactions
int_sp <- function(mat){
  d1 <- diag(mat)
  diag(mat) <- 0
  mm1 <- matrix(nrow = nrow(mat), ncol = 18)
  colnames(mm1) <- c("self", "nComp", "CompIn", "CompOut", "nMut", "MutIn", "MutOut", "nPred", "PredIn", "PredOut", "nAmens", "AmensIn", "AmensOut", "nComm", "CommIn", "CommOut", "allIn", "allOut")
  for(i in 1:nrow(mat)){
    i2 <- mat[,i]
    i1 <- mat[i,]
    
    comp <- which(i1 < 0 & i2 < 0)
    mut <- which(i1 > 0 & i2 > 0)
    pred <- which(i1 > 0 & i2 < 0 | i1 < 0 & i2 > 0)
    amens <- which(i1 < 0 & i2  == 0 | i1 == 0 & i2 < 0)
    comm <- which(i1 > 0 & i2  == 0 | i1 == 0 & i2 > 0)
    
    mm1[i, "nComp"] <- length(comp)
    mm1[i, "nMut"] <- length(mut)
    mm1[i, "nPred"] <- length(pred)
    mm1[i, "nAmens"] <- length(amens)
    mm1[i, "nComm"] <- length(comm)
    
    mm1[i, "CompIn"] <- ifelse(is.na(mean(i1[comp])), 0, mean(i1[comp]))
    mm1[i, "MutIn"] <- ifelse(is.na(mean(i1[mut])), 0, mean(i1[mut]))
    mm1[i, "PredIn"] <- ifelse(is.na(mean(i1[pred])), 0, mean(i1[pred]))
    mm1[i, "AmensIn"] <- ifelse(is.na(mean(i1[amens])), 0, mean(i1[amens]))
    mm1[i, "CommIn"] <- ifelse(is.na(mean(i1[comm])), 0, mean(i1[comm]))
    
    mm1[i, "CompOut"] <- ifelse(is.na(mean(i2[comp])), 0, mean(i2[comp]))
    mm1[i, "MutOut"] <- ifelse(is.na(mean(i2[mut])), 0, mean(i2[mut]))
    mm1[i, "PredOut"] <- ifelse(is.na(mean(i2[pred])), 0, mean(i2[pred]))
    mm1[i, "AmensOut"] <- ifelse(is.na(mean(i2[amens])), 0, mean(i2[amens]))
    mm1[i, "CommOut"] <- ifelse(is.na(mean(i2[comm])), 0, mean(i2[comm]))
    
    mm1[i, "self"] <- d1[i]
    mm1[i, "allIn"] <- ifelse(sum(i1 != 0) > 0, mean(i1[i1 != 0]), 0)
    mm1[i, "allOut"] <- ifelse(sum(i2 != 0) > 0, mean(i2[i2 != 0]), 0)
  }
  
  
  return(mm1)
}

# network level data

sp_role <- function(mat, dyn){
  # convert to unweighted and weighted graphs 
  g <- graph.adjacency(abs(sign(mat)))
  g2 <- graph.adjacency(abs(mat), weighted = T)
  # unweighted and weighted betweenness
  b.uw <- betweenness(g)
  b.w <- betweenness(g2)
  # degree
  d.in <- degree(g, mode = "in")
  d.out <- degree(g, mode = "out")
  d.tot <- degree(g, mode = "total")
  # clustering
  cc.uw <- transitivity(g2, "local")
  cc.w <- transitivity(g2, "weighted")
  # average path length
  apl.uw.mean <- colMeans(distances(g))
  #apl.uw.median <- apply(distances(g), 2, median)
  apl.w.mean <- colMeans(distances(g2))
  #apl.w.median <- apply(distances(g2), 2, median)
  # biomass
  bio <- dyn[nrow(dyn), -1][dyn[nrow(dyn), -1] > 0]
  cvbio <- apply(dyn[(nrow(dyn)-10):(nrow(dyn)),-1][,dyn[nrow(dyn), -1] > 0], 2, function(x) sd(x)/mean(x))
  # modularity
  mod <- (rnetcarto::netcarto(abs(sign(mat)))[[1]])
  mod <- mod[order(as.numeric(as.character(mod$name))),]
  wc.uw <- walktrap.community(g)
  wc.mod.uw <- wc.uw$modularity
  wc.mem.uw <- wc.uw$membership
  wc.w <- walktrap.community(g2)
  wc.mod.w <- wc.w$modularity
  wc.mem.w <- wc.w$membership
  
  #res <- matrix(c(b.uw, b.w, d.in, d.out, d.tot, cc.uw, cc.w, apl.uw.mean, apl.uw.median, apl.w.mean, apl.w.median),
  #              nrow = nrow(mat), ncol = 11)
  res <- matrix(c(b.uw, b.w, d.in, d.out, d.tot, cc.uw, cc.w, apl.uw.mean, apl.w.mean, bio, cvbio,
                  wc.mod.uw, wc.mem.uw, wc.mod.w, wc.mem.w),
                nrow = nrow(mat), ncol = 15)
  colnames(res) <- c("bet.uw", "bet.w", "d.in", "d.out", "d.tot", "cc.uw", "cc.w", 
                     "apl.uw.mu", "apl.w.mu", "bio", "cvbio", "mod.uw", "mem.uw", "mod.w", "mem.w")
  
  res <- cbind(res, mod)
  return(res)
}

# split network into component interaction webs
typ_bet <- function(mat){
  mat.alt <- mat
  
  prs <- t(combn(1:nrow(mat), 2))
  
  i2 <- c()
  i1 <- c()
  
  for(i in 1:nrow(prs)){
    i1[i] <- mat[prs[i, 1], prs[i, 2]]
    i2[i] <- mat[prs[i, 2], prs[i, 1]]
    if(i1[i] < 0 & i2[i] < 0){
      mat.alt[prs[i, 1], prs[i, 2]] <- "comp"
      mat.alt[prs[i, 2], prs[i, 1]] <- "comp"
    }
    if(i1[i] > 0 & i2[i] > 0){
      mat.alt[prs[i, 1], prs[i, 2]] <- "mut"
      mat.alt[prs[i, 2], prs[i, 1]] <- "mut"
    }
    if(i1[i] > 0 & i2[i] < 0 | i1[i] < 0 & i2[i] > 0){
      mat.alt[prs[i, 1], prs[i, 2]] <- "pred"
      mat.alt[prs[i, 2], prs[i, 1]] <- "pred"
    }
    if(i1[i] < 0 & i2[i]  == 0 | i1[i] == 0 & i2[i] < 0){
      mat.alt[prs[i, 1], prs[i, 2]] <- "amens"
      mat.alt[prs[i, 2], prs[i, 1]] <- "amens"
    }
    if(i1[i] > 0 & i2[i]  == 0 | i1[i] == 0 & i2[i] > 0){
      mat.alt[prs[i, 1], prs[i, 2]] <- "comm"
      mat.alt[prs[i, 2], prs[i, 1]] <- "comm"
    }
    
  }
  
  m <- ifelse(mat.alt == "comp", 1, 0)
  b.comp <- betweenness(graph.adjacency(m))
  
  m <- ifelse(mat.alt == "mut", 1, 0)
  b.mut <- betweenness(graph.adjacency(m))
  
  m <- ifelse(mat.alt == "pred", 1, 0)
  b.pred <- betweenness(graph.adjacency(m))
  
  m <- ifelse(mat.alt == "amens", 1, 0)
  b.amens <- betweenness(graph.adjacency(m))
  
  m <- ifelse(mat.alt == "comm", 1, 0)
  b.comm <- betweenness(graph.adjacency(m))
  
  m <- matrix(c(b.comp, b.mut, b.pred, b.amens, b.comm), nrow = nrow(mat), ncol = 5)
  colnames(m) <- c("comp", "mut", "pred", "amens", "comm")
  return(m)
}

###########################################
###########################################
###########################################
###########################################
tatoosh <- as.matrix(read.csv("C:/Users/jjborrelli/Desktop/GitHub/rKeystone/tatoosh.csv", header = F))

library(parallel)
library(doSNOW)
strt <- Sys.time()
cl <- makeCluster(detectCores() - 1)
clusterExport(cl, c("lvm", "ext1", "fill_mat", "isim", "persist", "biodiff", "cvar", "key_effect", "int_sp", "sp_role", "tatoosh"))
registerDoSNOW(cl)
iter = 1000

foreach(x = 1:iter, .packages = c("deSolve", "rnetcarto", "igraph", "rootSolve", "R.utils")) %dopar% {
  init <- isim(S = 50, tf = 2000, efun = ext1, idis = "beta", dp1 = 1, dp2 = 4, Rmax = 1, self = 1, plot = TRUE)
  #init <- isim_tatoosh(tf = 2000, efun = ext1, idis = "beta", dp1 = 1, dp2 = 4, Rmax = 1, self = 1, plot = TRUE, mats = tatoosh)
  #init <- isim_competition(S = 80, tf = 2000, efun = ext1, idis = "beta", dp1 = 1, dp2 = 4, Rmax = 1, self = 1, plot = TRUE)
  #init <- isim_predation(S = 50, tf = 2000, efun = ext1, idis = "beta", dp1 = 1, dp2 = 4, Rmax = 1, self = 1, plot = TRUE)
  
  mat <- init$m[init$dyn1[2000,-1] > 10^-10,init$dyn1[2000,-1] > 10^-10]
  diag(mat) <- 0
  
  nc <- rnetcarto::netcarto(abs(sign(mat)))[[1]]
  nc <- nc[order(as.numeric(as.character(nc$name))),]
  
  
  ke1 <- key_effect(init, nc$module)
  sr1 <- sp_role(mat, init$dyn1)
  
  mat <- init$m[init$dyn1[2000,-1] > 10^-10,init$dyn1[2000,-1] > 10^-10]
  
  isp <- int_sp(mat)
  
  keylist <- list(effect = ke1[[1]], modeffect = ke1[[2]], roles = sr1, ints = isp)
  
  saveRDS(keylist, paste("D:/jjborrelli/keystonePRED/", "key", x, ".rds", sep = ""))
  
  #return(keylist)
}

stopCluster(cl)
fin <- Sys.time()
fin - strt

iter = 1000
t1 <- Sys.time()
keylistP <- list()
keylistC <- list()
for(x in 1:iter){
  keylistP[[x]] <- readRDS(paste("D:/jjborrelli/keystonePRED/", "key", x, ".rds", sep = ""))
  keylistC[[x]] <- readRDS(paste("D:/jjborrelli/keystoneCOMP/", "key", x, ".rds", sep = ""))
  print(x)
}
iter = 5000
t1 <- Sys.time()
keylist <- list()
for(x in 1:iter){
  keylist[[x]] <- readRDS(paste("D:/jjborrelli/keystone/", "key", x, ".rds", sep = ""))
  print(x)
}
t2 <- Sys.time()
t2-t1

ke <- do.call(rbind, lapply(keylist, "[[", 1))
keP <- do.call(rbind, lapply(keylistP, "[[", 1))
keC <- do.call(rbind, lapply(keylistC, "[[", 1))
kem <- do.call(rbind, lapply(keylist,"[[", 2))
ke <- data.frame(ke, typ = "Mixed")
keP <- data.frame(keP, typ = "Predation")
keC <- data.frame(keC, typ = "Competition")

ro <- do.call(rbind, lapply(keylist, "[[", 3))
ro <- data.frame(ro, typ = "Mixed")
roP <- do.call(rbind, lapply(keylistP, "[[", 3))
roP <- data.frame(roP, typ = "Predation")
roC <- do.call(rbind, lapply(keylistC, "[[", 3))
roC <- data.frame(roC, typ = "Competition")
ints <- do.call(rbind, lapply(keylist, "[[", 4))
ints <- data.frame(ints, typ = "Mixed")
intsP <- do.call(rbind, lapply(keylistP, "[[", 4))
intsP <- data.frame(intsP, typ = "Predation")
intsC <- do.call(rbind, lapply(keylistC, "[[", 4))
intsC <- data.frame(intsC, typ = "Competition")

g.o <- colnames(intsC)[grep("Out",colnames(intsC))]
g.i <- colnames(intsC)[grep("In",colnames(intsC))]

colnames(intsC)[c(grep("In",colnames(intsC)), grep("Out",colnames(intsC)))] <- c(g.o, g.i)

g.o <- colnames(intsP)[grep("Out",colnames(intsP))]
g.i <- colnames(intsP)[grep("In",colnames(intsP))]

colnames(intsP)[c(grep("In",colnames(intsP)), grep("Out",colnames(intsP)))] <- c(g.o, g.i)

g.o <- colnames(ints)[grep("Out",colnames(ints))]
g.i <- colnames(ints)[grep("In",colnames(ints))]

colnames(ints)[c(grep("In",colnames(ints)), grep("Out",colnames(ints)))] <- c(g.o, g.i)


commNspp <- sapply(lapply(keylist, "[[", 1), nrow)
commID <- rep(1:iter, commNspp)
commnsp <- rep(commNspp, commNspp)

commNsppC <- sapply(lapply(keylistC, "[[", 1), nrow)
commIDC <- rep(1:1000, commNsppC)
commnspC <- rep(commNsppC, commNsppC)

commNsppP <- sapply(lapply(keylistP, "[[", 1), nrow)
commIDP <- rep(1:1000, commNsppP)
commnspP <- rep(commNsppP, commNsppP)



modNum <- lapply(lapply(keylist, "[[", 3), function(x){
  nmod <- as.vector(table(x$module))
  return(nmod[x$module+1])
})

modNumC <- lapply(lapply(keylistC, "[[", 3), function(x){
  nmod <- as.vector(table(x$module))
  return(nmod[x$module+1])
})

modNumP <- lapply(lapply(keylistP, "[[", 3), function(x){
  nmod <- as.vector(table(x$module))
  return(nmod[x$module+1])
})


#fz <- list()
#for(i in 1:100){
#  fz[[i]] <- fzmod(get_abundvec(ro$bio[commID == i], 2000))
#}
#fz1 <- do.call(rbind, fz)
#fz2 <- do.call(rbind, fz)
#fz3 <- do.call(rbind, fz)
#fz4 <- do.call(rbind, fz)
ggplot(reshape2::melt(ints[,c(2,5,8,11,14)]), aes(x = value, y = ..density..)) +
  geom_histogram() + facet_wrap(~Var2, scales = "free", nrow = 5) + 
  labs(x = "Number of Interactions per Species", y = "Density") + theme_bw(base_size = 20)

ggplot(reshape2::melt(ints[,-c(1,2,5,8,11,14,17,18)]), aes(x = value, y = ..density..)) +
  geom_histogram() + facet_wrap(~Var2, scales = "free", nrow = 5) + 
  labs(x = "", y = "Density") + theme_bw()

ggplot(reshape2::melt(ro[,-c(2,3,4,7,9,13,14,15,16,17,20)]), aes(x = value, y = ..density..)) +
  geom_histogram() + facet_wrap(~variable, scales = "free", nrow = 3) + 
  labs(x = "", y = "Density") + theme_bw(base_size = 20)

qplot(data = dfall, x = abs(ext), y = ..density.., geom = "histogram", binwidth = 1, xlab = "Number of Extinctions", ylab = "Density") + theme_bw(base_size = 20) + facet_wrap(~typ)

dfall <- (data.frame(ext = c(ke$ext, keP$ext, keC$ext), typ = c(as.character(ke$ctyp),as.character(keP$ctyp),as.character(keC$ctyp))))

ggplot(dfall, aes(x = abs(ext), y = ..density..)) + geom_histogram(binwidth = 1) +
  facet_wrap(~factor(typ, levels = c("Mixed", "Predation", "Competition")), scales = "free") +
  theme_bw(base_size = 20) + labs(x = "Number of Extinctions", y = "Density")

ggplot(dfall, aes(y = abs(ext), 
                  x = factor(typ, levels = c("Mixed", "Predation", "Competition")))) +
  geom_point(alpha = .1, position = "jitter", size = 1) +
  geom_violin() +
  theme_bw(base_size = 14) + labs(y = "Time to Stable Equilibrium", x = "Community Type")

ggsave("D:/key-images/allextinction.svg", width = 5, height = 4.5)

dfall <- (data.frame(ext = c(ke$tteq, keP$tteq, keC$tteq), typ = c(as.character(ke$typ),as.character(keP$typ),as.character(keC$typ))))

ggplot(dfall, aes(x = abs(ext), y = ..density..)) + geom_histogram(binwidth = 100) +
  facet_wrap(~factor(typ, levels = c("Mixed", "Predation", "Competition")), scales = "free") +
  theme_bw(base_size = 20) + labs(x = "Time to Stable Equilibrium", y = "Density")

ggsave("D:/key-images/alltteq1.png", width = 5, height = 4.5)

qplot(x = abs(ke[,"tteq"]), y = ..density.., geom = "histogram", binwidth = 100, xlab = "Time to Equilibrium", ylab = "Density") + theme_bw(base_size = 20)

chdf <- data.frame(PercentChange = c(ke$tot, keP$tot, keC$tot),
                   SignChange = sign(c(ke$tot, keP$tot, keC$tot)),
                   Type = c(as.character(ke$typ),as.character(keP$typ),as.character(keC$typ)))
chdf <- chdf[chdf$PercentChange != 0,]
chdf <- chdf[complete.cases(chdf),]
chdf <- chdf[chdf$PercentChange < 250,]
#chdf$PercentChange <- log10(abs(chdf$PercentChange))
chdf$SignChange <- ifelse(chdf$SignChange < 0, "Decrease", "Increase")

ggplot(chdf, aes(x = PercentChange, y = ..density..)) + geom_histogram(binwidth = 20) +
  facet_grid(~factor(Type, levels = c("Mixed", "Predation", "Competition")), scales = "free") + 
  labs(x = "Percent Change in Biomass", y = "Density") + theme_bw(base_size = 20)

aggregate(c(ke$tot, keP$tot, keC$tot), list(c(as.character(ke$ctyp),as.character(keP$ctyp),as.character(keC$ctyp))), function(x) sum(x > 250, na.rm= T))

ggsave("D:/key-images/allperchange.svg", width = 15, height = 10)

ggplot(chdf, aes(x = factor(Type, levels = c("Mixed", "Predation", "Competition")))) +
  geom_bar(aes(fill = SignChange), position = "fill") +
  theme_bw(base_size = 20) + labs(x = "Community Type", y = "Count", fill = "Change in \n Biomass") 

ggsave("D:/key-images/allsignchange.svg", width = 12, height = 10)

imps <- (reshape2::melt(dplyr::select(data.frame(ke), ext, tteq, ls, tot, npos, mpos, nneg, mneg, cvi, cvf, div)))
ggplot(imps, aes(x = value, y = ..density..)) + geom_histogram() + facet_wrap(~variable, scales = "free") + theme_bw()

hist(ke[,"tot"][ke[,"tot"]<0])
sum(ke[,"tot"] < 0, na.rm  = T)
dim(ke)

colnames(ro)
ggplot(kem, aes(x = factor(Group.1), y = log10(x), fill = factor(rmod))) + geom_boxplot()

df1 <- data.frame(ro, ints, ke)
df1$ext[df1$ext < 0] <- 0
df1$ext2 <- df1$ext/commnsp
df1$cID <- commID
df1$cnsp <- commnsp
df1$mnum <- unlist(modNum)
df1$stab <- ke[,"ls"] < 0
df1$pn <- ke[,"npos"]/ke[,"nneg"]
df1 <- df1[complete.cases(df1),]
#cls <- cls[complete.cases(df1)]

dfC <- data.frame(roC, intsC, keC)
dfC$ext[dfC$ext < 0] <- 0
dfC$ext2 <- dfC$ext/commnspC
dfC$cID <- commIDC
dfC$cnsp <- commnspC
dfC$mnum <- unlist(modNumC)
dfC$stab <- keC[,"ls"] < 0
dfC$pn <- keC[,"npos"]/keC[,"nneg"]
dfC <- dfC[complete.cases(dfC),]

dfP <- data.frame(roP, intsP, keP)
dfP$ext[dfP$ext < 0] <- 0
dfP$ext2 <- dfP$ext/commnspP
dfP$cID <- commIDP
dfP$cnsp <- commnspP
dfP$mnum <- unlist(modNumP)
dfP$stab <- keP[,"ls"] < 0
dfP$pn <- keP[,"npos"]/keP[,"nneg"]
dfP <- dfP[complete.cases(dfP),]

imps2 <- reshape2::melt(dplyr::select(df1, ext, pn))

ggplot(imps2, aes(x = value, y = ..density..)) + geom_histogram(bins = 40) + facet_wrap(~variable, scales = "free") + theme_bw()

uwfit <- glm(stab ~ bet.uw + d.tot + cc.uw + apl.uw.mu + bio + cvbio + mod.uw + connectivity + mnum,
            data = df1, family = "binomial")
summary(uwfit)

wfit <- glm(stab ~ bet.w + d.tot + cc.w + apl.w.mu + bio + cvbio + mod.w + connectivity + mnum,
           data = df1, family = "binomial")
summary(wfit)

library(rpart)
library(rpart.plot)
library(randomForest)
f1 <- formula((signch) ~ bet.uw + d.tot + cc.uw + apl.uw.mu + bio + cvbio + mod.uw + connectivity + participation + self + nComp + CompIn + CompOut + nMut + MutIn + MutOut + nPred + PredIn + PredOut + nAmens + AmensIn + AmensOut + nComm + CommIn + CommOut+ cID + cnsp + mnum)
f1 <- formula(factor(cls) ~ bet.w + d.tot + cc.w + apl.w.mu + bio + cvbio + connectivity + participation + self + nComp + CompIn + CompOut + nMut + MutIn + MutOut + nPred + PredIn + PredOut + nAmens + AmensIn + AmensOut + nComm + CommIn + CommOut+ cID + cnsp + mnum)

df1$signch <- (sign(df1$tot) == -1)
df1$signch <- ifelse(df1$signch, -1, 1)
dfC$signch <- (sign(dfC$tot) == -1)
dfP$signch <- (sign(dfP$tot) == -1)
library(FFTrees)
fftree1 <- FFTrees(f1, data = df1)
cart <- rpart(f1, data = df1[-which(df1[,"bio"] > 1000),], method = "class")
plotcp(cart)
prn <- prune(cart, cp = cart$cptable[which.min(cart$cptable[,"xerror"]),"CP"])
prp(prune(cart, cp = cart$cptable[which.min(cart$cptable[,"xerror"]),"CP"]), extra = 1)
rf1 <- randomForest(f1, data = dfP, mtry = 20, ntree = 500)

cart <- rpart(f1, data = dfP, method = "class")

f2 <- formula(sqrt(abs(tot))  ~ bet.uw + d.tot + cc.uw + apl.uw.mu + bio + cvbio + mod.uw + connectivity + participation + self + nComp + CompIn + CompOut + nMut + MutIn + MutOut + nPred + PredIn + PredOut + nAmens + AmensIn + AmensOut + nComm + CommIn + CommOut+ cID + cnsp + mnum)

cart2 <- rpart(f2, data = df1[-which(df1[,"bio"] > 1000),], method = "anova")
#plotcp(cart2)
prp(prune(cart2, cp = cart2$cptable[which.min(cart2$cptable[,"xerror"]),"CP"]), extra = 1)
cart2 <- rpart(f2, data = dfP, method = "anova")

f3 <- formula(ext ~ bet.w + d.tot + cc.w + apl.w.mu + cvbio + connectivity + participation + self + nComp + CompIn + CompOut + nMut + MutIn + MutOut + nPred + PredIn + PredOut + nAmens + AmensIn + AmensOut + nComm + CommIn + CommOut + cnsp + mnum )

cart3 <- rpart(f3, data = df1[-which(df1[,"bio"] > 1000),], method = "anova")
#plotcp(cart3)
prp(prune(cart3, cp = cart3$cptable[which.min(cart3$cptable[,"xerror"]),"CP"]), extra = 1)
cart3 <- rpart(f3, data = dfP, method = "anova")

gfit <- lme4::glmer(f3, data = df1[-1000,], family = "poisson")
summary(gfit)
predict(gfit, df1[1000,])


f4 <- formula(tteq ~ bet.uw + d.tot + cc.uw + apl.uw.mu + bio + cvbio + connectivity + participation + self + nComp + CompIn + CompOut + nMut + MutIn + MutOut + nPred + PredIn + PredOut + nAmens + AmensIn + AmensOut + nComm + CommIn + CommOut+ cID + cnsp + mnum)

cart4 <- rpart(f4, data = df1, method = "anova")
plotcp(cart4)
prn <- prune(cart4, cp = cart4$cptable[which.min(cart4$cptable[,"xerror"]),"CP"])
prp(prune(cart4, cp = cart4$cptable[which.min(cart4$cptable[,"xerror"]),"CP"]), extra = 1)
summary(cart4)
cart4 <- rpart(f4, data = dfC, method = "anova")

conf.matrix <- table(factor(sign(df1$tot)), predict(prn,type="class"))
rownames(conf.matrix) <- paste("Actual", rownames(conf.matrix), sep = ":")
colnames(conf.matrix) <- paste("Pred", colnames(conf.matrix), sep = ":")
print(conf.matrix)


f5 <- formula(log10(pn) ~ bet.w + d.tot + cc.w + apl.w.mu + cvbio + connectivity + participation + self + nComp + CompIn + CompOut + nMut + MutIn + MutOut + nPred + PredIn + PredOut + nAmens + AmensIn + AmensOut + nComm + CommIn + CommOut+ cID + cnsp + mnum)

cart5 <- rpart(f5, data = df1, method = "anova")
rsq.rpart(cart5)
prp(prune(cart5, cp = cart5$cptable[which.min(cart5$cptable[,"xerror"]),"CP"]), extra = 1)


# ext + tteq + ls + tot + npos + mpos + maxp + nneg + mneg + maxn + cvi + cvf + diff + div + ext2

# bet.w + d.in + d.out + d.tot + cc.w + apl.w.mu + bio + cvbio + mod.w + mem.w + name + module + connectivity + participation + role + self + nComp + CompIn + CompOut + nMut + MutIn + MutOut + nPred + PredIn + PredOut + nAmens + AmensIn + AmensOut + nComm + CommIn + CommOut + allIn + allOut  + cID + cnsp + mnum

# bet.uw + d.in + d.out + d.tot + cc.uw + apl.uw.mu + bio + cvbio + mod.uw + name + module + connectivity + participation + role + self + nComp + CompIn + CompOut + nMut + MutIn + MutOut + nPred + PredIn + PredOut + nAmens + AmensIn + AmensOut + nComm + CommIn + CommOut + allIn + allOut + cID + cnsp + mnum

plot(df1$ext, df1$diff)

form1 <- formula(ext ~ bet.w + d.in + d.out + d.tot + cc.w + apl.w.mu + bio + cvbio + mod.w + mem.w + connectivity + participation + self + nComp + CompIn + CompOut + nMut + MutIn + MutOut + nPred + PredIn + PredOut + nAmens + AmensIn + AmensOut + nComm + CommIn + CommOut + allIn + allOut  + cID + cnsp + mnum)

lin1 <- lm(form1, data = df1, x = F, y = F, model = F)
summary(lin1)

form2 <- formula(tteq ~ bet.w + d.in + d.out + d.tot + cc.w + apl.w.mu + bio + cvbio + mod.w + mem.w + connectivity + participation + self + nComp + CompIn + CompOut + nMut + MutIn + MutOut + nPred + PredIn + PredOut + nAmens + AmensIn + AmensOut + nComm + CommIn + CommOut + allIn + allOut  + cID + cnsp + mnum)

lin2 <- lm(form2, data = df1, x = F, y = F, model = F)
summary(lin2)

form3 <- formula(ls ~ bet.w + d.in + d.out + d.tot + cc.w + apl.w.mu + bio + cvbio + mod.w + mem.w + connectivity + participation + self + nComp + CompIn + CompOut + nMut + MutIn + MutOut + nPred + PredIn + PredOut + nAmens + AmensIn + AmensOut + nComm + CommIn + CommOut + allIn + allOut  + cID + cnsp + mnum)

lin3 <- lm(form3, data = df1, x = F, y = F, model = F)
summary(lin3)

form4 <- formula(tot ~ bet.w + d.in + d.out + d.tot + cc.w + apl.w.mu + bio + cvbio + mod.w + mem.w + connectivity + participation + self + nComp + CompIn + CompOut + nMut + MutIn + MutOut + nPred + PredIn + PredOut + nAmens + AmensIn + AmensOut + nComm + CommIn + CommOut + allIn + allOut  + cID + cnsp + mnum)

lin4 <- lm(form4, data = df1, x = F, y = F, model = F)
summary(lin4)

form5 <- formula(cvi ~ bet.w + d.in + d.out + d.tot + cc.w + apl.w.mu + bio + cvbio + mod.w + mem.w + connectivity + participation + self + nComp + CompIn + CompOut + nMut + MutIn + MutOut + nPred + PredIn + PredOut + nAmens + AmensIn + AmensOut + nComm + CommIn + CommOut + allIn + allOut  + cID + cnsp + mnum)

lin5 <- lm(form5, data = df1, x = F, y = F, model = F)
summary(lin5)

dfimp <- (dplyr::select(df1, ext:ext2))
pc1 <- princomp(dfimp)
summary(pc1)

#######################

df1$hasext <- df1$ext != 0
df1$bigneg <- df1$tot < -10

form6 <- formula(hasext ~ bet.w + d.in + d.out + d.tot + cc.w + apl.w.mu + bio + cvbio + mod.w + mem.w + connectivity + participation + self + nComp + CompIn + CompOut + nMut + MutIn + MutOut + nPred + PredIn + PredOut + nAmens + AmensIn + AmensOut + nComm + CommIn + CommOut + allIn + allOut  + cID + cnsp + mnum)

glin6 <- glm(form6, data = df1, x = F, y = F, model = F, family = "binomial")
summary(glin6)

form7 <- formula(bigneg ~ bet.w + d.in + d.out + d.tot + cc.w + apl.w.mu + bio + cvbio + mod.w + mem.w + connectivity + participation + self + nComp + CompIn + CompOut + nMut + MutIn + MutOut + nPred + PredIn + PredOut + nAmens + AmensIn + AmensOut + nComm + CommIn + CommOut + allIn + allOut  + cID + cnsp + mnum)

glin7 <- glm(form7, data = df1, x = F, y = F, model = F, family = "binomial")
summary(glin7)
DAAG::cv.binary(glin7)

df1$pc1 <- pc1$scores[,1]

form8 <- formula(pc1 ~ bet.w + d.in + d.out + d.tot + cc.w + apl.w.mu + bio + cvbio + mod.w + mem.w + connectivity + participation + self + nComp + CompIn + CompOut + nMut + MutIn + MutOut + nPred + PredIn + PredOut + nAmens + AmensIn + AmensOut + nComm + CommIn + CommOut + allIn + allOut  + cID + cnsp + mnum)

lin8 <- lm(form8, data = df1, x = F, y = F, model = F)
summary(lin8)
