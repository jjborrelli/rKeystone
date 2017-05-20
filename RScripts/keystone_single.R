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
    i1 <- mat[,i]
    i2 <- mat[i,]
    
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
  
  i1 <- c()
  i2 <- c()
  
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


library(parallel)
library(doSNOW)
strt <- Sys.time()
cl <- makeCluster(detectCores() - 1)
clusterExport(cl, c("lvm", "ext1", "fill_mat", "isim", "persist", "biodiff", "cvar", "key_effect", "int_sp", "sp_role"))
registerDoSNOW(cl)
iter = 200

key.res <- foreach(x = 1:iter, .packages = c("deSolve", "rnetcarto", "igraph", "rootSolve")) %dopar% {
  init <- isim(S = 50, tf = 2000, efun = ext1, idis = "beta", dp1 = 1, dp2 = 4, Rmax = 1, self = 1, plot = TRUE)

  mat <- init$m[init$dyn1[2000,-1] > 10^-10,init$dyn1[2000,-1] > 10^-10]
  diag(mat) <- 0
  
  nc <- rnetcarto::netcarto(abs(sign(mat)))[[1]]
  nc <- nc[order(as.numeric(as.character(nc$name))),]
  
  
  ke1 <- key_effect(init, nc$module)
  sr1 <- sp_role(mat, init$dyn1)
  
  mat <- init$m[init$dyn1[2000,-1] > 10^-10,init$dyn1[2000,-1] > 10^-10]
  
  isp <- int_sp(mat)
  
  keylist <- list(effect = ke1[[1]], modeffect = ke1[[2]], roles = sr1, ints = isp)
  
  saveRDS(keylist, paste("D:/jjborrelli/keystone/", "key", x, ".rds", sep = ""))
  
  return(keylist)
}

stopCluster(cl)
fin <- Sys.time()
fin - strt

ke <- do.call(rbind, lapply(key.res, "[[", 1))
ro <- do.call(rbind, lapply(key.res, "[[", 3))
ints <- do.call(rbind, lapply(key.res, "[[", 4))

hist(ke[,"tot"][ke[,"tot"]<0])
sum(ke[,"tot"] < -20, na.rm  = T)
dim(ke)
