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
fill_mat <- function(mat, sdevp = .5, sdevn = 1){
  t1 <- mat
  diag(t1) <- 0  #-rbeta(length(diag(t1)), 1.1, 5)*5
  #t1[t1 == 1] <- abs(rnorm(sum(t1 == 1), 0, sdevp))#
  t1[t1 == 1] <- runif(sum(t1 == 1), 0, 1) #
 
  #t1[t1 == -1] <- -abs(rnorm(sum(t1 == -1), 0, sdevn)) #
  t1[t1 == -1] <- runif(sum(t1 == -1), -1, 0) 
  return(t1)
}

# Simulate dynamics
## returns connected community and deterministic and stochastic dynamics
isim <- function(S, tf, efun = ext1, plot = FALSE){
  cond <- FALSE
  while(!cond){
    
    p1 <- runif(1,0,1)
    p2 <- runif(1, p1, 1)
    c1 <- runif(1, .1, .8)
    
    mats <- get.adjacency(erdos.renyi.game(S, c1, "gnp", directed = F), sparse = F)
    
    multityp <- mats*sample(c(-1,1,0), length(mats), replace = T, prob = c(p1,p2-p1,1-(p2)))
    
    multityp.fill <- fill_mat(multityp, sdevp = .5, sdevn = .5)
    #diag(multityp.fill) <- 0
    #self <- runif(length(diag(multityp.fill)), 0, 1)
    diag(multityp.fill) <- runif(length(diag(multityp.fill)), -2, 0)
    a.i <- runif(nrow(multityp.fill), .1, .5)
    
    par1 <- list(alpha = runif(nrow(multityp.fill), 0, .1), m = multityp.fill)
    dyn <-(ode(a.i, times = 1:tf, func = lvm, parms = par1, events = list(func = efun, time =  1:tf)))
    
    if(any(is.na(dyn))){cond <- FALSE;next}
    #if(nrow(dyn) == tf & nrow(dyn2) == tf){
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
  totalbd <- (sum(dyn[tf,-1]) - sum(dyn[1,-1]))/sum(dyn[1,-1])
  npos <- sum(ch > 0)
  mpos <- median(ch[ch > 0])
  mp <- max(ch[ch > 0])
  nneg <- sum(ch < 0)
  mneg <- median(ch[ch < 0])
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

key_effect <- function(init, plots = TRUE){ 
  nt <- nrow(init$dyn1)
  spp <- init$dyn1[nt,-1] > 10^-10
  mat <- init$m[spp,spp]
  alp <- init$alpha[spp]
  ia <- init$dyn1[nt,-1][init$dyn1[nt,-1] > 10^-10]
  pars <- list(alpha = alp, m = mat) 
  
  div <- c()
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
  }
  cova <- t(sapply(cvx, function(x) apply(x, 2, median, na.rm = T)))
  impact <- cbind(per, bdiff, cova)
  
  return(impact)
}


###########################################
###########################################
###########################################
###########################################

s1 <- Sys.time()
init <- isim(50, 1000, efun = ext1, plot = TRUE)
#matplot(init$dyn1[,-1], typ = 'l', main = sum(init$dyn1[1000,-1] > 10^-10))
s1.2 <- Sys.time()
mat <- init$m[init$dyn1[1000,-1] > 10^-10,init$dyn1[1000,-1] > 10^-10]
diag(mat) <- 0

g <- graph.adjacency(abs(sign(mat)))
lay <- layout_with_kk(g)

#plot(g, vertex.size = degree(g), layout = lay)

ra <- init$dyn1[1000,-1][init$dyn1[1000,-1]> 10^-10]/sum(init$dyn1[1000,-1][init$dyn1[1000,-1]> 10^-10])
plot(g, vertex.size = ra*100, layout = lay, main = sum(init$dyn1[1000,-1] > 10^-10))

s2.1 <- Sys.time()
ke1 <- key_effect(init)
s2 <- Sys.time()
s2 - s1

ke1
###########################################
###########################################
###########################################
###########################################

# info on species level interactions

# network level data

# split network into component interaction webs