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
    
    par1 <- list(alpha = runif(nrow(multityp.fill), 0, .1), m = multityp.fill, aii = self)
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


init <- isim(50, 1000, efun = ext1, plot = TRUE)
#matplot(init$dyn1[,-1], typ = 'l', main = sum(init$dyn1[1000,-1] > 10^-10))

mat <- init$m[init$dyn1[1000,-1] > 10^-10,init$dyn1[1000,-1] > 10^-10]
diag(mat) <- 0

g <- graph.adjacency(abs(sign(mat)))
lay <- layout_with_kk(g)

#plot(g, vertex.size = degree(g), layout = lay)

ra <- init$dyn1[1000,-1][init$dyn1[1000,-1]> 10^-10]/sum(init$dyn1[1000,-1][init$dyn1[1000,-1]> 10^-10])
plot(g, vertex.size = ra*100, layout = lay, main = sum(init$dyn1[1000,-1] > 10^-10))


###########################################
###########################################
###########################################
###########################################

spp <- init$dyn1[1000,-1] > 10^-10
mat <- init$m[spp,spp]
alp <- init$alpha[spp]
ia <- init$dyn1[1000,-1][init$dyn1[1000,-1] > 10^-10]
pars <- list(alpha = alp, m = mat)
ch <- c()
tteq <- c()
ei <- c()
ext <- c()
div <- c()
cvi <- matrix(nrow = nrow(mat), ncol = ncol(mat))
cvf <- matrix(nrow = nrow(mat), ncol = ncol(mat))
for(i in 1:nrow(mat)){
  ia <- init$dyn1[1000,-1][init$dyn1[1000,-1] > 10^-10]
  ia[i] <- 0
  
  dyn <-(ode(ia, times = 1:1000, func = lvm, parms = pars, events = list(func = ext1, time =  1:1000)))
  if(nrow(dyn) != 1000){next}
  
  e_tr <- apply(dyn[,-1], 1, function(x){
    p1 <- list(alpha = alp[x>10^-10], m = mat[x>10^-10,x>10^-10])
    return(max(Re(eigen(jacobian.full(x[x>10^-10], func = lvm, parms = p1))$values)))
  })
  
  matplot(dyn[,-1], typ = "l", main = i)
  
  ext[i] <- sum(dyn[1,-1] > 10^-10) - sum(dyn[1000,-1] > 10^-10)
  ch[i] <- (sum(dyn[1000,-1]) - sum(dyn[1,-1]))/sum(dyn[1,-1])
  tteq[i] <- sum(e_tr > 0)
  ei[i] <- e_tr[1000]
  div[i] <- vegan::diversity(dyn[1000,-1][dyn[1000,-1]>10^-10])
  cvi[i,] <- apply(dyn[1:100,-1], 2, function(x) sd(x)/mean(x))
  cvf[i,] <- apply(dyn[900:1000,-1], 2, function(x) sd(x)/mean(x))
}
ch*100
div
tteq
plot(ei~ch)
plot(tteq~ch)
