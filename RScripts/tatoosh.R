# save.image("~/Desktop/tatoosh.Rdata") 
# load("~/Desktop/tatoosh.Rdata")


library(igraph)
library(NetIndices)
library(deSolve)
library(ggplot2)
library(MuMIn)
library(rootSolve)
library(DAAG)
library(reshape2)


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


tatoosh <- as.matrix(read.csv("~/Desktop/GitHub/rKeystone/tatoosh.csv", header = F))
sum(tatoosh != 0)/(nrow(tatoosh)*(nrow(tatoosh) - 1))
itypes(tatoosh)
itypes.sp(tatoosh)

nINT <- sum(tatoosh != 0)
SteinInt <- read.csv("~/Desktop/GitHub/microbial-dyn/Data/ecomod-ints.csv", row.names = 1)
INTs <- c(SteinInt[upper.tri(SteinInt)],SteinInt[lower.tri(SteinInt)])
INTstr <- abs(rnorm(nINT, mean(INTs), sd(INTs)))
tat2 <- tatoosh
tat2[tat2 != 0] <- INTstr
tat2 <- tat2 * tatoosh
diag(tat2) <- -rbeta(nrow(tat2), 1.1, 5)*5

gr <- runif(nrow(tat2), .1, 1)
parms <- list(alpha = gr, m = tat2)

# numerical integration of ODE, simulates dynamics of local community
test <- ode(runif(nrow(tat2), .1, 10), 1:2000, parms = parms, func = lvmod, events = list(func = ext1, time =  1:2000))
matplot(test[,-1], type = "l")
which(test[2000,-1] > 0)



x = 0
eqst <- list()
parms1 <- list()
eqab <- list()
init.par <- list()
for(i in 1:500){
  INTstr <- abs(rnorm(nINT, mean(INTs), sd(INTs)))
  tat2 <- tatoosh
  tat2[tat2 != 0] <- INTstr
  tat2 <- tat2 * tatoosh
  diag(tat2) <- -rbeta(nrow(tat2), 1.1, 5)*5
  
  gr <- runif(nrow(tat2), .1, 1)
  parms <- list(alpha = gr, m = tat2)
  
  test <- ode(runif(nrow(tat2), .1, 10), 1:1000, parms = parms, func = lvmod, events = list(func = ext1, time =  1:1000))
  
  print(i)
  if(nrow(test) == 1000){
    x <- x+1
    init.par[[x]] <- parms
    eqab[[x]] <- test[1000,-1] 
    eqst[[x]] <- tat2[which(test[1000,-1] > 0),which(test[1000,-1] > 0)]
    parms1[[x]] <- list(alpha = gr[which(test[1000,-1] > 0)], m = eqst[[x]])
  }else{
    next
  }
}

removal <- function(iter, parms1, eqab, eqst){
  eqmat <- matrix(nrow = nrow(eqst[[iter]]), ncol = nrow(eqst[[iter]]))
  for(i in 1:nrow(eqst[[iter]])){
    parms <- parms1[[iter]]
    parms$alpha <- parms$alpha[-i]
    parms$m <- parms$m[-i,-i]
    test2 <- ode(eqab[[iter]][which(eqab[[iter]] > 0)][-i], 1:1000, parms = parms, func = lvmod, events = list(func = ext1, time =  1:1000))
    #matplot(test2[,-1], type = "l", main = i)
    #plot(test2[2000,-1],test[2000,-1][which(test[2000,-1] > 0)][-i], main = i)
    if(nrow(test2) == 1000){
      eqmat[i, -i] <- test2[1000,-1]-eqab[[iter]][which(eqab[[iter]] > 0)][-i]
    }else{
      eqmat[i, -i] <- NA
    }
    
  }
  return(eqmat)
}

removal.alt <- function(iter, parms1, eqab, eqst){
  eqmat <- matrix(nrow = nrow(eqst[[iter]]), ncol = nrow(eqst[[iter]]))
  for(i in 1:nrow(eqst[[iter]])){
    parms <- parms1[[iter]]
    eq1 <- eqab[[iter]][which(eqab[[iter]] > 0)]
    eq1[i] <- 10^-4
    
    test2 <- ode(eq1, 1:1000, parms = parms, func = lvmod, events = list(func = ext1, time =  1:1000))
    #matplot(test2[,-1], type = "l", main = i)
    #plot(test2[2000,-1],test[2000,-1][which(test[2000,-1] > 0)][-i], main = i)
    if(nrow(test2) == 1000){
      eqmat[i, ] <- test2[1000,-1]-eqab[[iter]][which(eqab[[iter]] > 0)]
    }else{
      eqmat[i, ] <- NA
    }
    
  }
  return(eqmat)
}

alleq <- lapply(1:length(eqst), function(x) removal(x, parms1, eqab, eqst))
alleq.alt <- lapply(1:length(eqst), function(x) removal.alt(x, parms1, eqab, eqst))

gs <- lapply(eqst, function(x) graph.adjacency(abs(x), weighted = T, diag = F))
gs2 <- lapply(eqst, function(x){x[x<0] <- 0; graph.adjacency(abs(x), weighted = T, diag = F)})
hist(unlist(sapply(gs[sapply(gs, is.connected)], betweenness)))
hist(unlist(sapply(gs, closeness)))

n.ext <- lapply(alleq[sapply(gs, is.connected)], function(x) apply(x, 1, function(y) cbind(length(y)-length(which(y == 0))-1, length(which(y == 0)))))
n.ext2 <- do.call(rbind, lapply(n.ext, t))
p.ext <- n.ext2[,2]/rowSums(n.ext2)

p.ext

summary(glm(n.ext2~unlist(sapply(gs[sapply(gs, is.connected)], betweenness))+unlist(sapply(gs[sapply(gs, is.connected)], closeness)), family = "binomial"))

#plot(apply(eqmat, 1, function(x) sum(x == 0, na.rm = T))~closeness(graph.adjacency(abs(eqst[[iter]]), weighted = T)), col = "blue", pch = 20)
#points(glm(apply(eqmat, 1, function(x) sum(x == 0, na.rm = T))~closeness(graph.adjacency(abs(eqst[[iter]]), weighted = T)), family = "poisson")$fitted.values~closeness(graph.adjacency(abs(eqst[[iter]]), weighted = T)))

#plot(rowMeans(eqmat, na.rm = T)~closeness(graph.adjacency(abs(eqst[[iter]]), weighted = T)), col = "blue", pch = 20)
#summary(glm(rowMeans(eqmat, na.rm = T)~closeness(graph.adjacency(abs(eqst[[iter]]), weighted = T)), family = "gaussian"))
#points(glm(rowMeans(eqmat, na.rm = T)~closeness(graph.adjacency(abs(eqst[[iter]]), weighted = T)), family = "gaussian")$fitted.values~closeness(graph.adjacency(abs(eqst[[iter]]), weighted = T)))


x = 0
eqst1 <- list()
parms2 <- list()
eqab1 <- list()
init.par1 <- list()
initAB <- runif(nrow(tat2), .1, 10)
for(i in 1:200){
  INTstr <- abs(rnorm(nINT, mean(INTs), sd(INTs)))
  tat2 <- tatoosh
  tat2[tat2 != 0] <- INTstr
  tat2 <- tat2 * tatoosh
  p1 <- runif(1, .05, .2)
  tat2 <- tat2 * sample(c(-1,1), length(tat2), prob = c(p1,1-p1), replace = T)
  diag(tat2) <- -rbeta(nrow(tat2), 1.1, 5)*5
  
  gr <- runif(nrow(tat2), .1, 1)
  parms <- list(alpha = gr, m = tat2)
  
  test <- ode(initAB, 1:1000, parms = parms, func = lvmod, events = list(func = ext1, time =  1:1000))
  
  print(i)
  if(nrow(test) == 1000){
    x <- x+1
    init.par1[[x]] <- parms
    eqab1[[x]] <- test[1000,-1] 
    eqst1[[x]] <- tat2[which(test[1000,-1] > 0),which(test[1000,-1] > 0)]
    parms2[[x]] <- list(alpha = gr[which(test[1000,-1] > 0)], m = eqst[[x]])
  }else{
    next
  }
}