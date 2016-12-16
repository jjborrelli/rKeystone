
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
    #states <- states*parms$r[ceiling(times)]
    states[states < 10^-5] <- 0 
    #if(sum(states >= 100) >= 1){states<-rep(0, length(states))} 
    return(c(states))
  })
}


ints1 <- read.csv("~/Desktop/GitHub/microbial-dyn/Data/ecomod-ints.csv", row.names = 1)
grow1 <- read.csv("~/Desktop/GitHub/microbial-dyn/Data/ecomod-Growth.csv")

parms <- list(alpha = unlist(grow1), m = as.matrix(ints1))

g1 <- graph.adjacency((parms$m), weighted = T, diag = F)
w1 <- abs(E(g1)$weight)
s1 <- (sign(E(g1)$weight))
s2 <- factor(rep("red", length(s1)), levels = c("red", "blue"))
s2[s1 == 1] <- factor("blue", levels = c("red", "blue"))

plot(g1, edge.width = w1*5, edge.color = as.character(s2), edge.arrow.size = 1.5, vertex.label.color = "black", vertex.color = "gray")

## #################################
## Stochasticity ###################
## parms$r <- rnorm(1000, 1, .1) ###
## #################################



res1 <- ode(runif(nrow(ints1),1,10), 1:1000, parms = parms, func = lvmod, events = list(func = ext1, time =  1:1000))
res1[1000,-1]
matplot(res1[,-1], typ = "l", lwd = 2)

g2 <- graph.adjacency((parms$m[which(res1[1000,-1] != 0),which(res1[1000,-1] != 0)]), weighted = T, diag = F)
w2 <- abs(E(g2)$weight)
sA <- (sign(E(g2)$weight))
sB <- factor(rep("red", length(s1)), levels = c("red", "blue"))
sB[sA == 1] <- factor("blue", levels = c("red", "blue"))

plot(g2, edge.width = w2*5, edge.color = as.character(sB), edge.arrow.size = 1.5, vertex.label.color = "black", vertex.color = "gray")


m2 <- parms$m[which(res1[1000,-1] != 0),which(res1[1000,-1] != 0)]
#m2[m2 < .2 & m2 > 0] <- 0
#m2[m2 > -.2 & m2 < 0] <- 0 
m2[m2 > 0] <- 0

g3 <- graph.adjacency(m2, weighted = T, diag = F)
w3 <- abs(E(g3)$weight)
sA1 <- (sign(E(g3)$weight))
sB1 <- factor(rep("red", length(sA1)), levels = c("red", "blue"))
sB1[sA1 == 1] <- factor("blue", levels = c("red", "blue"))

plot(g3, edge.width = w3*8, edge.color = as.character(sB1), edge.arrow.size = 1.5, vertex.label.color = "black", vertex.color = "gray")
betweenness(g3, weights = w3)

##################################################################
##################################################################
###
###     Species removals
###
##################################################################
##################################################################

iter <- 100
iAb <- sapply(1:iter, function(x) runif(11, 1, 10))

removal <- function(x, iter, iA, p = parms, lv = lvmod, ext = ext1){
  iA[x,] <- 0
  out <- list()
  for(i in 1:iter){
    out[[i]] <- ode(iA[,i], 1:1000, parms = p, func = lv, events = list(func = ext, time =  1:1000))
    print(i)
  }
  return(out)
}

iter <- 100
iAb <- sapply(1:iter, function(x) runif(11, 1, 10))
out <- removal(2, iter = 100, iA = iAb)
sapply(out, function(x) x[1000,-1])
