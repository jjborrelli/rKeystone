# save.image("~/Desktop/tatoosh.Rdata") 
# load("~/Desktop/tatoosh.Rdata")


library(igraph)
library(NetIndices)
library(deSolve)
library(ggplot2)
library(rootSolve)
library(reshape2)
library(boot)
library(data.table)


# Basic Lotka-Volterra model
lvmod <- function(times, state, parms){
  with(as.list(c(state, parms)), {
    dB <- state * parms$alpha + state * parms$m %*% state 
    
    list(dB)
  })
}

# Basic Lotka-Volterra model
lvmodK <- function(times, state, parms){
  with(as.list(c(state, parms)), {
    dB <- state * parms$alpha * (1 - state/(parms$K/sum(state > 0))) + state * parms$m %*% state 
    
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


tatoosh <- as.matrix(read.csv("~/Desktop/GitHub/rKeystone/tatoosh.csv", header = F))
sum(tatoosh != 0)/(nrow(tatoosh)*(nrow(tatoosh) - 1))
itypes(tatoosh)
itypes.sp(tatoosh)


p1 <- runif(1,0,1)
tat <- tatoosh*sample(c(1,-1), length(tatoosh), replace = T, prob = c(p1,1-p1))
itypes(tat)


#SteinInt <- read.csv("C:/Users/jjborrelli/Desktop/GitHub/microbial-dyn/Data/ecomod-ints.csv", row.names = 1)

SteinInt <- read.csv("~/Desktop/GitHub/microbial-dyn/Data/ecomod-ints.csv", row.names = 1)
INTs <- c(SteinInt[upper.tri(SteinInt)],SteinInt[lower.tri(SteinInt)])


tats <- lapply(1:10, function(x){
  p1 <- runif(1,0,1)
  tat <- tatoosh*sample(c(1,-1), length(tatoosh), replace = T, prob = c(p1,1-p1))
  return((tat))
})

ity <- sapply(tats, itypes)

par(mfrow = c(5,1))
apply(ity, 1, hist)
dev.off()

dyn <- list()
tats2 <- list()
grs <- list()
for(i in 1:length(tats)){
  t1 <- tats[[i]]
  diag(t1) <- 0#-rbeta(length(diag(t1)), 1.1, 5)*5
  t1[t1 == 1] <- abs(rnorm(sum(t1 == 1), mean(INTs), sd(INTs)))
  t1[t1 == -1] <- -abs(rnorm(sum(t1 == -1), mean(INTs), sd(INTs)))
  
  gr <- runif(nrow(t1), .1, 1)
  parms <- list(alpha = gr, m = t1, K = 20)
  
  # numerical integration of ODE, simulates dynamics of local community
  test <- ode(runif(nrow(t1), .1, 10), 1:2000, parms = parms, func = lvmodK, events = list(func = ext1, time =  1:2000))
  
  if(nrow(test) == 2000){
    dyn[[i]] <- test[,-1]
    tats2[[i]] <- t1
    grs[[i]] <- gr
  }else{
    dyn[[i]] <- NA
    tats2[[i]] <- NA
    grs[[i]] <- NA
  }
  
  matplot(test[,-1], type = "l", main = i)
  print(length(dyn[!is.na(dyn)]))
}

sum(is.na(dyn))

#initial interaction matrices for communities that worked
tats2 <- tats2[!is.na(dyn)]#[sapply(dyn[!is.na(dyn)], function(x) sum(is.na(x)) == 0) & !sapply(dyn, function(x) sum(x < 0) > 0)]
#growth rates of spp for communities that worked
grs <- grs[!is.na(dyn)]#[sapply(dyn[!is.na(dyn)], function(x) sum(is.na(x)) == 0)]
#dynamics of communities that worked
dyn <- dyn[!is.na(dyn)]#[sapply(dyn[!is.na(dyn)], function(x) sum(is.na(x)) == 0)]

#final networks of communities that worked
ags <- lapply(1:length(dyn), function(x) graph.adjacency(abs(sign(tats2[[x]][dyn[[x]][2000,] > 0,dyn[[x]][2000,] > 0])), diag = F))
#final networks of communities that worked, without isolated nodes
ags.rv <- lapply(ags, function(x) delete_vertices(x, which(degree(x) == 0)))
#get ids of connected nodes in ags.rv
conn.com <- lapply(ags.rv, function(x) colnames(tatoosh) %in% names(V(x)))
#get dynamics of connected nodes
condyn <- lapply(1:length(dyn), function(x) dyn[[x]][,conn.com[[x]]])
#get interaction matrices of connected nodes
conmat <- lapply(1:length(tats2), function(x) tats2[[x]][conn.com[[x]], conn.com[[x]]])

p1 <- list(alpha = grs[[13]][conn.com[[13]]], m = conmat[[13]], K = 30)
test1 <- ode(condyn[[13]][2000,], 1:100, parms = p1, func = lvmodK, events = list(func = ext1, time =  1:100))
matplot(test1[,-1], typ = "l")


ags.i <- lapply(tats2, function(x) graph.adjacency(abs(sign(x)), diag = F))
bet.i <- lapply(1:length(ags.i), function(x) betweenness(ags.i[[x]]))
bet.f <- lapply(ags.rv, betweenness)

ityspl <- rbindlist(lapply(tats2, function(x) as.data.frame(itypes.sp(x))))
ityspl2 <- as.data.frame(t(apply(ityspl, 1, function(x) x/sum(x))))

par(mfrow = c(1,2))
plot(betsur$S~ityspl2$V5)
fitA <- glm(betsur$S~ityspl2$V5, family = "binomial")
cbA <- predict(fitA, se.fit = T)
points(fitA$fitted.values~ityspl2$V5, pch = 20)
points(inv.logit(cbA$fit+1.96*cbA$se.fit)~ityspl2$V5, pch = 20, col = "blue")
points(inv.logit(cbA$fit-1.96*cbA$se.fit)~ityspl2$V5, pch = 20, col = "blue")


plot(betsur$S~ityspl2$V4)
fitB <- glm(betsur$S~ityspl2$V4, family = "binomial")
cbB <- predict(fitB, se.fit = T)
points(fitB$fitted.values~ityspl2$V4, pch = 20)
points(inv.logit(cbB$fit+1.96*cbB$se.fit)~ityspl2$V4, pch = 20, col = "blue")
points(inv.logit(cbB$fit-1.96*cbB$se.fit)~ityspl2$V4, pch = 20, col = "blue")

betsur <- rbindlist(lapply(1:length(ags.i), function(x) data.frame(B = bet.i[[x]], S = (dyn[[x]][2000,] != 0), A = dyn[[x]][2000,])))
betsur1 <- rbindlist(lapply(1:length(ags.i), function(x) data.frame(B = bet.i[[x]][conn.com[[x]]], B2 = bet.f[[x]], A = condyn[[x]][2000,])))
plot(betsur1$B, betsur1$B2)
abline(a = 0, b = 1, xpd = F)
plot(betsur1$B2, betsur1$A)


################################################################################################################
################################################################################################################
# GLM Predicting extinction risk as a funciton of variability

coefvar <- function(dyna){
  cv10 <- apply(dyna[1:5,], 2, function(x) sd(x)/mean(x))
  cv100 <- apply(dyna[1:100,], 2, function(x) sd(x)/mean(x))
  cv1k <- apply(dyna[1:1000,], 2, function(x) sd(x)/mean(x))
  cvfi <- apply(tail(dyna, 100), 2, function(x) sd(x)/mean(x))
  
  return(data.frame(cv10 = cv10, cv100 = cv100, cv1k = cv1k, cvfi100 = cvfi))
}

cvi <- lapply(dyn, coefvar)
cdat <- data.frame(S = betsur$S, cv10 = unlist(rbindlist(cvi)$cv10), cv100 = unlist(rbindlist(cvi)$cv100), iB = betsur$B)
fitcv <- glm(S~cv10, data = cdat, family = "binomial")
cb <- predict.glm(fitcv, se.fit = T)

plot(betsur$S~unlist(rbindlist(cvi)$cv10))
plot(inv.logit(cb$fit)~cdat$cv10, col = "blue")
lines(inv.logit(cb$fit+1.96*cb$se.fit)~seq(0, 5, .01), col = "darkgreen")
lines(inv.logit(cb$fit-1.96*cb$se.fit)~seq(0, 5, .01), col = "darkgreen")


plot(inv.logit(cb$fit)~cdat$cv10, col = "blue", pch = 20)
points(inv.logit(cb$fit+1.96*cb$se.fit)~cdat$cv10, col = "darkgreen", pch = 18)
points(inv.logit(cb$fit-1.96*cb$se.fit)~cdat$cv10, col = "darkgreen", pch = 18)

################################################################################################################
################################################################################################################
# not add log(deltaBIOplus/deltaBIOminus)

remove.sp <- function(sp, parms, states){
  states[sp] <- 0
  test <- ode(states, 1:2000, parms = parms, func = lvmodK, events = list(func = ext1, time =  1:2000))

  if(nrow(test) == 2000){
    bi <- betweenness(graph.adjacency(abs(sign(parms$m))))
    br <- rep(NA, length(bi))
    br[-sp] <- betweenness(graph.adjacency(abs(sign(parms$m[-sp,-sp]))))
    bf <- rep(NA, length(bi))
    bf[test[2000,-1]!=0] <- betweenness(graph.adjacency(abs(sign(parms$m[test[2000,-1]!=0,test[2000,-1]!=0]))))
    
    di <- degree(graph.adjacency(abs(sign(parms$m))))
    dr <- rep(NA, length(di))
    dr[-sp] <- degree(graph.adjacency(abs(sign(parms$m[-sp,-sp]))))
    df <- rep(NA, length(di))
    df[test[2000,-1]!=0] <- degree(graph.adjacency(abs(sign(parms$m[test[2000,-1]!=0,test[2000,-1]!=0]))))
    
    
    sppdf <- data.frame(spR = sp, spT = 1:nrow(parms$m), coefvar(test[,-1]), iBio = test[1,-1],
                        t10Bio = test[10,-1], t100Bio = test[100,-1], fBio = test[2000,-1],
                        iBet = bi, rBet = br, fBet = bf, iDeg = di, rDeg = dr, fDeg = df,
                        gr = parms$alpha, self = diag(parms$m),
                        exts = (test[2000,-1] != 0))
  }else{
    sppdf <- NA
  }
  
  return(sppdf)
}


spr <- list()
errs <- list()
ity <- list()
n.exts <- list()
for(comm in 1:length(dyn)){
  par <- list(alpha = grs[[comm]][conn.com[[comm]]], m = conmat[[comm]], K =20)
  st1 <- condyn[[comm]][2000,]
  
  dflist <- list()
  for(i in 1:nrow(par$m)){
    dflist[[i]] <- remove.sp(i, parms = par, states = st1)
    print(i)
  }
  
  errs[[comm]] <- is.na(dflist)
  ity[[comm]] <- itypes.sp(par$m)[!errs[[comm]],]
  dflist <- dflist[!is.na(dflist)]
  spr[[comm]] <- data.frame(comm = comm, do.call(rbind, dflist))
  n.exts[[comm]] <- sapply(dflist, function(x) sum(!x$exts)-1)
  
  cat("---------------------|| ", comm, " ||---------------------", "\n")
}

sprl <- rbindlist(spr)

ag1 <- aggregate(sprl$fBio-sprl$iBio, list(sprl$comm, sprl$spR), vegan::diversity)
aggregate(ag1$x, list(ag1$Group.1), median)

t(sapply(conmat, itypes))[,5]
i=1
mdiv <- c()
lRR <- list()
for(i in 1:length(spr)){
  spl1 <- split(spr[[i]], f = spr[[i]]$spR)
  mdiv[i] <- median(sapply(spl1, function(x) (vegan::diversity(x$fBio) - vegan::diversity(x$iBio))/vegan::diversity(x$iBio))*100)
  lRR[[i]] <- sapply(1:length(spl1), function(x){
    log(abs(mean((spl1[[x]]$fBio - spl1[[x]]$iBio)[-x][(spl1[[x]]$fBio - spl1[[x]]$iBio)[-x] < 0])/mean((spl1[[x]]$fBio - spl1[[x]]$iBio)[-x][(spl1[[x]]$fBio - spl1[[x]]$iBio)[-x] > 0])))
  })
}
mdiv

lRR <- sapply(1:length(spl1), function(x){
  log(abs(mean((spl1[[x]]$fBio - spl1[[x]]$iBio)[-x][(spl1[[x]]$fBio - spl1[[x]]$iBio)[-x] < 0])/mean((spl1[[x]]$fBio - spl1[[x]]$iBio)[-x][(spl1[[x]]$fBio - spl1[[x]]$iBio)[-x] > 0])))
})
spl1

posint <- rbindlist(lapply(ity, function(x) as.data.frame(x[,c(1,2,3,4,5)])))
ex1 <- unlist(n.exts)

cval <- c()
for(i in 1:22){
  co <- i
  ctest <- cor.test(abs(tapply(spr[[co]]$fBio - spr[[co]]$iBio, list(spr[[co]]$spR), median)), tapply(spr[[co]]$cv10, list(spr[[co]]$spR), median, na.rm = T))
  cval[i] <- ctest$estimate 
}
co <- 14
cor.test(unlist(lapply(spr, function(x) tapply(x$cv10, list(x$spR), median, na.rm = T))), ex1)

comm <- 3
par <- list(alpha = grs[[comm]][conn.com[[comm]]], m = conmat[[comm]])
st1 <- condyn[[comm]][2000,]

dflist <- list()
for(i in 1:nrow(par$m)){
  dflist[[i]] <- remove.sp(i, parms = par, states = st1)
  print(i)
}
errs <- is.na(dflist)
dflist <- dflist[!is.na(dflist)]
spr <- rbindlist(dflist)

sum(!rbindlist(spr)$exts) - sum(sapply(spr, function(x) max(x$spR)))


plot(spr$cv10,spr$exts)
fit1 <- glm(exts~cv10, data = spr, family = "binomial")
cb1 <- predict.glm(fit1, se.fit = T)
points(fit1$fitted.values~spr$cv10[!is.na(spr$cv10)], col = "darkgreen", pch = 20)
points(inv.logit(cb1$fit+1.96*cb1$se.fit)~spr$cv10[!is.na(spr$cv10)], col = "blue", pch = 18)
points(inv.logit(cb1$fit-1.96*cb1$se.fit)~spr$cv10[!is.na(spr$cv10)], col = "blue", pch = 18)
#

plot(spr$cv10, spr$t10Bio)

summary(lm(sapply(dflist, function(x) sum(!x$exts)-1)~(itypes.sp(par$m)[,5]+itypes.sp(par$m)[,2])))
plot(spr$cv10, spr$cv100)
abline(a = 0, b = 1, xpd = F)

################################################################################################################
################################################################################################################

cm1 <- conmat[[1]]
g1 <- grs[[1]][conn.com[[1]]]

eig.i <- max(Re(eigen(jacobian.full(condyn[[1]][2000,], func = lvmod, parms = list(alpha = g1, m = cm1)))$values))


qss <- function(cd1, cm1, g1){
  eig.c <- c()
  cm2 <- cm1
  adj.ints <- c(cm1[upper.tri(cm1)][cm1[upper.tri(cm1)]!=0],cm1[lower.tri(cm1)][cm1[lower.tri(cm1)]!=0])
  for(i in 1:1000){
    cm2[cm1 != 0] <- rnorm(sum(cm1 != 0), mean = cm1[cm1 != 0], sd = abs(cm1[cm1 != 0]*.1))
    diag(cm2) <- diag(cm1)
    eig.c[i] <- max(Re(eigen(jacobian.full(cd1, func = lvmod, parms = list(alpha = g1, m = cm2)))$values))
  }
  
  return(eig.c)
}

q2 <- c()
for(x in 1:length(dyn)){
  par1 <- condyn[[x]][2000,]
  par2 <- conmat[[x]]
  par3 <- grs[[x]][conn.com[[x]]]
  q1 <- qss(par1, par2, par3)
  q2[x] <- sum(q1 < 0)/1000
}
q2
inty <- t(apply(sapply(conmat, itypes), 2, function(x) x/sum(x)))
summary(lm(q2~inty[,3]))
################################################################################################################
################################################################################################################
################################################################################################################
################################################################################################################

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
