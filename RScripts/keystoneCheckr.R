# load("~/Desktop/simul-example3.Rdata")

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
###
###

kbp <- list()
ext2 <- list()
for(i in 1:sum(use)){
  spp1 <- eqcomm[[i]]
  kbio <- ks2[[i]] 
  eqa <- eq.abund2[[i]]
  kbioper <- matrix(nrow = nrow(kbio), ncol = ncol(kbio))
  sec.ext <- matrix(nrow = length(spp1), ncol = length(spp1))
  for(j in 1:nrow(kbio)){
    kbioper[j,] <- kbio[j,]/eqa[-j]
    
    sec.ext[j,] <- spp1 %in% spp1[-i][which(kbioper[j,] == -1)] 
  }
  
  ext2[[i]] <- sec.ext 
  kbp[[i]] <- kbioper 
}

lapply(ext2, function(x) sapply(apply(x, 1, which), length))

eqcomm2 <- rep(sapply(eqcomm, length), sapply(eqcomm, length))[ccak]
hist(mydat$pers/eqcomm2)

dab <- cbind(unlist(lapply(kbp, function(x) apply(x, 1, function(q) mean(q[q < 0])))),unlist(lapply(ks2, function(x) apply(x, 1, function(q) mean(q[q > 0])))))
head(dab)


###
###
###
mats2 <- mats
mats2[mats2!=0] <- 0 
diag(mats2) <- diag(mats)
mats2

r3 <- list()
for(i in 1:sum(use)){
  #spp[[i]] <- c(sample(71:200, 20), sample(1:70, 30))               # sample species from the pool (two samples ensure some overlap)
  isp <- spp[use][[i]]                                                   # local species community
  parms <- list(alpha = growth[isp], m = mats2[isp,isp], k = K[isp]) # named parameter list (growth rates and int mat; K not used in sim)
  
  # numerical integration of ODE, simulates dynamics of local community
  r3[[i]] <- ode(runif(length(isp), .1, 1), 1:1000, parms = parms, func = lvmod, events = list(func = ext1, time =  1:1000))
  matplot(r3[[i]][,-1], typ = "l")                                  # plot community dynamics
  print(i)                                                          # which iteration are we on again? 
}

eq.ni <- t(sapply(r3, function(x) x[1000,-1]))
eq.wi <- t(sapply(r2[use], function(x) x[1000,-1]))

eqdiff <- matrix(nrow = sum(use), ncol = 200)
for(i in 1:sum(use)){
  eqdiff[i,spp[use][[i]]] <- eq.ni[i,] - eq.wi[i,]
}

eqdiff


###
###
###

# initial interaction matrices
inmatuse <- lapply(1:sum(use), function(x) mats[spp[use][[x]],spp[use][[x]]])
# equilibrium interaction matrices for all iterations that didn't mess up
matuse <- lapply(1:sum(use), function(i) mats[eqcomm[[i]], eqcomm[[i]]])

in.itysp <- lapply(inmatuse, itypes.sp)
in.istrsp <- lapply(inmatuse, istr.sp)

gr1 <- unlist(lapply(1:sum(use), function(x) growth[spp[[x]]]))
dia1 <- unlist(lapply(inmatuse, diag))

in.itysp[[1]]
eq1 <- lapply(r2[use], function(x) x[1000,-1] > 0)
eq2 <- lapply(r2[use], function(x) x[1000,-1])

eqi <- cbind(do.call(rbind, in.itysp)[,1:3], unlist(eq1), unlist(eq2))
head(eqi)
colnames(eqi) <- c("comp", "mut", "pred", "eqcom", "eqabund")
summary(glm(eqi[,4]~eqi[,1:3], family = "binomial"))

eqi2 <- cbind(do.call(rbind, in.itysp)[,1:3], do.call(rbind, in.istrsp), unlist(eq1), unlist(eq2))
colnames(eqi2) <- c("comp", "mut", "pred", "compS", "mutS", "predS", "eqcom", "eqabund")
eqi2 <- as.data.frame(eqi2)
fitA <- glm(eqcom~comp+abs(compS)+mut+mutS+pred+predS, data = eqi2, family = "binomial", na.action = "na.fail")
fitB <- glm(eqabund~comp+abs(compS)+mut+mutS+pred+predS+gr1[which(eqi2$eqcom == 1)]+abs(dia1[which(eqi2$eqcom == 1)]), data = eqi2[which(eqi2$eqcom == 1),], family = "gaussian", na.action = "na.fail")
ciA <- confint(fitA)
ciB <- confint(fitB)
modat <- data.frame(ciA, coeff <- fitA$coefficients, rownames(ciA))
colnames(modat) <- c("lower", "upper", "coeff", "met")
summary(fitA)
confint(model.avg(dredge(fitB), subset = delta < 2))

ggplot(modat) + geom_segment(aes(x = lower, y = met, xend = upper, yend = met)) + geom_vline(aes(xintercept = 0)) + 
  geom_point(aes(x = coeff, y = met)) + xlab("Value") + ylab("Variable") + theme_bw()
