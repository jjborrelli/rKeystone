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
gr2 <- sapply(lapply(1:sum(use), function(x) growth[spp[[x]]]), mean)
dia2 <- sapply(lapply(inmatuse, diag), mean)


in.itysp[[1]]
eq1 <- lapply(r2[use], function(x) x[1000,-1] > 0)
eq2 <- lapply(r2[use], function(x) x[1000,-1])

eqi <- cbind(do.call(rbind, in.itysp), unlist(eq1), unlist(eq2))
head(eqi)
colnames(eqi) <- c("comp", "mut", "pred", "amens", "comm", "eqcom", "eqabund")
summary(glm(eqi[,6]~eqi[,1:5], family = "binomial"))

eqi2 <- cbind(do.call(rbind, in.itysp), do.call(rbind, in.istrsp), unlist(eq1), unlist(eq2))
colnames(eqi2) <- c("comp", "mut", "pred", "amens", "comm", "compS", "mutS", "predS", "eqcom", "eqabund")
eqi2 <- as.data.frame(eqi2)
fitA <- glm(eqcom~comp+abs(compS)+mut+mutS+pred+predS+amens+comm, data = eqi2, family = "binomial", na.action = "na.fail")
fitB <- glm(eqabund~comp+abs(compS)+mut+mutS+pred+predS+amens+comm+gr1[which(eqi2$eqcom == 1)]+abs(dia1[which(eqi2$eqcom == 1)]), data = eqi2[which(eqi2$eqcom == 1),], family = "gaussian", na.action = "na.fail")

plotCI(fitA)
plotCI(fitB)

ityA <- apply(sapply(inmatuse, itypes), 2, function(x) x/sum(x))
in.conn <- sapply(inmatuse, function(x) sum(x != 0)/(50*49))
sumneg <- sapply(inmatuse, function(X) sum(X < 0)/sum(X != 0))
imean <- sapply(inmatuse, function(x) mean(c(x[upper.tri(x)][x[upper.tri(x)]!=0],x[lower.tri(x)][x[lower.tri(x)]!=0])))

test <- cbind(t(ityA),t(sapply(in.istrsp, colMeans)))
summary(glm(sapply(eq1, function(x) sum(x)/length(x))~test+imean+abs(dia2), family = "quasibinomial"))
summary(glm(sapply(eq2, function(x) mean(x[x!=0]))~dia2, family = "gaussian"))
summary(glm(sapply(eq1, function(x) sum(x)/length(x))~imean+abs(dia2), family = "quasibinomial"))

spi <- lapply(inmatuse, spints)
spi2 <- do.call(rbind, spi)
colnames(spi2) <- c("pos1", "neg1", "non1", "spos1", "sneg1", "pos2", "neg2", "non2", "spos2", "sneg2")

spi3 <- data.frame(eqa = eqi2$eqabund, eco = eqi2$eqcom, spi2)
fitC <- glm(eco~pos1+neg1+spos1+sneg1+pos2+neg2+spos1+sneg2, data = abs(spi3), family = "binomial")
fitD <- glm(eqa~pos1+neg1+spos1+sneg1+pos2+neg2+spos1+sneg2, data = abs(spi3), family = "gaussian")
summary(fitC)
summary(fitD)

plotCI(fitC)
plotCI(fitD)
