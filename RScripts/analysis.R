### SAVED WORK
# last saved 10-25-16
# alt save 10-26-16
# ms save 11-8-16 == example3
# save.image("~/Desktop/simex2.Rdata") 
# load("~/Desktop/simexT-5.Rdata")


###
### run sim
###
simstart <- Sys.time()
source("~/Desktop/GitHub/rKeystone/RScripts/communitySIM.R")
simend1 <- Sys.time()

simend1 - simstart
source("~/Desktop/GitHub/rKeystone/RScripts/keystoneSIM.R")
simend2 <- Sys.time()
simend2 - simend1
save.image("~/Desktop/simexT-5.Rdata")
###
### ANALYSIS
###
inmats1 <- lapply(1:1000, function(x) mats[spp[[x]], spp[[x]]])
umat <- sapply(inmats1[use], itypes)
unmat <- sapply(inmats1[!use], itypes)


# how equilibrial are the communities?
dyn <- lapply(r2[use], function(x){x[x < 0] <- 0; x})               # get new object of dynamics list for runs that worked
cv.eq <- sapply(dyn, function(x) apply(x[990:1000,-1][,x[1000,-1] > 0], 2, sd)/(colMeans(x[990:1000, -1][,x[1000,-1] > 0])))
#cv.eq[is.nan(cv.eq)] <- 0
#range(colMeans(cv.eq))
hist(unlist(cv.eq))

eq.abund2 <- (lapply(dyn, function(x) x[1000,-1][x[1000,-1] > 0]))
# matrix of species found in each local equilibrium community
# can be used to determine compositional similarity of communities
eqmat <- matrix(0, nrow = sum(use), ncol = S)                       # initialize eqmat
eqmat1 <- matrix(0, nrow = sum(use), ncol = S)
for(i in 1:sum(use)){
  eqmat[i,eqcomm[[i]]] <- 1 
  eqmat1[i,eqcomm[[i]]] <- eq.abund2[[i]] # if the species is present in local comm i it gets a one, 0 otherwise
}

# initial interaction matrices
inmatuse <- lapply(1:sum(use), function(x) mats[spp[use][[x]],spp[use][[x]]])
# equilibrium interaction matrices for all iterations that didn't mess up
matuse <- lapply(1:sum(use), function(i) mats[eqcomm[[i]], eqcomm[[i]]])

###
###
###

# how each spp removal affects eigenval of comm
eigkey <- lapply(1:sum(use), function(x) eigenkey(mat = mats, growth = growth, isp = eqcomm[[x]], dyna = dyn[[x]]))        

allg <- lapply(matuse, getgraph)                                    # get the network for each local eq comm
betw <- lapply(allg, betweenness)                                   # get betweenness of each node
clocent <- lapply(allg, closeness)                                  # get closeness centrality
# get neighborhood of each spp going out 2 links
g.neighbors2 <- lapply(1:length(allg), function(x){sapply(graph.neighborhood(allg[[x]], 2), function(y) length(V(y)))})
ecent <- lapply(allg, function(x) eigen_centrality(x)$vector)       # get eigenvector centrality
hscore <- lapply(allg, function(x) hub_score(x)$vector)             # get hub score
p.rank <- lapply(allg, function(x) page_rank(x)$vector)             # get page rank algo



###
###
###
# get frequency of interaction types for each species in each equilibrial comm
itySP <- lapply(matuse, itypes.sp)
# get mean strength of each interaction type species participate in
istrSP <- lapply(matuse, istr.sp)

itySP2 <- do.call(rbind, itySP)                                     # get species level participation in interaction types
istrSP2 <- do.call(rbind, istrSP)                                   # get species level interaction type strengths

itySP3 <- t(apply(itySP2, 1, function(x) x/sum(x)))
itySP3[is.nan(itySP3)] <- 0

spi <- lapply(matuse, spints)
spin2 <- do.call(rbind, spi)
colnames(spin2) <- c("pos1", "neg1", "non1", "spos1", "sneg1", "pos2", "neg2", "non2", "spos2", "sneg2")
#es <- spi2[,5] + spi2[,4]
#eo <- spi2[,10] + spi2[,9]

ityinsp <- lapply(inmatuse, itypes.sp)
ityinsp2 <- do.call(rbind, ityinsp)


# put all data together in single matrix
allks <- do.call(rbind, lapply(ks1, function(x) x[,1:4]))           # all biomass, variation, and persistence
allks <- cbind(allks, eig = unlist(eigkey), sp.id = unlist(eqcomm), n.comp = itySP2[,1], n.mut = itySP2[,2], n.pred = itySP2[,3], 
               n.amen = itySP2[,4], n.com = itySP2[,5], 
               s.comp = istrSP2[,1], s.mut = istrSP2[,2], s.pred = istrSP2[,3], bet = unlist(betw), close = unlist(clocent),
               neigh = unlist(g.neighbors2),  ec = unlist(ecent), hub = unlist(hscore), pr = unlist(p.rank))
ccak <- complete.cases(allks)                                       # only use complete cases

spin3 <- spin2[ccak,]
eq.abund <- unlist(lapply(dyn, function(x) x[1000,-1][x[1000,-1] > 0]))[ccak]



########################
# 
ev.init <- c()
for(i in 1:sum(use)){
  dyna <- r2[use][[i]]
  isp <- eqcomm[[i]]
  eq.biom <- dyna[1000,-1][dyna[1000,-1] > 0]
  j1 <- jacobian.full(eq.biom, lvmod, parms = list(alpha = growth[isp], m = mats[isp,isp]))
  
  ev.init[i] <- max(Re(eigen(j1)$values))
}

evinit <- rep(ev.init, sapply(eqcomm, length))[ccak]
mydat <- as.data.frame(allks[ccak,])
destab <- (mydat$eig >= evinit)

eab <- unlist(lapply(eq.abund2, function(x) x/sum(x)))[ccak]
icv <- rep(sapply(cv.eq, mean), sapply(eqcomm, length))[ccak]
neq <- rep(sapply(eqcomm, length), sapply(eqcomm, length))[ccak]
evi <- rep(ev.init, sapply(eqcomm, length))[ccak]

dAbund1 <- lapply(lapply(ks2, function(x){apply(x, 1, function(y){c(mNeg = mean(y[y < 0]), nNeg = sum(y < 0)/length(y), mPos = mean(y[y > 0]), nPos = sum(y > 0)/length(y))})}),t)
dAbund2 <- do.call(rbind, dAbund1)[ccak,]

## Modeling
#### with MuMIn package

fit1A <- glm(delta.biom~n.comp*s.comp+n.mut*s.mut+n.pred*s.pred+bet+close+neigh+ec+hub+pr, family = "gaussian", data = mydat, na.action = "na.fail")
fit2A <- glm(mean.vary~n.comp*s.comp+n.mut*s.mut+n.pred*s.pred+bet+close+neigh+ec+hub+pr, family = "gaussian", data = mydat, na.action = "na.fail")
fit3A <- glm(m.init.vary~n.comp*s.comp+n.mut*s.mut+n.pred*s.pred+bet+close+neigh+ec+hub+pr, 
            family = "gaussian", data = mydat, na.action = "na.fail")
fit4A <- glm(cbind(pers,rep(sapply(eqcomm, length), sapply(eqcomm, length))[ccak])~n.comp*s.comp+n.mut*s.mut+n.pred*s.pred+bet+close+neigh+ec+hub+pr,
            family = "binomial", data = mydat, na.action = "na.fail")
fit5A <- glm(eig~n.comp*s.comp+n.mut*s.mut+n.pred*s.pred+bet+close+neigh+ec+hub+pr,
            family = "gaussian", data = mydat, na.action = "na.fail")


flist.A <- list(fit1A, fit2A, fit3A, fit4A, fit5A)
fAnames <- c("abund", "cvend", "cvinit", "pers", "eig")

dfall.A <- multimod(flist.A, fAnames)

ggplot(dfall.A) + geom_segment(aes(x = lower, y = met, xend = upper, yend = met, col = sig)) + geom_vline(aes(xintercept = 0)) + 
  geom_point(aes(x = coef, y = met, col = sig)) + facet_wrap(~mod, scales = "free_x") + 
  scale_color_manual(name = "Significant", values = c("grey", "blue")) + xlab("Value") + ylab("Variable") + theme_bw()

#####################################
#####################################
fit1B <- glm(delta.biom~n.comp+n.mut+n.pred+n.amen+n.com+s.comp+s.mut+s.pred, family = "gaussian", data = mydat, na.action = "na.fail")
fit2B <- glm(mean.vary~n.comp+n.mut+n.pred+n.amen+n.com+s.comp+s.mut+s.pred, family = "gaussian", data = mydat, na.action = "na.fail")
fit3B <- glm(m.init.vary~n.comp+n.mut+n.pred+n.amen+n.com+s.comp+s.mut+s.pred, 
            family = "gaussian", data = mydat, na.action = "na.fail")
fit4B <- glm(cbind(pers,rep(sapply(eqcomm, length), sapply(eqcomm, length))[ccak])~n.comp+n.mut+n.pred+n.amen+n.com+s.comp+s.mut+s.pred,
            family = "binomial", data = mydat, na.action = "na.fail")
fit5B <- glm(eig~n.comp+n.mut+n.pred+n.amen+n.com+s.comp+s.mut+s.pred, family = "gaussian", data = mydat, na.action = "na.fail")


flist.B <- list(fit1B, fit2B, fit3B, fit4B, fit5B)
fBnames <- c("abund", "cvend", "cvinit", "pers", "eig")

dfall.B <- multimod(flist.B, fBnames)

ggplot(dfall.B) + geom_segment(aes(x = lower, y = met, xend = upper, yend = met, col = sig)) + geom_vline(aes(xintercept = 0)) + 
  geom_point(aes(x = coef, y = met, col = sig)) + facet_wrap(~mod, scales = "free_x") + 
  scale_color_manual(name = "Significant", values = c("grey", "blue")) + xlab("Value") + ylab("Variable") + theme_bw()



########################
# 
CI.abund <- (mydat$delta.biom/rep(sapply(eq.abund2, mean), sapply(eqcomm, length))[ccak]) * (1/eab)
CI.ivary <- ((icv - mydat$m.init.vary)/icv) * (1/eab) 
CI.pers <- ((neq - mydat$pers)/neq) * (1/eab)
CI.eig <- ((evi - mydat$eig)/evi) * (1/eab)


fit1 <- glm(CI.abund~n.comp+abs(s.comp)+n.mut+s.mut+n.pred+s.pred+bet+close+neigh+ec+hub+pr, family = "gaussian", data = mydat, na.action = "na.fail")

fit3 <- glm(CI.ivary~n.comp+abs(s.comp)+n.mut+s.mut+n.pred+s.pred+bet+close+neigh+ec+hub+pr, 
            family = "gaussian", data = mydat, na.action = "na.fail")
fit4 <- glm(CI.pers~n.comp+abs(s.comp)+n.mut+s.mut+n.pred+s.pred+bet+close+neigh+ec+hub+pr,
            family = "gaussian", data = mydat, na.action = "na.fail")
fit5 <- glm(CI.eig~n.comp+abs(s.comp)+n.mut+s.mut+n.pred+s.pred+bet+close+neigh+ec+hub+pr, family = "gaussian", data = mydat, na.action = "na.fail")


flist.alt <- list(fit1, fit3, fit4, fit5)
fBnames <- c("abund", "cvinit", "pers", "eig")

dfall.alt <- multimod(flist.alt, fBnames)

ggplot(dfall.alt) + geom_segment(aes(x = lower, y = met, xend = upper, yend = met, col = sig)) + geom_vline(aes(xintercept = 0)) + 
  geom_point(aes(x = coef, y = met, col = sig)) + facet_wrap(~mod, scales = "free_x") + 
  scale_color_manual(name = "Significant", values = c("grey", "blue")) + xlab("Value") + ylab("Variable") + theme_bw()


########################
########################
################################
################################

eqsp <- rep(sapply(eqcomm, length), sapply(eqcomm, length))


eqs1 <- lapply(1:sum(use), function(x) cbind(sp = spp[[x]], eq = r2[use][[x]][1000,-1], x))
eqs2 <- do.call(rbind, eqs1)


####################################
####################################

CI.abund <- ((mydat$delta.biom/rep(sapply(eq.abund2, mean), sapply(eqcomm, length))[ccak]))
CI.ivary <- (((icv - mydat$m.init.vary)/icv))
CI.pers <- (((neq - mydat$pers)/neq) * (1/eab))
CI.eig <- (((evi - mydat$eig)/evi) * (1/eab))



subAb <- abs(mydat$delta.biom/rep(sapply(eq.abund2, mean), sapply(eqcomm, length))[ccak]) > 1
subIv <- ((icv - mydat$m.init.vary)/icv) > .1
subPe <- ((neq - mydat$pers)/neq) > 1
subEi <- ((evi - mydat$eig)/evi) > 1

newd2 <- mydat
newd2$G <- destab#G1
newd2$G2 <- CI.pers#G2
newd2$G3 <- CI.abund#G3
newd2$G4 <- CI.eig#G4
newd2$G5 <- CI.ivary#G5

fitCI <- glm(G~n.comp+abs(s.comp)+n.mut+s.mut+n.pred+s.pred+bet+close+neigh+ec+hub+pr, data = newd2[-newd2$pers == 0,], family = "binomial", na.action = "na.fail")
fitCI2 <- glm(G2~n.comp+abs(s.comp)+n.mut+s.mut+n.pred+s.pred+bet+close+neigh+ec+hub+pr, data = newd2[subPe,], family = "gaussian", na.action = "na.fail")
fitCI3 <- glm(G3~n.comp+abs(s.comp)+n.mut+s.mut+n.pred+s.pred+bet+close+neigh+ec+hub+pr, data = newd2[subAb,], family = "gaussian", na.action = "na.fail")
fitCI4 <- glm(G4~n.comp+abs(s.comp)+n.mut+s.mut+n.pred+s.pred+bet+close+neigh+ec+hub+pr, data = newd2[subEi,], family = "gaussian", na.action = "na.fail")
fitCI5 <- glm(G5~n.comp+abs(s.comp)+n.mut+s.mut+n.pred+s.pred+bet+close+neigh+ec+hub+pr, data = newd2[subIv,], family = "gaussian", na.action = "na.fail")

flist <- list(fitCI, fitCI2, fitCI3, fitCI4, fitCI5)
fnames <- c("destab", "pers", "abund", "eig", "ivary")

dfall <- multimod(flist, fnames)

ggplot(dfall) + geom_segment(aes(x = lower, y = met, xend = upper, yend = met, col = sig)) + geom_vline(aes(xintercept = 0)) + 
  geom_point(aes(x = coef, y = met, col = sig)) + facet_wrap(~mod, scales = "free_x") + 
  scale_color_manual(name = "Significant", values = c("grey", "blue")) + xlab("Value") + ylab("Variable") + theme_bw()
ggsave(filename = "~/Desktop/modelpar.jpeg", width = 7, height = 5)

ggplot(dfall) + geom_segment(aes(x = lower, y = met, xend = upper, yend = met, col = sig)) + geom_vline(aes(xintercept = 0)) + 
  geom_point(aes(x = coef, y = met, col = sig)) + facet_wrap(~mod) + 
  scale_color_manual(name = "Significant", values = c("grey", "blue")) + scale_x_continuous(limits = c(-5, 5)) +
  xlab("Value") + ylab("Variable") + theme_bw()
ggsave(filename = "~/Desktop/modelparZOOM.jpeg", width = 7, height = 5)


ggplot(dfall, aes(x = met, y = impt, fill = sig)) + geom_bar(stat = "identity") + scale_fill_manual(name = "Significant", values = c("grey", "blue")) + facet_wrap(~mod) + ylab("Importance") + xlab("Variable") + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(filename = "~/Desktop/parimpt.jpeg", width = 7, height = 5)
####################################
####################################
spi2 <- data.frame(spin2[ccak,], G = newd2$G, mydat)
fitCI.1 <- glm(G~pos1+neg1+spos1+sneg1+pos2+neg2+spos1+sneg2, data = spi2, family = "binomial", na.action = "na.fail")
fitCI2.1 <- glm(pers~pos1+neg1+spos1+sneg1+pos2+neg2+spos1+sneg2, data = spi2, family = "poisson", na.action = "na.fail")
fitCI3.1 <- glm(delta.biom~pos1+neg1+spos1+sneg1+pos2+neg2+spos1+sneg2, data = spi2, family = "gaussian", na.action = "na.fail")
fitCI4.1 <- glm(eig~pos1+neg1+spos1+sneg1+pos2+neg2+spos1+sneg2, data = spi2, family = "gaussian", na.action = "na.fail")
fitCI5.1 <- glm(m.init.vary~pos1+neg1+spos1+sneg1+pos2+neg2+spos1+sneg2, data = spi2, family = "gaussian", na.action = "na.fail")

flist.1 <- list( fitCI2.1, fitCI3.1, fitCI4.1, fitCI5.1)
fnames <- c("destab", "pers", "abund", "eig", "ivary")

dfall.1 <- multimod(flist.1, fnames)

ggplot(dfall.1) + geom_segment(aes(x = lower, y = met, xend = upper, yend = met, col = sig)) + geom_vline(aes(xintercept = 0)) + 
  geom_point(aes(x = coef, y = met, col = sig)) + facet_wrap(~mod, scales = "free_x") + 
  scale_color_manual(name = "Significant", values = c("grey", "blue")) + xlab("Value") + ylab("Variable") + theme_bw()

#####################################
#####################################
#####################################


head(mydat)



#####################################
#####################################
#####################################
# Alternatives for keystones
# keystone as species whose removal leads to at more than 1 extinction
eqnum <- rep(sapply(eqcomm, length), sapply(eqcomm, length))[ccak]
mydat$e <- eqnum-1-mydat$pers > 2
gfit <-glm(e~n.comp+n.mut+n.pred+n.amen+n.com+s.comp+s.mut+bet+close+neigh+ec+pr, family = "binomial", data = mydat, na.action = "na.fail")
gfit <-glm(e~n.comp+n.mut+n.pred+bet, family = "binomial", data = mydat, na.action = "na.fail")
cv.binary(gfit)
confint(model.avg(dredge(gfit)))

ks2[[1]]
cost1 <- function(r, pi = 0) mean(abs(r-pi) > 0.5)


mydat1 <- data.frame(e = mydat$e, deg = unlist(lapply(allg, degree))[ccak])
fitlda1 <- lda(e~deg, data = data.frame(e = mydat$e, deg = unlist(lapply(allg, degree))[ccak]))
fitlda2 <- lda(e~n.comp+n.mut+n.pred+n.amen+n.com+s.comp+s.mut+bet+close+neigh+ec+pr, data = mydat)

sam <- c(sample(which(mydat$e), sum(mydat$e)/2), sample(which(!mydat$e), sum(!mydat$e)/2))
traindat <- mydat[sam,]
testdat <- mydat[-sam,]

ldatrain <- lda(e~n.comp+n.mut+n.pred+n.amen+n.com+s.comp+s.mut+bet+close+neigh+ec+pr, data = traindat)
plda <- predict(ldatrain, newdata = testdat)

confusion <- function(actual, predicted, names = NULL, printit = TRUE,
                       prior = NULL) {
   if (is.null(names))
     names <- levels(actual)
     tab <- table(actual, predicted)
     acctab <- t(apply(tab, 1, function(x) x/sum(x)))
     dimnames(acctab) <- list(Actual = names, "Predicted (cv)" = names)
     if (is.null(prior)) {
       relnum <- table(actual)
       prior <- relnum/sum(relnum)
       acc <- sum(tab[row(tab) == col(tab)])/sum(tab)
       }
     else {
       acc <- sum(prior * diag(acctab))
       names(prior) <- names
       }
     if (printit)
       print(round(c("Overall accuracy" = acc, "Prior frequency" = prior),
                     4))
     if (printit) {
       cat("\nConfusion matrix", "\n")
       print(round(acctab, 4))
       }
     invisible(acctab)
}

confusion(testdat$e, plda$class)
#####################################
#####################################
#####################################
library(reshape2)

getShPath <- function(mat, ks, ksA, eqa){
  test <- mat
  diag(test) <- 0
  test <- melt(test)
  test <- test[test[,3] != 0,]
  g1 <- graph.edgelist(as.matrix(test)[,1:2])
  
  ks <- lapply(ks, function(x) if(sum(is.na(x))>0){rep(NA, length(eqa))}else{x})
  
  spath1 <- sapply(1:length(V(g1)), function(x) sapply(get.shortest.paths(g1, x, mode = "out", output = "epath")$epath, length))
  spath2 <- sapply(1:length(V(g1)), function(x) sapply(get.shortest.paths(g1, x, mode = "out", weights = abs(test[,3]), output = "epath")$epath, length))
  dat1 <- lapply(1:nrow(spath1), function(x) cbind(from = rep(x, length(eqa)-1), to = (1:nrow(spath1))[-x], spuw = spath1[x,-x], spw = spath2[x,-x], d1e = ksA[x,-c(1,(x+1))], dstrt = eqa[-x], tte= ks[[x]][-1]))
  dat2 <- do.call(rbind, dat1)
  
  return(dat2)
}
x=2
mat = matuse[[x]]
ks = ks4[[x]]
ksA = ks5[[x]]
eqa = eq.abund2[[x]]

spaths <- lapply(1:sum(use), function(x) getShPath(matuse[[x]], ks4[[x]], ks5[[x]], eq.abund2[[x]]))
spaths1 <- lapply(1:sum(use), function(x) cbind(spaths[[x]], comm = x))
spaths2 <- do.call(rbind, spaths1)
plot(aggregate(spaths2[,"spw"], list(spaths2[,"tte"]), mean))
plot(spaths2[,"spw"]~spaths2[,"tte"])


mat = matuse[[1]]
ks = ks3[[1]]
ksB = ks4[[9]]
getShPath.p <- function(mat, ks, ksB){
  test <- mat
  diag(test) <- 0
  test <- melt(test)
  test <- test[test[,3] != 0,]
  g1 <- graph.edgelist(as.matrix(test)[,1:2])
  
  spath1 <- sapply(1:length(V(g1)), function(x) sapply(get.shortest.paths(g1, x, mode = "out", output = "epath")$epath, length))
  ksB <-lapply(ksB, function(x) if(sum(is.na(x)) >= 1){rep(NA, nrow(spath1))}else{x})
  dat1 <- lapply(1:nrow(spath1), function(x) cbind(from = rep(x, length(spath1[x,-x])), to = (1:nrow(spath1))[-x], path = spath1[x,-x], dB = ks[x,-x],
                                                   tte = ksB[[x]][-1]))
  dat2 <- do.call(rbind, dat1)
  
  return(dat2)
}


spaths.p <- lapply(1:sum(use), function(x) getShPath.p(matuse[[x]], ks3[[x]], ks4[[x]]))
spaths.p <- lapply(1:length(spaths.p), function(x) cbind(comm = x, spaths.p[[x]]))
spathsp <- do.call(rbind, spaths.p)
hist(spathsp[spathsp[,3] == 0,2])


path2first <- function(mat, t2e, ks){
  test <- mat
  diag(test) <- 0
  test <- melt(test)
  test <- test[test[,3] != 0,]
  g1 <- graph.edgelist(as.matrix(test)[,1:2])
  
  extinctions <- which(apply(ks, 1, function(x) sum(!x, na.rm = T))!=0)
  firstE <- sapply(t2e[extinctions], which.min)
  
  pl <- c()
  for(i in 1:length(extinctions)){
    spa1 <- get.shortest.paths(g1, extinctions[i], firstE[i], mode = "out", output = "epath")$epath[[1]]
    pl[i] <- length(spa1)
  }
  return(pl)
}

hist(unlist(lapply(1:sum(use), function(x) path2first(matuse[[x]], ks4[[x]], ks3[[x]]))))

path <- spa1[[1]]
mat <- matuse[[5]]
ipath <- function(mat, ks){
  test <- mat
  #diag(test) <- 0
  test <- melt(test)
  test <- test[test[,3] != 0,]
  g1 <- graph.edgelist(as.matrix(test)[,1:2])
  
  pathints <- list()
  #extinctions <- which(apply(ks, 1, function(x) sum(!x, na.rm = T))!=0)
  for(i in 1:length(V(g1))){
    
    #spa1 <- get.shortest.paths(g1, extinctions[i], which(!ks[extinctions[i],])[which(!ks[extinctions[i],]) %in% V(g1)],
    #                           mode = "out", output = "epath")$epath
    spa1 <- get.shortest.paths(g1, i, mode = "out", output = "epath")$epath

    vspa <- lapply(spa1, as.vector)
    pathints[[i]] <- lapply(vspa, function(x) test[x,3])
    #print(i)

  }
  return(pathints)
}

pfun <- function(g, x){
  p1 <- c()
  for(i in 1:(length(x)-1)){
    p1[i] <- g[g$Var1 == x[i] & g$Var2 == x[i+1], 3]
    #print(i)
  }
  return(p1)
}

ipathAll <- function(mat, ks){
  test <- mat
  diag(test) <- 0
  test <- melt(test)
  test <- test[test[,3] != 0,]
  g1 <- graph.edgelist(as.matrix(test)[,1:2])
  
  pathints <- list()
  extinctions <- which(apply(ks, 1, function(x) sum(!x, na.rm = T))!=0)
  for(i in 1:length(extinctions)){
    spa1 <- all_simple_paths(g1, extinctions[i], which(!ks[extinctions[i],]), mode = "out")
    vspa <- lapply(spa1, as.vector)
    pathints[[i]] <- lapply(vspa, function(x) pfun(test, x))
  }
  return(pathints)
}

mpath <- list()
mpath2 <- list()
for(i in 1:sum(use)){
  #mpath[[i]] <- do.call(rbind,lapply(ipath(matuse[[i]], ks3[[i]]), function(x) cbind(sapply(x, mean),sapply(x, length))))
  
  pth <- ipath(matuse[[i]], ks3[[i]])
  
  dat <- lapply(pth, function(x) t(sapply(x, function(y){c(nNeg = sum(y < 0), mNeg = mean(y[y < 0]), nPos = sum(y > 0), mPos = mean(y[y > 0]))})))
  
  dat2 <- lapply(1:length(dat), function(x) data.frame(dat[[x]][-x,], from = rep(x, nrow(dat[[x]])-1), 
                                                  to = (1:nrow(dat[[x]]))[-x], ext = ks3[[i]][x,-x], abund = ks2[[i]][x,]))
  
  mpath2[[i]] <- do.call(rbind, dat2)
  print(i)
}

a1 <- sapply(allg, function(x) length(V(x)))
eg1 <- list()
for(i in 1:length(allg)){
  eg1[[i]] <- expand.grid((1:a1[i])[-1], 1:a1[i])
}
mp1 <- do.call(rbind, mpath)
mp2 <- do.call(rbind, mpath2)
mp2[is.nan(mp2)] <- 0

head(mp2)
plot(mp2[,c(1,3)])

plot(mp2[,2:1])
#####################################
ks3[[2]][walktrap.community(g1, weights = abs(test[,3]))[3][[1]],walktrap.community(g1, weights = abs(test[,3]))[3][[1]]]
ks3[[2]][walktrap.community(g1, weights = abs(test[,3]))[3][[1]],-walktrap.community(g1, weights = abs(test[,3]))[3][[1]]]
plot(induced_subgraph(g1, walktrap.community(g1, weights = abs(test[,3]))[3][[1]]))
#####################################
#####################################