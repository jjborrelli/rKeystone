### SAVED WORK
# last saved 10-25-16
# alt save 10-26-16
# ms save 11-8-16 == example3
# save.image("~/Desktop/simul-example4.Rdata") 
# load("~/Desktop/simul-example3.Rdata")

###
### ANALYSIS
###


itySP2 <- do.call(rbind, itySP)                                     # get species level participation in interaction types
istrSP2 <- do.call(rbind, istrSP)                                   # get species level interaction type strengths
allks <- do.call(rbind, lapply(ks1, function(x) x[,1:4]))           # all biomass, variation, and persistence

itySP3 <- t(apply(itySP2, 1, function(x) x/sum(x)))
itySP3[is.nan(itySP3)] <- 0

# put all data together in single matrix
allks <- cbind(allks, eig = unlist(eigkey), sp.id = unlist(eqcomm), n.comp = itySP2[,1], n.mut = itySP2[,2], n.pred = itySP2[,3], 
               n.amen = itySP2[,4], n.com = itySP2[,5], 
               s.comp = istrSP2[,1], s.mut = istrSP2[,2], s.pred = istrSP2[,3], bet = unlist(betw), close = unlist(clocent),
               neigh = unlist(g.neighbors2),  ec = unlist(ecent), hub = unlist(hscore), pr = unlist(p.rank))
ccak <- complete.cases(allks)                                       # only use complete cases

dim(allks[ccak,])

eq.abund <- unlist(lapply(dyn, function(x) x[1000,-1][x[1000,-1] > 0]))[ccak]
eq.abund2 <- (lapply(dyn, function(x) x[1000,-1][x[1000,-1] > 0]))

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
destab <- (mydat$eig >= evinit)

eab <- unlist(lapply(eq.abund2, function(x) x/sum(x)))[ccak]
icv <- rep(sapply(cv.eq, mean), sapply(eqcomm, length))[ccak]
neq <- rep(sapply(eqcomm, length), sapply(eqcomm, length))[ccak]
evi <- rep(ev.init, sapply(eqcomm, length))[ccak]



## Modeling
#### with MuMIn package


mydat <- as.data.frame(allks[ccak,])

fit1 <- glm(delta.biom~n.comp+n.mut+n.pred+s.comp+s.mut+s.pred+bet+close+neigh+ec+hub+pr, family = "gaussian", data = mydat, na.action = "na.fail")
fit2 <- glm(mean.vary~n.comp+n.mut+n.pred+s.comp+s.mut+s.pred+bet+close+neigh+ec+hub+pr, family = "gaussian", data = mydat, na.action = "na.fail")
fit3 <- glm(m.init.vary~n.comp+n.mut+n.pred+s.comp+s.mut+s.pred+bet+close+neigh+ec+hub+pr, 
            family = "gaussian", data = mydat, na.action = "na.fail")
fit4 <- glm(cbind(pers,rep(sapply(eqcomm, length), sapply(eqcomm, length))[ccak])~n.comp+n.mut+n.pred+s.comp+s.mut+s.pred+
              bet+close+neigh+ec+hub+pr,
            family = "binomial", data = mydat, na.action = "na.fail")
fit5 <- glm(eig~n.comp+n.mut+n.pred+s.comp+s.mut+s.pred+bet+close+neigh+ec+hub+pr,
            family = "gaussian", data = mydat, na.action = "na.fail")

d1.fit <- dredge(fit1)
d2.fit <- dredge(fit2)
d3.fit <- dredge(fit3)
d4.fit <- dredge(fit4)
d5.fit <- dredge(fit5)


head(d1.fit)
head(d2.fit)
head(d3.fit)
head(d4.fit)
head(d5.fit)


dmat1 <- matrix(c(colMeans(d1.fit[d1.fit$delta < 2,], na.rm = T),colMeans(d2.fit[d2.fit$delta < 2,], na.rm = T),colMeans(d3.fit[d3.fit$delta < 2,], na.rm = T),colMeans(d4.fit[d4.fit$delta < 2,], na.rm = T),colMeans(d5.fit[d5.fit$delta < 2,], na.rm = T)), nrow = 5, byrow = T)
colnames(dmat1) <- names(colMeans(d4.fit[d4.fit$delta < 2,]))
rownames(dmat1) <- c("biomass", "meanvary", "initvary", "persist", "eigen")
dmat1
dmat2 <- matrix(c((d1.fit[1,]),(d2.fit[1,]),(d3.fit[1,]),(d4.fit[1,]),(d5.fit[1,])), nrow = 5, byrow = T)
colnames(dmat2) <- names(colMeans(d4.fit[d4.fit$delta < 2,]))
rownames(dmat2) <- c("biomass", "meanvary", "initvary", "persist", "eigen")
dmat2

#####################################
#####################################
fit1 <- glm(delta.biom~n.comp+n.mut+n.pred+n.amen+n.com+s.comp+s.mut+s.pred, family = "gaussian", data = mydat, na.action = "na.fail")
fit2 <- glm(mean.vary~n.comp+n.mut+n.pred+n.amen+n.com+s.comp+s.mut+s.pred, family = "gaussian", data = mydat, na.action = "na.fail")
fit3 <- glm(m.init.vary~n.comp+n.mut+n.pred+n.amen+n.com+s.comp+s.mut+s.pred, 
            family = "gaussian", data = mydat, na.action = "na.fail")
fit4 <- glm(cbind(pers,rep(sapply(eqcomm, length), sapply(eqcomm, length))[ccak])~n.comp+n.mut+n.pred+n.amen+n.com+s.comp+s.mut+s.pred,
            family = "binomial", data = mydat, na.action = "na.fail")
fit5 <- glm(eig~n.comp+n.mut+n.pred+n.amen+n.com+s.comp+s.mut+s.pred, family = "gaussian", data = mydat, na.action = "na.fail")


d1.fit <- dredge(fit1)
d2.fit <- dredge(fit2)
d3.fit <- dredge(fit3)
d4.fit <- dredge(fit4)
d5.fit <- dredge(fit5)


model.avg(d1.fit, subset = delta < 5)
model.avg(d2.fit, subset = delta < 5)
model.avg(d3.fit, subset = delta < 5)
model.avg(d4.fit, subset = delta < 5)
model.avg(d5.fit, subset = delta < 5)

dmat1 <- matrix(c(colMeans(d1.fit[d1.fit$delta < 2,], na.rm = T),colMeans(d2.fit[d2.fit$delta < 2,], na.rm = T),colMeans(d3.fit[d3.fit$delta < 2,], na.rm = T),colMeans(d4.fit[d4.fit$delta < 2,], na.rm = T),colMeans(d5.fit[d5.fit$delta < 2,], na.rm = T)), nrow = 5, byrow = T)
colnames(dmat1) <- names(colMeans(d4.fit[d4.fit$delta < 2,]))
rownames(dmat1) <- c("biomass", "meanvary", "initvary", "persist", "eigen")
dmat1
dmat2 <- matrix(c((d1.fit[1,]),(d2.fit[1,]),(d3.fit[1,]),(d4.fit[1,]),(d5.fit[1,])), nrow = 5, byrow = T)
colnames(dmat2) <- names(colMeans(d4.fit[d4.fit$delta < 2,]))
rownames(dmat2) <- c("biomass", "meanvary", "initvary", "persist", "eigen")
dmat2



########################
# 
CI.abund <- (mydat$delta.biom/rep(sapply(eq.abund2, mean), sapply(eqcomm, length))[ccak]) * (1/eab)
CI.ivary <- ((icv - mydat$m.init.vary)/icv) * (1/eab) 
CI.pers <- ((neq - mydat$pers)/neq) * (1/eab)
CI.eig <- ((evi - mydat$eig)/evi) * (1/eab)


fit1 <- glm(CI.abund~n.comp+n.mut+n.pred+n.amen+n.com+s.comp+s.mut+s.pred, family = "gaussian", data = mydat, na.action = "na.fail")

fit3 <- glm(CI.ivary~n.comp+n.mut+n.pred+n.amen+n.com+s.comp+s.mut+s.pred, 
            family = "gaussian", data = mydat, na.action = "na.fail")
fit4 <- glm(CI.pers~n.comp+n.mut+n.pred+n.amen+n.com+s.comp+s.mut+s.pred,
            family = "gaussian", data = mydat, na.action = "na.fail")
fit5 <- glm(CI.eig~n.comp+n.mut+n.pred+n.amen+n.com+s.comp+s.mut+s.pred, family = "gaussian", data = mydat, na.action = "na.fail")


d1.fit <- dredge(fit1)
d3.fit <- dredge(fit3)
d4.fit <- dredge(fit4)
d5.fit <- dredge(fit5)


head(d1.fit)
head(d3.fit)
head(d4.fit)
head(d5.fit)

########################
########################
################################
################################

eqsp <- rep(sapply(eqcomm, length), sapply(eqcomm, length))


eqs1 <- lapply(1:sum(use), function(x) cbind(sp = spp[[x]], eq = r2[use][[x]][1000,-1], x))
eqs2 <- do.call(rbind, eqs1)


####################################
####################################

CI.abund <- ((mydat$delta.biom/rep(sapply(eq.abund2, mean), sapply(eqcomm, length))[ccak]) * (1/eab))[(mydat$delta.biom/rep(sapply(eq.abund2, mean), sapply(eqcomm, length))[ccak]) < .1]
CI.ivary <- (((icv - mydat$m.init.vary)/icv) * (1/eab))[((icv - mydat$m.init.vary)/icv) < .1]
CI.pers <- (((neq - mydat$pers)/neq) * (1/eab))[((neq - mydat$pers)/neq) < .1]
CI.eig <- (((evi - mydat$eig)/evi) * (1/eab))[((evi - mydat$eig)/evi) < .1]

plot(((neq - mydat$pers)/neq)[((evi - mydat$eig)/evi) < .1]~eab[((evi - mydat$eig)/evi) < .1])


quant1 <- .9
G1 <- (abs(CI.pers) > quantile(abs(CI.pers), probs = quant1) & abs(CI.abund) > quantile(abs(CI.abund), probs = quant1) & (CI.eig) > quantile((CI.eig), probs = quant1) & abs(CI.ivary) > quantile(abs(CI.ivary), probs = quant1))*1
sum(G1)
G2 <- (abs(CI.pers) > quantile(abs(CI.pers), probs = quant1))*1
sum(G2)
G3 <- (abs(CI.abund) > quantile(abs(CI.abund), probs = quant1))*1
sum(G3)
G4 <- ((CI.eig) > quantile((CI.eig), probs = quant1))*1
sum(G4)
G5 <- (abs(CI.ivary) > quantile(abs(CI.ivary), probs = quant1))*1
sum(G5)


G1 <- (abs(CI.pers) > 100 & abs(CI.abund) > 100 & (CI.eig) > 100 & abs(CI.ivary) > 100)
sum(G1)
G2 <- (abs(CI.pers) > 100)*1
sum(G2)
G3 <- (abs(CI.abund) > 100)*1
sum(G3)
G4 <- ((CI.eig) > 100)*1
sum(G4)
G5 <- (abs(CI.ivary) > 100)*1
sum(G5)

t1 <- unlist(lapply(matuse, colSums))[ccak]
t2 <- unlist(lapply(matuse, rowSums))[ccak]
summary(glm(G~t1+t2, data = newd2, family = "binomial", na.action = "na.fail"))


newd2 <- mydat[eab > 10e-5,]
newd2$G <- destab[eab > 10e-5]#G1
newd2$G2 <- CI.pers#G2
newd2$G3 <- CI.abund#G3
newd2$G4 <- CI.eig#G4
newd2$G5 <- CI.ivary#G5

fitCI <- glm(G~n.comp+abs(s.comp)+n.mut+s.mut+n.pred+s.pred, data = newd2, family = "binomial", na.action = "na.fail")
fitCI2 <- glm(G2~n.comp+abs(s.comp)+n.mut+s.mut+n.pred+s.pred, data = newd2[((neq - mydat$pers)/neq) < .1,], family = "gaussian", na.action = "na.fail")
fitCI3 <- glm(G3~n.comp+abs(s.comp)+n.mut+s.mut+n.pred+s.pred, data = newd2[(mydat$delta.biom/rep(sapply(eq.abund2, mean), sapply(eqcomm, length))[ccak]) < .1,], family = "gaussian", na.action = "na.fail")
fitCI4 <- glm(G4~n.comp+abs(s.comp)+n.mut+s.mut+n.pred+s.pred, data = newd2[((evi - mydat$eig)/evi) < .1,], family = "gaussian", na.action = "na.fail")
fitCI5 <- glm(G5~n.comp+abs(s.comp)+n.mut+s.mut+n.pred+s.pred, data = newd2[((icv - mydat$m.init.vary)/icv) < .1,], family = "gaussian", na.action = "na.fail")

dCI <- dredge(fitCI)
dCI2 <- dredge(fitCI2)
dCI3 <- dredge(fitCI3)
dCI4 <- dredge(fitCI4)
dCI5 <- dredge(fitCI5)

ma1 <- model.avg(dCI, subset = delta < 2)
summary(ma1)

# pers
ma2 <- model.avg(dCI2, subset = delta < 2)
summary(ma2)

# abund
ma3 <- model.avg(dCI3, subset = delta < 2)
summary(ma3)

# eig
ma4 <- model.avg(dCI4, subset = delta < 2)
summary(ma4)

# init vary
ma5 <- model.avg(dCI5, subset = delta < 2)
summary(ma5)



length(unlist(eqcomm)[ccak][newd2$G2 == 1])
df0 <- data.frame(confint(ma1, level = .95), rownames(confint(ma1)), ma1$coefficients[1,], "Stab")
df1 <- data.frame(confint(ma2, level = .95), rownames(confint(ma2)), ma2$coefficients[1,], "Persistence")
df2 <- data.frame(confint(ma3, level = .95), rownames(confint(ma3)), ma3$coefficients[1,], "Change in Abundance")
df3 <- data.frame(confint(ma4, level = .95), rownames(confint(ma4)), ma4$coefficients[1,], "Local Stability")
df4 <- data.frame(confint(ma5, level = .95), rownames(confint(ma5)), ma5$coefficients[1,], "Initial Variation")

colnames(df1) <- c("lower", "upper", "met", "coef", "mod")
colnames(df2) <- c("lower", "upper", "met", "coef", "mod")
colnames(df3) <- c("lower", "upper", "met", "coef", "mod")
colnames(df4) <- c("lower", "upper", "met", "coef", "mod")

df1$sig <- df1$lower < 0 & df1$upper < 0 | df1$lower > 0 & df1$upper > 0
df2$sig <- df2$lower < 0 & df2$upper < 0 | df2$lower > 0 & df2$upper > 0
df3$sig <- df3$lower < 0 & df3$upper < 0 | df3$lower > 0 & df3$upper > 0
df4$sig <- df4$lower < 0 & df4$upper < 0 | df4$lower > 0 & df4$upper > 0

df1$impt <- ma2$importance[rownames(df1)]
df2$impt <- ma3$importance[rownames(df2)]
df3$impt <- ma4$importance[rownames(df3)]
df4$impt <- ma5$importance[rownames(df4)]

dfall <- rbind(df1, df2, df3, df4)

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

####################################
####################################
##  Universal Interactions



eqsp <- rep(sapply(eqcomm, length), sapply(eqcomm, length))


eqs1 <- lapply(1:sum(use), function(x) cbind(sp = spp[use][[x]], eq = r2[use][[x]][1000,-1], x))
eqs2 <- do.call(rbind, eqs1)

c1 <- combn(1:181, 2)
nshare <- c()
nshare2 <- c()
d1 <- c()
d2 <- c()
d3 <- c()
d4 <- c()
d5 <- c()
for(i in 1:ncol(c1)){
  nshare[i] <- sum(colSums(eqmat[c1[,i],]) == 2)/sum(colSums(eqmat[c1[,i],]) != 0)
  shsp <- which(colSums(eqmat[c1[,i],]) == 2)
  l1 <- eqs1[c1[,i]]
  e1 <- l1[[1]][l1[[1]][,"sp"] %in% shsp, "eq"]
  e2 <- l1[[2]][l1[[2]][,"sp"] %in% shsp, "eq"]
  
  d1[i] <- dist(rbind(e1, e2))
  d2[i] <- dist(rbind(e1/sum(e1), e2/sum(e2)))
  
  l2 <- lapply(l1, function(x) x[,2] <- x[,2]/sum(x[,2]))
  e1 <- l2[[1]][l1[[1]][,"sp"] %in% shsp]
  e2 <- l2[[2]][l1[[2]][,"sp"] %in% shsp]
  
  nshare2[i] <- sum((e1 + e2)/2)
  
  d3[i] <- dist(rbind(e1, e2))
  d4[i] <- dist(rbind(e1/sum(e1), e2/sum(e2)))
  
  eA <- e1/sum(e1)
  eB <- e2/sum(e2)
  m <- (eA + eB)/2
  
  d5[i] <- sqrt((sum(eA * log10(eA/m)) + sum(eB * log10(eB/m)))/2)
}
range(nshare)
hist(nshare)
par(mfrow = c(2,2))
plot(d1~nshare, main = "Abs Abund")
plot(d2~nshare, main = "Scaled Abs Abund")
plot(d3~nshare, main = "Rel Abund")
plot(d4~nshare, main = "Scaled Rel Abund")

plot(d1~nshare2, main = "Abs Abund")
plot(d2~nshare2, main = "Scaled Abs Abund")
plot(d3~nshare2, main = "Rel Abund")
plot(d4~nshare2, main = "Scaled Rel Abund")

ggplot(data.frame(nshare2, d5), aes(x = nshare2, y = d5)) + geom_point(alpha = 0.1) + geom_smooth()
