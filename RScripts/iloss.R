library(vegan)
library(data.table)
library(igraph)
library(deSolve)
library(dplyr)
library(ggplot2)


SteinInt <- read.csv("~/Desktop/GitHub/microbial-dyn/Data/ecomod-ints.csv", row.names = 1)
#SteinInt <- read.csv("C:/Users/jjborrelli/Desktop/GitHub/microbial-dyn/Data/ecomod-ints.csv", row.names = 1)
INTs <- c(SteinInt[upper.tri(SteinInt)],SteinInt[lower.tri(SteinInt)])


t0 <- Sys.time()
Ns <- rep(c(10, 30, 50, 70, 90), each = 10)
Cs <- rep(c(.05, .15, .25), each = 10)
t.final <- 4000

nc <- expand.grid(Ns, Cs)
res <- list()
inds <- list()
for(i in 1:nrow(nc)){
  eqc <- get_eqcomm(nc[i,1], nc[i,2], INTs, 1:t.final)
  
  cv10 <- apply(eqc$comm.dyn[1:10,], 2, function(x) sd(x)/mean(x))
  cv100 <- apply(eqc$comm.dyn[1:100,], 2, function(x) sd(x)/mean(x))
  cv1000 <- apply(eqc$comm.dyn[1:1000,], 2, function(x) sd(x)/mean(x))
  
  upptr.i <- mean(eqc$comm.mat[upper.tri(eqc$comm.mat)][eqc$comm.mat[lower.tri(eqc$comm.mat)] != 0])
  lowtr.i <- mean(eqc$comm.mat[lower.tri(eqc$comm.mat)][eqc$comm.mat[lower.tri(eqc$comm.mat)] != 0])
  
  gr <- eqc$init.parms$alpha
  iD <- vegan::diversity(eqc$comm.dyn[1,])
  eq.exp <- -eqc$init.parms$alpha/diag(eqc$comm.mat)
  
  if(any(is.nan(eqc$comm.dyn))){
    res[[i]] <- data.frame(co = i, N.ext = NA, t.bio = NA, exp.bio = eq.exp, ShD = NA, iShD = iD, eig.i = NA, eig.f = NA,
                           N = NA, C = NA, Lin = NA, Lfi = NA, gr.in = mean(gr), di.in = mean(diag(eqc$comm.mat)), (itypes(eqc$comm.mat)),
                           gr.fi = NA, di.fi = NA, (itypes(eqc$comm.mat)), ut.in = upptr.i, ut.fi = NA, lt.in = lowtr.i, lt.fi = NA, 
                           cv10 = median(cv10), cv100 = median(cv100), cv1k = median(cv1000))
    next
  }
  
  eqc2 <- eqc$comm.mat[eqc$comm.dyn[t.final,] > 0, eqc$comm.dyn[t.final,] > 0]
  
  abund <- eqc$comm.dyn[t.final,]
  abund[abund < 0] <- 0
  
  ex <- sum(abund == 0)
  wex <- 1*(abund != 0)
  
  D <- vegan::diversity(abund)
  
  Lin <- (sum(eqc$comm.mat != 0) - nrow(eqc$comm.mat))
  Lfi <- (sum(eqc2 != 0) - nrow(eqc2))
  
  jf1 <- rootSolve::jacobian.full(y = eqc$comm.dyn[1,], func = lvmod, parms = eqc$init.parms)
  jf2 <- rootSolve::jacobian.full(y = eqc$comm.dyn[t.final,eqc$comm.dyn[t.final,] != 0], 
                                  func = lvmod, parms = list(alpha = gr[eqc$comm.dyn[t.final,] != 0], m = eqc2))
  
  eig1 <- max(Re(eigen(jf1)$values))
  eig2 <- max(Re(eigen(jf2)$values))
  
  upptr.f <- mean(eqc2[upper.tri(eqc2)][eqc2[upper.tri(eqc2)] != 0])
  lowtr.f <- mean(eqc2[lower.tri(eqc2)][eqc2[lower.tri(eqc2)] != 0])
  
  
  res[[i]] <- data.frame(co = i, N.ext = ex, t.bio = mean(abund), exp.bio = mean(eq.exp), ShD = D, iShD = iD, eig.i = eig1, eig.f = eig2,
                         N = nc[i,1], C = nc[i,2], Lin, Lfi, 
                         gr.in = mean(gr), di.in = mean(diag(eqc$comm.mat)), (itypes(eqc$comm.mat)),
                         gr.fi = mean(gr[eqc$comm.dyn[t.final,] != 0]), di.fi = mean(diag(eqc2)), (itypes(eqc2)),
                         ut.in = upptr.i, ut.fi = upptr.f, lt.in = lowtr.i, lt.fi = lowtr.f, cv10 = median(cv10), cv100 = median(cv100), cv1k = median(cv1000))
  
  i1 <- itypes.sp(eqc$comm.mat)
  im <- matrix(0, nrow = nc[i,1], ncol = 5)
  im[as.logical(wex), ] <- itypes.sp(eqc2)
  colnames(im) <- colnames(i1)
  
  g1 <- graph.adjacency(sign(abs(eqc$comm.mat)), diag = F)
  g2 <- graph.adjacency(sign(abs(eqc2)), diag = F)
  b1 <- betweenness(g1)
  b2 <- vector("numeric", length = nc[i,1])
  b2[as.logical(wex)] <- betweenness(g2)
  
  is.out <- colMeans(eqc$comm.mat)
  is.in <- rowMeans(eqc$comm.mat)
  is.out2 <- vector("numeric", nc[i,1])
  is.in2 <- vector("numeric", nc[i,1])
  is.out2[as.logical(wex)] <- colMeans(eqc2)
  is.in2[as.logical(wex)] <- rowMeans(eqc2)
  
  inds[[i]] <- data.frame(co = i, n = nc[i,1], C = nc[i,2], ext = wex, bio = abund, ebio = eq.exp, self = diag(eqc$comm.mat), gro = gr, 
                          bet.i = b1, bet.f = b2, i1, im, is.out1 = is.out, is.in1 = is.in, is.out2 = is.out2, is.in2 = is.in2, 
                          cv10 = cv10, cv100 = cv100, cv1k = cv1000)
  
  cat(rep(c(i, "\n"), 10))
}

webres <- rbindlist(res)
sppres <- rbindlist(inds)

t1 <- Sys.time()
t1-t0

save.image("~/Desktop/simstab.Rdata")

plot(ext~bet.i, data = sppres)
cartA <- rpart::rpart(ext~n+C+self+gro+bet.i+comp+mut+pred+amens+comm+is.out1+is.in1, data = sppres, method = "class")
plot(cartA)
text(cartA, use.n=TRUE, all=TRUE, cex=.8)

cartB <- rpart::rpart(N.ext~N+C+Lin+iShD+di.in+gr.in+eig.i+comp+mut+pred+amens+comm, data = webres, method = "anova")
plot(cartB)
text(cartB, use.n=TRUE, all=TRUE, cex=.8)


dfA <- data.frame(P = webres$N.ext/webres$N, C = webres$C, iShD = webres$iShD,t(apply(select(webres, comp:comm),1,function(x) x/sum(x))))
dfB <- data.frame(t(apply(select(webres, comp:comm),1,function(x) x/sum(x))), t(apply(select(webres, comp.1:comm.1),1,function(x) x/sum(x))))
head(melt(t(apply(select(webres, comp:comm),1,function(x) x/sum(x)))))
test <- (data.frame(melt(t(apply(select(webres, comp.1:comm.1),1,function(x) x/sum(x)))), melt(t(apply(select(webres, comp:comm),1,function(x) x/sum(x))))))
test2 <- (data.frame(melt(t(apply(select(sppres[as.logical(sppres$ext)], comp.1:comm.1),1,function(x) x/sum(x)))), melt(t(apply(select(sppres[as.logical(sppres$ext)], comp:comm),1,function(x) x/sum(x))))))

p1 <- ggplot(test, aes(x = value.1, y = value)) + geom_point() + geom_abline(slope = 1, intercept = 0) 
p1 <- p1 + facet_wrap(~Var2.1) + xlab("Initial Proportion") + ylab("Final Proportion") + theme_bw()
p1

s1 <- sample(1:nrow(dfA), ceiling(nrow(dfA)/2))
fit <- (lm(P~C+iShD+mut, data = dfA[s1,]))
summary(fit)

mean((dfA$P[-s1] - predict(fit, newdata = dfA[(1:nrow(dfA))[-s1],]))^2)


p2 <- ggplot(test2, aes(x = value.1, y = value)) + geom_point() + geom_abline(slope = 1, intercept = 0) 
p2 <- p2 + facet_wrap(~Var2.1) + xlab("Initial Proportion") + ylab("Final Proportion") + theme_bw()
p2
plot(sppres$comp, sppres$ext)
points(glm(sppres$ext~sppres$comp)$fitted.values~sppres$comp, col = "blue", pch = 20)
