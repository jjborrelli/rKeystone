library(vegan)


Ns <- rep(c(10, 20, 30, 40, 50, 60), each = 10)
Cs <- rep(c(.05, .1, .15, .2, .25, .3), each = 10)


nc <- expand.grid(Ns, Cs)
res <- list()
for(i in 1:nrow(nc)){
  eqc <- get_eqcomm(nc[i,1], nc[i,2], INTs, 1:2000)
  eqc2 <- eqc$comm.mat[eqc$comm.dyn[2000,] != 0, eqc$comm.dyn[2000,] != 0]
  gr <- eqc$init.parms$alpha
  
  abund <- eqc$comm.dyn[2000,]
  ex <- sum(abund == 0)
  eq.exp <- -eqc$init.parms$alpha/diag(eqc$comm.mat)
  iD <- vegan::diversity(eqc$comm.dyn[1,])
  D <- vegan::diversity(abund)
  
  Lin <- (sum(eqc$comm.mat != 0) - nrow(eqc$comm.mat))
  Lfi <- (sum(eqc2 != 0) - nrow(eqc2))
  
  jf1 <- rootSolve::jacobian.full(y = eqc$comm.dyn[1,], func = lvmod, parms = eqc$init.parms)
  jf2 <- rootSolve::jacobian.full(y = eqc$comm.dyn[2000,eqc$comm.dyn[2000,] != 0], func = lvmod, parms = list(alpha = gr[eqc$comm.dyn[2000,] != 0], m = eqc2))
  eig1 <- max(Re(eigen(jf1)$values))
  eig2 <- max(Re(eigen(jf2)$values))
  
  res[[i]] <- data.frame(N.ext = ex, t.bio = mean(abund), exp.bio = mean(eq.exp), ShD = D, iShD = iD, eig.i = eig1, eig.f = eig2,
                         N = nc[i,1], C = nc[i,2], Lin, Lfi, 
                         gr.in = mean(gr), di.in = mean(diag(eqc$comm.mat)), (itypes(eqc$comm.mat)),
                         gr.fi = mean(gr[eqc$comm.dyn[2000,] != 0]), di.fi = mean(diag(eqc2)), (itypes(eqc2)))
  
  cat(rep(c(i, "\n"), 10))
}

rbindlist(res)
