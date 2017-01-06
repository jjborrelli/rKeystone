

Bbef <- r2[use][[1]][1000,-1][r2[use][[1]][1000,-1] != 0]
sapply(1:25, function(x) log(KS[[2]][x,]/Bbef[-x]))


# grow tree 
fit <- rpart(Mileage~Price + Country + Reliability + Type, 
             method="anova", data=cu.summary)

fit <- rpart(pers~n.comp+abs(s.comp)+n.mut+s.mut+n.pred+s.pred+bet+close+neigh+ec+hub+pr, 
             method="anova", data=mydat)

printcp(fit) # display the results 
plotcp(fit) # visualize cross-validation results 
summary(fit) # detailed summary of splits

# create additional plots 
par(mfrow=c(1,2)) # two plots on one page 
rsq.rpart(fit) # visualize cross-validation results  	

# plot tree 
plot(fit, uniform=TRUE, 
     main="Regression Tree")
text(fit, use.n=TRUE, all=TRUE, cex=.8)


head(allks)

ks2
eq.abs <- lapply(r2[use], function(x) x[1000,-1][x[1000,-1] != 0])

ri <- list()
for(i in 1:22){
  test <- sapply(1:nrow(ks2[[i]]), function(x) log(ks2[[i]][x,]/eq.abs[[i]][-x]))
  rownames(test) <- 1:nrow(test)
  mt <- melt(test)
  colnames(mt) <- c("Target", "Removed", "Bafter")
  
  #eqa <- c()
  #for(j in 1:nrow(mt)){
  #  eqa[j] <- eq.abs[[i]][-mt[j,2]][mt[j,1]]
  #}
  #ri[[i]] <- cbind(mt, Bbef = eqa)
  ri[[i]] <- mt
}