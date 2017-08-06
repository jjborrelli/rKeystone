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


key_effect2 <- function(init, i, mod = NULL,  plots = TRUE){ 
  nt <- nrow(init$dyn1)
  spp <- init$dyn1[nt,-1] > 10^-10
  mat <- init$m[spp,spp]
  alp <- init$alpha[spp]
  ia <- init$dyn1[nt,-1][init$dyn1[nt,-1] > 10^-10]
  pars <- list(alpha = alp, m = mat) 
  
  div <- c()
  mmch <- list()
  per <- matrix(nrow = nrow(mat), ncol = 3)
  colnames(per) <- c("ext", "tteq", "ls")
  bdiff <- matrix(nrow = nrow(mat), ncol = 7)
  colnames(bdiff) <- c("tot", "npos", "mpos", "maxp", "nneg", "mneg", "maxn")
  cvx <- list()
  
  #for(i in 1:nrow(mat)){
    ia <- init$dyn1[nt,-1][spp]
    ia[i] <- 0
    
    dyn <-(ode(ia, times = 1:nt, func = lvm, parms = pars, events = list(func = ext1, time =  1:nt)))
    
    if(nrow(dyn) != nt){
      cvx[[i]] <- rbind(rep(NA, 3), rep(NA, 3))
      next
    }
    
    if(plots){
      matplot(dyn[1:100,-1], typ = "l", main = i)
    }
    
    per <- persist(dyn, nt, alp, mat)
    bdiff <- biodiff(dyn, nt)
    cvx <- cvar(dyn, 100)
    div <- vegan::diversity(dyn[nt,-1][dyn[nt,-1] > 10^-10])
    
    if(!is.null(mod)){
      ch <- dyn[nt,-1][dyn[1,-1] > 0] - dyn[1,-1][dyn[1,-1] > 0]
      mmch[[i]] <- cbind(aggregate(ch, list(mod[-i]), mean), rmod = mod[i], node = i)
    }
  #}
  cova <- (apply(cvx, 2, median, na.rm = T))
  impact <- c(per, bdiff, cova, div = div)
  
  if(!is.null(mod)){
    mod.imp <- do.call(rbind, mmch)
    return(list(impact, mod.imp, dyn))
  }
  
  return(impact)
}



ilist <- list()
for(i in 1:20){
  init <- isim(S = 50, tf = 2000, efun = ext1, idis = "beta", dp1 = 1, dp2 = 4, Rmax = 1, self = 1, plot = TRUE)
  ilist[[i]] <- init
}

matplot(ilist[[9]]$dyn1[,-1], typ = "l")
matplot(ilist[[12]]$dyn1[,-1], typ = "l")
matplot(ilist[[20]]$dyn1[,-1], typ = "l")

init <- ilist[[12]]
mat <- init$m[init$dyn1[2000,-1] > 10^-10,init$dyn1[2000,-1] > 10^-10]
diag(mat) <- 0

nc <- rnetcarto::netcarto(abs(sign(mat)))[[1]]
nc <- nc[order(as.numeric(as.character(nc$name))),]


ke1 <- key_effect2(init, i = 35, mod = nc$module)
#15, 24, 26, 31, 33, 34, 37
# 27 - compensatory
# 30, 35, 36 mostly positive

library(ggraph)


nsize <- ceiling(init$dyn1[1,-1][init$dyn1[1,-1] > 10^-10])+1
ncolor <- rep("black", 50)
g <- graph.adjacency(abs((init$m)),mode = "directed", weighted = T)
elist <- get.edgelist(g)
ecolor <- c();for(i in 1:nrow(elist)){ecolor[i] <- ifelse(init$m[elist[i,1], elist[i,2]] > 0, "blue", "red")}
lay <- create_layout(g, "nicely")
pres <- as.numeric(init$dyn1[2000,-1] > 10^-10)
ggraph(lay) + 
  geom_edge_fan(arrow = arrow(angle = 20,length = unit(.15, "inches"),
                               type = "closed"),  aes(color = ecolor)) + 
  geom_node_point(size = nsize) + theme_void() + theme(legend.position="none")  

ndat <- data.frame(get.edgelist(g)-1)
colnames(ndat) <- c("src", "target")

ggraph(lay, layout = 'kk') + 
  geom_edge_fan(aes(alpha = ..index..), show.legend = FALSE) + 
  geom_node_point()  



nsize <- ceiling(init$dyn1[2000,-1][init$dyn1[2000,-1] > 10^-10])+1

nsizeF <- ceiling(ke1[[3]][2000,-1])+1
ncolor <- rep("black", 37)
nfill <- rep("black", 37)
nfill[26] <- "white"

nfillF <- rep("black", 37)
nfillF[ke1[[3]][2000,-1] < 10^-10] <- rep("white", sum(ke1[[3]][2000,-1] < 10^-10))
g <- graph.adjacency(abs((mat)),mode = "directed", weighted = T)
elist <- get.edgelist(g)
ecolor <- c();for(i in 1:nrow(elist)){ecolor[i] <- ifelse(mat[elist[i,1], elist[i,2]] > 0, "blue", "red")}

#lay2 <- create_layout(g, "nicely")
ggraph(lay2) + 
  geom_edge_link(arrow = arrow(angle = 20,length = unit(.15, "inches"),
                               type = "closed"), edge_width = 1, aes(color = ecolor)) + 
  geom_node_point(size = nsize, shape = 21, fill = nfill, colour = ncolor) + theme_void() + theme(legend.position="none")  
ggsave("D:/key-images/eqnet-rem-sp26-I.svg")

ggraph(lay2) + 
  geom_edge_link(arrow = arrow(angle = 20,length = unit(.15, "inches"),
                               type = "closed"), edge_width = 1, aes(color = ecolor)) + 
  geom_node_point(size = nsizeF, shape = 21, fill = nfillF, colour = ncolor) + theme_void() + theme(legend.position="none")  
ggsave("D:/key-images/eqnet-rem-sp35-F.svg")
#saveRDS(init, "D:/keystone_ex.rds")


spp1 <- as.numeric(colnames(ke1[[3]][,-1]))
dim(init$dyn1)
tm1 <- matrix(0, nrow = 2000, ncol = 50)
tm1[,spp1] <- ke1[[3]][,-1]
alldyn <- rbind(init$dyn1[,-1], tm1)
matplot(alldyn, typ = "l", ylab = "Biomass", xlab = "Time")
matplot(alldyn[2000:2100,], typ = "l", ylab = "Biomass", xlab = "Time")


ke1 <- key_effect2(init, i = 35, mod = nc$module)
spp1 <- as.numeric(colnames(ke1[[3]][,-1]))#[ke1[[3]][2000,-1] > 10^-10]
ity <- itypes(init$m[spp1, spp1])
df2 <- data.frame(typ = names(ity),num = ity/sum(ity))
df35F <- data.frame(typ = names(ity),num = ity/sum(ity))

df.all <- rbind(cbind(df1, Net = "Initial"), cbind(df2, Net = "Final"), cbind(df26I, Net = "Initial"), cbind(df26F, Net = "Final"), cbind(df35I, Net = "Initial"), cbind(df35F, Net = "Final"))
df.all$Rem <- rep(c("Full", "Removal 1", "Removal 2"), each = 10)
ggplot(df.all, aes(x = Net, y = num, fill = typ)) + geom_bar(stat = "identity") + facet_wrap(~Rem) + 
  scale_fill_manual(labels = c("Amensalism", "Commensalism", "Competition", "Mutualism", "Predation"), values =   c("#999999", "#56B4E9", "#009E73", "#0072B2", "#E69F00")) + 
  labs(x = "", y = "Relative Frequency", fill = "Interaction") + theme_bw()
ggsave("D:/key-images/ints.svg")
