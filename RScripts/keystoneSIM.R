###
### CHECK FOR KEYSTONE SPECIES
###

t.key <- Sys.time()                                                 # note time keystone simm starts



# quick test of function
#system.time(
#ks2 <- keystone(2, dyn = r2[use], eqcomm, mats)
#)


# Simulated removal for each species in each equilibrium community
ks1 <- list()
ks2 <- list()
ks3 <- list()
for(i in 1:sum(use)){
  KS <- keystone(i, dyn = r2[use], eqcomm, mats, growth)            # keystone species simulation
  ks1[[i]] <- KS[[1]]                                               # biomass variability and persistence
  ks2[[i]] <- KS[[2]]                                               # change in spp biomass with removal
  ks3[[i]] <- KS[[3]]                                               # who went extinct
  
  cat(paste("\n ------------------|   ", i, "   |------------------ \n"))
}


t.end <- Sys.time()                                                 # what time does the simulation end
t.end - t.start                                                     # total time spent simulating from start to finish
