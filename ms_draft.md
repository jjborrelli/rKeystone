# Interactions and impacts of keystone species   
**J. J. Borrelli**

## Introduction

Keystone species are those that have a disproportional impact on the communtiy relative to thier abundance. The textbook case of a keystone species is the sea star _Pisaster_ in tidal communities of the Pacific Northwest (Paine 1966, 1969). _Pisaster_ feeds on highly competitive mussels, clearing precious real-estate in the inter-tidal zone for multiple sedentary species to colonize. Removing _Pisaster_ from these communities leads to the competitive exclusion of most of these species, and the domination of the community by a single mussel. Other examples of keystone species include otters, whose predatory interaction with sea urchins helps stabilize the kelp forest (Estes et al. 1998), and the removal of bass in lakes leads to a trophic cascade affecting smaller fish, grazers and algae abundances (Carpenter 1985, Power 1985, Mittelbach 2006). 

Most studies involving keystone species are focused on a single interaction type, predator-prey. For example, Paine (1980) highlights that _Pisaster_ is effective at reducing competition among mussels because it is a generalist that consumes prey at a range of sizes. The impact that different species may have on their community, however, may be mediated by different types of interaction including competition and mutualism. Moreover, keystone interactions of different types may be variable in the form of the impact to the community.   

Microbial communities offer an exciting avenue for research into the effects of keystone species. With their short generation times, long term community dynamics may be observed over comparatively short human time scales. Advances in metagenomics also allow us to identify microbial taxa and quantify their abundance over time _in vivo_. 

In the human gut there is a thriving microbial community with hundreds of coexisting species. Recent advances in metagenomic sequencing have allowed us to catalogue these species and describe the variation in community structure and composition across individuals and within individuals over time. Several factors likely impact the composition of the community. Variation in diet has been shown to influence microbial community composition. This is likely because different microbes are able to more efficiently use different resources (macronutrients). Addtionally, either by hitchhiking on food or other means of transfer, new microbes may invade the community. For example, _Eschericha coli_ may be introduced to the gut via ingestion of contaminated or undercooked foods. The host enviroment can also influence which species are able to coexist in the gut community either through an immune response or mediated through some kind of niche-selection. Finally, the interactions among microbes may both set the boundaries for community composition and drive the response of the community to the other external impacts. 

Despite the many ways the gut microbiome may be affected, longitudinal studies have revealed that their composition tends to be stable over long periods of time. An understanding of how microbes interact may allow us to understand why the human gut microbiome is able to remain stable. Specifically we may want to know whether the role a species plays in the community can tell us about how important they are for the stability of the system. Becasue it is also likely that the interactions among microbes are universal (REF) being able to identify the impact species with different roles have may be generalizable to the larger population.  

Microbes can compete with one another directly for limited resources in the gut (both food and space). Alternatively they may interact via the production of metabolites. These compounds can be either beneficial or detrimental to the growth of other microbial populations. Thus, in microbial communities we may expect to see all five major interaction types: competition (-,-), mutualism (+,+), parasitism/predation (+,-), amensalism (-,0), and commensalism (+,0). While prevailing wisdom on microbial communties suggested that their stability and function resulted from mutually beneficial relationships among taxa, recent evidence suggests the opposite. Coyte et al. (2015) found that increasing the number of competitive interactions increased the stability of simulated microbial communities and the opposite for mutualistic interactions.      

Recent efforts to infer microbial interaction from time series data have further demonstrated that competition is indeed more common than mutualism (Stein et al. 2013; Marino et al 2014). Coyte et al. (2015) further demonstrated that their results for simulated microbial communities held true for the empirical result of Stein et al. (2013). Inferring interactions is challenging, however, and limitations in the methods meant that the communities in these studies were reduced from several hundred species to either 11 (Stein) or 17 (Marino). Each study focused on the most abundant taxa throughout the course of the time series and lumped the less abundant taxa into a single category ("Other").   


Typically our understanding of keystone species' impacts have come from species removal experiments. In this paper we describe an _in silico_ species removal experiment. We determine the impact that each species has on a simulated community and how that impact is mediated by the different ways species interact with one another.  


## Methods

Simulations began by generating a species pool whose interactions were known. Interactions among species were assigned using the Erdos-Renyi model for random directed networks with 200 nodes and a connectance of 0.2. Each species was assigned a self-damping term drawn from a beta distributions (a = 1.1, b = 5) scaled between -5 and 0. Off-diagonal interaction strengths were drawn from a normal distribution based on the interactions derived by Stein et al. (2013). This gave a normal distribution with a mean of -0.07 and standard deviation of 0.354. Growth rates for all species were positive and drawn from a uniform distribution between 0.1 and 1. 

Each individual community was created by sampling 50 species from the pool. All interactions defined above were assumed to be universal, so individual communities represented subsets of the initial species pool. Dynamics of each community were simulated using the Lotka-Volterra model of species interactions,

$$
\begin{equation}
\frac{dN_i}{dt} = rN_i + \sum(a_{ij}*N_j)   
\end{equation}
$$

where _r_ is the species specific growth rate, _aij_ is the effect of species _j_ on species _i_. The simulation was run for 1000 time steps, which was long enough for most communities to reach equilibrium (constant abundances) or steady-state (stable attractor). Note the difference is that at steady state the abundances of populations may be fluctuating, but they will remain that way unless perturbed. 

In order to identify which species were important to the stability of the community we systematically removed each species (one at a time) and measured the changes in the community. The starting point for species removals were the equilibrial/steady-state communities. Following the removal of a species, the resulting community dynamics were simulated for 1000 time steps using the same model and parameters as the initial community. At the end of each simulation, the impact on the community was measured using four metrics: (1) the mean change in abundance, (2) mean coefficient of variation in the first 50 time steps following removal, (3) persistence (the fraction of species with positive abundance), and (4) eigenvalue of the resulting Jacobian matrix with the largest real part.


The types and strengths of the interactions each species participated in were identified for every community. In addition we identified the structural roles of each species in every community. The structural properties we measured were the betweenness, closeness, eigenvector centrality, hub score, page rank, and the number of species within 2 links of the target species (neighborhood size). Betweenness is the number of shortest paths along which the target species lies. Closeness is the number of steps required to access every other species in the interaction network. The eigenvector centrality, hub score, and page rank are three methods to compute the importance of species in the network. 

A species was considered to be keystone based on its level of community importance (Power et al. 1996) with respece to four measures of impact on the community following species removal. Community importance is measured as

$$
\begin{equation}
CI_i = \frac{(t_N - t_D)}{t_N} * \frac{1}{p_i} 
\end{equation}
$$

where _t_N_ is the quantity of interest in the initial community, _t_D_ is the quantity of interest following the removal of the species _i_, and _p_i_ is the relative abundance of species _i_. Keystone species were those whose community importance was in the top 10th percentile for all four metrics.   

To determine what makes a species a keystone in the community, we used a generalized linear model to identify the effect of the measured species properties on keystone status. Keystoneness was modeled as a binomial variable. All combinations of predictor variables were assessed using the _dredge_ function of the __MuMIn__ R package, and ranked according to AIC. All models with deltaAIC < 2 were averaged together. In addition to defining keystone species as those in the top 10% of community importance for all four metrics (K_full), we created additional generalized linear models for keystone species defined as the top 10% in each community metric (K_persist, K_abund, K_eigen, K_var). To evaluate the predictive accuracy of the models we used k-fold cross-validation on the best model (lowest AICc) in each case. Cross-validation was done using the __DAAG__ package's _cv.binary_ function.   
  

## Results 

Equilibrium local communities ranged from 16 to 34 species (median = 25), with connectances between 0.15 and 0.3 (median = 0.215). 

The K_full model response variable inlcuded 42 instances of a keystone species out of 4537. In each individual metric model there were 454 instances. In the K_full averaged model all predictor variables were included. The strongest effects on keystoneness were the types and strengths of the target species' interactions. In particular a species was more likely to be identified as keystone if it participated in more competition links and fewer mutualistic links. The effect of the strengths of those interactions, however, is the opposite. Mutualism strength had a positive effect on keystoneness while competition strength had a negative effect. These results indicate that keystone species are those that compete weakly with many species but are also strongly mutualistic with few. The number of predator/parasitic links had a weak positive effect on keystoneness and the strength of those interactions had a large negative effect. Topological predictor variables had relatively small effects on keystoneness. However, eigenvalue centrality and Page Rank (two measures of a species topological importance) had a positive and negative effect on keystoneness respectively. 
     
By examining the effects of each predictor variable on each individual metric of keystoneness, we are able to parse out why each predictor may be influencing keystoneness in the full model. The averaged K_persist model included all predictor variables except Page Rank. The largest negative effects on keystoneness based on persistence were the strengths of competition and predation. The number of competitive links had a positive effect in this model as well, although that effect was smaller than in K_full. For mutualism both the number and strength of the interactions had a negative effect on keystoneness. Both closeness centrality and eigenvector centrality had positive effects on keystoneness in K_persist as well. This suggests that keystone species identified by persistence are generalist weak competitors. 

In the K_abund average model all predictor variables were included except the number of competitive links and the neighborhood size. The strongest positive effects on keystoneness defined as change in abundance were by the topological measures of closeness and Page Rank. Keystoneness was also positively related to the strength of mutualistic interactions. The number of mutualistic and predator/parasitic interactions and the strength of competition and predation/parasitism were negatively related to keystoneness. These results indicate that keystone species influencing abundance are those that are specialized mutualists that interact with species who interact with many species. 

All predictor variables were included in the averaged K_eigen model. The largest effects on keystoneness defined by having a positive effect on the eigenvalue of the Jacobian matrix were competition and predation/parasitism strength (negative), and eigenvalue centrality (positive). We found a weak positive effect of the number of competitive and predation/parasitic links as well. These results suggest that keystone species defined by their effect on local stability are those that compete weakly with many species, or with species that interact with many species. 

The K_var model, where keystone species were defined by the coefficient of variation in species abundances following removal, included all predictor variables. The largest positive effects on keystoneness were closeness centrality and Page Rank. Mutualistic interaction strength had a positive effect on keystoneness but the number of mutualistic interactions had a negative effect. Competitive and predation/parasitic interaction strength had negative effects as well. Thus species that have a large impact on the variation in abundance are those that are specialized mutualists interacting with species that interact with many species. 

![Figure 1](/Users/jjborrelli/Desktop/modelpar.jpeg "Figure 1")  

**Figure 1: Model averaged parameters with 95% confidence intervals for the four individual metric based impact models. Models with delta AICc < 2 were included in the averaging. Blue points indicate significantly different from 0. Note changes to the scale of the x-axis.**     


## Discussion

A species is more likely to be a keystone, important for the stability of the community, if it is a weak generalist competitor and a strong specialized mutualist. Each of these interaction types is allowing the prospective keystone species to be influencing the community in different ways. When it comes to the persistence and local stability of the community, removing species that have many weak competitive interactions has a large effect. Removing strong specialized mutualists leads to large impacts on equilibrium abundances and variability in the community.    

Our results corroborate the recent result of Coyte et al. (2015). Coyte and colleagues demonstrated that competition increases the stability of microbial networks, while mutualism decreases it. The number of competitive links a species participated in was positively related to keystoneness in all models, while the number of mutualistic links was negative in all models. Interaction strength, however, is also important to consider when assessing the stability of a community. While Coyte et al. (2015) focused only on the number of interactions, we found that the strength of competition is negatively related to keystoneness, but mutualistic interaction strength can positively influence keystoneness. Therefore it may not just be the number of competitive interactions that increase stability, but rather the number of weak competitive interactions.  

The position of a species in the network tended to have less of an effect, although it was dependent on the measure used. Metrics like closeness (how many steps it takes to get to every other node) and Page Rank (measure of importance) Had strong positive relationships with keystoneness as measured by change in abundance and initial variation. These effects disappeared, however, when keystoneness was defined by all impact metrics. 

In most macrobiological studies the species that have been identified as keystone are those that prey upon competitively dominant taxa. Here we have shown that when multiple interaction types and impact is broadly defined, keystone species are those that positively affect taxa that compete weakly with many others. In this context, having more weak and fewer strong links can, when defining impact as persistence or local stability, make a species more likely to be a keytsone. Species with strong predatory/parasitic links were less likely to be keystone species, despite what we may have expected from other experimentally determined keystones like _Pisaster_. 


It is clear that our definition of keystone species depends on how we define "large impact on the commmunity relative to their abundance." Paine (1969) suggested that the impact of a keystone species is on the stability of the community. Power et al. (1996) extended that definition to mean a change in any measurable community trait. Which community trait we pick is therefore going to determine what species are likely to be identified as keystones. We also need to consider what magnitude of impact relative to abundance makes for keystone status. In this paper we called any species whose community importance was in the top 10% a keystone. 

With our quantitative definition, not all communities will have a keystone species, and some may have multiple. Our strictest definition, combining all four impact metrics, only has 42 of the possible 4537 species-community combinations were assigned keystone status. For the single impact measures, 454 species-community combinations were considered keystones.  
   
These results will help us as we move forward with new methods to infer species interactions from time series data. With this new data we will now be able to identify putatively important species. These results can also guide future studies into probiotics. When attempting to introduce new species into the human gut for the purposes of manipulating the microbiome, we need to be able to identify those taxa that will have the largest impact on the community.   

Our study was nonetheless limited in its scope. The communities simulated here are orders of magnitude smaller thant those that exist in the gut. Furthermore, the species in these communities interacted through competition, mutualism, and predations/parasitism. There were no commensal or amensalistic links. 

## References


## Supplemental

![SuppFig 1](/Users/jjborrelli/Desktop/parimpt.jpeg "Supplemental Figure 1")

**Supplemental Figure 1: The importance (fraction of models containing that parameter) of each parameter in the averaged model for each of the four cases. All models with delta AICc < 2 were averaged. Blue bars indicate statistical significance of the parameter in the averaged model (p < 0.05).**