###############################################################
## Title: Impute missing trait data using Rphylopars and MICE 
## User: Lisbeth Hordley 
## email: l.hordley@pgr.reading.ac.uk
## Date: October 2019
##############################################################

rm(list=ls()) # clear R

library(ape)
library(phangorn)
library(Rphylip)
library(phytools)
library(Rphylopars)
library(mice)
library(evobiR)
library(factoextra)
library(PVR)
library(MPSEM)

## read in the subset of species trees
bird_subset <- read.nexus("../Data/Phylogeny/Bird_subset.nex") ## 1000 trees (Hackett family-level backbone)
bird_traits <- read.csv("../Data/Trait data/ALL_trait_data.csv", header=TRUE)

bird_phylo <- maxCladeCred(bird_subset)
plotTree(bird_phylo,ftype="i",fsize=0.8,lwd=1) ## 81 species!

## match order of species in phylo tree with trait dataframe
correct_order <- bird_phylo$tip.label
correct_order <- as.data.frame(correct_order)
correct_order <- gsub("_", " ", correct_order$correct_order)
bird_traits <- bird_traits[ match(correct_order, bird_traits$Scientific_name), ]
bird_phylo$tip.label = gsub("_", " ", bird_phylo$tip.label)

## check each trait with missing data for phylogenetic signal (9 traits with missing data)
gapewidth_K <- phylosig(bird_phylo, bird_traits$Gape_width, method="K", test=TRUE) ## K=1.3 - STRONG phylogenetic signal
latitude_K <- phylosig(bird_phylo, bird_traits$Latitude_mean, method="K", test=TRUE) ## K=0.64 - weak phylogenetic signal
SSI_K <- phylosig(bird_phylo, bird_traits$SSI, method="K", test=TRUE) ## K=0.26 - weak phylogenetic signal
maxlong_K <- phylosig(bird_phylo, bird_traits$maximum_longevity_y, method="K", test=TRUE) ## K=0.49 - weak phylogenetic signal
thermmax_K <- phylosig(bird_phylo, bird_traits$thermal_max, method="K", test=TRUE) ## K=0.27 - weak phylogenetic signal

###########################################################################################################
## First try RPhylopars to impute data

## get trait dataset in right format
bird_traits <- bird_traits[,-c(1,3,5,22)] ## remove common name, CBC code and duplicated clutch size and lifespan traits
colnames(bird_traits)[1] <- "species" ## change name of first column to species
bird_traits <- bird_traits[,c(1:2,13:20)] ## keep the 9 traits with missing data

## impute traits using Rphylopars
bird_imp_pylo <- phylopars(trait_data = bird_traits, tree = bird_phylo) 
bird_imp_pylo$anc_recon ## imputed species means

bird_imp1 <- as.data.frame(bird_imp_pylo$anc_recon)
## remove the random extra 80 rows of species
bird_imp1 <- bird_imp1[-c(82:161),]
bird_imp1$species <- row.names(bird_imp1)
row.names(bird_imp1) <- 1:nrow(bird_imp1)
bird_imp1 <- bird_imp1[,c(10,1,2,3,4,5,6,7,8,9)]

write.tree(bird_phylo, file="../Data/Phylogeny/bird_max_clade_cred_tree.tre")
write.csv(bird_imp1, file="../Data/Trait data/bird_traits_imp_phylopars.csv", row.names=FALSE)

###########################################################################################################
## try MICE to impute data

## calculate pairwise distances from the phylogenetic tree
bird_pvr <- PVRdecomp(bird_phylo)
bird_traits <- bird_traits[,2:10]
row.names(bird_traits) <- correct_order
bird_traits <- as.matrix(bird_traits)
trait1 <- bird_traits$Body_Mass
bird_eigenv <- PVR(bird_pvr, bird_phylo, bird_traits, method="moran")


TreeGraph = Phylo2DirectedGraph(bird_phylo)

#Option 1: Try to fit a model to estimate steepness in the tree - cant get this working for your data for some reason
BasePEM = PEM.fitSimple(
  y = bird_traits$Gape_width,
  x = NULL,
  w = TreeGraph,
  d = "distance", sp="species")
BasePEM_DF = as.data.frame(BasePEM)

#Option 2: Define the steepness e.g. 0.2. Read the vignette to decide which option
BasePEM = PEM.build(TreeGraph, a = 1)

BasePEM_LM = lmforwardsequentialsidak(
  y = bird_traits$Gape_width,
  object = BasePEM, 
  alpha = 0.2)

BaseEigens = as.data.frame(BasePEM_LM$qr$qr)
TreeTips = as.character(rownames(BaseEigens))
BaseEigens = cbind(TreeTips, BaseEigens) #remove intercept









bird_imp <- mice(bird_traits[,2:10], m=20, maxit=100, meth="rf")
## m = 20 (number of imputed datasets)
## maxit = 100 (number of iterations taken to impute missing values)
## method = predictive mean matching
completedData <- complete(bird_imp,1) ## save the first imputed dataset
completedData$species <- correct_order ## put in species
gapewidth <- hist(bird_traits$Gape_width) ## distribution of trait before imputation
gapewidth2 <- hist(completedData$Gape_width) ## distribution of trait before imputation
densityplot(bird_imp)
bird_imp$imp$Latitude_mean
bird_imp$imp$SSI
bird_imp$imp$clutch_size
