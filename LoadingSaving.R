
library(phyloseq)
library(dada2)
library(tidyverse)
library(vegan)
library(plotrix)
library(fossil) #for sorenson()
library(nlme)
library(cowplot)
library(car)
library(MASS)
library(pscl)
library(multcomp)
library(Hmisc)
library(deming)#for deming regression with x and y error
library(emmeans)#for least squares means, new
#library(lsmeans) #for least squares means, old
#library(ggforce)#for some ellipses on ordintions

options(contrasts=c("contr.helmert","contr.poly"));options("contrasts")

save.image("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Stats/Gradient/workspace5.Rdata")  # 
load("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Stats/Gradient/workspace5.Rdata")  # 


#Bioinformatics on fungi and bacteria
workspacebioinformatics.Rdata

#After bioinformatics only necessary dataframes
workspacecleandata.Rdata
#rm(list=setdiff(ls(), c("datITSS4","datITSS5", "datITSS5c","datITSS5otu","datITSS5cotu","datBac3","datBac4","datBac4c","datBac4otu","datBac4cotu","dat17")))
datITSS4
datITSS5
datITSS5c
datITSS5otu
datITSS5cotu
datBac3
datBac4
datBac4c
datBac4otu
datBac4cotu
dat17 #has microbial richness added to it


#Analysis workspaces 
workspace5.Rdata





#List of files to export for Christina (workspaceCB.Rdata)
#rm(list=setdiff(ls(), c("spcomp4", "dat17","datBac3rc","richbac2","datITS2rcsoil","richITS2")))
spcomp4
dat17
datBac3rc
richbac2
datITS2rcsoil
richITS2


