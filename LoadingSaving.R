
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

save.image("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Stats/Gradient/workspacecleandata.Rdata")  # 
load("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Stats/Gradient/workspacebioinformatics.Rdata")  # 


#my redo of the fungi bioinformatics, in qiime2 folder
Roots_SoilFungi2017EF2.RData

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


#Analysis workspaces old fungi/bac data
workspace1.Rdata
workspace2.Rdata

#Analysis new data
workspace3.Rdata
workspace4.Rdata





#List of files to export for Christina (workspaceCB.Rdata)
#rm(list=setdiff(ls(), c("spcomp4", "dat17","datBac3rc","richbac2","datITS2rcsoil","richITS2")))
spcomp4
dat17
datBac3rc
richbac2
datITS2rcsoil
richITS2


