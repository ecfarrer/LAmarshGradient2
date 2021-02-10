# exporting for Dryad

#plot info, biomass, diversity measures

which(colnames(dat17)=="Phragmites australis")
which(colnames(dat17)=="Vigna luteola")

PlotData<-dat17[,-c(5:51)]
PlotData<-PlotData%>%
  dplyr::select(-phraus,-MarshClassS,-MarshClassVTransect,-MarshClass_Transect,-Transect2,-MarshClassV.Transect,-Calciumppm,-Copperppm,-Magnesiumppm,-Phosphorusppm,-Potassiumppm,-Sodiumppm,-Sulfurppm,-Zincppm,-Pathotroph,-Saprotroph,-Symbiotroph,-Endophyte,-Arbuscular.Mycorrhizalhighlyprobable,-Endophytehighlyprobable,-Plant.Pathogenhighlyprobable,-Ectomycorrhizal,-PlantMutualist,-Plant.Pathogenanywhere,-Ectomycorrhizal.nooutliers,-Salinity15cmppt,-Salinity45cmppt)%>%
  rename(PercentCarbon=Carbon.,PercentNitrogen=Nitrogen.)
PlotData<-cbind(PlotData[,1:31],RichnessBac.nooutlier=PlotData$RichnessBac,PlotData[,32:42])
PlotData<-cbind(PlotData[,1:33],Chao1Bac.nooutlier=PlotData$Chao1Bac,PlotData[,34:43])

PlotData$RichnessBac.nooutlier[167]<-NA
PlotData$Chao1Bac.nooutlier[167]<-NA

head(PlotData)

write.csv(PlotData, "/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Manuscripts/Gradientms/Dryad/PlotData.csv",row.names=F)



#plant composition
#from dat17

PlantComp<-cbind(Plot=dat17$Plot,dat17[,5:51])
head(PlantComp)
write.csv(PlantComp, "/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Manuscripts/Gradientms/Dryad/PlantComp.csv",row.names=F)


#fungi otu table, need to merge with taxonomy
#datITSS5cotu[1:5,1:5]
taxtempf<-tax_table(datITSS5c)
otutempf<-t(otu_table(datITSS5c))
FungiOTUtable<-cbind(otutempf,taxtempf)
head(FungiOTUtable)
write.csv(FungiOTUtable, "/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Manuscripts/Gradientms/Dryad/FungiASVtable.csv",row.names=T)


#bacteria otu table, need to merge with taxonomy
#datBac4cotu
taxtempb<-tax_table(datBac4c)
otutempb<-otu_table(datBac4c)
BacteriaOTUtable<-cbind(otutempb,taxtempb)
head(BacteriaOTUtable)
dim(BacteriaOTUtable)
colnames(BacteriaOTUtable)[162]<-"Kingdom"
colnames(BacteriaOTUtable)[163]<-"Phylum"
colnames(BacteriaOTUtable)[164]<-"Class"
colnames(BacteriaOTUtable)[165]<-"Order"
colnames(BacteriaOTUtable)[166]<-"Family"
colnames(BacteriaOTUtable)[167]<-"Genus"
colnames(BacteriaOTUtable)[168]<-"Species"
write.csv(BacteriaOTUtable, "/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Manuscripts/Gradientms/Dryad/BacteriaASVtable.csv",row.names=T)



