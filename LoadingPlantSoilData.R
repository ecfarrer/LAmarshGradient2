#Loading the salinity and vegetation data


##### CRMS data #####
# CRMS daily salinity data (hydrology) downloaded from https://lacoast.gov/chart/Charting.aspx?laf=crms&tab=2 (the crms site not cims). 
# Monthly hydrographic data (soil porewater) downloaded from cims  https://cims.coastal.louisiana.gov/DataDownload/DataDownload.aspx?type=hydro_monthly. this is what Christina/Pawel used for their ms. you need to go in and remove degree symbols from the column names for the file to read properly.

#H data are hydrology - Water salinity, temperature, and water level (in the water/canal). there is 1 H plot per station
#S are Soil- Bulk Density, % Organic Matter, Salinity, pH, & Moisture. 6 S plots per station
#P are Soil Porewater, Salinity, and Temperature. There are 3 P plots per station
#A are accretion. 3 A plots
#V are veg, 10 veg plots per station
#M hydrography - Water Level, Water Temperature. 1 M plot


##### Water #####
#I did not use this data 

#sal<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Data/CRMS/Hydrology/8148_CalendarYear/8148.csv")

#head(sal)

#Barataria CRMS0188-H01
#Pearl River	CRMS4110-H01
#Turtle cove	CRMS0030-H01
#Fontainblebleau	CRMS2854-H01
#Big Branch	CRMS0006-H01
#Bayou Sauvage	CRMS4107-H01
#LUMCON1	CRMS0311-H01
#LUMCON2	CRMS0311-H01


#sal2<-sal%>%
#  filter(Station_id=="CRMS0188-H01")
#sal2



###### Soil #####
#I will use this

sals<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Data/CRMS/Porewater/HYDROGRAPHIC_MONTHLY.csv")

#note that the latest measurement is 10/2019. Should download updated data in a few months to get nov and dec 2019. For this analysis I used 2010-2017. I could go earlier if I wanted, but Turtle Cove increases.
sals2<-sals%>%
  separate(CPRA.Station.ID,c("StationFront","ID"),"-",remove=F)%>%
  separate(ID,c("StationType",NA),1,remove=F)%>%
  filter(StationFront%in%c("CRMS0188","CRMS0030","CRMS4110","CRMS2854","CRMS0006","CRMS4107","CRMS0311"))%>%
  filter(Measurement.Depth..ft.==.328)%>%
  dplyr::select(Station=CPRA.Station.ID,StationFront,ID,StationType,Date=Date..mm.dd.yyyy.,Salinityppt=Soil.Porewater.Salinity..ppt.)%>%
  separate(Date,c("Month","Day","Year"),"/",remove=F)%>%
  filter(Year%in%c("10","11","12","13","14","15","16","17"),StationType=="P")%>%#originally I used 2010-19
  group_by(StationFront)%>%
  summarise(AnnualSalinityppt=mean(Salinityppt, na.rm=T))
head(sals2)
sals2$Site<-c("Big Branch","Turtle Cove","Barataria","LUMCON","Fontainebleau","Bayou Sauvage","Pearl River")
#sals2[order(sals2$MeanSalinityppt),]

sals2$Site[sals2$Site=="LUMCON"]<-"LUMCON 1"
sals2[8,]<-sals2[4,]
sals2[8,3]<-"LUMCON 2"
sals2

sals2[order(sals2$AnnualSalinityppt),]



#### GPS coordinates####
gps<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Data/2017/Phragmites_plots_ and_site_crds_2017.csv")

gps2<-gps%>%
  dplyr::select(Plot,Line,Transect,Lat,Long)


  

##### Vegetation #####
spcomp<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Data/LAmarsh_Survey_2017to2019_SppComp_cleaned.csv")

spcomp2<-spcomp%>%
  dplyr::mutate(Species=dplyr::recode(Species,"Phragmites australis (L)"="Phragmites australis"))%>%
  filter(Species!="Phragmites australis (D)")
spcomp2$Species<-factor(spcomp2$Species)
sort(levels(spcomp2$Species))

spcomp3<-spcomp2%>%
  spread(Species,StemCount,fill=0)%>%
  arrange(Year,Plot)
spcomp4<-cbind(spcomp3[,1:4],"Phragmites australis"=spcomp3[,"Phragmites australis"],spcomp3[,5:28],spcomp3[30:51])
spcomp4$Richness<-rowSums(spcomp4[,5:51]>0)
spcomp4$Shannon<-diversity(spcomp4[,5:51], index="shannon")
spcomp4$Evenness<-spcomp4$Shannon/log(spcomp4$Richness)
spcomp4$NatAbun<-rowSums(spcomp4[,6:51])
spcomp4$NatRichness<-rowSums(spcomp4[,6:51]>0)
spcomp4$NatShannon<-diversity(spcomp4[,6:51], index="shannon")
spcomp4$NatEvenness<-spcomp4$NatShannon/log(spcomp4$NatRichness)
spcomp4$Simpson<-diversity(spcomp4[,5:51], index="simpson")
spcomp4$NatSimpson<-diversity(spcomp4[,6:51], index="simpson")
spcomp4$InvSimpson<-diversity(spcomp4[,5:51], index="invsimpson")
spcomp4$NatInvSimpson<-diversity(spcomp4[,6:51], index="invsimpson")




##### Biomass litter ######

biomass<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Data/2017/LAmarsh_Survey_2017_biomasslitter.csv")

biomass2 <- biomass%>%
  filter(is.na(Mass.g.20x20cm.area)==T|Mass.g.20x20cm.area<400)%>% #there are two plots in BS one litter, one biomass with super high values, they are possibly correct but they are just super huge so I'm removing them
  left_join(gps2)%>%
  dplyr::select(Plot,Site,Line,Transect,L.B,Mass=Mass.g.20x20cm.area)%>%
  mutate(Mass=Mass*25)
biomass2$Site<-factor(biomass2$Site,levels=c("Barataria","Turtle Cove","Pearl River","Fontainebleau","Big Branch","Bayou Sauvage","LUMCON 2","LUMCON 1"))
biomass2$Transect=factor(biomass2$Transect,levels=c("Native","Transition","Phragmites"))

biomass3<-biomass2%>%
  spread(L.B,Mass,fill=NA)




###### Field salinity measurements ######
fieldsal<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Data/2017/LAmarsh_Survey_2017_WaterChemistry.csv")
fieldsal2<-fieldsal%>%
  dplyr::select(Plot,Salinity15cmppt,Salinity45cmppt, WaterDepthcm)

aggregate.data.frame(fieldsal$Salinity15cmppt,by=list(fieldsal$Site),mean,na.rm=T)
aggregate.data.frame(fieldsal$Salinity15cmppt,by=list(fieldsal$Site,fieldsal$Transect),mean,na.rm=T)



###### Getting one data file together ######

dat17<-spcomp4%>%
  filter(Year==2017)%>%
  full_join(sals2)%>%
  full_join(gps2)%>%
  full_join(biomass3)%>%
  full_join(fieldsal2)%>%
  arrange(Plot)

dat17$Site<-factor(dat17$Site,levels=c("Barataria","Turtle Cove","Pearl River","Fontainebleau","Big Branch","Bayou Sauvage","LUMCON 2","LUMCON 1"))
dat17$Transect<-factor(dat17$Transect,levels=c("Native","Transition","Phragmites"))

#Make an additional phragmites column with no space
dat17$phraus<-dat17$`Phragmites australis`

#Salinity classes based on salinity levels
dat17$MarshClassS<-ifelse(dat17$Site%in%c("Barataria","Turtle Cove","Pearl River"),"Fresh",ifelse(dat17$Site%in%c("Fontainebleau","Big Branch","Bayou Sauvage"),"Brackish","Saline"))
dat17$MarshClassS<-factor(dat17$MarshClassS,levels=c("Fresh","Brackish","Saline"))

#Salinity classes based on vegetation composition
dat17$MarshClassV<-ifelse(dat17$Site%in%c("Barataria","Turtle Cove"),"Fresh",ifelse(dat17$Site%in%c("Pearl River","Fontainebleau","Big Branch","Bayou Sauvage"),"Brackish","Saline"))
dat17$MarshClassV<-factor(dat17$MarshClassV,levels=c("Fresh","Brackish","Saline"))

#One variable with levels for marshclassv*transect
dat17$MarshClassVTransect<-factor(paste(dat17$MarshClassV,dat17$Transect,sep=" "),levels=c("Fresh Native","Fresh Transition","Fresh Phragmites","Brackish Native","Brackish Transition","Brackish Phragmites","Saline Native","Saline Transition","Saline Phragmites"))

dat17$MarshClass_Transect<-plyr::revalue(dat17$MarshClassVTransect,c("Fresh Phragmites"="Fresh Monoculture","Brackish Phragmites" = "Brackish Monoculture","Saline Phragmites"="Saline Monoculture"))

dat17$MarshClass_Transect<-factor(dat17$MarshClass_Transect,levels=c("Fresh Native","Fresh Transition","Fresh Monoculture","Brackish Native","Brackish Transition","Brackish Monoculture","Saline Native","Saline Transition","Saline Monoculture"))

dat17$Transect2<-plyr::revalue(dat17$Transect,c("Phragmites"="Monoculture"))

dat17$MarshClassV.Transect<-interaction(dat17$MarshClassV,dat17$Transect)
dat17$MarshClassV.Transect<-factor(dat17$MarshClassV.Transect,levels=c("Fresh.Native","Fresh.Transition","Fresh.Phragmites","Brackish.Native","Brackish.Transition","Brackish.Phragmites","Saline.Native","Saline.Transition","Saline.Phragmites"))





###### run above, then add richess info from fungi and bacteria and add to dat17, then run below to mke supplemntary table #####

##### Soil and environemtnal data #####
soilenv<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Data/2017/Combined_soil_results.csv")
ph<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Data/2017/Soil_pH_PhragmitesProject_2018.csv")

ph2<-ph%>%
  dplyr::select(Plot,pH)
ph2$Plot<-as.numeric(ph2$Plot)

soilenv2<-soilenv%>%
  separate(SampleID,into=c(NA,"Plot"),sep=2)%>%
  mutate(Plot=as.numeric(Plot))%>%
  dplyr::select(-Site,-SiteID,-Transect)%>%
  full_join(ph2)

dat17b<-dat17%>%
  full_join(soilenv2)
dat17<-dat17b
  



#only haplotype I (run after doing fungal/bacterial bioinformatics)
dat17I<-dat17%>%
  filter(Site!="LUMCON 1",Site!="LUMCON 2")
dat17aI<-dat17a%>%
  filter(Site!="LUMCON 1",Site!="LUMCON 2")

#errors in nat/trans/phrag designation in the sp comp files
#2019
#37 should be native
#58 is native
#79 native
#155 phrag
#2018
#35 36 should be transition






