#Phrag abundance, biomass, and litter analyses

options(contrasts=c("contr.treatment","contr.poly"));options("contrasts")
options(contrasts=c("contr.helmert","contr.poly"));options("contrasts")

#for glht/multcomp
model.matrix.gls <- function(object, ...) {
  model.matrix(terms(object), data = getData(object), ...)}
model.frame.gls <- function(object, ...) {
  model.frame(formula(object), data = getData(object), ...)}
terms.gls <- function(object, ...) {
  terms(model.frame(object), ...)}


#### Phrag abundance #####

dat17phrag<-dat17%>%
  filter(Transect%in%c("Transition","Phragmites"))
dat17phrag$Haplotype<-factor(ifelse(dat17phrag$Site%in%c("LUMCON 1","LUMCON 2"),"M1","I"))

#Final stats model for manuscript
#m1<-lme(phraus~MarshClassV*Transect2,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17phrag)
m2<-lme(phraus~MarshClassV*Transect2,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17phrag)#
anova(m1,m2) #var ident not significant
anova(m2,type="margin")
m1m<-as.data.frame(summary(emmeans(m2,~Transect2|MarshClassV)))

m1<-lme(phraus~MarshClassV.Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17phrag)
summary(glht(m1, linfct = mcp(MarshClassV.Transect = "Tukey")))

#only haplotype I
dat17phragI<-dat17phrag%>%
  filter(Site!="LUMCON 1",Site!="LUMCON 2")
dat17phragI$MarshClassV.Transect<-factor(dat17phragI$MarshClassV.Transect, levels=c("Fresh.Transition","Fresh.Phragmites","Brackish.Transition","Brackish.Phragmites"))
m3<-lme(phraus~MarshClassV*Transect2,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17phragI)#
#m4<-lme(phraus~MarshClassV*Transect2,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17phragI)#
anova(m3,m4)
anova(m3,type="margin")
m3<-lme(phraus~MarshClassV.Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17phragI)
summary(glht(m3, linfct = mcp(MarshClassV.Transect = "Tukey")))

#Final fig for manuscript
pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Manuscripts/Gradientms/Figs/phragabun.pdf",width=2.2,height=2.2)
ggplot(m1m,aes(x=Transect2,y=emmean,color=Transect2,group=MarshClassV))+
  #labs(x = "",y=expression(paste(italic("Phragmites")," density (stems"~m^-2~")"))) +
  labs(x = "",y="Phragmites stems m-2") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside")+
  geom_line(stat = "identity", position = "identity",size=.5,col="black")+
  geom_point(size=1.8)+
  geom_errorbar(aes(ymax = emmean+SE, ymin=emmean-SE),width=.25,size=.5) +
  scale_color_manual(values = c("gray50", "gray30"))+
  facet_wrap(vars(MarshClassV),strip.position = "bottom")
dev.off()


#old analyses/figs

dat17p<-dat17%>%
  #filter(Transect%in%c("Phragmites","Transition"))%>%
  dplyr::select(Line,Transect,Site,nat=NatAbun,AnnualSalinityppt,Phragmites="Phragmites australis")%>%
  group_by(Site,Transect,AnnualSalinityppt)%>%
  summarise(mean=mean(Phragmites),std=std.error(Phragmites,na.rm=T),sal=mean(AnnualSalinityppt))

phragabunfig<-ggplot(dat17p,aes(x=Transect,y=mean,color=Transect,group=Site))+
  labs(x = "",y="Phrag abundance")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5))+
  geom_line(stat = "identity", position = "identity",size=.8,col="black")+
  geom_point(size=3)+
  geom_errorbar(aes(ymax = mean+std, ymin=mean-std),width=.15,size=.8) +
  facet_wrap(vars(Site),nrow=1)#,scales="free"

##Trying a figure average per marshclassV, to do this correctly I would use least squares means via glht to incorporate error/random effects.
dat17p<-dat17%>%
  filter(Transect%in%c("Phragmites","Transition"))%>%
  dplyr::select(Line,Transect,Site,MarshClassV,nat=NatAbun,AnnualSalinityppt,Phragmites="Phragmites australis")%>%
  group_by(MarshClassV,Transect)%>%
  summarise(mean=mean(Phragmites),std=std.error(Phragmites,na.rm=T),sal=mean(AnnualSalinityppt))

ggplot(dat17p,aes(x=Transect,y=mean,color=Transect,group=MarshClassV))+
  labs(x = "",y="Phrag abundance")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5))+
  geom_line(stat = "identity", position = "identity",size=.8,col="black")+
  geom_point(size=3)+
  geom_errorbar(aes(ymax = mean+std, ymin=mean-std),width=.15,size=.8) +
  facet_wrap(vars(MarshClassV),nrow=1)#,scales="free"




#Stats

#There is zero spatial autocorrelation, interesting, maybe b/c of using the site effect?
#m1<-lme(phraus~MarshClassS*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17phrag)
#anova(m1,type="marginal")
m2<-lme(phraus~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17phrag)
anova(m2,type="marginal")
#m3<-lme(phraus~log(AnnualSalinityppt)*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17phrag)
#anova(m3,type="marginal")
#anova(m1,m2)

plot(log(dat17phrag$AnnualSalinityppt),dat17phrag$phraus,col=dat17phrag$Transect)


#Is site or haplotype more important random effect. Site is very slightly better
m1<-lme(phraus~MarshClassS*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17phrag)
m2<-lme(phraus~MarshClassS*Transect,random=~1|Haplotype,correlation=corSpher(form = ~ Lat+Long),data=dat17phrag)
anova(m1,m2)

#Can I put haplotype in fixed effects? NO
m2<-lme(phraus~MarshClassV*Transect+Haplotype,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17phrag)
anova(m2,type="marginal")

#Can I put both haplotype and site as random effect? NO
library(lme4)
m1<-lmer(phraus~MarshClassV*Transect+(1|Site)+(1|Haplotype),data=dat17phrag)
summary(m1)


#Doing analysis only on haplotype I
dat17phragI<-dat17phrag%>%
  filter(Haplotype=="I")

m1<-lme(phraus~MarshClassS*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17phragI)
anova(m1,type="marginal")
m2<-lme(phraus~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17phragI)
anova(m2,type="marginal")
m3<-lme(phraus~log(AnnualSalinityppt)*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17phragI)
anova(m3,type="marginal")
anova(m1,m2)


#Getting least squares means
dat17phrag
dat17phrag$MarshClassV.Transect<-interaction(dat17phrag$MarshClassV,dat17phrag$Transect)
dat17phrag$MarshClassV.Transect<-factor(dat17phrag$MarshClassV.Transect,levels=c("Fresh.Transition","Fresh.Phragmites","Brackish.Transition","Brackish.Phragmites","Saline.Transition","Saline.Phragmites"))

m1<-lme(phraus~MarshClassV.Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17phrag)
anova(m1)
summary(glht(m1, linfct = mcp(MarshClassV.Transect = "Tukey")))
emmeans(m1,~MarshClassV.Transect)






##### Biomass and litter #####

dat17$biomasskg<-dat17$biomass/1000
dat17$litterkg<-dat17$litter/1000

#Live biomass
m1<-lme(biomasskg~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17,na.action=na.omit,control = lmeControl(maxIter=100,msMaxIter = 100))
#m2<-lme(biomasskg~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17,na.action=na.omit)
anova(m1,m2) #het var is sig
anova(m1,type="marginal")

m1biomass<-as.data.frame(summary(emmeans(m1,~Transect|MarshClassV)))

m1<-lme(biomasskg~MarshClassVTransect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17,na.action=na.omit)
summary(glht(m1, linfct = mcp(MarshClassVTransect = "Tukey")))

#Only on haplotype I
m1<-lme(biomasskg~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17I,na.action=na.omit,control = lmeControl(maxIter=100,msMaxIter = 100))
m2<-lme(biomasskg~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17I,na.action=na.omit)
anova(m1,m2) #het var is sig
anova(m1,type="marginal")



#Litter
m1<-lme(litterkg~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17,na.action=na.omit,control = lmeControl(maxIter=100,msMaxIter = 100))
#m2<-lme(litterkg~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17,na.action=na.omit)
anova(m1,m2) #het var very significant
anova(m1,type="marginal")

m1litter<-as.data.frame(summary(emmeans(m1,~Transect|MarshClassV)))

#m1<-lme(litterkg~MarshClassVTransect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17,na.action=na.omit)
#summary(glht(m1, linfct = mcp(MarshClassVTransect = "Tukey")))

#multiple comparisons within marsh type
m2<-lme(litterkg~MarshClassV.Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17,na.action=na.omit,control = lmeControl(maxIter=100,msMaxIter = 100))
#m2<-lme(litterkg~MarshClassV.Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17,na.action=na.omit)
summary(glht(m2, linfct = mcp(MarshClassV.Transect=c("Fresh.Native-Fresh.Transition=0","Fresh.Native-Fresh.Phragmites=0","Fresh.Transition-Fresh.Phragmites=0","Brackish.Native-Brackish.Transition=0","Brackish.Native-Brackish.Phragmites=0","Brackish.Transition-Brackish.Phragmites=0","Saline.Native-Saline.Transition=0","Saline.Native-Saline.Phragmites=0","Saline.Transition-Saline.Phragmites=0"))))

#Note: when doing Tukey contrasts, if you have interaction terms you need to use helmert contrasts. it should be fine to have interactions as long as they are not significant.

#Only on haplotype I
m1<-lme(litterkg~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17I,na.action=na.omit,control = lmeControl(maxIter=100,msMaxIter = 100))
m2<-lme(litterkg~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17I,na.action=na.omit)
anova(m1,m2) #het var is sig
anova(m1,type="marginal")


#Figures by marsh class

#Biomass (live)

pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Manuscripts/Gradientms/Figs/biomasshetvar.pdf",width=2.2,height=2.2)#2.1
ggplot(m1biomass,aes(x=Transect,y=emmean,color=Transect,group=MarshClassV))+
  labs(x = "",y="Live biomass (kg m-2)")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside")+
  geom_point(size=1.8)+
  geom_line(stat = "identity", position = "identity",size=.5,color="black")+
  geom_errorbar(aes(ymax = emmean+SE, ymin=emmean-SE),width=.25,size=.5) +
  scale_color_manual(values = c("gray70","gray50", "gray30"))+
  facet_wrap(vars(MarshClassV),nrow=1,strip.position = "bottom")#,scales="free"
dev.off()
#,axis.text.x = element_text(angle = -45)


#Litter
pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Manuscripts/Gradientms/Figs/litterhtvar.pdf",width=2.2,height=2.2)#2.1
ggplot(m1litter,aes(x=Transect,y=emmean,color=Transect,group=MarshClassV))+
  labs(x = "",y="Litter mass (kg m-2)")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside")+
  geom_point(size=1.8)+
  geom_line(stat = "identity", position = "identity",size=.5,color="black")+
  geom_errorbar(aes(ymax = emmean+SE, ymin=emmean-SE),width=.25,size=.5) +
  scale_color_manual(values = c("gray70","gray50", "gray30"))+
  facet_wrap(vars(MarshClassV),nrow=1,strip.position = "bottom")#,scales="free"
dev.off()



#Figures by site (old)
#Biomass (live)
biomass1<-dat17%>%
  group_by(Site, Transect)%>%
  summarise(mean=mean(Mass,na.rm=T),std=std.error(Mass,na.rm=T))

livefig<-ggplot(biomass4,aes(x=Transect,y=mean,color=Transect,group=Site))+
  labs(x = "",y="Live mass (g/m2)")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5))+
  geom_point(size=3)+
  geom_line(stat = "identity", position = "identity",size=.8,color="black")+
  geom_errorbar(aes(ymax = mean+std, ymin=mean-std),width=.15,size=.8) +
  facet_wrap(vars(Site),nrow=1)#,scales="free"

#litter
dat17%>%
  group_by(MarshClassV, Transect)%>%
  summarise(mean=mean(litterkg,na.rm=T),se=std.error(litterkg,na.rm=T))

biomass2<-dat17%>%
  #filter(litterkg<6.6)%>%
  group_by(Site, Transect)%>%
  summarise(mean=mean(litterkg,na.rm=T),se=std.error(litterkg,na.rm=T))

#litterfig<-
  ggplot(biomass2,aes(x=Transect,y=mean,color=Transect,group=Site))+
  labs(x = "",y="Litter mass (g/m2)")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5))+
  #scale_color_manual(values=c("#C0C0C0", "#808080", "#000000"))+ #for black and white
  geom_point(size=3)+
  geom_line(stat = "identity", position = "identity",size=.8,color="black")+
  geom_errorbar(aes(ymax = mean+se, ymin=mean-se),width=.15,size=.8) +
  facet_wrap(vars(Site),nrow=1)#,scales="free"


pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Stats/Gradient/Figs/livelitter.pdf")
plot_grid(livefig, litterfig,nrow = 2)
dev.off()








##### Native abundance #####

#Figure wth least squares means and grouped by marsh class

#Final stats model for manuscript
m1<-lme(NatAbun~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassVTransect),data=dat17)
#m2<-lme(NatAbun~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17) #,weights=varIdent(form=~1|MarshClassVTransect), i tried heterogeneous errors using varident to see if abundance in phrag plot at lumcon which only varied from 0-2 native stems would ever be different from fresh or brackish. is it not. i was confused about this, but I thnk it is still correct b/c the ls means are taking into accout that we only have 2 saline sites so we dont know much about those areas.
anova(m1,m2) #het var significant
anova(m1,type="margin")
m1m<-as.data.frame(summary(emmeans(m1,~Transect|MarshClassV)))

m1<-lme(NatAbun~MarshClassV.Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17)
#m1<-lme(NatAbun~MarshClassV.Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17)
#summary(glht(m1, linfct = mcp(MarshClassV.Transect = "Tukey")))
summary(glht(m1, linfct = mcp(MarshClassV.Transect=c("Fresh.Native-Fresh.Transition=0","Fresh.Native-Fresh.Phragmites=0","Fresh.Transition-Fresh.Phragmites=0","Brackish.Native-Brackish.Transition=0","Brackish.Native-Brackish.Phragmites=0","Brackish.Transition-Brackish.Phragmites=0","Saline.Native-Saline.Transition=0","Saline.Native-Saline.Phragmites=0","Saline.Transition-Saline.Phragmites=0"))))

#how to do least squares means using glht()
#summary(glht(m1, linfct = lsm(~ MarshClassV.Transect)))

#Only on haplotype I
m1<-lme(NatAbun~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17I,na.action=na.omit,control = lmeControl(maxIter=100,msMaxIter = 100))
m2<-lme(NatAbun~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17I,na.action=na.omit)
anova(m1,m2) #het var is sig
anova(m1,type="marginal")


#Final fig for manuscript
pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Manuscripts/Gradientms/Figs/natabunhetvar.pdf",width=2.2,height=2.2)
ggplot(m1m,aes(x=Transect,y=emmean,color=Transect,group=MarshClassV))+
  labs(x = "",y="Native stems m-2") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside")+
  geom_line(stat = "identity", position = "identity",size=.5,col="black")+
  geom_point(size=1.8)+
  geom_errorbar(aes(ymax = emmean+SE, ymin=emmean-SE),width=.25,size=.5) +
  scale_color_manual(values = c("gray70","gray50", "gray30"))+
  facet_wrap(vars(MarshClassV),strip.position = "bottom")
dev.off()


#testing adding litter or phrag abun to the model
m2<-lme(NatAbun~MarshClassV*phraus,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17,na.action = na.omit)
anova(m2,type="margin")
m2<-lme(NatAbun~MarshClassV+phraus+litterkg,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17,na.action = na.omit)
anova(m2,type="margin")
plot(dat17$phraus,dat17$NatAbun,col=dat17$MarshClassV)


# Old stats and figs by site
dat17b <- dat17%>%
  dplyr::select(MarshClassV,Transect,Site,NatAbun)%>%
  group_by(MarshClassV,Transect,)%>%
  summarise(mean=mean(NatAbun),se=std.error(NatAbun,na.rm=T))#,sal=mean(AnnualSalinityppt)

pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Manuscripts/Gradientms/Figs/phragabun.pdf",width=2.2,height=2.2)
ggplot(dat17b,aes(x=Transect,y=mean,color=Transect,group=MarshClassV))+
  labs(x = "",y="Phragmites stems m-2") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside")+
  geom_line(stat = "identity", position = "identity",size=.5,col="black")+
  geom_point(size=1.8)+
  geom_errorbar(aes(ymax = mean+se, ymin=mean-se),width=.25,size=.5) +
  facet_wrap(vars(MarshClassV),strip.position = "bottom")
dev.off()

pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Stats/Gradient/Figs/phragnatabun.pdf")
plot_grid(phragabunfig, natabunfig,nrow = 2)
dev.off()



dat17%>%
  filter(Site=="LUMCON 1"|Site=="LUMCON 2")%>%
  dplyr::select(Site,Transect, NatAbun)

dat17%>%
  filter(MarshClassV=="Brackish")%>%
  filter(Transect=="Phragmites")%>%
  dplyr::select(Site,Transect, NatAbun)

dat17%>%
  filter(MarshClassV=="Fresh")%>%
  filter(Transect=="Phragmites")%>%
  dplyr::select(Site,Transect, NatAbun)

















###### LRR ######

#In each site take the mean nat abun of each transect. 
dat2 <- dat%>%
  filter(Year==2017)%>%
  group_by(Transect,Site)%>%
  summarise(nat=mean(NatRichness,na.rm=T),sal=mean(MeanSalinityppt))%>%
  spread(Transect,nat)%>%
  group_by(Site)%>%
  summarise(nat=mean(Native,na.rm=T),phrag=mean(Phragmites,na.rm=T),sal=mean(sal,na.rm=T))%>%
  mutate(lrr=log((phrag)/(nat)))
#  mutate(lrr=log((phrag+1)/(nat+1)))

plot(dat2$sal,dat2$lrr)

dat2 <- dat%>%
  filter(Year==2017)%>%
  dplyr::select(Line,Transect,Site,nat=NatAbun,MeanSalinityppt)%>%#NatRichness
  spread(Transect,nat)%>%
  mutate(lrr=log((Phragmites+1)/(Native+1)))%>%
  group_by(Site,MeanSalinityppt)%>%
  summarise(stdlrr=std.error(lrr),lrr=mean(lrr,na.rm=T),sal=mean(MeanSalinityppt))
dat2

dat2$Site<-factor(dat2$Site,levels=c("Barataria","Turtle Cove","Pearl River","Fontainebleau","Big Branch","Bayou Sauvage","LUMCON 2","LUMCON 1"))
#plot(dat2$sal,dat2$lrr)

ggplot(dat2,aes(x=Site,y=lrr))+
  labs(x = "",y="lrr abundance")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5))+
  #geom_line(stat = "identity", position = "identity",size=.8)+
  geom_point(size=3)+
  geom_errorbar(aes(ymax = lrr+stdlrr, ymin=lrr-stdlrr),width=.15,size=.8) 
#facet_wrap(~type,scales="free")

