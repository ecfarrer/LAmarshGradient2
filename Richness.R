# Richness analyses

options(contrasts=c("contr.treatment","contr.poly"));options("contrasts")
options(contrasts=c("contr.helmert","contr.poly"));options("contrasts")


#### Plants ####

###### Plant Richness #####

#Final stats model for manuscript. adding heteroscedasticity is significant and helps the tukey results

m1<-lme(NatRichness~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17)
#m2<-lme(NatRichness~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17)
anova(m1,m2) #het var sig
anova(m1,type="marginal")
summary(glht(m1,linfct=mcp(Transect="Tukey")))
summary(glht(m1,linfct=mcp(MarshClassV="Tukey")))
hist(resid(m1,type="normalized"))
plot(fitted(m1),resid(m1,type="normalized"))

m1P<-as.data.frame(summary(emmeans(m1,~Transect|MarshClassV)))

m2<-lme(NatRichness~MarshClassV.Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17)#
#summary(glht(m1, linfct = mcp(Transect = "Tukey")))
summary(glht(m2, linfct = mcp(MarshClassV.Transect=c("Fresh.Native-Fresh.Transition=0","Fresh.Native-Fresh.Phragmites=0","Fresh.Transition-Fresh.Phragmites=0","Brackish.Native-Brackish.Transition=0","Brackish.Native-Brackish.Phragmites=0","Brackish.Transition-Brackish.Phragmites=0","Saline.Native-Saline.Transition=0","Saline.Native-Saline.Phragmites=0","Saline.Transition-Saline.Phragmites=0"))))

#Only on haplotype I
m1<-lme(NatRichness~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17I,na.action=na.omit,control = lmeControl(maxIter=100,msMaxIter = 100))
m2<-lme(NatRichness~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17I,na.action=na.omit)
anova(m1,m2) #het var is sig
anova(m1,type="marginal")

#Final fig for manuscript
pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Manuscripts/Gradientms/Figs/natrichhetvar.pdf",width=2.2,height=2.2)
ggplot(m1P,aes(x=Transect,y=emmean,color=Transect,group=MarshClassV))+
  labs(x = "",y="Plant richness") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside")+
  geom_line(stat = "identity", position = "identity",size=.5,col="black")+
  geom_point(size=1.8)+
  geom_errorbar(aes(ymax = emmean+SE, ymin=emmean-SE),width=.25,size=.5) +
  scale_color_manual(values = c("gray70","gray50", "gray30"))+
  facet_wrap(vars(MarshClassV),strip.position = "bottom")
dev.off()


### phrag abundance vs richness
#litterkg
ggplot(dat17a,aes(x=phraus,y=Chao1Bac,color=Transect,group=MarshClassV))+
  labs(x = "",y="") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),panel.spacing=unit(0,"cm"),strip.placement = "outside")+
  geom_point(size=1.8)+
  geom_smooth(method="lm",se=F)+
  facet_wrap(vars(MarshClassV),strip.position = "bottom")

temp<-dat17a%>%
  filter(MarshClassV=="Saline")
rcorr(cbind(phraus=temp$phraus,litterkg=temp$litter,NatAbun=temp$NatAbun,Chao1ITS=temp$Chao1ITS,Chao1Bac=temp$Chao1Bac))


### Old fig
forplot<-dat17%>%
  group_by(MarshClassV,Transect)%>%
  summarise(mean=mean(NatRichness),se=std.error(NatRichness))

forplot<-dat17%>%
  group_by(Site,Transect)%>%
  summarise(mean=mean(Shannon),std=std.error(Shannon))

ggplot(forplot,aes(x=Transect,y=mean,group=Site,color=Transect))+
  labs(x = "",y="Richness")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5))+
  geom_point(size=3)+
  geom_line(stat = "identity", position = "identity",size=.8,color="black")+
  geom_errorbar(aes(ymax = mean+std, ymin=mean-std),width=.15,size=.8) +
  facet_wrap(~Site,nrow=1)#,scales="free"




###### Plant Shannon #####

#Final stats model for manuscript. adding heteroscedasticity is significant and helps the tukey results

m1<-lme(Shannon~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17)
#m2<-lme(Shannon~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17)
anova(m1,m2) #het var not sig
anova(m1,type="marginal")
summary(glht(m1,linfct=mcp(Transect="Tukey")))
summary(glht(m1,linfct=mcp(MarshClassV="Tukey")))
hist(resid(m1,type="normalized"))
plot(fitted(m1),resid(m1,type="normalized"))

m1P<-as.data.frame(summary(emmeans(m1,~Transect|MarshClassV)))

m2<-lme(Shannon~MarshClassV.Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17)#,weights=varIdent(form=~1|MarshClassV.Transect)
#summary(glht(m1, linfct = mcp(Transect = "Tukey")))
summary(glht(m2, linfct = mcp(MarshClassV.Transect=c("Fresh.Native-Fresh.Transition=0","Fresh.Native-Fresh.Phragmites=0","Fresh.Transition-Fresh.Phragmites=0","Brackish.Native-Brackish.Transition=0","Brackish.Native-Brackish.Phragmites=0","Brackish.Transition-Brackish.Phragmites=0","Saline.Native-Saline.Transition=0","Saline.Native-Saline.Phragmites=0","Saline.Transition-Saline.Phragmites=0"))))

#Only on haplotype I
m1<-lme(NatShannon~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17I,na.action=na.omit,control = lmeControl(maxIter=100,msMaxIter = 100))
m2<-lme(NatShannon~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17I,na.action=na.omit)
anova(m1,m2) #het var is sig
anova(m1,type="marginal")

#Final fig for manuscript
pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Manuscripts/Gradientms/Figs/shannonhetvar.pdf",width=2.2,height=2.2)
ggplot(m1P,aes(x=Transect,y=emmean,color=Transect,group=MarshClassV))+
  labs(x = "",y="Plant Shannon diversity") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside")+
  geom_line(stat = "identity", position = "identity",size=.5,col="black")+
  geom_point(size=1.8)+
  geom_errorbar(aes(ymax = emmean+SE, ymin=emmean-SE),width=.25,size=.5) +
  scale_color_manual(values = c("gray70","gray50", "gray30"))+
  facet_wrap(vars(MarshClassV),strip.position = "bottom")
dev.off()



###### Plant Simpson #####

#Final stats model for manuscript. adding heteroscedasticity is significant and helps the tukey results

m1<-lme(Simpson~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17)
#m2<-lme(Simpson~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17,control = lmeControl(maxIter=100,msMaxIter = 100))
anova(m1,m2) #het var not sig
anova(m1,type="marginal")
summary(glht(m1,linfct=mcp(MarshClassV="Tukey")))
hist(resid(m1,type="normalized"))
plot(fitted(m1),resid(m1,type="normalized"))

m1P<-as.data.frame(summary(emmeans(m1,~Transect|MarshClassV)))

m2<-lme(Simpson~MarshClassV.Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17)#,weights=varIdent(form=~1|MarshClassV.Transect)
#summary(glht(m1, linfct = mcp(Transect = "Tukey")))
summary(glht(m2, linfct = mcp(MarshClassV.Transect=c("Fresh.Native-Fresh.Transition=0","Fresh.Native-Fresh.Phragmites=0","Fresh.Transition-Fresh.Phragmites=0","Brackish.Native-Brackish.Transition=0","Brackish.Native-Brackish.Phragmites=0","Brackish.Transition-Brackish.Phragmites=0","Saline.Native-Saline.Transition=0","Saline.Native-Saline.Phragmites=0","Saline.Transition-Saline.Phragmites=0"))))

#Only on haplotype I
m1<-lme(NatShannon~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17I,na.action=na.omit,control = lmeControl(maxIter=100,msMaxIter = 100))
m2<-lme(NatShannon~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17I,na.action=na.omit)
anova(m1,m2) #het var is sig
anova(m1,type="marginal")

#Final fig for manuscript
pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Manuscripts/Gradientms/Figs/simpsonhetvar.pdf",width=2.2,height=2.2)
ggplot(m1P,aes(x=Transect,y=emmean,color=Transect,group=MarshClassV))+
  labs(x = "",y="Plant Simpson diversity") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside")+
  geom_line(stat = "identity", position = "identity",size=.5,col="black")+
  geom_point(size=1.8)+
  geom_errorbar(aes(ymax = emmean+SE, ymin=emmean-SE),width=.25,size=.5) +
  scale_color_manual(values = c("gray70","gray50", "gray30"))+
  facet_wrap(vars(MarshClassV),strip.position = "bottom")
dev.off()




###### Plant InvSimpson #####

#Final stats model for manuscript. adding heteroscedasticity is significant and helps the tukey results

#m1<-lme(InvSimpson~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17)
m1<-lme(InvSimpson~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17)
anova(m1,m2) #het var sig
anova(m1,type="marginal")
summary(glht(m1,linfct=mcp(Transect="Tukey")))
hist(resid(m1,type="normalized"))
plot(fitted(m1),resid(m1,type="normalized"))

m1P<-as.data.frame(summary(emmeans(m1,~Transect|MarshClassV)))

m2<-lme(InvSimpson~MarshClassV.Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17)#
#summary(glht(m1, linfct = mcp(Transect = "Tukey")))
summary(glht(m2, linfct = mcp(MarshClassV.Transect=c("Fresh.Native-Fresh.Transition=0","Fresh.Native-Fresh.Phragmites=0","Fresh.Transition-Fresh.Phragmites=0","Brackish.Native-Brackish.Transition=0","Brackish.Native-Brackish.Phragmites=0","Brackish.Transition-Brackish.Phragmites=0","Saline.Native-Saline.Transition=0","Saline.Native-Saline.Phragmites=0","Saline.Transition-Saline.Phragmites=0"))))

#Only on haplotype I
m1<-lme(NatShannon~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17I,na.action=na.omit,control = lmeControl(maxIter=100,msMaxIter = 100))
m2<-lme(NatShannon~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17I,na.action=na.omit)
anova(m1,m2) #het var is sig
anova(m1,type="marginal")

#Final fig for manuscript
pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Manuscripts/Gradientms/Figs/invsimpsonhetvar.pdf",width=2.2,height=2.2)
ggplot(m1P,aes(x=Transect,y=emmean,color=Transect,group=MarshClassV))+
  labs(x = "",y="Plant inverse Simpson") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside")+
  geom_line(stat = "identity", position = "identity",size=.5,col="black")+
  geom_point(size=1.8)+
  geom_errorbar(aes(ymax = emmean+SE, ymin=emmean-SE),width=.25,size=.5) +
  scale_color_manual(values = c("gray70","gray50", "gray30"))+
  facet_wrap(vars(MarshClassV),strip.position = "bottom")
dev.off()



##### Fungi ####

###### Fungal Richness ######

#m1<-lme(Chao1ITS~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17,na.action=na.omit)
m1<-lme(Chao1ITS~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17,na.action=na.omit)
anova(m1,m2) #het var not sig
anova(m1,type="marginal")
hist(resid(m1,type="normalized"))
plot(fitted(m1),resid(m1,type="normalized"))

m1F<-as.data.frame(summary(emmeans(m1,~Transect|MarshClassV)))

m2<-lme(Chao1ITS~MarshClassV.Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17,na.action=na.omit)
#summary(glht(m1, linfct = mcp(Transect = "Tukey")))
summary(glht(m2, linfct = mcp(MarshClassV.Transect=c("Fresh.Native-Fresh.Transition=0","Fresh.Native-Fresh.Phragmites=0","Fresh.Transition-Fresh.Phragmites=0","Brackish.Native-Brackish.Transition=0","Brackish.Native-Brackish.Phragmites=0","Brackish.Transition-Brackish.Phragmites=0","Saline.Native-Saline.Transition=0","Saline.Native-Saline.Phragmites=0","Saline.Transition-Saline.Phragmites=0"))))

#Only on haplotype I
m1<-lme(Chao1ITS~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17I,na.action=na.omit)
m2<-lme(Chao1ITS~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17I,na.action=na.omit)
anova(m1,m2) #het var is not sig
anova(m2,type="marginal")

#Final fig for manuscript
pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Manuscripts/Gradientms/Figs/fungirich.pdf",width=2.2,height=2.2)
ggplot(m1F,aes(x=Transect,y=emmean,color=Transect,group=MarshClassV))+
  labs(x = "",y="Fungal richness") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside")+
  geom_line(stat = "identity", position = "identity",size=.5,col="black")+
  geom_point(size=1.8)+
  geom_errorbar(aes(ymax = emmean+SE, ymin=emmean-SE),width=.25,size=.5) +
  scale_color_manual(values = c("gray70","gray50", "gray30"))+
  facet_wrap(vars(MarshClassV),strip.position = "bottom")
dev.off()


###Old analyses/figs

richITS2means<-dat17%>%
  group_by(Site,Transect)%>%
  summarise(mean=mean(Chao1ITS,na.rm=T),std=std.error(Chao1ITS,na.rm=T))

ggplot(richITS2means,aes(x=Transect,y=mean,group=Site,color=Transect))+
  labs(x = "",y="Richness")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5))+
  geom_point(size=3)+
  geom_line(stat = "identity", position = "identity",size=.8,color="black")+
  geom_errorbar(aes(ymax = mean+std, ymin=mean-std),width=.15,size=.8) +
  facet_wrap(~Site,nrow=1)#,scales="free"

m1<-lme(Chao1ITS~Transect*MarshClassV,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17,na.action = na.omit)
anova(m1,type="marginal")
summary(glht(m1,linfct=mcp(Transect="Tukey")))




###### Fungal Shannon ######

m1<-lme(ShannonITS~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17,na.action=na.omit)
#m2<-lme(ShannonITS~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17,na.action=na.omit)
anova(m1,m2) #het var not sig
anova(m1,type="marginal")
hist(resid(m1,type="normalized"))
plot(fitted(m1),resid(m1,type="normalized"))

m1F<-as.data.frame(summary(emmeans(m1,~Transect|MarshClassV)))

m2<-lme(ShannonITS~MarshClassV.Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17,na.action=na.omit)
#summary(glht(m1, linfct = mcp(Transect = "Tukey")))
summary(glht(m2, linfct = mcp(MarshClassV.Transect=c("Fresh.Native-Fresh.Transition=0","Fresh.Native-Fresh.Phragmites=0","Fresh.Transition-Fresh.Phragmites=0","Brackish.Native-Brackish.Transition=0","Brackish.Native-Brackish.Phragmites=0","Brackish.Transition-Brackish.Phragmites=0","Saline.Native-Saline.Transition=0","Saline.Native-Saline.Phragmites=0","Saline.Transition-Saline.Phragmites=0"))))

#Only on haplotype I
m1<-lme(Chao1ITS~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17I,na.action=na.omit)
m2<-lme(Chao1ITS~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17I,na.action=na.omit)
anova(m1,m2) #het var is not sig
anova(m2,type="marginal")

#Final fig for manuscript
pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Manuscripts/Gradientms/Figs/fungishannon.pdf",width=2.2,height=2.2)
ggplot(m1F,aes(x=Transect,y=emmean,color=Transect,group=MarshClassV))+
  labs(x = "",y="Fungal Shannon diversity") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside")+
  geom_line(stat = "identity", position = "identity",size=.5,col="black")+
  geom_point(size=1.8)+
  geom_errorbar(aes(ymax = emmean+SE, ymin=emmean-SE),width=.25,size=.5) +
  scale_color_manual(values = c("gray70","gray50", "gray30"))+
  facet_wrap(vars(MarshClassV),strip.position = "bottom")
dev.off()



###### Fungal Simpson ######

#m1<-lme(SimpsonITS~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17,na.action=na.omit)
m1<-lme(SimpsonITS~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17,na.action=na.omit)
anova(m1,m2) #het var sig
anova(m1,type="marginal")
hist(resid(m1,type="normalized"))
plot(fitted(m1),resid(m1,type="normalized"))

m1F<-as.data.frame(summary(emmeans(m1,~Transect|MarshClassV)))

m2<-lme(SimpsonITS~MarshClassV.Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17,na.action=na.omit)
#summary(glht(m1, linfct = mcp(Transect = "Tukey")))
summary(glht(m2, linfct = mcp(MarshClassV.Transect=c("Fresh.Native-Fresh.Transition=0","Fresh.Native-Fresh.Phragmites=0","Fresh.Transition-Fresh.Phragmites=0","Brackish.Native-Brackish.Transition=0","Brackish.Native-Brackish.Phragmites=0","Brackish.Transition-Brackish.Phragmites=0","Saline.Native-Saline.Transition=0","Saline.Native-Saline.Phragmites=0","Saline.Transition-Saline.Phragmites=0"))))

#Only on haplotype I
m1<-lme(Chao1ITS~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17I,na.action=na.omit)
m2<-lme(Chao1ITS~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17I,na.action=na.omit)
anova(m1,m2) #het var is not sig
anova(m2,type="marginal")

#Final fig for manuscript
pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Manuscripts/Gradientms/Figs/fungisimpsonhetvar.pdf",width=2.2,height=2.2)
ggplot(m1F,aes(x=Transect,y=emmean,color=Transect,group=MarshClassV))+
  labs(x = "",y="Fungal Simpson diversity") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside")+
  geom_line(stat = "identity", position = "identity",size=.5,col="black")+
  geom_point(size=1.8)+
  geom_errorbar(aes(ymax = emmean+SE, ymin=emmean-SE),width=.25,size=.5) +
  scale_color_manual(values = c("gray70","gray50", "gray30"))+
  facet_wrap(vars(MarshClassV),strip.position = "bottom")
dev.off()




###### Fungal InvSimpson ######

#m1<-lme(InvSimpsonITS~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17,na.action=na.omit)
m1<-lme(InvSimpsonITS~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17,na.action=na.omit,control = lmeControl(msMaxIter = 10000,tolerance=1e-06,warnOnly= T))
anova(m1,m2) #het var sig
anova(m1,type="marginal")
hist(resid(m1,type="normalized"))
plot(fitted(m1),resid(m1,type="normalized"))

m1F<-as.data.frame(summary(emmeans(m1,~Transect|MarshClassV)))

m2<-lme(InvSimpsonITS~MarshClassV.Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17,na.action=na.omit,control = lmeControl(msMaxIter = 10000,tolerance=1e-06,warnOnly= T))
#summary(glht(m1, linfct = mcp(Transect = "Tukey")))
summary(glht(m2, linfct = mcp(MarshClassV.Transect=c("Fresh.Native-Fresh.Transition=0","Fresh.Native-Fresh.Phragmites=0","Fresh.Transition-Fresh.Phragmites=0","Brackish.Native-Brackish.Transition=0","Brackish.Native-Brackish.Phragmites=0","Brackish.Transition-Brackish.Phragmites=0","Saline.Native-Saline.Transition=0","Saline.Native-Saline.Phragmites=0","Saline.Transition-Saline.Phragmites=0"))))

#Only on haplotype I
m1<-lme(Chao1ITS~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17I,na.action=na.omit)
m2<-lme(Chao1ITS~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17I,na.action=na.omit)
anova(m1,m2) #het var is not sig
anova(m2,type="marginal")

#Final fig for manuscript
pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Manuscripts/Gradientms/Figs/fungiinvsimpsonhetvar.pdf",width=2.2,height=2.2)
ggplot(m1F,aes(x=Transect,y=emmean,color=Transect,group=MarshClassV))+
  labs(x = "",y="Fungal inverse Simpson") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside")+
  geom_line(stat = "identity", position = "identity",size=.5,col="black")+
  geom_point(size=1.8)+
  geom_errorbar(aes(ymax = emmean+SE, ymin=emmean-SE),width=.25,size=.5) +
  scale_color_manual(values = c("gray70","gray50", "gray30"))+
  facet_wrap(vars(MarshClassV),strip.position = "bottom")
dev.off()







##### Bacteria #####

###### Bacteria richness #####

#dat17a (use ths b/c the 9000 richness plot is removed)
dat17a$MarshClassV.Transect<-interaction(dat17a$MarshClassV,dat17a$Transect)
dat17a$MarshClassV.Transect<-factor(dat17a$MarshClassV.Transect,levels=c("Fresh.Native","Fresh.Transition","Fresh.Phragmites","Brackish.Native","Brackish.Transition","Brackish.Phragmites","Saline.Native","Saline.Transition","Saline.Phragmites"))

#m1<-lme(Chao1Bac~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17a,na.action=na.omit)
m2<-lme(Chao1Bac~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17a,na.action=na.omit)
anova(m1,m2) #sig.
anova(m2,type="marginal")
summary(glht(m2,linfct=mcp(MarshClassV="Tukey")))
hist(resid(m2,type="normalized"))
plot(fitted(m2),resid(m2,type="normalized"))

m1B<-as.data.frame(summary(emmeans(m2,~Transect|MarshClassV)))

m2<-lme(Chao1Bac~MarshClassV.Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17a,na.action=na.omit)
#summary(glht(m1, linfct = mcp(Transect = "Tukey")))
summary(glht(m2, linfct = mcp(MarshClassV.Transect=c("Fresh.Native-Fresh.Transition=0","Fresh.Native-Fresh.Phragmites=0","Fresh.Transition-Fresh.Phragmites=0","Brackish.Native-Brackish.Transition=0","Brackish.Native-Brackish.Phragmites=0","Brackish.Transition-Brackish.Phragmites=0","Saline.Native-Saline.Transition=0","Saline.Native-Saline.Phragmites=0","Saline.Transition-Saline.Phragmites=0"))))

#Only on haplotype I
m1<-lme(Chao1Bac~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17aI,na.action=na.omit)
m2<-lme(Chao1Bac~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17aI,na.action=na.omit)
anova(m1,m2) #het var is sig
anova(m1,type="marginal")

#Final fig for manuscript
pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Manuscripts/Gradientms/Figs/bacteriarichhetvar.pdf",width=2.2,height=2.2)
ggplot(m1B,aes(x=Transect,y=emmean,color=Transect,group=MarshClassV))+
  labs(x = "",y="Bacterial richness") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside")+
  geom_line(stat = "identity", position = "identity",size=.5,col="black")+
  geom_point(size=1.8)+
  geom_errorbar(aes(ymax = emmean+SE, ymin=emmean-SE),width=.25,size=.5) +
  scale_color_manual(values = c("gray70","gray50", "gray30"))+
  facet_wrap(vars(MarshClassV),strip.position = "bottom")
dev.off()



## Old analyses/figs
richbac2means<-dat17%>%
  filter(Chao1Bac<9000)%>% #one native smaple from LUMCON 2 (168) is huge (it was the samples with 2 million reads)
  group_by(Site,Transect)%>%
  summarise(mean=mean(Chao1Bac),std=std.error(Chao1Bac))

ggplot(richbac2means,aes(x=Transect,y=mean,group=Site,color=Transect))+
  labs(x = "",y="Richness")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5))+
  geom_point(size=3)+
  geom_line(stat = "identity", position = "identity",size=.8,color="black")+
  geom_errorbar(aes(ymax = mean+std, ymin=mean-std),width=.15,size=.8) +
  facet_wrap(~Site,nrow=1)#,scales="free"

m1<-lme(Chao1Bac~Transect*MarshClassS,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17,na.action = na.omit)
anova(m1,type="marginal")
summary(glht(m1,linfct=mcp(MarshClassS="Tukey")))

m1<-lme(Chao1Bac~Transect*MarshClassV,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17,na.action = na.omit)
anova(m1,type="marginal")
summary(glht(m1,linfct=mcp(MarshClassS="Tukey")))



dat17a$Taxon<-"Plant"
richITS2means$Taxon<-"Fungi"
richbac2means$Taxon<-"Bacteria"


plotrichplants <- 
  ggplot(dat17a,aes(x=Transect,y=mean,group=Site,color=Transect))+
  labs(x = "",y="Plant Richness")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5))+
  geom_point(size=3)+
  geom_line(stat = "identity", position = "identity",size=.8,color="black")+
  geom_errorbar(aes(ymax = mean+std, ymin=mean-std),width=.15,size=.8) +
  facet_wrap(~Site,nrow=1)

plotrichITS <- 
  ggplot(richITS2means,aes(x=Transect,y=mean,group=Site,color=Transect))+
  labs(x = "",y="Fungi Richness")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5))+
  geom_point(size=3)+
  geom_line(stat = "identity", position = "identity",size=.8,color="black")+
  geom_errorbar(aes(ymax = mean+std, ymin=mean-std),width=.15,size=.8) +
  facet_wrap(~Site,nrow=1)

plotrichbact <- 
  ggplot(richbac2means,aes(x=Transect,y=mean,group=Site,color=Transect))+
  labs(x = "",y="Bacteria Richness")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5))+
  geom_point(size=3)+
  geom_line(stat = "identity", position = "identity",size=.8,color="black")+
  geom_errorbar(aes(ymax = mean+std, ymin=mean-std),width=.15,size=.8) +
  facet_wrap(~Site,nrow=1)

pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Stats/Gradient/Figs/richnesspfb.pdf")
plot_grid(plotrichplants, plotrichITS, plotrichbact, nrow = 3)
dev.off()


#old plot trying to not use cowplot but it doesn't work b/c you need scales ="free" for between taxa but not within taxa
richmeans<-rbind(dat17a,richbac2means,richITS2means)
ggplot(richmeans,aes(x=Transect,y=mean,group=Site,color=Transect))+
  labs(x = "",y="Richness")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5))+
  geom_point(size=3)+
  geom_line(stat = "identity", position = "identity",size=.8,color="black")+
  geom_errorbar(aes(ymax = mean+std, ymin=mean-std),width=.15,size=.8) +
  facet_wrap(~Taxon*Site,nrow=3)#,scales="free"




###### Bacteria Shannon #####

#dat17a (use ths b/c the 9000 richness plot is removed)

m1<-lme(ShannonBac~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17a,na.action=na.omit)
#m2<-lme(ShannonBac~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17a,na.action=na.omit)
anova(m1,m2) #not sig.
anova(m1,type="marginal")
summary(glht(m1,linfct=mcp(MarshClassV="Tukey")))
hist(resid(m1,type="normalized"))
plot(fitted(m1),resid(m1,type="normalized"))

m1B<-as.data.frame(summary(emmeans(m1,~Transect|MarshClassV)))

m2<-lme(ShannonBac~MarshClassV.Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17a,na.action=na.omit)#,weights=varIdent(form=~1|MarshClassV.Transect)
#summary(glht(m1, linfct = mcp(Transect = "Tukey")))
summary(glht(m2, linfct = mcp(MarshClassV.Transect=c("Fresh.Native-Fresh.Transition=0","Fresh.Native-Fresh.Phragmites=0","Fresh.Transition-Fresh.Phragmites=0","Brackish.Native-Brackish.Transition=0","Brackish.Native-Brackish.Phragmites=0","Brackish.Transition-Brackish.Phragmites=0","Saline.Native-Saline.Transition=0","Saline.Native-Saline.Phragmites=0","Saline.Transition-Saline.Phragmites=0"))))

#Only on haplotype I
m1<-lme(Chao1Bac~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17aI,na.action=na.omit)
m2<-lme(Chao1Bac~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17aI,na.action=na.omit)
anova(m1,m2) #het var is sig
anova(m1,type="marginal")

#Final fig for manuscript
pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Manuscripts/Gradientms/Figs/bacteriashannonhetvar.pdf",width=2.2,height=2.2)
ggplot(m1B,aes(x=Transect,y=emmean,color=Transect,group=MarshClassV))+
  labs(x = "",y="Bacterial Shannon diversity") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside")+
  geom_line(stat = "identity", position = "identity",size=.5,col="black")+
  geom_point(size=1.8)+
  geom_errorbar(aes(ymax = emmean+SE, ymin=emmean-SE),width=.25,size=.5) +
  scale_color_manual(values = c("gray70","gray50", "gray30"))+
  facet_wrap(vars(MarshClassV),strip.position = "bottom")
dev.off()



###### Bacteria Simpson #####

#dat17a (use ths b/c the 9000 richness plot is removed)

#m1<-lme(SimpsonBac~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17a,na.action=na.omit)
m1<-lme(SimpsonBac~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17a,na.action=na.omit)
anova(m1,m2) # sig.
anova(m1,type="marginal")
summary(glht(m1,linfct=mcp(Transect="Tukey")))
hist(resid(m1,type="normalized"))
plot(fitted(m1),resid(m1,type="normalized"))

m1B<-as.data.frame(summary(emmeans(m1,~Transect|MarshClassV)))

m2<-lme(SimpsonBac~MarshClassV.Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17a,na.action=na.omit,control = lmeControl(msMaxIter = 10000,tolerance=1e-06,warnOnly= T))#
AIC(m1,m2)
#summary(glht(m1, linfct = mcp(Transect = "Tukey")))
summary(glht(m2, linfct = mcp(MarshClassV.Transect=c("Fresh.Native-Fresh.Transition=0","Fresh.Native-Fresh.Phragmites=0","Fresh.Transition-Fresh.Phragmites=0","Brackish.Native-Brackish.Transition=0","Brackish.Native-Brackish.Phragmites=0","Brackish.Transition-Brackish.Phragmites=0","Saline.Native-Saline.Transition=0","Saline.Native-Saline.Phragmites=0","Saline.Transition-Saline.Phragmites=0"))))

#Only on haplotype I
m1<-lme(Chao1Bac~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17aI,na.action=na.omit)
m2<-lme(Chao1Bac~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17aI,na.action=na.omit)
anova(m1,m2) #het var is sig
anova(m1,type="marginal")

#Final fig for manuscript
pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Manuscripts/Gradientms/Figs/bacteriasimpsonhetvar.pdf",width=2.2,height=2.2)
ggplot(m1B,aes(x=Transect,y=emmean,color=Transect,group=MarshClassV))+
  labs(x = "",y="Bacterial Simpson diversity") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside")+
  geom_line(stat = "identity", position = "identity",size=.5,col="black")+
  geom_point(size=1.8)+
  geom_errorbar(aes(ymax = emmean+SE, ymin=emmean-SE),width=.25,size=.5) +
  scale_color_manual(values = c("gray70","gray50", "gray30"))+
  facet_wrap(vars(MarshClassV),strip.position = "bottom")
dev.off()




###### Bacteria InvSimpson #####

#dat17a (use ths b/c the 9000 richness plot is removed)

#m1<-lme(InvSimpsonBac~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17a,na.action=na.omit)
m1<-lme(InvSimpsonBac~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17a,na.action=na.omit)
anova(m1,m2) # sig.
anova(m1,type="marginal")
summary(glht(m1,linfct=mcp(Transect="Tukey")))
hist(resid(m1,type="normalized"))
plot(fitted(m1),resid(m1,type="normalized"))

m1B<-as.data.frame(summary(emmeans(m1,~Transect|MarshClassV)))

m2<-lme(InvSimpsonBac~MarshClassV.Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17a,na.action=na.omit)#,control = lmeControl(msMaxIter = 10000,tolerance=1e-06,warnOnly= T)
AIC(m1,m2)
#summary(glht(m1, linfct = mcp(Transect = "Tukey")))
summary(glht(m2, linfct = mcp(MarshClassV.Transect=c("Fresh.Native-Fresh.Transition=0","Fresh.Native-Fresh.Phragmites=0","Fresh.Transition-Fresh.Phragmites=0","Brackish.Native-Brackish.Transition=0","Brackish.Native-Brackish.Phragmites=0","Brackish.Transition-Brackish.Phragmites=0","Saline.Native-Saline.Transition=0","Saline.Native-Saline.Phragmites=0","Saline.Transition-Saline.Phragmites=0"))))

#Only on haplotype I
m1<-lme(Chao1Bac~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),weights=varIdent(form=~1|MarshClassV.Transect),data=dat17aI,na.action=na.omit)
m2<-lme(Chao1Bac~MarshClassV*Transect,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17aI,na.action=na.omit)
anova(m1,m2) #het var is sig
anova(m1,type="marginal")

#Final fig for manuscript
pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/LAmarsh/Survey/Manuscripts/Gradientms/Figs/bacteriainvsimpsonhetvar.pdf",width=2.2,height=2.2)
ggplot(m1B,aes(x=Transect,y=emmean,color=Transect,group=MarshClassV))+
  labs(x = "",y="Bacterial inverse Simpson") +
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5),legend.position = "none",panel.spacing=unit(0,"cm"),strip.placement = "outside")+
  geom_line(stat = "identity", position = "identity",size=.5,col="black")+
  geom_point(size=1.8)+
  geom_errorbar(aes(ymax = emmean+SE, ymin=emmean-SE),width=.25,size=.5) +
  scale_color_manual(values = c("gray70","gray50", "gray30"))+
  facet_wrap(vars(MarshClassV),strip.position = "bottom")
dev.off()









##### Regressions of plants, bacteria, fungi richness ####

dat17a<-dat17
ind<-which(dat17a$Chao1Bac>9000)
dat17a$Chao1Bac[ind]<-NA
dat17a$RichnessBac[ind]<-NA

plot(jitter(dat17a$NatRichness),dat17a$Chao1Bac)
plot(jitter(dat17a$NatRichness),dat17a$Chao1ITS)
plot(dat17a$Chao1ITS,dat17a$Chao1Bac)

ggplot(dat17a,aes(x=NatRichness,y=Chao1ITS,group=Site,color=Transect))+
  labs(x = "",y="Richness")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5))+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  facet_wrap(~Site,nrow=3)#,scales="free"

rcorr(dat17a$NatRichness,dat17$Chao1ITS)
rcorr(dat17a$NatRichness,dat17$Chao1Bac)
rcorr(dat17a$Chao1ITS,dat17$Chao1Bac)

m1<-lme(Chao1Bac~Chao1ITS,random=~1|Site,correlation=corSpher(form = ~ Lat+Long),data=dat17,na.action = na.omit)
anova(m1,type="marginal")





##### Coordinated Response Analysis - Richness #####
#Options for coordinated response analysis: 

#ordinations
#sum species comp for all native and all phrag plots per site and calculate sorensons dissimilarity.
#do all pairs of native-phrg plot dissimilarity and average
#calculate bray-curtis dissimilarity so that it is more congruent with ordination results
#try the ordinations just with native and phrag plots to see if it is more consistent with the sorensons results

#richness
#difference in means (or log response ratio) between native and phrag plots (somehow calculate an error bar on the difference)

#I like using bray dist rather than sorensons (at least for microbes), even with the lines it matched up better with the ordination results


##### Rich dif on line and then averaging, "a mean of the difference" #####
#Use this for now to be consistent with the composition analysis

#use dat17a b/c the large bacteria richness sample (168) is NA

rich<-dat17a%>%
  filter(Transect!="Transition")%>%
  dplyr::select(Line,Transect, Site, NatRichness,Chao1ITS,Chao1Bac)%>%
  gather(key="Taxon",value="Richness",NatRichness,Chao1ITS,Chao1Bac)%>%
  arrange(Line)%>%
  spread(Transect,Richness)%>%
  #mutate(RichDif=log((Phragmites+1)/(Native+1)))%>%
  mutate(RichDif=Phragmites-Native)%>%
  dplyr::select(Line,Site,Taxon,RichDif)%>%
  spread(Taxon,RichDif)%>%
  dplyr::select(Site,NatRichness,Chao1Bac,Chao1ITS)%>%
  group_by(Site)%>%
  summarise_all(list(~mean(.,na.rm=T),~sd(.,na.rm=T),~std.error(.,na.rm=T)))

data.frame(rich)  
plot(rich$RichnessBac,rich$RichnessITS,col=rich$Site)
plot(rich$RichnessBac_mean,rich$RichnessITS_mean)

ggplot(rich,aes(x=NatRichness_mean,y=Chao1ITS_mean))+
  labs(x = "Richness Difference Plants",y="Richness Difference Fungi")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5))+
  geom_point(size=2)+
  geom_errorbar(aes(ymax = Chao1ITS_mean+Chao1ITS_std.error, ymin=Chao1ITS_mean-Chao1ITS_std.error),width=.15,size=.8)+#.5
  geom_errorbarh(aes(xmax=NatRichness_mean+NatRichness_std.error,xmin=NatRichness_mean-NatRichness_std.error,height=.15),size=.8)+#.5
  guides(col = guide_legend(ncol = 1))


m1<-deming(rich$Chao1ITS_mean~rich$NatRichness_mean,xstd=rich$NatRichness_sd,ystd = rich$Chao1ITS_sd)#,cv=T or F gives same exact results
m1
m2<-deming(rich$Chao1Bac_mean~rich$NatRichness_mean,xstd=rich$NatRichness_sd,ystd = rich$Chao1Bac_sd)#,cv=T or F gives same exact results
m2
m3<-deming(rich$Chao1Bac_mean~rich$Chao1ITS_mean,xstd=rich$Chao1ITS_sd,ystd = rich$Chao1Bac_sd)#,cv=T or F gives same exact results
m3



#Only haplotype I
richI<-rich[1:6,]
m1<-deming(richI$Chao1ITS_mean~richI$NatRichness_mean,xstd=richI$NatRichness_sd,ystd = richI$Chao1ITS_sd)#,cv=T or F gives same exact results
m1
m2<-deming(richI$Chao1Bac_mean~richI$NatRichness_mean,xstd=richI$NatRichness_sd,ystd = richI$Chao1Bac_sd)#,cv=T or F gives same exact results
m2
m3<-deming(richI$Chao1Bac_mean~richI$Chao1ITS_mean,xstd=richI$Chao1ITS_sd,ystd = richI$Chao1Bac_sd)#,cv=T or F gives same exact results
m3



##### Averaging by transect then difference, "a difference of the means" #####
richb<-dat17a%>%
  filter(Transect!="Transition")%>%
  dplyr::select(Transect,Site,NatRichness,Chao1ITS,Chao1Bac)%>%
  gather(Taxon,Richness,NatRichness,Chao1ITS,Chao1Bac)%>%
  group_by(Taxon,Site,Transect)%>%
  summarise(mean=mean(Richness,na.rm=T),sd=sd(Richness,na.rm=T),count=sum(!is.na(Richness)))%>%#,se=std.error(Richness,na.rm=T)
  gather(variable,value,mean,sd,count)%>%
  unite(temp,Transect,variable)%>%
  spread(temp,value)%>%
  mutate(dif=Phragmites_mean-Native_mean,se=sqrt(Phragmites_sd^2/Phragmites_count+Native_sd^2/Native_count),sd=sqrt(((Native_count-1)*Native_sd^2+(Phragmites_count-1)*Phragmites_sd^2)/(Native_count+Phragmites_count-2)))%>%
  ungroup()%>%
  dplyr::select(Taxon,Site,dif,sd,se)%>%
  mutate(Taxon = plyr::revalue(Taxon,c(NatRichness="Plant",Chao1ITS="Fungi",Chao1Bac="Bacteria")))%>%
  gather(variable,value,dif,sd,se)%>%
  unite(temp,Taxon,variable)%>%
  spread(temp,value)

as.data.frame(richb)

plot(richb$Bacteria_dif,richb$Fungi_dif)

ggplot(richb,aes(x=Plant_dif,y=Fungi_dif))+
  labs(x = "Richness Difference Plants",y="Richness Difference Fungi")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5))+
  geom_point(size=2)+
  geom_errorbar(aes(ymax = Fungi_dif+Fungi_se, ymin=Fungi_dif-Fungi_se),width=.15,size=.8)+#.5
  geom_errorbarh(aes(xmax=Plant_dif+Plant_se,xmin=Plant_dif-Plant_se,height=.15),size=.8)+#.5
  guides(col = guide_legend(ncol = 1))

ggplot(richb,aes(x=Fungi_dif,y=Bacteria_dif))+
  labs(x = "Richness Difference Fungi",y="Richness Difference Bacteria")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5))+
  geom_point(size=2)+
  geom_errorbar(aes(ymax = Bacteria_dif+Bacteria_se, ymin=Bacteria_dif-Bacteria_se),width=11,size=.8)+#
  geom_errorbarh(aes(xmax=Fungi_dif+Fungi_se,xmin=Fungi_dif-Fungi_se),height=45,size=.8)+#
  guides(col = guide_legend(ncol = 1))














####### FUNGI - Looking at composition at the phylum level #####
datITS2rcsoilphylum <- tax_glom(datITS2rcsoil, taxrank="Rank2")
tax_table(datITS2rcsoilphylum)

tempphy<-datITS2rcsoilphylum%>%
  #subset_samples(Site==sites[i])%>%
  filter_taxa(function(x) sum(x>0) >2, prune=T)
mynmds <- ordinate(tempphy, "CAP", "bray",formula=as.formula(~Transect+Condition(Site)))
mynmds <- ordinate(tempphy, "CAP", "bray",formula=as.formula(~Transect*Site))
plot_ordination(tempphy, mynmds, type="samples", color="Transect")+
  theme_classic()+
  #theme(legend.position = "none")+
  stat_ellipse(geom = "polygon", type="t",alpha=0.0,level=.95)# , aes(fill=Site),

anova(mynmds,by = "terms",permutations = how(nperm=1000))


datITS2rcsoilphylum2 <- tax_glom(datITS2rcsoil, taxrank="Rank2")%>%
  filter_taxa(function(x) sum(x>0) >10, prune=T)


options(contrasts=c("contr.helmert","contr.poly"));options("contrasts")

anovaoutput<-data.frame(taxon=rep(NA,ntaxa(datITS2rcsoilphylum2)),site=rep(NA,ntaxa(datITS2rcsoilphylum2)),transect=rep(NA,ntaxa(datITS2rcsoilphylum2)),interaction=rep(NA,ntaxa(datITS2rcsoilphylum2)))
for(i in 1:ntaxa(datITS2rcsoilphylum2)){
  temptaxa<-tax_table(datITS2rcsoilphylum2)[i,"Rank2"]
  tempdata<-data.frame(otu_table(datITS2rcsoilphylum2)[,i],sample_data(datITS2rcsoilorder2))
  colnames(tempdata)[1]<-"OTU"
  m1<-gls(OTU~Site*Transect,data=tempdata)
  anovaoutput[i,1]<-temptaxa
  anovaoutput[i,2]<-anova(m1,type="marginal")[2,"p-value"]
  anovaoutput[i,3]<-anova(m1,type="marginal")[3,"p-value"]
  anovaoutput[i,4]<-anova(m1,type="marginal")[4,"p-value"]
  #hist(resid(m1))
}
anovaoutput$siteq<-p.adjust(anovaoutput$site,method="fdr")
anovaoutput$transectq<-p.adjust(anovaoutput$transect,method="fdr")
anovaoutput$interactionq<-p.adjust(anovaoutput$interaction,method="fdr")

anovaoutput[which(anovaoutput$transect<0.05),]
anovaoutput[which(anovaoutput$interaction<0.05),]

temp<-otu_table(datITS2rcsoilphylum2)
colnames(temp)<-tax_table(datITS2rcsoilphylum2)[,"Rank2"]
datITS2rcsoilphylum3<-data.frame(temp,sample_data(datITS2rcsoilphylum2))

means<-datITS2rcsoilphylum3%>%
  group_by(Site,Transect)%>%
  summarize(mean=mean(p__Glomeromycota),std=std.error(p__Glomeromycota))
means$Site<-factor(means$Site,levels=c("Barataria","Turtle Cove","Pearl River","Fontainebleau","Big Branch","Bayou Sauvage","LUMCON 1","LUMCON 2"))

ggplot(means,aes(x=Transect,y=mean,group=Site,color=Transect))+
  labs(x = "",y="Relative abundance")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5))+
  geom_point(size=3)+
  geom_line(stat = "identity", position = "identity",size=.8,color="black")+
  geom_errorbar(aes(ymax = mean+std, ymin=mean-std),width=.15,size=.8) +
  facet_wrap(~Site,nrow=1)#,scales="free"

means<-datITS2rcsoilorder3%>%
  group_by(Transect)%>%
  summarize(mean=mean(o__Cystobasidiales),std=std.error(o__Cystobasidiales))

ggplot(means,aes(x=Transect,y=mean,group=Transect,color=Transect))+
  labs(x = "",y="Relative abundance")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5))+
  geom_point(size=3)+
  geom_line(stat = "identity", position = "identity",size=.8,color="black")+
  geom_errorbar(aes(ymax = mean+std, ymin=mean-std),width=.15,size=.8) 




####### Looking at composition at the class level #####
temp<-left_join(sample_data(datITS2rcsoil),sals2)
sample_data(datITS2rcsoil)$Salinity<-temp$MeanSalinityppt
sample_data(datITS2rcsoil)$Salinitylog<-log(temp$MeanSalinityppt)

datITS2rcsoilclass <- tax_glom(datITS2rcsoil, taxrank="Rank3")
tax_table(datITS2rcsoilclass)

tempphy<-datITS2rcsoilclass%>%
  #subset_samples(Site==sites[i])%>%
  filter_taxa(function(x) sum(x>0) >2, prune=T)
mynmds <- ordinate(tempphy, "CAP", "bray",formula=as.formula(~Transect+Condition(Site)))
mynmds <- ordinate(tempphy, "CAP", "bray",formula=as.formula(~Transect*Salinitylog))
plot_ordination(tempphy, mynmds, type="samples", color="Transect")+
  theme_classic()+
  #theme(legend.position = "none")+
  stat_ellipse(geom = "polygon", type="t",alpha=0.0,level=.95)# , aes(fill=Site),

anova(mynmds,by = "terms",permutations = how(nperm=10000))


datITS2rcsoilorder2 <- tax_glom(datITS2rcsoil, taxrank="Rank4")%>%
  filter_taxa(function(x) sum(x>0) >10, prune=T)


options(contrasts=c("contr.helmert","contr.poly"));options("contrasts")

anovaoutput<-data.frame(taxon=rep(NA,ntaxa(datITS2rcsoilorder2)),site=rep(NA,ntaxa(datITS2rcsoilorder2)),transect=rep(NA,ntaxa(datITS2rcsoilorder2)),interaction=rep(NA,ntaxa(datITS2rcsoilorder2)))
for(i in 1:ntaxa(datITS2rcsoilorder2)){
  temptaxa<-tax_table(datITS2rcsoilorder2)[i,"Rank4"]
  tempdata<-data.frame(otu_table(datITS2rcsoilorder2)[,i],sample_data(datITS2rcsoilorder2))
  colnames(tempdata)[1]<-"OTU"
  m1<-gls(OTU~Site*Transect,data=tempdata)
  anovaoutput[i,1]<-temptaxa
  anovaoutput[i,2]<-anova(m1,type="marginal")[2,"p-value"]
  anovaoutput[i,3]<-anova(m1,type="marginal")[3,"p-value"]
  anovaoutput[i,4]<-anova(m1,type="marginal")[4,"p-value"]
  #hist(resid(m1))
}
anovaoutput$siteq<-p.adjust(anovaoutput$site,method="fdr")
anovaoutput$transectq<-p.adjust(anovaoutput$transect,method="fdr")
anovaoutput$interactionq<-p.adjust(anovaoutput$interaction,method="fdr")

anovaoutput[which(anovaoutput$transect<0.05),]
anovaoutput[which(anovaoutput$interaction<0.05),]

temp<-otu_table(datITS2rcsoilorder2)
colnames(temp)<-tax_table(datITS2rcsoilorder2)[,"Rank4"]
datITS2rcsoilorder3<-data.frame(temp,sample_data(datITS2rcsoilorder2))

means<-datITS2rcsoilorder3%>%
  group_by(Site,Transect)%>%
  summarize(mean=mean(o__Cystobasidiales),std=std.error(o__Cystobasidiales))
means$Site<-factor(means$Site,levels=c("Barataria","Turtle Cove","Pearl River","Fontainebleau","Big Branch","Bayou Sauvage","LUMCON 1","LUMCON 2"))

ggplot(means,aes(x=Transect,y=mean,group=Site,color=Transect))+
  labs(x = "",y="Relative abundance")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5))+
  geom_point(size=3)+
  geom_line(stat = "identity", position = "identity",size=.8,color="black")+
  geom_errorbar(aes(ymax = mean+std, ymin=mean-std),width=.15,size=.8) +
  facet_wrap(~Site,nrow=1)#,scales="free"

means<-datITS2rcsoilorder3%>%
  group_by(Transect)%>%
  summarize(mean=mean(o__Cystobasidiales),std=std.error(o__Cystobasidiales))

ggplot(means,aes(x=Transect,y=mean,group=Transect,color=Transect))+
  labs(x = "",y="Relative abundance")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5))+
  geom_point(size=3)+
  geom_line(stat = "identity", position = "identity",size=.8,color="black")+
  geom_errorbar(aes(ymax = mean+std, ymin=mean-std),width=.15,size=.8) 



####### Looking at composition at the order level #####

datITS2rcsoilorder <- tax_glom(datITS2rcsoil, taxrank="Rank4")
tax_table(datITS2rcsoilorder)
ntaxa(datITS2rcsoilorder)

tempphy<-datITS2rcsoilorder%>%
  #subset_samples(Site==sites[i])%>%
  filter_taxa(function(x) sum(x>0) >2, prune=T)
mynmds <- ordinate(tempphy, "CAP", "bray",formula=as.formula(~Transect+Condition(Site)))
mynmds <- ordinate(tempphy, "CAP", "bray",formula=as.formula(~Transect+SalinityClass))
plot_ordination(tempphy, mynmds, type="samples", color="Transect")+
  theme_classic()+
  #theme(legend.position = "none")+
  stat_ellipse(geom = "polygon", type="t",alpha=0.0,level=.95)# , aes(fill=Site),

anova(mynmds,by = "margin",permutations = how(nperm=1000))


datITS2rcsoilorder2 <- tax_glom(datITS2rcsoil, taxrank="Rank4")%>%
  filter_taxa(function(x) sum(x>0) >10, prune=T)

options(contrasts=c("contr.helmert","contr.poly"));options("contrasts")

anovaoutput<-data.frame(taxon=rep(NA,ntaxa(datITS2rcsoilorder2)),site=rep(NA,ntaxa(datITS2rcsoilorder2)),transect=rep(NA,ntaxa(datITS2rcsoilorder2)),interaction=rep(NA,ntaxa(datITS2rcsoilorder2)))
for(i in 1:ntaxa(datITS2rcsoilorder2)){
  temptaxa<-tax_table(datITS2rcsoilorder2)[i,"Rank4"]
  tempdata<-data.frame(otu_table(datITS2rcsoilorder2)[,i],sample_data(datITS2rcsoilorder2))
  colnames(tempdata)[1]<-"OTU"
  m1<-gls(log(OTU+1)~Site*Transect,data=tempdata)
  anovaoutput[i,1]<-temptaxa
  anovaoutput[i,2]<-anova(m1,type="marginal")[2,"p-value"]
  anovaoutput[i,3]<-anova(m1,type="marginal")[3,"p-value"]
  anovaoutput[i,4]<-anova(m1,type="marginal")[4,"p-value"]
  hist(resid(m1))
}
anovaoutput$siteq<-p.adjust(anovaoutput$site,method="fdr")
anovaoutput$transectq<-p.adjust(anovaoutput$transect,method="fdr")
anovaoutput$interactionq<-p.adjust(anovaoutput$interaction,method="fdr")

anovaoutput[which(anovaoutput$transectq<0.05),]
anovaoutput[which(anovaoutput$interactionq<0.05),]

temp<-otu_table(datITS2rcsoilorder2)
colnames(temp)<-tax_table(datITS2rcsoilorder2)[,"Rank4"]
datITS2rcsoilorder3<-data.frame(temp,sample_data(datITS2rcsoilorder2))

means<-datITS2rcsoilorder3%>%
  group_by(SalinityClass,Transect)%>%
  summarize(mean=mean(o__Erythrobasidiales),std=std.error(o__Erythrobasidiales))
means$Site<-factor(means$Site,levels=c("Barataria","Turtle Cove","Pearl River","Fontainebleau","Big Branch","Bayou Sauvage","LUMCON 1","LUMCON 2"))

ggplot(means,aes(x=Transect,y=mean,group=SalinityClass,color=Transect))+
  labs(x = "",y="Relative abundance")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5))+
  geom_point(size=3)+
  geom_line(stat = "identity", position = "identity",size=.8,color="black")+
  geom_errorbar(aes(ymax = mean+std, ymin=mean-std),width=.15,size=.8) +
  facet_wrap(~SalinityClass,nrow=1)#,scales="free"

means<-datITS2rcsoilorder3%>%
  group_by(Transect)%>%
  summarize(mean=mean(o__Mytilinidiales),std=std.error(o__Mytilinidiales))

ggplot(means,aes(x=Transect,y=mean,group=Transect,color=Transect))+
  labs(x = "",y="Relative abundance")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5))+
  geom_point(size=3)+
  geom_line(stat = "identity", position = "identity",size=.8,color="black")+
  geom_errorbar(aes(ymax = mean+std, ymin=mean-std),width=.15,size=.8) 







####### Looking at composition at the family level #####

datITS2rcsoilfamily <- tax_glom(datITS2rcsoil, taxrank="Rank5")
tax_table(datITS2rcsoilfamily)

tempphy<-datITS2rcsoilfamily%>%
  #subset_samples(Site==sites[i])%>%
  filter_taxa(function(x) sum(x>0) >2, prune=T)
mynmds <- ordinate(tempphy, "CAP", "bray",formula=as.formula(~Transect+Condition(Site)))
mynmds <- ordinate(tempphy, "CAP", "bray",formula=as.formula(~Transect+Salinitylog))
plot_ordination(tempphy, mynmds, type="samples", color="Transect")+
  theme_classic()+
  #theme(legend.position = "none")+
  stat_ellipse(geom = "polygon", type="t",alpha=0.0,level=.95)# , aes(fill=Site),

anova(mynmds,by = "margin",permutations = how(nperm=1000))


datITS2rcsoilfamily2 <- tax_glom(datITS2rcsoil, taxrank="Rank5")%>%
  filter_taxa(function(x) sum(x>0) >10, prune=T)


options(contrasts=c("contr.helmert","contr.poly"));options("contrasts")

anovaoutput<-data.frame(taxon=rep(NA,ntaxa(datITS2rcsoilfamily2)),site=rep(NA,ntaxa(datITS2rcsoilfamily2)),transect=rep(NA,ntaxa(datITS2rcsoilfamily2)),interaction=rep(NA,ntaxa(datITS2rcsoilfamily2)))
for(i in 1:ntaxa(datITS2rcsoilfamily2)){
  temptaxa<-tax_table(datITS2rcsoilfamily2)[i,"Rank5"]
  tempdata<-data.frame(otu_table(datITS2rcsoilfamily2)[,i],sample_data(datITS2rcsoilorder2))
  colnames(tempdata)[1]<-"OTU"
  m1<-gls(log(OTU+1)~SalinityClass*Transect,data=tempdata)
  anovaoutput[i,1]<-temptaxa
  anovaoutput[i,2]<-anova(m1,type="marginal")[2,"p-value"]
  anovaoutput[i,3]<-anova(m1,type="marginal")[3,"p-value"]
  anovaoutput[i,4]<-anova(m1,type="marginal")[4,"p-value"]
  hist(resid(m1))
}
anovaoutput$siteq<-p.adjust(anovaoutput$site,method="fdr")
anovaoutput$transectq<-p.adjust(anovaoutput$transect,method="fdr")
anovaoutput$interactionq<-p.adjust(anovaoutput$interaction,method="fdr")

anovaoutput[which(anovaoutput$transectq<0.05),]
anovaoutput[which(anovaoutput$interactionq<0.05),]

temp<-otu_table(datITS2rcsoilfamily2)
colnames(temp)<-tax_table(datITS2rcsoilfamily2)[,"Rank5"]
datITS2rcsoilfamily3<-data.frame(temp,sample_data(datITS2rcsoilfamily2))

means<-datITS2rcsoilfamily3%>%
  group_by(Site,Transect)%>%
  summarize(mean=mean(f__Glomeraceae),std=std.error(f__Glomeraceae))
means$Site<-factor(means$Site,levels=c("Barataria","Turtle Cove","Pearl River","Fontainebleau","Big Branch","Bayou Sauvage","LUMCON 1","LUMCON 2"))

ggplot(means,aes(x=Transect,y=mean,group=Site,color=Transect))+
  labs(x = "",y="Relative abundance")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5))+
  geom_point(size=3)+
  geom_line(stat = "identity", position = "identity",size=.8,color="black")+
  geom_errorbar(aes(ymax = mean+std, ymin=mean-std),width=.15,size=.8) +
  facet_wrap(~Site,nrow=1)#,scales="free"

means<-datITS2rcsoilorder3%>%
  group_by(Transect)%>%
  summarize(mean=mean(o__Mytilinidiales),std=std.error(o__Mytilinidiales))

ggplot(means,aes(x=Transect,y=mean,group=Transect,color=Transect))+
  labs(x = "",y="Relative abundance")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5))+
  geom_point(size=3)+
  geom_line(stat = "identity", position = "identity",size=.8,color="black")+
  geom_errorbar(aes(ymax = mean+std, ymin=mean-std),width=.15,size=.8) 








####### BACTERIA - Looking at composition at the phylum level #####

datBac3rcphylum <- tax_glom(datBac3rc, taxrank="Rank2")%>%
  subset_taxa(Rank2!="p__")
unique(tax_table(datBac3rcphylum))[,"Rank2"]

tempphy<-datBac3rcphylum%>%
  #subset_samples(Site==sites[i])%>%
  filter_taxa(function(x) sum(x>0) >2, prune=T)
mynmds <- ordinate(tempphy, "CAP", "bray",formula=as.formula(~Transect+Condition(Site)))
mynmds <- ordinate(tempphy, "CAP", "bray",formula=as.formula(~Transect+SalinityClass))
plot_ordination(tempphy, mynmds, type="samples", color="Transect")+
  theme_classic()+
  #theme(legend.position = "none")+
  stat_ellipse(geom = "polygon", type="t",alpha=0.0,level=.95)# , aes(fill=Site),

anova(mynmds,by = "margin",permutations = how(nperm=1000))


datBac3rcphylum2 <- tax_glom(datBac3rc, taxrank="Rank2")%>%
  filter_taxa(function(x) sum(x>0) >20, prune=T)


#options(contrasts=c("contr.treatment","contr.poly"));options("contrasts")
options(contrasts=c("contr.helmert","contr.poly"));options("contrasts")
#notes: when you do poisson or negbin regression use Chi2 stats test, when you use quasipoisson use F test. the glm.nb and Anova were giving errors in fitting  b/c there are site*transect combinations that are all zero, the Anova often gave a warning about convergence. even if models would fit, nbinom was making things really significant when there were site*transect combinations that were all zero. The residual plots on the link scale are weird when there are sites*transect combinations that are all zero. taking out the site*transect combinations that are all zero typically fixes the zero inflated issue. and I cant figure out how to do an anova with the zero inflated model fit. I can't take out site*transect combinations and then fit a site*transect interactions b/c then the model is singular, so I have to take out hte whole site for any instance when site*transect combo is zero. I can still use the quasipoisson or normal when i have site*transect combos that are all 0

#m2 <- zeroinfl(OTU ~ Site*Transect | Site*Transect,data = tempdata, dist = "negbin")#, EM = TRUE
#vuong(m1, m2)

anovaoutput<-data.frame(taxon=rep(NA,ntaxa(datBac3rcphylum2)),site=rep(NA,ntaxa(datBac3rcphylum2)),transect=rep(NA,ntaxa(datBac3rcphylum2)),interaction=rep(NA,ntaxa(datBac3rcphylum2)))
for(i in 1:ntaxa(datBac3rcphylum2)){
  temptaxa<-tax_table(datBac3rcphylum2)[i,"Rank2"]
  tempdata<-data.frame(t(otu_table(datBac3rcphylum2))[,i],sample_data(datBac3rcphylum2))
  colnames(tempdata)[1]<-"OTU"
  #temp1<-aggregate.data.frame(tempdata$OTU,by=list(tempdata$Site,tempdata$Transect),sum);temp1
  #temp2<-temp1$Group.1[which(temp1$x==0)]
  #tempdata<-tempdata%>%filter(!Site%in%temp2)
  #m1<-gls(log(OTU+1)~Site*Transect,data=tempdata)
  #anova(m1,type="marginal")
  #m1<-gls(log(OTU+1)~Site*Transect,correlation=corSpher(form = ~ Lat+ Long),data=tempdata)
  m1<-glm(OTU~Site*Transect,family=quasipoisson,data=tempdata)
  #Anova(m1,type="III",test.statistic ="F")
  #m1<-glm.nb(OTU~Site*Transect, maxit = 1000,data=tempdata)
  #Anova(m1,type="III")
  #summary(m1)
  #plot(predict(m1,type="response"),resid(m1,type="deviance"))
  anovaoutput[i,1]<-temptaxa
  #anovaoutput[i,2]<-anova(m1,type="marginal")[2,"p-value"]
  #anovaoutput[i,3]<-anova(m1,type="marginal")[3,"p-value"]
  #anovaoutput[i,4]<-anova(m1,type="marginal")[4,"p-value"]
  anovaoutput[i,2]<-Anova(m1,type="III", test.statistic = "F")[1,"Pr(>F)"]
  anovaoutput[i,3]<-Anova(m1,type="III", test.statistic = "F")[2,"Pr(>F)"]
  anovaoutput[i,4]<-Anova(m1,type="III", test.statistic = "F")[3,"Pr(>F)"]
  #anovaoutput[i,2]<-Anova(m1,type="III")[1,"Pr(>Chisq)"]
  #anovaoutput[i,3]<-Anova(m1,type="III")[2,"Pr(>Chisq)"]
  #anovaoutput[i,4]<-Anova(m1,type="III")[3,"Pr(>Chisq)"]
  #hist(resid(m1))
}
anovaoutput$siteq<-p.adjust(anovaoutput$site,method="fdr")
anovaoutput$transectq<-p.adjust(anovaoutput$transect,method="fdr")
anovaoutput$interactionq<-p.adjust(anovaoutput$interaction,method="fdr")

#anovaoutput[which(anovaoutput$siteq<0.05),]
anovaoutput[which(anovaoutput$transectq<0.05),]
anovaoutput[which(anovaoutput$interactionq<0.05),]

temp<-t(otu_table(datBac3rcphylum2))
colnames(temp)<-tax_table(datBac3rcphylum2)[,"Rank2"]
datBac3rcphylum3<-data.frame(temp,sample_data(datBac3rcphylum2))


means<-datBac3rcphylum3%>%
  group_by(Site,Transect)%>%
  summarize(mean=mean(p__Euryarchaeota),std=std.error(p__Euryarchaeota))

ggplot(means,aes(x=Transect,y=mean,group=Site,color=Transect))+
  labs(x = "",y="Relative abundance")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5))+
  geom_point(size=3)+
  geom_line(stat = "identity", position = "identity",size=.8,color="black")+
  geom_errorbar(aes(ymax = mean+std, ymin=mean-std),width=.15,size=.8) +
  facet_wrap(~Site,nrow=1)#,scales="free"


means<-datBac3rcphylum3%>%
  group_by(SalinityClass,Transect)%>%
  summarize(mean=mean(p__Crenarchaeota),std=std.error(p__Crenarchaeota))

ggplot(means,aes(x=Transect,y=mean,group=SalinityClass,color=Transect))+
  labs(x = "",y="Relative abundance")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5))+
  geom_point(size=3)+
  geom_line(stat = "identity", position = "identity",size=.8,color="black")+
  geom_errorbar(aes(ymax = mean+std, ymin=mean-std),width=.15,size=.8) +  
  facet_wrap(~SalinityClass,nrow=1)#,scales="free"


means<-datBac3rcphylum3%>%
  group_by(Transect)%>%
  summarize(mean=mean(p__Acidobacteria),std=std.error(p__Acidobacteria))

ggplot(means,aes(x=Transect,y=mean,group=Transect,color=Transect))+
  labs(x = "",y="Relative abundance")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5))+
  geom_point(size=3)+
  geom_line(stat = "identity", position = "identity",size=.8,color="black")+
  geom_errorbar(aes(ymax = mean+std, ymin=mean-std),width=.15,size=.8) 
