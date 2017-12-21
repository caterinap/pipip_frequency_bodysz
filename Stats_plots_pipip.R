# Analyse the effects of climatic variables on body size (mist-net dataset)
# and call frequency (acoustic dataset) 
# Sensitivity analyses
# Plot results (fig. 1 and 2)

# R script Caterina Penone, University of Bern
# caterina.penone [at] gmail.com
# Under Git repository - pipip_frequency_bodysz

# Load required libraries
library(nlme)
library(ggplot2)


#### Mist-net dataset analyses ####

# Read data
capt<-read.table("./data_pipip_capture.txt",h=T)
capt$ANNEE<-as.factor(capt$ANNEE)
capt<-capt[!"Age" %in% "Imm"]
#remove 2 outliers
capt<-subset(capt,Ab<40)
names(capt)[14]<-"Bio7"

# PCA of environmental variables from Worldclim
PCAcapt<-prcomp(capt  [,c("X","Y","Bio1","Bio7","Bio10","Bio11","Bio12","Bio16")],center=T,scale=T)
capt$Prin_comp_1<-as.data.frame(PCAcapt$x)[,1]
capt$Prin_comp_2<-as.data.frame(PCAcapt$x)[,2]
capt$Prin_comp_3<-as.data.frame(PCAcapt$x)[,3]
summary(PCAcapt,scale=T)

# Extract loadings
write.csv(PCAcapt$rotation,"PCAloadCAPT.csv")

# Model
lmeCapt<-lme(Ab~Sexe+Prin_comp_1+Prin_comp_2,random=~ 1|CLE_Comm/ANNEE,na.action=na.omit, data=capt)
lmeCaptRE<-lm(Ab~Sexe+Prin_comp_1+Prin_comp_2,na.action=na.omit, data=capt)
anova(lmeCapt,lmeCaptRE) #ok, with random effect model is better
summary(lmeCapt)
anova(lmeCapt)

# Sensitivity analysis: analyse separately East and west
captW<-subset(capt,X<349356)
lmeCaptW<-lme(Ab~Sexe+Prin_comp_1+Prin_comp_2,random=~ 1|CLE_Comm/ANNEE,na.action=na.omit, data=captW)
summary(lmeCaptW)
anova(lmeCaptW)

captE<-subset(capt,X>349356)
lmeCaptE<-lme(Ab~Sexe+Prin_comp_1+Prin_comp_2,random=~ 1|CLE_Comm/ANNEE,na.action=na.omit, data=captE)
summary(lmeCaptE)
anova(lmeCaptE)



#### Acoustic dataset analyses ####

# Read and prepare data
sound<-read.csv('./data_pipip_sound.csv',h=T)
# Pipip freq is between 41 and 52 KHz
sound<-subset(sound, LowFc<52)
# Remove island
sound<-subset(sound, id_cir!=48 & id_cir!=49)
# Filter out "extreme" slope QCF
sound<-subset(sound, -5<SlopeQCF)
sound<-subset(sound,SlopeQCF<1)
# Remove 3rd survey (very few)
sound<-subset(sound, passage!=3 & passage!=4)
# Restruc to altitude<500 m
sound<-subset(sound, alti<500)
# Remove other species  
sound<-subset(sound,rf.pred!='Pipnat')
sound<-subset(sound,rf.pred!='Pippyg')
# Tranform into factor and remove missing data
sound$id_cir<-as.factor(sound$id_cir)
sound$id_tron<-as.factor(sound$id_tron)
sound$passage<-as.factor(sound$passage)
sound<-sound[!is.na(sound$b1MN),]

# Rename in order to be the same as mist-net data
names(sound)[c(128,137:144)]<-c("altitude","Bio10","Bio11","Bio12","Bio1","Bio7","Bio15","Bio16","Bio18")

# PCA of environmental variables from Worldclim
PCAsound<-prcomp(sound[,c("X","Y","Bio1","Bio7","Bio10","Bio11","Bio12","Bio16")],center=T,scale=T,na.rm=T)
sound$Prin_comp_1<-as.data.frame(PCAsound$x)[,1]
sound$Prin_comp_2<-as.data.frame(PCAsound$x)[,2]
summary(PCAsound,scale=T)

# Extract loadings
write.csv(PCAsound$rotation,"PCAloadSOUND.csv")

# Model
lmeSound<-lme(LowFc~SlopeQCF+passage+Prin_comp_1+Prin_comp_2,random=~ 1|id_cir,na.action=na.omit, data=sound)
lmeSoundRE<-lm(LowFc~SlopeQCF+jour+Prin_comp_1+Prin_comp_2,na.action=na.omit, data=sound)
anova(lmeSound,lmeSoundRE) #ok, with random effect model is better
summary(lmeSound)
anova(lmeSound)


# Sensitivity analysis: the results hold with/without local temp in the model?
sound2<-sound[!is.na(sound$temp),]
lmeSoundTEMP<-lme(LowFc~SlopeQCF+temp+passage+Prin_comp_1+Prin_comp_2,random=~ 1|id_cir,na.action=na.omit, data=sound2)
lmeSoundNOTEMP<-lme(LowFc~SlopeQCF+passage+Prin_comp_1+Prin_comp_2,random=~ 1|id_cir,na.action=na.omit, data=sound2)
anova(lmeSoundTEMP,lmeSoundNOTEMP) #with temperature it is better
summary(lmeSoundTEMP)
summary(lmeSoundNOTEMP)

# Sensitivity analysis: the results hold with/without SlopeQCF?
lmeSoundSlope<-lme(LowFc~passage+Prin_comp_1+Prin_comp_2,random=~ 1|id_cir,na.action=na.omit, data=sound)
anova(lmeSoundSlope)
summary(lmeSoundSlope)


# Sensitivity analysis: the results hold when south circuits are removed from the analysis?
soundN<-subset(sound,Y>2155220)
lmeSoundN<-lme(LowFc~SlopeQCF+passage+Prin_comp_1+Prin_comp_2,random=~ 1|id_cir,na.action=na.omit, data=soundN)
summary(lmeSoundN)
anova(lmeSoundN)

soundS<-subset(sound,Y<2155220)
lmeSoundS<-lme(LowFc~SlopeQCF+passage+Prin_comp_1+Prin_comp_2,random=~ 1|id_cir,na.action=na.omit, data=soundS)
summary(lmeSoundS)
anova(lmeSoundS)



#################################### FIGURES #############################################

#### Fig1: PCA plots #####
PCAsound<-prcomp(sound[,c("X","Y","Bio1","Bio7","Bio10","Bio11","Bio12","Bio16")],center=T,scale=T,na.rm=T)
PCAcapt<-prcomp(capt[,c("X","Y","Bio1","Bio7","Bio10","Bio11","Bio12","Bio16")],center=T,scale=T)

##biplots
pdf("PCAs.pdf", height = 7.5, width = 14)
par(mfrow=c(1,2),mar = c(1,5,1,2))
biplot(PCAcapt,scale=T,xlabs=rep("o", nrow(capt)),expand=0.9,xlab="PC1: 58%",ylab="PC2: 29%",cex.lab=1.5)
text(-40, 42, "a",cex=1.5)
biplot(PCAsound,scale=T,xlabs=rep("o", nrow(sound)),expand=1.2,xlim=c(-0.02,0.032),xlab="PC1: 59%",ylab="PC2: 34%",cex.lab=1.5)
text(-65, 175, "b",cex=1.5)
dev.off()


##### Fig2: Scatterplots #####

# Function to create multiplots with ggplot2 
# (source: http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/)
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                    layout.pos.col = matchidx$col))
      }
  }
}

# Plot of capture dataset
p1<-ggplot(capt, aes(x = Prin_comp_1, y = Ab)) +
  theme_bw(base_size=20) +
  xlab("PC1")  +      ylab("Forearm length (mm)")  +
  #ggtitle("Capture dataset") +
  geom_point() +
  stat_smooth(method = "lm", col = "black")+
  theme(axis.title.x=element_text(vjust=-1.5)) +
  theme(axis.title.y=element_text(vjust=1.5)) +
  theme(plot.title=element_text(vjust=1.5)) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  theme(panel.grid.major=element_blank())+
  geom_text(data = NULL, x = -7.2, y = 37.5, label = "a",size=10)+
  scale_y_continuous(breaks=seq(28.5, 37.5, 2))+
  scale_x_continuous(breaks=seq(-7, 4, 2))
p1
ggsave("bd_sz.pdf", width=16, height=16, units="cm")

# Plot of acoustic dataset
p2<- ggplot(sound, aes(x = Prin_comp_2, y = LowFc)) + 
  theme_bw(base_size=20) +
  xlab("PC2")  +      ylab("Characteristic frequency (MHz)")  +
  #ggtitle("Acoustic dataset")+
  geom_point() + 
  stat_smooth(method = glm, col = "black") +
  stat_smooth(method = "lm", col = "black")+
  theme(axis.title.x=element_text(vjust=-1.5)) +
  theme(axis.title.y=element_text(vjust=1.5)) +
  theme(plot.title=element_text(vjust=1.5)) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))+
  theme(panel.grid.major=element_blank())+
  geom_text(data = NULL, x = -2.7, y = 52.1, label = "b",size=10)+
  scale_y_continuous(breaks=seq(42, 52.5, 3))
p2
ggsave("freq.pdf", width=16, height=16, units="cm")


# Plot the two together
png("Plots_pipgeo.png", width=320, height=160, units="cm", res=600)
multiplot(p1, p2, cols=2)
dev.off()