# Analyse the effects of climatic variables on body size (mist-net dataset)
# and call frequency (acoustic dataset) 
# Sensitivity analyses
# Plot results (fig. 1 and 2)

# R script Caterina Penone, University of Bern
# caterina.penone [at] gmail.com
# Under Git repository - pipip_frequency_bodysz

# Load required libraries
library(lme4)
library(lmerTest)
library(ggplot2)

# Set working directory
setwd("path")

#### Mist-net dataset analyses ####

# Read data
capt<-read.table("/data_pipip_capture.txt",h=T)
capt$Year<-as.factor(capt$Year)
#remove Immature individuals
capt<-capt[!"Age" %in% "Immature"]
#remove 2 outliers
capt<-subset(capt,Forearm_length_mm<40)
names(capt)[names(capt) %in% c("Worldclim1","Worldclim7","Worldclim10",
                               "Worldclim11","Worldclim12","Worldclim16")]<-c("Bio1","Bio7","Bio10","Bio11","Bio12","Bio16")

# PCA of environmental variables from Worldclim
PCAcapt<-prcomp(capt[,c("X","Y","Bio1","Bio7","Bio10","Bio11","Bio12","Bio16")],center=T,scale=T)
capt$Prin_comp_1<-as.data.frame(PCAcapt$x)[,1]
capt$Prin_comp_2<-as.data.frame(PCAcapt$x)[,2]
capt$Prin_comp_3<-as.data.frame(PCAcapt$x)[,3]
summary(PCAcapt)

# Extract loadings
#write.csv(PCAcapt$rotation,"PCAloadCAPT.csv")

# Model
lmeCapt<-lmer(Forearm_length_mm~Sex+day_of_year+Prin_comp_1+Prin_comp_2 +(1|Year/Locality_id),na.action=na.omit, data=capt)
summary(lmeCapt)

# check model assumptions
plot(lmeCapt) #heteroscedasticity
hist(residuals(lmeCapt)) #normality of residuals
qqnorm(resid(lmeCapt))
qqline(resid(lmeCapt))

# Sensitivity analysis: subsampling of the dense areas
captW<-subset(capt,(X<349356 & Y>2155220))
lmeCaptW<-lmer(Forearm_length_mm~Sex+day_of_year+Prin_comp_1+Prin_comp_2+(1|Year/Locality_id),na.action=na.omit, data=captW)
summary(lmeCaptW)
plot(lmeCaptW) #heteroscedasticity
hist(residuals(lmeCaptW)) #normality of residuals
qqnorm(resid(lmeCaptW))
qqline(resid(lmeCaptW))



#### Acoustic dataset analyses ####

# Read and prepare data
sound<-read.csv('/data_pipip_sound.csv',h=T)
# Pipip freq is between 41 and 52 KHz
sound<-subset(sound, LowFc<52)
# Remove island
sound<-subset(sound, circuit_id!=48 & circuit_id!=49)
# Filter out "extreme" slope QCF
sound<-subset(sound, -5<SlopeQCF)
sound<-subset(sound, SlopeQCF<1)
# Remove 3rd survey (very few)
sound<-subset(sound, survey %in% c(1,2))
# Restruc to altitudetude<500 m
sound<-subset(sound, altitude<500)
# Tranform into factor and remove missing data
sound$circuit_id<-as.factor(sound$circuit_id)
sound$segment_id<-as.factor(sound$segment_id)
sound$survey<-as.factor(sound$survey)
#remove NA's in dataset
sound<-sound[!is.na(sound$Worldclim1),]

# Rename in order to be the same as mist-net data
names(sound)[names(sound) %in% c("Worldclim10","Worldclim11","Worldclim12","Worldclim1",
                               "Worldclim7","Worldclim16")]<-c("Bio10","Bio11","Bio12","Bio1","Bio7","Bio16")

# PCA of environmental variables from Worldclim
PCAsound<-prcomp(sound[,c("X","Y","Bio1","Bio7","Bio10","Bio11","Bio12","Bio16")],center=T,scale=T)
sound$Prin_comp_1<-as.data.frame(PCAsound$x)[,1]
sound$Prin_comp_2<-as.data.frame(PCAsound$x)[,2]
summary(PCAsound)


# Extract loadings
#write.csv(PCAsound$rotation,"PCAloadSOUND.csv")

# Model
lmeSound<-lmer(LowFc~SlopeQCF+survey+call_abundance+Prin_comp_1+Prin_comp_2+(1|circuit_id/segment_id),na.action=na.omit, data=sound)
summary(lmeSound)


# check model assumptions
plot(lmeSound) #heteroscedasticity
hist(residuals(lmeSound)) #normality of residuals
qqnorm(resid(lmeSound))
qqline(resid(lmeSound))

# analysis using only the June survey (survey 1)
lmeSound1<-lmer(LowFc~SlopeQCF+day_of_year+call_abundance+Prin_comp_1+Prin_comp_2+(1|circuit_id/segment_id),na.action=na.omit, data=subset(sound,survey==1))
summary(lmeSound1)
plot(lmeSound1) #heteroscedasticity
hist(residuals(lmeSound1)) #normality of residuals

# analysis using information on the day of the year (instead of survey)
lmeSound2<-lmer(LowFc~SlopeQCF+day_of_year+call_abundance+Prin_comp_1+Prin_comp_2+(1|circuit_id/segment_id),na.action=na.omit, data=subset(sound,survey==1))
summary(lmeSound2)
plot(lmeSound2) #heteroscedasticity
hist(residuals(lmeSound2)) #normality of residuals

# Sensitivity analysis: the results hold with/without local temp in the model?
sound2<-sound[!is.na(sound$temperature_degreeC),]
lmeSoundTEMP<-lmer(LowFc~temperature_degreeC+SlopeQCF+survey+call_abundance+Prin_comp_1+Prin_comp_2+(1|circuit_id/segment_id),na.action=na.omit, data=sound2)
lmeSoundNOTEMP<-lmer(LowFc~SlopeQCF+survey+call_abundance+Prin_comp_1+Prin_comp_2+(1|circuit_id/segment_id),na.action=na.omit, data=sound2)
anova(lmeSoundTEMP,lmeSoundNOTEMP) #without temperature it is better
summary(lmeSoundTEMP)
summary(lmeSoundNOTEMP)

# Sensitivity analysis: subsampling of the dense areas
plot(Y~X,data=sound)
abline(v=410000)
abline(h=2155220)
soundN<-subset(sound,(X>410000 & Y>2155220))

lmeSoundN<-lmer(LowFc~SlopeQCF+survey+call_abundance+Prin_comp_1+Prin_comp_2+ (1|circuit_id/segment_id),na.action=na.omit, data=soundN)
summary(lmeSoundN)


# Sensitivity analysis: potential confusion between P. pipistrellus and P. pygmaeus - analyse only South ####REMOVE?
soundS<-subset(sound,(Y<2155220))
lmeSoundS<-lmer(LowFc~SlopeQCF+survey+call_abundance+Prin_comp_1+Prin_comp_2+ (1|circuit_id/segment_id),na.action=na.omit, data=soundS)
summary(lmeSoundS)




#################################### FIGURES #############################################

#### Fig1: PCA plots #####
PCAsound<-prcomp(sound[,c("X","Y","Bio1","Bio7","Bio10","Bio11","Bio12","Bio16")],center=T,scale=T)
PCAcapt<-prcomp(capt[,c("X","Y","Bio1","Bio7","Bio10","Bio11","Bio12","Bio16")],center=T,scale=T)

##biplots
pdf("PCAs.pdf", height = 7.5, width = 14)
par(mfrow=c(1,2),mar = c(1,5,1,2))
biplot(PCAcapt,scale=T,xlabs=rep("o", nrow(capt)),expand=0.9,xlab="M.PC1: 58%",ylab="M.PC2: 29%",cex.lab=1.5)
text(-40, 42, "a",cex=1.5)
biplot(PCAsound,scale=T,xlabs=rep("o", nrow(sound)),expand=1.2,xlim=c(-0.02,0.032),xlab="A.PC1: 59%",ylab="A.PC2: 34%",cex.lab=1.5)
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
p1<-ggplot(capt, aes(x = Prin_comp_1, y = Forearm_length_mm)) +
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
#ggsave("bd_sz.pdf", width=16, height=16, units="cm")

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
#ggsave("freq.pdf", width=16, height=16, units="cm")


# Plot the two together
pdf("Plots_pipgeo.pdf", height = 7.5, width = 14)
multiplot(p1, p2, cols=2)
dev.off()