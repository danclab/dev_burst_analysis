library("lme4")
library("car")
library("doBy")
library(tidyverse)
library(emmeans)
library(Rmisc)
library(cowplot)
source('R_rainclouds.R')

infant_data<-read.csv('../../data/UMD/12m/deriv/bursts_adult_thresh.csv')
infant_data$Age<-'12m'

adult_data<-read.csv('../../data/UMD/adult/deriv/bursts.csv')
adult_data$Subject<-paste0(adult_data$Subject,'A')
adult_data$Age<-'adult'

data<-rbind(infant_data,adult_data)

data<-data[data$Condition=='execute',]
data$Age<-as.factor(data$Age)

clusters<-unique(data$Cluster)

# Analyze amplitude in each cluster
for(cluster in clusters) {
  cdata<-data[data$Cluster==cluster,]
  print(cluster)
  model <- lmer(Amp ~ Age+(1|Subject), data = cdata, control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
  results<-Anova(model,type=3)
  print(results)
  
  print(paste0('adult M=',mean(cdata$Amp[cdata$Age=='adult'],na.rm=TRUE),' SD=',sd(cdata$Amp[cdata$Age=='adult'],na.rm=TRUE)))
  print(paste0('12m M=',mean(cdata$Amp[cdata$Age=='12m'],na.rm=TRUE),' SD=',sd(cdata$Amp[cdata$Age=='12m'],na.rm=TRUE)))
}

sum_dat<-summarySE(data=data, 'Amp', groupvars=c('Subject','Age','Cluster'),na.rm=TRUE)
summary_data<-summarySE(data=sum_dat, 'Amp', groupvars=c('Age', 'Cluster'), na.rm=TRUE)
dev.new()
p <- ggplot(sum_dat, aes_string(x = 'Age', y = 'Amp', fill = 'Age')) +
  geom_flat_violin(data=data, position = position_nudge(x = .1, y = 0), adjust = 2, trim = FALSE, alpha = 1, colour = NA)+
  geom_point(data=data, aes(x = as.numeric(Age)-.15, color=Subject),position = position_jitter(width = .05), size = 2, shape = 20)+
  geom_boxplot(outlier.shape = NA, alpha = 1, width = .1, colour = "black")+
  geom_line(data = summary_data, aes(x = as.numeric(Age)+.075), linetype = 3)+
  geom_point(data = summary_data, aes(x = as.numeric(Age)+.075), shape = 18) +
  geom_errorbar(data = summary_data, aes(x = as.numeric(Age)+.075, colour = Age, ymin = Amp-se, ymax = Amp+se), width = .025)+
  facet_grid(~Cluster)+
  ylim(0,80)+
  theme_cowplot()
print(p)
