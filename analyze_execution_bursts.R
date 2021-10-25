library("lme4")
library("car")
library("doBy")
library(tidyverse)
library(emmeans)
library(Rmisc)
library(cowplot)
source('R_rainclouds.R')

infant_data<-read.csv('../../data/UMD/12m/deriv/nbursts.csv')
infant_data$Age<-'12m'
adult_data<-read.csv('../../data/UMD/adult/deriv/nbursts.csv')
adult_data$Subject<-paste0(adult_data$Subject,'A')
adult_data$Age<-'adult'
data<-rbind(infant_data,adult_data)
data<-data[data$Condition=='execute',]

clusters<-unique(data$Cluster)
ages<-unique(data$Age)

# Run analysis of burst count for each age and cluster
for(age in ages) {
  for(cluster in clusters) {
    c_data<-data[data$Cluster==cluster & data$Age==age,]
    
    print(cluster)
    model <- glmer(Count ~ Epoch+(1|Subject), data = c_data, family=poisson(link='log'), control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
    results<-Anova(model,type=3)
    print(results)
    
    print(paste0(age,': baseline M=',mean(c_data$Count[c_data$Epoch=='baseline']),' SD=',sd(c_data$Count[c_data$Epoch=='baseline'])))
    print(paste0(age,': exp M=',mean(c_data$Count[c_data$Epoch=='exp']),' SD=',sd(c_data$Count[c_data$Epoch=='exp'])))
  }
}

sum_dat<-summarySE(data=data,'Count',groupvars=c('Subject','Age','Cluster','Epoch'),na.rm=TRUE)
summary_data<-summarySE(data =sum_dat, 'Count', groupvars=c('Age', 'Cluster','Epoch'), na.rm=TRUE)
p <- ggplot(sum_dat, aes_string(x = 'Epoch', y = 'Count', fill = 'Epoch')) +
  geom_flat_violin(data=data, position = position_nudge(x = .1, y = 0), adjust = 2, trim = FALSE, alpha = 1, colour = NA)+
  geom_point(data = data, aes(x = as.numeric(Epoch)-.15, color=Subject),position = position_jitter(width = .05), size = 2, shape = 20, show.legend = FALSE)+
  geom_boxplot(outlier.shape = NA, alpha = 1, width = .1, colour = "black")+
  geom_line(data = summary_data, aes(x = as.numeric(Epoch)+.075), linetype = 3)+
  geom_point(data = summary_data, aes(x = as.numeric(Epoch)+.075), shape = 18) +
  geom_errorbar(data = summary_data, aes(x = as.numeric(Epoch)+.075, colour = Epoch, ymin = Count-se, ymax = Count+se), width = .025)+
  facet_grid(Age~Cluster)+
  theme_cowplot()
print(p)
ggsave('../../manuscript/figures/count.eps')
ggsave('../../manuscript/figures/count.png')


infant_data<-read.csv('../../data/UMD/12m/deriv/bursts.csv')
infant_data$Age<-'12m'
# Compute duration
infant_data$Duration<-infant_data$Offset-infant_data$Onset
# Convert duration to cycles
infant_data$DurationCycles<-0
mean_betas=c(15.5,16)
infant_data$DurationCycles[infant_data$Cluster=='C3']<-infant_data$Duration[infant_data$Cluster=='C3']/(1000/mean_betas[1])
infant_data$DurationCycles[infant_data$Cluster=='C4']<-infant_data$Duration[infant_data$Cluster=='C4']/(1000/mean_betas[2])

adult_data<-read.csv('../../data/UMD/adult/deriv/bursts.csv')
adult_data$Subject<-paste0(adult_data$Subject,'A')
adult_data$Age<-'adult'
# Compute duration
adult_data$Duration<-adult_data$Offset-adult_data$Onset
# Convert duration to cycles
adult_data$DurationCycles<-0
mean_betas=c(22.5,21.5)
adult_data$DurationCycles[adult_data$Cluster=='C3']<-adult_data$Duration[adult_data$Cluster=='C3']/(1000/mean_betas[1])
adult_data$DurationCycles[adult_data$Cluster=='C4']<-adult_data$Duration[adult_data$Cluster=='C4']/(1000/mean_betas[2])

data<-rbind(infant_data,adult_data)

data<-data[data$Condition=='execute',]
data$Age<-as.factor(data$Age)

clusters<-unique(data$Cluster)

# Analyze duration in each cluster
for(cluster in clusters) {
  cdata<-data[data$Cluster==cluster,]
  
  print(cluster)
  model <- lmer(Duration ~ Age+(1|Subject), data = cdata, control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
  results<-Anova(model,type=3)
  print(results)
  
  print(paste0('adult M=',mean(cdata$Duration[cdata$Age=='adult'],na.rm=TRUE),' SD=',sd(cdata$Duration[cdata$Age=='adult'],na.rm=TRUE)))
  print(paste0('12m M=',mean(cdata$Duration[cdata$Age=='12m'],na.rm=TRUE),' SD=',sd(cdata$Duration[cdata$Age=='12m'],na.rm=TRUE)))
}

# Analyze duration (cycles) in each cluster
for(cluster in clusters) {
  cdata<-data[data$Cluster==cluster,]
  
  print(cluster)
  model <- lmer(DurationCycles ~ Age+(1|Subject), data = cdata, control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
  results<-Anova(model,type=3)
  print(results)
  
  print(paste0('adult M=',mean(cdata$DurationCycles[cdata$Age=='adult'],na.rm=TRUE),' SD=',sd(cdata$DurationCycles[cdata$Age=='adult'],na.rm=TRUE)))
  print(paste0('12m M=',mean(cdata$DurationCycles[cdata$Age=='12m'],na.rm=TRUE),' SD=',sd(cdata$DurationCycles[cdata$Age=='12m'],na.rm=TRUE)))
}

# Analyze amplitude in eahc cluster
for(cluster in clusters) {
  cdata<-data[data$Cluster==cluster,]
  print(cluster)
  model <- lmer(Amp ~ Age+(1|Subject), data = cdata, control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
  results<-Anova(model,type=3)
  print(results)

  print(paste0('adult M=',mean(cdata$Amp[cdata$Age=='adult'],na.rm=TRUE),' SD=',sd(cdata$Amp[cdata$Age=='adult'],na.rm=TRUE)))
  print(paste0('12m M=',mean(cdata$Amp[cdata$Age=='12m'],na.rm=TRUE),' SD=',sd(cdata$Amp[cdata$Age=='12m'],na.rm=TRUE)))
}

sum_dat<-summarySE(data=data,'Duration',groupvars=c('Subject','Age','Cluster'),na.rm=TRUE)
summary_data<-summarySE(data =sum_dat, 'Duration', groupvars=c('Age', 'Cluster'), na.rm=TRUE)
dev.new()
p <- ggplot(sum_dat, aes_string(x = 'Age', y = 'Duration', fill = 'Age')) +
  geom_flat_violin(data=data, position = position_nudge(x = .1, y = 0), adjust = 2, trim = FALSE, alpha = 1, colour = NA)+
  geom_point(data = data, aes(x = as.numeric(Age)-.15, color=Subject),position = position_jitter(width = .05), size = 2, shape = 20)+
  geom_boxplot(outlier.shape = NA, alpha = 1, width = .1, colour = "black")+
  geom_line(data = summary_data, aes(x = as.numeric(Age)+.075), linetype = 3)+
  geom_point(data = summary_data, aes(x = as.numeric(Age)+.075), shape = 18) +
  geom_errorbar(data = summary_data, aes(x = as.numeric(Age)+.075, colour = Age, ymin = Duration-se, ymax = Duration+se), width = .025)+
  facet_grid(~Cluster)+
  ylim(0,750)+
  theme_cowplot()
print(p)

sum_dat<-summarySE(data=data,'DurationCycles',groupvars=c('Subject','Age','Cluster'),na.rm=TRUE)
summary_data<-summarySE(data =sum_dat, 'DurationCycles', groupvars=c('Age', 'Cluster'), na.rm=TRUE)
dev.new()
p <- ggplot(sum_dat, aes_string(x = 'Age', y = 'DurationCycles', fill = 'Age')) +
  geom_flat_violin(data=data, position = position_nudge(x = .1, y = 0), adjust = 2, trim = FALSE, alpha = 1, colour = NA)+
  geom_point(data = data, aes(x = as.numeric(Age)-.15, color=Subject),position = position_jitter(width = .05), size = 2, shape = 20)+
  geom_boxplot(outlier.shape = NA, alpha = 1, width = .1, colour = "black")+
  geom_line(data = summary_data, aes(x = as.numeric(Age)+.075), linetype = 3)+
  geom_point(data = summary_data, aes(x = as.numeric(Age)+.075), shape = 18) +
  geom_errorbar(data = summary_data, aes(x = as.numeric(Age)+.075, colour = Age, ymin = DurationCycles-se, ymax = DurationCycles+se), width = .025)+
  facet_grid(~Cluster)+
  ylim(0,21)+
  theme_cowplot()
print(p)


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

