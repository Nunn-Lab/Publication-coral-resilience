#lipid analysis for M cap  
library(ggplot2)
setwd("~/Lipid Analyses")
lipid.dat<-read.csv('lipid analysis.csv', na.strings='na')
View(lipid.dat)

#remove one sample with na in Tolerance
lipid.dat2<-lipid.dat[!is.na(lipid.dat$Tolerance),]
lipid.dat2$Time<-as.factor(lipid.dat2$Time)

#just T1 and T2
lipid.dat2<-lipid.dat[!is.na(lipid.dat$Tolerance),]
lipid.dat3<-subset(lipid.dat2, Time<3)
lipid.dat3$Time<-as.factor(lipid.dat3$Time)
View(lipid.dat3)
#Susceptible p<0.01 = hex 762a83  ; border =hex 8857A8  **note that -log10 (.o1) would be >2
#Susceptible p<0.05 = hex af8dc3  ; border =hex 8857A8     **note that -log10 (.05) would be >1.301
#Susceptible p<0.10 = hex E7D4E8  ; border =hex 8857A8     **note that -log10 (0.10) would be > 1.0
#all other Susceptible no fill  or hex FFFFFF ;   ; border =hex 8857A8

########Greens####
#Resilient  p<0.01 = hex 1b7837 ; border =hex 95C38B   **note that -log10 (.o1) would be >2
#Resilient  p<0.05 = hex 7fbf7b ; border =hex 95C38B     **note that -log10 (.05) would be >1.301
#Resilient  p<0.10 = hex d9f0d3 ; border =hex 95C38B     **note that -log10 (0.10) would be > 1.0
#all other Resilient  no fill or hex FFFFFF;  ; border =hex 95C38B

lipid.dat3$Tolerance.Tim <-paste(lipid.dat3$Tolerance,lipid.dat3$Time)
lipid.dat3$Toler.Tim.Stat <-paste(lipid.dat3$Tolerance,lipid.dat3$Time, lipid.dat3$Status)

#I want RT1='#1b7837','ST1#762a83','Rt2=#7fbf7b','St2=#af8dc3' and NB samples to be behind the othe T2 samples in light grey
#also want x labels to be T1 and T2 (to match others)
#order of labels etc= RT1,Rt2,RT2NB,ST1,ST2,ST2NB
ggplot(data=lipid.dat3,aes(x=Time, y=Lipid..Biomass..g.gdw., fill=Toler.Tim.Stat, color=Toler.Tim.Stat)) +
  geom_boxplot() +
  theme_classic(base_size = 20) +
  scale_color_manual(values=c('#95C38B','#1B7837','#1B7837','#8857A8','#762A83','#762A83'), labels=c('Bleached', 'Non-bleached','bb','aa','cc','dd'))+
  scale_fill_manual(values=c('#1b7837','grey90','#7fbf7b','#762a83','#FAFAFA','#af8dc3'), labels=c('RT1', 'RT2','ST1', 'ST2')) +
  ylab('Lipid Biomass (g/g dry weight)') +
  scale_x_discrete(labels=c('Timepoint 1', 'Timepoint 2', 'Timepoint 2NB')) +
  facet_wrap(~Time*Status, scales='free_x') +
  theme(strip.background=element_blank(), strip.text.x=element_blank())

#c('#1b7837','#7fbf7b','#762a83','#af8dc3','#af8dc3','grey')
png("~/Documents/genome_sciences_postdoc/corals/time series/Mcap T1T2 lipids.png", width=7, height=7, units = "in", res = 1200)
T1T2
dev.off()

#2-way ANOVA of Tolerance and Status
res.aov3<-aov(Lipid..Biomass..g.gdw. ~ Tolerance*Status, data=lipid.dat3)
summary(res.aov3)

Df   Sum Sq  Mean Sq F value   Pr(>F)    
Tolerance         1 0.019907 0.019907  15.007 0.000945 ***
  Status            1 0.020698 0.020698  15.604 0.000790 ***
  Tolerance:Status  1 0.006106 0.006106   4.603 0.044374 *  
  Residuals        20 0.026529 0.001326                     
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
12 observations deleted due to missingness






