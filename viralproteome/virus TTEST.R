#load biostats package

library(dplyr)
library(vegan)
setwd("D:/Dropbox/Brook_Nunn_WORK/Collaborators/Jackie_Padilla_Gamino/Tanya Brown/ViralProteins_")
#readin tab seperated values
dat<-read.delim('ABACUS_T1T2_5Virus_output.tsv', header=T, row.names=1)
View(dat)


specstot<-select(dat, contains('NUMSPECSTOT'))

#keep only proteins that are viral
specstot$protein<-row.names(specstot)
virus.prot<-subset(specstot, grepl(c('NP_|YP_|A0|K7'), specstot$protein))
virus.prot = virus.prot[-2,] # removes the second row where there is trypsin
#remove 1st column & 2nd row (trypsin)
virus.prot<-virus.prot[-c(1)]
#also remove columns 5, 13 corresponding to 20BT2 and 28BT2; remove protein name column 25
virus.prot<-virus.prot[-c(5,13,25)]
 

#Are there different numbers of viral spectral counts in resilient vs. susceptible corals at T1 and T2?
#two-tailed t-test to test if there is a difference
#do not do paired because different corals are in R and S groups
RT1<-subset(virus.prot, select=c(X2019_APRIL_16_CORALS_TANYABROWN_10NBT0_22_RT1_NUMSPECSTOT,X2019_APRIL_16_CORALS_TANYABROWN_13NBT0_112_RT1_NUMSPECSTOT,X2019_APRIL_16_CORALS_TANYABROWN_21NBT0_35_RT1_NUMSPECSTOT,X2019_APRIL_16_CORALS_TANYABROWN_26NBT0_81_RT1_NUMSPECSTOT,X2019_APRIL_16_CORALS_TANYABROWN_66NBT0_129_RT1_NUMSPECSTOT,X2019_APRIL_16_CORALS_TANYABROWN_74NBT0_64_RT1_NUMSPECSTOT))
ST1<-subset(virus.prot, select=c(X2019_APRIL_16_CORALS_TANYABROWN_20NBT0_73_ST1_NUMSPECSTOT,X2019_APRIL_16_CORALS_TANYABROWN_25NBT0_102_ST1_NUMSPECSTOT,X2019_APRIL_16_CORALS_TANYABROWN_28NBT0_120_ST1_NUMSPECSTOT,X2019_APRIL_16_CORALS_TANYABROWN_55NBT0_19_ST1_NUMSPECSTOT,X2019_APRIL_16_CORALS_TANYABROWN_5NBT0_39_ST1_NUMSPECSTOT,X2019_APRIL_16_CORALS_TANYABROWN_73NBT0_91_ST1_NUMSPECSTOT))
RT2<-subset(virus.prot, select=c(X2019_APRIL_16_CORALS_TANYABROWN_10BT2_43_RT2_NUMSPECSTOT,X2019_APRIL_16_CORALS_TANYABROWN_13BT2_71_RT2_NUMSPECSTOT,X2019_APRIL_16_CORALS_TANYABROWN_21BT2_40_RT2_NUMSPECSTOT,X2019_APRIL_16_CORALS_TANYABROWN_26BT2_127_RT2_NUMSPECSTOT,X2019_APRIL_16_CORALS_TANYABROWN_66BT2_17_RT2_NUMSPECSTOT,X2019_APRIL_16_CORALS_TANYABROWN_74BT2_88_RT2_NUMSPECSTOT))
ST2<-subset(virus.prot, select=c(X2019_APRIL_16_CORALS_TANYABROWN_25NBT0_32_ST2_NUMSPECSTOT,X2019_APRIL_16_CORALS_TANYABROWN_55BT2_96_ST2_NUMSPECSTOT,X2019_APRIL_16_CORALS_TANYABROWN_5BT2_110_ST2_NUMSPECSTOT,X2019_APRIL_16_CORALS_TANYABROWN_73BT2_21_ST2_NUMSPECSTOT))
ST1.4<-subset(virus.prot, select=c(X2019_APRIL_16_CORALS_TANYABROWN_25NBT0_102_ST1_NUMSPECSTOT,X2019_APRIL_16_CORALS_TANYABROWN_55NBT0_19_ST1_NUMSPECSTOT,X2019_APRIL_16_CORALS_TANYABROWN_5NBT0_39_ST1_NUMSPECSTOT,X2019_APRIL_16_CORALS_TANYABROWN_73NBT0_91_ST1_NUMSPECSTOT))

sumRT1<-colSums(RT1)
sumST1<-colSums(ST1)
sumRT2<-colSums(RT2)
sumST2<-colSums(ST2)
sumST1.4<-colSums(ST1.4)




#compare VALUES ACROSS ROWS R and S T1
RvST1<-t.test((RT1), (ST1), alternative='two.sided', paired=F, conf.level=0.95)
RvST1
Welch Two Sample t-test
data:  (RT1) and (ST1)
t = 0.21904, df = 1773.6, p-value = 0.8266
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
  -0.1522775  0.1905658
sample estimates:
  mean of x mean of y 
0.5315315 0.5123874 


#compare VALUES ACROSS ROWS R and S T2
RvST2<-t.test((RT2), (ST2), alternative='two.sided', paired=F, conf.level=0.95)
RvST2
Welch Two Sample t-test

data:  (RT2) and (ST2)
t = -0.65347, df = 1021.9, p-value = 0.5136
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
  -0.3065249  0.1533717
sample estimates:
  mean of x mean of y 
0.5281532 0.6047297 



#compare R and S T1
resT1<-t.test((sumRT1), (sumST1), alternative='two.sided', paired=F, conf.level=0.95)
resT1
Welch Two Sample t-test

data:  (sumRT1) and (sumST1)
t = 1.0808, df = 7.7056, p-value = 0.3124
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
  -3.252281  8.918948
sample estimates:
  mean of x mean of y 
78.66667  75.83333 

#compare R and S T2
resT2<-t.test((sumRT2), (sumST2), alternative='two.sided', paired=F, conf.level=0.95)
resT2

Welch Two Sample t-test

data:  (sumRT2) and (sumST2)
t = -1.2294, df = 4.4351, p-value = 0.2801
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
  -35.96702  13.30036
sample estimates:
  mean of x mean of y 
78.16667  89.50000 

#compare T1 and T2 R
resR<-t.test((sumRT1), (sumRT2), alternative='two.sided', paired=T, conf.level=0.95)
resR


Paired t-test

data:  (sumRT1) and (sumRT2)
t = 0.12039, df = 5, p-value = 0.9089
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
  -10.17643  11.17643
sample estimates:
  mean of the differences 
0.5 

#compare T1 and T2 S
resS<-t.test((sumST1.4), (sumST2), alternative='two.sided', paired=T, conf.level=0.95)
resS

Paired t-test

data:  (sumST1.4) and (sumST2)
t = -1.5476, df = 3, p-value = 0.2195
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
  -42.78898  14.78898
sample estimates:
  mean of the differences 
-14 

#compare T1 and T2 S unpaired out of curriousity

resSunp<-t.test((sumST1.4), (sumST2), alternative='two.sided', paired=F, conf.level=0.95)
resSunp
Welch Two Sample t-test

data:  (sumST1.4) and (sumST2)
t = -1.6442, df = 3.3263, p-value = 0.1897
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
  -39.6541  11.6541
sample estimates:
  mean of x mean of y 
75.5      89.5 