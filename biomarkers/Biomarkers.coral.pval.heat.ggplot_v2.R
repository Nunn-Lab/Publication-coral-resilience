

setwd("~/Biomarkers")
coral.nsaf<-read.csv('NSAF tanya coral T0 T2.csv', header=T, row.names=1)


library(pheatmap)
library(RColorBrewer)
library(viridis)
library(dplyr)

#average NSAF for R and S at each time point, for bleached and NB
coral.R<-subset(coral.nsaf, select=c(X10NBT0_22, X13NBT0_112, X21NBT0_35,X26NBT0_81,X66NBT0_129,X74NBT0_64, X10BT2_43, X13BT2_71,X21BT2_40,X26BT2_127,X66BT2_17,X74BT2_88, X10NBT2_90, X13NBT2_101, X21NBT2_26,X26NBT2_49,X66NBT2_119,X74NBT2_42))
coral.S<-subset(coral.nsaf, select=c(X20NBT0_73,X25NBT0_102,X28NBT0_120,X55NBT0_19,X5NBT0_39,X73NBT0_91 ,X55BT2_96,X5BT2_110,X73BT2_21, X20NBT2_123, X25NBT0_32,X28NBT2_24,X55NBT2_79,X5NBT2_66,X73NBT2_108))

coral.R.T1<-rowMeans(coral.R[,1:6])
coral.R.T2B<-rowMeans(coral.R[,7:12])
coral.R.T2NB<-rowMeans(coral.R[,13:18])

coral.S.T1<-rowMeans(coral.S[,1:6])
coral.S.T2B<-rowMeans(coral.R[,7:9])
coral.S.T2NB<-rowMeans(coral.R[,10:15])

coral.avg<-cbind(coral.R.T1, coral.R.T2B, coral.R.T2NB, coral.S.T1, coral.S.T2B, coral.S.T2NB)


# check to see if data are normally distributed with shapiro wilk test before T-test
  #shapiro.test(coral.R.T1)
  #shapiro.test(coral.R.T2B)
  #shapiro.test(coral.S.T1)
  #shapiro.test(coral.R.T2B)
  #hist(coral.R.T1,xlim = c(10,500),breaks=100)



#1) do t-test to see if T1 R > T1 S, put result in new column
coral.T1<-subset(coral.nsaf, select=c(X10NBT0_22, X13NBT0_112, X21NBT0_35,X26NBT0_81,X66NBT0_129,X74NBT0_64,X20NBT0_73,X25NBT0_102,X28NBT0_120,X55NBT0_19,X5NBT0_39,X73NBT0_91))

coral.T1$protein<-row.names(coral.T1)
col1 = protein
col2= NSAF
col3=R or S
col4=coral ID
coral10<-subset(coral.T1, select=c(protein, X10NBT0_22))
coral13<-subset(coral.T1, select=c(protein, X13NBT0_112))
coral21<-subset(coral.T1, select=c(protein, X21NBT0_35))
coral26<-subset(coral.T1, select=c(protein, X26NBT0_81))
coral66<-subset(coral.T1, select=c(protein, X66NBT0_129))
coral74<-subset(coral.T1, select=c(protein, X74NBT0_64))
coral20<-subset(coral.T1, select=c(protein, X20NBT0_73))
coral25<-subset(coral.T1, select=c(protein, X25NBT0_102))
coral28<-subset(coral.T1, select=c(protein, X28NBT0_120))
coral55<-subset(coral.T1, select=c(protein, X55NBT0_19))
coral5<-subset(coral.T1, select=c(protein, X5NBT0_39))
coral73<-subset(coral.T1, select=c(protein, X73NBT0_91))

#not sure why next line there writtenby emma-Brook doesn't know why here
names(coral73)[names(coral73)=='X73NBT0_91']<-'NSAF'

#Does T test on whole dataset to see if R>S
coral.ttest<-rbind(coral10$X10NBT0_22, coral13$X13NBT0_112,coral13$X13NBT0_112, coral21$X21NBT0_35, coral26$X26NBT0_81, coral66$X66NBT0_129, coral74$X74NBT0_64, coral20$X20NBT0_73, coral25$X25NBT0_102, coral28$X28NBT0_120, coral55$X55NBT0_19, coral5$X5NBT0_39, coral73$NSAF)
RorS<-c(rep('R', 15636), rep('S', 15636))
coralID<-c(rep('10T1',2606), rep('13T1', 2606), rep('21T1', 2606),rep('26T1', 2606), rep('66T1', 2606), rep('74T1', 2606), rep('20T1', 2606), rep('25T1', 2606), rep('28T1', 2606), rep('55T1', 2606), rep('5T1', 2606), rep('73T1', 2606))
coral.ttest$RorS<-RorS
coral.ttest$coralID<-coralID

coral.T1.t<-t.test(NSAF~RorS, data=coral.ttest, paired=T, alternative="two.sided")
#overall, no difference p=0.5709

#t-test on each protein/row- important test
prot.pval<-apply(coral.T1[-13], 1, function(x) t.test(x[1:6], x[7:12])$p.value)


#add rows to pvlaue sheet that have T1 avg and T2 rates
coral.T1.pval<-cbind(coral.T1[-13], prot.pval)
coral.T1.pval$T1S.avg<-coral.S.T1
coral.T1.pval$T1R.avg<-coral.R.T1

coralnsaf.R.T2B<-(coral.R[,7:12])
coralnsaf.S.T2B<-(coral.R[,7:9])
coral.T1T2.pval<-cbind(coral.T1.pval,coralnsaf.S.T2B)
coral.T1T2.pval<-cbind(coral.T1T2.pval,coralnsaf.R.T2B)
coral.T1T2.pval$coralT2B.R.avg<-coral.R.T2B
coral.T1T2.pval$coralT2B.S.avg<-coral.S.T2B


#This is all of data 
write.csv(coral.T1T2.pval, 'proteins with pvalue.nsaf.rates.csv', quote=F)

# using subset function -you can add select only colums e.g. select=c(ID, Weight)
T1.p.val.01 <- subset(coral.T1T2.pval, prot.pval <= .01)
write.csv(T1.p.val.01, 'p.value.01.rates.nsaf.csv')

#keep all proteins with p-value < or = 0.05
T1.p.val.05<-subset(coral.T1T2.pval, prot.pval<=0.05)
write.csv(T1.p.val.05, 'p.value.05.rates.nsaf', quote=F)




#HEATMAPS
# data 0.01 is df= T1.p.val.01
pval.01.heat.data <-subset(T1.p.val.01,
                           select=c(X10NBT0_22, X13NBT0_112, X21NBT0_35,X26NBT0_81,X66NBT0_129,X74NBT0_64,
                                    X20NBT0_73,X25NBT0_102,X28NBT0_120,X55NBT0_19,X5NBT0_39,X73NBT0_91,T1S.avg,T1R.avg))
#normalize resilient by row mean
pval.01.heat.data$rw.mean<- ((pval.01.heat.data$T1R.avg + pval.01.heat.data$T1S.avg)/2)
pval.01.heat.data$normX10NBT0_22 <-(pval.01.heat.data$X10NBT0_22 / pval.01.heat.data$rw.mean)
pval.01.heat.data$normX13NBT0_112 <-(pval.01.heat.data$X13NBT0_112 / pval.01.heat.data$rw.mean)
pval.01.heat.data$normX21NBT0_35 <-(pval.01.heat.data$X21NBT0_35 / pval.01.heat.data$rw.mean)
pval.01.heat.data$normX26NBT0_81 <-(pval.01.heat.data$X26NBT0_81 / pval.01.heat.data$rw.mean)
pval.01.heat.data$normX66NBT0_129 <-(pval.01.heat.data$X66NBT0_129 / pval.01.heat.data$rw.mean)
pval.01.heat.data$normX74NBT0_64 <-(pval.01.heat.data$X74NBT0_64 / pval.01.heat.data$rw.mean)

#normalize susceptible by row
pval.01.heat.data$normX20NBT0_73 <-(pval.01.heat.data$X20NBT0_73 / pval.01.heat.data$rw.mean)
pval.01.heat.data$normX25NBT0_102 <-(pval.01.heat.data$X25NBT0_102 / pval.01.heat.data$rw.mean)
pval.01.heat.data$normX28NBT0_120 <-(pval.01.heat.data$X28NBT0_120 / pval.01.heat.data$rw.mean)
pval.01.heat.data$normX55NBT0_19 <-(pval.01.heat.data$X55NBT0_19 / pval.01.heat.data$rw.mean)
pval.01.heat.data$normX5NBT0_39 <-(pval.01.heat.data$X5NBT0_39 / pval.01.heat.data$rw.mean)
pval.01.heat.data$normX73NBT0_91 <-(pval.01.heat.data$X73NBT0_91 / pval.01.heat.data$rw.mean)

#replace non exist values NaN from 0/0 with "0"
pval.01.heat.data[is.na(pval.01.heat.data)] <- 0

#subsetting NSAF data for mapping  
pval.01.heat.nsafdata <-subset(pval.01.heat.data,
                               select=c(X10NBT0_22, X13NBT0_112, X21NBT0_35, X26NBT0_81, X66NBT0_129, X74NBT0_64,
                                        X20NBT0_73,X25NBT0_102,X28NBT0_120, X55NBT0_19, X5NBT0_39, X73NBT0_91))

##NORMALIZED NSAF per protein
#subsetting normalized data for mapping  
pval.01.heat.normdata <-subset(pval.01.heat.data,
                           select=c(normX10NBT0_22, normX13NBT0_112, normX21NBT0_35, normX26NBT0_81, normX66NBT0_129,normX74NBT0_64,
                                    normX20NBT0_73,normX25NBT0_102,normX28NBT0_120, normX55NBT0_19, normX5NBT0_39, normX73NBT0_91))
#change col names
colnames(pval.01.heat.normdata) <- c("T1R_10", "T1R_13","T1R_21", "T1R_26", "T1R_66","T1R_74",
                                 "T1S_20","T1S_25","T1S_28","T1S_55","T1S_5","T1S_73")
#change row names
pval.01.row.names <- read.csv("~/Biomarkers/pval.01.row.names.csv", header=FALSE)




#make heatmap and extracting the list of elements after clustering
p.val.heat<-pheatmap(pval.01.heat.normdata, 
              color = colorRampPalette((brewer.pal(n = 5, name ="BuPu")))(100),
              border_color = NA,
              cellwidth = 20, 
              cellheight = 20,
              angle_col="45",
              show_rownames=T, 
              show_colnames=T,
              clustering_distance_rows="correlation",
              cluster_rows=T, 
              cluster_cols= T, 
              fontsize=12, 
              treeheight_row = 30,
              cutree_rows = 2,
              cutree_cols = 2)



library(ggplot2)
#Make small ggplot facet that show T1 and T2 avg and error bars.
#subsetting NSAF data for mapping  dataset pval.01.heat.nsafdata has 

#grab NSAF for R and S at each time point, for bleached and NB
linePlotData.1 <- read.csv("linePlotData.1.csv", header=T)

p<- ggplot(linePlotData.1, aes(x=Timepoint, y=R.avg, group=Gene, color=R.S)) + 
  geom_line() + 
  geom_point(size=3) + 
  scale_fill_brewer(palette="Set1") +
  geom_errorbar(aes(ymin=R.avg-R.std, ymax=R.avg+R.std), width=.01, position=position_dodge(1))

#final labels and colors
p+labs(title="Time Dependent NSAF", x="Timepoint", y = "NSAF")+ 
  facet_grid(Gene ~ R.S, scales="free")+
  theme_classic()
  #facet_wrap(~ protein+ R.S,scales="free_y",ncol=2) #+theme_classic()
          
