install.packages("viridis")

library(pheatmap)
library(RColorBrewer)
library(viridis)
library(dplyr)

setwd("~/Heatmap")
abacusoutput <- read.csv("ABACUS_03092020_livevdead_output.csv")
View(abacusoutput)


#PREPARE NSAF file for NMDS:  next subset will keep any column with header PROT, DEFLINE or ADJNSAF
subset_abacusADJNSAF<-abacusoutput[,grepl('PROT|DEFLINE|ADJNSAF', names(abacusoutput))]
View(subset_abacusADJNSAF)


#across all columns replace "X2019_APRIL_16_CORALS_" with ""
names(subset_abacusADJNSAF) <- gsub("X2019_APRIL_16_CORALS_", "", names(subset_abacusADJNSAF))




#Samples 10, 13, 21, 26, 66, and 74 lived. Samples 5, 20, 25, 28, 55, and 73 eventually died
#sample 25NBT0_32 is an 25NBT2_32. The 102 is NBT0. The original NB T2 was not a good run (it was the bleached version)

names(subset_abacusADJNSAF) <- gsub("25NBT0_32", "25NBT2_32", names(subset_abacusADJNSAF))
View(subset_abacusADJNSAF)


#Subset all T0 and T2 samples
subset_abacusADJNSAFT0T2<-subset_abacusADJNSAF[,grepl('PROTID|T0|T2', names(subset_abacusADJNSAF),)]
View(subset_abacusADJNSAFT0T2)

#make csv- to create average NSAF values
write.csv(subset_abacusADJNSAFT0T2, file = "SUBSET.ABACUS.ADJNSAF.T0.T2.csv")

####SHORT LIST of most SIGNIFICANT- low p values####
#Make new list with the IDs that are logfold <1.5 and Zscore >3 in T0 or T2 in excel- save inheatmap folder
logfold_1.5_z_3 <- read.csv("logfold_1.5_z_3.csv")
View(logfold_1.5_z_3)


#Make new list with the IDs that are logfold 1 or <1 in T0 or T2 in excel- save inheatmap folder
logfold_1_T0T2 <- read.csv("logfold_1_T0T2.csv")
View(logfold_1_T0T2)

#join choice "left" chooses 1st listed file and adds to it
Sig1_NSAF <- left_join(logfold_1_T0T2,subset_abacusADJNSAFT0T2, by = c("PROTID"="PROTID"))
View(Sig1_NSAF)

#Average NSAF for following samples
#Samples 10, 13, 21, 26, 66, and 74 lived. Samples 5, 20, 25, 28, 55, and 73 died

#View(T2Sus)
#T0Sus Row means
Sig1_NSAF$AvgT0Sus = rowMeans(Sig1_NSAF[,grepl('5NBT0|20NBT0|25NBT0|28NBT0|55NBT0|73NBT0', names(Sig1_NSAF),)],na.rm=FALSE)
#T0Res Row means
Sig1_NSAF$AvgT0Res = rowMeans(Sig1_NSAF[,grepl('10NBT0|13NBT0|21NBT0|26NBT0|66NBT0|74NBT0', names(Sig1_NSAF),)],na.rm=FALSE)
#T2Sus Row means
Sig1_NSAF$AvgT2Sus = rowMeans(Sig1_NSAF[,grepl('5BT2|20BT2|25BT2|28BT2|55BT2|73BT2', names(Sig1_NSAF),)],na.rm=FALSE)
#T2Res Row means
Sig1_NSAF$AvgT2Res = rowMeans(Sig1_NSAF[,grepl('10BT2|13BT2|21BT2|26BT2|66BT2|74BT2', names(Sig1_NSAF),)],na.rm=FALSE)



#Subset the NSAF averages
Sig1_AVG_NSAF<-Sig1_NSAF[,grepl('PROTID|AvgT0Sus|AvgT0Res|AvgT2Sus|AvgT2Res', names(Sig1_NSAF),)]
View(Sig1_AVG_NSAF)
#Normalize each value by the row mean
Sig1_AVG_NSAF$normAvgT0Sus = (Sig1_AVG_NSAF$AvgT0Sus)/rowMeans(Sig1_AVG_NSAF[,grepl('AvgT0Sus|AvgT0Res|AvgT2Sus|AvgT2Res', 
                             names(Sig1_AVG_NSAF),)], 
                             na.rm =FALSE)
Sig1_AVG_NSAF$normAvgT0Res = (Sig1_AVG_NSAF$AvgT0Res)/rowMeans(Sig1_AVG_NSAF[,grepl('AvgT0Sus|AvgT0Res|AvgT2Sus|AvgT2Res', 
                             names(Sig1_AVG_NSAF),)], 
                             na.rm =FALSE)
Sig1_AVG_NSAF$normAvgT2Sus = (Sig1_AVG_NSAF$AvgT2Sus)/rowMeans(Sig1_AVG_NSAF[,grepl('AvgT0Sus|AvgT0Res|AvgT2Sus|AvgT2Res', 
                             names(Sig1_AVG_NSAF),)], 
                             na.rm =FALSE)
Sig1_AVG_NSAF$normAvgT2Res = (Sig1_AVG_NSAF$AvgT2Res)/rowMeans(Sig1_AVG_NSAF[,grepl('AvgT0Sus|AvgT0Res|AvgT2Sus|AvgT2Res', 
                             names(Sig1_AVG_NSAF),)], 
                             na.rm =FALSE)



#pheatmap only takes row name and then datamatrix.-change row name to what I want to view on figure
#rownames(annotated.norm.AVG.NSAFplot) = sapply(annotated.norm.AVG.NSAFplot$m.ID,function(x) strsplit(as.character(x),split = "\\\\")[[1]][1])
#norm.AVG_NSAFplot.1<-annotated.norm.AVG.NSAFplot[,grepl('normAvgT0Sus|normAvgT0Res|normAvgT2Sus|normAvgT2Res', names(annotated.norm.AVG.NSAFplot),)]
Sig1_AVG_NSAF.1 <- Sig1_AVG_NSAF[-1]
row.names(Sig1_AVG_NSAF.1) <- Sig1_AVG_NSAF$PROTID
View(Sig1_AVG_NSAF.1)
norm.AVG_NSAFplot<-Sig1_AVG_NSAF.1[,grepl('normAvgT0Sus|normAvgT0Res|normAvgT2Sus|normAvgT2Res', names(Sig1_AVG_NSAF.1),)]
write.csv(norm.AVG_NSAFplot,'norm.AVG_NSAFplot')
View(norm.AVG_NSAFplot)


#make heatmap and extracting the list of elements after clustering
out<-pheatmap(norm.AVG_NSAFplot, 
              color = colorRampPalette((brewer.pal(n = 5, name ="Blues")))(100),
              border_color = NA,
              cellwidth = 35, 
              cellheight = 6,
              angle_col="0",
              show_rownames=T, 
              show_colnames=T,
              clustering_distance_rows="correlation",
              cluster_rows=T, 
              cluster_cols= T, 
              fontsize=5, 
              treeheight_row = 30,
              cutree_rows =12)
#generate file of cluter list
Clust_List<-norm.AVG_NSAFplot[c(out$tree_row[['order']]), out$tree_col[['order']]]
View(Clust_List)
Clust_List$PROTID<-row.names(Clust_List)

#Annotate final file with all options
MC.MASTER.ANNOTATE <- read.delim("MC.MASTER.ANNOTATE.txt")
Clust_list.annotate <- left_join(Clust_List,MC.MASTER.ANNOTATE, by = c("PROTID"="PROTID"))
View(Clust_list.annotate)
write.csv(Clust_list.annotate, file = "Clust.List.annotate.correlation.csv")




#
