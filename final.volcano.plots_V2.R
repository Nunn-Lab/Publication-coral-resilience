##########################################################################################
##########################################################################################

library(dplyr)
library(ggplot2)
library(ggrepel)

#read in the new file
setwd("~/volcano")
QSpecProt.dat<-read.csv('QSPEC.T1.R.v.S.pvalue.log.axes.csv', header=T, na.strings='null')

#QSpecProt.dat_MT0<-read.csv('T0 R vs S QSpecProtmics.csv', header=T, na.strings='null')



T1Select_genes <- read.csv("~/volcano/T1Select.genes.volcano.csv")

T1QSpecProt.dat <- merge(QSpecProt.dat, T1Select_genes, by = "PROTEIN.ID", all.x = TRUE)


#write.csv(QSpecProt.dat, "check.csv")

#create limits for colors for T1
QSpecProt.dat["group"] <- "NotSignificant"
QSpecProt.dat[which(QSpecProt.dat['axes.LogFoldChange'] >= 0.5 & between(QSpecProt.dat$axes.log.p.value, 0, 1.0)),"group"] <- "S_none"
QSpecProt.dat[which(QSpecProt.dat['axes.LogFoldChange'] >= 0.5 & between (QSpecProt.dat$axes.log.p.value, 1.0, 1.301)),"group"] <- "S_0.1"
QSpecProt.dat[which(QSpecProt.dat['axes.LogFoldChange'] >= 0.5 & between (QSpecProt.dat$axes.log.p.value, 1.301, 2)),"group"] <- "S_0.05"
QSpecProt.dat[which(QSpecProt.dat['axes.LogFoldChange'] >= 0.5 & between (QSpecProt.dat$axes.log.p.value, 2, 10)),"group"] <- "S_0.01"

QSpecProt.dat[which(QSpecProt.dat['axes.LogFoldChange'] <= -0.5 & between (QSpecProt.dat$axes.log.p.value, 0, 1.0)),"group"] <- "R_none"
QSpecProt.dat[which(QSpecProt.dat['axes.LogFoldChange'] <= -0.5 & between (QSpecProt.dat$axes.log.p.value, 1.0, 1.301)),"group"] <- "R_0.1"
QSpecProt.dat[which(QSpecProt.dat['axes.LogFoldChange'] <= -0.5 & between (QSpecProt.dat$axes.log.p.value, 1.301, 2)),"group"] <- "R_0.05"
QSpecProt.dat[which(QSpecProt.dat['axes.LogFoldChange'] <= -0.5 & between (QSpecProt.dat$axes.log.p.value, 2, 10)),"group"] <- "R_0.01"

#write.csv(QSpecProt.dat, "color.check.csv")
########Purples####
#Susceptible p<0.01 = hex 762a83  ; border =hex 8857A8  **note that -log10 (.o1) would be >2
#Susceptible p<0.05 = hex af8dc3  ; border =hex 8857A8     **note that -log10 (.05) would be >1.301
#Susceptible p<0.10 = hex E7D4E8  ; border =hex 8857A8     **note that -log10 (0.10) would be > 1.0
#all other Susceptible no fill  or hex FFFFFF ;   ; border =hex 8857A8

########Greens####
#Resilient  p<0.01 = hex 1b7837 ; border =hex 95C38B   **note that -log10 (.o1) would be >2
#Resilient  p<0.05 = hex 7fbf7b ; border =hex 95C38B     **note that -log10 (.05) would be >1.301
#Resilient  p<0.10 = hex d9f0d3 ; border =hex 95C38B     **note that -log10 (0.10) would be > 1.0
#all other Resilient  no fill or hex FFFFFF;  ; border =hex 95C38B

#dark purple #762a83 outline #8857A8
#med purple #af8dc3 outline #8857A8
#light purple #E7D4E8 outline #C69AC9

#dark green #1b7837 outline #95C38B
#med green #7fbf7b outline #2E9337
#light green #d9f0d3 outline #004E00

#create colors
cols = c("#FFFFFF", "#FFFFFF", "#E7D4E8", "#af8dc3", "#762a83", 
         "#FFFFFF", "#d9f0d3", "#7fbf7b", "#1b7837")
names(cols) = c("NotSignificant","S_none","S_0.1", "S_0.05", "S_0.01",
                "R_none","R_0.1", "R_0.05", "R_0.01")

cols<- c(NotSignificant="#FFFFFF", S_none="#FFFFFF", S_0.1="#E7D4E8", S_0.05="#af8dc3", S_0.01="#762a83", 
         R_none="#FFFFFF", R_0.1="#d9f0d3", R_0.05="#7fbf7b", R_0.01="#1b7837")

outline <- c("gray", "#8857A8", "#8857A8", "#8857A8", "#8857A8", 
             "#95C38B", "#95C38B", "#95C38B", "#95C38B")
names(outline) = c("NotSignificant","S_none","S_0.1", "S_0.05", "S_0.01",
                "R_none","R_0.1", "R_0.05", "R_0.01")
outline<- c(NotSignificant="gray", S_none="#8857A8", S_0.1="#8857A8", S_0.05="#8857A8", S_0.01="#8857A8", 
         R_none="#95C38B", R_0.1="#95C38B", R_0.05="#95C38B", R_0.01="#95C38B")


##################

#Set label names T1
T1label_limit_S <- T1QSpecProt.dat[ which(T1QSpecProt.dat$axes.LogFoldChange >= 0.5 & T1QSpecProt.dat$axes.log.p.value > 0.5), ]
T1label_limit_R <- T1QSpecProt.dat[ which(T1QSpecProt.dat$axes.LogFoldChange <= -0.5 & T1QSpecProt.dat$axes.log.p.value > -0.5), ]
 
#Qspec - T1
ggplot(QSpecProt.dat) + 
  geom_point(aes(x=axes.LogFoldChange, y=axes.log.p.value, size=2,
                 color = factor(group),fill = factor(group)), 
             shape=21)+
  scale_colour_manual(values = outline) +
  scale_fill_manual(values = cols)+
  scale_shape_identity()+
  scale_x_continuous(breaks=seq(-2.5,2.5,1), limits = c(-2.5,2.5))+
  scale_y_continuous(breaks=seq(0,10,2.5), limits = c(0, 10.47))+
  xlab("Log Fold Change") +
  ylab("Log(p-value)")+
  ggtitle("Timepoint 1")+
  theme_classic(base_size = 10) +
  theme(legend.position = "none", text = element_text(size=15)) + 
  geom_hline(yintercept = 0, linetype= "dashed", colour="gray")+
  geom_vline (xintercept = 0.5, linetype='dashed', colour="gray")+
  geom_vline (xintercept = -0.5, linetype='dashed', colour="gray")+  
  geom_text(x=-2.3, y=10.65, label="Resilient", size = 7.5,colour ='#1B7837') +
  geom_text(x=2.13, y=10.65, label="Susceptible", size = 7.5, colour ='#762A83')+  
 
   geom_text_repel(data = T1label_limit_R,
                  min.segment.length = 2, 
                  seed=42,
                  box.padding = 0.25,
                  aes(x = axes.LogFoldChange , y = axes.log.p.value),
                  label=T1label_limit_R$GENE.ID, size=5) +
  geom_text_repel(data = T1label_limit_S,
                  min.segment.length = 2, 
                  seed=42,
                  box.padding = 0.25,
                  aes(x = axes.LogFoldChange , y = axes.log.p.value ),
                  label=T1label_limit_S$GENE.ID, size=5) 
#sometimes in code above we get error Aes must be length one or same as data--- look at the column name for gene.ID- it appears differently deending on mac vs pc
######Timepoint 2#######


#read in the new file
T2QSpecProt.dat<-read.csv('QSPEC.T2.R.v.S.pvalue.log.axes.csv', header=T, na.strings='null')
 

T2Select_genes <- read.csv("T2Select.genes.volcano.csv")
T2QSpecProt.dat <- merge(T2QSpecProt.dat, T2Select_genes, by = "PROTEIN.ID", all.x = TRUE)




#create limits for colors for T2
T2QSpecProt.dat["group"] <- "NotSignificant"
T2QSpecProt.dat[which(T2QSpecProt.dat['axes.LogFoldChange'] >= 0.5 & between(T2QSpecProt.dat$axes.log.p.value, 0, 1.0)),"group"] <- "S_none"
T2QSpecProt.dat[which(T2QSpecProt.dat['axes.LogFoldChange'] >= 0.5 & between (T2QSpecProt.dat$axes.log.p.value, 1.0, 1.301)),"group"] <- "S_0.1"
T2QSpecProt.dat[which(T2QSpecProt.dat['axes.LogFoldChange'] >= 0.5 & between (T2QSpecProt.dat$axes.log.p.value, 1.301, 2)),"group"] <- "S_0.05"
T2QSpecProt.dat[which(T2QSpecProt.dat['axes.LogFoldChange'] >= 0.5 & between (T2QSpecProt.dat$axes.log.p.value, 2, 10)),"group"] <- "S_0.01"

T2QSpecProt.dat[which(T2QSpecProt.dat['axes.LogFoldChange'] <= -0.5 & between (T2QSpecProt.dat$axes.log.p.value, 0, 1.0)),"group"] <- "R_none"
T2QSpecProt.dat[which(T2QSpecProt.dat['axes.LogFoldChange'] <= -0.5 & between (T2QSpecProt.dat$axes.log.p.value, 1.0, 1.301)),"group"] <- "R_0.1"
T2QSpecProt.dat[which(T2QSpecProt.dat['axes.LogFoldChange'] <= -0.5 & between (T2QSpecProt.dat$axes.log.p.value, 1.301, 2)),"group"] <- "R_0.05"
T2QSpecProt.dat[which(T2QSpecProt.dat['axes.LogFoldChange'] <= -0.5 & between (T2QSpecProt.dat$axes.log.p.value, 2, 10)),"group"] <- "R_0.01"

#write.csv(T2QSpecProt.dat, "color.check.csv")
########Purples####
#Susceptible p<0.01 = hex 762a83  ; border =hex 8857A8  **note that -log10 (.o1) would be >2
#Susceptible p<0.05 = hex af8dc3  ; border =hex 8857A8     **note that -log10 (.05) would be >1.301
#Susceptible p<0.10 = hex E7D4E8  ; border =hex 8857A8     **note that -log10 (0.10) would be > 1.0
#all other Susceptible no fill  or hex FFFFFF ;   ; border =hex 8857A8

########Greens####
#Resilient  p<0.01 = hex 1b7837 ; border =hex 95C38B   **note that -log10 (.o1) would be >2
#Resilient  p<0.05 = hex 7fbf7b ; border =hex 95C38B     **note that -log10 (.05) would be >1.301
#Resilient  p<0.10 = hex d9f0d3 ; border =hex 95C38B     **note that -log10 (0.10) would be > 1.0
#all other Resilient  no fill or hex FFFFFF;  ; border =hex 95C38B

#dark purple #762a83 outline #8857A8
#med purple #af8dc3 outline #8857A8
#light purple #E7D4E8 outline #C69AC9

#dark green #1b7837 outline #95C38B
#med green #7fbf7b outline #2E9337
#light green #d9f0d3 outline #004E00

#create colors
#fill colors in order
cols = c("grey90","grey90","grey90","grey90","grey90",
         "grey90","grey90","grey90","grey90")

#assign fill colors to names
names(cols) = c("NotSignificant","S_none","S_0.1", "S_0.05", "S_0.01",
                "R_none","R_0.1", "R_0.05", "R_0.01")

cols<- c(NotSignificant="white", S_none="white", S_0.1="grey100", S_0.05="grey88", S_0.01="grey75", 
         R_none="white", R_0.1="grey100", R_0.05="grey88", R_0.01="grey75")
#outline colors in order
outline <- c("gray", "#8857A8", "#8857A8", "#8857A8", "#8857A8", 
             "#95C38B", "#D9F0D3", "#7fbf7b", "#1B7837")

#assign fill colors to names
names(outline) = c("NotSignificant","S_none","S_0.1", "S_0.05", "S_0.01",
                   "R_none","R_0.1", "R_0.05", "R_0.01")
outline<- c(NotSignificant="gray", S_none="#E7D4E8", S_0.1="#E7D4E8", S_0.05="#AF8DC3", S_0.01="#762A83", 
            R_none="#D9F0D3", R_0.1="#95C38B", R_0.05="#7FBF7B", R_0.01="#1B7837")



#Set label names T2
T2label_limit_S <- T2QSpecProt.dat[ which(T2QSpecProt.dat$axes.LogFoldChange >= 0.5 & T2QSpecProt.dat$axes.log.p.value > .05), ]
T2label_limit_R <- T2QSpecProt.dat[ which(T2QSpecProt.dat$axes.LogFoldChange <= -0.5 & T2QSpecProt.dat$axes.log.p.value > -.05), ]


#View(T2label_limit_R)
# to thicken outline line add stroke = 3 next to shape
#Qspec - T2
ggplot(T2QSpecProt.dat) + 
  geom_point(aes(x=axes.LogFoldChange, y=axes.log.p.value, size=2,
                 color = factor(group),fill = factor(group)), 
             shape=21)+
  scale_colour_manual(values = outline) +
  scale_fill_manual(values = cols)+
  scale_shape_identity()+
  scale_x_continuous(breaks=seq(-2.5,2.5,1), limits = c(-2.5,2.5))+
  scale_y_continuous(breaks=seq(0,10,2.5), limits = c(0, 10.47))+
  xlab("Log Fold Change") +
  ylab("Log(p-value)")+
  ggtitle("Timepoint 2")+
  theme_classic(base_size = 10) +
  theme(legend.position = "none", text = element_text(size=15)) + 
  geom_hline(yintercept = 0, linetype= "dashed", colour="gray")+
  geom_vline (xintercept = 0.5, linetype='dashed', colour="gray")+
  geom_vline (xintercept = -0.5, linetype='dashed', colour="gray") +
  geom_text(x=-2.3, y=10.65, label="Resilient", size = 7.5,colour ='#1B7837') +
  geom_text(x=2.13, y=10.65, label="Susceptible", size = 7.5, colour ='#762A83')+ 
  geom_text_repel(data = T2label_limit_R,
                  min.segment.length = 2, 
                  seed=42,
                  box.padding = 0.25,
                  aes(x = axes.LogFoldChange , y = axes.log.p.value),
                  label=T2label_limit_R$GENE.ID, size=5)+
  geom_text_repel(data = T2label_limit_S,
                  min.segment.length = 2, 
                  seed=42,
                  box.padding = 0.25,
                  aes(x = axes.LogFoldChange , y = axes.log.p.value ),
                  label=T2label_limit_S$GENE.ID, size=5) 
 
                