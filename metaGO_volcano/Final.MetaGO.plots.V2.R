##########################################################################################
##########################################################################################
#MetaGoMic Volcano
library(dplyr)
library(ggplot2)
library(ggrepel)

#read in the new file

#T0
setwd("~/metagomics")

metaGO.T1.dat<-read.csv("T0 R vs S metaGOmics.csv", header=T, na.strings='null')
#resilient = run 1, susceptible = run 2 
#higher PSM in R = -LFC 

#calculate log p values
metaGO.T1.dat$logpval.all <- -log(metaGO.T1.dat$Laplace.corr..Bonf..corr.p.value,10)
 
  

#create limits for colors for T1
metaGO.T1.dat["group"] <- "NotSignificant"
metaGO.T1.dat[which(metaGO.T1.dat['Laplace.corr..Log.2..fold.change'] >= 0.5 & between(metaGO.T1.dat$logpval.all, 0, 2)),"group"] <- "S_none"
 
metaGO.T1.dat[which(metaGO.T1.dat['Laplace.corr..Log.2..fold.change'] >= 0.5 & between (metaGO.T1.dat$logpval.all, 2, 100)),"group"] <- "S_0.01"

metaGO.T1.dat[which(metaGO.T1.dat['Laplace.corr..Log.2..fold.change'] <= -0.5 & between (metaGO.T1.dat$logpval.all, 0, 2)),"group"] <- "R_none"
 
metaGO.T1.dat[which(metaGO.T1.dat['Laplace.corr..Log.2..fold.change'] <= -0.5 & between (metaGO.T1.dat$logpval.all, 2, 100)),"group"] <- "R_0.01"

#write.csv(metaGO.dat, "color.check.csv")
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
names(cols) = c("NotSignificant","S_none", "S_0.1","S_0.05","S_0.01",
                "R_none","R_0.1", "R_0.05", "R_0.01")

cols<- c(NotSignificant="#FFFFFF", S_none="#FFFFFF", S_0.1="#E7D4E8", S_0.05="#af8dc3", S_0.01="#762a83", 
         R_none="#FFFFFF", R_0.1="#d9f0d3", R_0.05="#7fbf7b", R_0.01="#1b7837")

outline <- c("gray", "gray", "#8857A8", "#8857A8", "#8857A8", 
             "gray", "#95C38B", "#95C38B", "#95C38B")
names(outline) = c("NotSignificant","S_none","S_0.1", "S_0.05", "S_0.01",
                "R_none","R_0.1", "R_0.05", "R_0.01")
outline<- c(NotSignificant="gray", S_none="gray", S_0.1="#8857A8", S_0.05="#8857A8", S_0.01="#8857A8", 
         R_none="gray", R_0.1="#95C38B", R_0.05="#95C38B", R_0.01="#95C38B")


##################

#Set label names T1
T1label_limit_S <- metaGO.T1.dat[ which(metaGO.T1.dat$Laplace.corr..Log.2..fold.change >= 0.5 & metaGO.T1.dat$logpval.all >=2), ]
T1label_limit_R <- metaGO.T1.dat[ which(metaGO.T1.dat$Laplace.corr..Log.2..fold.change <= -0.5 & metaGO.T1.dat$logpval.all >= 2), ]
#label_limit_P <- metaGO.dat[ which(metaGO.dat$logpval.all >= 5 ), ]
#label_limit_P10 <- metaGO.dat[ which(metaGO.dat$logpval.all >= 10 ), ]
#can nudge_y label and keep_y p-value to add layers
 
#MetaGOmics - T1
PlotT1 <-ggplot(metaGO.T1.dat) + 
  geom_point(aes(x=Laplace.corr..Log.2..fold.change, y=logpval.all, 
                 color = factor(group),fill = factor(group)), size=4,
             shape=21, stroke = 1.2)+
  scale_colour_manual(values = outline) +
  scale_fill_manual(values = cols)+
  scale_shape_identity()+
  scale_x_continuous(breaks=seq(-5,5,1), limits = c(-5,5))+
  scale_y_continuous(breaks=seq(0,110,25), limits = c(-1.0, 110))+
  xlab("Laplace Corrected Log Fold Change") +
  ylab("-Log(p-value)")+
  ggtitle("Timepoint 1")+
  theme_classic(base_size = 10) +
  theme(legend.position = "none", text = element_text(size=15)) + 
  geom_hline(yintercept = 2, linetype= "dashed", colour="gray")+
  geom_hline(yintercept = 0, linetype= "solid", colour="black")+
  geom_vline (xintercept = 0.5, linetype='dashed', colour="gray")+
  geom_vline (xintercept = -0.5, linetype='dashed', colour="gray") +
  geom_text(x=-1.3, y=110, label="Resilient", size = 7.5,colour ='#1B7837') +
  geom_text(x=4.13, y=110, label="Susceptible", size = 7.5, colour ='#762A83') 


 PlotT1+ geom_text_repel(data = T1label_limit_R,
                       min.segment.length = 8, 
                       segment.color = "gray",
                       nudge_x = -10,
                      
                       hjust =0,
                       seed=42,
                       box.padding = 0.25,
                       aes(x = Laplace.corr..Log.2..fold.change,
                           y = logpval.all),
                       label=T1label_limit_R$GO.name, size=5 )+
   geom_text_repel(data = T1label_limit_S,
                   min.segment.length = 11, 
                   segment.color = "gray", 
                   nudge_x = 10,
                   vjust =-22,
                   hjust =1,
                   seed=42,
                   box.padding = 0.25,
                   aes(x = Laplace.corr..Log.2..fold.change , y = logpval.all),
                   label=T1label_limit_S$GO.name, size=5 )

#########################Timepoint 2##############################
  #T2
  setwd("~/metagomics")
  metaGO.T2.dat<-read.csv("T2 R vs S metaGOmics.csv", header=T, na.strings='null')
  #resilient = run 1, susceptible = run 2 
  #higher PSM in R = -LFC 
  
  #calculate log p values
  metaGO.T2.dat$logpval.all <- -log(metaGO.T2.dat$Laplace.corr..Bonf..corr.p.value,10)
  #write.csv(metaGO.T2.dat,file = "~/Dropbox/Brook_Nunn_WORK/Collaborators/Jackie_Padilla_Gamino/Tanya Brown/metagomics/T2metaGO.dat.csv")
  
  
  #T2Select_genes <- read.csv("T2Select.genes.volcano.csv")
  #T2metaGO.dat <- merge(metaGO.dat, T2Select_genes, by = "PROTEIN.ID", all.x = TRUE)
  
  #View(metaGO.T2.dat)
  #write.csv(metaGO.dat, "check.csv")
  
  #create limits for colors for T2
  metaGO.T2.dat["group"] <- "NotSignificant"
  metaGO.T2.dat[which(metaGO.T2.dat['Laplace.corr..Log.2..fold.change'] >= 0.5 & between(metaGO.T2.dat$logpval.all, 0, 2)),"group"] <- "S_none"
  metaGO.T2.dat[which(metaGO.T2.dat['Laplace.corr..Log.2..fold.change'] >= 0.5 & between (metaGO.T2.dat$logpval.all, 2, 100)),"group"] <- "S_0.01"
  
  metaGO.T2.dat[which(metaGO.T2.dat['Laplace.corr..Log.2..fold.change'] <= -0.5 & between (metaGO.T2.dat$logpval.all, 0, 2)),"group"] <- "R_none"
  metaGO.T2.dat[which(metaGO.T2.dat['Laplace.corr..Log.2..fold.change'] <= -0.5 & between (metaGO.T2.dat$logpval.all, 2, 100)),"group"] <- "R_0.01"
  
  #write.csv(metaGO.dat, "color.check.csv")
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
  cols = c("#FFFFFF", "#FFFFFF", "#FAFAFA", 
           "#FFFFFF", "grey90")
  names(cols) = c("NotSignificant","S_none", "S_0.01",
                  "R_none","R_0.01")
  
  cols<- c(NotSignificant="#FFFFFF", S_none="#FFFFFF", S_0.01="#FAFAFA", 
           R_none="#FFFFFF", R_0.01="grey90")
  
  outline <- c("gray", "gray", "#762A83", 
               "gray", "#1b7837")
  names(outline) = c("NotSignificant","S_none","S_0.01",
                     "R_none","R_0.01")
  outline<- c(NotSignificant="gray", S_none="gray", S_0.01="#762A83", 
              R_none="gray", R_0.01="#1b7837")
  
  
  ##################
  
  #Set label names T2- since so many names, 2nd file wiht NR names 
  metaGO.T2.dat.names<-read.csv("T2metaGO.dat.names.csv", header=T, na.strings='null')
      #View(metaGO.T2.dat.names)
  T2label_limit_S <- metaGO.T2.dat.names[ which(metaGO.T2.dat.names$Laplace.corr..Log.2..fold.change >= 0.5 & metaGO.T2.dat.names$logpval.all >=2), ]
  T2label_limit_R <- metaGO.T2.dat.names[ which(metaGO.T2.dat.names$Laplace.corr..Log.2..fold.change <= -0.5 & metaGO.T2.dat.names$logpval.all >= 2), ]
  #label_limit_P <- metaGO.dat[ which(metaGO.dat$logpval.all >= 5 ), ]
  #label_limit_P10 <- metaGO.dat[ which(metaGO.dat$logpval.all >= 10 ), ]
  #can nudge_y label and keep_y p-value to add layers
  #View(T2label_limit_R)
  #MetaGOmics - T2
  Plot <-ggplot(metaGO.T2.dat) + 
    geom_point(aes(x=Laplace.corr..Log.2..fold.change, y=logpval.all, color = factor(group),
                   fill = factor(group)),size=4, 
               shape=21, stroke=1.2)+
    scale_colour_manual(values = outline) +
    scale_fill_manual(values = cols)+
    scale_shape_identity()+
    scale_x_continuous(breaks=seq(-4,4,1), limits = c(-4,4))+
    scale_y_continuous(breaks=seq(0,110,25), limits = c(-1.0, 110))+
    xlab("Laplace Corrected Log Fold Change") +
    ylab("-Log(p-value)")+
    ggtitle("Timepoint 2")+
    theme_classic(base_size = 10) +
    theme(legend.position = "none", text = element_text(size=15)) + 
    geom_hline(yintercept = 2, linetype= "dashed", colour="gray")+
    geom_hline(yintercept = 0, linetype= "solid", colour="black")+
    geom_vline (xintercept = 0.5, linetype='dashed', colour="gray")+
    geom_vline (xintercept = -0.5, linetype='dashed', colour="gray") +
    geom_text(x=-3.3, y=110, label="Resilient", size = 7.5,colour ='#1B7837') +
    geom_text(x=3.13, y=110, label="Susceptible", size = 7.5, colour ='#762A83') 
  
  Plot+ geom_text_repel(data = T2label_limit_R,
                        min.segment.length = 11, 
                        segment.color = "gray",
                        nudge_x = -10,
                        vjust =-22,
                        hjust =0,
                        seed=42,
                        box.padding = 0.25,
                        aes(x = Laplace.corr..Log.2..fold.change,
                            y = logpval.all),
                        label=T2label_limit_R$GO.name, size=5 )+
    geom_text_repel(data = T2label_limit_S,
                    min.segment.length = 11, 
                    segment.color = "gray",
                    nudge_x = 10,
                    vjust =-22,
                    hjust =1,
                    seed=42,
                    box.padding = 0.25,
                    aes(x = Laplace.corr..Log.2..fold.change , y = logpval.all),
                    label=T2label_limit_S$GO.name, size=5 )
  
  