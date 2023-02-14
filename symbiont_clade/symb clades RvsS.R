library(ggplot2)

symb.dat<-read.csv('symbiont data RvsS.csv')
symb.dat$Time.point<-as.factor(symb.dat$Time.point)

symb3<-ggplot(data=symb.dat, aes(x=factor(interaction(Time.point, Bleaching), levels=c('1.T1', '2.NB', '2.B')), y=Abundance, fill=Symbiont, color=Bleaching)) +
  geom_bar(position='fill',stat='identity') +
  facet_wrap(~Sample..) +
  xlab('Time Point and Bleaching') +
  ylab('Proportion C and D Clades') +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90,size=10)) +
  scale_x_discrete(labels=c('T1', 'T2 NB', 'T2 B')) +
  scale_fill_manual('Clade',values=c('deepskyblue',  'lightgoldenrod1'), labels=c('C', 'D')) +
  scale_color_manual(values=c('white', 'black','grey'), labels=c('B', 'NB', 'T1'))

png("symbiont v3.png", width=7, height=7, units = "in", res = 1200)
symb3
dev.off()

#ANOVA
symb.Donly<-subset(symb.dat, Symbiont=="D")

res.aov2<-aov(Abundance ~ Bleaching*Time.point*Tolerance*Sample.., data=symb.Donly)
summary(res.aov2)

Df Sum Sq Mean Sq F value Pr(>F)  
Bleaching                     2 0.2803 0.14013   2.560 0.0983 .
Tolerance                     1 0.1961 0.19606   3.581 0.0706 .
Sample..                      1 0.0430 0.04299   0.785 0.3843  
Bleaching:Tolerance           2 0.1597 0.07985   1.458 0.2525  
Bleaching:Sample..            2 0.0152 0.00762   0.139 0.8708  
Tolerance:Sample..            1 0.0737 0.07371   1.346 0.2574  
Bleaching:Tolerance:Sample..  2 0.2711 0.13555   2.476 0.1053  
Residuals                    24 1.3140 0.05475  