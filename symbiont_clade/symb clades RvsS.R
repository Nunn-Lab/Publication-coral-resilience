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


# Logistic Regressin Analysis to determine if the proportion of Durusdinium:Cladocopium influence survivability
# Input the data from symbiont proportions into dataframe
data <- data.frame(
  Coral_Genet_ID = c(5, 20, 25, 28, 55, 73, 10, 13, 21, 26, 66, 74),
  Susceptibility = c("S", "S", "S", "S", "S", "S", "R", "R", "R", "R", "R", "R"),
  Proportion_Cladocopium = c(0.196, 0.142, 0, 0.081, 0, NA, 0, 1, 0.843, 0, 0, 0),
  Proportion_Durusdinium = c(0.804, 0.858, 1, 0.919, 1, NA, 1, 0, 0.157, 1, 1, 1)
)

# Remove rows with missing values
data <- na.omit(data)

# Convert 'Susceptibility' to a binary variable (1 for 'S', 0 for 'R')
data$Susceptibility_Binary <- ifelse(data$Susceptibility == "S", 1, 0)

# Perform logistic regression
model <- glm(Susceptibility_Binary ~ Proportion_Durusdinium, data = data, family = binomial)

# Display the summary of the model
summary(model)


Call:
glm(formula = Susceptibility_Binary ~ Proportion_Durusdinium, family = binomial, data = data)

Deviance Residuals: 
   Min      1Q  Median      3Q     Max  
-1.278  -0.886  -0.533   1.051   1.210  

Coefficients:
                        Estimate Std. Error z value Pr(>|z|)  
(Intercept)              -2.1581     2.2870  -0.944   0.345  
Proportion_Durusdinium    2.3817     2.5040   0.951   0.342  

(Dispersion parameter for binomial family taken to be 1)

    Null deviance: 15.162  on 10  degrees of freedom
Residual deviance: 13.873  on  9  degrees of freedom
AIC: 17.873

Number of Fisher Scoring iterations: 6
