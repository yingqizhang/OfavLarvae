library(survival) # for survival analysis 
library(vegan) # for plotting ellipses in principal coordinates analysis plots
library(corrr) # for correlations
library(MCMCglmm) # for mcmcglmm stats
library(car) # for data transformations 
library(nlme) #for lme
library(MASS) # for stepAIC
library(coxme) # for mixed effects cox model
library(tidyverse) # for data wrangling and visualization
library(ggridges) # for ridge plots
library(reshape2) # for melt
library(corrplot) # for correlations
library(summarytools) # for dfSummary
library(ggplot2)
library(RColorBrewer)
#install.packages("survminer",dependencies=TRUE)
library(survminer)
#install.packages("coxme")
library(coxme)

source('~/summarySE.R')


################### Kaplan-Meier survival curves

tod=read.csv("~/TOD_2019.csv")

str(tod)	

tod$Bulk=as.factor(tod$Bulk)
tod$Rep=as.factor(tod$Rep)
tod$Origin=as.factor(tod$Origin)
tod$Trmt=as.factor(tod$Trmt)

summary(tod)


todCC<-tod[complete.cases(tod$Dead),] #Remove rows that have NA for the column Dead

todCC$Origin<-relevel(todCC$Origin,ref="Offshore")## picked origin with best survival as reference (want to know how much worse others are)
todCC$Trmt<-relevel(todCC$Trmt,ref="Control")## picked trmt with best survival as reference

died = Surv(todCC$TOD, todCC$Dead)

#stats

#nw1k$Family <- relevel(nw1k$Family,ref="27")

surT = survfit(died ~ Trmt, data=todCC)
surT

1-(119/510) # 77% survival of ctrl
1-(274/520) # 47% survival of heats

surO = survfit(died ~ Origin, data=todCC)
surO

1-(122/360) # 66.1% survival of offshore
1-(118/350) # 66.3% survival of Cross
1-(153/320) # 52.2% survival of Inshore

mAll<-coxme(died ~ Origin*Trmt+(1|Rep),data=todCC) #does survival differ among inshore/offshore
(print(mAll))

# Cox mixed-effects model fit by maximum likelihood
# Data: todCC
# events, n = 393, 1030
# Iterations= 19 100 
# NULL Integrated    Fitted
# Log-likelihood -2639.711  -2563.713 -2500.116
# 
# Chisq    df p    AIC    BIC
# Integrated loglik 152.00  6.00 0 140.00 116.15
# Penalized loglik 279.19 54.56 0 170.08 -46.73
# 
# Model:  died ~ Origin * Trmt + (1 | Rep) 
# Fixed coefficients
# coef exp(coef)  se(coef)     z       p
# OriginCross            -0.01923847 0.9809454 0.3110761 -0.06 9.5e-01
# OriginInshore           0.63663217 1.8901046 0.2919973  2.18 2.9e-02 *
# TrmtHeat                1.13397202 3.1079770 0.2728321  4.16 3.2e-05 *
# OriginCross:TrmtHeat    0.07240922 1.0750952 0.3939036  0.18 8.5e-01
# OriginInshore:TrmtHeat -0.10260012 0.9024878 0.3799777 -0.27 7.9e-01
# 
# Random effects
# Group Variable  Std Dev   Variance 
# Rep   Intercept 0.5540231 0.3069416
# Cox mixed-effects model fit by maximum likelihood
# Data: todCC
# events, n = 393, 1030
# Iterations= 19 100 
# NULL Integrated    Fitted
# Log-likelihood -2639.711  -2563.713 -2500.116
# 
# Chisq    df p    AIC    BIC
# Integrated loglik 152.00  6.00 0 140.00 116.15
# Penalized loglik 279.19 54.56 0 170.08 -46.73
# 
# Model:  died ~ Origin * Trmt + (1 | Rep) 
# Fixed coefficients
# coef exp(coef)  se(coef)     z       p
# OriginCross            -0.01923847 0.9809454 0.3110761 -0.06 9.5e-01
# OriginInshore           0.63663217 1.8901046 0.2919973  2.18 2.9e-02 ***
# TrmtHeat                1.13397202 3.1079770 0.2728321  4.16 3.2e-05 ***
# OriginCross:TrmtHeat    0.07240922 1.0750952 0.3939036  0.18 8.5e-01
# OriginInshore:TrmtHeat -0.10260012 0.9024878 0.3799777 -0.27 7.9e-01
# 
# Random effects
# Group Variable  Std Dev   Variance 
# Rep   Intercept 0.5540231 0.3069416

exp(ranef(mAll)[[1]])

# Hazard ratio?

# 1         2         3         4         5         6         7         8         9 
# 1.6859013 0.9370296 1.0094936 1.2070739 1.2508946 0.6578139 0.7644960 1.0094936 0.6395354 
# 10        11        12        13        14        16        17        19        20 
# 0.6395354 3.0246452 1.0094936 1.6175692 1.4469134 0.5290256 0.5290256 1.2149646 1.4567126 
# 21        22        23        24        25        26        27        28        29 
# 1.0175415 1.0019871 0.8349139 1.0553126 1.0124158 1.0019871 0.8167130 0.8349139 0.6577916 
# 30        31        32        33        34        35        36        37        38 
# 0.8349139 2.7274904 2.0298605 0.6577916 0.8277227 0.8349139 0.6577916 1.4696224 0.6615201 
# 39        40        41        42        43        44        45        46        47 
# 0.8401707 1.2987027 0.8329965 0.8329965 0.8329965 1.2391676 0.8220104 1.0461951 0.6615201 
# 48        50        51        52        53        54        55        57        58 
# 0.6615201 1.0093042 0.8220104 6.1011140 0.8293822 0.8329965 1.3260481 0.9814756 1.3805398 
# 59        61        62        63        64        65        66        67        68 
# 1.6283405 0.6259328 0.8203859 1.3295454 0.8298509 0.6451636 1.0101479 0.7253408 0.6491093 
# 69        70        71        72        73        74        75        76        77 
# 1.2093241 0.9012644 1.0574856 1.7059673 1.0552329 1.3183834 1.1017214 0.7640321 0.9607236 
# 78        79        80        81        82        83        84        85        86 
# 0.6319571 1.0456293 0.8878122 0.6414252 1.2906345 1.2419159 1.6171652 1.1419022 1.3701537 
# 87        88        89        90        91        92        93        94        95 
# 1.4344916 0.9142867 0.8573554 0.5179911 1.1218766 1.1243052 1.0127461 0.5968296 1.0568414 
# 96        97        98        99       100       101       102       103       104 
# 1.0221274 0.9708839 0.6725686 0.7742733 0.6725686 1.1473106 0.7876128 1.0127461 0.9537494 
# 105       106       107       108 
# 1.1201712 0.7290084 1.7051084 2.9382632 

sr=exp(ranef(mAll)[[1]])
sr=stack(sr)
names(sr)=c("HR","ID")

survival=read.csv("~/Survival_2019.csv")
info=survival[,c(1:4)]

sr=merge(sr,info,by="ID")

setwd("~/Acute_stress")
write.table(x=sr, file='Hazard_ratio.csv', append=FALSE,  quote=FALSE, row.names=FALSE, col.names=TRUE, sep=',')
