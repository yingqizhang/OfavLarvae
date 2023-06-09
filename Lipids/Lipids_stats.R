##2019 Florida Spawning - Gametes

library(car)
library(ggplot2)
library(cowplot)

gametes=read.csv("2019_gametes.csv")
gametes

hist(gametes$Lipid_gamete)

#test for normality with Shapiro-Wilks test to have official consideration of normal
#p<0.05 = not a normal dist.
shapiro.test(gametes$Phospholip)
gametes$LogWax = log(gametes$WAX)
# only WAX needs to be log transformed

#to test for factors one-way ANOVA
gametes_anova = lm(gametes$Total~gametes$site)
anova(gametes_anova)
# Analysis of Variance Table
# 
# Response: gametes$Total
# Df    Sum Sq    Mean Sq F value Pr(>F)
# gametes$site  1 0.0015884 0.00158839  2.6649 0.1412
# Residuals     8 0.0047683 0.00059604


#Testing Lipid Classes

gametes_anova = lm(gametes$LogWax~gametes$site)
anova(gametes_anova)
# Analysis of Variance Table
# 
# Response: gametes$LogWax
# Df  Sum Sq  Mean Sq F value  Pr(>F)  
# gametes$site  1 0.29339 0.293385   6.266 0.03676 *
#   Residuals     8 0.37458 0.046822                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


gametes_anova = lm(gametes$TAG~gametes$site)
anova(gametes_anova)
# Analysis of Variance Table
# 
# Response: gametes$TAG
# Df  Sum Sq Mean Sq F value Pr(>F)
# gametes$site  1  0.6409 0.64088  0.4323 0.5293
# Residuals     8 11.8598 1.48248 


gametes_anova = lm(gametes$Phospholip~gametes$site)
anova(gametes_anova)
# Analysis of Variance Table
# 
# Response: gametes$Phospholip
# Df Sum Sq Mean Sq F value  Pr(>F)  
# gametes$site  1 44.321  44.321  6.0686 0.03911 *
#   Residuals     8 58.427   7.303                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##2021 Florida Spawning - Gametes and Larvae

gametes=read.csv("2021_gametes.csv")
gametes

hist(gametes$Lipid_gamete)

#test for normality with Shapiro-Wilks test to have official consideration of normal
#p<0.05 = not a normal dist.
shapiro.test(gametes$Lipid_gamete)
## p-value = 0.07303, so normal (yayyy!)

#to test for factors two-way ANOVA
gametes_anova = lm(gametes$Lipid_gamete~gametes$Date*gametes$Site)
anova(gametes_anova)
##not significantly different Response: gametes$Lipid_gamete
#Df     Sum Sq    Mean Sq F value Pr(>F)
#gametes$Date  1 0.00000671 0.00000671  0.0101 0.9246
#gametes$Site  1 0.00015734 0.00015734  0.2378 0.6513
#Residuals     4 0.00264619 0.00066155

#gamete size
hist(gametes$Avg..gamete.volume)
#test for normality with Shapiro-Wilks test to have official consideration of normal
#p<0.05 = not a normal dist.
shapiro.test(gametes$Avg..gamete.volume)
##  p-value = 0.06932, normal!

#to test for factors two-way ANOVA  -> are there differences between the batches?
gametes_anova = lm(gametes$Avg..gamete.volume~gametes$Date*gametes$Site)
anova(gametes_anova)
# sign difference by date
#Analysis of Variance Table
#Response: gametes$Avg..gamete.volume
#              Df   Sum Sq    Mean Sq  F value  Pr(>F)  
#gametes$Date  1 7.7657e-05 7.7657e-05 14.9122 0.01812 *
#gametes$Site  1 4.9910e-06 4.9910e-06  0.9584 0.38303  
#Residuals     4 2.0830e-05 5.2080e-06            


#Testing Lipid Classes

hist(gametes$WAX)
#test for normality with Shapiro-Wilks test to have official consideration of normal
shapiro.test(gametes$WAX)
##  p-value = 0.8741, normal!

#to test for factors two-way ANOVA  -> are there differences between the batches?
gametes_anova = lm(gametes$WAX~gametes$Date*gametes$Site)
anova(gametes_anova)
# Analysis of Variance Table
# 
# Response: gametes$WAX
# Df     Sum Sq    Mean Sq F value Pr(>F)
# gametes$Date  1 1.1048e-05 1.1048e-05  2.3706 0.1985
# gametes$Site  1 7.0974e-06 7.0974e-06  1.5229 0.2847
# Residuals     4 1.8642e-05 4.6606e-06 

hist(gametes$TAG)
#test for normality with Shapiro-Wilks test to have official consideration of normal
#p<0.05 = not a normal dist.
shapiro.test(gametes$TAG)
##  p-value = 0.1012, normal!

#to test for factors two-way ANOVA  -> are there differences between the batches?
gametes_anova = lm(gametes$TAG~gametes$Date*gametes$Site)
anova(gametes_anova)
# Analysis of Variance Table
# 
# Response: gametes$TAG
# Df     Sum Sq    Mean Sq F value  Pr(>F)  
# gametes$Date  1 1.8423e-10 1.8423e-10  0.8504 0.40863  
# gametes$Site  1 1.1310e-09 1.1310e-09  5.2206 0.08433 .
# Residuals     4 8.6655e-10 2.1664e-10 

hist(gametes$ST)
#test for normality with Shapiro-Wilks test to have official consideration of normal
#p<0.05 = not a normal dist.
shapiro.test(gametes$ST)
##  p-value = 0.2916, normal!

#to test for factors two-way ANOVA  -> are there differences between the batches?
gametes_anova = lm(gametes$ST~gametes$Date*gametes$Site)
anova(gametes_anova)
# not sign.

hist(gametes$AMPL)
#test for normality with Shapiro-Wilks test to have official consideration of normal
#p<0.05 = not a normal dist.
shapiro.test(gametes$AMPL)
##  p-value = 0.06655, normal!!

#to test for factors two-way ANOVA  -> are there differences between the batches?
gametes_anova = lm(gametes$AMPL~gametes$Date*gametes$Site)
anova(gametes_anova)
# not sign.

shapiro.test(gametes$PL_ug)
#p-value = 0.6341
gametes_anova = lm(gametes$PL_ug~gametes$Date*gametes$Site)
anova(gametes_anova)
# Analysis of Variance Table
# 
# Response: gametes$PL_ug
# Df Sum Sq Mean Sq F value Pr(>F)
# gametes$Date  1 0.4069 0.40690  0.3065 0.6093
# gametes$Site  1 0.4410 0.44096  0.3321 0.5953
# Residuals     4 5.3107 1.32767

#####################
larvae=read.csv("Larvae_YZ.csv")
str(larvae)

hist(larvae$Total_lipid)

#test for normality with Shapiro-Wilks test to have official consideration of normal
#p<0.05 = not a normal dist.
shapiro.test(larvae$Total_lipid)
##  p-value = 2.197e-06, need transformation!
larvae["Total_lipid"][larvae["Total_lipid"] == 0] <- NA
larvae$Total_lipid
larvae <- larvae[complete.cases(larvae), ]
larvae
larvae$Log_Total_lipid = log(larvae$Total_lipid)
larvae$Log_Total_lipid
shapiro.test(larvae$Log_Total_lipid)
# p-value = 0.2316

#to test for factors two-way ANOVA  -> are there differences between the batches?
larvae_anova = lm(larvae$Log_Total_lipid~larvae$Site*larvae$Temperature)
anova(larvae_anova)
# Response: larvae$Log_Total_lipid
# Df  Sum Sq Mean Sq F value Pr(>F)
# larvae$Site                     2  0.9649 0.48244  0.4417 0.6491
# larvae$Temperature              1  0.0002 0.00017  0.0002 0.9902
# larvae$Site:larvae$Temperature  2  4.0395 2.01977  1.8492 0.1833
# Residuals                      20 21.8445 1.09223               

#Testing Lipid Classes

larv = read.csv("Larvae_YZ.csv")
hist(larv$WAX)
hist(larv$Phospholipids)
hist(larv$TAG)
hist(larv$ST)
hist(larv$AMPL)

#test for normality with Shapiro-Wilks test to have official consideration of normal
#p<0.05 = not a normal dist.
shapiro.test(larv$WAX)
##  p-value = 0.1024

shapiro.test(larv$Phospholipids)
## p-value=0.0004027, needs transformation
larv$Log_phospholipid = log(larv$Phospholipids)
shapiro.test(larv$Log_phospholipid)
# p = 0.3459

shapiro.test(larv$TAG)
larv$Log_TAG = log(larv$TAG)
shapiro.test(larv$Log_TAG)
# Too many zero values... nvmd
# Same with ST and AMPL

# WAX
WAX_anova = lm(larv$WAX~larv$Site*larv$Temperature)
anova(WAX_anova)
# Response: larv$WAX
# Df  Sum Sq Mean Sq F value Pr(>F)
# larv$Site                   2   4.114  2.0568  0.2213 0.8031
# larv$Temperature            1   3.465  3.4646  0.3728 0.5472
# larv$Site:larv$Temperature  2   6.741  3.3704  0.3626 0.6996
# Residuals                  24 223.073  9.2947 

#Phospholipids
Phos_anova = lm(larv$Log_phospholipid~larv$Site*larv$Temperature)
anova(Phos_anova)
# Response: larv$Log_phospholipid
# Df  Sum Sq Mean Sq F value Pr(>F)
# larv$Site                   2  2.1669 1.08347  0.9841 0.3884
# larv$Temperature            1  0.5934 0.59339  0.5389 0.4700
# larv$Site:larv$Temperature  2  0.2033 0.10166  0.0923 0.9121
# Residuals                  24 26.4247 1.10103
