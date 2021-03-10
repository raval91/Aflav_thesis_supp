#### FULL ANALYSIS OF HETEROTHERMY INDEX ####

rm(list = ls())

#setwd("~/Documents/Thesis_analysis")

library(reshape2)
library(lme4)
library(aod)
library(bestNormalize)
library(SNPRelate)
library(ggplot2)
library(ggpubr)
library(MASS)
library(car)
library(rstatix)
library(dplyr)

#### DATA ORGANISATION ####
#read in the data
phenotypic_data_full <- read.table("data/phenotypic_data/A.flavicollis_phenotypic_data_clean.txt", header = T, stringsAsFactors = F, sep = "\t")
phenotypic_data_full$Trapping_year <- as.factor(phenotypic_data_full$Trapping_year)

###create long format table for HI data
#melt df by season
HI_data <- melt(phenotypic_data_full, id.vars = c("ID", "merged_samples", "Sex", "Age", "mb1", "mb2", "Trapping_year", "Birth_year", "HI1f", "TS1"), 
                 measure.vars = c("BMR_autumn", "minBMR_winter.spring"))
HI_data$variable <- as.character(HI_data$variable)
HI_data$variable[HI_data$variable == "BMR_autumn"] <- "Autumn"
HI_data$variable[HI_data$variable == "minBMR_winter.spring"] <- "Winter/Spring"
HI_data$Sex[HI_data$Sex == "m"] <- "Male"
HI_data$Sex[HI_data$Sex == "f"] <- "Female"
colnames(HI_data)[c(11:12)] <- c("Season", "BMR_min")

BMR_data <- melt(phenotypic_data_full, id.vars = c("merged_samples"), measure.vars = c("mb_BMR_autumn", "mb_BMR_winter.spring_minBMR"))
colnames(BMR_data)[2:3] <- c("Season", "mb")
BMR_data$Season <- as.character(BMR_data$Season)
BMR_data$Season[BMR_data$Season == unique(BMR_data$Season)[1]] <- "Autumn"
BMR_data$Season[BMR_data$Season == unique(BMR_data$Season)[2]] <- "Winter/Spring"
HI_data$mb_bmr <- BMR_data$mb
HI_data$ID <- as.character(HI_data$ID)
rm(phenotypic_data_full)

#turn questionable ages into proper ages
#perhaps remove these from the analysis altogether??
#how much will it affect the analysis? best to remove as unsure of their age!!
HI_data$Age[HI_data$Age == "ad?"] <- NA
HI_data$Age[HI_data$Age == "juv?"] <- NA

#TRANSFORM THE HETERITHERMY INDICES (HI1f) 
bestNormalize(HI_data$HI1f)#which transfrmation? orderNorm transformation
#trans <- orderNorm(jitter(HI_data$HI1f))#random noise added to ensure no ties in data resulting in true gaussian distribution
trans <- orderNorm(HI_data$HI1f)
HI_data$HI_orderNorm <- trans$x.t

#calculate the proportion of mass lost during torpor
HI_data$mass_loss <- HI_data$mb1 - HI_data$mb2#mass lost
HI_data$proportion_mass_loss <- HI_data$mass_loss/HI_data$mb1#proportion mass lost

#adjust HI for mb
HI_data$adjusted_HI <- predict(lmer(HI_orderNorm ~ mb1 + (1 | ID), data = HI_data))

#calculate the mean for HI1f (untransformed HI data)
mean_HI1f <- melt(tapply(HI_data$HI1f, HI_data$ID, mean), value.name = "mean_HI1f")
colnames(mean_HI1f)[1] <- "ID"
HI_data <- merge(HI_data, mean_HI1f, by = "ID")
rm(mean_HI1f)

#function to calculate standard error
std <- function(x){
  sd(x)/sqrt(length(x))
}

#calculate std err
std_errHIorderNorm <- data.frame(ID = names(tapply(HI_data$HI_orderNorm, HI_data$ID, std)),
                                 std_errHIorderNorm = tapply(HI_data$HI_orderNorm, HI_data$ID, std))

std_errHI1f <- data.frame(ID = names(tapply(HI_data$HI1f, HI_data$ID, std)),
                          std_errHI1f = tapply(HI_data$HI1f, HI_data$ID, std))

#merge phenotypic data and standard error data by ID
# HI_data <- Reduce(function(x,y) merge(x = x, y = y, by = "ID"), list(HI_data, std_errHIorderNorm))
HI_data <- Reduce(function(x,y) merge(x = x, y = y, by = "ID"), list(HI_data, std_errHI1f))
rm(std_errHIorderNorm, std_errHI1f)

##order the data by box cox transformed HI and add new column for the ordered index
HI_data <- HI_data[order(HI_data$mean_HI1f),]

#function to greate ordereg group IDs
grpid <- function(x) {
  match(x, unique(x))
}

#create unique group IDs
HI_data <- HI_data %>% mutate(ordered_sample = group_indices(., ID) %>% grpid())

# AVERAGE DATA TO ACCOUNT FOR PSEUDO-SAMPLING 
#this is to check the differences between sex/years/age for averaged data 
#as the data is pseudoreplicated. some individuals have been masured more than once.
#the same trends are found BUT the differences are less significant as the degrees
#of freedom are no longer inflated. this is still a problem for the heritability models 
#though. one solution is to work only with the first "pseudo-replicate" to ensure the true 
#degrees of freedom are taken into account.

#average the data
ave_HI <- data.frame(HI_ave = tapply(HI_data$HI1f, HI_data$merged_samples, mean))

#order bmr data by sample
HI_data <- HI_data[order(HI_data$merged_samples),]

#add relevant data to df
ave_HI$Sex <- HI_data[match(unique(HI_data$merged_samples), HI_data$merged_samples),]$Sex #sex
ave_HI$merged_samples <- row.names(ave_HI) #sample id
ave_HI$Age <- HI_data[match(unique(HI_data$merged_samples), HI_data$merged_samples),]$Age #Age
ave_HI$Trapping_year <- HI_data[!duplicated(HI_data$merged_samples),]$Trapping_year#Trapping year

trans <- orderNorm(ave_HI$HI_ave)#transform averaged HI data
ave_HI$HI_ave_gaus <- trans$x.t

ave_HI$mean_mb1 <- as.numeric(tapply(HI_data$mb1, HI_data$merged_samples, mean, na.rm = T))#mb1. no need to transform as conforms gaus dist
ave_HI$mean_mb2 <- as.numeric(tapply(HI_data$mb2, HI_data$merged_samples, mean, na.rm = T))#mb2. no need to transform as conforms gaus dist
ave_HI$mean_mass_loss <- ave_HI$mean_mb1 - ave_HI$mean_mb2#mass lost
ave_HI$proportion_mean_mass_loss <- ave_HI$mean_mass_loss/ave_HI$mean_mb1#proportion mass lost

#### SUMMARY AND ANALYSIS OF HETEROTHERMY DATA ####

#how meany measures of each mouse?
replicates <- data.frame(table(HI_data[!is.na(HI_data$BMR_min),]$merged_samples))
hist(replicates$Freq)
summary(replicates$Freq)

#how many mice phenotyped by...
table(HI_data[which(!duplicated(HI_data$ID) == T), ]$Sex)#sex?
table(HI_data[which(!duplicated(HI_data$ID) == T), ]$Age)#age?
table(HI_data[which(!duplicated(HI_data$ID) == T), ]$Trapping_year)#trapping year?

#Summary of...
summary(HI_data$HI1f)#untransformed Hi
sd(HI_data$HI1f)
summary(HI_data$TS1)#Tb_mode - Tb_torpor
sd(HI_data$TS1)

#hi vs mb
summary(lm(HI_data$HI_orderNorm ~ HI_data$mb1))
cor.test(ave_HI$HI_ave, ave_HI$mean_mb1)

#adjust HI for mb in averaged data
ave_HI$adjusted_HI <- predict(lm(ave_HI$HI_ave_gaus ~ ave_HI$mean_mb1))

#adjust bmr for mb (full data frame)
ave_HI$mean_bmr_min <- tapply(HI_data$BMR_min, HI_data$merged_samples, mean, na.rm = T)#calculate mean bmr_min
ave_HI$mean_mb_bmr <- tapply(HI_data$mb_bmr, HI_data$merged_samples, mean, na.rm = T)#calculate mean mb_bmr
ave_HI$adjusted_bmr <- predict(lm(ave_HI$mean_bmr_min ~ ave_HI$mean_mb_bmr))#adjust

#any difference in adjusted HI by sex and age or trapping year???
#age was removed due to non-significance. will not improve model fits!!!
t.test(ave_HI$adjusted_HI ~ ave_HI$Sex)# YES... BUT SEX ONLY. MALES HAVE LOWER HI.
summary(lm(ave_HI$HI_ave_gaus ~ ave_HI$Trapping_year))#NO ... not mb adjusted
summary(lm(ave_HI$adjusted_HI ~ ave_HI$Trapping_year))#NO ... mb adjusted

##HOW DOES THE PROPORTION OF MASS LOST DURING FASTING DIFFER BETWEEN SEXES? 
#higher HI means less mass lost (2nd and 3rd models). makes sense... 
#males need to conserve more mass than females as have a higher cost to need a higher HI. 
#this makes no comment on body condition however. 
#averaged data
# summary(glm(ave_HI$proportion_mean_mass_loss ~ ave_HI$HI_ave + ave_HI$Sex, family = "binomial"))#NO DIFFERENCE IN COST FOR MALES OR FEMALES!!!
wilcox.test(ave_HI$proportion_mean_mass_loss ~ ave_HI$Sex)#NO DIFFERENCE IN COST FOR MALES OR FEMALES!!! glm is overkill... simple wilcox is better.

#DO MALES HAVE HIGHER MAINTENANCE COSTS THAN FEMALES?

trans <- orderNorm(ave_HI$mean_bmr_min)#transform first as BMR_min is not normal
ave_HI$mean_bmr_min_gaus <- trans$x.t
summary(lm(ave_HI$adjusted_bmr ~ ave_HI$Sex))#males have higher maintenance costs even when mass adjusted


#WHAT IS THE RELATIONSHIP BETWEEN HI AND BMR?
#seems as though we have the opposite trend to boratynski et al. HI and BMR are inversly correlated.
#summer
hi_reduced <- HI_data[!is.na(HI_data$BMR_min),]
# trans <- orderNorm(jitter(hi_reduced$BMR_min))#transform first as BMR_min is not normal
trans <- orderNorm(hi_reduced$BMR_min)#transform first as BMR_min is not normal
hi_reduced$mean_bmr_min_gaus <- trans$x.t
#autumn
# HI_autumn <- HI_data[HI_data$Season == "Autumn",]
# trans <- orderNorm(jitter(HI_autumn$BMR_min))
# HI_autumn$mean_bmr_min_gaus <- trans$x.t

#create a df of mass adjusted data
adjusted_HI <- predict(lmer(HI_orderNorm ~ mb1 + Sex + (1|ID), data = hi_reduced))
adjusted_bmr <- predict(lm(mean_bmr_min_gaus ~ mb1, data = hi_reduced))
adjusted_data <- data.frame(ID = hi_reduced$ID, Sex = hi_reduced$Sex, hi_reduced$mb1, adjusted_HI, adjusted_bmr)
rm(adjusted_HI, adjusted_bmr)

#mice with higher HI have a lower BMR. these are primarily females!
summary(lm(adjusted_data$adjusted_HI ~ adjusted_data$adjusted_bmr))#
cor.test(adjusted_data$adjusted_HI, adjusted_data$adjusted_bmr)

#### PLOTTING ####
#plot the variation in average HI
ggplot(data = HI_data, aes(x = ordered_sample , y = mean_HI1f, colour = TS1)) +
  geom_point() +
  scale_colour_gradient(low = "#6495ED85", high = "#C1403D85",
                        guide = guide_colourbar(title.position = "top",
                                                ticks = T, frame.colour = "black",
                                                title.hjust = -0.5, ticks.colour = "black")) +
  geom_pointrange(data = HI_data, aes(ymin = mean_HI1f - std_errHI1f,
                                       ymax = mean_HI1f + std_errHI1f)) +
  labs(colour = expression(paste("Deviation from normothermy", " ("^"o", "C)")), 
       y = expression(paste("Hi"[mu], " ("^"o", "C)")),
       x = "Sample") +
  theme_bw() + 
  theme(
    legend.position = c(0.25, 0.66), 
    legend.direction = "horizontal", 
    axis.text.x = element_blank())

#relationship between HI and mb
ggplot(HI_data, aes(x = mb1, y = HI_orderNorm, colour = Sex)) +
  geom_point() +
  scale_colour_manual(values = c("#6495ED85", "#C1403D85")) +
  ylab(expression(paste("Order-Norm Hi", " ("^"o", "C)"))) +
  xlab(expression(paste(m[b], " (g)"))) +
  stat_smooth(method = "lm", colour = "#660000", se = F) +
  theme_bw() +
  theme(legend.position = c(0.75, 0.8))

##relationship between adjusted HI and adjusted BMRmin
ggplot(adjusted_data, aes(x = adjusted_bmr, y = adjusted_HI, colour = Sex)) +
  geom_point() +
  scale_colour_manual(values = c("#6495ED85", "#C1403D85")) +
  ylab(expression(paste("Residual Order-Norm Hi", " ("^"o", "C)"))) +
  xlab(expression(paste("Residual ", BMR[min], " (mW)"))) +
  stat_smooth(method = "lm", colour = "#660000", se = F) +
  theme_bw() +
  theme(legend.position = c(0.75, 0.8))

#### IS HI REPEATABLE? #### 
#Age and BMR were nonsignificant terms so removed. 
#Trapping year and Birth year caused singilar fits in the model as not enough variation 
#in the data so these terms were removed. Although sex is non-significant... keep it in
#and report repeatability and heritability for both. 

#create male and female dataframes
HI_male <- HI_data[HI_data$Sex == "Male",]
HI_female <- HI_data[HI_data$Sex == "Female",]
#repeatability using REML
repeatability_mod <- lmer(HI_orderNorm ~ Sex + mb1 + BMR_min + 
                            (1 | Birth_year) + 
                            (1 | ID), 
                          data = HI_data)
summary(repeatability_mod)
Anova(repeatability_mod)

repeatability_mod1 <- lmer(HI_orderNorm ~ Sex + mb1 + BMR_min +
                            (1 | ID), 
                          data = HI_data)
summary(repeatability_mod1)
Anova(repeatability_mod1)

repeatability_mod2 <- lmer(HI_orderNorm ~ Sex  + mb1 + 
                             (1 | ID), 
                           data = HI_data)
summary(repeatability_mod2)
Anova(repeatability_mod2)

repeatability_mod3 <- lmer(HI_orderNorm ~ mb1 + 
                             (1 | ID), 
                           data = HI_data)
summary(repeatability_mod3)
Anova(repeatability_mod3)

repeatability_mod_males <- lmer(HI_orderNorm ~ mb1 + BMR_min +
                                                        (1 | ID), 
                                                      data = HI_male)
summary(repeatability_mod_males)
Anova(repeatability_mod_males)

repeatability_mod_males1 <- lmer(HI_orderNorm ~ mb1 +
                                                        (1 | ID), 
                                                      data = HI_male)
summary(repeatability_mod_males1)
Anova(repeatability_mod_males1)

repeatability_mod_females <- lmer(HI_orderNorm ~ mb1 + BMR_min +
                                  (1 | ID), 
                                data = HI_female)
summary(repeatability_mod_females)
Anova(repeatability_mod_females)

repeatability_mod_females1 <- lmer(HI_orderNorm ~ mb1 +
                                  (1 | ID), 
                                data = HI_female)
summary(repeatability_mod_females1)
Anova(repeatability_mod_females1)

#anova(repeatability_mod1, repeatability_mod2, repeatability_mod3)

#model diagnostics
plot(repeatability_mod1)
plot(repeatability_mod2)
plot(repeatability_mod3)
plot(repeatability_mod_males)
plot(repeatability_mod_males1)
plot(repeatability_mod_females)
plot(repeatability_mod_females1)
qqPlot(HI_data$HI_orderNorm)
qqPlot(HI_male$HI_orderNorm)
qqPlot(HI_female$HI_orderNorm)
####calculate repeatability
##predcted values from the model with fixed effects
#overall
repeatability_predicted1 <- predict(lm(HI_data$HI_orderNorm ~ HI_data$mb1 + HI_data$Sex + HI_data$BMR_min), interval = "confidence")# mb + sex + bmr
repeatability_predicted2 <- predict(lm(HI_data$HI_orderNorm ~ HI_data$mb1 + HI_data$Sex), interval = "confidence")# sex + mb
repeatability_predicted3 <- predict(lm(HI_data$HI_orderNorm ~ HI_data$mb1), interval = "confidence")# mb

#males
repeatability_predicted_m <- predict(lm(HI_male$HI_orderNorm ~ HI_male$mb1 + HI_male$BMR_min), interval = "confidence")# mb + bmr
repeatability_predicted_m1 <- predict(lm(HI_male$HI_orderNorm ~ HI_male$mb1), interval = "confidence")# mb

#females
repeatability_predicted_f <- predict(lm(HI_female$HI_orderNorm ~ HI_female$mb1 + HI_female$BMR_min), interval = "confidence")# mb + bmr
repeatability_predicted_f1 <- predict(lm(HI_female$HI_orderNorm ~ HI_female$mb1), interval = "confidence")# mb


##extract variance components
#overall
repeatability_var1 <- as.data.frame(VarCorr(repeatability_mod1), comp = ("Variance"))# sex + mb + bmr
repeatability_var2 <- as.data.frame(VarCorr(repeatability_mod2), comp = ("Variance"))# sex + mb
repeatability_var3 <- as.data.frame(VarCorr(repeatability_mod3), comp = ("Variance"))# mb

#males
repeatability_var_m <- as.data.frame(VarCorr(repeatability_mod_males), comp = ("Variance"))# mb + bmr
repeatability_var_m1 <- as.data.frame(VarCorr(repeatability_mod_males1), comp = ("Variance"))# mb

#females
repeatability_var_f <- as.data.frame(VarCorr(repeatability_mod_females), comp = ("Variance"))# mb + bmr
repeatability_var_f1 <- as.data.frame(VarCorr(repeatability_mod_females1), comp = ("Variance"))# mb

##estimate repeatability
#overall
repeatability1 <- repeatability_var1$vcov[1]/(sum(repeatability_var1$vcov) + var(repeatability_predicted1[,1]))# sex + mb + bmr
repeatability2 <- repeatability_var2$vcov[1]/(sum(repeatability_var2$vcov) + var(repeatability_predicted2[,1]))# sex + mb
repeatability3 <- repeatability_var3$vcov[1]/(sum(repeatability_var3$vcov) + var(repeatability_predicted3[,1]))# mb

repeatability1# sex + mb + bmr... HI is highly repeatable!!! 0.75
repeatability2# sex + mb... HI is highly repeatable!!! 0.75
repeatability3# mb

#males
repeatability_males <- repeatability_var_m$vcov[1]/(sum(repeatability_var_m$vcov) + var(repeatability_predicted_m[,1]))# mb + bmr
repeatability_males1 <- repeatability_var_m1$vcov[1]/(sum(repeatability_var_m1$vcov) + var(repeatability_predicted_m1[,1]))# mb

repeatability_males
repeatability_males1

#females
repeatability_females <- repeatability_var_f$vcov[1]/(sum(repeatability_var_f$vcov) + var(repeatability_predicted_f[,1]))# mb + bmr
repeatability_females1 <- repeatability_var_f1$vcov[1]/(sum(repeatability_var_f1$vcov) + var(repeatability_predicted_f1[,1]))# mb

repeatability_females
repeatability_females1

##calculate confidence intervals of the fixed effects variance components (sd^2)
#overall
rep_confint1 <- confint(repeatability_mod1, oldNames = F, method = "boot")^2# sex + mb + bmr
rep_confint2 <- confint(repeatability_mod2, oldNames = F, method = "boot")^2# sex + mb
rep_confint3 <- confint(repeatability_mod3, oldNames = F, method = "boot")^2# mb

#males
rep_confint_m <- confint(repeatability_mod_males, oldNames = F, method = "boot")^2# mb + bmr
rep_confint_m1 <- confint(repeatability_mod_males1, oldNames = F, method = "boot")^2# mb

#females
rep_confint_f <- confint(repeatability_mod_females, oldNames = F, method = "boot")^2# mb + bmr
rep_confint_f1 <- confint(repeatability_mod_females1, oldNames = F, method = "boot")^2# mb

###calculate the confidence intervals for heritability
#uses confidence intervals... see how it differs with prediction intervals instead
##overall
# sex + mb + bmr
rep_confint1[1,1]/(sum(rep_confint1[1:2,1]) + var(repeatability_predicted1[,2]))#lower
rep_confint1[1,2]/(sum(rep_confint1[1:2,2]) + var(repeatability_predicted1[,3]))#upper

#sex + mb
rep_confint2[1,1]/(sum(rep_confint2[1:2,1]) + var(repeatability_predicted2[,2]))#lower
rep_confint2[1,2]/(sum(rep_confint2[1:2,2]) + var(repeatability_predicted2[,3]))#upper

#mb
rep_confint3[1,1]/(sum(rep_confint3[1:2,1]) + var(repeatability_predicted3[,2]))#lower
rep_confint3[1,2]/(sum(rep_confint3[1:2,2]) + var(repeatability_predicted3[,3]))#upper

##males
# mb + bmr
rep_confint_m[1,1]/(sum(rep_confint_m[1:2,1]) + var(repeatability_predicted_m[,2]))#lower
rep_confint_m[1,2]/(sum(rep_confint_m[1:2,2]) + var(repeatability_predicted_m[,3]))#upper

#mb
rep_confint_m1[1,1]/(sum(rep_confint_m1[1:2,1]) + var(repeatability_predicted_m1[,2]))#lower
rep_confint_m1[1,2]/(sum(rep_confint_m1[1:2,2]) + var(repeatability_predicted_m1[,3]))#upper

##females
#mb + bmr
rep_confint_f[1,1]/(sum(rep_confint_f[1:2,1]) + var(repeatability_predicted_f[,2]))#lower
rep_confint_f[1,2]/(sum(rep_confint_f[1:2,2]) + var(repeatability_predicted_f[,3]))#upper

#mb
rep_confint_f1[1,1]/(sum(rep_confint_f1[1:2,1])+ var(repeatability_predicted_f1[,2]))#lower
rep_confint_f1[1,2]/(sum(rep_confint_f1[1:2,2]) + var(repeatability_predicted_f1[,3]))#upper

#### CREATE GRM ####

##generate a grm of only phenotyped mice
genofile0.05 <- snpgdsOpen("data/populations_out/0.05/r2_prune_filtered.gds")#read in the binary data

#calculate IBD for all SNPs
start_time <- Sys.time()
IBD_0.05 <- snpgdsIBDMLE(genofile0.05, sample.id = unique(BMR_data$merged_samples), kinship = T, autosome.only = F, num.thread = 8)
end_time <- Sys.time()
end_time - start_time

#turn the kinship matrtix into a standardised GRM
grm_phenotyped <- IBD_0.05$kinship*2
colnames(grm_phenotyped) <- IBD_0.05$sample.id#ensure cols named as samples
row.names(grm_phenotyped) <- IBD_0.05$sample.id#ensure rows are names as samples
grm_phenotyped <- scale(grm_phenotyped, scale = F)

#create a data framw from the relatedness matrix to plot a heatmap
upper_triangle <- upper.tri(IBD_0.05$kinship)#new logical matrix of only the upper triangle
relatedness.upperTriangle <- grm_phenotyped #take a copy of the original relatedness matrix
grm_phenotyped[!upper_triangle]<-NA#set everything not in upper triangle to NA
relatedness<-na.omit(melt(grm_phenotyped, value.name ="relatedness")) #use melt to reshape the matrix into a along format df
relatedness <- relatedness[order(relatedness$relatedness, decreasing = T),] #sort by descending relatedness
names(relatedness) <- c("IID1", "IID2", "genomic_relatedness")

#plot the relatedness matrix. ONLY WORKS WHEN LINES 428-430 ARE RUN ALSO
ggplot(relatedness, aes(x = IID1, y = IID2, fill = genomic_relatedness)) + 
  geom_tile() +  
  theme(axis.text = element_blank(),
        legend.position = c(0.8, 0.4)) +
  labs(fill = "Genomic relatedness") +
  xlab("Sample") +
  ylab("Sample")

mean(relatedness$genomic_relatedness)
sd(relatedness$genomic_relatedness)
sum(relatedness$genomic_relatedness >= 0.025)

#### CALCULATE HERITABILITY ON MICE WHERE HIGHLY RELATED INDIVIDUALS ARE REMOVED (N = 24) ####
##remove all the individuals related by more than 0.025 (as suggested by GCTA)? 

##which mice have 0.025 >= r < 1?
#subset for mice related by < 0.025
unrelated_mice <- relatedness[relatedness$genomic_relatedness <= 0.025 |
                                relatedness$genomic_relatedness >=0.6,]

#turn the molten df back into a relatedness matrix
grm_unrelated <- dcast(unrelated_mice, IID1~IID2)
row.names(grm_unrelated) <- as.character(grm_unrelated[,1])
grm_unrelated <- na.omit(grm_unrelated[,-1])
grm_unrelated <- data.matrix(grm_unrelated[,colnames(grm_unrelated) %in% row.names(grm_unrelated)])

### CALCULATE HERITABILITY ACCORDING TO VILLERMEREUIL ET AL. 2017 ###

#create a df of unrelated samples only
heritability_data <- HI_data[HI_data$merged_samples %in% colnames(grm_unrelated),]#overall
hist(heritability_data$HI1f)#HIGHLY SKEWED!!
# trans <- orderNorm(jitter(heritability_data$HI1f))#OrderNorm transform
trans <- orderNorm(heritability_data$HI1f)#OrderNorm transform
hist(trans$x.t)
heritability_data$HI_orderNorm <- trans$x.t #replace with transformed data

heritability_data_m <- heritability_data[heritability_data$Sex == "Male",]#males
heritability_data_f <- heritability_data[heritability_data$Sex == "Female",]#females

#hertiability models
heritability_mod1 <- lmer(HI_orderNorm ~ Sex + mb1 + BMR_min +
                             (1 | ID), 
                           data = heritability_data)
summary(heritability_mod1)
Anova(heritability_mod1)

heritability_mod2 <- lmer(HI_orderNorm ~ Sex  + mb1 + 
                             (1 | ID), 
                           data = heritability_data)
summary(heritability_mod2)
Anova(heritability_mod2)

heritability_mod3 <- lmer(HI_orderNorm ~ mb1 + 
                             (1 | ID), 
                           data = heritability_data)
summary(heritability_mod3)
Anova(heritability_mod3)

heritability_mod_males <- lmer(HI_orderNorm ~ mb1 + BMR_min +
                                  (1 | ID), 
                                data = heritability_data_m)
summary(heritability_mod_males)
Anova(heritability_mod_males)

heritability_mod_males1 <- lmer(HI_orderNorm ~ mb1 +
                                   (1 | ID), 
                                 data = heritability_data_m)
summary(heritability_mod_males1)
Anova(heritability_mod_males1)

heritability_mod_females <- lmer(HI_orderNorm ~ mb1 + BMR_min +
                                    (1 | ID), 
                                  data = heritability_data_f)
summary(heritability_mod_females)
Anova(heritability_mod_females)

heritability_mod_females1 <- lmer(HI_orderNorm ~ mb1 +
                                     (1 | ID), 
                                   data = heritability_data_f)
summary(heritability_mod_females1)
Anova(heritability_mod_females1)

#anova(heritability_mod1, heritability_mod2, heritability_mod3)

#model diagnostics
plot(heritability_mod1)
plot(heritability_mod2)
plot(heritability_mod3)
plot(heritability_mod_males)
plot(heritability_mod_males1)
plot(heritability_mod_females)
plot(heritability_mod_females1)
qqPlot(heritability_data$HI_orderNorm)
qqPlot(heritability_data_m$HI_orderNorm)
qqPlot(heritability_data_f$HI_orderNorm)

##predcted values from the model with fixed effects
#overall
heritability_predicted1 <- predict(lm(heritability_data$HI_orderNorm ~ 
                                         heritability_data$mb1 + heritability_data$Sex + heritability_data$BMR_min), 
                                    interval = "confidence")# mb + sex + bmr
heritability_predicted2 <- predict(lm(heritability_data$HI_orderNorm ~ heritability_data$mb1 + heritability_data$Sex), 
                                    interval = "confidence")# sex + mb
heritability_predicted3 <- predict(lm(heritability_data$HI_orderNorm ~ heritability_data$mb1), interval = "confidence")# mb

#males
heritability_predicted_m <- predict(lm(heritability_data_m$HI_orderNorm ~ heritability_data_m$mb1 + heritability_data_m$BMR_min), 
                                     interval = "confidence")# mb + bmr
heritability_predicted_m1 <- predict(lm(heritability_data_m$HI_orderNorm ~ heritability_data_m$mb1), 
                                      interval = "confidence")# mb

#females
heritability_predicted_f <- predict(lm(heritability_data_f$HI_orderNorm ~ heritability_data_f$mb1 + heritability_data_f$BMR_min), 
                                     interval = "confidence")# mb + bmr
heritability_predicted_f1 <- predict(lm(heritability_data_f$HI_orderNorm ~ heritability_data_f$mb1), 
                                      interval = "confidence")# mb

##extract variance components
#overall
heritability_var1 <- as.data.frame(VarCorr(heritability_mod1), comp = ("Variance"))# sex + mb + bmr
heritability_var2 <- as.data.frame(VarCorr(heritability_mod2), comp = ("Variance"))# sex + mb
heritability_var3 <- as.data.frame(VarCorr(heritability_mod3), comp = ("Variance"))# mb

#males
heritability_var_m <- as.data.frame(VarCorr(heritability_mod_males), comp = ("Variance"))# mb + bmr
heritability_var_m1 <- as.data.frame(VarCorr(heritability_mod_males1), comp = ("Variance"))# mb

#females
heritability_var_f <- as.data.frame(VarCorr(heritability_mod_females), comp = ("Variance"))# mb + bmr
heritability_var_f1 <- as.data.frame(VarCorr(heritability_mod_females1), comp = ("Variance"))# mb

##estimate heritability
#overall
heritability1 <- heritability_var1$vcov[1]/(sum(heritability_var1$vcov) + var(heritability_predicted1[,1]))# sex + mb + bmr
heritability2 <- heritability_var2$vcov[1]/(sum(heritability_var2$vcov) + var(heritability_predicted2[,1]))# sex + mb
heritability3 <- heritability_var3$vcov[1]/(sum(heritability_var3$vcov) + var(heritability_predicted3[,1]))# mb

heritability1# sex + mb + bmr... HI is highly repeatable!!! 0.75
heritability2# sex + mb... HI is highly repeatable!!! 0.75
heritability3# mb

#males
heritability_males <- heritability_var_m$vcov[1]/(sum(heritability_var_m$vcov) + var(heritability_predicted_m[,1]))# mb + bmr
heritability_males1 <- heritability_var_m1$vcov[1]/(sum(heritability_var_m1$vcov) + var(heritability_predicted_m1[,1]))# mb

heritability_males
heritability_males1

#females
heritability_females <- heritability_var_f$vcov[1]/(sum(heritability_var_f$vcov) + var(heritability_predicted_f[,1]))# mb + bmr
heritability_females1 <- heritability_var_f1$vcov[1]/(sum(heritability_var_f1$vcov) + var(heritability_predicted_f1[,1]))# mb

heritability_females
heritability_females1

##calculate confidence intervals of the fixed effects variance components (sd^2)
#overall
rep_confint1 <- confint(heritability_mod1, oldNames = F, method = "boot")^2# sex + mb + bmr
rep_confint2 <- confint(heritability_mod2, oldNames = F, method = "boot")^2# sex + mb
rep_confint3 <- confint(heritability_mod3, oldNames = F, method = "boot")^2# mb

#males
rep_confint_m <- confint(heritability_mod_males, oldNames = F, method = "boot")^2# mb + bmr
rep_confint_m1 <- confint(heritability_mod_males1, oldNames = F, method = "boot")^2# mb

#females
rep_confint_f <- confint(heritability_mod_females, oldNames = F, method = "boot")^2# mb + bmr
rep_confint_f1 <- confint(heritability_mod_females1, oldNames = F, method = "boot")^2# mb

###calculate the confidence intervals for heritability
#uses confidence intervals... see how it differs with prediction intervals instead
##overall
# sex + mb + bmr
rep_confint1[1,1]/(sum(rep_confint1[1:2,1]) + var(heritability_predicted1[,2]))#lower
rep_confint1[1,2]/(sum(rep_confint1[1:2,2]) + var(heritability_predicted1[,3]))#upper

#sex + mb
rep_confint2[1,1]/(sum(rep_confint2[1:2,1]) + var(heritability_predicted2[,2]))#lower
rep_confint2[1,2]/(sum(rep_confint2[1:2,2]) + var(heritability_predicted2[,3]))#upper

#mb
rep_confint3[1,1]/(sum(rep_confint3[1:2,1]) + var(heritability_predicted3[,2]))#lower
rep_confint3[1,2]/(sum(rep_confint3[1:2,2]) + var(heritability_predicted3[,3]))#upper

##males
# mb + bmr
rep_confint_m[1,1]/(sum(rep_confint_m[1:2,1]) + var(heritability_predicted_m[,2]))#lower
rep_confint_m[1,2]/(sum(rep_confint_m[1:2,2]) + var(heritability_predicted_m[,3]))#upper

#mb
rep_confint_m1[1,1]/(sum(rep_confint_m1[1:2,1]) + var(heritability_predicted_m1[,2]))#lower
rep_confint_m1[1,2]/(sum(rep_confint_m1[1:2,2]) + var(heritability_predicted_m1[,3]))#upper

##females
#mb + bmr
rep_confint_f[1,1]/(sum(rep_confint_f[1:2,1]) + var(heritability_predicted_f[,2]))#lower
rep_confint_f[1,2]/(sum(rep_confint_f[1:2,2]) + var(heritability_predicted_f[,3]))#upper

#mb
rep_confint_f1[1,1]/(sum(rep_confint_f1[1:2,1])+ var(heritability_predicted_f1[,2]))#lower
rep_confint_f1[1,2]/(sum(rep_confint_f1[1:2,2]) + var(heritability_predicted_f1[,3]))#upper

