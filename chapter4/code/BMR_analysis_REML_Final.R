#### FULL ANALYSIS OF BMRmin ####

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

#### data organisation ####
#read in the data
phenotypic_data_full <- read.table("data/phenotypic_data/A.flavicollis_phenotypic_data_clean.txt", header = T, stringsAsFactors = F, sep = "\t")
phenotypic_data_full$Trapping_year <- as.factor(phenotypic_data_full$Trapping_year)

###create long format table for BMR data
#transform phenotypic data into long format by season
BMR_data <- melt(phenotypic_data_full, id.vars = c("ID", "merged_samples", "Sex", "Age","Birth_year", "Trapping_year"), 
                 measure.vars = c("BMR_autumn", "minBMR_winter.spring"))
BMR_data$variable <- as.character(BMR_data$variable)
BMR_data$variable[BMR_data$variable == "BMR_autumn"] <- "Autumn"
BMR_data$variable[BMR_data$variable == "minBMR_winter.spring"] <- "Winter/Spring"
BMR_data$Sex[BMR_data$Sex == "m"] <- "Male"
BMR_data$Sex[BMR_data$Sex == "f"] <- "Female"
colnames(BMR_data)[c(7:8)] <- c("Season", "BMR_min")

#add the masses associated with bmr measurements
BMR_mass <- melt(phenotypic_data_full, id.vars = c("merged_samples"), measure.vars = c("mb_BMR_autumn", "mb_BMR_winter.spring_minBMR"))
colnames(BMR_mass)[2:3] <- c("Season", "mb")
BMR_mass$Season <- as.character(BMR_mass$Season)
BMR_mass$Season[BMR_mass$Season == unique(BMR_mass$Season)[1]] <- "Autumn"
BMR_mass$Season[BMR_mass$Season == unique(BMR_mass$Season)[2]] <- "Winter/Spring"
BMR_data$mb_bmr <- BMR_mass$mb

#combine both df and remomve duplicates
BMR_data <- BMR_data[!is.na(BMR_data$BMR_min),]
BMR_data <- BMR_data[!duplicated(BMR_data$BMR_min),]
rm(phenotypic_data_full, BMR_mass)

#turn questionable ages into NA
BMR_data$Age[BMR_data$Age == "ad?"] <- NA
BMR_data$Age[BMR_data$Age == "juv?"] <- NA

##is BMR normal?
#seasons together
hist(BMR_data$BMR_min)
shapiro.test(BMR_data$BMR_min)#definitely not normal

##normalise data for when needed as a response variable
trans <- orderNorm(BMR_data$BMR_min)
shapiro.test(trans$x.t)
BMR_data$BMR_min_gaus <- trans$x.t
hist(BMR_data$BMR_min_gaus)

#calculate the mean for BMR_min
mean_BMRmin <- melt(tapply(BMR_data$BMR_min, BMR_data$ID, mean, na.rm = T), value.name = "mean_BMR_min")
colnames(mean_BMRmin)[1] <- "ID"
BMR_data <- merge(BMR_data, mean_BMRmin, by = "ID")
rm(mean_BMRmin)

#function to calculate standard error
std <- function(x){
  sd(x, na.rm = T)/sqrt(length(x))
}

#calculate std err
std_errBMR_min <- data.frame(ID = names(tapply(BMR_data$BMR_min, BMR_data$ID, std)),
                             std_errBMR_min = tapply(BMR_data$BMR_min, BMR_data$ID, std))


#merge phenotypic data and standard error data by ID
BMR_data <- Reduce(function(x,y) merge(x = x, y = y, by = "ID"), list(BMR_data, std_errBMR_min))
rm(std_errBMR_min)

#### DATA SUMMARY AND ANALYSIS ####
#data summary
summary(BMR_data$BMR_min)#overall
summary(BMR_data[BMR_data$Season == "Autumn",]$BMR_min)#autumn
summary(BMR_data[BMR_data$Season == "Winter/Spring",]$BMR_min)#winer/spring

sd(BMR_data$BMR_min)#overall
sd(BMR_data[BMR_data$Season == "Autumn",]$BMR_min, na.rm = T)
sd(BMR_data[BMR_data$Season == "Winter/Spring",]$BMR_min, na.rm = T)

#how many mice measured in...
length(unique(BMR_data[BMR_data$Season == "Autumn" & !is.na(BMR_data$BMR_min),]$ID))#autumn?
length(unique(BMR_data[BMR_data$Season == "Winter/Spring" & !is.na(BMR_data$BMR_min),]$ID))#winter/spring?

#how many mice phenotyped by...
table(BMR_data[which(!duplicated(BMR_data$ID) == T),]$Sex)#sex?
table(BMR_data[which(!duplicated(BMR_data$ID) == T),]$Age)#age?
table(BMR_data[which(!duplicated(BMR_data$ID) == T),]$Trapping_year)#trapping year?

##correct BMR_min by mb
#model summary
summary(lm(BMR_min_gaus ~ mb_bmr, data = BMR_data))
cor.test(BMR_data$BMR_min_gaus, BMR_data$mb_bmr)

#data frame of predicted values
mass_corrected_bmr <- data.frame(rows = names(predict(lm(BMR_min_gaus ~ mb_bmr, data = BMR_data))), 
                                 mass_corrected_BMR_min_gaus = predict(lm(BMR_min_gaus ~ mb_bmr, data = BMR_data)))

#merge mass corrected bmr with bmr data
BMR_data$rows <- row.names(BMR_data)
BMR_data <- merge(mass_corrected_bmr, BMR_data, by = "rows", all = T)
BMR_data <- BMR_data[,-1]
rm(mass_corrected_bmr)

##Is there a difference between...
#Sex? NOT MASS CORRECTED
t.test(BMR_data$BMR_min_gaus ~ BMR_data$Sex)#YES overall. pattern may be different by season
sd(BMR_data[BMR_data$Sex == "Female",]$BMR_min_gaus, na.rm = T)
sd(BMR_data[BMR_data$Sex == "Male",]$BMR_min_gaus, na.rm = T)

#SeX? MASS CORRCTED
#diff between sexes?
t.test(BMR_data$mass_corrected_BMR_min_gaus ~ BMR_data$Sex)#YES overall. pattern may be different by season
sd(BMR_data[BMR_data$Sex == "Female",]$mass_corrected_BMR_min_gaus, na.rm = T)
sd(BMR_data[BMR_data$Sex == "Male",]$mass_corrected_BMR_min_gaus, na.rm = T)

#Age? NOT MASS CORRECTED
t.test(BMR_data$BMR_min_gaus ~ BMR_data$Age)#NO

#Age? MASS CORRECTED
t.test(BMR_data$mass_corrected_BMR_min_gaus ~ BMR_data$Age)#NO

#mass corrected bmr vs trapping year and season and sex 
summary(lm(BMR_data$mass_corrected_BMR_min_gaus ~ BMR_data$Trapping_year + BMR_data$Season))#NO for year but YES for season
plot(lm(BMR_data$mass_corrected_BMR_min_gaus ~ BMR_data$Trapping_year + BMR_data$Season))

#does mb differ between seasons?
t.test(BMR_data$mb_bmr ~ BMR_data$Season)#YES... mice measured in autumn are lighter/smaller
tapply(BMR_data$mb_bmr, BMR_data$Season, sd)

#### SUMMARY PLOTTING ####
##plot the distribution of BMR 
#plots by season and year
ggplot(BMR_data, aes(x = Trapping_year, y = mass_corrected_BMR_min_gaus, fill = Season)) +
  geom_boxplot(position = position_dodge(preserve = "single")) +
  ylab(expression(paste("Order-Norm ", BMR[min], " (mW)"))) + #(Order-NormBMR[min] (mW))) +
  xlab("Trapping year") +
  scale_fill_manual(values = c("lightblue", "#FFCC99")) +
  #facet_grid(. ~ Sex) +
  theme_bw() +
  theme(legend.position = c(0.8, 0.15))

### plot individual variation in BMRmin
##order the data by mean_BMR_min
BMR_data <- BMR_data[order(BMR_data$mean_BMR_min),]

#function to greate ordereg group IDs
grpid <- function(x) {
  match(x, unique(x))
}

#create unique group IDs
BMR_data <- BMR_data %>% mutate(ordered_sample = group_indices(., ID) %>% grpid())

ggplot(data = BMR_data, aes(x = ordered_sample , y = mean_BMR_min, colour = Sex)) +
  geom_point() +
  geom_pointrange(data = BMR_data, aes(ymin = mean_BMR_min - std_errBMR_min,
                                       ymax = mean_BMR_min + std_errBMR_min)) +
  scale_colour_manual(values = c("#6495ED85", "#C1403D85")) +
  ylab(expression(paste(BMR[min], " (mW)"))) +
  xlab("Sample") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        legend.position = c(0.2, 0.8))

##plot BMR vs mb
#plot
ggplot(BMR_data, aes(x = mb_bmr, y = BMR_min_gaus, colour = Sex)) +
  geom_point() +#colour = "midnightblue") +
  scale_colour_manual(values = c("#6495ED85", "#C1403D85")) +
  ylab(expression(paste("Order-Norm ", BMR[min], " (mW)"))) +
  xlab(expression(paste(m[b], "(g)"))) +
  stat_smooth(method = "lm", colour = "#660000", se = F) +
  #annotate(geom = "text", x = 40, y = 2.5, label = expression(paste("y = -1.74 + 0.046x, ", "r"^2, " = 0.102")), fontface = "italic") +
  theme_bw() +
  theme(legend.position = c(0.8, 0.2))

#### CALCULATE REPEATABILITY ####

## REPEATABILITY OF BMR CAN ONLY BE ESTIMATED IN MICE MEASURED IN BOTH SEASONS AS THERE IS NO WITHIN 
## INDIVIDUAL VARIATION FOR MICE ONLY MEASURED IN WINTER/SPRING (ONLY BMRmin MEASURED). 

## THERE IS NOT ENOUGH WITHIN INDIVIDUAL VARIANCE (ONLY 2 MEASURES PER INDIVIDUAL - ONE FOR EACH SEASON)
## AND BMR BETWEEN SEASONS VARIES TOO MUCH SO BMR IS NOT REPEATABLE!!! CANNOT CALCULATE HERITABILITY!!!

#model for repeatability including sex, mb, trapping year, season
#age not included as initial analysis showed no difference anyway!
repeatability_mod <- lmer(BMR_min_gaus ~ Sex + mb_bmr +#sex insignificant
                            (1 | ID) + #causing singular fit as not enough data
                            (1 | Trapping_year) + 
                            (1 | Season), 
                          data = BMR_data)

summary(repeatability_mod)
Anova(repeatability_mod)






