###################################################
#### MEASURE ALLELE FREQUENCY CHANGE OVER TIME ####
###################################################

# rm(list = ls())

#set working directory
setwd("~/Documents/Thesis_analysis/")

#load required libraries
library(ggplot2)
library(car)
#library(bestNormalize)
library(reshape2)
library(parallel)
library(foreach)
library(MASS)
library(snpStats)

#### DATA CLEANING ####
##read in the required data
#vcf file of SNPs filtered for MAF, HWE and linkage disequilibrium and missingness
#with hash removed from col headers
variants <- read.table("data/populations_out/0.05/plink_out/filtered_data/r2_pruned_filtered.txt", 
                       header = T, stringsAsFactors = F)

#demographic data (BY, sex)
Af_demography <- read.table("~/Documents/Data/A.flavicollis_demographic_data_duplicates_merged.txt", 
                            header = T, sep = "\t", stringsAsFactors = F)

##standardise all sample IDs for merged samples
#for demographic data, remove dashes
Af_demography$Sample <- sapply(strsplit(Af_demography$Sample, "-"), `[`, 1)

#for variant data, remove periods
colnames(variants)[10:ncol(variants)] <- sapply(strsplit(colnames(variants[,10:ncol(variants)]), "[.]"), `[`, 1)

#pull demographic data for samples that remain after filtering
Af_demography <- Af_demography[which(Af_demography$Sample %in% colnames(variants[,10:ncol(variants)]) == T),]

#make df of demographic data: ID Sex BirthYear and format for visped
LifeHist <- Af_demography[,c(1,3,8)]

#how many birth years are there?
unique(LifeHist$Birth_year)#4, 2015-2018

#remove NA birth years
LifeHist <- LifeHist[!is.na(LifeHist$Birth_year),]#366 samples remaining

#which mice were born when?
born_2015 <- LifeHist[LifeHist$Birth_year == 2015,]$Sample
born_2016 <- LifeHist[LifeHist$Birth_year == 2016,]$Sample
born_2017 <- LifeHist[LifeHist$Birth_year == 2017,]$Sample
born_2018 <- LifeHist[LifeHist$Birth_year == 2018,]$Sample

#how many in each?
length(born_2015)
length(born_2016)
length(born_2017)
length(born_2018)

#### CALCULATE ALLELE FREQUENCIES BY YEAR ####

#function calculate minor allele frequencies
MAF <- function(x) {
  minor_allele <- length(which(x == "1/0" | x == "0/1" | x == "1/1"))
  major_allele <- length(which(x == "1/0" | x == "0/1" | x == "0/0"))
  maf <- (minor_allele)/(minor_allele+major_allele)
  return(maf)
}

#calculate minor allele frequencies for each birth year and arrange into a single df
MAF_2015 <- data.frame(SNP_ID = variants$ID, BY = rep(2015, nrow(variants)), 
                       allele_frequency =  apply(variants[,which(colnames(variants) %in% born_2015)], 1, FUN = MAF)
)
MAF_2016 <- data.frame(SNP_ID = variants$ID, BY = rep(2016, nrow(variants)), 
                       allele_frequency =  apply(variants[,which(colnames(variants) %in% born_2016)], 1, FUN = MAF)
)
MAF_2017 <- data.frame(SNP_ID = variants$ID, BY = rep(2017, nrow(variants)), 
                       allele_frequency =  apply(variants[,which(colnames(variants) %in% born_2017)], 1, FUN = MAF)
)

allele_frequencies <- rbind(MAF_2015, MAF_2016, MAF_2017)
rm(list = ls(pattern = "MAF_"))
rm(list = ls(pattern = "born_"))

#how many alleles are missing in each year?
#most alleles are present every year except for 2018 where almost half the SNPs are missing.
#only 2 mice were caught in that year!
#change all 0 frequencies to NA
allele_frequencies$allele_frequency[allele_frequencies$allele_frequency == 0] <- NA

sum(is.na(allele_frequencies))
tapply(allele_frequencies$allele_frequency, allele_frequencies$BY, function (x) sum(is.na(x)))#21,189 missing total
allele_frequencies[which(is.na(allele_frequencies$allele_frequency)),]

#as only one SNP is missing in one year, remove it from the analysis
allele_frequencies <- na.omit(allele_frequencies)
allele_frequencies <- droplevels(allele_frequencies)

#### CAN WE PREDICT ALLELE FREQUENCIES FOR EACH SNP FOR A GIVEN YEAR? ####

#set up the dataframe in wide format to run lm by year for each SNP
freqs <- dcast(allele_frequencies[,1:3], BY ~ SNP_ID)
which(is.na(freqs[3,]))
freqs <- freqs[,-12612]#dcast still passes this NA but dont know why. remove it!
sum(is.na(freqs))

#in wide format in case you need it
allele_frequencies_wide <- data.frame(SNP_ID = colnames(freqs)[2:ncol(freqs)],
                                      by1 = as.numeric(freqs[1,2:ncol(freqs)]),
                                      by2 = as.numeric(freqs[2,2:ncol(freqs)]),
                                      by3 = as.numeric(freqs[3,2:ncol(freqs)]))

#this model predicts surprisingly well! will go with this one!
#run in parallel
cl <- makeCluster(4)#initialise a cluster with 4 cores (one for each row of df)
#untransformed data
clusterExport(cl, "freqs")#export df to each core in the cluster object
system.time(mod_out <- parApply(cl, freqs[,2:ncol(freqs)], 2, function(x) lm(x ~ BY, data = freqs)))#run model and time
stopCluster(cl)

### how well do observed and predicted values agree?
##extract model predictions for each SNP
mod_predictions <- list()

for (i in 1:length(mod_out)) {
  mod_predictions[[i]] <- predict(mod_out[[i]])
}

##add to a df
n.obs <- sapply(mod_predictions, length)#lengths f model predictions
seq.max <- seq_len(max(n.obs))#generate a sequence for max no. of predictions
mod_predictions <- data.frame(t(sapply(mod_predictions, "[", i = seq.max)))#combine predictions and fill with NAs
mod_predictions <- data.frame(SNP_ID = allele_frequencies_wide$SNP_ID, mod_predictions)
colnames(mod_predictions) <- colnames(allele_frequencies_wide)

##correlate model predictions with observed allele frequencies
#model predicts with 90% accuracy when untransformed!!
cor.test(c(allele_frequencies_wide$by1, allele_frequencies_wide$by2, allele_frequencies_wide$by3, allele_frequencies_wide$by4), 
         c(mod_predictions$by1, mod_predictions$by2, mod_predictions$by3, mod_predictions$by4))

plot(c(allele_frequencies_wide$by1, allele_frequencies_wide$by2, allele_frequencies_wide$by3, allele_frequencies_wide$by4) ~ 
       c(mod_predictions$by1, mod_predictions$by2, mod_predictions$by3, mod_predictions$by4))

abline(lm(c(allele_frequencies_wide$by1, allele_frequencies_wide$by2, allele_frequencies_wide$by3) ~ 
            c(mod_predictions$by1, mod_predictions$by2, mod_predictions$by3)))


#### CAN WE PREDICT ALLELE FREQUENCY CHANGE FOR EACH SNP #### 
#df to write allele frequency change to
allele_frequency_change_wide <- allele_frequencies_wide[,1:4]

#what is the allele frequency change per SNP?
allele_frequency_change_wide$by1 <- allele_frequencies_wide$by2 - allele_frequencies_wide$by1
allele_frequency_change_wide$by2 <- allele_frequencies_wide$by3 - allele_frequencies_wide$by2
#allele_frequency_change_wide$by3 <- allele_frequencies_wide$by4 - allele_frequencies_wide$by3
allele_frequency_change_wide$by3 <- allele_frequencies_wide$by3 - allele_frequencies_wide$by1

colnames(allele_frequency_change_wide)[2:4] <- c("delta2_1", "delta3_2", "delta3_1")

freqs1 <- matrix(nrow = 3, ncol = ncol(freqs))
freqs1[1,] <- c(1, as.numeric(allele_frequency_change_wide$delta2_1))
freqs1[2,] <- c(2, as.numeric(allele_frequency_change_wide$delta3_2))
freqs1[3,] <- c(3, as.numeric(allele_frequency_change_wide$delta3_1))

freqs1 <- data.frame(freqs1)
colnames(freqs1) <- c("delta_yr", colnames(freqs[2:ncol(freqs)]))

#year on year changes
cl <- makeCluster(4)
clusterExport(cl, "freqs1")
system.time(mod_out1 <- parApply(cl, freqs1[,2:ncol(freqs1)], 2, function(x) lm(x ~ delta_yr, data = freqs1)))#run model and time
stopCluster(cl)#close cluster object

### how well do observed and predicted values agree?
##extract model predictions for each SNP
mod_predictions1 <- list()

for (i in 1:length(mod_out1)) {
  mod_predictions1[[i]] <- predict(mod_out1[[i]])
}

##add to a df
n.obs <- sapply(mod_predictions1, length)#lengths f model predictions
seq.max <- seq_len(max(n.obs))#generate a sequence for max no. of predictions
mod_predictions1 <- data.frame(t(sapply(mod_predictions1, "[", i = seq.max)))#combine predictions and fill with NAs
mod_predictions1 <- data.frame(SNP_ID = allele_frequency_change_wide$SNP_ID, mod_predictions1)
colnames(mod_predictions1) <- colnames(allele_frequency_change_wide[1:3])

##correlate model predictions with observed allele frequencies
#model predicts with 88%% accuracy when untransformed!! predicts allele frequency change well!
cor.test(c(allele_frequency_change_wide$delta2_1, allele_frequency_change_wide$delta3_2), 
         c(mod_predictions1$delta2_1, mod_predictions1$delta3_2))

plot(c(allele_frequency_change_wide$delta2_1, allele_frequency_change_wide$delta3_2, allele_frequency_change_wide$delta4_3) ~ 
       c(mod_predictions1$delta2_1, mod_predictions1$delta3_2, mod_predictions1$delta4_3))

abline(lm(c(allele_frequency_change_wide$delta2_1, allele_frequency_change_wide$delta3_2, allele_frequency_change_wide$delta4_3) ~ 
            c(mod_predictions1$delta2_1, mod_predictions1$delta3_2, mod_predictions1$delta4_3)))

#how many snps increase and how many snps decrease in the population overall?
length(allele_frequency_change_wide$delta3_1[allele_frequency_change_wide$delta3_1 == 0]) #no net change
length(allele_frequency_change_wide$delta3_1[allele_frequency_change_wide$delta3_1 < 0]) #decrease
length(allele_frequency_change_wide$delta3_1[allele_frequency_change_wide$delta3_1 > 0]) #increase
#### PLOTTING ####

#function for estimating point density for plotting
# Get density of points in 2 dimensions.
# @param x A numeric vector.
# @param y A numeric vector.
# @param n Create a square n by n grid to compute density.
# @return The density within each square.
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

## PREDICTED AND OBSERVED ALLELE FREQUENCY CORRELATIONS
#make a df
cor.test(c(allele_frequencies_wide$by1, allele_frequencies_wide$by2, 
           allele_frequencies_wide$by3, allele_frequencies_wide$by4), 
         c(mod_predictions$by1, mod_predictions$by2, 
           mod_predictions$by3, mod_predictions$by4))

AF_correlations <- na.omit(data.frame(obs_AF = c(allele_frequencies_wide$by1, allele_frequencies_wide$by2, 
                                                 allele_frequencies_wide$by3, allele_frequencies_wide$by4),
                                      pred_AF = c(mod_predictions$by1, mod_predictions$by2, 
                                                  mod_predictions$by3, mod_predictions$by4)))

AF_correlations$density <- get_density(AF_correlations$obs_AF, AF_correlations$pred_AF, n = 1000)

#png("results/allele_frequency_change/obs_predAF.png", units = "in", width = 7.15, height = 6.38, res = 200)
ggplot(AF_correlations, aes(x = obs_AF, y = pred_AF, colour = density)) +
  geom_point(alpha = 0.8) +
  stat_smooth(method = "lm", se = F, colour = "darkred") +
  # geom_abline(col = "darkred") +
  xlab("Observed allele frequency") +
  ylab("Predicted allele frequency") +
  scale_color_continuous(low = "midnightblue", high = "dodgerblue4") +
  theme_bw() +
  theme(legend.position = "none")
#dev.off()

##plot for predicted and observed correlations for allele frequency change per year
#make a df
AFchange_correlations <- na.omit(data.frame(obs_AF = c(allele_frequency_change_wide$delta2_1,
                                                       allele_frequency_change_wide$delta3_2,
                                                       allele_frequency_change_wide$delta4_3),
                                            pred_AF = c(mod_predictions1$delta2_1, 
                                                        mod_predictions1$delta3_2, 
                                                        mod_predictions1$delta4_3)))
AFchange_correlations$density <- get_density(AFchange_correlations$obs_AF, 
                                             AFchange_correlations$pred_AF, n = 1000)

# png("results/allele_frequency_change/obs_predAFchange.png", units = "in", width = 7.15, height = 6.38, res = 200)
ggplot(AFchange_correlations, aes(x = obs_AF, y = pred_AF, colour = density)) +
  geom_point(alpha = 0.8) +
  stat_smooth(method = "lm", se = F, colour = "darkred") +
  #geom_abline(col = "darkred") +
  xlab("Observed allele frequency change") +
  ylab("Predicted allele frequency change") +
  scale_color_continuous(low = "midnightblue", high = "dodgerblue4") +
  theme_bw() +
  theme(legend.position = "none")
# dev.off()

# EXAMPLE OF LARGEST INCREASE IN ALLELE FREQUENCIES 2015-2017

# Function to calculate the mean and the standard deviation
# for each group
# @data : a data frame
# @varname : the name of a column containing the variable
#to be summariezed
# @groupnames : vector of column names to be used as
# grouping variables
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

# function to calculate 95% CI of null allele frequencies
CI <- function(x){
  upper <- sort(x)[0.975*length(x)]
  lower <- sort(x)[0.025*length(x)]
  return(data.frame(upper = upper, lower = lower))
}

#greatest increase allele
hist(BY3_1_null_frequency_change[12389,])
range(BY3_1_null_frequency_change[1,])
allele_frequency_change_wide$delta3_1[12389]

#convert by3-by1 matrix into a df for ggplot
BY3_1_null_frequency_change_gg <- as.data.frame(BY3_1_null_frequency_change)
BY3_1_null_frequency_change_gg$SNP_ID <- allele_frequencies_wide$SNP_ID
BY3_1_null_frequency_change_gg <- melt(BY3_1_null_frequency_change_gg[BY3_1_null_frequency_change_gg$SNP_ID == "un_26255",], 
                                       id.vars = c("SNP_ID"))
colnames(BY3_1_null_frequency_change_gg)[2:3] <- c("replicate", "AF_change")
un_26255 <- allele_frequency_change_wide[allele_frequency_change_wide$SNP_ID == "un_26255",]$delta3_1

# pdf("results/allele_frequency_change/un_26266_AFchange.pdf", width = 7.15, height = 6.38)
ggplot(BY3_1_null_frequency_change_gg, aes(x = AF_change)) + 
  geom_histogram(alpha = 0.8, size = 0) +
  geom_vline(aes(xintercept = mean(AF_change)), 
             colour = "blue", linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = un_26255), colour = "darkred", 
             linetype = "dashed", size = 1) + 
  annotate(geom = "text", x = 0.14, y = 400, label = "Observed change", 
           colour = "darkred") +
  annotate(geom = "text", x = 0.24, y = 400, label = "Expected change", 
           colour = "blue") +
  xlab("Allele frequency change (SNP un_26255)") +
  ylab("Count") +
  theme_bw()
# dev.off()

un_26255_observed_freq <- as.numeric(allele_frequencies_wide[allele_frequencies_wide$SNP_ID == "un_26255",][2:4])

un_26255_trace <- data.frame(SNP_ID = rep("un_26255", 15003), 
                             study_year = c(rep("2015", 5000),
                                            rep("2016", 5000),
                                            rep("2017", 5000),
                                            # rep("2018", 5000),
                                            as.character(2015:2017)),
                             type = c(rep("Expected", 15000), rep("Observed", 3)),
                             allele_frequency = c(BY1_null_frequencies[12389,],
                                                  BY2_null_frequencies[12389,],
                                                  BY3_null_frequencies[12389,],
                                                  # BY4_null_frequencies[12389,],
                                                  un_26255_observed_freq))

#calculate mean null allele frequencies, sd and CI for each year
un_26255_summary <- data_summary(un_26255_trace, varname = "allele_frequency",
                                 groupnames = c("study_year", "type"))

#calculate where 95% of null points lie
confidence <- tapply(un_26255_trace[1:15000,]$allele_frequency, un_26255_trace[1:15000,]$study_year, CI)
un_26255_summary$upper <- c(confidence$`2015`$upper, NA, confidence$`2016`$upper, 
                            NA, confidence$`2017`$upper, NA)
un_26255_summary$lower <- c(confidence$`2015`$lower, NA, confidence$`2016`$lower, 
                            NA, confidence$`2017`$lower, NA)
un_26255_summary

# pdf("results/allele_frequency_change/un_26266_AFchange_perYear.pdf", width = 7.15, height = 6.38)
ggplot(un_26255_summary, aes(x = study_year, y = allele_frequency, group = type, colour = type)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4, colour = NA) +
  geom_line() +
  scale_colour_manual(values = c("blue", "darkred")) +
  xlab(NULL) +
  ylab("Allele frequency (SNP un_26255)") +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = c(0.8, 0.25),
        legend.text = element_text(size = 11))
# dev.off()

# EXAMPLE OF LARGEST DECREASE IN ALLELE FREQUENCIES 2015-2017

#greatest decrease allele
hist(BY3_1_null_frequency_change[7572,])
range(BY3_1_null_frequency_change[1,])
allele_frequency_change_wide$delta3_1[7572]

#convert by3-by1 matrix into a df for ggplot
BY3_1_null_frequency_change_gg <- as.data.frame(BY3_1_null_frequency_change)
BY3_1_null_frequency_change_gg$SNP_ID <- allele_frequencies_wide$SNP_ID
BY3_1_null_frequency_change_gg <- melt(BY3_1_null_frequency_change_gg[BY3_1_null_frequency_change_gg$SNP_ID == "un_19825",], 
                                       id.vars = c("SNP_ID"))
colnames(BY3_1_null_frequency_change_gg)[2:3] <- c("replicate", "AF_change")
un_19825 <- allele_frequency_change_wide[allele_frequency_change_wide$SNP_ID == "un_19825",]$delta3_1

# pdf("results/allele_frequency_change/un_19825_AFchange.pdf", width = 7.15, height = 6.38)
ggplot(BY3_1_null_frequency_change_gg, aes(x = AF_change)) + 
  geom_histogram(alpha = 0.8, size = 0) +
  geom_vline(aes(xintercept = mean(AF_change)), 
             colour = "blue", linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = un_19825), colour = "darkred", 
             linetype = "dashed", size = 1) + 
  annotate(geom = "text", x = -0.24, y = 450, label = "Observed change", 
           colour = "darkred") +
  annotate(geom = "text", x = -0.325, y = 450, label = "Expected change", 
           colour = "blue") +
  xlab("Allele frequency change (SNP un_19825)") +
  ylab("Count") +
  theme_bw()
# dev.off()

un_19825_observed_freq <- as.numeric(allele_frequencies_wide[allele_frequencies_wide$SNP_ID == "un_19825",][2:4])

un_19825_trace <- data.frame(SNP_ID = rep("un_19825", 15003), 
                             study_year = c(rep("2015", 5000),
                                            rep("2016", 5000),
                                            rep("2017", 5000),
                                            # rep("2018", 5000),
                                            as.character(2015:2017)),
                             type = c(rep("Expected", 15000), rep("Observed", 3)),
                             allele_frequency = c(BY1_null_frequencies[7572,],
                                                  BY2_null_frequencies[7572,],
                                                  BY3_null_frequencies[7572,],
                                                  # BY4_null_frequencies[12389,],
                                                  un_19825_observed_freq))

un_19825_summary <- data_summary(un_19825_trace, varname = "allele_frequency",
                                 groupnames = c("study_year", "type"))

confidence <- tapply(un_19825_trace[1:15000,]$allele_frequency, un_19825_trace[1:15000,]$study_year, CI)
un_19825_summary$upper <- c(confidence$`2015`$upper, NA, confidence$`2016`$upper, 
                            NA, confidence$`2017`$upper, NA)
un_19825_summary$lower <- c(confidence$`2015`$lower, NA, confidence$`2016`$lower, 
                            NA, confidence$`2017`$lower, NA)
un_19825_summary

# pdf("results/allele_frequency_change/un_19825_AFchange_perYear.pdf", width = 7.15, height = 6.38)
ggplot(un_19825_summary, aes(x = study_year, y = allele_frequency, group = type, colour = type)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4, colour = NA) +
  geom_line() +
  scale_colour_manual(values = c("blue", "darkred")) +
  xlab(NULL) +
  ylab("Allele frequency (SNP un_19825)") +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = c(0.8, 0.8),
        legend.text = element_text(size = 11))
# dev.off()

