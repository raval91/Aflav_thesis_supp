#### NULL ALLELE FREQUENCY SIMULATIONS ####

#rm(list=ls())

#set the wd
setwd("~/Documents/Thesis_analysis/")

#load the required libraries
library(parallel)
library(doParallel)
library(foreach)
library(ggplot2)

# READ IN THE REQUIRED DATA
# MAF for each SNP in each study/birth year
allele_frequencies <- read.table("data/allele_frequency_change/allele_frequencies.txt", 
                                      header = T, sep = "\t", stringsAsFactors = F)

# MAF change for each year and SNP
allele_frequency_change <- read.table("data/allele_frequency_change/allele_frequency_change.txt", 
                                      header = T, sep = "\t", stringsAsFactors = F)

#### GENERATE NULL ALLELE FREQUENCY DISTRIBUTIONS FOR COMPARISON  WITH OBSERVED ALLELE FREQUENCIES ####

# function to calculate MAF
MAF <- function(x) {
  minor_allele <- length(which(x == "1/0" | x == "0/1" | x == "1/1"))
  major_allele <- length(which(x == "1/0" | x == "0/1" | x == "0/0"))
  maf <- (minor_allele)/(minor_allele+major_allele)
  return(maf)
}

# function to calculate 95% CI of null allele frequencies
CI <- function(x){
  upper <- sort(x, method = "quick")[0.975*length(x)]
  lower <- sort(x, method = "quick")[0.025*length(x)]
  return(list(upper = upper, lower = lower))
}

# function to deregister any foreach clusters
unregister <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}

#  function to simulate null MAF data sampling 500 genotypes per SNP
# @samplings: the number of times to draw alleles to form a distribution of null frequencies
simulations <- function(samplings) {  
  
  #initialise vectors and matrices to write null allele frequencies to
  null_frequencies <- vector(mode = "list", length = 5)#list to add matrices of null MAF to
  x <- matrix(ncol = ncol(allele_frequencies)-1, nrow = nrow(allele_frequencies))#matrix to write MAF  for every replicate to
  null_matrices <- list(matrix(nrow = nrow(allele_frequencies), ncol = 593),
                        matrix(nrow = nrow(allele_frequencies), ncol = 593),
                        matrix(nrow = nrow(allele_frequencies), ncol = 593))
  
  #for n replicates to generate a distribution
  for(i in 1:samplings) {
    print(paste("replicate", i))
    
    #for each study year
    for (k in 1:length(null_matrices)){
      print(paste("Study year ", k))
      
      if (k == 1){#if study year is 1
        #for each SNP
        for (l in 1:nrow(null_matrices[[k]])) {
          
          #sample alleles with probability = frequency of observed MAF in year 1
          null_matrices[[k]][l,] <- sample(alleles, ncol(null_matrices[[k]]), replace = T,
                                           prob = c((1-allele_frequencies$by1[l])^2,
                                                    2*(allele_frequencies$by1[l] * (1-allele_frequencies$by1[l])),
                                                    (allele_frequencies$by1[l])^2))
        }
      } else if (k == 2) {# else if study year is 2
        #for each SNP
        for (l in 1:nrow(null_matrices[[k]])) {
          
          #sample alleles with probability = frequency of observed MAF  in year 2
          null_matrices[[k]][l,] <- sample(alleles, ncol(null_matrices[[k]]), replace = T,
                                           prob = c((1-allele_frequencies$by2[l])^2,
                                                    2*(allele_frequencies$by2[l] * (1-allele_frequencies$by2[l])),
                                                    (allele_frequencies$by2[l])^2))
        }
      } else if (k == 3) {# else if study year is 3
        #for each SNP
        for (l in 1:nrow(null_matrices[[k]])) {
          
          #sample alleles with probability = frequency of observed MAF  in year 3
          null_matrices[[k]][l,] <- sample(alleles, ncol(null_matrices[[k]]), replace = T,
                                           prob = c((1-allele_frequencies$by3[l]^2),
                                                    2*(allele_frequencies$by3[l] * (1-allele_frequencies$by3[l])),
                                                    (allele_frequencies$by3[l])^2))
        }
      } else {# else sample from study year 4
        #for each SNP
        for (l in 1:nrow(null_matrices[[k]])) {
          
          #sample alleles with probability = frequency of observed MAF  in year 4
          null_matrices[[k]][l,] <- sample(alleles, ncol(null_matrices[[k]]), replace = T,
                                           prob = c((1-allele_frequencies$by4[l])^2,
                                                    2*(allele_frequencies$by4[l] * (1-allele_frequencies$by4[l])),
                                                    (allele_frequencies$by4[l])^2))
        }
      }
      
      # calculate MAF for each study year and write to matrix
      x[,k] <- apply(null_matrices[[k]], 1, MAF)
    }
    #save matrix
    null_frequencies[[i]] <- x
  }
  return(null_frequencies)
}


# SIMULATE NULL MAF DISTRIBUTIONS IN PARALLEL

# alleles to sample from
alleles <- c("0/0","0/1","1/1")#homozygous ref, heterozygous, homozygous alt

#run simulation
cl <- parallel::makeForkCluster(detectCores())
doParallel::registerDoParallel(cl)
system.time(
  null_MAF_sims <- foreach(i=1:8) %dopar% {
    set.seed(i)
    simulations(samplings = 5000)
  }
)
# parallel::stopCluster(cl)
# registerDoSEQ()
unregister()

#### CALCULATE 95% RANGE OF AF CHANGE FOR EVERY SNP (ROW) IN EVERY MATRIX  FOR EACH SIMULATION ####
#empty matrices to write simulation output to for each study/birth year
BY1_2_simOut <- vector(mode = "list", length = length(null_MAF_sims))
BY2_3_simOut <- vector(mode = "list", length = length(null_MAF_sims))
BY3_1_simOut <- vector(mode = "list", length = length(null_MAF_sims))

#empty matrices to write 95% CI for AF change to
CI_1_2 <- vector(mode = "list", length = length(null_MAF_sims))
CI_2_3 <- vector(mode = "list", length = length(null_MAF_sims))
CI_3_1 <- vector(mode = "list", length = length(null_MAF_sims))

# empty matrices to write matrix output for each year
BY1_null_frequencies <- matrix(nrow = nrow(null_MAF_sims[[1]][[1]]), ncol = length(null_MAF_sims[[1]]))
BY2_null_frequencies <- matrix(nrow = nrow(null_MAF_sims[[1]][[1]]), ncol = length(null_MAF_sims[[1]]))
BY3_null_frequencies <- matrix(nrow = nrow(null_MAF_sims[[1]][[1]]), ncol = length(null_MAF_sims[[1]]))

system.time(
  #for every simulation  
  for(i in 1:length(null_MAF_sims)){
    #for every replicate
    for(j in 1:length(null_MAF_sims[[i]])){
      
      #extract null frequencies for each year and write their own matrix
      BY1_null_frequencies[,j] <- null_MAF_sims[[i]][[j]][,1]
      BY2_null_frequencies[,j] <- null_MAF_sims[[i]][[j]][,2]
      BY3_null_frequencies[,j] <- null_MAF_sims[[i]][[j]][,3]
    }
    
    #calculate allele frequency change
    BY1_2_simOut[[i]] <- BY2_null_frequencies - BY1_null_frequencies
    BY2_3_simOut[[i]] <- BY3_null_frequencies - BY2_null_frequencies
    BY3_1_simOut[[i]] <- BY3_null_frequencies - BY1_null_frequencies
    
    #calculate 95% confidence limits
    CI_1_2[[i]] <- apply(BY1_2_simOut[[i]], 1, CI)
    CI_2_3[[i]] <- apply(BY2_3_simOut[[i]], 1, CI)
    CI_3_1[[i]] <- apply(BY3_1_simOut[[i]], 1, CI)
  }
)
# rm(list = ls(pattern = "simOut"))
# rm(list = ls(pattern = "null_frequencies"))

#### CALCULATE P-VALUES FOR ALLELE FREQUENCY CHANGE ####

# make empty matrices to write pvalues for each SNP from each simulation to
sig_counts_BY1_2 <- matrix(nrow = length(CI_1_2[[i]]), ncol = length(CI_1_2))
sig_counts_BY2_3 <- matrix(nrow = length(CI_1_2[[i]]), ncol = length(CI_1_2))
sig_counts_BY3_1 <- matrix(nrow = length(CI_1_2[[i]]), ncol = length(CI_1_2))

system.time(
  
  # for each simulation
  for (i in 1:length(CI_1_2)){
  
    # for each SNP
    for (j in 1:length(CI_1_2[[i]])){
    
      # is the change in MAF for each SNP outside the 95% confidence range of null distributions?
      sig_counts_BY1_2[j,i] <- allele_frequency_change[j,2] >= CI_1_2[[i]][[j]]$upper | 
        allele_frequency_change[j,2] <= CI_1_2[[i]][[j]]$lower
      
      sig_counts_BY2_3[j,i] <- allele_frequency_change[j,3] >= CI_2_3[[i]][[j]]$upper | 
        allele_frequency_change[j,3] <= CI_2_3[[i]][[j]]$lower
      
      sig_counts_BY3_1[j,i] <- allele_frequency_change[j,4] >= CI_3_1[[i]][[j]]$upper | 
        allele_frequency_change[j,4] <= CI_3_1[[i]][[j]]$lower
    }
    # calculate the p vals
    pvals_1_2 <- 1-(rowSums(sig_counts_BY1_2)/length(null_MAF_sims))
    pvals_2_3 <- 1-(rowSums(sig_counts_BY2_3)/length(null_MAF_sims))
    pvals_3_1 <- 1-(rowSums(sig_counts_BY3_1)/length(null_MAF_sims))
  }
)
rm(list = ls(pattern = "sig_counts"))

p_values <- data.frame(SNP_ID = allele_frequencies$SNP_ID, 
                       pvals_1_2 = pvals_1_2, 
                       pvals_2_3 = pvals_2_3, 
                       pvals_3_1 = pvals_3_1)

hist(p_values$pvals_3_1)
hist(p_values$pvals_1_2)
hist(p_values$pvals_2_3)
# write.table(p_values, "results/allele_frequency_change/test/p_values.txt", col.names = T, 
#             row.names = F, quote = F, sep = "\t")  

#### PLOTTING ####

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

# DF of null change from 1 simulation
null_change_3_1 <- apply(BY3_1_simOut[[1]], 1, median)
null_change_3_1 <- data.frame(SNP = allele_frequencies$SNP_ID, AF_change = null_change_3_1)

# which allele showed the greatest change in frequency from the sims throughout the study?
hist(BY3_1_simOut[[1]][which.max(null_change_3_1$AF_change),])# MAX, SNP un_16920
hist(BY3_1_simOut[[1]][which.min(null_change_3_1$AF_change),])# MIN, SNP un_19825

un_16920 <- data.frame(un_16920 = BY3_1_simOut[[1]][which.max(null_change_3_1$AF_change),])
un_19825 <- data.frame(un_19825 = BY3_1_simOut[[1]][which.min(null_change_3_1$AF_change),])

# greatest increase in allele frequency overall
ggplot(un_16920, aes(x = un_16920)) + 
  geom_histogram(alpha = 0.8, size = 0) +
  geom_vline(aes(xintercept = median(un_16920)), 
             colour = "blue", linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = allele_frequency_change[which(allele_frequency_change$SNP_ID == "un_16920"),]$delta3_1), 
             colour = "darkred", 
             linetype = "dashed", size = 1) + 
  annotate(geom = "text", x = 0.18, y = 410, label = "Observed change", 
           colour = "darkred") +
  annotate(geom = "text", x = 0.12, y = 410, label = "Expected change", 
           colour = "blue") +
  xlab("Allele frequency change (SNP un_16920)") +
  ylab("Count") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14),
        axis.text = element_text(size = 12))


# greatest decrease in allele frequency overall
ggplot(un_19825, aes(x = un_19825)) + 
  geom_histogram(alpha = 0.8, size = 0) +
  geom_vline(aes(xintercept = median(un_19825)), 
             colour = "blue", linetype = "dashed", size = 1) +
  geom_vline(aes(xintercept = allele_frequency_change[which(allele_frequency_change$SNP_ID == "un_19825"),]$delta3_1), 
             colour = "darkred", 
             linetype = "dashed", size = 1) + 
  annotate(geom = "text", x = -0.24, y = 450, label = "Observed change", 
           colour = "darkred") +
  annotate(geom = "text", x = -0.3, y = 450, label = "Expected change", 
           colour = "blue") +
  xlab("Allele frequency change (SNP un_19825)") +
  ylab("Count") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14),
        axis.text = element_text(size = 12))

## plots by year
# greatest increase
un_16920_observed <- as.numeric(allele_frequencies[which(allele_frequency_change$SNP_ID == "un_16920"),2:4])

un_16920_trace <- data.frame(SNP_ID = rep("un_16920", 15003), 
                             study_year = c(rep("2015", 5000),
                                            rep("2016", 5000),
                                            rep("2017", 5000),
                                            as.character(2015:2017)),
                             type = c(rep("Expected", 15000), rep("Observed", 3)),
                             allele_frequency = c(BY1_null_frequencies[which(allele_frequencies$SNP_ID == "un_16920"),],
                                                  BY2_null_frequencies[which(allele_frequencies$SNP_ID == "un_16920"),],
                                                  BY3_null_frequencies[which(allele_frequencies$SNP_ID == "un_16920"),],
                                                  un_16920_observed))

#calculate mean null allele frequencies, sd and CI for each year
un_16920_summary <- data_summary(un_16920_trace, varname = "allele_frequency",
                                 groupnames = c("study_year", "type"))

#calculate where 95% of null points lie
confidence <- tapply(un_16920_trace[1:15000,]$allele_frequency, un_16920_trace[1:15000,]$study_year, CI)
un_16920_summary$upper <- c(confidence$`2015`$upper, NA, confidence$`2016`$upper, 
                            NA, confidence$`2017`$upper, NA)
un_16920_summary$lower <- c(confidence$`2015`$lower, NA, confidence$`2016`$lower, 
                            NA, confidence$`2017`$lower, NA)
un_16920_summary

ggplot(un_16920_summary, aes(x = study_year, y = allele_frequency, group = type, colour = type)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4, colour = NA) +
  geom_line() +
  scale_colour_manual(values = c("blue", "darkred")) +
  xlab(NULL) +
  ylab("Allele frequency (SNP un_16920)") +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = c(0.2, 0.8),
        legend.text = element_text(size = 12),
        axis.title= element_text(size = 14),
        axis.text = element_text(size = 12))

# greatest decrease
un_19825_observed <- as.numeric(allele_frequencies[which(allele_frequencies$SNP_ID == "un_19825"),2:4])

un_19825_trace <- data.frame(SNP_ID = rep("un_19825", 15003), 
                             study_year = c(rep("2015", 5000),
                                            rep("2016", 5000),
                                            rep("2017", 5000),
                                            as.character(2015:2017)),
                             type = c(rep("Expected", 15000), rep("Observed", 3)),
                             allele_frequency = c(BY1_null_frequencies[which(allele_frequencies$SNP_ID == "un_19825"),],
                                                  BY2_null_frequencies[which(allele_frequencies$SNP_ID == "un_19825"),],
                                                  BY3_null_frequencies[which(allele_frequencies$SNP_ID == "un_19825"),],
                                                  un_19825_observed))

#calculate mean null allele frequencies, sd and CI for each year
un_19825_summary <- data_summary(un_19825_trace, varname = "allele_frequency",
                                 groupnames = c("study_year", "type"))

#calculate where 95% of null points lie
confidence <- tapply(un_19825_trace[1:15000,]$allele_frequency, un_19825_trace[1:15000,]$study_year, CI)
un_19825_summary$upper <- c(confidence$`2015`$upper, NA, confidence$`2016`$upper, 
                            NA, confidence$`2017`$upper, NA)
un_19825_summary$lower <- c(confidence$`2015`$lower, NA, confidence$`2016`$lower, 
                            NA, confidence$`2017`$lower, NA)
un_19825_summary

ggplot(un_19825_summary, aes(x = study_year, y = allele_frequency, group = type, colour = type)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4, colour = NA) +
  geom_line() +
  scale_colour_manual(values = c("blue", "darkred")) +
  xlab(NULL) +
  ylab("Allele frequency (SNP un_19825)") +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = c(0.8, 0.8),
        legend.text = element_text(size = 12),
        axis.title= element_text(size = 14),
        axis.text = element_text(size = 12))

