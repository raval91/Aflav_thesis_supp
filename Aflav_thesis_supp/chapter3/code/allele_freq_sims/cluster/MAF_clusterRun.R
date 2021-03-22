#### NULL ALLELE FREQUENCY SIMULATIONS ####

#rm(list=ls())

#set the wd
#setwd("~/Documents/Thesis_analysis/")


# READ IN THE REQUIRED DATA
# MAF for each SNP in each study/birth year
 allele_frequencies <- read.table("/local/u1767986/data/allele_frequencies.txt",
                                  header = T, sep = "\t", stringsAsFactors = F)

# MAF change for each year and SNP
 allele_frequency_change <- read.table("/local/u1767986/data/allele_frequency_change.txt",
                                       header = T, sep = "\t", stringsAsFactors = F)


#### GENERATE NULL ALLELE FREQUENCY DISTRIBUTIONS FOR COMPARISON  WITH OBSERVED ALLELE FREQUENCIES ####

#set iteration variable to automatically generate a random seed
batch <- as.numeric(Sys.getenv("PBS_ARRAYID"))
cline_args <- commandArgs()
iter <- batch*(as.numeric(cline_args[4]) + 20)

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
set.seed(iter)
system.time(null_frequencies <- simulations(samplings = 1000))

#### CALCULATE CALCULATE AF CHANGE BETWEEN EACH YEAR ####

# function to calculate null AF change and extract confidence intervals
confidence_intervals <- function(){
  # empty matrices to write matrix output for each year
  BY1_null_frequencies <- matrix(nrow = nrow(null_frequencies[[1]]), ncol = length(null_frequencies))
  BY2_null_frequencies <- matrix(nrow = nrow(null_frequencies[[1]]), ncol = length(null_frequencies))
  BY3_null_frequencies <- matrix(nrow = nrow(null_frequencies[[1]]), ncol = length(null_frequencies))
  
  CI_1_2 <- vector(mode = "list", length = nrow(allele_frequency_change))
  CI_2_3 <- vector(mode = "list", length = nrow(allele_frequency_change))
  CI_3_1 <- vector(mode = "list", length = nrow(allele_frequency_change))
  
  # for every sampling
  for(i in 1:length(null_frequencies)){
    
    # extract null frequencies for each year and write their own matrix
    BY1_null_frequencies[,i] <- null_frequencies[[i]][,1]
    BY2_null_frequencies[,i] <- null_frequencies[[i]][,2]
    BY3_null_frequencies[,i] <- null_frequencies[[i]][,3]
  }
  # calculate allele frequency change between years
  BY1_2AFchange <- BY2_null_frequencies - BY1_null_frequencies
  BY2_3AFchange <- BY3_null_frequencies - BY2_null_frequencies
  BY3_1AFchange <- BY3_null_frequencies - BY1_null_frequencies
  
  # calculate CI
  for(i in 1:nrow(BY1_2AFchange)) {
    CI_1_2[[i]] <- CI(BY1_2AFchange[i,])
    CI_2_3[[i]] <- CI(BY2_3AFchange[i,])
    CI_3_1[[i]] <- CI(BY3_1AFchange[i,])
  }
  
  # write CI to df 
  CI <- data.frame(SNP_ID = rep(allele_frequencies$SNP_ID, 3), 
                   years = c(rep("1_2", nrow(allele_frequencies)),
                             rep("2_3", nrow(allele_frequencies)),
                             rep("3_1", nrow(allele_frequencies))),
                   upper = c(sapply(CI_1_2, "[[", 1), 
                             sapply(CI_2_3, "[[", 1), 
                             sapply(CI_3_1, "[[", 1)), 
                   lower = c(sapply(CI_1_2, "[[", 2),
                             sapply(CI_2_3, "[[", 2),
                             sapply(CI_3_1, "[[", 2)))
  
  return(CI)
}

system.time(confidence <- confidence_intervals())

write.table(confidence, 
            paste0("/local/u1767986/results/confint_batch", batch, "_sim", cline_args[4], ".txt"),
            col.names = T, row.names = F, quote = F, sep = "\t")

# write.table(confidence, paste0("results/allele_frequency_change/test/confint_", iter, ".txt"),
#             col.names = T, row.names = F, quote = F, sep = "\t")

