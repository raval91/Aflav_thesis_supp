#### EXTRACT DATA FROM CONFINT FILES AND CALCULATE PVALS FOR NULL MAF SIMS ####
rm(list=ls())

#setwd("~/Huddersfield/christmas/Untitled Folder/")

#list of confint files

confint_files <- list.files("results/nullMAFsim_out/", pattern = "confint", full.names = T)

#read in observed MAF
allele_frequency_change <- read.table("data/allele_frequency_change/allele_frequency_change.txt", header = T, stringsAsFactors = F)


# FUNCTION TO CALCULATE P VALUES
pvals <- function(){
  
  # list to write pvals to
  pvals <- vector(mode = "list", length = 3)
  
  # temp matrices to write to for each sim and year
  sig_counts_BY1_2 <- matrix(nrow = nrow(allele_frequency_change), ncol = length(confint_files))
  sig_counts_BY2_3 <- matrix(nrow = nrow(allele_frequency_change), ncol = length(confint_files))
  sig_counts_BY3_1 <- matrix(nrow = nrow(allele_frequency_change), ncol = length(confint_files))
  
  
  # for every confint file from null MAF simulations
  for (i in 1:length(confint_files)){
    
    # read file
    print(paste("Working on simulation ", i))
    confint_file <- read.table(confint_files[i], header = T, stringsAsFactors = F)
    confint_file1 <- split(confint_file[confint_file$years == "1_2",3:4], 
                           seq(nrow(confint_file[confint_file$years == "1_2",])))
    confint_file2 <- split(confint_file[confint_file$years == "2_3",3:4], 
                           seq(nrow(confint_file[confint_file$years == "2_3",])))
    confint_file3 <- split(confint_file[confint_file$years == "3_1",3:4], 
                           seq(nrow(confint_file[confint_file$years == "3_1",])))
    
    #for each simulation
    for (j in 1:nrow(allele_frequency_change)){
      
      # print(paste("SNP ", allele_frequency_change[j,]$SNP_ID))
      ## is the observed AF outside 95th percentiles of the simulations? 
      # between year 1-2
      sig_counts_BY1_2[j,i] <- allele_frequency_change[j,2] >= confint_file1[[j]]$upper | 
        allele_frequency_change[j,2] <= confint_file1[[j]]$lower
      
      sig_counts_BY2_3[j,i] <- allele_frequency_change[j,3] >= confint_file2[[j]]$upper | 
        allele_frequency_change[j,3] <= confint_file2[[j]]$lower
      
      sig_counts_BY3_1[j,i] <- allele_frequency_change[j,4] >= confint_file3[[j]]$upper | 
        allele_frequency_change[j,4] <= confint_file3[[j]]$lower
    }
  }
  # calculate the p vals
  p_values <- data.frame(SNP_ID = allele_frequency_change$SNP_ID, 
                         pvals_1_2 = 1-(rowSums(sig_counts_BY1_2)/length(confint_files)), 
                         pvals_2_3 = 1-(rowSums(sig_counts_BY2_3)/length(confint_files)), 
                         pvals_3_1 = 1-(rowSums(sig_counts_BY3_1)/length(confint_files)))
  return(p_values)
}

# CALCULATE P VALUES FOR EACH YEAR
system.time(p_values <- pvals())



