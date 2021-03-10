rm(list=ls())

library(sequoia)
library(kinship2)
library(visPedigree)
library(SNPRelate)
library(reshape2)
library(ggplot2)
library(MASS)
#setwd("~/Documents/Thesis_analysis/")

#### ESTIMARE LINKAGE DISEQUILIBRIUM ####

##read in the data and inspect
linkage_disequilibrium0.45 <- read.table("data/populations_out/duplicates_merged/0.45/plink_out/filtered_snps.ld", header = T, stringsAsFactors = F)
head(linkage_disequilibrium0.45)

linkage_disequilibrium0.05 <- read.table("data/populations_out/duplicates_merged/0.05/plink_out/filtered_snps.ld", header = T, stringsAsFactors = F)
head(linkage_disequilibrium0.05)
#sort by r2 values and inspect to make sure what you expect is there
linkage_disequilibrium0.45 <- linkage_disequilibrium0.45[order(-linkage_disequilibrium0.45$R2),]
head(linkage_disequilibrium0.45)
tail(linkage_disequilibrium0.45)

linkage_disequilibrium0.05 <- linkage_disequilibrium0.05[order(-linkage_disequilibrium0.05$R2),]
head(linkage_disequilibrium0.05)
tail(linkage_disequilibrium0.05)

length(which(linkage_disequilibrium0.45$R2 > 0.035))
length(which(linkage_disequilibrium0.05$R2 > 0.035))
##plot a histogram of r2 values to decide on a threshold to exclude snps in the tails of the distribution
#MAF = 0.45
ggplot(linkage_disequilibrium0.45, aes(x = R2)) +
  geom_histogram(binwidth = 0.005, colour = "darkblue", fill = "lightblue") + 
  #scale_y_log10() + #scale the y axis optional
  labs(x = expression(paste("r"^"2")), y = "count") +
  xlim(c(0,0.25))#limit the x axis to 0.25 optional

#MAF = 0.05
ggplot(linkage_disequilibrium0.05, aes(x = R2)) +
  geom_histogram(binwidth = 0.005, colour = "darkblue", fill = "lightblue") + 
  #scale_y_log10() + #scale the y axis optional
  labs(x = expression(paste("r"^"2")), y = "count") +
  xlim(c(0,0.25))#limit the x axis to 0.25 optional

#from all histograms, r2 of 0.035 chosen to exclude SNP pairs in the tail of the distribution
#how many SNP pairs have r2 > 0.035?
length(which(linkage_disequilibrium0.45$R2 > 0.03)) #71 pairs - MAF = 0.45
length(which(linkage_disequilibrium0.05$R2 > 0.03)) #71 pairs - MAF = 0.05

#write a list of SNP pairs to a new .ld file for plink to exclude SNPS above the chosen threshold
linkage_excluded0.45 <- linkage_disequilibrium0.45[which(linkage_disequilibrium0.45$R2 > 0.03),]
linkage_included0.45 <- linkage_disequilibrium0.45[which(linkage_disequilibrium0.45$R2 <= 0.03),]

linkage_excluded0.05 <- linkage_disequilibrium0.05[which(linkage_disequilibrium0.05$R2 > 0.03),]
linkage_included0.05 <- linkage_disequilibrium0.05[which(linkage_disequilibrium0.05$R2 <= 0.03),]

#n\how many unique snps are there if we take the last 400 pairs after calculating r2?
linkage_included0.45_reduced <- tail(linkage_disequilibrium0.45, 400) #gives 701 unique SNPS
linkage_included0.05_reduced <- tail(linkage_disequilibrium0.05, 400) #gives 790 unique SNPS
#create a list of SNP IDs to keep in and remove
prune.out0.45 <- data.frame(SNP_IDs = c(as.character(linkage_excluded0.45$SNP_A), as.character(linkage_excluded0.45$SNP_B)))#to exclude using PLINK
prune.in0.45 <- data.frame(SNP_IDs = c(as.character(linkage_included0.45$SNP_A), as.character(linkage_included0.45$SNP_B)))#to include using PLINK
prune.in0.45_reduced <- data.frame(SNP_IDs = c(as.character(linkage_included0.45_reduced$SNP_A), as.character(linkage_included0.45_reduced$SNP_B)))#to include using PLINK

prune.out0.05 <- data.frame(SNP_IDs = c(as.character(linkage_excluded0.05$SNP_A), as.character(linkage_excluded0.05$SNP_B)))#to exclude using PLINK
prune.in0.05 <- data.frame(SNP_IDs = c(as.character(linkage_included0.05$SNP_A), as.character(linkage_included0.05$SNP_B)))#to include using PLINK
prune.in0.05_reduced <- data.frame(SNP_IDs = c(as.character(linkage_included0.05_reduced$SNP_A), as.character(linkage_included0.05_reduced$SNP_B)))#to include using PLINK
#write a prune.out file to exclude SNPS using PLINK
#write.table(prune.in, "results/r2_list.prune.in", quote = F, col.names = F, row.names = F, sep = "\t")
#write.table(prune.out0.45, "data/populations_out/duplicates_merged/0.45/plink_out/r2_list.prune.out", quote = F, col.names = F, row.names = F, sep = "\t")
#write.table(prune.in0.45_reduced, "data/populations_out/duplicates_merged/0.45/plink_out/filtered_data/r2_list_reduced.prune.in", quote = F, col.names = F, 
#            row.names = F, sep = "\t")

#write a prune.out file to exclude SNPS using PLINK
# write.table(prune.out0.05, "data/populations_out/duplicates_merged/0.05/plink_out/r2_list.prune.out", quote = F, col.names = F, row.names = F, sep = "\t")
# write.table(prune.in0.05_reduced, "data/populations_out/duplicates_merged/0.05/plink_out/filtered_data/r2_list_reduced.prune.in", quote = F, col.names = F,
#            row.names = F, sep = "\t")

#######create another vcf file of just the reduced pruned SNP list in plink

#### ESTIMATE IBD ####

#create a binary .GDS file from the pruned and filtered .VCF file
#MAF 0.45
# snpgdsVCF2GDS("data/populations_out/duplicates_merged/0.45/plink_out/filtered_data/r2_pruned_filtered.vcf", 
#               "data/populations_out/duplicates_merged/0.45/plink_out/filtered_data/r2_prune_filtered.gds")

#MAF 0.05
# snpgdsVCF2GDS("data/populations_out/duplicates_merged/0.05/plink_out/filtered_data/r2_pruned_filtered.vcf",
#               "data/populations_out/duplicates_merged/0.05/plink_out/filtered_data/r2_prune_filtered.gds")

genofile0.45 <- snpgdsOpen("data/populations_out/duplicates_merged/0.45/plink_out/filtered_data/r2_prune_filtered.gds")
genofile0.05 <- snpgdsOpen("data/populations_out/duplicates_merged/0.05/plink_out/filtered_data/r2_prune_filtered.gds")

start_time <- Sys.time()
IBD_0.45 <- snpgdsIBDMLE(genofile0.45, kinship = T, autosome.only = F, num.thread = 8)
end_time <- Sys.time()
end_time - start_time

names(IBD_0.45)
head(IBD_0.45$kinship)

start_time <- Sys.time()
IBD_0.05 <- snpgdsIBDMLE(genofile0.05, kinship = T, autosome.only = F, num.thread = 8)
end_time <- Sys.time()
end_time - start_time

names(IBD_0.05)
head(IBD_0.05$kinship)
#turn the kinship matrix into a relatedness matrix (kinship coeff to relatedness coff)

relatedness_matrix0.45 <- IBD_0.45$kinship*2
colnames(relatedness_matrix0.45) <- IBD_0.45$sample.id
row.names(relatedness_matrix0.45) <- IBD_0.45$sample.id

relatedness_matrix0.05 <- IBD_0.05$kinship*2
colnames(relatedness_matrix0.05) <- IBD_0.05$sample.id
row.names(relatedness_matrix0.05) <- IBD_0.05$sample.id

#create a data framw from the relatedness matrix to plot a heatmap
upper_triangle <- upper.tri(relatedness_matrix0.05)#new logical matrix of only the upper triangle
#relatedness.upperTriangle <- relatedness_matrix #take a copy of the original relatedness matrix
relatedness_matrix0.05[!upper_triangle]<-NA#set everything not in upper triangle to NA
relatedness_melted<-na.omit(melt(relatedness_matrix0.05, value.name ="relatednessCoeff")) #use melt to reshape the matrix into a along format df
names(relatedness_melted)[3] <- "relatednessCoeff"
relatedness_melted <- relatedness_melted[order(relatedness_melted$relatednessCoeff, decreasing = T),] #sort by descending relatedness
names(relatedness_melted) <- c("IID1", "IID2", "genomic_relatedness")
#identify and save duplicates as a separate df
duplicated_samples <- relatedness_melted[relatedness_melted[,3] >= 0.9,]
colnames(duplicated_samples) <- c("sampleID_1", "sampleID_2", "relatednessCoefficient")

p <- ggplot(relatedness_melted, aes(x = IID1, y = IID2, fill = genomic_relatedness)) + 
  geom_tile() +
  theme(axis.title = element_blank())
p


# write.table(relatedness_melted, "results/relatedness_results_duplicatesMerged.txt", sep = "\t", col.names = T, row.names = F, quote = F)
# ggsave("results/relatedness_duplicatesMerged.pdf", p, device = "pdf", width = 500, height = 500, units = "cm", limitsize = F)


#### PEDIGREE RECONSTRUCTION ####

#obtain a list of samples that remain after filtering
samples <- read.table("data/populations_out/duplicates_merged/0.45/plink_out/filtered_data/filtered_samples.csv", header = F, sep = "\t", stringsAsFactors = F)
names(samples) <- "sample_ID"
#read in the demographic data
Af_demography <- read.table("~/Documents/Data/A.flavicollis_demographic_data_duplicates_merged.txt", header = T, sep = "\t", stringsAsFactors = F)

#pull demographic data for samples that remain after filtring
Af_demography <- Af_demography[which(Af_demography$Sample %in% samples$sample_ID == T),]

#make df of demographic data: ID Sex BirthYear
LifeHist <- Af_demography[,c(1,3,8)]
LifeHist$Sex[LifeHist$Sex == "f"] <- 1
LifeHist$Sex[LifeHist$Sex == "m"] <- 2
LifeHist$Sex[is.na(LifeHist$Sex)] <- "other"
head(LifeHist)
nrow(LifeHist)

#convert genotyping data from the PLINK .raw format to a format seqquoia can understand
GenoM0.45 <- GenoConvert("data/populations_out/duplicates_merged/0.45/plink_out/filtered_data/Af_filtered4sequoia_Reduced.raw")


#make sure the sample names match the GenoM matrix for sequoia
#row.names(GenoM0.45) <- paste0("Af_", row.names(GenoM0.45))

#GENERATE THE PEDIGREE
##reduced SNP set - 701 SNPs

#check for obvious errors in the data (impossible assignments), findmayberel has been deprecated, as has maxmismatch.
start_time <- Sys.time()
ped_out0.45_parentage <- sequoia(GenoM = GenoM0.45, LifeHistData = LifeHist, MaxSibIter = 0, Err = 0.0196836)
end_time <- Sys.time()
end_time - start_time


#generate full pedigree
start_time <- Sys.time()
ped_out0.45_full <- sequoia(GenoM = GenoM0.45, LifeHistData = LifeHist, MaxSibIter = 40, Err = 0.0196836, Plot = T)
end_time <- Sys.time()
end_time - start_time

#write the ourput to .txt files
# for (i in 1:length(ped_out0.45_full)){
#   write.table(ped_out0.45_full[i], paste0("results/ped_out0.45_duplicates_merged/", names(ped_out0.45_full[i]), ".txt"), sep = "\t", quote = F, col.names = T, row.names = F)
# }



#### PEDIGREE SUMMARY OUTPUTS ####

#summry SNP data - MAF distribution, missingness rate, medelian error count
snp_summary <- SnpStats(GenoM = GenoM0.45, Ped = ped_out0.45_full$Pedigree)
#number of parents and grandparents assigned to genotyped individuals
ped_summary <- SummarySeq(SeqList = ped_out0.45_full, Ped = ped_out0.45_full$Pedigree)

#summary of maybe relatives
table(ped_out0.45_full$MaybeRel$TopRel)
table(ped_out0.45_full$MaybePar$TopRel)

#count number of offspring for each individual
n_offspring <- data.frame(sample = names(c(table(ped_out0.45_full$Pedigree$dam), table(ped_out0.45_full$Pedigree$sire))), 
                          offspring = c(table(ped_out0.45_full$Pedigree$dam), table(ped_out0.45_full$Pedigree$sire)))

zero_offspring <- ped_out0.45_full$Pedigree$id[!(ped_out0.45_full$Pedigree$id %in% ped_out0.45_full$Pedigree$dam) & !(ped_out0.45_full$Pedigree$id %in% ped_out0.45_full$Pedigree$sire)]
n_offspring <- rbind(n_offspring, data.frame(sample = zero_offspring, offspring = rep(0, length(zero_offspring))))
rm(zero_offspring)

#mean offspring overall
mean(n_offspring$offspring)
sd(n_offspring$offspring)

#mean offspring when excluding 0
mean(n_offspring[n_offspring$offspring > 0,]$offspring)
sd(n_offspring[n_offspring$offspring > 0,]$offspring)

##does lifetime reproductive success differ between the sexes??
#make a df by merging the full and tidy pedigree with n_offspring data
n_offspring1 <- n_offspring
names(n_offspring1)[1] <- "Ind"
n_offspring1 <- merge(ped_full_tidy, n_offspring1, by = "Ind")

#male ped
male_ped <- n_offspring1[n_offspring1$Sex == "male",]
female_ped <- n_offspring1[n_offspring1$Sex == "female",]

#wilcoxon sum rank test
wilcox.test(male_ped[male_ped$offspring > 0,]$offspring,
            female_ped[female_ped$offspring > 0,]$offspring)

#merge n_offspring with 
#order by no. of offspring for a nice plot
n_offspring <- n_offspring[order(n_offspring$offspring, decreasing = T),]
n_offspring$sample_no <- c(1:nrow(n_offspring))

ggplot(n_offspring, aes(x = offspring)) +
  geom_histogram(colour = "darkblue", fill = "lightblue", binwidth = 1) +
  xlab("Lifetime reproductive success") +
  ylab("Count") +
  xlim(c(-1,30)) +
  theme_bw()

ggplot(n_offspring, aes(x = sample_no, y = offspring)) +
  geom_bar(colour = "midnightblue", fill = "midnightblue", stat = "identity") +
  xlab("Sample") +
  ylab("Lifetime reproductive success") +
  theme_bw() +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  scale_y_continuous(breaks = seq(0,30, by = 5))

#### CONFIDENCE IN PEDIGREES ####

##simulate pedigrees to estimate assignment error rate

start_time <- Sys.time()
set.seed(1)
confidence <- EstConf(Ped = ped_out0.45_full$PedigreePar[,1:3], LifeHistData = LifeHist, nSim = 60, 
                       args.sim = list(nSnp = 701, ParMis = c(0.4,0.4), SnpError = 0.0196836), args.seq = list(Err = 0.0196836), return.PC = T)
end_time <- Sys.time()
end_time - start_time

##calculate the mean assignment error rate after 50 iterations of simulated pedigrees##
##calculated as Number of Mismatches + P2 only / 2N##
##this is inherently biased as is calculated by simulating genotypes absed on the newly generated pedigree.##
##if this new pedigree isnt corret, the second pedigree that compares against it will be inherently biased.##
#overall mean
overall_SimCounts <- confidence$SimCounts[,"TT", , "sire"] + confidence$SimCounts[,"TT", , "dam"]
mean((overall_SimCounts[,3] + overall_SimCounts[,5]) / (2*nrow(ped_out0.45_full$Pedigree)))

#overall mean for dams
mean((confidence$SimCounts[,"TT", , "dam"][,3] + confidence$SimCounts[,"TT", , "dam"][,5]) / (2*nrow(ped_out0.45_full$Pedigree)))

#overall mean for sires
mean((confidence$SimCounts[,"TT", , "sire"][,3] + confidence$SimCounts[,"TT", , "sire"][,5]) / (2*nrow(ped_out0.45_full$Pedigree)))

##CALCULATE THE PAIRWISE PEDIGREE RELATEDNESS
#clean the pedigree to make it compatible with kinship2
#code for the PedPolish function was obtained from https://github.com/cran/sequoia/blob/master/R/PedPolish.R
#as it is not available in the version of sequoia ued for this analysis
Ped_kinship <- PedPolish(ped_out0.45_full$Pedigree, FillParents = T) #adds fake parents when only when has been assigned in the ped
names(LifeHist)[1] <- "id" #change nme of id column to match kinships requirements
Ped_kinship <- merge(Ped_kinship, LifeHist, by.x = "id", all.x = T)
Ped_kinship$Sex[which(Ped_kinship$Sex == 1)] <- "female" #change sed id for females
Ped_kinship$Sex[which(Ped_kinship$Sex == 2)] <-  "male" #change sex id for males
#Ped_kinship$Sex[which(Ped_kinship$Sex == "female")] <- 2 #change sex id for females to 2
#Ped_kinship$Sex[which(Ped_kinship$Sex == "male")] <- 1 #change sex id for females to 1

#fix sexes of parents for kinship2
Ped_kinship <- fixParents(id = Ped_kinship$id, dadid = Ped_kinship$sire, momid = Ped_kinship$dam, sex = Ped_kinship$Sex)

#create a pedigree object for kinship2
Ped_kinship <- pedigree(id = Ped_kinship$id, dadid = Ped_kinship$dadid, momid = Ped_kinship$momid, sex = Ped_kinship$sex)

#CALCULATE THE PAIRWISE PEDIGREE RELATEDNESS MATRIX
kinship_matrix <- kinship(Ped_kinship)#calculate the kinship matrix first
dim(kinship_matrix)#stupidity check

#turn matrix into a df
#create a data framw from the relatedness matrix
upper_triangle <- upper.tri(kinship_matrix)#new logical matrix of only the upper triangle
kinship_matrix[!upper_triangle]<-NA#set everything not in upper triangle to NA
kinship_melted<-na.omit(melt(kinship_matrix, value.name ="relatednessCoeff")) #use melt to reshape the matrix into a along format df
names(kinship_melted)[3] <- "relatednessCoeff"
kinship_melted <- kinship_melted[order(kinship_melted$relatednessCoeff, decreasing = T),] #sort by descending relatedness
names(kinship_melted) <- c("IID1", "IID2", "kinship_ped")

#CALCULATE PEDIGREE RELATEDNESS
kinship_melted$relatedness_ped <- kinship_melted$kinship_ped*2

#merge together the pedigree relatedness and genomic relatedness (obtained from a full set of snps in snprelate) data
kinship_grm <- merge(kinship_melted, relatedness_melted, by = c("IID1", "IID2"))

##PLOT THE DATA
#rough plotting
plot(kinship_grm$genomic_relatedness ~ kinship_grm$relatedness_ped)
grm_ped <- lm(kinship_grm$genomic_relatedness ~ kinship_grm$relatedness_ped)
summary(grm_ped)
abline(grm_ped)

##scale the colour of the points according to the density
#define a function to calculate the kernel density
# Get density of points in 2 dimensions.
# @param x A numeric vector.
# @param y A numeric vector.
# @param n Create a square n by n grid to compute density.
# @return The density within each square.
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, h = c(ifelse(bandwidth.nrd(x) == 0, 0.0001, bandwidth.nrd(x)),
                                  ifelse(bandwidth.nrd(y) == 0, 0.0001, bandwidth.nrd(y))), ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

#estimate kernal density
#kinship_grm$density <- get_density(kinship_grm$relatedness_ped, kinship_grm$relatedness_ped, n = 5000)

#PLOT RESULTS
#pdf("results/grm_kinship_corr.pdf")
ggplot(data = kinship_grm, aes(x = relatedness_ped, y = genomic_relatedness)) +
  geom_point(colour = "midnightblue") +
  #geom_jitter(alpha = 0.1, size = 0.8, width = 0.02, height = 0.02) +
  xlab("Pedigree relatedness") +
  ylab("Genomic relatedness") +
  theme_bw() +
  stat_smooth(method = "lm", colour = "darkred")#maybe remove the abline as it can make the plot misleading... it represents a hypothetical prfect correlation between ped and genomic relatedness
#dev.off()
#check correlation between grm and pairwise pedigree relatedness
cor.test(kinship_grm$relatedness_ped, kinship_grm$genomic_relatedness, method = "pearson")#when all individuals with ped relatednes = 0 are removed, correlation is 92%

#estimate the amount of pedigree error using rped - rgrm > 0.2
length(which(kinship_grm$relatedness_ped - kinship_grm$genomic_relatedness >0.2))/length(kinship_grm$relatedness_ped - kinship_grm$genomic_relatedness)