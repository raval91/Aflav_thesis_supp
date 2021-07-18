#### CALCULATING GENEALOGICAL AND GENETIC CONTRIBUTIONS OF EACH MOUSE IN A PEDIGREE ####

# setwd("~/Documents/Thesis_analysis/")
rm(list = ls())

library(visPedigree)
library(kinship2)
library(optiSel)
library(ggplot2)
library(sequoia)
library(reshape2)

#read in the pedigree
AF_pedigree <- read.table("results/pedigrees/R_v3.6/ped_out0.45_duplicates_merged/pedigree_tidy.txt", 
                       header = T, stringsAsFactors = F, sep = "\t")
AF_pedigree <- AF_pedigree[order(AF_pedigree$Ind),]

demographic_data <- read.table("data/A.flavicollis_demographic_data_duplicates_merged.txt", 
                             header = T, stringsAsFactors = F, sep = "\t")
demographic_data <- demographic_data[,c(1,6)]

#format pedigree for kinship2
pedigree1 <- PedPolish(AF_pedigree[,c(1:3)], GenoNames = F, ZeroToNA = F, NAToZero = F, DropNonSNPd = F, FillParents = T)

names(pedigree1) <- names(AF_pedigree[c(1,3,2)])
rownames(pedigree1) <- NULL
pedigree1 <- pedigree1[-which(pedigree1 == F),]
pedigree1$Sex <- c(AF_pedigree$Sex, rep("female", 23), rep("male", 15))
pedigree1$fID <- makefamid(pedigree1$Ind, pedigree1$Sire, pedigree1$Dam)

final_ped <- pedigree(id = pedigree1$Ind, dadid = pedigree1$Sire, momid = pedigree1$Dam, sex = pedigree1$Sex)

# format the pedigree for visped
pedigree2 <- as.data.frame(final_ped)#pedigree1[,c(1,3,2,4)]
final_visped <- tidyped(pedigree2, trace = "down")#tidyped(pedigree1)
# final_visped <- final_visped[,c(1,3,2,4)]
# names(final_visped)[1:4] <- c("Ind", "Sire", "Dam", "Sex")

#PLOT PEDIGREE
plot.pedigree(final_ped, id = rep(NA, length(final_ped$id)))#kinship2
visped(final_visped, cex = 0.5, compact = T)#visped

#### USE VISPED TO ESTIMATE THE NUMBER OF GENEALOGICAL DESCENDENTS FOR FOUNDERS ONLY ####

#create column to indicate if an individual is a founder
founders_ped <- data.frame(final_visped)
founders_ped$founder <- NA
founders_ped$founder[which(is.na(founders_ped$Sire) == T & is.na(founders_ped$Dam) == T)] <- T
founders_ped$founder[which(is.na(founders_ped$founder == T))] <- F

#list to write descendents of founders to
gen_contributions <- vector(mode = "list", length = 1)

#list of founders from full pedigree
founders_list_full <- founders_ped[founders_ped$founder == T,]$Ind

#extract founders pedigrees
for(i in 1:length(founders_list_full)){
  gen_contributions[[i]] <- try(tidyped(founders_ped, cand = founders_list_full[i], trace = "down"))
}

# #extract founders pedigrees
# for(i in 1:nrow(founders_ped)){
#   gen_contributions[[i]] <- try(tidyped(founders_ped, cand = founders_ped$Ind[i], trace = "down"))
# }

##remove elements of the list with errors due to no descendents
#list of objects with length 1 (just error messages)
failed_peds <- sapply(gen_contributions, function(x) class(x) == c("data.table", "data.frame"))

#subset for only objects length >1 (vector of descendents in each element, first element of each is the founder)
gen_contributions <- gen_contributions[failed_peds[1,]]

#list of founders
founders_list <- as.character(sapply(gen_contributions, function(x) x[1,1]))

#length of each list element -1 (as first of each vector is the founder)
contributions <- sapply(gen_contributions, function(x)  nrow(x) - 1)

#create df of founders genealogical contributions
descendants <- data.frame(founders = founders_list, descendants = contributions)
descendants <- descendants[order(descendants$descendants, decreasing = T),]
descendants$sample <- c(1:nrow(descendants))

#add generation number to descendants
founders <- founders_ped[founders_ped$Ind %in% descendants$founders,]#founders only data
names(founders)[1] <- "founders"
descendants <- merge(descendants, founders, by = "founders")


#### GENETIC CONTRIBUTIONS OF FOUNDERS ####

# genetic contribution to a given generation is the number of alleles 
# inherited by all descendants from the focal individual as a proportion 
# of all alleles in the population
# overall genetic contribution is the total contribution over all generations
# from the focal individual

# table of when founders first appeared in dataset
founders <- read.table("results/pedigrees/R_v3.6/ped_out0.45_duplicates_merged/founders.txt", 
                       header = T, sep = "\t", stringsAsFactors = F)

# compute a relatedness matrix based on the pedigree
relatedness <- 2*kinship(final_ped)

# list to write the relatedness of descendants to founders
relatedness2founder <- vector(mode = "list", length = 1)

## extract the relatedness of each descendant to the founder in the pedigree 
#for every pedigree from every founder
for(i in 1:length(gen_contributions)){
  # vector to write to
  r <- c()
  # for each descendant
  for(j in 2:nrow(gen_contributions[[i]])){
    # extract the relatedness from the relatedness matrix
    r[j] <- relatedness[which(rownames(relatedness) == gen_contributions[[i]]$Ind[1]),
                                            which(colnames(relatedness) == gen_contributions[[i]]$Ind[j])]
  }
  # write the vector of relatedness to a list
  relatedness2founder[[i]] <- as.numeric(na.omit(r))
  # relatedness2founder[[i]] <- as.numeric(r)
}

# assign founder names
names(relatedness2founder) <- founders_list

## calculate expected genetic contributions of founders
# list to write contributions for each founder to
founder_genetic_contrib <- vector("list")

# for every founders descendants
for(i in 1:length(relatedness2founder)){
  #check the generation the founder was assigned (generation used as many birth years are missing)
  generation <- founders[which(founders$founders == names(relatedness2founder)[i]),]$Gen# generation <- founders[which(founders$founders == names(relatedness2founder)[i]),]$First_trapped
  contribution <- c()# contribution <- vector(mode = "list", length = 3)
  if(generation == 1){
    # contribution[1] <- sum(as.numeric(relatedness2founder[[i]]))/334
    contribution[2] <- sum(as.numeric(relatedness2founder[[i]]))/98
    contribution[3] <- sum(as.numeric(relatedness2founder[[i]]))/173
    contribution[4] <- sum(as.numeric(relatedness2founder[[i]]))/90
  } else if(generation == 2) {
    # contribution[2] <- sum(as.numeric(relatedness2founder[[i]]))/98
    contribution[3] <- sum(as.numeric(relatedness2founder[[i]]))/173
    contribution[4] <- sum(as.numeric(relatedness2founder[[i]]))/90
  } else {
    contribution[4] <- sum(as.numeric(relatedness2founder[[i]]))/90
  }
  founder_genetic_contrib[[i]] <- contribution
}

# write genetic contributions to dataframe and tidy
genetic_contributions <- data.frame(do.call(rbind, founder_genetic_contrib))
genetic_contributions$founder <- founders_list
genetic_contributions <- genetic_contributions[,c(5,1,2,3,4)]
names(genetic_contributions) <- c("founders", "gen1", "gen2", "gen3", "gen4")
genetic_contributions[genetic_contributions == "NULL"] <- NA
genetic_contributions$gen1[1:10] <- rep(0, 10)#add 0s for gen 1 mice whcih contribute to gen 2 so plots start from gen 1

#### GENEALOGICAL CONTRIBUTIONS OF FOUNDERS ####
## given as a prortion of the subsequent generation

# list to write contributions for each founder to
founder_genealogical_contrib <- vector("list")

# for every founders descendants
for(i in 1:length(relatedness2founder)){
  #check the generation the founder was assigned (generation used as many birth years are missing)
  generation <- founders[which(founders$founders == names(relatedness2founder)[i]),]$Gen
  contribution <- c()
  if(generation == 1){#if first generation...
    #calculate the cumulative genealogical contribution as a proportion of the total births in the next generation
    contribution[1] <- nrow(gen_contributions[[i]][which(gen_contributions[[i]]$Gen == 1 & gen_contributions[[i]]$Cand == F),])/334
    contribution[2] <- (nrow(gen_contributions[[i]][which(gen_contributions[[i]]$Gen == 2 & gen_contributions[[i]]$Cand == F),])/98) + contribution[1]
    contribution[3] <- (nrow(gen_contributions[[i]][which(gen_contributions[[i]]$Gen == 3 & gen_contributions[[i]]$Cand == F),])/173) + contribution[2]
    contribution[4] <- (nrow(gen_contributions[[i]][which(gen_contributions[[i]]$Gen == 3 & gen_contributions[[i]]$Cand == F),])/90) + contribution[3]
  } else if(generation == 2) {
    contribution[2] <- nrow(gen_contributions[[i]][which(gen_contributions[[i]]$Gen == 2 & gen_contributions[[i]]$Cand == F),])/98
    contribution[3] <- (nrow(gen_contributions[[i]][which(gen_contributions[[i]]$Gen == 3 & gen_contributions[[i]]$Cand == F),])/173) + contribution[2]
    contribution[4] <- (nrow(gen_contributions[[i]][which(gen_contributions[[i]]$Gen == 3 & gen_contributions[[i]]$Cand == F),])/90) + contribution[3]
  } else if(generation == 3) {
    contribution[3] <- nrow(gen_contributions[[i]][which(gen_contributions[[i]]$Gen == 3 & gen_contributions[[i]]$Cand == F),])/173
    contribution[4] <- (nrow(gen_contributions[[i]][which(gen_contributions[[i]]$Gen == 3 & gen_contributions[[i]]$Cand == F),])/90) + contribution[3]
  } else {
    contribution[4] <- (nrow(gen_contributions[[i]][which(gen_contributions[[i]]$Gen == 3 & gen_contributions[[i]]$Cand == F),])/90)
  }
  founder_genealogical_contrib[[i]] <- contribution
}

#write genealogical contributions to dataframe and tidy
genealogical_contributions <- data.frame(do.call(rbind, founder_genealogical_contrib))
genealogical_contributions$founder <- founders_list
genealogical_contributions <- genealogical_contributions[,c(5,1,2,3,4)]
names(genealogical_contributions) <- c("founders", "gen1", "gen2", "gen3", "gen4")


## tidy and merge to single df
# change to long format for easy analysis, and combine dfs
genealogical_contributions <- melt(genealogical_contributions, id.vars = "founders")
genetic_contributions <- melt(genetic_contributions, id.vars = "founders")

genealogical_contributions$type <- rep("genealogical", nrow(genealogical_contributions))
genetic_contributions$type <- rep("genetic", nrow(genetic_contributions))
contribution <- rbind(genealogical_contributions, genetic_contributions)


#tidy df
contribution <- contribution[,c(1,4,2,3)]
colnames(contribution)[3:4] <- c("generation", "contribution")
contribution <- merge(contribution, descendants[,c("founders", "sex")], by = "founders") 


### #write results to text files ####
# write.table(contribution, "results/pedigrees/R_v3.6/ped_out0.45_duplicates_merged/genealogical_genetic_contributions.txt",
#             col.names = T, row.names = F, quote = F)
# write.table(genetic_contributions, "results/pedigrees/R_v3.6/ped_out0.45_duplicates_merged/genetic_contributions.txt",
#             col.names = T, row.names = F, quote = F)
# write.table(genealogical_contributions, "results/pedigrees/R_v3.6/ped_out0.45_duplicates_merged/genealogical_contributions.txt",
#             col.names = T, row.names = F, quote = F)

#### DATA SUMMARY AND ANALYSIS ####

#number of genealogical descendants for each founder
summary(descendants$descendants)
sd(descendants$descendants)

#how many founders caught in each generation
table(descendants$Gen)

#mean genealogical contributions for founders in each generation
tapply(descendants$descendants, descendants$Gen, mean)
tapply(descendants$descendants, descendants$Gen, sd)

#compare fitness as no. of descendants of males and females
wilcox.test(descendants[descendants$sex == "male",]$descendants, 
            descendants[descendants$sex == "female",]$descendants)

## compare fitness as proportion of genetic and genealogical contribution of males and females to the population
#genealogical contribution
wilcox.test(contribution[contribution$type == "genealogical" & contribution$sex == "male",]$contribution,
            contribution[contribution$type == "genealogical" & contribution$sex == "female",]$contribution)

#genetic contribution
wilcox.test(contribution[contribution$type == "genetic" & contribution$sex == "male",]$contribution,
            contribution[contribution$type == "genetic" & contribution$sex == "female",]$contribution)

#Are genealogical and genetic contributions to future generations correlated as in Chen et al (2019)?
cor.test(genealogical_contributions$value, genetic_contributions$value)# YES... highly correlated

## is genealogical or genetic contribution to the population larger?
#this is a paired test (signed rank test) as it is the same individuals with 2 differenr measures of fitness
wilcox.test(contribution[contribution$type == "genealogical",]$contribution,
            contribution[contribution$type == "genetic",]$contribution,
            paired = T)#YES... genealogical contributions are greater than genetic

boxplot(contribution$contribution ~ contribution$type)
summary(contribution[contribution$type == "genealogical",]$contribution)
summary(contribution[contribution$type == "genetic",]$contribution)

#### PLOTTING ####

#number of genealogical descendants
ggplot(descendants, aes(x = sample, y = descendants)) +
  geom_bar(colour = "midnightblue", fill = "midnightblue", stat = "identity") +
  xlab("Sample") +
  ylab("Genealogical descendants") +
  theme_bw() +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  scale_y_continuous(breaks = seq(0,40, by = 5))

#subset genealogical and genetic contributions for founder Af_237 and Af_258
Af_contribs <- contribution[contribution$founders == "Af_237" | contribution$founders == "Af_258",]

#genealogical and genetic contributions of 2 founders to each generation
ggplot(Af_contribs, aes(x = generation, y = contribution, group = type)) +
  geom_point(colour = "midnightblue") +
  geom_line(aes(linetype = type), colour = "midnightblue") +
  xlab("Generation") +
  ylab("Contribution") +
  facet_wrap( ~ founders, nrow = 1) +
  theme_bw() +
  scale_x_discrete(labels = c("1", "2", "3", "4")) +
  theme(legend.position = c(0.8, 0.8),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text.x = element_text(size = 14))

#subset the pedigrees for Af_237 and Af_258, and plot using the optiSel package
#reorder pedigree1 to have the columns in the following format id,sire,dam
new_pedigree1 <- pedigree1[,c(1,3,2,4,5)]

#ensure the format is correcct for the whole pedigree
new_pedigree1 <- subPed(new_pedigree1, succGen = 4, keep = new_pedigree1$Ind)

#subset the pedigree for Af_237 and its descendants
Af_237Kinship <- subPed(Pedig = new_pedigree1, keep = "Af_237", succGen = 4)

subset_pedigree <- subPed(Pedig = new_pedigree1, keep = c("Af_237", "Af_258"), succGen = 4)

#plot the pedigree using optiSel
pedplot(Af_237Kinship)#Af_237
pedplot(subset_pedigree, label = NULL)
