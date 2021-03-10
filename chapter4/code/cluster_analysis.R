#### CLUSTER ANALYSIS OF HI DATA ####


rm(list = ls())

#setwd("~/Documents/Thesis_analysis")

library(reshape2)
library(ggplot2)
library(ggpubr)
library(car)
library(dplyr)
library(cluster)
library(factoextra)
library(FactoMineR)
library(NbClust)
library(bestNormalize)

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

#transform HI1f
bestNormalize(HI_data$HI1f)#which transfrmation? orderNorm transformation
# trans <- orderNorm(jitter(HI_data$HI1f))#random noise added to ensure no ties in data resulting in true gaussian distribution
trans <- orderNorm(HI_data$HI1f)#random noise added to ensure no ties in data resulting in true gaussian distribution
HI_data$HI_orderNorm <- trans$x.t

#create a df with no missing values and is scaled to all variables can be directly compared
cluster_data <- HI_data[!is.na(HI_data$BMR_min),c(1:3,5,9,12,14)]#remove nas
cluster_data[,4:7] <- scale(cluster_data[,4:7])#scale the data
rm(BMR_data, trans)

#### ASSESS CLUSTERING TENDENCY ####

#create a null dataset with no structure generated from our data to compare against
#data should alredy be scaled/standardised!
set.seed(1)
random_df <- apply(cluster_data[,c(4:6)], 2, function(x){runif(length(x), min(x), (max(x)))})
random_df <- as.data.frame(random_df)

#plot the data
fviz_pca_ind(prcomp(cluster_data[,c(4:6)]), 
             habillage = cluster_data$Sex,
             geom = "point",
             ggtheme = theme_bw())

fviz_pca_ind(prcomp(random_df), 
             geom = "point",
             ggtheme = theme_bw())

#compute the clustering tendency - threshold to reject null = 0.5
HI_tendency <- get_clust_tendency(cluster_data[,c(4:6)], 
                                  n = nrow(cluster_data[,c(4:6)])-1, graph = T)#real data
rand_tendency <- get_clust_tendency(random_df, 
                                    n = nrow(random_df)-1, graph = T)#random data
HI_tendency$hopkins_stat# H = 0.869... high clustering tendency
rand_tendency$hopkins_stat# H = 0.496... low clustering tendency

#compute euclidean distance metric
dist_euc <- dist(cluster_data[,c(4:6)], method = "euclidean")

#visualise the dissimilarity matrix based on euclidean distance
fviz_dist(dist_euc, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))

#### K-MEANS CLUSTERING ####

##estimate the optimal number of clusters
set.seed(1)

#elbow method
fviz_nbclust(cluster_data[,c(4:6)], kmeans, method = "wss") +#not clear what is best
  theme_bw()

#silhouette method
fviz_nbclust(cluster_data[,c(4:6)], kmeans, method = "silhouette") +#not clear what is best
  theme_bw()

#gap statistic method
gap_stat <- clusGap(cluster_data[,c(4:6)], FUN = kmeans, nstart = 25,
                    K.max = 10, B = 50)
# Print the result
print(gap_stat, method = "firstmax")
fviz_gap_stat(gap_stat)

# estimate the optimal number of clusters based on 30 different metrics
#3 clusters chosen by the function BUT PLOTTING IT LOOKS STUPID. like no
#structure exists at all... with k=2, it looks liek there is some structure.
n_clust <- NbClust(data = cluster_data[,c(4:6)], distance = "euclidean", min.nc = 2, max.nc = 10, method = "kmeans")

#k-means clustering
k2 <- kmeans(cluster_data[,c(4:6)], centers = 2, nstart = 25)
k3 <- kmeans(cluster_data[,c(4:6)], centers = 3, nstart = 25)
k4 <- kmeans(cluster_data[,c(4:6)], centers = 4, nstart = 25)
cluster_data$k2 <- k2$cluster
cluster_data$k3 <- k3$cluster
cluster_data$k4 <- k4$cluster
  
#plot
fviz_cluster(k2, data = cluster_data[,c(4:6)],
             main = NULL,
             geom = "point",
             shape = 16,
             #alpha = 0,
             show.clust.cent = F,
             ellipse.type = "norm",
             ggtheme = theme_bw()
             ) +
  scale_colour_discrete(name = "Thermal strategy", labels = c("Generalists", "Specialists")) +
  scale_fill_discrete(guide = F) +
  geom_hline(yintercept = 0, lty = "dotted") +
  geom_vline(xintercept = 0, lty = "dotted") +
  theme(legend.position = c(0.8, 0.1)) #+
  #geom_point(data = cluster_data, aes(x = k2.pc1, y = k2.pc2, colour = k2, shape = Sex)) 

fviz_pca_biplot(PCA(cluster_data[,4:6]))

#how many males and females belong to each cluster?
table(cluster_data[cluster_data$Sex == "Male",]$k2)#k2 males
table(cluster_data[cluster_data$Sex == "Female",]$k2)#k2 females

#how many males and females were specialists or generalists?
generalists <- names(table(cluster_data[cluster_data$k2 == 1,]$merged_samples))#sample names of speciaists
specialists <- names(table(cluster_data[cluster_data$k2 == 2,]$merged_samples))#sample names of generalists

#vector of sexes for specialist and generalist mice
specialist_count <- c()
generalist_count <- c()

#populate the vectors
for(i in 1:length(specialists)){
  specialist_count[i] <- unique(cluster_data[cluster_data$merged_samples == specialists[i], ]$Sex)#specialists
}

for(i in 1:length(generalists)){
  generalist_count[i] <- unique(cluster_data[cluster_data$merged_samples == generalists[i], ]$Sex)#generalists
}

table(specialist_count)
table(generalist_count)

table(cluster_data[cluster_data$Sex == "Male",]$k3)#k3 males
table(cluster_data[cluster_data$Sex == "Female",]$k3)#k3 females

#what is the proportion of variance explained by each principle component?
(eigen(cov(cluster_data[,4:6]))$vectors)^2
get_eigenvalue(PCA(cluster_data[,4:6]))

#correlation of variables
pca <- prcomp(cluster_data[,4:6], scale = F)
pca$rotation
#### k - means clustering on averaged data ####

#order data to make sure nothing is jumbled
HI_data <- HI_data[order(HI_data$merged_samples),]

#average the data
cluster_ave <- tapply(HI_data$HI_orderNorm, HI_data$merged_samples, mean)# HI order norm
mb1 <- tapply(HI_data$mb1, HI_data$merged_samples, mean)# mb
HI1f <- tapply(HI_data$HI1f, HI_data$merged_samples, mean)# HI
BMR_ave <- tapply(HI_data$BMR_min, HI_data$merged_samples, mean, na.rm = T)# BMR_min

#make df
cluster_ave <- data.frame(merged_samples = names(cluster_ave), mean_mb = mb1, HI1f = HI1f, BMR_min = BMR_ave, HI_ave_gaus = cluster_ave)

#add sex and bmr
cluster_ave$Sex <- HI_data[match(unique(HI_data$merged_samples), HI_data$merged_samples),]$Sex 

#scale the data
cluster_ave[,2:5] <- scale(cluster_ave[,2:5])

#k-means clustering
k2_ave <- kmeans(cluster_ave[,2:4], centers = 2, nstart = 25) 
k3_ave <- kmeans(cluster_ave[,2:4], centers = 3, nstart = 25) 
#plot
fviz_cluster(k2_ave, data = cluster_ave[,2:4],
             main = NULL,
             geom = "point",
             ellipse.type = "norm",
             ggtheme = theme_bw()
)

cluster_ave$kmeans_k2 <- k2_ave$cluster

table(cluster_ave[cluster_ave$Sex == "Male",]$kmeans_k2)#k2 males
table(cluster_ave[cluster_ave$Sex == "Female",]$kmeans_k2)#k2 females


#### HIERARCHICAL CLUSTERING ####
## alternative to k-means clustering which does not require the user
## to pre-specify the number of clusters.

#k3 makes no sensse here either. the third group that pops up does not fit with the
#results of the linear regressions showing male and female mb/HI. k2 using the ward 
#method males most biological sense. more males show in cluster 1 than 2. but females 
#are not that dissimilar in which group they belong. males are more specialist, and 
#females are more generalist.

## compute dissimilarity matrices
#euclidean
dist_euc <- dist(cluster_data[,c(4:6)], method = "euclidean")

##hierarchical clustering
#assess clustering strength between different methods
m <- c( "average", "single", "complete", "ward")#clustering methods
names(m) <- c( "average", "single", "complete", "ward")

# function to compute coefficient/strength of clustering
ac <- function(x) {
  agnes(dist_euc, method = x)$ac
}

#calculate the agglomerative coefficient for each method
lapply(m, ac)#Ward found as the most strongly clustered. agg. coeff = 0.99
hc_euc_ward <- agnes(dist_euc, method = "ward")

#plot a dentrogram
plot(hc_euc_ward)
rect.hclust(hclust(dist_euc, method = "ward.D"), k = 2, border = 2:4)

#cut the tree into the groups and add the grouping info to cluster data
sub_grp <- cutree(hc_euc_ward, k = 2)
cluster_data$cluster <- sub_grp

#visualise
fviz_cluster(list(data = cluster_data[,c(4:6)], cluster = sub_grp), 
             main = NULL,
             geom = "point",
             ellipse.type = "norm",
             ggtheme = theme_bw()
             )# +
  #scale_colour_manual(values = c("#6495ED85", "#C1403D85")) +
  #scale_fill_manual(values = c("#6495ED85", "#C1403D85")) +
  #ylab(expression(paste("H"[i]))) +
  #xlab(expression(paste(m[b])))
  

#how many males and females belong to each cluster?
table(cluster_data[cluster_data$Sex == "Male",]$cluster)# males
table(cluster_data[cluster_data$Sex == "Female",]$cluster)# females
             
#as above with averaged data

dist_euc_ave <- dist(cluster_ave[,2:4], method = "euclidean")

# function to compute coefficient/strength of clustering
ac1 <- function(x) {
  agnes(dist_euc_ave, method = x)$ac
}

lapply(m, ac1)#Ward found as the most strongly clustered. agg. coeff = 0.99
hc_euc_ward_ave <- agnes(dist_euc_ave, method = "ward")

#plot a dentrogram
plot(hc_euc_ward_ave)
rect.hclust(hclust(dist_euc_ave, method = "ward.D"), k = 2, border = 2:4)

#cut the tree into the groups and add the grouping info to cluster data
sub_grp_ave <- cutree(hc_euc_ward_ave, k = 2)
cluster_ave$ward_k2 <- sub_grp_ave

#visualise
fviz_cluster(list(data = cluster_ave[,2:3], cluster = sub_grp_ave), 
             main = NULL,
             geom = "point",
             ellipse.type = "norm",
             ggtheme = theme_bw()
)

table(cluster_ave[cluster_data$Sex == "Male",]$ward_k2)# males
table(cluster_ave[cluster_data$Sex == "Female",]$ward_k2)# females
