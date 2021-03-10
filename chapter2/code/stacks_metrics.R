#### PLOTTING METRICS FOR PARAMETER OPTIMISATION ####
rm(list = ls())

library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(scales)
#setwd("/media/rohan/Maxtor/rohan_parameter_optimisation/results/")

#read in the data
optimisation_metrics <- read.table("parameter_optimisation_metrics.txt", header = T, sep = "\t")#iteration over all parameters
param <- as.character(optimisation_metrics$parameter)
optimisation_metrics$params <- param
coverage <- read.table("../data/m_coverage.txt", header = T, sep = "\t")#coverage data
head(coverage)
r_80 <- read.table("populations/r80_metrics.txt", header = T, sep = "\t")#r80 metrics from populations

#split the parameter column into the parameter and value
optimisation_metrics$parameter <- unlist(lapply(strsplit(param, "", perl = T), `[[`, 1))
optimisation_metrics$parameter_value <- unlist(lapply(strsplit(param, "[mMn]", perl = T), `[[`, 2))
optimisation_metrics <- optimisation_metrics[,c(1,6,2,7,3,4,5)]#reorder the df columns
optimisation_metrics$n_assembled_loci_K <- optimisation_metrics$n_assembled_loci/1000
optimisation_metrics$n_polymorphic_loci_K <- optimisation_metrics$n_polymorphic_loci/1000
optimisation_metrics$n_snps_K <- optimisation_metrics$n_snps/1000

#split the df into 3 for each paramter to make plotting easier
m_optimisation <- optimisation_metrics[optimisation_metrics$parameter == "m",]
M_optimisation <- optimisation_metrics[optimisation_metrics$parameter == "M",]
n_optimisation <- optimisation_metrics[optimisation_metrics$parameter == "n",]

#order the n dataframe so it is in order of values of n
n_optimisation$parameter_value <- as.numeric(n_optimisation$parameter_value)
n_optimisation <- n_optimisation[order(n_optimisation$parameter_value),]
n_optimisation$parameter_value <- as.factor(n_optimisation$parameter_value)#converts the ordered variable into a factor where the levels are in the correct order


#organise the r_80 df so you can plot the values above each boxplot
r80_m <- r_80[r_80$parameter == "m",]
r80_m$parameter_value <- as.factor(c(2:6))
r80_m$assembled_r80_loci_k <- r80_m$assembled_r80_loci/1000
r80_m$polymorphic_r80_loci_k <- r80_m$polymorphic_r80_loci/1000
r80_m$r80_snps_k <- r80_m$r80_snps/1000

r80_M <- r_80[r_80$parameter == "M",]
r80_M$parameter_value <- as.factor(c(0:8))
r80_M$assembled_r80_loci_k <- r80_M$assembled_r80_loci/1000
r80_M$polymorphic_r80_loci_k <- r80_M$polymorphic_r80_loci/1000
r80_M$r80_snps_k <- r80_M$r80_snps/1000

r80_n <- r_80[r_80$parameter == "n",]
r80_n$parameter_value <- as.factor(c(0:10))
r80_n$assembled_r80_loci_k <- r80_n$assembled_r80_loci/1000
r80_n$polymorphic_r80_loci_k <- r80_n$polymorphic_r80_loci/1000
r80_n$r80_snps_k <- r80_n$r80_snps/1000

r80_nM <- r_80[r_80$parameter == "n_M",]
r80_nM$parameter_value <- as.factor(c(0:8))
r80_nM$assembled_r80_loci_k <- r80_nM$assembled_r80_loci/1000
r80_nM$polymorphic_r80_loci_k <- r80_nM$polymorphic_r80_loci/1000
r80_nM$r80_snps_k <- r80_nM$r80_snps/1000
r80_nM <- r80_nM[order(r80_nM$parameter_value),]
####plot the data

##m
#number of assembled loci
m_assembled_loci <- ggplot(m_optimisation, aes(x = parameter_value, y = n_assembled_loci_K, fill = parameter_value)) + 
  geom_boxplot() + 
  theme_bw() +
  xlab("m") +
  ylab("Number of assembled loci (K)") + 
  theme(legend.position = "none", text = element_text(size = 19), axis.text.x = element_text(size = 19), axis.text.y = element_text(size = 19)) + 
  scale_fill_brewer(palette = "Greens") + 
  geom_point(data = r80_m, aes(x = parameter_value, y = assembled_r80_loci_k), shape = 23, fill = "blue", size = 4) + 
  scale_y_continuous(breaks = seq(0,200, by = 10))

#number of polymorphic loci
m_polymorphic_loci <- ggplot(m_optimisation, aes(x = parameter_value, y = n_polymorphic_loci_K, fill = parameter_value)) + 
  geom_boxplot() + 
  theme_bw() +
  xlab("m") +
  ylab("Number of polymorphic loci (K)") + 
  theme(legend.position = "none", text = element_text(size = 19), axis.text.x = element_text(size = 19), axis.text.y = element_text(size = 19)) + 
  scale_fill_brewer(palette = "Greens") + 
  geom_point(data = r80_m, aes(x = parameter_value, y = polymorphic_r80_loci_k), shape = 23, fill = "blue", size = 4) + 
  scale_y_continuous(breaks = seq(0,200, by = 5))

#number of SNPs
m_n_snps <- ggplot(m_optimisation, aes(x = parameter_value, y = n_snps_K, fill = parameter_value)) + 
  geom_boxplot() + 
  theme_bw() +
  xlab("m") +
  ylab("Number of SNPs (K)") + 
  theme(legend.position = "none", text = element_text(size = 19), axis.text.x = element_text(size = 19), axis.text.y = element_text(size = 19)) + 
  scale_fill_brewer(palette = "Greens") + 
  geom_point(data = r80_m, aes(x = parameter_value, y = r80_snps_k), shape = 23, fill = "blue", size = 4) + 
  scale_y_continuous(breaks = seq(0,200, by = 10))

##M
#number of assembled loci
M_assembled_loci <- ggplot(M_optimisation, aes(x = parameter_value, y = n_assembled_loci_K, fill = parameter_value)) + 
  geom_boxplot() + 
  theme_bw() +
  xlab("M") +
  ylab("Number of assembled loci (K)") + 
  theme(legend.position = "none", text = element_text(size = 19), axis.text.x = element_text(size = 19), axis.text.y = element_text(size = 19)) + 
  scale_fill_brewer(palette = "Oranges") + 
  geom_point(data = r80_M, aes(x = parameter_value, y = assembled_r80_loci_k), shape = 23, fill = "blue", size = 4) + 
  scale_y_continuous(breaks = seq(0,200, by = 10))

#number of polymorphic loci
M_polymorphic_loci <- ggplot(M_optimisation, aes(x = parameter_value, y = n_polymorphic_loci_K, fill = parameter_value)) + 
  geom_boxplot() + 
  theme_bw() +
  xlab("M") +
  ylab("Number of polymorphic loci (K)") + 
  theme(legend.position = "none", text = element_text(size = 19), axis.text.x = element_text(size = 19), axis.text.y = element_text(size = 19)) + 
  scale_fill_brewer(palette = "Oranges") + 
  geom_point(data = r80_M, aes(x = parameter_value, y = polymorphic_r80_loci_k), shape = 23, fill = "blue", size = 4) + 
  scale_y_continuous(breaks = seq(0,200, by = 5))

#number of SNPs
M_n_snps <- ggplot(M_optimisation, aes(x = parameter_value, y = n_snps_K, fill = parameter_value)) + 
  geom_boxplot() + 
  theme_bw() +
  xlab("M") +
  ylab("Number of SNPs (K)") + 
  theme(legend.position = "none", text = element_text(size = 19), axis.text.x = element_text(size = 19), axis.text.y = element_text(size = 19)) + 
  scale_fill_brewer(palette = "Oranges") + 
  geom_point(data = r80_M, aes(x = parameter_value, y = r80_snps_k), shape = 23, fill = "blue", size = 4) + 
  scale_y_continuous(breaks = seq(0,200, by = 10))

##n
colourCount <- length(unique(n_optimisation$parameter_value))
getPalette <- colorRampPalette(brewer.pal(9, "Purples"))(colourCount)
#number of assembled loci
n_assembled_loci <- ggplot(n_optimisation, aes(x = parameter_value, y = n_assembled_loci_K, fill = parameter_value)) + 
  geom_boxplot() + 
  theme_bw() +
  xlab("n") +
  ylab("Number of assembled loci (K)") + 
  theme(legend.position = "none", text = element_text(size = 19), axis.text.x = element_text(size = 19), axis.text.y = element_text(size = 19)) + 
  scale_fill_manual(values = getPalette) + 
  geom_point(data = r80_n, aes(x = parameter_value, y = assembled_r80_loci_k), shape = 23, fill = "blue", size = 4) + 
  scale_y_continuous(breaks = seq(0,200, by = 10))

#number of polymorphic loci
n_polymorphic_loci <- ggplot(n_optimisation, aes(x = parameter_value, y = n_polymorphic_loci_K, fill = parameter_value)) + 
  geom_boxplot() + 
  theme_bw() +
  xlab("n") +
  ylab("Number of polymorphic loci (K)") + 
  theme(legend.position = "none", text = element_text(size = 19), axis.text.x = element_text(size = 19), axis.text.y = element_text(size = 19)) + 
  scale_fill_manual(values = getPalette) + 
  geom_point(data = r80_n, aes(x = parameter_value, y = polymorphic_r80_loci_k), shape = 23, fill = "blue", size = 4) + 
  scale_y_continuous(breaks = seq(0,200, by = 5))

#number of SNPs
n_n_snps <- ggplot(n_optimisation, aes(x = parameter_value, y = n_snps_K, fill = parameter_value)) + 
  geom_boxplot() + 
  theme_bw() +
  xlab("n") +
  ylab("Number of SNPs (K)") + 
  theme(legend.position = "none", text = element_text(size = 19), axis.text.x = element_text(size = 19), axis.text.y = element_text(size = 19)) + 
  scale_fill_manual(values = getPalette) + 
  geom_point(data = r80_n, aes(x = parameter_value, y = r80_snps_k), shape = 23, fill = "blue", size = 4) + 
  scale_y_continuous(breaks = seq(0,200, by = 10))

ggarrange(m_assembled_loci + rremove("xlab") + rremove("x.text"), M_assembled_loci + rremove("ylab") + rremove("xlab") + rremove("x.text"), 
          n_assembled_loci + rremove("ylab") + rremove("xlab") + rremove("x.text"),
          m_polymorphic_loci + rremove("xlab") + rremove("x.text"), M_polymorphic_loci + rremove("ylab") + rremove("xlab") + rremove("x.text"), 
          n_polymorphic_loci+ rremove("ylab") + rremove("xlab") + rremove("x.text"),
          m_n_snps, M_n_snps + rremove("ylab"), n_n_snps + rremove("ylab"),
          ncol = 3, nrow = 3)

#ggsave("parameter_optimisation_plots.pdf", width = 15, height = 15, units = "in")

##effect of m on mean coverage and mean merged coverage
ggplot(coverage, aes(x = parameter_value, y = coverage, fill = measure)) + 
  geom_boxplot(outlier.size = 0) + 
  #coord_cartesian(ylim = c(0,70)) +
  ylab("Coverage (x)") + 
  theme(legend.position = "top", axis.text.x = element_text(size = 4),
        axis.title.y = element_text(size = 4), axis.text.y = element_text(size = 4),
        legend.title = element_blank(), legend.text = element_text(size = 4.5)) +
  theme_bw() +
  theme(legend.position = "top", legend.title = element_blank(), axis.title.x = element_blank()) +
  scale_fill_manual(values = c("#339933", "#993399"), labels = c("Mean coverage", "Mean merged coverage")) +
  scale_y_continuous(breaks = seq(0,150, by = 10))

#### PROBE EACH PARAMETER VALUE TO DECIDE WHAT IS THE OPTIMUM ####

## m
#mean coverage for each parameter value
mean_coverage <- coverage[coverage$measure == "mean",]#new df for mean coverage
mean_merged_coverage <- coverage[coverage$measure == "merged",]#new df for mean merged coverage

tapply(mean_coverage$coverage, mean_coverage$parameter_value, FUN = mean, na.rm = T)#stabilises after m3. mean coverage @ m3=21.9, m4=23.06, m5=23.9, m6=24.7
tapply(mean_merged_coverage$coverage, mean_merged_coverage$parameter_value, FUN = mean, na.rm = T)#stabilises after m3. mean merged coverage @ m3=28.5, m4=29.9, m5=31.1, m6=31.9

#number of assembled and polymorphic loci, and number of SNPs
tapply(m_optimisation$n_assembled_loci_K, m_optimisation$params, FUN = mean, na.rm = T)#stabilises after m3. assembled loci @ m3=57.5k, m4=54.8k, m5=52.8k, m6=51.0k
tapply(m_optimisation$n_polymorphic_loci_K, m_optimisation$params, FUN = mean, na.rm = T)
tapply(m_optimisation$n_snps_K, m_optimisation$params, FUN = mean, na.rm = T)

##M
#number of assembled and polymorphic loci, and number of SNPs
tapply(M_optimisation$n_assembled_loci_K, M_optimisation$params, FUN = mean, na.rm = T)#stabilises after m3. assembled loci @ m3=57.5k, m4=54.8k, m5=52.8k, m6=51.0k
tapply(M_optimisation$n_polymorphic_loci_K, M_optimisation$params, FUN = mean, na.rm = T)
tapply(M_optimisation$n_snps_K, M_optimisation$params, FUN = mean, na.rm = T)

#number of new polymorphic r80 loci added at each iteration of M
new_pol_loci <- c()
difference_per_M <- c("M0/M1", "M1/M2", "M2/M3", "M3/M4", "M4/M5", "M5/M6", "M6/M7", "M7/M8")
for (i in 1:length(r80_nM$polymorphic_r80_loci)) {
  new_pol_loci[i] <- r80_nM$polymorphic_r80_loci[i+1] - r80_nM$polymorphic_r80_loci[i]
}


#plot the data
new_poly_loci <- data.frame(param = difference_per_M, new_polymorphic_loci = new_pol_loci[1:8])#, log_mod_transf = log_mod_transformation[1:8])
ggplot(data = new_poly_loci, aes(x = param, y = new_polymorphic_loci)) + 
  theme_bw() +
  geom_point(size = 3, colour = "midnightblue") + 
  geom_line(aes(group = 1), colour = "midnightblue") + 
  scale_y_continuous(trans =  modulus_trans(0.25), breaks = c(-1000, -100, -10, -1, 1, 10, 100, 1000, 10000) ) +  
  ylab("Number of new polymorphic loci") + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 14),
        legend.position = "none", panel.grid.minor.y = element_blank())


#ggsave("new_polymorphic_loci_nM.pdf", width = 10, height = 7, units = "in")




