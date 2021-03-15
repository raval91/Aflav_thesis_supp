THE ECOLOGICAL AND POPULATION GENOMICS OF THE
WILD YELLOW-NECKED MOUSE APODEMUS FLAVICOLLIS


This repository contains the supplementary materials for the above titled thesis.

Included within are the following:

CHAPTER2

code
- denovo_map_parameter_optimisation.sh 
	Optimises parameters to genotype samples using the denovo_map.pl pipeline of Stacks 2.3d.

- stacks_data_extraction.sh
	Extracts the relevant metrics from Stacks 2.3d to optimise parameters. 

- denovo_pipeline.sh
	runs denovo_map.pl on the optimised parameter combination and prepares data for pedigree construction

- M.musculus_insilico_digestion_parallel_SbfI_MseI.R
	Runs an in silico restriction digest of a reference genome and outputs various plots

- M.musculus_insilico_digestion_parallel_SbfI_MseI.R
	Runs an in silico restriction digest of a reference genome and outputs various plots

 tables_and_figs
- A.flavicollis_demographic_data.txt - Life history and sample ID data on samples
	- Sample - sample name
	- ID - ID number
	- Sex - m for male, f for female
	- Age - estimated life stage
	- Trapping_month - month the mouse was trapped
	- Trapping_year - year the mouse was trapped
	- Birth_year - estimated birth year 
	- Dupliated - independently replicated libraries prepared for genotyping error rate estimation. 0 for not duplicated, 1 for duplicated sample 

- coverage_duplicates_merged.txt - mean sequencing coverage per sample after duplicate samples have been merged
	- sample - sample name
	- coverage_type - primary or secondary read
	- mean - mean sequencing coverage
	- stdev - standard deviation of sequencing coverage
	- max - max sequencing coverage
	- n_reads - number of reads per sample

- parameter_optimisation_metrics.txt - metrics extracted from stacks 2.3d to optimise parameters
	- sample - sample name
	- parameter - stacks parameter and parameter value (e.g. M0 means stacks parameter M = 0)
	- n_assembled_loci - number of assembled loci
	- n_polymorphic_loci - number of polymorphic loci
	- n_snps - number of SNPs per sample

- relatedness_results.txt - Pairwise relatedness estimated by max-likelihood, GCTA, KING and PLINK method of moments to identify duplicate samples
	- sampleID_1 - sample name 1
	- sampleID_2 - sample name 2
	- MLE - maximum likelihood estimate of relatedness between sampleID_1 and sampleID_2
	- KING_MoM - KING method of moments estimate of relatedness between sampleID_1 and sampleID_2
	- PLINK_MoM - PLINK method of moments estimate of relatedness between sampleID_1 and sampleID_2
	- GCTA - GCTA maximum likelihood estimate of relatedness between sampleID_1 and sampleID_2
	
- relatedness_results_duplicatesMerged.txt - Pairwise relatedness estimated by max-likelihood estimation after duplicate samples had been merged
	- column names as relatedness_results.txt

- genotyping_error_rates.txt - genotyping error rates as extimated using TIGER
	- Batch - sequencing lane
	- ErrorRate - genotyping error rate

- sample_groups.txt - list of independently sequenced samples across multiple sequencing lanes
	- Sample - sample name
	- Group - unique id for groups of replicate samples
	- Batch - sequencing lane the replicates were sequenced on



CHAPTER3

 code 
- allele_freq_sims/cluster/  
	Scripts to simulate null allele frequency distributions on a high performance cluster

- allele_freq_sims/cluster/MAF_clusterRun_control.sh 
	instructs the cluster to run the simulations as a batch job
	
- allele_freq_sims/cluster/nullMAFsimulations.sh 
	parallelises the simulations
	
- allele_freq_sims/cluster/MAF_clusterRun.R 
 	runs the simulations in R

- allele_freq_sims/cluster/extract_confint.R 
	extracts the confidence intervals from the output of the simulations
	
- allele_freq_sims/nullMAFsimulations.R
	runs the above simulations in parallel on a desktop computer instead of on a cluster

- allele_frequency_change.R
	Analyses and plots allele frequency change over time

- pedigree_construction_R3.6.R
	Pedigree construction and analysis using sequoia in R v3.6

- genetic_contributions.R
	Calculates the genetic and genealogical contributions of founders to a population
	
- genetic_diversity.het
	estimates of genetic diversity (inbreeding coefficients) as calculated by PLINK


 tables_and_figs

- maybe_relatives.txt
	List of possible relationships but excluded from the pedigree

- pedigree.txt
	The pedigree generated by sequoia in R v3.6
	
- genealogical_genetic_contributions.txt
	the genealogical and genetic contributions of founders

- allele_frequencies.txt
	observed allele frequencies between 2015-2017
	
- allele_frequency_change.txt
	change in allele frequencies between 2015-2017

- pvals.txt
	p values of observed vs expected allele frequency change between 2015-2017

CHAPTER4

 code
- cluster_analysis.R
	Performs a K-means cluster analysis on the heterothermy data

- torpor_analysis.R
	Runs a full analysis of heterothermic responses including estimating repeatability, calculating relatedness and estimating heritability

- BMR_analysis_REML_Final
	Runs a full analysis of BMR data from _A. flavicollis_.

 tables_and_figs	
	
- cluster_analysis/cluster_indices.txt
	Table of indices to determine optimal number of clusters in a K-means cluster analysis

- cluster_analysis/k3.pdf - cluster_analysis/k10.pdf
	figures of cluster analyses for heterothermy data when K=3-10

- diagnostic_plots/repeatability/
	Residuals vs fitted values plots of repeatability models including sex+mb+bmr, sex+mb and mb as variables in the models. Also included are diagnostic plots for models of each individual sex
		
- diagnostic_plots/heritability/
	Residuals vs fitted values plots of heritability models including sex+mb+bmr, sex+mb and mb as variables in the models. Also included are diagnostic plots for for models of each individual sex
	


