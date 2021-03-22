#!/bin/bash
# Author: Rohan Raval
# Script: denovo_pipeline.sh
# Desc: Runs denovo pipeline with optimal parameters
# Arguments: NA
# Date: Apr 2020

#run denovo_map with the optimal parameters
#	denovo_map.pl --samples ~/Documents/Data/demultiplexed_sequencing_data/all_samples/ \
#	--popmap ~/Documents/Data/demultiplexed_sequencing_data/popmap.txt \
#	-o ./ \
#	-T 8 \
#	-M 2 \
#	-m 3 \
#	-n 2 \
#	--paired \
#	-X "populations:--vcf"

#maf=(0.05 0.3 0.4)
#r=$(seq 0.5 0.1 0.9)
#run populations to filter for MAF and output a new .vcf file
#for i in ${maf[@]}
#do
#	populations -P ~/Documents/Data/denovo_map_out/m3_M2_n2 \
#	-M ~/Documents/Data/demultiplexed_sequencing_data/popmap.txt \
#	-O ~/Documents/Thesis_analysis/data/populations_out/$i \
#	-r 0.8 \
#	--min_maf $i \
#	--write_single_snp \
#	--plink \
#	-t 8 \
#	--vcf
#	echo "Populations for MAF = $i done" 
#done

##use plink to filter SNPs for HWE, >70% missing data and calculate r2 for all snp pairs
#recode to plink files
#for i in $(ls ~/Documents/Thesis_analysis/data/populations_out)
#do
#~/Documents/Thesis_analysis/./plink \
#--vcf ~/Documents/Thesis_analysis/data/populations_out/$i/pops.snps.vcf \
#--out ~/Documents/Thesis_analysis/data/populations_out/$i/plink_out/recoded \
#--allow-extra-chr \
#--recode
#done

#filter with plink and the recoded files
#for i in $(ls ~/Documents/Thesis_analysis/data/populations_out)
#do
#~/Documents/Thesis_analysis/./plink \
#--file ~/Documents/Thesis_analysis/data/populations_out/$i/plink_out/recoded \
#--out ~/Documents/Thesis_analysis/data/populations_out/$i/plink_out/filtered_snps \
#--allow-extra-chr \
#--allow-no-sex \
#--mind 0.7 \
#--hwe 0.05 'midp' \
#--r2 \ 
#--ld-window-r2 0 \
#--recode vcf
#done

#exclude SNP pairs estimated to be in linkage disequilibrium from the unfiltered data (calculated in R)
#for i in $(ls ~/Documents/Thesis_analysis/data/populations_out)
#do
#~/Documents/Thesis_analysis/./plink \
#--vcf ~/Documents/Thesis_analysis/data/populations_out/$i/pops.snps.vcf \
#--out ~/Documents/Thesis_analysis/data/populations_out/$i/plink_out/filtered_data/r2_pruned_filtered \
#--allow-extra-chr \
#--recode
#done

#for i in $(ls ~/Documents/Thesis_analysis/data/populations_out)
#do
#~/Documents/Thesis_analysis/./plink \
#--file ~/Documents/Thesis_analysis/data/populations_out/$i/plink_out/filtered_data/r2_pruned_filtered \
#--out ~/Documents/Thesis_analysis/data/populations_out/$i/plink_out/filtered_data/r2_pruned_filtered \
#--allow-extra-chr \
#--allow-no-sex \
#--exclude ~/Documents/Thesis_analysis/data/populations_out/$i/plink_out/r2_list.prune.out \
#--mind 0.7 \
#--hwe 0.05 'midp' \
#--recode vcf
#done

#convert the final filtered data into .raw format using plink for sequoia
#for i in $(ls ~/Documents/Thesis_analysis/data/populations_out)
#do
#~/Documents/Thesis_analysis/./plink \
#--file ~/Documents/Thesis_analysis/data/populations_out/$i/plink_out/filtered_data/r2_pruned_filtered \
#--out ~/Documents/Thesis_analysis/data/populations_out/$i/plink_out/filtered_data/Af_filtered4sequoia \
#--allow-extra-chr \
#--recodeA 
#done

#for MAF = 0.45
#~/Documents/Thesis_analysis/./plink \
#--file ~/Documents/Thesis_analysis/data/populations_out/0.4/plink_out/filtered_data/r2_pruned_filtered \
#--out ~/Documents/Thesis_analysis/data/populations_out/0.4/plink_out/filtered_data/r2_pruned_filtered0.45 \
#--allow-extra-chr \
#--allow-no-sex \
#--mind 0.7 \
#--hwe 0.05 'midp' \
#--maf 0.45 
#--recode vcf

#~/Documents/Thesis_analysis/./plink \
#--file ~/Documents/Thesis_analysis/data/populations_out/0.4/plink_out/filtered_data/r2_pruned_filtered \
#--out ~/Documents/Thesis_analysis/data/populations_out/0.4/plink_out/filtered_data/Af_filtered4sequoia_reduced1 \
#--allow-extra-chr \
#--allow-no-sex \
#--extract ~/Documents/Thesis_analysis/data/populations_out/0.4/plink_out/filtered_data/r2_list_reduced.prune1.in \
#--recodeA

##############################################################
###### Pipeline once duplicate samples have been merged ######
##############################################################

#run denovo map for samples after duplicates have been merged
#denovo_map.pl --samples ~/Documents/Data/demultiplexed_sequencing_data/all_samples/ \
#	--popmap ~/Documents/Data/demultiplexed_sequencing_data/popmap_duplicates_merged.txt \
#	-o ./ \
#	-T 8 \
#	-M 2 \
#	-m 3 \
#	-n 2 \
#	--paired \
#	-X "populations:--vcf" \
#	-d

#run populations to filter for MAF (> 0.45) and output a new .vcf file (for generation of pedigree)
populations -P ~/Documents/Data/denovo_map_out/m3_M2_n2_duplicates_merged \
	-M ~/Documents/Data/demultiplexed_sequencing_data/popmap_duplicates_merged.txt \
	-O ~/Documents/Thesis_analysis/data/populations_out/duplicates_merged/0.45 \
	-r 0.8 \
	--min_maf 0.45 \
	--write_single_snp \
	--plink \
	-t 8 \
	--vcf
	echo "Populations for MAF = 0.45 done" 

#run populations to filter for MAF > 0.05 and output a new .vcf	to generate a full panel of snps
#used to correlate pedigree relatedness and genomic relatedness
populations -P ~/Documents/Data/denovo_map_out/m3_M2_n2_duplicates_merged \
	-M ~/Documents/Data/demultiplexed_sequencing_data/popmap_duplicates_merged.txt \
	-O ~/Documents/Thesis_analysis/data/populations_out/duplicates_merged/0.05 \
	-r 0.8 \
	--min_maf 0.05 \
	--write_single_snp \
	--plink \
	-t 8 \
	--vcf
	echo "Populations for MAF = 0.05 done" 

##use plink to filter SNPs for HWE, >70% missing data and calculate r2 for all snp pairs

#recode to plink files
#pops.snps.vcf is an edited version of populations.snps.vcf with the chromosome column edited 
#to "un" and unique IDs in the format "un_##". script in thesis_analysis/code directory

#MAF 0.45
~/Documents/Thesis_analysis/./plink \
--vcf ~/Documents/Thesis_analysis/data/populations_out/duplicates_merged/0.45/pops.snps.vcf \
--out ~/Documents/Thesis_analysis/data/populations_out/duplicates_merged/0.45/plink_out/recoded \
--allow-extra-chr \
--double-id \
--recode

#MAF 0.05
~/Documents/Thesis_analysis/./plink \
--vcf ~/Documents/Thesis_analysis/data/populations_out/duplicates_merged/0.05/pops.snps.vcf \
--out ~/Documents/Thesis_analysis/data/populations_out/duplicates_merged/0.05/plink_out/recoded \
--allow-extra-chr \
--double-id \
--recode

#filter with plink and the recoded files

#MAF 0.45
~/Documents/Thesis_analysis/./plink \
--file ~/Documents/Thesis_analysis/data/populations_out/duplicates_merged/0.45/plink_out/recoded \
--out ~/Documents/Thesis_analysis/data/populations_out/duplicates_merged/0.45/plink_out/filtered_snps \
--allow-extra-chr \
--allow-no-sex \
--mind 0.7 \
--hwe 0.05 'midp' \
--r2 \
--ld-window-r2 0 \
--recode vcf

#MAF 0.05
~/Documents/Thesis_analysis/./plink \
--file ~/Documents/Thesis_analysis/data/populations_out/duplicates_merged/0.05/plink_out/recoded \
--out ~/Documents/Thesis_analysis/data/populations_out/duplicates_merged/0.05/plink_out/filtered_snps \
--allow-extra-chr \
--allow-no-sex \
--mind 0.7 \
--hwe 0.05 'midp' \
--r2 \
--ld-window-r2 0 \
--recode vcf

#exclude SNP pairs estimated to be in linkage disequilibrium from the unfiltered data (calculated in R)
#MAF 0.45
~/Documents/Thesis_analysis/./plink \
--vcf ~/Documents/Thesis_analysis/data/populations_out/duplicates_merged/0.45/pops.snps.vcf \
--out ~/Documents/Thesis_analysis/data/populations_out/duplicates_merged/0.45/plink_out/filtered_data/r2_pruned_filtered \
--allow-extra-chr \
--recode

~/Documents/Thesis_analysis/./plink \
--file ~/Documents/Thesis_analysis/data/populations_out/duplicates_merged/0.45/plink_out/filtered_data/r2_pruned_filtered \
--out ~/Documents/Thesis_analysis/data/populations_out/duplicates_merged/0.45/plink_out/filtered_data/r2_pruned_filtered \
--allow-extra-chr \
--allow-no-sex \
--exclude ~/Documents/Thesis_analysis/data/populations_out/duplicates_merged/0.45/plink_out/r2_list.prune.out \
--mind 0.7 \
--hwe 0.05 'midp' \
--recode vcf

#MAF 0.05
~/Documents/Thesis_analysis/./plink \
--vcf ~/Documents/Thesis_analysis/data/populations_out/duplicates_merged/0.05/pops.snps.vcf \
--out ~/Documents/Thesis_analysis/data/populations_out/duplicates_merged/0.05/plink_out/filtered_data/r2_pruned_filtered \
--allow-extra-chr \
--recode

~/Documents/Thesis_analysis/./plink \
--file ~/Documents/Thesis_analysis/data/populations_out/duplicates_merged/0.05/plink_out/filtered_data/r2_pruned_filtered \
--out ~/Documents/Thesis_analysis/data/populations_out/duplicates_merged/0.05/plink_out/filtered_data/r2_pruned_filtered \
--allow-extra-chr \
--allow-no-sex \
--exclude ~/Documents/Thesis_analysis/data/populations_out/duplicates_merged/0.05/plink_out/r2_list.prune.out \
--mind 0.7 \
--hwe 0.05 'midp' \
--recode vcf

#create another vcf file of only the unique snps in the last 400 snp pairs from calculating r2 (pruned snplist created by r)
#MAF 0.45
~/Documents/Thesis_analysis/./plink \
--file ~/Documents/Thesis_analysis/data/populations_out/duplicates_merged/0.45/plink_out/filtered_data/r2_pruned_filtered \
--out ~/Documents/Thesis_analysis/data/populations_out/duplicates_merged/0.45/plink_out/filtered_data/r2_pruned_filteredReduced \
--allow-extra-chr \
--allow-no-sex \
--extract ~/Documents/Thesis_analysis/data/populations_out/duplicates_merged/0.45/plink_out/filtered_data/r2_list_reduced.prune.in \
--mind 0.7 \
--hwe 0.05 'midp' \
--recode vcf

#MAF 0.05
~/Documents/Thesis_analysis/./plink \
--file ~/Documents/Thesis_analysis/data/populations_out/duplicates_merged/0.05/plink_out/filtered_data/r2_pruned_filtered \
--out ~/Documents/Thesis_analysis/data/populations_out/duplicates_merged/0.05/plink_out/filtered_data/r2_pruned_filteredReduced \
--allow-extra-chr \
--allow-no-sex \
--extract ~/Documents/Thesis_analysis/data/populations_out/duplicates_merged/0.05/plink_out/filtered_data/r2_list_reduced.prune.in \
--mind 0.7 \
--hwe 0.05 'midp' \
--recode vcf
##################################################################################### not done for 0.05 ###################################
##convert the final filtered data into .raw format using plink for sequoia (main and reduced data)
#main
~/Documents/Thesis_analysis/./plink \
--vcf ~/Documents/Thesis_analysis/data/populations_out/duplicates_merged/0.45/plink_out/filtered_data/r2_pruned_filtered.vcf \
--out ~/Documents/Thesis_analysis/data/populations_out/duplicates_merged/0.45/plink_out/filtered_data/Af_filtered4sequoia \
--allow-extra-chr \
--double-id \
--recodeA 
#reduced
~/Documents/Thesis_analysis/./plink \
--vcf ~/Documents/Thesis_analysis/data/populations_out/duplicates_merged/0.45/plink_out/filtered_data/r2_pruned_filteredReduced.vcf \
--out ~/Documents/Thesis_analysis/data/populations_out/duplicates_merged/0.45/plink_out/filtered_data/Af_filtered4sequoia_Reduced \
--allow-extra-chr \
--double-id \
--recodeA 
