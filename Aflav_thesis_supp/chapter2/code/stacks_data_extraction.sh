#!/bin/bash
# Author: Rohan Raval
# Script: parameter_optimisation_data_extraction.sh
# Desc: Data extraction for chosing optimal Stacks parameters
# Arguments: user promted
# Date: Mar 2020

read -p "Provide a results directory: " r
read -p "Insert the full path to the m directory: " m
read -p "Insert the full path to the M directory: " M
read -p "Insert the full path to the n directory: " n
read -p "Insert the full path to the n=M directory: " nM

#create a list of filenames
results_files=$(ls $r/stacks_param_optimisation*.txt)
file_names_m=$(ls $m/m2 | grep -P 'Af_\d+\.alleles\.tsv\.gz' | grep -Po 'Af_\d+')
file_names_M=$(ls $M/M0 | grep -P 'Af_\d+\.alleles\.tsv\.gz' | grep -Po 'Af_\d+')
file_names_n=$(ls $n/N0 | grep -P 'Af_\d+\.alleles\.tsv\.gz' | grep -Po 'Af_\d+')


#extract metrics for m2-m6
cd $m
pwd

for i in {2..6}
do
	for file in $file_names_m
	do
		cd $m/m$i/
		echo "Working on sample: $file"
		echo "number of assembled loci: "
		zcat $file.tags.tsv.gz | grep model | wc -l
		echo "Number of polymorphic loci: "
		zcat $file.tags.tsv.gz | grep E | wc -l
		echo "Number of SNPs: "
		zcat $file.snps.tsv.gz | grep E | wc -l 
	done >> ../../results/stacks_param_optimisation_m$i.txt
cd ../ 
done

#extract metrics for M0-M8
cd $M
pwd



for i in {0..8}
do
	for file in $file_names_M
	do
		cd $M/M$i/
		echo "Working on sample: $file"
		echo "number of assembled loci: "
		zcat $file.tags.tsv.gz | grep model | wc -l
		echo "Number of polymorphic loci: "
		zcat $file.tags.tsv.gz | grep E | wc -l
		echo "Number of SNPs: "
		zcat $file.snps.tsv.gz | grep E | wc -l 
	done >> ../../results/stacks_param_optimisation_M$i.txt
cd ../ 
done



#extract metrics for n0-n10
cd $n
pwd



for i in {0..10}
do
	for file in $file_names_n
	do
		cd $n/N$i/
		echo "Working on sample: $file"
		echo "number of assembled loci: "
		zcat $file.tags.tsv.gz | grep model | wc -l
		echo "Number of polymorphic loci: "
		zcat $file.tags.tsv.gz | grep E | wc -l
		echo "Number of SNPs: "
		zcat $file.snps.tsv.gz | grep E | wc -l 
	done >> ../../results/stacks_param_optimisation_n$i.txt
cd ../ 
done


#pull the data from each file output file to get only a single file of each 
#assembled and polymorphic loci, and number of snps
cd $r
for i in $results_files
do
	param=$(echo $i | grep -Po '[Mmn]\d+')
	grep -Po 'Af_\d+' $i > samples.txt
	grep -E 'assembled' $i -A 1 | grep -P '\d' > assembled_loci.txt
	grep -E 'polymorphic' $i -A 1 | grep -P '\d' > polymorphic_loci.txt
	grep -E 'SNPs' $i -A 1 | grep -P '\d' > SNPs.txt
	n_samples=$(wc -l assembled_loci.txt | grep -Po '\d+')
	yes $param | head -n $n_samples > parameter.txt
	paste samples.txt parameter.txt assembled_loci.txt polymorphic_loci.txt SNPs.txt >> metrics_$param.txt
	rm assembled_loci.txt polymorphic_loci.txt SNPs.txt parameter.txt samples.txt
done
cat metrics*.txt > parameter_optimisation_metrics.txt
sed -i '1 i\sample\tn_assembled_loci\tn_polymorphic_loci\tn_snps' parameter_optimisation_metrics.txt
rm metrics_*

##run populations with r=0.4, 0.6 and 0.8
##for m parameters
for i in $(ls $m)
do
	mkdir $m/$i/populations/
	mkdir $m/$i/populations/p0.4/ $m/$i/populations/p0.6/ $m/$i/populations/p0.8/
	populations -P $i -O $m/$i/populations/p0.4 -t $(nproc --all) -r 40 --vcf --vcf_haplotypes
	populations -P $i -O $m/$i/populations/p0.6 -t $(nproc --all) -r 60 --vcf --vcf_haplotypes
	populations -P $i -O $m/$i/populations/p0.8 -t $(nproc --all) -r 80 --vcf --vcf_haplotypes
	pop_dir=$m/$i/populations/
	stacks-dist-extract $pop_dir/p0.4/populations.log.distribs snps_per_loc_postfilters > $pop_dir/p0.4/locus_distribs_$i.txt
	stacks-dist-extract $pop_dir/p0.6/populations.log.distribs snps_per_loc_postfilters > $pop_dir/p0.6/locus_distribs_$i.txt
	stacks-dist-extract $pop_dir/p0.8/populations.log.distribs snps_per_loc_postfilters > $pop_dir/p0.8/locus_distribs_$i.txt
done

##for M parameters
for i in $(ls $M)
do
	mkdir $M/$i/populations/
	mkdir $M/$i/populations/p0.4/ $M/$i/populations/p0.6/ $M/$i/populations/p0.8/
	populations -P $i -O $M/$i/populations/p0.4 -t $(nproc --all) -r 40 --vcf --vcf_haplotypes
	populations -P $i -O $M/$i/populations/p0.6 -t $(nproc --all) -r 60 --vcf --vcf_haplotypes
	populations -P $i -O $M/$i/populations/p0.8 -t $(nproc --all) -r 80 --vcf --vcf_haplotypes
	pop_dir=$M/$i/populations/
	stacks-dist-extract $pop_dir/p0.4/populations.log.distribs snps_per_loc_postfilters > $pop_dir/p0.4/locus_distribs_$i.txt
	stacks-dist-extract $pop_dir/p0.6/populations.log.distribs snps_per_loc_postfilters > $pop_dir/p0.6/locus_distribs_$i.txt
	stacks-dist-extract $pop_dir/p0.8/populations.log.distribs snps_per_loc_postfilters > $pop_dir/p0.8/locus_distribs_$i.txt 
done

##for n parameters
for i in $(ls $n)
do
	mkdir $n/$i/populations/
	mkdir $n/$i/populations/p0.4/ $n/$i/populations/p0.6/ $n/$i/populations/p0.8/
	populations -P $i -O $n/$i/populations/p0.4 -t $(nproc --all) -r 40 --vcf --vcf_haplotypes
	echo "$i p0.4 done"
	populations -P $i -O $n/$i/populations/p0.6 -t $(nproc --all) -r 60 --vcf --vcf_haplotypes
	echo "$i p0.6 done"
	populations -P $i -O $n/$i/populations/p0.8 -t $(nproc --all) -r 80 --vcf --vcf_haplotypes
	echo "$i p0.8 done"
	pop_dir=$n/$i/populations/
	stacks-dist-extract $pop_dir/p0.4/populations.log.distribs snps_per_loc_postfilters > $pop_dir/p0.4/locus_distribs_$i.txt
	stacks-dist-extract $pop_dir/p0.6/populations.log.distribs snps_per_loc_postfilters > $pop_dir/p0.6/locus_distribs_$i.txt
	stacks-dist-extract $pop_dir/p0.8/populations.log.distribs snps_per_loc_postfilters > $pop_dir/p0.8/locus_distribs_$i.txt 
	echo "$i ALL done"
done


#organise the distribution files for r 0.4, 0.6 and 0.8
for i in $(ls $m)
do
	mv $m/$i/populations/p0.4/locus_distribs_$i.txt $r/locus_distribs_p0.4_$i.txt
	mv $m/$i/populations/p0.6/locus_distribs_$i.txt $r/locus_distribs_p0.6_$i.txt
	mv $m/$i/populations/p0.8/locus_distribs_$i.txt $r/locus_distribs_p0.8_$i.txt
done

for i in $(ls $M)
do
	mv $M/$i/populations/p0.4/locus_distribs_$i.txt $r/locus_distribs_p0.4_$i.txt
	mv $M/$i/populations/p0.6/locus_distribs_$i.txt $r/locus_distribs_p0.6_$i.txt
	mv $M/$i/populations/p0.8/locus_distribs_$i.txt $r/locus_distribs_p0.8_$i.txt
done

for i in $(ls $N)
do
	mv $n/$i/populations/p0.4/locus_distribs_$i.txt $r/locus_distribs_p0.4_$i.txt
	mv $n/$i/populations/p0.6/locus_distribs_$i.txt $r/locus_distribs_p0.6_$i.txt
	mv $n/$i/populations/p0.8/locus_distribs_$i.txt $r/locus_distribs_p0.8_$i.txt
done

##extract the metrics from the populations outputs for r80 loci only
#m
for i in $(ls $m)
do
	echo $i 'r80 assembled loci'
	awk '{if(NR>1)print}' $m/$i/populations/p0.8/populations.haplotypes.tsv | wc -l
	echo $i 'r80 polymorphic loci'
	awk '{if(NR>1)print}' $m/$i/populations/p0.8/populations.haplotypes.tsv | grep -v consensus | wc -l
	echo $i 'r80 snps'
	grep -v "^#" $m/$i/populations/p0.8/populations.snps.vcf | wc -l
done >> $r/populations/m_r80_metrics.txt

#M
for i in $(ls $M)
do
	echo $i 'r80 assembled loci'
	awk '{if(NR>1)print}' $M/$i/populations/p0.8/populations.haplotypes.tsv | wc -l
	echo $i 'r80 polymorphic loci'
	awk '{if(NR>1)print}' $M/$i/populations/p0.8/populations.haplotypes.tsv | grep -v consensus | wc -l
	echo $i 'r80 snps'
	grep -v "^#" $M/$i/populations/p0.8/populations.snps.vcf | wc -l
done >> $r/populations/M_r80_metrics.txt

#n
for i in $(ls $n)
do
	echo $i 'r80 assembled loci'
	awk '{if(NR>1)print}' $n/$i/populations/p0.8/populations.haplotypes.tsv | wc -l
	echo $i 'r80 polymorphic loci'
	awk '{if(NR>1)print}' $n/$i/populations/p0.8/populations.haplotypes.tsv | grep -v consensus | wc -l
	echo $i 'r80 snps'
	grep -v "^#" $n/$i/populations/p0.8/populations.snps.vcf | wc -l
done >> $r/populations/n_r80_metrics.txt

#n=M
for i in $(ls $nM)
do
	echo $i 'r80 assembled loci'
	awk '{if(NR>1)print}' $nM/$i/populations.haplotypes.tsv | wc -l
	echo $i 'r80 polymorphic loci'
	awk '{if(NR>1)print}' $nM/$i/populations.haplotypes.tsv | grep -v consensus | wc -l
	echo $i 'r80 snps'
	grep -v "^#" $nM/$i/populations.snps.vcf | wc -l
done >> $r/populations/nM_r80_metrics.txt

for i in $(ls $r/populations)
do
	param=$(grep -Po '[mMN]\d' $r/populations/$i | head -n 1 | grep -Eo [mMN])
	grep -E 'assembled' $r/populations/$i -A 1 | grep -P '\d\d\d+' > $r/populations/assembled_loci_$param.txt
	grep -E 'polymorphic' $r/populations/$i -A 1 | grep -P '\d\d\d+' > $r/populations/polymorphic_loci_$param.txt
	grep -E 'snps' $r/populations/$i -A 1 | grep -P '\d\d\d+' > $r/populations/SNPs_$param.txt
	#the following was done in excel as it was wasting too much time to figure out why it wasnt working
	#it combines the r80 metrics for each parameter combination into its relevant text file
	#paste $r/populations/assembled_loci_$param.txt $r/populations/polymorphic_loci_$param.txt $r/populations/SNPs_$param.txt >> r80_metrics_$param.txt
	#sed -i '1 i\assembled_r80_loci\tpolymorphic_r80_loci\tr80_snps' r80_metrics_$param.txt
	#rm assembled_loci_$param.txt polymorphic_loci_$param.txt SNPs_$param.txt
done



#for i in $(ls)
#do
#	cp $i/populations.log ../../results/populations/populations_$i.log
#done
