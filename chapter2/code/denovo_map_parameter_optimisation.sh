#!/bin/bash
# Author: Rohan Raval
# Script: denovo_map_parameter_optimisation.sh
# Desc: performs the denovo_map.pl and iterates over parameters for optimisation
# Arguments: NA
# Date: Nov 2019

M_params=$(seq 1 8)
m_params=$(seq 1 6)
n_params=$(seq 0 10)

#to vary the -M parameter of ustacks
for i in $M_params
	do denovo_map.pl --samples ~/Documents/Data/demultiplexed_sequencing_data/all_samples/ \
	--popmap ~/Documents/Data/demultiplexed_sequencing_data/popmap.txt \
	-o ~/Documents/Data/denovo_map_out/M/M$i/ \
	-T 8 \
	-M $i \
	-m 3 \
	-n 0 \
	--paired \
	-d

done


#to vary the -m parameter of ustacks
for i in $m_params

	do denovo_map.pl --samples ~/Documents/Data/demultiplexed_sequencing_data/all_samples/ \
	--popmap ~/Documents/Data/demultiplexed_sequencing_data/popmap.txt \
	-o ~/Documents/Data/denovo_map_out/m/m$i/ \
	-T 8 \
	-M 2 \
	-m $i \
	-n 0 \
	--paired \
	-d
done


#to vary the -n parameter of cstacks
for i in $n_params

	do denovo_map.pl --samples ~/Documents/Data/demultiplexed_sequencing_data/all_samples/ \
	--popmap ~/Documents/Data/demultiplexed_sequencing_data/popmap.txt \
	-o ~/Documents/Data/denovo_map_out/N/N$i/ \
	-T 8 \
	-M 2 \
	-m 3 \
	-n $i \
	--paired \
	-d

done

#to vary the n parameter of cstacks where n=M
for i in $n_params

	do denovo_map.pl --samples ~/Documents/Data/demultiplexed_sequencing_data/all_samples/ \
	--popmap ~/Documents/Data/demultiplexed_sequencing_data/popmap.txt \
	-o /media/rohan/Maxtor/rohan_parameter_optimisation/denovo_map_out/n_M/N$i/ \
	-T 8 \
	-M $i \
	-m 3 \
	-n $i \
	--paired \
	-r 0.8 \
	-X "populations:--vcf" \
#	-d
done


#to vary the n parameter of cstacks where n=M-1
for i in $n_params

	do denovo_map.pl --samples ~/Documents/Data/demultiplexed_sequencing_data/all_samples/ \
	--popmap ~/Documents/Data/demultiplexed_sequencing_data/popmap.txt \
	-o /media/rohan/Maxtor/rohan_parameter_optimisation/denovo_map_out/n_Mminus1/N$i/ \
	-T 8 \
	-M $(expr $i + 1) \
	-m 3 \
	-n $i \
	--paired \
	-r 0.8 \
	-X "populations:--vcf" \
	-d
done

#to vary the n parameter of cstacks where n=M+1
for i in $n_params

	do denovo_map.pl --samples ~/Documents/Data/demultiplexed_sequencing_data/all_samples/ \
	--popmap ~/Documents/Data/demultiplexed_sequencing_data/popmap.txt \
	-o /media/rohan/Maxtor/rohan_parameter_optimisation/denovo_map_out/n_Mplus1/N$i/ \
	-T 8 \
	-M $(expr $i - 1) \
	-m 3 \
	-n $i \
	--paired \
	-r 0.8 \
	-X "populations:--vcf" \
	-d
done
