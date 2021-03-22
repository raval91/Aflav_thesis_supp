#!/bin/bash

# version 1.0 | 13-Feb-2020_9-31 | contact: j.bryk@hud.ac.uk

# READ THIS FIRST TO LEARN ABOUT THE DIFFERENT PARAMETERS
# https://www.osc.edu/supercomputing/batch-processing-at-osc/pbs-directives-summary

# ADJUST MEMORY (MAX 376gb FOR jarek NODES OR 126gb FOR marco NODES)
#PBS -l mem=200gb

# ADJUST WALLTIME AFTER WHICH YOUR JOB WILL BE KILLED EVEN IF NOT FINISHED hh:mm:ss
#PBS -l walltime=10:00:00 

# IF YOU WANT CLUSTER TO SEND YOU EMAIL WHEN THE JOB ABORTS, BEGINS AND/OR ENDS, USE THIS 
#PBS -m abe
#PBS -M rohan.raval@hud.ac.uk

# ADJUST NODES AND CORES
# IMPORTANT: MAKE SURE YOUR EXECUTABLE CAN BE RUN ON MULTIPLE CORES OTHERWISE YOU WILL KEEP THE ENTIRE NODE OCCUPIED WHILE ONLY USING A SINGLE CORE
# EACH jarek NODE HAS 20 CORES (2 PROCESSORS @ 10 CORES EACH), EACH marco NODE HAS 48 (2 PROCESSORS @ 24 CORES EACH)

# ONLY nodes=1 ARE ALLOWED
# SET ppn TO BETWEEN 1-20 FOR jarek OR BETWEEN 1-48 for marco (see below)
# SET -q TO jarek FOR BIOLOGY NODES OR marco for CHEMISTRY NODES
# SET -j TO oe TO REPORT OUTPUTS AND ERRORS TO THE SAME FILE job_id.oe
# -t n-m to submit an array of jobs with a unique identifier

#PBS -l nodes=1:ppn=20
#PBS -q jarek
#PBS -V
#PBS -j oe
#PBS -t 9-20%1

# NAME OF YOUR SCRIPT TO RUN
PROG="/home/u1767986/code/nullMAFsimulations.sh"

# SHORT NAME OF YOUR JOB (WILL APPEAR IN qstat)
JOB="MAFsims"

##############################################

# YOUR SCRIPT CALL GOES HERE
# Do not put the entire script content here, just call the script from here and keep the script in a separate file.

bash $PROG


# IMPORTANT:
# /local is the hard drive on the node itself 
# /home is the drive on the head node (there is no /home folder on the node directly)
