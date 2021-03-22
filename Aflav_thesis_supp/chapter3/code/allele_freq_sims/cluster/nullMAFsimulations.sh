#!/bin/bash

#### THIS SCRIPT RUNS NULL MAF SIMULATIONS ON THE HUDDERSFIELD UNI CLUSTER ####

#trap to kill parallel backends when loop is complete
trap "exit" INT TERM ERR
trap "kill 0" EXIT

#run the batch of 20 simulations with each loop as a parallel backend
for i in {1..20}
do
R --vanilla --args $i < /home/u1767986/code/MAF_clusterRun.R &
done

wait

#scp output to the home directory
scp /local/u1767986/results/* /home/u1767986/results/

