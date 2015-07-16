#!/bin/bash

# Submission script for sqtls.R

#$ -N sqtls  
#$ -cwd
#$ -M diego.garrido@crg.eu -m abe 
#$ -q short-sl65 
#$ -o ../../../clust_out/$JOB_NAME-$JOB_ID.out 
#$ -e ../../../clust_err/$JOB_NAME-$JOB_ID.err
#$ -l virtual_free=20G 
#$ -l h_rt=01:00:00

Rscript ../../../bin/sqtls.R 2>>stderr.txt >>stdout.txt
