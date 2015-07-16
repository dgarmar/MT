#!/bin/bash

# Submission script for index_geno.R

#$ -N index.geno  
#$ -cwd
#$ -M diego.garrido@crg.eu -m abe 
#$ -q rg-el6 
#$ -o ../clust_out/$JOB_NAME-$JOB_ID.out 
#$ -e ../clust_err/$JOB_NAME-$JOB_ID.err
#$ -l virtual_free=4G 
#$ -l h_rt=01:00:00

GENO="../input/snps.tsv"
Rscript index_geno.R $GENO
