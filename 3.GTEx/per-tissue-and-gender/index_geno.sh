
## Submission script for index_geno.R

#!/bin/bash

#$ -N index_geno  
#$ -cwd
#$ -M diego.garrido@crg.eu -m abe 
#$ -q rg-el6 
#$ -o ../clust_out/$JOB_NAME-$JOB_ID.out 
#$ -e ../clust_err/$JOB_NAME-$JOB_ID.err
#$ -l virtual_free=20G 
#$ -l h_rt=03:00:00

GENO="../input/snps.tsv"
Rscript index_geno.R $GENO
