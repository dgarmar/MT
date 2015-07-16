#!/bin/bash

# Submission script. Run cufflinks on T-cells' BAM files (mRNA)
#$ -N cuff-tcell-mrna
#$ -cwd
#$ -M diego.garrido@crg.eu -m abe
#$ -q short-sl65,long-sl65,rg-el6 
#$ -o clust_out/$JOB_NAME-$JOB_ID.out 
#$ -e clust_err/$JOB_NAME-$JOB_ID.err
#$ -l virtual_free=10G,h_rt=06:00:00
#$ -t 1-40
#$ -pe smp 4

export WD='/no_backup_isis/rg/projects/Blueprint/WP10/alignments/homo_sapiens/Venous_blood'
export ANNO='/users/rg/dgarrido/run-blueprint/annotation/gencode.v15.annotation.nochr.gtf'

module load Cufflinks/2.2.1-goolf-1.4.10-no-OFED

INPUT=$( /bin/find *.bam $WD 2>/tmp/flux.err | /bin/grep '.bam' | /bin/grep -e 'CD4.*mRNA' | /bin/sort | /bin/sed -n ${SGE_TASK_ID}p )
OUTPUT=("output.names" $( /bin/find *.bam $WD 2>/tmp/flux.err | /bin/grep '.bam' | /bin/grep -e 'CD4.*mRNA' | /bin/sort | /bin/sed 's/^.*McGill\/\(.*\)\.polyA_RNA.*/\1/' )) 

cufflinks -q -G $ANNO -u --library-type fr-firststrand --max-bundle-length 1000000000 -p 4 -o ../cufflinks/tcells/mrna/${OUTPUT[ $SGE_TASK_ID ]} $INPUT
