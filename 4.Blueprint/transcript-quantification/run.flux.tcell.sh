#!/bin/bash

# Submission script. Run Flux-Capacitor on T-cells' BAM files (total RNA)

#$ -N flux-tcell-totalrna
#$ -cwd
#$ -M diego.garrido@crg.eu -m abe
#$ -q short-sl65,long-sl65,rg-el6 
#$ -o clust_out/$JOB_NAME-$JOB_ID.out 
#$ -e clust_err/$JOB_NAME-$JOB_ID.err
#$ -l virtual_free=50G,h_rt=06:00:00
#$ -t 1-40
#$ -pe smp 2 

if [ "$SGE_TASK_ID" -eq "$SGE_TASK_LAST" ]
then
    echo "Job Array: $JOB_ID ($JOB_NAME) ran $SGE_TASK_LAST tasks and ::::::    finished" | /bin/mail -s "Last Task: $SGE_TASK_LAST, current task $SGE_TASK_ID" diego.garrido@crg.es
fi

export WD='/no_backup_isis/rg/projects/Blueprint/WP10/alignments/homo_sapiens/Venous_blood'
export ANNO='/users/rg/dgarrido/run-blueprint/annotation/gencode.v15.annotation.nochr.gtf'

module load flux-capacitor/1.6.1-Java-1.7.0_21

INPUT=$( /bin/find *.bam $WD 2>/tmp/flux.err | /bin/grep '.bam' | /bin/grep -e 'CD4.*total_RNA' | /bin/sort | /bin/sed -n ${SGE_TASK_ID}p )
OUTPUT=("output.names" $( /bin/find *.bam $WD 2>/tmp/flux.err | /bin/grep '.bam' | /bin/grep -e 'CD4.*total_RNA' | /bin/sort | /bin/sed 's/^.*McGill\/\(.*\)\.total.*/\1/' ))

flux-capacitor -i $INPUT -a $ANNO -m PAIRED_STRANDED --read-strand MATE2_SENSE -o ../flux/tcells/total_rna/gtf/${OUTPUT[ $SGE_TASK_ID ]}.gtf


