#!/bin/bash

# Per cell type exon closeness computation

# Define all environment variables
export WD='/users/rg/dgarrido/run-blueprint/run/results'
export EXONS_PROTCOD='/users/rg/dgarrido/enrichments/blueprint/files/protcod.exons.bed'

cd $WD

# Loop over tissues

for i in $( ls -d */ ); do
	echo $i
	cd $WD/$i/cufflinks/enrichments # or flux
	mkdir closeness
	cd closeness
	
	bedtools closest -a ../exons/sqtls.5fdr.bed -b $EXONS_PROTCOD -t first -d | cut -f4,9 > closeness.5fdr.txt
	bedtools closest -a ../exons/non-sqtls.5fdr.bed -b $EXONS_PROTCOD -t first -d | cut -f4,9 > closeness.control.5fdr.txt
	
	Rscript /users/rg/dgarrido/enrichments/blueprint/scripts/closeness.R closeness.5fdr.txt closeness.control.5fdr.txt closeness.5fdr.png


        bedtools closest -a ../exons/sqtls.1fdr.bed -b $EXONS_PROTCOD -t first -d | cut -f4,9 > closeness.1fdr.txt
        bedtools closest -a ../exons/non-sqtls.1fdr.bed -b $EXONS_PROTCOD -t first -d | cut -f4,9 > closeness.control.1fdr.txt

        Rscript /users/rg/dgarrido/enrichments/blueprint/scripts/closeness.R closeness.1fdr.txt closeness.control.1fdr.txt closeness.1fdr.png


done
