#!/bin/bash

# Tissue by tissue exon closeness values

# Define all environment variables

export WD='/users/rg/dgarrido/run-GTEx1/results'
export PROTCOD='/users/rg/dgarrido/enrichments/gtex/all/protcod.list'
export SNPPOS='/users/rg/dgarrido/enrichments/gtex/all/snps.positions.bed'
export REF_ANNO='/users/rg/dgarrido/sqtlseeker/GTEx/ref_annot/gencode.v19.annotation.gtf'
export EXONS_PROTCOD='/users/rg/dgarrido/enrichments/gtex/all/protcod.exons.bed'

cd $WD

# Loop over tissues

for i in $( ls -d */ ); do
	echo $i
	cd $WD/$i/enrichments
	mkdir closeness
	cd closeness
	
	bedtools closest -a ../exons/sqtls.bed -b $EXONS_PROTCOD -t first -d | cut -f4,9 > closeness.txt
	bedtools closest -a ../exons/non-sqtls.bed -b $EXONS_PROTCOD -t first -d | cut -f4,9 > closeness.control.txt
	
	Rscript /users/rg/dgarrido/enrichments/gtex/all/closeness.R closeness.txt closeness.control.txt closeness.png	

done
