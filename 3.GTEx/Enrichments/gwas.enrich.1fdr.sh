#!/bin/bash

# Tissue by tissue GWAS (1Kb) enrichments

# Define all environment variables

export WD='/users/rg/dgarrido/run-GTEx1/results'
export GWAS_BED='/users/rg/dgarrido/enrichments/gtex/all/gwas.bed'

cd $WD

# Loop over tissues

for i in $( ls -d */ ); do
	echo $i
	cd $WD/$i/enrichments/gwas

	SQTLS_IN="$(bedtools intersect -a ../exons_1fdr/sqtls.bed -b $GWAS_BED -c | awk '$5>0 {print $4"\t"$5}' | wc -l)"
	NON_SQTLS_IN="$(bedtools intersect -a ../exons_1fdr/non-sqtls.bed -b $GWAS_BED -c | awk '$5>0 {print $4"\t"$5}' | wc -l)"

	SQTLS_NB="$(cat ../exons_1fdr/sqtls.bed | wc -l)"
	NON_SQTLS_NB="$(cat ../exons_1fdr/non-sqtls.bed | wc -l)"

	SQTLS_PERC_IN=$(echo "100*$SQTLS_IN/$SQTLS_NB" | bc -l)
	NON_SQTLS_PERC_IN=$(echo "100*$NON_SQTLS_IN/$NON_SQTLS_NB" | bc -l)

	RATIO=$(echo "$SQTLS_PERC_IN/$NON_SQTLS_PERC_IN" | bc -l)	
	echo $RATIO

done
