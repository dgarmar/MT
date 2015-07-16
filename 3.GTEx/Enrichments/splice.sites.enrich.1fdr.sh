#!/bin/bash

# Tissue by tissue enrichment of sQTLs in splice sites

# Define all environment variables

export WD='/users/rg/dgarrido/run-GTEx1/results'
export SSITES='/users/rg/dgarrido/enrichments/gtex/all/splice_sites/splice.sites.bed'

cd $WD

# Loop over tissues

for i in $( ls -d * ); do
	echo $i

	SQTLS_IN="$(bedtools intersect -a $i/enrichments/exons_1fdr/sqtls.bed -b $SSITES -c | awk '$5>0 {print $4"\t"$5}' | wc -l)"
	NON_SQTLS_IN="$(bedtools intersect -a $i/enrichments/exons_1fdr/non-sqtls.bed -b $SSITES -c | awk '$5>0 {print $4"\t"$5}' | wc -l)"

	SQTLS_NB="$(cat $i/enrichments/exons_1fdr/sqtls.bed | wc -l)"
	NON_SQTLS_NB="$(cat $i/enrichments/exons_1fdr/non-sqtls.bed | wc -l)"

	SQTLS_PERC_IN=$(echo "100*$SQTLS_IN/$SQTLS_NB" | bc -l)
	NON_SQTLS_PERC_IN=$(echo "100*$NON_SQTLS_IN/$NON_SQTLS_NB" | bc -l)

	RATIO=$(echo "$SQTLS_PERC_IN/$NON_SQTLS_PERC_IN" | bc -l)	
	echo $RATIO

done
