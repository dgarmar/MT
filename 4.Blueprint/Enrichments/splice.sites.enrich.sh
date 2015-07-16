#!/bin/bash

# Compute sQTLs enrichment in splice sites, per cell type

# Define all environment variables

export WD='/users/rg/dgarrido/run-blueprint/run/results'
export SSITES='/users/rg/dgarrido/enrichments/blueprint/files/splice.sites.bed'

cd $WD

# Loop over cell types

for i in $( ls -d * ); do
	echo $i
	cd $WD/$i/cufflinks/enrichments/exons # or flux

	# 5%FDR
	SQTLS_IN="$(bedtools intersect -a sqtls.5fdr.bed -b $SSITES -c | awk '$5>0 {print $4"\t"$5}' | wc -l)"
	NON_SQTLS_IN="$(bedtools intersect -a non-sqtls.5fdr.bed -b $SSITES -c | awk '$5>0 {print $4"\t"$5}' | wc -l)"

	SQTLS_NB="$(cat sqtls.5fdr.bed | wc -l)"
	NON_SQTLS_NB="$(cat non-sqtls.5fdr.bed | wc -l)"

	SQTLS_PERC_IN=$(echo "100*$SQTLS_IN/$SQTLS_NB" | bc -l)
	NON_SQTLS_PERC_IN=$(echo "100*$NON_SQTLS_IN/$NON_SQTLS_NB" | bc -l)

	RATIO=$(echo "$SQTLS_PERC_IN/$NON_SQTLS_PERC_IN" | bc -l)	
	echo $RATIO # Fold enrichment in splice sites (%sQTLs/%non-sQTLs) at 5% FDR

	# 1%FDR
        SQTLS_IN="$(bedtools intersect -a sqtls.1fdr.bed -b $SSITES -c | awk '$5>0 {print $4"\t"$5}' | wc -l)"
        NON_SQTLS_IN="$(bedtools intersect -a non-sqtls.1fdr.bed -b $SSITES -c | awk '$5>0 {print $4"\t"$5}' | wc -l)"

        SQTLS_NB="$(cat sqtls.1fdr.bed | wc -l)"
        NON_SQTLS_NB="$(cat non-sqtls.1fdr.bed | wc -l)"

        SQTLS_PERC_IN=$(echo "100*$SQTLS_IN/$SQTLS_NB" | bc -l)
        NON_SQTLS_PERC_IN=$(echo "100*$NON_SQTLS_IN/$NON_SQTLS_NB" | bc -l)

        RATIO=$(echo "$SQTLS_PERC_IN/$NON_SQTLS_PERC_IN" | bc -l)
        echo $RATIO # Fold enrichment in splice sites (%sQTLs/%non-sQTLs) at 1% FDR


done
