#!/bin/bash

# Compute enrichments in GWAS hits and their vicinity (1Kb), per cell type.

# Define all environment variables

export WD='/users/rg/dgarrido/run-blueprint/run/results'
export GWAS_BED='/users/rg/dgarrido/enrichments/blueprint/files/gwas.bed'

cd $WD

# Loop over cell types

for i in $( ls -d */ ); do
	echo $i
	cd $WD/$i/flux/enrichments/exons # Or flux

	SQTLS_IN="$(bedtools intersect -a sqtls.5fdr.bed -b $GWAS_BED -c | awk '$5>0 {print $4"\t"$5}' | wc -l)"
	NON_SQTLS_IN="$(bedtools intersect -a non-sqtls.5fdr.bed -b $GWAS_BED -c | awk '$5>0 {print $4"\t"$5}' | wc -l)"

	SQTLS_NB="$(cat sqtls.5fdr.bed | wc -l)"
	NON_SQTLS_NB="$(cat non-sqtls.5fdr.bed | wc -l)"

	SQTLS_PERC_IN=$(echo "100*$SQTLS_IN/$SQTLS_NB" | bc -l)
	NON_SQTLS_PERC_IN=$(echo "100*$NON_SQTLS_IN/$NON_SQTLS_NB" | bc -l)

	RATIO5=$(echo "$SQTLS_PERC_IN/$NON_SQTLS_PERC_IN" | bc -l)	
	echo $RATIO5 # Fold enrichment in GWAS-1Kb (%sQTLs/%non-sQTLs) at 5% FDR

	SQTLS_IN="$(bedtools intersect -a sqtls.1fdr.bed -b $GWAS_BED -c | awk '$5>0 {print $4"\t"$5}' | wc -l)"
        NON_SQTLS_IN="$(bedtools intersect -a non-sqtls.1fdr.bed -b $GWAS_BED -c | awk '$5>0 {print $4"\t"$5}' | wc -l)"

        SQTLS_NB="$(cat sqtls.1fdr.bed | wc -l)"
        NON_SQTLS_NB="$(cat non-sqtls.1fdr.bed | wc -l)"

        SQTLS_PERC_IN=$(echo "100*$SQTLS_IN/$SQTLS_NB" | bc -l)
        NON_SQTLS_PERC_IN=$(echo "100*$NON_SQTLS_IN/$NON_SQTLS_NB" | bc -l)

        RATIO1=$(echo "$SQTLS_PERC_IN/$NON_SQTLS_PERC_IN" | bc -l)
        echo $RATIO1 # Fold enrichment in GWAS-1Kb (%sQTLs/%non-sQTLs) at 1% FDR


