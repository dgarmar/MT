#!/bin/bash

# Compute enrichment of sQTLs (pooled) in splice sites.

# Define all environment variables

export WD='/users/rg/dgarrido/run-blueprint/run/results'
export SSITES='/users/rg/dgarrido/enrichments/blueprint/files/splice.sites.bed'

cd $WD/globalenrich

# Re-using the previous so only cufflinks and total RNA

SQTLS_IN="$(bedtools intersect -a pooled.sqtls.5fdr.bed -b $SSITES -c | awk '$5>0 {print $4"\t"$5}' | wc -l)"
NON_SQTLS_IN="$(bedtools intersect -a pooled.non-sqtls.5fdr.bed -b $SSITES -c | awk '$5>0 {print $4"\t"$5}' | wc -l)"

SQTLS_NB="$(cat pooled.sqtls.5fdr.bed | wc -l)"
NON_SQTLS_NB="$(cat pooled.non-sqtls.5fdr.bed | wc -l)"

SQTLS_PERC_IN=$(echo "100*$SQTLS_IN/$SQTLS_NB" | bc -l)
NON_SQTLS_PERC_IN=$(echo "100*$NON_SQTLS_IN/$NON_SQTLS_NB" | bc -l)

RATIO5=$(echo "$SQTLS_PERC_IN/$NON_SQTLS_PERC_IN" | bc -l)
echo $RATIO5 # Fold enrichment in GWAS-1Kb (%sQTLs/%non-sQTLs) at 5% FDR

SQTLS_IN="$(bedtools intersect -a pooled.sqtls.1fdr.bed -b $SSITES -c | awk '$5>0 {print $4"\t"$5}' | wc -l)"
NON_SQTLS_IN="$(bedtools intersect -a pooled.non-sqtls.1fdr.bed -b $SSITES -c | awk '$5>0 {print $4"\t"$5}' | wc -l)"

SQTLS_NB="$(cat pooled.sqtls.1fdr.bed | wc -l)"
NON_SQTLS_NB="$(cat pooled.non-sqtls.1fdr.bed | wc -l)"

SQTLS_PERC_IN=$(echo "100*$SQTLS_IN/$SQTLS_NB" | bc -l)
NON_SQTLS_PERC_IN=$(echo "100*$NON_SQTLS_IN/$NON_SQTLS_NB" | bc -l)

RATIO1=$(echo "$SQTLS_PERC_IN/$NON_SQTLS_PERC_IN" | bc -l)
echo $RATIO1 # Fold enrichment in GWAS-1Kb (%sQTLs/%non-sQTLs) at 1% FDR

