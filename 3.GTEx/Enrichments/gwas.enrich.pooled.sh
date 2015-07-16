#!/bin/bash

# All sQTLs pooled - enrichments in GWAS hits (1Kb)

# Define all environment variables

export WD='/users/rg/dgarrido/run-GTEx1/results'
export GWAS_BED='/users/rg/dgarrido/enrichments/gtex/all/gwas.bed'

cd $WD/globalenrich

SQTLS_IN="$(bedtools intersect -a pooled.sqtls.bed -b $GWAS_BED -c | awk '$5>0 {print $4"\t"$5}' | wc -l)"
NON_SQTLS_IN="$(bedtools intersect -a pooled.non-sqtls.bed -b $GWAS_BED -c | awk '$5>0 {print $4"\t"$5}' | wc -l)"

SQTLS_NB="$(cat pooled.sqtls.bed | wc -l)"
NON_SQTLS_NB="$(cat pooled.non-sqtls.bed | wc -l)"

SQTLS_PERC_IN=$(echo "100*$SQTLS_IN/$SQTLS_NB" | bc -l)
NON_SQTLS_PERC_IN=$(echo "100*$NON_SQTLS_IN/$NON_SQTLS_NB" | bc -l)

RATIO5=$(echo "$SQTLS_PERC_IN/$NON_SQTLS_PERC_IN" | bc -l)
echo $RATIO5 # Enrichment in exons of sQTLs at 5% FDR

SQTLS_IN="$(bedtools intersect -a pooled.sqtls.1fdr.bed -b $GWAS_BED -c | awk '$5>0 {print $4"\t"$5}' | wc -l)"
NON_SQTLS_IN="$(bedtools intersect -a pooled.non-sqtls.1fdr.bed -b $GWAS_BED -c | awk '$5>0 {print $4"\t"$5}' | wc -l)"

SQTLS_NB="$(cat pooled.sqtls.1fdr.bed | wc -l)"
NON_SQTLS_NB="$(cat pooled.non-sqtls.1fdr.bed | wc -l)"

SQTLS_PERC_IN=$(echo "100*$SQTLS_IN/$SQTLS_NB" | bc -l)
NON_SQTLS_PERC_IN=$(echo "100*$NON_SQTLS_IN/$NON_SQTLS_NB" | bc -l)

RATIO1=$(echo "$SQTLS_PERC_IN/$NON_SQTLS_PERC_IN" | bc -l)
echo $RATIO1 # Enrichment in exons of sQTLs at 1% FDR

