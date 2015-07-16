#!/bin/bash

# Pooled sQTLs - enrichment in exons

# Define all environment variables

export WD='/users/rg/dgarrido/run-GTEx1/results'
export PROTCOD='/users/rg/dgarrido/enrichments/gtex/all/protcod.list'
export SNPPOS='/users/rg/dgarrido/enrichments/gtex/all/snps.positions.bed'
export REF_ANNO='/users/rg/dgarrido/sqtlseeker/GTEx/ref_annot/gencode.v19.annotation.gtf'
export EXONS_PROTCOD='/users/rg/dgarrido/enrichments/gtex/all/protcod.exons.bed'

cd $WD/globalenrich

cut -f1-2 ../*/sQTLs-sig5.tsv | grep -F -w -f $PROTCOD | cut -f 2 | sort | uniq > pooled.sqtls.list # Pooled sqtls for protein coding genes
cut -f1-2 ../*/sQTLs-sig1.tsv | grep -F -w -f $PROTCOD | cut -f 2 | sort | uniq > pooled.sqtls.1fdr.list # Pooled sqtls for protein coding genes at 1%fdr
cut -f1-2 ../*/sQTLs-all.tsv | grep -F -w -f $PROTCOD | cut -f 2 | sort | uniq > pooled.allsnps.list # Pooled snps for protein coding genes

comm -23 pooled.allsnps.list pooled.sqtls.list > pooled.non-sqtls.list
comm -23 pooled.allsnps.list pooled.sqtls.1fdr.list > pooled.non-sqtls.1fdr.list

grep -F -w -f pooled.sqtls.list $SNPPOS > pooled.sqtls.bed
grep -F -w -f pooled.sqtls.1fdr.list $SNPPOS > pooled.sqtls.1fdr.bed


grep -F -w -f pooled.non-sqtls.list $SNPPOS > pooled.non-sqtls.bed
grep -F -w -f pooled.non-sqtls.1fdr.list $SNPPOS > pooled.non-sqtls.1fdr.bed

SQTLS_IN="$(bedtools intersect -a pooled.sqtls.bed -b $EXONS_PROTCOD -c | awk '$5>0 {print $4"\t"$5}' | wc -l)"
NON_SQTLS_IN="$(bedtools intersect -a pooled.non-sqtls.bed -b $EXONS_PROTCOD -c | awk '$5>0 {print $4"\t"$5}' | wc -l)"

SQTLS_NB="$(cat pooled.sqtls.bed | wc -l)"
NON_SQTLS_NB="$(cat pooled.non-sqtls.bed | wc -l)"

SQTLS_PERC_IN=$(echo "100*$SQTLS_IN/$SQTLS_NB" | bc -l)
NON_SQTLS_PERC_IN=$(echo "100*$NON_SQTLS_IN/$NON_SQTLS_NB" | bc -l)

RATIO5=$(echo "$SQTLS_PERC_IN/$NON_SQTLS_PERC_IN" | bc -l)
echo $RATIO5

SQTLS_IN="$(bedtools intersect -a pooled.sqtls.1fdr.bed -b $EXONS_PROTCOD -c | awk '$5>0 {print $4"\t"$5}' | wc -l)"
NON_SQTLS_IN="$(bedtools intersect -a pooled.non-sqtls.1fdr.bed -b $EXONS_PROTCOD -c | awk '$5>0 {print $4"\t"$5}' | wc -l)"

SQTLS_NB="$(cat pooled.sqtls.1fdr.bed | wc -l)"
NON_SQTLS_NB="$(cat pooled.non-sqtls.1fdr.bed | wc -l)"

SQTLS_PERC_IN=$(echo "100*$SQTLS_IN/$SQTLS_NB" | bc -l)
NON_SQTLS_PERC_IN=$(echo "100*$NON_SQTLS_IN/$NON_SQTLS_NB" | bc -l)

RATIO1=$(echo "$SQTLS_PERC_IN/$NON_SQTLS_PERC_IN" | bc -l)
echo $RATIO1
