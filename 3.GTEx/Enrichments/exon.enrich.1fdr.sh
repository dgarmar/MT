#!/bin/bash

# Tissue-by-tissue enrichments in exons at 1% FDR

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
	cd $WD/$i

	mkdir enrichments/exons_1fdr
	
	grep -F -w -f $PROTCOD sQTLs-sig1.tsv | cut -f2 | sort | uniq > enrichments/exons_1fdr/sqtls.list # Take all the sqtls in protein coding genes
	grep -F -w -f $PROTCOD sQTLs-all.tsv | cut -f2 | sort | uniq > enrichments/exons_1fdr/allsnps.list # Take all the SNPs

	cd enrichments/exons_1fdr

	comm -23 allsnps.list sqtls.list > non-sqtls.list
	grep -F -w -f sqtls.list $SNPPOS > sqtls.bed
	grep -F -w -f non-sqtls.list $SNPPOS > non-sqtls.bed

	SQTLS_IN="$(bedtools intersect -a sqtls.bed -b $EXONS_PROTCOD -c | awk '$5>0 {print $4"\t"$5}' | wc -l)"
	NON_SQTLS_IN="$(bedtools intersect -a non-sqtls.bed -b $EXONS_PROTCOD -c | awk '$5>0 {print $4"\t"$5}' | wc -l)"

	SQTLS_NB="$(cat sqtls.bed | wc -l)"
	NON_SQTLS_NB="$(cat non-sqtls.bed | wc -l)"

	SQTLS_PERC_IN=$(echo "100*$SQTLS_IN/$SQTLS_NB" | bc -l)
	NON_SQTLS_PERC_IN=$(echo "100*$NON_SQTLS_IN/$NON_SQTLS_NB" | bc -l)

	RATIO=$(echo "$SQTLS_PERC_IN/$NON_SQTLS_PERC_IN" | bc -l)	
	echo $RATIO

done
