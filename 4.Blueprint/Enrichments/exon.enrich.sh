#!/bin/bash

# Compute per cell type exon enrichments

# Define all environment variables

export WD='/users/rg/dgarrido/run-blueprint/run/results'
# cut -f4 /users/rg/dgarrido/run-blueprint/run/input/genes.bed | sort | uniq > /users/rg/dgarrido/enrichments/blueprint/files/protcod.list
export PROTCOD='/users/rg/dgarrido/enrichments/blueprint/files/protcod.list'
# cut -f1-4 /users/rg/dgarrido/run-blueprint/run/input/snps.tsv | sed '1d' > /users/rg/dgarrido/enrichments/blueprint/files/snps.positions.bed
export SNPPOS='/users/rg/dgarrido/enrichments/blueprint/files/snps.positions.bed'
# grep '\sexon\s' /users/rg/dgarrido/run-blueprint/annotation/gencode.v15.annotation.nochr.gtf | grep 'gene_type "protein_coding' | 
# perl -ne 'my @l=split(/\t/,$_);if($l[2] eq "exon"){$l[0]=~s/chr//;$_ =~/gene_id "([^"]+)"/;print $l[0]."\t".$l[3]."\t".$l[4]."\t".$1."\n";}' > /users/rg/dgarrido/enrichments/blueprint/files/protcod.exons.bed 
export EXONS_PROTCOD='/users/rg/dgarrido/enrichments/blueprint/files/protcod.exons.bed'


cd $WD

# Loop over cell types

for i in $( ls -d */ ); do
	
	echo $i
		
	cd $WD/$i/cufflinks # Also flux
	
	grep -F -w -f $PROTCOD sQTLs-sig5.tsv | cut -f2 | sort | uniq > enrichments/exons/sqtls.5fdr.list # Take all the sqtls in protein coding genes 5% fdr
	grep -F -w -f $PROTCOD sQTLs-sig1.tsv | cut -f2 | sort | uniq > enrichments/exons/sqtls.1fdr.list # Take all the sqtls in protein coding genes 1% fdr
	grep -F -w -f $PROTCOD sQTLs-all.tsv | cut -f2 | sort | uniq > enrichments/exons/allsnps.list # Take all the SNPs

	cd enrichments/exons

	comm -23 allsnps.list sqtls.5fdr.list > non-sqtls.5fdr.list
	grep -F -w -f sqtls.5fdr.list $SNPPOS > sqtls.5fdr.bed
	grep -F -w -f non-sqtls.5fdr.list $SNPPOS > non-sqtls.5fdr.bed

	SQTLS_IN="$(bedtools intersect -a sqtls.5fdr.bed -b $EXONS_PROTCOD -c | awk '$5>0 {print $4"\t"$5}' | wc -l)"
	NON_SQTLS_IN="$(bedtools intersect -a non-sqtls.5fdr.bed -b $EXONS_PROTCOD -c | awk '$5>0 {print $4"\t"$5}' | wc -l)"

	SQTLS_NB="$(cat sqtls.5fdr.bed | wc -l)"
	NON_SQTLS_NB="$(cat non-sqtls.5fdr.bed | wc -l)"

	SQTLS_PERC_IN=$(echo "100*$SQTLS_IN/$SQTLS_NB" | bc -l)
	NON_SQTLS_PERC_IN=$(echo "100*$NON_SQTLS_IN/$NON_SQTLS_NB" | bc -l)

	RATIO5=$(echo "$SQTLS_PERC_IN/$NON_SQTLS_PERC_IN" | bc -l)	
	echo $RATIO5  # Fold enrichment in exons (%sQTLs/%non-sQTLs) at 5% FDR


	comm -23 allsnps.list sqtls.1fdr.list > non-sqtls.1fdr.list
        grep -F -w -f sqtls.1fdr.list $SNPPOS > sqtls.1fdr.bed
	grep -F -w -f non-sqtls.1fdr.list $SNPPOS > non-sqtls.1fdr.bed

	SQTLS_IN="$(bedtools intersect -a sqtls.1fdr.bed -b $EXONS_PROTCOD -c | awk '$5>0 {print $4"\t"$5}' | wc -l)"
	NON_SQTLS_IN="$(bedtools intersect -a non-sqtls.1fdr.bed -b $EXONS_PROTCOD -c | awk '$5>0 {print $4"\t"$5}' | wc -l)"

	SQTLS_NB="$(cat sqtls.1fdr.bed | wc -l)"
	NON_SQTLS_NB="$(cat non-sqtls.1fdr.bed | wc -l)"

  	SQTLS_PERC_IN=$(echo "100*$SQTLS_IN/$SQTLS_NB" | bc -l)
	NON_SQTLS_PERC_IN=$(echo "100*$NON_SQTLS_IN/$NON_SQTLS_NB" | bc -l)

	RATIO1=$(echo "$SQTLS_PERC_IN/$NON_SQTLS_PERC_IN" | bc -l)
  	echo $RATIO1  # Fold enrichment in exons (%sQTLs/%non-sQTLs) at 1% FDR

done
