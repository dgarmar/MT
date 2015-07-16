#!/bin/bash

# Compute enrichment in exons of sQTLs 

# Define all environment variables

export WD='/users/rg/dgarrido/run-blueprint/run/results'
# cut -f4 /users/rg/dgarrido/run-blueprint/run/input/genes.bed | sort | uniq > /users/rg/dgarrido/enrichments/blueprint/files/protcod.list
export PROTCOD='/users/rg/dgarrido/enrichments/blueprint/files/protcod.list'
# cut -f1-4 /users/rg/dgarrido/run-blueprint/run/input/snps.tsv | sed '1d' > /users/rg/dgarrido/enrichments/blueprint/files/snps.positions.bed
export SNPPOS='/users/rg/dgarrido/enrichments/blueprint/files/snps.positions.bed'
# grep '\sexon\s' /users/rg/dgarrido/run-blueprint/annotation/gencode.v15.annotation.nochr.gtf | grep 'gene_type "protein_coding' | 
# perl -ne 'my @l=split(/\t/,$_);if($l[2] eq "exon"){$l[0]=~s/chr//;$_ =~/gene_id "([^"]+)"/;print $l[0]."\t".$l[3]."\t".$l[4]."\t".$1."\n";}' > /users/rg/dgarrido/enrichments/blueprint/files/protcod.exons.bed 
export EXONS_PROTCOD='/users/rg/dgarrido/enrichments/blueprint/files/protcod.exons.bed'

cd $WD/globalenrich

# ATTENTION Only for total RNA samples, only for cufflinks

cut -f1-2 ../*o*/cufflinks/sQTLs-sig5.tsv | grep -F -w -f $PROTCOD | cut -f 2 | sort | uniq > pooled.sqtls.5fdr.list # Pooled sqtls for protein coding genes at 5% fdr
cut -f1-2 ../*o*/cufflinks/sQTLs-sig1.tsv | grep -F -w -f $PROTCOD | cut -f 2 | sort | uniq > pooled.sqtls.1fdr.list # Pooled sqtls for protein coding genes at 1%fdr
cut -f1-2 ../*o*/cufflinks/sQTLs-all.tsv | grep -F -w -f $PROTCOD | cut -f 2 | sort | uniq > pooled.allsnps.list # Pooled snps for protein coding genes

comm -23 pooled.allsnps.list pooled.sqtls.5fdr.list > pooled.non-sqtls.5fdr.list
comm -23 pooled.allsnps.list pooled.sqtls.1fdr.list > pooled.non-sqtls.1fdr.list

grep -F -w -f pooled.sqtls.5fdr.list $SNPPOS > pooled.sqtls.5fdr.bed
grep -F -w -f pooled.sqtls.1fdr.list $SNPPOS > pooled.sqtls.1fdr.bed


grep -F -w -f pooled.non-sqtls.5fdr.list $SNPPOS > pooled.non-sqtls.5fdr.bed
grep -F -w -f pooled.non-sqtls.1fdr.list $SNPPOS > pooled.non-sqtls.1fdr.bed

SQTLS_IN="$(bedtools intersect -a pooled.sqtls.5fdr.bed -b $EXONS_PROTCOD -c | awk '$5>0 {print $4"\t"$5}' | wc -l)"
NON_SQTLS_IN="$(bedtools intersect -a pooled.non-sqtls.5fdr.bed -b $EXONS_PROTCOD -c | awk '$5>0 {print $4"\t"$5}' | wc -l)"

SQTLS_NB="$(cat pooled.sqtls.5fdr.bed | wc -l)"
NON_SQTLS_NB="$(cat pooled.non-sqtls.5fdr.bed | wc -l)"

SQTLS_PERC_IN=$(echo "100*$SQTLS_IN/$SQTLS_NB" | bc -l)
NON_SQTLS_PERC_IN=$(echo "100*$NON_SQTLS_IN/$NON_SQTLS_NB" | bc -l)

RATIO5=$(echo "$SQTLS_PERC_IN/$NON_SQTLS_PERC_IN" | bc -l)
echo $RATIO5 # Fold enrichment in exons (%sQTLs/%non-sQTLs) at 5% FDR

SQTLS_IN="$(bedtools intersect -a pooled.sqtls.1fdr.bed -b $EXONS_PROTCOD -c | awk '$5>0 {print $4"\t"$5}' | wc -l)"
NON_SQTLS_IN="$(bedtools intersect -a pooled.non-sqtls.1fdr.bed -b $EXONS_PROTCOD -c | awk '$5>0 {print $4"\t"$5}' | wc -l)"

SQTLS_NB="$(cat pooled.sqtls.1fdr.bed | wc -l)"
NON_SQTLS_NB="$(cat pooled.non-sqtls.1fdr.bed | wc -l)"

SQTLS_PERC_IN=$(echo "100*$SQTLS_IN/$SQTLS_NB" | bc -l)
NON_SQTLS_PERC_IN=$(echo "100*$NON_SQTLS_IN/$NON_SQTLS_NB" | bc -l)

RATIO1=$(echo "$SQTLS_PERC_IN/$NON_SQTLS_PERC_IN" | bc -l)
echo $RATIO1 # Fold enrichment in exons (%sQTLs/%non-sQTLs) at 1% FDR


