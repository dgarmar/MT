#!/bin/bash

# It runs parse.geno.pl

perl parse.gwas.gen.pl genotype.gwas.gen.s > genotype.tsv
awk '{ print $1"\t"$3"\t"$3+1"\t"$2 }' genotype.tsv > start
cut -f 4- genotype.tsv > end
paste start end > genotype.tsv
rm start end 
exit 0
