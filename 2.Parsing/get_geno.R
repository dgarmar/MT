
# Obtain the genotype from GTEx in sQTLseekeR required format

## previous sort: sort -n -k1 -k2 snps.tsv > snps_s.tsv
## previous label repetitions: sed..

setwd("/nfs/users/rg/dgarrido/run-GTEx/input/")

# Read genotype file, content and header (previously separated)
snps.samples= read.table(file="geno_parsed.ok.noh.nodup.sorted.vcf", header=FALSE, as.is=TRUE, sep="\t")
snps.header = read.table(file="snps.tsv.header",as.is=TRUE,sep="\t") # Obtained: head -n1 geno_parsed.ok.noh.nodup.sorted.vcf | sed 's/\t/\n/g' > snps.tsv.header CAREFUL WITH HYPHEN/DOTS AFTER

# Remove two first colnames
snps.header = snps.header[-c(1:2),1] 

# Put all colnames, except "end"
colnames(snps.samples) = c("chr","start","snpId",snps.header) 

# Adding the "end" column and colname
# spot = which(names(snps.samples)=="start") = 2L
snps.samples <- data.frame(snps.samples[1:2L],end=snps.samples$start+1,snps.samples[(2L+1):ncol(snps.samples)])

# Change "." by "-" in sample names 
colnames(snps.samples) = chartr(".","-",colnames(snps.samples))

# Leave only the individual name
colnames(snps.samples) = gsub(pattern="(GTEX-[^-]+)-.*", replacement="\\1", perl=TRUE, colnames(snps.samples))

# Sample:
snps.samples[,1:8]

# Save
#colnames(snps.samples)<-chartr("-",".",colnames(snps.samples)) HYPHEN VS DOT PROBLEM
save(snps.samples, file="snps.samples.RData")
write.table(snps.samples, file="snps2.tsv", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
