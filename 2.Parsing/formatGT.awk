
##
##  Parse VCF genotype files, adapting the content to sQTLseekeR input requirements.
##

BEGIN{OFS="\t";FS="\t"}
{
    l=$1"\t"$2"\t"$3;
    for( i=10;i<=NF;i++ ){
	    if( $1=="#CHROM" ){
	       g=$i
	    } else {
	       g=-1
	       if($i~":0/0:"){g=0}   # Homozygous ref/ref
	       if($i~":1/0:"){g=1}   # Heterozygous 
	       if($i~":0/1:"){g=1}   # Heterozygous
	       if($i~":1/1:"){g=2}   # Homozygous alt/alt
    	}
	  l=l"\t"g
    }
    print l;
}

