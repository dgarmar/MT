#!/bin/bash

# Compute male and female specificities

# Define variables
export WD='/users/rg/dgarrido/run-GTEx/results'
export PI='/users/rg/dgarrido/sharing/pi1.txt'
export SHARINGR='/users/rg/dgarrido/sharing/sharing.R'

# Loop over tissues and run sharing.R (see 3.GTEx/sharing)
cd $WD

for i in $( ls -d *); do
	echo "$i" 
	cd $WD/$i/female
	echo "female"
	cut -f1,2 sQTLs-sig1-common.tsv | sort | uniq > /tmp/temp_pairs
	pvals="$(grep -F -w -f /tmp/temp_pairs ../male/sQTLs-all-common2.tsv | sed '1d' | cut -f10)"
	pi1="$(Rscript $SHARINGR $pvals)"
	echo "->"$pi1

	cd $WD/$i/male
	echo "male"
	cut -f1,2 sQTLs-sig1-common.tsv | sort | uniq > /tmp/temp_pairs
        pvals="$(grep -F -w -f /tmp/temp_pairs ../female/sQTLs-all-common2.tsv | sed '1d' | cut -f10)"
	pi1="$(Rscript $SHARINGR $pvals)"
	echo "->"$pi1 

done
