#!/bin/bash

# Iterate over tissues to compute sharing (pi1 values) through sharing.R (qvalue R package)

# Define varibles
export WD='/users/rg/dgarrido/run-GTEx1/results'
export PI='/users/rg/dgarrido/sharing/pi1.txt'
export SHARINGR='/users/rg/dgarrido/sharing/sharing.R'

# Loop over tissues. Keep track in the console
cd $WD

for i in $( ls -d *); do
	echo "$i" 
	echo "$i" > /users/rg/dgarrido/sharing/results/$i.$1.txt
	for j in $( ls -d *); do
		if [ "$j" != "$i" ];then
			echo "	$j"
			pvals="$(grep -F -w -f $i/sharing/pairs-sig$1.list $j/sharing/pairs-all.list | cut -f3)"
			pi1="$(Rscript $SHARINGR $pvals)"
			echo "	$pi1" 
			echo "$pi1" >> /users/rg/dgarrido/sharing/results/$i.$1.txt
			#	if [ "$j" == "Adrenal_Gland" ]; then	## For debugging ##
			#		exit 0
			#	fi
		else
			echo "-1" >> /users/rg/dgarrido/sharing/results/$i.$1.txt		
		fi
	done
done
