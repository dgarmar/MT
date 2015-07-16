
#!/bin/bash

# Iterate over tissues to compute sharing (pi1 values) through sharing.R (qvalue R package)

# Define varibles
export WD='/users/rg/dgarrido/run-blueprint/run/results'
export SHARINGR='/users/rg/dgarrido/sharing/sharing.R'

# Loop over cell types. Keep track in the console
cd $WD

for i in $( ls -d *); do
	echo "$i" 
	echo "$i" > /users/rg/dgarrido/sharing/blueprint/results/$i.$1.txt
	for j in $( ls -d *); do
		if [ "$j" != "$i" ];then
		
			echo "	$j"
			pvals="$(grep -F -w -f $i/cufflinks/sharing/pairs-sig$1.list $j/cufflinks/sharing/pairs-all.list | cut -f3)" # Or flux
			pi1="$(Rscript $SHARINGR $pvals)"
			echo "	$pi1" 
			echo "$pi1" >> /users/rg/dgarrido/sharing/blueprint/results/$i.$1.txt
		
		else
			echo "-1" >> /users/rg/dgarrido/sharing/blueprint/results/$i.$1.txt		
		fi
	
