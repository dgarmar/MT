
#!/bin/bash

# Prepare files for sQTL sharing estimation

# Define environmental variables

export WD='/users/rg/dgarrido/run-blueprint/run/results'
export PROTCOD='/users/rg/dgarrido/enrichments/blueprint/files/protcod.list'

cd $WD

# Loop over cell types and get pairs sQTL-sGene

for i in $( ls -d */ ); do
	echo $i
	cd $WD/$i/cufflinks # Also flux
	mkdir sharing
	cd sharing	

	grep -F -w -f $PROTCOD ../sQTLs-sig5.tsv | cut -f1-2 | sort | uniq > pairs-sig5.list
	grep -F -w -f $PROTCOD ../sQTLs-sig1.tsv | cut -f1-2 | sort | uniq > pairs-sig1.list
	grep -F -w -f $PROTCOD ../sQTLs-all.tsv | cut -f1-2,10 | sort | uniq > pairs-all.list
	
done
