
#!/bin/bash

## Iterate over all tissues and plot results (sQTls-sGenes)

export WD='/users/rg/dgarrido/run-GTEx1/results'

cd $WD

## Loop over tissues

for i in $( ls -d * ); do
	echo $i
	Rscript /users/rg/dgarrido/autoplot/autoplot.R $i
done	
