#!/bin/bash

#Converts Socrates output to a bed file:
#each breakpoint is converted to two bed lines with the single coordinate and support

if [ -z $1 ]; then
	echo "Usage: <socrates results file (paired)>"
	exit
fi

awk -F"\t" '{OFS=":"; print $1, $4, $7, $19}' $1 | awk -F":" 'BEGIN{print "#chr", "pos", "pos", "support"} NR>1{print $1, $2, $2, $5; print $3, $4, $4, $6}' > $1".bed"


echo "Done -- bed file generated for input with .bed extension"
