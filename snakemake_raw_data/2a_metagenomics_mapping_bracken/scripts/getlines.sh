#!/bin/bash
#Set input file
while IFS=, read -r Path	Sample	Reference_genome	ProviderName	Subject	sample1	plate1	spacer	sample2	sample2
do
	zgrep -c "^+$" "$Path/$ProviderName" + "R1*" | wc -l 
done



