#!/bin/bash

input=$1

while IFS= read -r line
do
	#echo "$line"

	if [[ $line == *">"* ]]; then
		peptide="${line:1}"
		
	else
		echo $line
		count=$(grep $line /u/corona/public_html/61bacc01525a0c787b964849d7eb6a27/HVYW7BGXB/star/1_0/Aligned.sortedByCoord.out.bam)
		echo $peptide $count
	fi
done < "$input"