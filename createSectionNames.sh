#!/bin/bash

while IFS= read -r line; do    
	section=$(echo "${line}" | cut -d',' -f1)
    outputName=$(printf '%s.csv' "${section}")    
    echo "WOCE Section $section: Creating Expocode file." 
	touch $outputName
    printf "Adding header info\n"    
	echo "Year,  GLODAP Expocode,GO-SHIP Expocode" > $outputName
done < $1
