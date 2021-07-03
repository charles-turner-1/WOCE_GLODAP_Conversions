#!/bin/bash

while IFS= read -r line; do    
    outputName=$(printf '%s.txt' "${line}")    
    printf "WOCE Section %s: Creating Expocode file.\n" "${line}"    
	touch $outputName
    printf "Adding header info\n"    
	echo "Year,  GLODAP Expocode,GO-SHIP Expocode" > $outputName
done < $1
