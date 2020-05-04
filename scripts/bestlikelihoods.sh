#! /usr/bin/env bash

LG=$(find ordermarkers -maxdepth 1 -name "ordered.*.*.txt" | cut -d "." -f2 | sort -V | uniq)
NUMITER=$(find ordermarkers -maxdepth 1 -name "ordered.*.*.txt" | cut -d "." -f3 | sort -V | uniq)
TOTALMAPS=$(find ordermarkers -maxdepth 1 -name "ordered.*.*.txt" | wc -l) 

# pull out best maps for each linkage group
echo "Best ordered maps:"
for i in $(seq 1 $NUMITER $TOTALMAPS); do
    LIKELYMAP=$(sed -n ${i}p ordermarkers/likelihoods.sorted.txt | cut -f1,2 | awk '{print $0, $1 "." $NF}' | cut -d ' ' -f2)
    echo "ordermarkers/$LIKELYMAP.txt" > bestlikelihoods.txt
    #cp ordermarkers/$LIKELYMAP.txt ordermarkers/bestlikelihoods
done
