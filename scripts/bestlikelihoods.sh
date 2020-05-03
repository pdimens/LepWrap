#! /usr/bin/env bash

LG=$(find ordermarkers -maxdepth 1 -name "ordered.*.*.txt" | cut -d "." -f3 | sort -V | uniq)
NUMITER=$(find ordermarkers -maxdepth 1 -name "ordered.*.*.txt" | cut -d "." -f4 | sort -V | uniq)
TOTALMAPS=$(find ordermarkers -maxdepth 1 -name "ordered.*.*.txt" | wc -l) 

for i in $LG ; do
    for j in $(seq 1 $NUMITER) ; do
        LG="ordered.$i"
        ITERUN="$j"
        LIKELIHOOD=$(head -1 ordermarkers/ordered.$i.$j.txt | tail -1 | cut -c 27-) 
        echo -e "$LG\t$ITERUN\t$LIKELIHOOD" >> ordermarkers/likelihoods.txt
    done
done
sort ordermarkers/likelihoods.txt -k1,1V -k3,3nr > ordermarkers/likelihoods.sorted.txt

# pull out best maps for each linkage group
echo "Best ordered maps:"
for i in $(seq 1 $NUMITER $TOTALMAPS); do
    LIKELYMAP=$(sed -n ${i}p ordermarkers/likelihoods.sorted.txt | cut -f1,2 | awk '{print $0, $1 "." $NF}' | cut -d ' ' -f2)
    echo "$LIKELYMAP.txt"
    cp ordermarkers/$LIKELYMAP.txt ordermarkers/bestlikelihoods
done
