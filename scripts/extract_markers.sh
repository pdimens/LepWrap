#! /usr/bin/env bash

file=$(echo $1)
echo "Extracting the marker names and positions from the LepMap3 data file ($file)"

zcat $1 | awk '(NR>=7){print $1"\t"$2}' > snps.txt

