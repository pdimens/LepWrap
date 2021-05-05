#! /usr/bin/env bash

# use this script to refine a particular map in an attempt to split out
# markers from a particular linkage group. The loop iterators over LOD

mkdir -p refine_map

# the map you want to refine
TARGETMAP=map.26
# the linkage group you want to try to split
TARGETLG=1
# minimum LOD score for iterating
LODMIN=22
# maximum LOD score for iterating
LODMAX=70
# remove linkage groups with markers less than this amount
SIZELIM=30

for i in $(seq $LODMIN $LODMAX); do
    zcat data_f.call.gz | java -cp LM3 SeparateChromosomes2 data=- map=maps.splitchrom/$TARGETMAP lg=$TARGETLG sizeLimit=$SIZELIM lodLimit=$i distortionLod=1 numThreads=10 > refine_map/map.$i
done

# generate a summary of the results
scripts/map_summary_refine.sh refine_map
