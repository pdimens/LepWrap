#! /usr/bin/env bash

# use this script to refine a particular map in an attempt to split out
# markers from a particular linkage group. The loop iterators over LOD

if [[ -z "$1" ]]; then
cat <<EOF
Attempt to refine a map from SeparateChromosomes2 by splitting out markers from an overclustered linkage group.
Uses 4 threads and distortionLod=1. Requires an input map, which linkage group to modify, LOD start:end values to iterate over, and a minimum size for new clusters.

[usage]:	refinemap.sh <mapfile> <linkage group> <LOD start> <LOD end> <size limit>
[example]: 	scripts/refinemap.sh 3_SeparateChromosomes/map.31 1 22 70 30
EOF
  exit 1
fi

mkdir -p $1.refine

echo "# The maps in this folder were refined based on the map below" > $1.refine/source.map
echo $1 >> $1.refine/source.map

# the map you want to refine
TARGETMAP=$1
# the linkage group you want to try to split
TARGETLG=$2
# minimum LOD score for iterating
LODMIN=$3
# maximum LOD score for iterating
LODMAX=$4
# remove linkage groups with markers less than this amount
SIZELIM=$5

for i in $(seq $LODMIN $LODMAX); do
    zcat 2_Filtering/data.filtered.lepmap3.gz | java -cp $CONDA_PREFIX/bin/lepmap3 SeparateChromosomes2 data=- map=$TARGETMAP lg=$TARGETLG sizeLimit=$SIZELIM lodLimit=$i distortionLod=1 numThreads=4 > $1.refine/map.$i
done

# generate a summary of the results
scripts/MapSummary.r $1.refine
