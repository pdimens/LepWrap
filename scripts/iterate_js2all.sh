#! /usr/bin/env bash

# use this script to iterate JoinSingles2All for a range of LOD limits

if [[ -z "$1" ]]; then
cat <<EOF
Iterate JoinSingles2All over a range of LODlimit values.
[usage]:	iterate_js2all.sh <mapfile> <LOD start> <LOD end> <LOD difference>
[example]: 	scripts/iterate_js2all.sh 3_SeparateChromosomes/map.31 22 31 5
EOF
  exit 1
fi

mkdir -p JoinSingles2All_iter

echo "# The maps in this folder were created based on the map below" > JoinSingles2All_iter/chosen.map
echo $1 >> JoinSingles2All_iter/chosen.map

# the map you want to refine
#TARGETMAP=map.31
TARGETMAP=$1
# the linkage group you want to try to split
# minimum LOD score for iterating
LODMIN=$2
# maximum LOD score for iterating
#LODMAX=70
LODMAX=$3
# LOD difference cutoff
LODDIFF=$4

for i in $(seq $LODMIN $LODMAX); do
    zcat 2_Filtering/data_f.call.gz | java -cp LM3 JoinSingles2All map=$TARGETMAP data=- lodLimit=$i lodDifference=4 iterate=1 distortionLod=1 numThreads=10 > JoinSingles2All_iter/map.$i.$4.js2all
done

echo "The generated maps are named map.LODlim.LODdiff.js2all"

# generate a summary of the results
scripts/map_summary.r JoinSingles2All_iter

