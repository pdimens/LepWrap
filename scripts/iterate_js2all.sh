#! /usr/bin/env bash

# use this script to iterate JoinSingles2All for a range of LOD limits

if [[ -z "$1" ]]; then
cat <<EOF
Iterate JoinSingles2All over a range of LODlimit values.
[usage]:	iterate_js2all.sh <mapfile> <LOD start> <LOD end> <LOD difference> <informativemask OPTIONAL>
[example]: 	scripts/iterate_js2all.sh 3_SeparateChromosomes/map.31 22 31 5 123
EOF
  exit 1
fi

mkdir -p JoinSingles2All_iter/logs

echo "# The maps in this folder were created based on the map below" > JoinSingles2All_iter/source.map
echo $1 >> JoinSingles2All_iter/source.map

# the map you want to refine
TARGETMAP=$1
# minimum LOD score for iterating
LODMIN=$2
# maximum LOD score for iterating
LODMAX=$3
# LOD difference cutoff
LODDIFF=$4

if [[ -z "$5" ]]; then
  INFMASK=123
else
  INFMASK=$5
fi

for i in $(seq $LODMIN $LODMAX); do
    zcat 2_Filtering/data.filtered.lepmap3.gz | java -cp $CONDA_PREFIX/bin/lepmap3lepmap3 JoinSingles2All map=$TARGETMAP data=- lodLimit=$i lodDifference=4 iterate=1 distortionLod=1 numThreads=10 informativeMask=$INFMASK > JoinSingles2All_iter/logs/map.$i.$4.js2all
    cut -f1 JoinSingles2All_iter/logs/map.$i.$4.js2all > JoinSingles2All_iter/LOD.$i.$4.js2all
done

echo "The generated maps are named LOD.LODlim.LODdiff.js2all"

# generate a summary of the results
scripts/MapSummary.r JoinSingles2All_iter

