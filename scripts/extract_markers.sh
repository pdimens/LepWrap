#! /usr/bin/env bash

if [[ -z "$1" ]]; then
cat <<EOF

Extracts marker information (contig and position) from a Lep-Map3 data file.

[usage]:	extract_markers.sh <data_file.gz>
[example]: 	scripts/extract_markers.sh data_f.call.gz

EOF
  exit 1
fi

file=$(echo $1)
echo "Extracting the marker names and positions from the LepMap3 data file ($file)"

zcat $1 | awk '(NR>=7){print $1"\t"$2}' > snps.txt
