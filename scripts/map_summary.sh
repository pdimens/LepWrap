#! /usr/bin/env bash

if [[ -z "$1" ]]; then
cat <<EOF

Print a quick summary of a map produced by SeparateChromosomes2.
Use scripts/map_summary.r to summarize many maps at once.

[usage]:	map_summary.sh <map.file>
[example]: 	scripts/map_summary.sh 3_SeparateChromosomes/map.13

EOF
  exit 1
fi

sort $1 | uniq -c | sort -k2n | tail -n +2