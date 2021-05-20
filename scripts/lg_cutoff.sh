#! /usr/bin/env bash

if [[ -z "$1" ]]; then
cat <<EOF

Unassign markers from linkage groups if their LG is greater than a specified value.
Essentially, set an LG maximum, and markers above it will be assigned to group 0 (singles).

[usage]:	lg_cutoff.sh <map.file> <LG max>
# unassign all markers in LG's 25+
[example]: 	scripts/lg_cutoff.sh 3_SeparateChromosomes/map.13 24

EOF
  exit 1
fi

echo "Converting markers from $1 in linkage groups >$2 to singles (0)"

awk -v lim=$2 '{if (NR==1 || $1<=lim) print; else print 0}' $1 > $1.cutoff

echo "Done! See the file $1.cutoff for the results."
echo "Here is a summary of the outcome:"

cat $1.cutoff | tail -n +2 | sort | uniq -c | sort -k2n | tail -n +2