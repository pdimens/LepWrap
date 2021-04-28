#! /usr/bin/env bash

#$1 = LOD end 

for i in refine_map/map.*; do
    FNAME=$(basename $i)
    # summarizes the maps, removes leading whitespaces, and sorts by LG
    sed '1,1d' $i | sort | uniq -c | sed 's/^[[:space:]]*//' | sort -k2n > refine_map/.$FNAME.summary.txt
    # prepend column names
    sed  -i "1i $FNAME LG" refine_map/.$FNAME.summary.txt
done

# find the map file with the largest number of linkage groups
MAXLG=$(find refine_map -type f -name ".map.*.summary.txt" -exec grep -H -c '[^[:space:]]' {} \; | sort -nr -t":" -k2 | awk -F: '{print $1; exit;}')
# initialize max number of LG's in summary file
cut -f2 -d " " $MAXLG > refine_map/maps.summary.txt

# merge summaries into one file
for summfile in $(find ./refine_map -maxdepth 1 -name ".map.*.summary.txt" | sort -V) ; do
    # append the marker numbers onto the summary file
    cut -f1 -d " " $summfile | paste refine_map/maps.summary.txt -  > refine_map/summtemp && mv refine_map/summtemp refine_map/maps.summary.txt && rm $summfile
done
# replace all spaces with tabs
sed -i 's/ /\t/g' refine_map/maps.summary.txt

cat refine_map/maps.summary.txt | column -t > refine_map/summary.txt
rm refine_map/maps.summary.txt && mv refine_map/summary.txt refine_map/maps.summary.txt

echo -e "\nExamine the maps produced ("refine_map/maps.summary.txt") and decide on the best map before proceeding"