#! /usr/bin/env bash

#$1 = LOD end 

for i in maps.splitchrom/map.*; do
    FNAME=$(basename $i)
    # summarizes the maps, removes leading whitespaces, and sorts by LG
    sed '1,1d' $i | sort | uniq -c | sed 's/^[[:space:]]*//' | sort -k2n > maps.splitchrom/.$FNAME.summary.txt
    # prepend column names
    sed  -i "1i $FNAME LG" maps.splitchrom/.$FNAME.summary.txt
done

# initialize max number of LG's in summary file
cut -f2 -d " " maps.splitchrom/.map.$1.summary.txt > maps.splitchrom/maps.summary.txt

# merge summaries into one file
for summfile in $(find ./maps.splitchrom -maxdepth 1 -name ".map.*.summary.txt" | sort -V) ; do
    # append the marker numbers onto the summary file
    cut -f1 -d " " $summfile | paste maps.splitchrom/maps.summary.txt -  > maps.splitchrom/summtemp && mv maps.splitchrom/summtemp maps.splitchrom/maps.summary.txt && rm $summfile
done
# replace all spaces with tabs
sed -i 's/ /\t/g' maps.splitchrom/maps.summary.txt