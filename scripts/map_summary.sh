#! /usr/bin/env bash

#$1 = LOD end 

for i in maps.splitchrom/map.*; do
    FNAME=$(basename $i)
    # summarizes the maps, removes leading whitespaces, and sorts by LG
    sed '1,1d' $i | sort | uniq -c | sed 's/^[[:space:]]*//' | sort -k2n > maps.splitchrom/.$FNAME.summary.txt
    # prepend column names
    sed  -i "1i $FNAME LG" maps.splitchrom/.$FNAME.summary.txt
done

# find the map file with the largest number of linkage groups
MAXLG=$(find maps.splitchrom -type f -name ".map.*.summary.txt" -exec grep -H -c '[^[:space:]]' {} \; | sort -nr -t":" -k2 | awk -F: '{print $1; exit;}')
# initialize max number of LG's in summary file
cut -f2 -d " " $MAXLG > maps.splitchrom/maps.summary.txt

# merge summaries into one file
for summfile in $(find ./maps.splitchrom -maxdepth 1 -name ".map.*.summary.txt" | sort -V) ; do
    # append the marker numbers onto the summary file
    cut -f1 -d " " $summfile | paste maps.splitchrom/maps.summary.txt -  > maps.splitchrom/summtemp && mv maps.splitchrom/summtemp maps.splitchrom/maps.summary.txt && rm $summfile
done
# replace all spaces with tabs
sed -i 's/ /\t/g' maps.splitchrom/maps.summary.txt

cat maps.splitchrom/maps.summary.txt | column -t > maps.splitchrom/summary.txt
rm maps.splitchrom/maps.summary.txt && mv maps.splitchrom/summary.txt maps.splitchrom/maps.summary.txt

echo -e "\nExamine the maps produced ("maps.splitchrom/maps.summary.txt") and decide on the best map before proceeding"
echo "if using a screen/tmux environment, detach this session and return to it with the appropriate command when ready"
echo "When you have decided on a map, press Enter to proceed"
read