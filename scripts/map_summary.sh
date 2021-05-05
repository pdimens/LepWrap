#! /usr/bin/env bash

#$1 = target folder 
MAPFLDR=$1

for i in $MAPFLDER/map.*; do
    FNAME=$(basename $i)
    # summarizes the maps, removes leading whitespaces, and sorts by LG
    sed '1,1d' $i | sort | uniq -c | sed 's/^[[:space:]]*//' | sort -k2n > $MAPFLDER/.$FNAME.summary.txt
    # prepend column names
    sed  -i "1i $FNAME LG" $MAPFLDER/.$FNAME.summary.txt
done

# find the map file with the largest number of linkage groups
MAXLG=$(find $MAPFLDER -type f -name ".map.*.summary.txt" -exec grep -H -c '[^[:space:]]' {} \; | sort -nr -t":" -k2 | awk -F: '{print $1; exit;}')
# initialize max number of LG's in summary file
cut -f2 -d " " $MAXLG > $MAPFLDER/maps.summary.txt

# merge summaries into one file
for summfile in $(find ./$MAPFLDER -maxdepth 1 -name ".map.*.summary.txt" | sort -V) ; do
    # append the marker numbers onto the summary file
    cut -f1 -d " " $summfile | paste $MAPFLDER/maps.summary.txt -  > $MAPFLDER/summtemp && mv $MAPFLDER/summtemp $MAPFLDER/maps.summary.txt && rm $summfile
done
# replace all spaces with tabs
sed -i 's/ /\t/g' $MAPFLDER/maps.summary.txt

cat $MAPFLDER/maps.summary.txt | column -t > $MAPFLDER/summary.txt
rm $MAPFLDER/maps.summary.txt && mv $MAPFLDER/summary.txt $MAPFLDER/maps.summary.txt

echo -e "\nExamine the maps produced ("$MAPFLDER/maps.summary.txt") and decide on the best map before proceeding"
echo "if using a screen/tmux environment, detach this session and return to it with the appropriate command when ready"
echo "When you have decided on a map, press Enter to proceed"
read
