#! /usr/bin/env bash

#$1 = target folder 
MAPFLDR=$1

for i in $MAPFLDR/map.*; do
    FNAME=$(basename $i)
    # summarizes the maps, removes leading whitespaces, and sorts by LG
    sed '1,1d' $i | sort | uniq -c | sed 's/^[[:space:]]*//' | sort -k2n > $MAPFLDR/.$FNAME.summary.txt
    # prepend column names
    sed  -i "1i $FNAME LG" $MAPFLDR/.$FNAME.summary.txt
done

# find the map file with the largest number of linkage groups
MAXLG=$(find $MAPFLDR -type f -name ".map.*.summary.txt" -exec grep -H -c '[^[:space:]]' {} \; | sort -nr -t":" -k2 | awk -F: '{print $1; exit;}')
# initialize max number of LG's in summary file
cut -f2 -d " " $MAXLG > $MAPFLDR/maps.summary.txt

# merge summaries into one file
for summfile in $(find ./$MAPFLDR -maxdepth 1 -name ".map.*.summary.txt" | sort -V) ; do
    # append the marker numbers onto the summary file
    cut -f1 -d " " $summfile | paste $MAPFLDR/maps.summary.txt -  > $MAPFLDR/summtemp && mv $MAPFLDR/summtemp $MAPFLDR/maps.summary.txt && rm $summfile
done
# replace all spaces with tabs
sed -i 's/ /\t/g' $MAPFLDR/maps.summary.txt

cat $MAPFLDR/maps.summary.txt | column -t > $MAPFLDR/summary.txt
rm $MAPFLDR/maps.summary.txt && mv $MAPFLDR/summary.txt $MAPFLDR/maps.summary.txt

echo -e "\nExamine the maps produced ("$MAPFLDR/maps.summary.txt") and decide on the best map before proceeding"
echo "if using a screen/tmux environment, detach this session and return to it with the appropriate command when ready"
echo "When you have decided on a map, press Enter to proceed"
read
