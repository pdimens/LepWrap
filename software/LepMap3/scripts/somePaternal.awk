#flips informative markers based on a map
#awk -f somePaternal map.txt data.call
BEGIN{
	FS="\t"
	OFS="\t"
}

(NR==FNR && /^[^#]/) {
	map[++markers] = $1
}

(NR!=FNR && FNR<=7){print}

(NR!=FNR && FNR==2){for (i = 3; i<=NF; ++i) {if (!($i in d)) d[$i] = ++p; f[i]=d[$i]}}
(NR!=FNR && FNR==4){for (i = 3; i<=NF; ++i) if ($i==0) pa[f[i], ++count[f[i]]]=i}

(NR!=FNR && FNR>7){
	for (j = 1; j <= p; ++j) {
		if (map[FNR-7]>0) {
			tmp = $(pa[j, 1])
			$(pa[j, 1]) = $(pa[j, 2])
			$(pa[j, 2]) = tmp
		}
 	}
	print
}
