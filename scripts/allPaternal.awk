#flips maternally informative markers to paternally inf
BEGIN{
	FS="\t"
	OFS="\t"
	inf["0 1.0 0 0 0 0 0 0 0 0"]=1
	inf["0 0 1.0 0 0 0 0 0 0 0"]=1
	inf["0 0 0 1.0 0 0 0 0 0 0"]=1
	inf["0 0 0 0 0 1.0 0 0 0 0"]=1
	inf["0 0 0 0 0 0 1.0 0 0 0"]=1
	inf["0 0 0 0 0 0 0 0 1.0 0"]=1

	inf["0 1 0 0 0 0 0 0 0 0"]=1
	inf["0 0 1 0 0 0 0 0 0 0"]=1
	inf["0 0 0 1 0 0 0 0 0 0"]=1
	inf["0 0 0 0 0 1 0 0 0 0"]=1
	inf["0 0 0 0 0 0 1 0 0 0"]=1
	inf["0 0 0 0 0 0 0 0 1 0"]=1

}

(NR<=7){print}

(NR==2){for (i = 3; i<=NF; ++i) {if (!($i in d)) d[$i] = ++p; f[i]=d[$i]}}

(NR==4){for (i = 3; i<=NF; ++i) if ($i==0) pa[f[i], ++count[f[i]]]=i}

(NR==6){for (i = 3; i<=NF; ++i) sex[i]=$i}

#(NR==7){print "#"}

(NR>7){
	for (j = 1; j <= p; ++j) {
		i = 0
		if (inf[$(pa[j, 1])]==1)
			i += sex[pa[j, 1]]
		if (inf[$(pa[j, 2])]==1)
			i += sex[pa[j, 2]]
		if (i == 2) {
			tmp = $(pa[j, 1])
			$(pa[j, 1]) = $(pa[j, 2])
			$(pa[j, 2]) = tmp
		}
 	}
	print
}
