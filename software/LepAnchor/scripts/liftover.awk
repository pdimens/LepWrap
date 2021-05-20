#awk [-vinverse=1] -f liftover.awk ref.agp chr_pos_file >chr_pos_file.liftover
#liftover genomic coordinates based on an agp file
#if using vinverse=1, coordinates are mapped backwards
#coordinates not found from the agp are kept untouched
BEGIN{
FS="\t"
OFS="\t"
}
(NR==FNR && /^[^#]/) {
	if (inverse) {
		tmp = $1
		$1 = $6
		$6 = tmp
		tmp = $2
		$2 = $7
		$7 = tmp
		tmp = $3
		$3 = $8
		$8 = tmp
	}
	++ni[$6]
	start1[$6, ni[$6]] = $7
	stop1[$6, ni[$6]] = $8
	start2[$6, ni[$6]] = $2
	stop2[$6, ni[$6]] = $3
	chr[$6, ni[$6]] = $1
	orientation[$6, ni[$6]] = $9
}
(NR!=FNR && /^[^#]/) {
	if ($1 in ni) {
		if ($2+0 != 0)
			$2=$2+0 #postion*
		for (i = 1; i <= ni[$1]; ++i)
			if ($2 >= start1[$1, i] && $2 <= stop1[$1, i]) {
				if (orientation[$1, i] == "-") {
					$2 = stop2[$1, i] - ($2 - start1[$1, i])
					$1 = chr[$1, i]
				} else {
					$2 = start2[$1, i] + $2 - start1[$1, i]
					$1 = chr[$1, i]
				}
				print
				next
			}
		print
		warning[$1 "\t" $2]
	} else {
		print
		warning[$1]
	}
}
END {
	for (i in warning)
		print "Warning: scaffold " i " not in agp file2" >"/dev/stderr"
}
