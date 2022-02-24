#awk [-vprefix=CHR -vlg=1 -vgapLength=100] -f makeagp_full.awk chr1.la
BEGIN{
	map["+"]="+"
	map["++"]="+"
	map["+++"]="+"
	map["-"]="-"
	map["--"]="-"
	map["---"]="-"
	map[""]="+"
	map["?"]="+"
	if (lg == "")
		lg = 1
	if (prefix == "")
		prefix = "CHR"
	if (gapLength == "")
		gapLength = 100
}

(/^[^#]/ && /^[^$]/ && $3-$2>=0){
	if (n != "") {
		if (gapLength > 0 && ($9=="null" || ($10!="null" && substr($10, 1, 5) != "chain"))) { #do not add gap for partial haplotypes
			if ($9 != $1 || prevO != $4 || ((($4 ~ /+/) || $4=="?") && $2 != prevE + 1) || (($4 ~ /-/) && $3 != prevS - 1)) { #do not add gap if we have cut the same contig...
				print prefix lg "\t" pos + 1 "\t" pos + gapLength "\t" ++n "\tU\t" gapLength "\tcontig\tno\tna"
				pos += gapLength
			}
		}
	}
	print prefix lg "\t" pos + 1 "\t" pos + ($3-$2+1)  "\t" ++n "\tW\t" $1 "\t" $2 "\t" $3 "\t" map[$4]
	pos += ($3-$2+1)
	prevO = $4
	prevS = $2
	prevE = $3
}

