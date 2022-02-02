#awk -f cutBed_fixN.awk contigs.length cut.bed
BEGIN{
#	FS="\t"
	OFS="\t"
}
(NR==FNR){
	l[$1]=$2
}
(NR!=FNR){
	if ($3 == "N*") {
		if ($1 in l)
			$3 = l[$1]"*"
		else {
			print "error: length for contig " $1 " not found!"
			exit
		}
	}
	d[++line] = $0
}
END{
	for (i = 1; i <= line; ++i)
		print d[i]
}
