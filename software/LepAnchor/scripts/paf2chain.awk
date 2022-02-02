#makes a chain file from minimap2 (-c) paf
#minimap2 --no-long-join -c -PD -g 20000 -r 1000,40000 -k15 -w8 -A3 -B6 -O12 -E1 -s200 -z600 -N50 --min-occ-floor=100 -t 20 ref.fa ref.fa >ref-self.paf
#awk -f chainpaf.awk ref-self.paf|awk -f sortpaf.awk| sort -T . -n -r -k 1,1 | cut -f 2- > ref-self_sorted.paf
#awk -f paf2chain.awk rf-self_sorted.paf | gzip > paf.chain.gz
BEGIN{
	OFS="\t"
	if (scale == "")
		scale = 33
}
($1!=$6){
	score = -1
	cigar = ""
	for (i=13; i<=NF; ++i) {
		if ($i ~ /^AS:i:/)
			score = substr($i, 6)+0
		if ($i ~ /^cg:Z:/)
			cigar = substr($i, 6)
	}
	if ($5 == "+")
		print "chain\t" scale * score, $1, $2, "+", $3, $4, $6, $7, "+" ,$8, $9, NR
	else
		print "chain\t" scale * score, $1, $2, "+", $3, $4, $6, $7, "-" ,$7-$9, $7-$8, NR

	c = cigar

	d = 0
	i = 0
	prev = ""
	do {
		n = c+0
		o = substr(c, length(n)+1, 1)
		c = substr(c, length(n)+2)
#		print n " " o
		if (o == "M") {
			if (prev != "")
				print prev, i, d 
			prev = n
			i = 0
			d = 0
		} else if (o == "I") {
			i+=n
		} else if (o == "D") {
			d+=n
		} else {
			print "error: cannot parse cigar operation " n " " o
			exit
		}
	} while (c != "")
	print prev
#	print cigar
#"\t" cigar
}
