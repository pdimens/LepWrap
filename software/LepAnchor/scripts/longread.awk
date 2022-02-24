#awk -f longread.awk aln.paf
#aln.paf: minimap2
#process long read alignment data for Lep-Anchor
BEGIN{
	#mapping quality limit
	if (qLimit == "")
		qLimit=1
	FS="\t"
	rev["+"] = "-"
	rev["-"] = "+"
}
(NF>=12 && $12>=qLimit){
	if ($1 == prevR) {
		for (i = 1; i <= n; ++i) { 
			processOnePair($6,$8,$9,$3,$4,$5,c[i],start[i],end[i],rstart[i],rend[i],o[i],min(q[i], $12))
		}
		++n
	} else
		n = 1
	#print within contig part
	print $6"\t"$8"\t"$6"\t"$9"\t++\t"$12
	print $6"\t"$9"\t"$6"\t"$8"\t--\t"$12
	o[n]=$5 #orientation 
	c[n]=$6 #contig

	start[n]=$8 #start
	end[n]=$9 #end

	rstart[n]=$3 #read start
	rend[n]=$4 #read end
	
	q[n]=$12 #quality
	
	prevR = $1

}
function min(a,b){
	if (a <= b)
		return a
	else
		return b
}
function processOnePair(c1, start1, end1, rstart1, rend1, o1, c2, start2, end2, rstart2, rend2, o2, quality   ,lp1,lp2,overlap)
{
	if (rstart1 > rstart2) {
		processOnePair(c2, start2, end2, rstart2, rend2, o2, c1, start1, end1, rstart1, rend1, o1, quality)
		return;
	}
	#now rstart1 <= rstart2
#	print "process"	
	overlap = 0 
	if (rend1 > rstart2)
		overlap = rend1 - rstart2
	overlap = int(overlap * 0.5)
	
	if (o1 == "+")
		lp1 = end1 + 1 - overlap
	else
		lp1 = start1 + overlap

	if (o2 == "+")
		lp2 = start2 + overlap
	else
		lp2 = end2 + 1 - overlap
		
	#print both ways
	print c1"\t"lp1"\t"c2"\t"lp2"\t"o1 o2"\t"quality
	print c2"\t"lp2"\t"c1"\t"lp1"\t"rev[o2] rev[o1]"\t"quality
#	print "end"	
}
END{
}
