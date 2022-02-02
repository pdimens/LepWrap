#samtools view hic.bam|awk -f hic.awk
#or samtools view -h hic.bam|samblaster/samblaster -r|awk -f hic.awk
#hic.bam: bwa mem -t16 -5SPM REF R1.fq.gz R2.fq.gz|samtools view -b - >hic.bam
#process hic data for Lep-Anchor
BEGIN{
	#mapping quality limit
	if (qLimit == "")
		qLimit=20
	FS="\t"
}
(NF>=7 && $5>=qLimit){
	if ($1 == prevR) {
		for (i = 1; i <= n; ++i) { #print both ways
			print $3"\t"$4"\t"c[i]"\t"p[i]"\t?\t"min(q[i], $5)
			print c[i]"\t"p[i]"\t"$3"\t"$4"\t?\t"min(q[i], $5)
		}
		++n
	} else
		n = 1
	c[n]=$3
	p[n]=$4
	q[n]=$5
	prevR = $1

}
function min(a,b){
	if (a <= b)
		return a
	else
		return b
}
END{
}
