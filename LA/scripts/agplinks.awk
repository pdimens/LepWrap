#awk -f agplinks.awk contigs.agp >contigs.paf
#contigs.agp will be from scaffolds => contigs direction, columns 1-3 will be contigs, 6-8 scaffolds
#contigs must be in the same order as in scaffolds, i.e. adjacent contigs should be adjacent (often sort -V will ensure this)
($9=="+"||$9=="-"){
	if ($6==prev) {
		if (prevO == "-")
			print "agp" ++n "\t200\t0\t100\t-\t" prevC "\t" (prevE-prevS+1) "\t" (prevS-1) "\t" (prevS+99) "\t60"
		else
			print "agp" ++n "\t200\t0\t100\t+\t" prevC "\t" (prevE-prevS+1) "\t" (prevE-100) "\t" (prevE) "\t60"
		if ($9 == "-")
			print "agp" n "\t200\t100\t200\t-\t" $1 "\t" ($3-$2+1) "\t" ($3-100) "\t" ($3) "\t60"
		else
			print "agp" n "\t200\t100\t200\t+\t" $1 "\t" ($3-$2+1) "\t" ($2-1) "\t" ($2+99) "\t60"
	}
	prev = $6
	prevC = $1
	prevS = $2
	prevE = $3
	prevO = $9
}
