#samtools mpileup -q 10 -Q 10 -s `cat sorted_bams`|awk -f pileupParser2.awk|awk -f pileup2posterior.awk
#needs mapping.txt
#Part of Lep-MAP3
BEGIN{
        "cat mapping.txt"|getline
        if (NF == 0) {
                print "Error: file mapping.txt not found!" > "/dev/stderr"
                exit
        }
	numIndividuals = 0
        for (i = 1; i <= NF; ++i) {
                mapping[i] = $i
		if (!($i in ind))
			++numIndividuals
		ind[$i]
	}
        close("cat mapping.txt")


	if (limit1 == "")
		limit1 = 3 # coverage per individual
	if (limit2 == "")
		limit2 = 0.3 * numIndividuals # number of allowed individuals with lower coverage
	if (limit3 == "")
		limit3 = 0.1 * numIndividuals # minimum allele coverage
	if (limit4 == "")
		limit4 = limit1 * (numIndividuals - limit2)  # minimum total counts

	print "PileupParser2 parameters: limit1=" limit1 " limit2=" limit2 " limit3=" limit3 " limit4="  limit4 > "/dev/stderr"

	FS="\t"
	OFS="\t"
}

function abs(a)
{
	if (a < 0)
		return -a
	else
		return a
}

{
	sum = 0
	for (i = 4; i <= NF; i+=4)
		sum += $i;

	if (sum < limit4)
		next

	delete count
	for (i = 4; i <= NF; i+=4)
		count[mapping[int(i / 4)]]+=$i

	missing = 0
	for (ci in count)
		if (count[ci] < limit1)
			++missing

	if (missing > limit2)
		next

	delete c
	for (i = 5; i <= NF; i+=4) {
		if ($(i-1) == 0)
			$i = ""
		gsub(/\$/,"",$i)  #remove end of reads
		gsub(/\^./,"",$i) #remove quality
		while (match($i, /[+-][1-9][0-9]*/) > 0) { #remove indels
			$i = substr($i, 1, RSTART - 1) substr($i, RSTART + RLENGTH + substr($i, RSTART + 1, RLENGTH - 1))
		}
		tmp = $i
		c[0] += gsub(/[Aa]/, "", tmp)
		c[1] += gsub(/[Cc]/, "", tmp)
		c[2] += gsub(/[Gg]/, "", tmp)
		c[3] += gsub(/[Tt]/, "", tmp)
	}
	alleles = 0
	if (c[0] >= limit3)
		++alleles;
	if (c[1] >= limit3)
		++alleles;
	if (c[2] >= limit3)
		++alleles;
	if (c[3] >= limit3)
		++alleles;

	if (alleles >= 2)
		print
}
