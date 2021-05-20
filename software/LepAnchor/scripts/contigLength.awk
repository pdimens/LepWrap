#calculate contig lengths from a fasta
#usage: awk -f contigLength file.fasta >contig.sizes
{
	if (/>/)
		name=substr($1, 2)
	else {
		len[name] += length()
	}
}
END{
	for (i in len)
		print i "\t" len[i]
}
