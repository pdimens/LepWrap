#converts a genotype file (SNP_name + AA/AB/BB for each individual) to posterior
#awk [-verror=0.001] -f affx2post.awk genotypes.txt >genotypes.post
function alleles2Index(a, b)
{
	if (b == 1)
		return a
	if (b == 2)
		return 3 + a
	if (b == 3)
		return 5 + a
	return 10
}
function allele1(code)
{
	if (code <= 4)
		return 1
	if (code <= 7)
		return 2
	if (code <= 9)
		return 3
	return 4
}
function allele2(code)
{
	if (code <= 4)
		return code
	if (code <= 7)
		return code - 3
	if (code <= 9)
		return code - 5
	return 4
}

function distance(code1, code2       ,a1,a2,b1,b2)
{
	if (code1 == code2)
		return 0
	a1 = allele1(code1)
	a2 = allele2(code1)

	b1 = allele1(code2)
	b2 = allele2(code2)

	if (a1 == b1 || a1 == b2 || a2 == b1 || a2 == b2)
		return 1
	return 2
}



BEGIN{
	FS="\t"
	OFS="\t"
	if (error == "")
		error = 0.001

	code = 1
	for (i = 1; i <= 4; ++i)
		for (j = i; j <= 4; ++j) {
			s = ""
			for (k = 1; k <= 10; ++k)
				s = s " " error ^ distance(code, k)
			map[i " " j] = substr(s, 2)
			map[j " " i] = map[i " " j]
			++code
		}

	map["0 0"] = "1 1 1 1 1 1 1 1 1 1"
	map2["AA"]="1 1"
	map2["AB"]="1 2"
	map2["BB"]="2 2"
	map2["NoCall"]="0 0"
}

($1 !~ /^#/) {
	$1=$1"\t"$1
	if (++line > 1)
		for (i = 2; i <= NF; ++i)
			$i = map[map2[$i]]
	print
} 


