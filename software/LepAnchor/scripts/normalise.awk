#normalise proximity data based on within bin proximity
#awk -f normalise prox.data prox.data
BEGIN{
	OFS="\t"
}
(NR==FNR && $1==$3 && $2==$4){ #same position
	sum2[$1,$2]+=$5
}
(NR!=FNR) {
	if (!first) {
		for (i in sum2) {
			s = sum2[i]
			if (max == "" || max <= s)
				max = s
			if (min == "" || min >= s)
				min = s
			sum  += s
			sumS += s**2
			++n
		}

		avg = sum / n
		iavg = 100.0 / avg
		sd = sqrt((sumS + n * avg**2 - 2 * avg * sum)/(n - 1))
		print avg "\tsd=\t" sd "\tmin=\t" min "\tmax=\t" max  "\tn=" n > "/dev/stderr"
		first = 1
	}
	s1 = sum2[$1,$2] + 0
	s2 = sum2[$3,$4] + 0
	$5 *= normalise(s1, s2)
	print
}
function maxf(s1, s2)
{
	if (s1 >= s2)
		return s1
	return s2
}
function minf(s1, s2)
{
        if (s1 < s2)
                return s1
        return s2

}
#maxf
function normalise(s1, s2)
{
	if (s1 <= 0.5 * avg && s2 <= 0.5 * avg)
		return 2 * iavg
	else
		return 100.0 / maxf(s1, s2)
}
#minf
function normalise2(s1, s2)
{
	if (s1 <= avg || s2 <= avg)
		return iavg
	else
		return 100.0 / minf(s1, s2)
}

