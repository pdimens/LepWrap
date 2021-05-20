#awk -f phasematch.awk order_reference.txt order_mapped.txt >order_mapped_in_reference_phase.txt
#
BEGIN{
}
(NR==FNR && /^[^#]/){
	for (f = 7; f < NF; f+=3)
		refdata[$1, f] = $f	
	if (numF != "" && numF != NF) {
		print "Error: different number of columns in the input orders" > "/dev/stderr"
		exit 1
	}
	numF = NF
}

(NR!=FNR){
	if (/^[#]/)
		;
	else {
		if (numF != NF) {
			print "Error: different number of columns in the input orders" > "/dev/stderr"
			exit 1
		}
	}
	data[FNR]=$0
}

END{
	for (i = 1; i <= FNR; ++i) {
		$0 = data[i]
		if (/^[#]/)
			;
		else {
			for (f = 7; f < NF; f+=3) {
				if (($1 SUBSEP f) in refdata) {
					ham1[f] += hamming1($f, refdata[$1, f])
					maxham1[f] += maxh
					ham2[f] += hamming2($f, refdata[$1, f])
					maxham2[f] += maxh
				}
			}
		}
		
	}
for (f = 7; f < NF; f+=3) {
		print "***" > "/dev/stderr"
		print "hamming distance1 is " ham1[f] " of " maxham1[f] " (" abs(ham1[f])/(maxham1[f]+0.000000000000000001) ") for family " ++family > "/dev/stderr"
		print "hamming distance2 is " ham2[f] " of " maxham2[f] " (" abs(ham2[f])/(maxham2[f]+0.000000000000000001) ") for family " family > "/dev/stderr"
	}

	for (i = 1; i <= FNR; ++i) {
		$0 = data[i]
		if (/^[#]/)
			print
		else {
			for (f = 7; f < NF; f+=3) {
				n = length($f) / 2
				p1 = substr($f, 1, n)
				p2 = substr($f, n + 1)
				if (ham1[f] < 0) 
					p1 = flip(p1)
				if (ham2[f] < 0) 
					p2 = flip(p2)
				$f = p1 p2
			}
			

			s = $1 "\t" $2 "\t" $3 "\t" $4 " " $5 " " $6
			for (f = 7; f < NF; f+=3)
				s = s "\t" $f " " $(f+1) " " $(f+2)
			print s
		}
		
	}
}
function abs(x) {
	if (x < 0)
		return -x
	return x
}

function flip(x) {
	gsub(/0/, "x", x)
	gsub(/1/, "0", x)
	gsub(/x/, "1", x)
	return x
}

function hamming1(x, y  ,i,xi,yi,ret, n)
{
	n = length(y)
	if (length(x) < n)
		n = length(x)
#	print x " " y
	ret = 0
	maxh = 0
	n = n / 2
	for (i = 1; i <= n; ++i) {
		xi = substr(x, i, 1)
		yi = substr(y, i, 1)
		if (yi != "-" && xi != "-") {
			++maxh
			if (xi == yi) {
				++ret
			}
			else
				--ret
		}
	}
	return ret
}
function hamming2(x, y  ,i,xi,yi,ret, n)
{
	n = length(y)
	if (length(x) < n)
		n = length(x)
#	print x " " y
	ret = 0
	maxh = 0
	for (i = n / 2 + 1; i <= n; ++i) {
		xi = substr(x, i, 1)
		yi = substr(y, i, 1)
		if (yi != "-" && xi != "-") {
			++maxh
			if (xi == yi) {
				++ret
			}
			else
				--ret
		}
	}
	return ret
}

