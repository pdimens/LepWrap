#zcat sorted_ld.gz|awk -f perm|awk -f maxmatch3
#awk -vOFS="\t" '(!(($1"\t"$2"\t"$3) in d) && !(($3"\t"$4"\t"$1) in d) ){d[$1"\t"$2"\t"$3];d[$3"\t"$4"\t"$1]; print;t=$1;$1=$3;$3=t;t=$2;$2=$4;$4=t;print}'
BEGIN{
	if (window == "")
		window = 10000
}
function round(p)
{
	return 1 + window * int((p-1)/window)
}
{
	r1 = round($2)
	r2 = round($4)
	i1 = $1"\t"$2"\t"$3"\t"r2"\t"r1
	i2 = $3"\t"$4"\t"$1"\t"r1"\t"r2
	if (d[i1] == "" && d[i2] == "") {
		++d[i1]
		++d[i2]
		print $1"\t"r1"\t"$3"\t"r2"\t"$5"\t"$6"\t"$2"\t"$4
		print $3"\t"r2"\t"$1"\t"r1"\t"$5"\t"$6"\t"$4"\t"$2
	}
}
