#ouputs a final fasta from an agp and a raw (contig) fasta
#awk [-vunUsed=1] -f makefasta.awk file.fasta file_V2.agp >file_V2.fasta

function reverseComplement(s ,i,n,map,ret,nc) {
        n = length(s)
        map["A"]="T"
        map["C"]="G"
        map["G"]="C"
        map["T"]="A"
        map["Y"]="R"
        map["R"]="Y"
        map["S"]="S"
        map["W"]="W"
        map["K"]="M"
        map["M"]="K"
        map["B"]="V"
        map["V"]="B"
        map["D"]="H"
        map["H"]="D"
        map["N"]="N"

        map["a"]="t"
        map["c"]="g"
        map["g"]="c"
        map["t"]="a"
        map["y"]="r"
        map["r"]="y"
        map["s"]="s"
        map["w"]="w"
        map["k"]="m"
        map["m"]="k"
        map["b"]="v"
        map["v"]="b"
        map["d"]="h"
        map["h"]="d"
        map["n"]="n"

        ret = ""
        for (i = n; i >= 1; --i) {
                nc = map[substr(s, i, 1)]
                if (nc == "") {
                        print "Error, letter " substr(s, i, 1) " in fasta"
                        exit
                }
                ret = ret nc
        }
        return ret
}
BEGIN{
}

#store raw-fasta
(NR==FNR){
	if ($1 ~ /^>/) {
		if (str != "") {
			d[name] = str
			str = ""
		}
		name = substr($1, 2)
		names[++numContigs]=name
	} else {
		str = str $0
	}
}
#print fasta according to agp
(NR!=FNR){
	if (str != "") {
		d[name] = str
		str = ""
	}
	name = $1
	if (name != oldName) {
		if (oldName != "")
			print ""
		print ">" name
	}
	oldName = name
	if ($7 == "contig") {
		for (i=1;i<=$6;++i)
			printf "N"
	} else {
		if (!($6 in d)) {
			print "Error, contig " $6 " not found"
			exit
		}
		if ($9 == "+") {
			printf substr(d[$6], $7, $8-$7+1)
		} else {
			printf reverseComplement(substr(d[$6], $7, $8-$7+1))
		}
		used[$6]
	}
}

#print unused contigs
END{
	if (unUsed) {
		print ""
		for (i = 1; i <= numContigs; ++i)
			if (!(names[i] in used)) {
				print ">" names[i]
				print d[names[i]]
			}
	}
}

