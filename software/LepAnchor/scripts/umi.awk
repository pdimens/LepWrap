#samtools view -f 0x40 10x.bam|awk -f umi.awk
#10x.bam: bwa mem -t 16 REF R1.fq.gz R2.fq.gz|samtools view -b -|samtools sort - -@4 -T tmp1 -o 10x.bam
#get 10x barcodes from bam aligning raw reads (without cutting barcode off)
#process 10x data for Lep-Anchor
BEGIN{
	#store possible cigar codes for soft clipping of barcode (16bp + 7bp), 16..24 should work
	for (i = 16; i<=24; ++i)
		s[i "S"]
}
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


(and($2, 0x40)!=0){
	if (and($2, 0x10)==0) { #+ orientation
		if (substr($6,1,3) in s)
			print $3"\t"$4"\t"substr($10,1,16) "\t+"
	}
	else if (substr($6,length($6)-2) in s) #-orientation
		print $3"\t"$4"\t"reverseComplement(substr($10,length($10)-15)) "\t-"
}
