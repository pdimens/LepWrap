#zcat all.chain.gz|awk [-vminScore=50 -vmaxGap=10000] -f findFullHaplotypes
#find haplotypes with score of at least score of minScore*length and at most maxGap unaligned
BEGIN{
	#use d as delimiter
	d = " "
	if (minScore == "")
		minScore=50
	if (maxGap == "")
		maxGap=10000
}
/chain/{
        numD = 0
        for (dummy = 1; dummy <= 2; ++dummy) { #iterate twice (by swapping alignment)
                if (($2 >= minScore * $4) && ($4 - ($7 - $6) <= maxGap) && !($8 in h) && !($3 in h)){
                        ++numD
                        gap[dummy]= $4 - ($7 - $6) #store gap length
                }

                #swap aligment (sometimes the chain has alignment only one way...)
                d = " "
                if ($5=="+"&&$10=="+") 
                        $0=$1d$2d$8d$9d$10d$11d$12d$3d$4d$5d$6d$7d$13
                if ($5=="+"&&$10=="-") 
                        $0=$1d$2d$8d$9d"+"d$9-$12d$9-$11d$3d$4d"-"d$4-$7d$4-$6d$13
        }
        if (numD == 0) # no haplotype found
                next
        for (dummy = 1; dummy <= 2; ++dummy) { #iterate twice (by swapping alignment)
                if (($2 >= minScore * $4) && ($4 - ($7 - $6) <= maxGap) && !($8 in h) && !($3 in h)){
                        if (numD == 1 || gap[dummy] <= gap[3 - dummy]) { #if either of $3 or $8 is haplotypes, pick one with less gap
		                if ($5=="+"&&$10=="+") 
                                	print $2"\t"$3"\t1\t"$4"\t"$8"\t1\t"$9"\t"$6+1"\t"$7"\t"$10"\t"$11+1"\t"$12
				else
                                	print $2"\t"$3"\t1\t"$4"\t"$8"\t1\t"$9"\t"$6+1"\t"$7"\t"$10"\t"$9-$12+1"\t"$9-$11
                                h[$3]
                        }
                }

                #swap aligment (sometimes the chain has alignment only one way...)
                d = " "
                if ($5=="+"&&$10=="+") 
                        $0=$1d$2d$8d$9d$10d$11d$12d$3d$4d$5d$6d$7d$13
                if ($5=="+"&&$10=="-") 
                        $0=$1d$2d$8d$9d"+"d$9-$12d$9-$11d$3d$4d"-"d$4-$7d$4-$6d$13
        }

}
