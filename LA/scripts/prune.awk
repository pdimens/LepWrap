#remove blocks of contigs without map support
#awk -f prune ichr1.la >ichr1_pruned.la 2>pruned.la
{
	if (/^#/)
		print
	else {
		if (prevO == "" || $4=="?" || (($4 == "+" || $4 == "-") && (prevO != "+" && prevO != "-")) || (($4 == "++" || $4 == "--") && (prevO != "++" && prevO != "--")) || (($4 == "+++" || $4 == "---") && (prevO != "+++" && prevO != "---")))
			++numBlocks
		data[++line] = $0
		block[line] = numBlocks
		prevO = $4
	}
}
END{
	for (i = 1; i <= line; ++i) {
		$0 = data[i]
		mapinfo[block[i]] += $15
	}
	for (i = 1; i <= line; ++i) {
		$0 = data[i]
		if (mapinfo[block[i]]+0 > 0)
			print 
		else
			print $0 >"/dev/stderr"
	}

}
