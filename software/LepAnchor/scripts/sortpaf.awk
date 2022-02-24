{
	for (j = 13; j <= NF; ++j)
		if ($j ~ /^AS:i:/) {
			print (substr($j, 6)+0) "\t" $0
			break
		}
}
