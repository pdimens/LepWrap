#puts loc files on single line
{
	i = index($0, ";")
	if (i > 0)
		printf(substr($0, 1, i-1))
	else
		print
}
