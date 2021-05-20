/#found orientation=/{
	printf $2
	i=3
	while ($i == "+" || $i=="-" ) {
		printf " " $i
		++i
	}
}
