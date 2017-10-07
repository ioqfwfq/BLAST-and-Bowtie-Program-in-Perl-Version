my $F = $ARGV[0]; # read $readFile all sam
my $segment;
	open(SEG, '<', "$F") or die "Can't open < $F!";
	while (my $line = <SEG>) {
		if ($line =~ /^@/) {
			$sam_head .= $line;
			next;
		}
		my @temp = split("\t", $line);
		my $direction = '+';
		if ($temp[1] eq '4') {
			next;
		}
		elsif ($temp[1] eq '16' || $temp[1] eq '272') {
			$direction = '-';
		}
		$temp[0] =~ m/(\w+\.\w+)-([0-9]+)/;  # Extract seg number
		my $name = $1;
		my $i = int($2);
		$segment{$name}[$i] .= "$temp[2]:$direction:$temp[3],";
	}
	close SEG;
print $segment;