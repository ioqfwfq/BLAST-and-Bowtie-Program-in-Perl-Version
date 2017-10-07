my $F = $ARGV[0];
open(SEG, '<', $F) or die "Can't open < $F!";
my %flag;
while (my $line = <SEG>) {
	if ($line =~ /^@/) {
			next;
	}
	my @temp = split(/\t/, $line);
	$flag{$temp[1]} ++;
	#print "$temp[0]";
}
for my $k(sort keys(%flag)) {
	print "flag $k has ", $flag{$k}, " entries.\n";
}
