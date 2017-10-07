$F = $ARGV[0];
$c1 = $ARGV[1]-1;
$c2 = $ARGV[2]-1;
open($inF,'<',$F) or die "Can't open < $F!";
<$inF>;
print "gene_id,RPKM\n";
while(<$inF>) {
	my @temp = split("\t");
	print "$temp[$c1],$temp[$c2]\n";
}
close $inF;
