use strict;

my $readFile = "$ARGV[0].un";

print "Splice file $readFile into segments...";
open(RDF, '<', $readFile) or die "Can't open < $readFile!";
open(SEG, '>', "$readFile.seg") or die "Can't open > $readFile.seg!";
my $Num_Of_Segment = 4;
my $Length_Of_Segment = 25;
while (<RDF>) {
	my @temp = split(" ");
	my $read_name = $temp[0];
	$read_name = substr($read_name, 1);
	my $read = <RDF>;
	chomp($read);
	my $middle_bar = <RDF>;
	my $quality = <RDF>;
	chomp($quality);
	for (my $i = 0; $i < $Num_Of_Segment; $i++) {
		print SEG "\@$read_name-$i\n";
		print SEG substr($read, $i*$Length_Of_Segment, $Length_Of_Segment), "\n";
		print SEG "+\n";
		print SEG substr($quality, $i*$Length_Of_Segment, $Length_Of_Segment), "\n";
	}
}
print "Complete.\n";
close(SEG);
close(RDF);
