#!/user/bin/perl
use warnings;

my ( $samFile, $outFile ) = @ARGV;
sortSAM( $samFile, $outFile );

sub sortSAM {
	my ( $samFile, $outFile ) = @_;
	open( SAM, '<', $samFile )
	  or die "cannot open < $samFile: $!";
	open( OUT, '>', $outFile )
	  or die "cannot open > $outFile: $!";
	my @samRec;
	while (<SAM>) {
		if (/\A@/) {

			#Match and output the header of SAM files
			print OUT $_;
		}
		else {
			chomp;
			push @samRec, [ split("\t") ];
		}
	}
	my @sorted =
	  @samRec[ sort { $samRec[$a][3] <=> $samRec[$b][3] } 0 .. $#samRec ];
	for $A (@sorted) {
		print OUT join( "\t", @$A ), "\n";
	}
	close SAM;
	close OUT;
}