#!/user/bin/perl
use warnings;

my ( $refSeq, $samFile, $outFile ) = @ARGV;
pileup( $refSeq, $samFile, $outFile );

sub pileup {
	my ( $refSeq, $samFile, $outFile ) = @_;
	open( REF, '<', $refSeq )
	  or die "cannot open < $refSeq: $!";
	open( SAM, '<', $samFile )
	  or die "cannot open < $samFile: $!";
	open( OUT, '>', $outFile )
	  or die "cannot open > $outFile: $!";
	my ( @sam, %pile, %quali );
	my ( $j, $k, $p, $seq );
	while (<SAM>) {
		if (/\A@/) {
			#print "SAM header!\n";
			next;
		}
		chomp;
		@sam = split("\t");
		$p   = $sam[3];
		$seq = $sam[-2];
		for ( $j = 0 ; $j < length($seq) ; $j++ ) {
			$pile{ $p + $j }  .= substr( $seq,     $j, 1 );
			$quali{ $p + $j } .= substr( $sam[-1], $j, 1 );
		}
	}
	my $curPos = 1;
	my $chr = substr( <REF>, 1 ) =~ /\s/;
	$chr = $`;
	my $ref = '';
	while (<REF>) {
		chomp;
		$ref .= $_;
	}
	for ( $j = 1 ; $j <= length($ref) ; $j++ ) {
		if ( !exists( $pile{$j} ) ) {
			next;
		}
		my $c = substr( $ref, $j - 1, 1 );
		my $pileRead = $pile{$j};
		print OUT $chr, "\t", $j, "\t", $c, "\t", length($pileRead), "\t";
		for ( $k = 0 ; $k < length($pileRead) ; $k++ ) {
			my $cc = substr( $pileRead, $k, 1 );
			if ( uc($cc) eq uc($c) ) {
				print OUT '.';
			}
			else {
				print OUT uc($cc);
			}
		}
		print OUT "\t", $quali{$j}, "\n";
	}
	close REF, SAM, OUT;
}