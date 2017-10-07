#!/user/bin/perl
use warnings;

open( my $inF, '<', "$ARGV[0].F" )
  or die "cannot open < $ARGV[0].F: $!";
open( my $inL, '<', "$ARGV[0].L" )
  or die "cannot open < $ARGV[0].L: $!";
open( my $inTally, '<', "$ARGV[0].tly" )
  or die "cannot open < $ARGV[0].tly: $!";
open( my $inQuery, '<', $ARGV[1] )
  or die "cannot open < $ARGV[1]: $!";
open( my $outFile, '>', $ARGV[2] )
  or die "cannot open > $ARGV[2]: $!";

my $MAXMISMATCH = 1;
my $MAXMULTIMATCH = 0;

#Read F line
my $genome_name = <$inF>;
chomp($genome_name);
my $n           = <$inF>;
chomp($n);
my ( $c, $t );
my ( @charNum, @charset, %charRange );
my $p = 0;
for ( my $i = 0 ; $i < $n ; $i++ ) {
	( $c, $t ) = split( ' ', <$inF> );
	$charNum{$c} = $t;
	$charset[$i] = $c;
	@{ $charRange{$c} } = ( $p, $p + $t - 1 );
	$p += $t;
}
shift @charset;    #Eliminate $ sign
close $inF;

#Read L line
$n = <$inL>;
chomp $n;
my $genome_length = $n-1;
my ( @L, @P );
for ( my $i = 0 ; $i < $n ; $i++ ) {
	( $L[$i], $P[$i] ) = split( ' ', <$inL> );
}
close $inL;

#Read Tally matrix
my %tally;
my @t;
for ( my $i = 0 ; $i < $n ; $i++ ) {
	@t = split( ' ', <$inTally> );
	for ( my $j = 0 ; $j < scalar(@charset) ; $j++ ) {
		$tally{ $charset[$j] }[$i] = $t[$j];
	}
}
close $inTally;

#Read query and do BWA
out2samHead($outFile, $genome_name, $genome_length);
while (<$inQuery>) {
	if (/@(\S+)/) {
		my $seqName = $1;
		chomp($seqName);
		my $seq = <$inQuery>;
		chomp($seq);
		#print $logFile "$seqName, seq = $seq\n";
		my $seqDirection = <$inQuery>;
		chomp($seqDirection);
		my $seqQuality = <$inQuery>;
		chomp($seqQuality);

		my ( $flag, $pos, $mapq, $cigar, $rnext, $pnext, $tlen ) = bwa($seq);
		foreach my $p (@$pos) {
			#print "$p\n";
			out2sam(
				$outFile, $seqName, $flag,  $genome_name,
				$p,       $mapq,    $cigar, $rnext,
				$pnext,   $tlen,    $seq,   $seqQuality
			);
		}
	}
}
close $inQuery;
close $outFile;

sub bwa {
	my $query = shift @_;
	my ( @pos, @r, $i, $c, $mismatch );
	if ( scalar(@_) == 0 ) {

		#Initialization
		$i        = length($query) - 1;
		$c        = substr( $query, $i, 1 );
		@r        = ( 1, $charNum{$c} );
		$mismatch = 0;
	}
	else {
		$i        = shift @_;
		$c        = shift @_;
		$r[0]     = shift @_;
		$r[1]     = shift @_;
		$mismatch = shift @_;
	}
	
	while ( $i > 1 ) {
		@r = ( $charRange{$c}[0] + $r[0] - 1, $charRange{$c}[0] + $r[1] - 1 )
		  ;    #Find the range in F line
		#print $logFile "i = $i, c = $c, r = $r[0], $r[1]\n";
		$i--;
		$c = substr( $query, $i, 1 );
		my @charset2;
		if ($c ne 'N') {
			push @charset2, $c;
		}
		foreach my $k (@charset) {
			if ($k ne $c) {
				push @charset2, $k;	
			} 
		}
		foreach my $k (@charset2) {
			my @r2 = ( $tally{$k}[ $r[0] - 1 ], $tally{$k}[ $r[1] ] );
			#print $logFile "Try i = $i, c = $k, r2 = $r2[0], $r2[1]\n";
			if ( $r2[1] != $r2[0] ) {
				if ( $c eq $k ) {
					push @pos, bwa( $query, $i, $k, $r2[0]+1, $r2[1], $mismatch );
				}
				else {
					if ( $mismatch < $MAXMISMATCH ) {
						# $c = 'N' or mismatch
						push @pos,
						  bwa( $query, $i, $k, $r2[0]+1, $r2[1], $mismatch + 1 );
					}
				}
			}
			if (scalar @pos > $MAXMULTIMATCH) {
				last;
			}
		}
		return @pos;
	}
	#print $logFile "r = $r[0], $r[1]\n";
	@r = ( $charRange{$c}[0] + $r[0] - 1, $charRange{$c}[0] + $r[1] - 1 )
	  ;    #Find the range in F line
	#print $logFile "r = $r[0], $r[1]\n";
	$c = substr( $query, 0, 1 );
	if ( $tally{$c}[ $r[0] - 1 ] < $tally{$c}[ $r[1] ] ) {

		#Find the position of the query
		for ( my $i = $r[0] ; $i <= $r[1] ; $i++ ) {
			if ( $L[$i] eq $c ) {
	
				#Find the position!
				#print "$P[$i]\n";
				push @pos, $P[$i];
				if (scalar @pos > $MAXMULTIMATCH) {
					last;
				}
			}
		}
	}

	my $flag  = 0;
	my $mapq  = 42;
	my $cigar = length($query) . 'M';
	my $rnext = '*';
	my $pnext = '0';
	my $tlen  = '0';
	if (@pos) {
	}
	else {
		#push @pos, 0;
		$flag  = 0x4;
		$mapq  = 255;
		$cigar = '*';
	}
	@pos = sort(@pos);

	return $flag, \@pos, $mapq, $cigar, $rnext, $pnext, $tlen;
}

sub out2samHead {
	my ($outFile, $genome_name, $genome_length) = @_;
	print $outFile "\@HD\tVN:1.0\tSO:unsorted\n";
	print $outFile "\@SQ\tSN:$genome_name\tLN:$genome_length\n";
	print $outFile "\@PG\tID:mybowtie\tPN:mybowtie\tVN:1.2.0\n";
}

sub out2sam {
	my (
		$outFile, $qName, $flag,  $rName, $pos, $mapq,
		$cigar,   $rnext, $pnext, $tlen,  $seq, $seqQuality
	) = @_;

	print $outFile
"$qName\t$flag\t$rName\t$pos\t$mapq\t$cigar\t$rnext\t$pnext\t$tlen\t$seq\t$seqQuality\n";
}
