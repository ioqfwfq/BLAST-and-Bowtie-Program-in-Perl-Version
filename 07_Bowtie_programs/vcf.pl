#!/user/bin/perl
use warnings;

my ( $pileupFile, $outFile ) = @ARGV;
vcf( $pileupFile, $outFile );

sub vcf {
	my ( $pileupFile, $outFile ) = @_;
	open( PI, '<', $pileupFile )
	  or die "cannot open < $pileupFile: $!";
	open( OUT, '>', $outFile )
	  or die "cannot open > $outFile: $!";
	print OUT "CHROM\tPOS\tREF\tALT\tQUAL\n"."\n";

	my $minReads = 3;
	my $minPurity = 0.5;
	my %count;
	my $num = 0;
	while (<PI>) {
		chomp;
		
		my ( $chr, $pos, $ref, $readsNum, $reads, $quali) = split("\t");
		undef %count;
		if ($readsNum >= $minReads) {
			for (my $i = 0; $i < $readsNum; $i++) {
				my $c = substr($reads, $i, 1);
				if ($c ne '.' && $c ne 'N') {
					if (exists($count{$c})) {
						$count{$c} ++;
					}
					else {
						$count{$c} = 1;
					}
				}
			}
			for my $k (keys %count) {
				my $purity = $count{$k} / $readsNum; 
				if ( $purity >= $minPurity) {
					#Find a SNP
					if ($purity == 1) {
						$purity = 0.9999;
					}
					$num ++;
					my $PhredScore = -10 * log(1-$purity)/log(10);
					my $printPScore = sprintf("%.2f", $PhredScore);
					print OUT "\n", $chr, "\t", $pos, "\t", $ref, "\t", $k, "\t", $printPScore, "\n"; 
				}
			}
		}	
	}
	print OUT $num;
}
