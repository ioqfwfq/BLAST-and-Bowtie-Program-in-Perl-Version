my $gtfFile = $ARGV[0];
my $readFileName = $ARGV[1];

my (%gtf, %gene_count, %exon_length);
my $total_mapped_reads = 0;

print "Read GTF files\n";
readGTF($gtfFile);

print "Count tophat sam result\n";
do_the_count($readFileName.'.un.th.sam');

print "Count bowtie sam result\n";
do_the_count($readFileName.'.sam');

print "Calculate and output RPKM\n";
output_rpkm($readFileName.'_rpkm.csv');


sub output_rpkm {
	my $F = shift;
	open ($ouF, '>', $F) or die "Can't open > $F!";
	print $ouF "gene_id,RPKM\n";
	for my $k (keys %gene_count) {
		my $rpkm = $gene_count{$k} * 10 ** 9 / $exon_length{$k} / $total_mapped_reads;
		print $ouF "$k,$rpkm\n";
	}
	close $ouF;
}

sub do_the_count {
	my $F = shift;
	open(SEG, '<', "$F") or die "Can't open < $F!";
	while (my $line = <SEG>) {
		if ($line =~ /^@/) {
			next;
		}
		my @temp = split("\t", $line);
		my $flag = int($temp[1]);
		my $mapq = int($temp[4]);
		if ($flag == 4 || $mapq < 10) {
			next;
		}
		my $chr = $temp[2];
		my $startp = int($temp[3]);
		my $geneid = $gtf{$chr}{$startp};
		my ($num_of_N, $len);
		if (exists($gene_count{$geneid})) {
			$gene_count{$geneid} ++; 
			$total_mapped_reads ++;
		}
		else {
			$num_of_N = $temp[5] =~ /M([0-9]+)N/;
			$endp = $startp + length($temp[9]) - 1 + $num_of_N;
			$geneid = $gtf{$chr}{$endp};
			if (exists($gene_count{$geneid})) {
				$gene_count{$geneid} ++;
				$total_mapped_reads ++;
			}
		}
	}
	close SEG;
}

sub readGTF {
	my $F = shift;
	open ($inF, '<', $F) or die "Can't open < $F!";
	while (<$inF>) {
		my @temp = split("\t");
		my ($chr, $startp, $endp, $geneid, $p);
		# print "$temp[2]\n";
		if ($temp[2] =~ /exon/) {
			$chr = $temp[0];
			$startp = int($temp[3]);
			$endp = int($temp[4]);
			$temp[8] =~ /gene_id \"([A-Za-z0-9]+)\"/;
			$geneid = $1;
			# print "$geneid, $chr, $startp, $endp\n";
			for ($p = $startp; $p <= $endp; $p++) {
				$gtf{$chr}{$p} = $geneid;
			}
			$gene_count{$geneid} = 0;			
			$exon_length{$geneid} += $endp - $startp + 1;
		}
	}
	close $inF;
}