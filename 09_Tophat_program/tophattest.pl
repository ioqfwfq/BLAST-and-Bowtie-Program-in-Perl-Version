## Parameters Handling
my $genomeFile = $ARGV[0];
my $readFileName = $ARGV[1];
my $logFile = $readFileName.'_tophat.log';
open(LOG, '>', $logFile) or die "Can't open < $logFile!";
## Initialization
my $Num_Of_Segment = 4;
my $Length_Of_Segment = 25;
my $Max_Intron_Length = 60000;
my $Max_Splice_Site_Shift = 5;
my @Intron_Site = (["GT", "AG"], ["GC", "AG"], ["AT", "AC"], ["CT", "AC"]);
my (%readSeq, %readQual);
my %genome;
my %mapHash;
my %segment;
my $sam_head;

print ("Read Segment files into Hash table\n");
read_segment_sam("$readFileName.un.seg.sam");

print ("Read reads\n");
read_reads("$readFileName.un");

print "Load genome data\n";
read_reference_genome($genomeFile);

print "Find intron-spanning reads\n";
find_intron_spanning_reads();

sub find_intron_spanning_reads {
	my $ouF;
	open($ouF, '>', "$readFileName.un.th.sam") or die "Can't open < $readFileName.un.th.sam!";
	print $ouF $sam_head;
	
	my $hit;
	my $flag;
	for my $k (sort keys %readSeq) {
		$hit = $segment{$k};
		$flag = 0;
		if (defined($hit->[0]) && defined($hit->[3]) &&
			(defined($hit->[1]) || defined($hit->[2])) ) {
			# Match the pattern X X _ X
			for my $temp (split(',', $hit->[0])) {
				if ($flag) {
					last;
				}
				my ($chr, $dir, $pos) = split(':', $temp);
				$pos = int($pos);
				my $sign = 1;
				if ($dir eq '-') {
					$sign = -1;
				}
				my $test = "$chr\\:\\$dir\\:";
				my $pos1 = $pos + $sign * $Length_Of_Segment;
				my $test1 = $test . $pos1;
				if ($hit->[1] =~ m/$test1/) {
					for my $m ($hit->[3] =~ m/$test([0-9]+)/g) {
						$m = int($m);
						if (($sign > 0 && $m < $pos1) || ($sign < 0 && $m > $pos1)) {
							next;
						}
						my $diff = abs($m - $pos1);
						if ($diff < $Max_Intron_Length && $diff > 2*$Length_Of_Segment) {
							my $splice_pos = find_junction($readSeq{$k}, $chr, $sign, $pos, $m, 1, 3);
							if ($splice_pos > 0) {
								print LOG "$k: $dir segment 2, splice position: $splice_pos\n";
								to_sam($ouF, $k, $chr, $sign, $pos, $m, 2, $splice_pos);
								$flag = 1;
								last;
							}
						}
					}
				}
			}
			# Match the pattern X _ X X
			for my $temp ( split(',', $hit->[3]) ) {
				if ($flag) {
					last;
				}
				my ($chr, $dir, $pos) = split(':', $temp);
				$pos = int($pos);
				my $sign = 1;
				if ($dir eq '-') {
					$sign = -1;
				}
				my $test = "$chr\\:\\$dir\\:";
				my $pos1 = $pos - $sign * $Length_Of_Segment;
				my $test1 = $test . $pos1;
				if ($hit->[2] =~ m/$test1/) {
					for my $m ($hit->[0] =~ m/$test([0-9]+)/g) {
						$m = int($m);
						if (($sign > 0 && $m > $pos1) || ($sign < 0 && $m < $pos1)) {
							next;
						}
						my $diff = abs($m - $pos1);
						if ($diff < $Max_Intron_Length && $diff > 2*$Length_Of_Segment) {
							my $splice_pos = find_junction($readSeq{$k}, $chr, $sign, $m, $pos, 0, 2);
							if ($splice_pos > 0) {
								print LOG "$k: $dir segment 1, splice position: $splice_pos\n";
								to_sam($ouF, $k, $chr, $sign, $m, $pos, 1, $splice_pos);
								$flag = 1;
								last;
							}
						}
					}
				}
			}					
		}	
	}
	close(ouF);
}

sub find_junction{
	my $MAX_MISMATCH = 3;
	my ($seq, $chr, $sign, $pos, $pos1, $ni, $nj) = @_;
	my $junc = uc(substr($seq, ($ni+1) * $Length_Of_Segment, $Length_Of_Segment));
	my $conc;
	if ($sign > 0) {
		$conc = uc(substr($genome{$chr}, $pos - 1 + ($ni+1) * $Length_Of_Segment, $Length_Of_Segment) . 
			substr($genome{$chr}, $pos1 - 1 - ($Num_Of_Segment - $nj) * $Length_Of_Segment, $Length_Of_Segment));		
	}
	else {
		$conc = uc(substr($genome{$chr}, $pos1 - 1 + ($Num_Of_Segment - $nj) * $Length_Of_Segment, $Length_Of_Segment) .
			substr($genome{$chr}, $pos - 1 - ($ni+1) * $Length_Of_Segment, $Length_Of_Segment));		
		# Reverse complement
		$junc =~ tr/ACGT/TGCA/;
		$junc = reverse $junc;
	}
	print LOG "\nSeq = $seq, Chr = $chr, Sign = $sign, Pos = $pos, Pos1 = $pos1, ni = $ni, nj = $nj\nJunc = $junc\nConc = $conc\n";
	my $j = length($junc) - 1;
	my $k = length($conc) - 1;
	my $i = 0;
	my $mismatch = 0;
	my $flag = 0;
	while ($j - $i > 1) {
		while ($j - $i > 1 && substr($junc, $i, 1) eq substr($conc, $i, 1)) {
			$i ++;
		}
		while ($j - $i > 1 && substr($junc, $j, 1) eq substr($conc, $k, 1)) {
			$j --;
			$k --;
		}
		if ($j - $i > 1 && $mismatch < $MAX_MISMATCH) {
			if ($flag) {
				$i ++;
				$flag = 0;
			}
			else {
				$j ++;
				$flag = 1;
			}
			$mismatch ++;
		}
		if ($mismatch == $MAX_MISMATCH) {
			# Can't find splicing pattern, return with -1
			print LOG "Fail: i = $i, j = $j, k = $k, mismatch = $mismatch\n";
			return -1;
		}
	}
	print LOG "i = $i, j = $j, k = $k, mismatch = $mismatch\n";

	# Find Intro splicing tag and correct the current result	
	$i = splice_shift($junc, $conc, $i, $j);

	return $i;
}

sub splice_shift{
	my ($junc, $conc, $i, $j) = @_;
	my ($p1, $p2, $s1, $s2, $ss1, $ss2);
	my $flag = 0;
	for (my $t = 0; $t <= $Max_Splice_Site_Shift; $t++ ){
		#Test shift left
		$p1 = $i - $t;
		$p2 = $Length_Of_Segment + $j - $t - 1;
		if ($p1 >= 0) {
			$s1 = substr($junc, $p1, $t);
			$s2 = substr($conc, $p2, $t);
			print LOG "left t = $t, p1 = $p1, p2 = $p2, s1 = $s1, s2 = $s2\n";			
			$ss1 = substr($conc, $p1, 2);
			$ss2 = substr($conc, $p2-2, 2);
			print LOG "ss1 = $ss1, ss2 = $ss2\n";
			if ($s1 eq $s2) {
				for my $intron (@Intron_Site) {
					if ($ss1 eq $intron->[0] && $ss2 eq $intron->[1]) {
						print LOG "--> intron0 = $intron->[0], intron1 = $intron->[1]\n";
						$i -= $t ;
						$j -= $t;
						$flag = 1;
						last;
					}
				}
			}
			if ($flag) {
				last;
			}
		}
		
		# Test shift right
		$p1 = $i;
		$p2 = $j - 1;
		if ($p2 + $t < $Length_Of_Segment) {
			$s1 = substr($conc, $p1, $t);
			$s2 = substr($junc, $p2, $t);
			print LOG "right t = $t, p1 = $p1, p2 = $p2, s1 = $s1, s2 = $s2\n";			
			$ss1 = substr($conc, $p1 + $t, 2);
			$ss2 = substr($conc, $Length_Of_Segment + $p2 + $t - 2, 2);
			print LOG "ss1 = $ss1, ss2 = $ss2\n";
			if ($s1 eq $s2) {
				for my $intron (@Intron_Site) {
					if ($ss1 eq $intron->[0] && $ss2 eq $intron->[1]) {
						print LOG "--> intron0 = $intron->[0], intron1 = $intron->[1]\n";
						$i += $t ;
						$j += $t;
						$flag = 1;
						last;
					}
				}
			}
			if ($flag) {
				last;
			}
		}
	}
	return $i;
} 


sub to_sam {
	my ($outFile, $qName, $rName, $sign, $pos, $pos1, $seg_num, $splice_pos) = @_;
	my ($seq, $qual, $mapq, $rnext, $pnext, $tlen) = ($readSeq{$qName}, $readQual{$qName}, 43, '*', 0, 0);
	my ($cut_pos, $pos_donar, $pos_acceptor, $num_of_N, $cigar);

	if ($sign > 0) {
		$cut_pos = $seg_num * $Length_Of_Segment + $splice_pos;
		$pos_donar = $pos;
		$pos_acceptor = $pos1 - ($Num_Of_Segment - $seg_num - 2) * $Length_Of_Segment - ($Length_Of_Segment - $splice_pos);
		$num_of_N = $pos_acceptor - $pos_donar - $cut_pos;
		$flag = 0;
		$cigar = $cut_pos. 'M' . $num_of_N . 'N' . (length($seq) - $cut_pos) . 'M';
		print $outFile "$qName\t$flag\t$rName\t$pos_donar\t$mapq\t$cigar\t$rnext\t$pnext\t$tlen\t$seq\t$qual\n";
	}
	elsif ($sign < 0) {
		$seq =~ tr/ACGT/TGCA/;
		$seq = reverse $seq;
		$qual = reverse $qual;
		$cut_pos = ($Num_Of_Segment - $seg_num - 1) * $Length_Of_Segment + $splice_pos;
		$pos_donar = $pos - ($seg_num - 1) * $Length_Of_Segment - ($Length_Of_Segment - $splice_pos);
		$pos_acceptor = $pos1;
		$num_of_N = $pos_donar - $pos_acceptor - $cut_pos;
		$flag = 16;
		$cigar = $cut_pos. 'M' . $num_of_N . 'N' . (length($seq) - $cut_pos) . 'M';
		print $outFile "$qName\t$flag\t$rName\t$pos_acceptor\t$mapq\t$cigar\t$rnext\t$pnext\t$tlen\t$seq\t$qual\n";
	}
}

sub read_reference_genome{
	my $F = shift;
	open(GEO, '<', $F) or die "Can't open < $F!"; 
	my $curGeo = "";
	my $line;
	while ($line = <GEO>) {
		chomp($line);
		if ($line =~ m/>/) {
			$curGeo = substr($line, 1);
			print "$curGeo\n";
			$genome{$curGeo} = "";	
			next;
		}
		$genome{$curGeo} .= $line;
	}
	close GEO;
}

sub read_segment_sam {
	my $F = shift; # read $readFile all sam
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
}

sub read_reads {
	my $F = shift;
	open(RDF, '<', $F) or die "Can't open < $F!";
	while (<RDF>) {
		my @temp = split(" ");
		my $read_name = $temp[0];
		$read_name = substr($read_name, 1);
		my $read = <RDF>;
		chomp($read);
		my $middle_bar = <RDF>;
		my $quality = <RDF>;
		chomp($quality);
		$readSeq{$read_name} = $read;
		$readQual{$read_name} = $quality;
	}
	close(RDF);
}