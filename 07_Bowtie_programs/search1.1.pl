#!/user/bin/perl
use warnings;

open (file1,"./lc.txt")||die "cannot open last_colume_file: $!";
open (file2,"./suffix.txt")||die "cannot open suffix_file: $!";
open (file3,"./totaltally.txt")||die "cannot open totaltally_file: $!";

my $lc="";
my @S;
my @a;

while($lcline=<file1>){
    chomp $lcline;
    $lc.=$lcline;
}
$si=0;

<file2>;
while($sline=<file2>){
    chomp $sline;
    $S[$si]=$sline;
    $si++;
}

$tti=0;
while($ttline=<file3>){
    chomp $ttline;
    if ($ttline eq "#"){
	$ttline = <file3>;
	chomp $ttline;
	@a=split(/,/,$ttline);
	$ai = $a[0];
	$ci = $a[1];
	$gi = $a[2];
        $ni = $a[3];
        $ti = $a[4];
	$reflen = $a[5];
	$refname = <file3>;
	chomp $refname;
}
    @a=split(/,/,$ttline);
    $tally{"A"}[$tti]=$a[0];
    $tally{"C"}[$tti]=$a[1];
    $tally{"G"}[$tti]=$a[2];
    $tally{"T"}[$tti]=$a[3];
    $tti++;
}

$left{"A"}=1;
$right{"A"}=$ai;
$left{"C"}=$ai+1;
$right{"C"}=$ai+$ci;
$left{"G"}=$ai+$ci+1;
$right{"G"}=$ai+$ci+$gi;
$left{"T"}=$ai+$ci+$gi+$ni+1;
$right{"T"}=$ai+$ci+$gi+$ni+$ti; 

##################get query
open (file5,"./reads_1.fq")||die "cannot open reads_1.fq_file: $!";
open (file7,"./lambda_chomp.fa")||die "cannot open ref_file: $!";
$ref=<file7>;
@ref=split(//,$ref);

my %results;
while($readname = <file5>){
    $readname=substr($readname,1);
    $readseq=<file5>;
    <file5>;
    $readquality=<file5>;
    chomp $readname;
    chomp $readseq;
    chomp $readquality;
    my $rlen=length($readseq);
    my $output;
    my $site=0;
    my $sfs;;
    my $sead;
    my $starts;
    my $tmm;
    $sfs=reverse($readseq);
    while ($sfs=~m/[ATCG]{24}/g){
	$starts=pos($sfs)-24;
	$sead=substr($sfs,$starts);
	$site = &single_read_search($sead);
    }
    $site=0 if $site<0;

    if ($site==0){
	$sfs=$readseq;
	$sfs=~tr/ATCG/TAGC/; 
	while ($sfs=~m/[ATCG]{24}/g){
	$starts=pos($sfs)-24;
        $sead=substr($sfs,$starts);
        $site = &single_read_search($sead);
	}
	$site=0 if $site<0;
	if ($site==0){$output="$readname\t4\t*\t0\t0\t*\t*\t0\t0\t$readseq\t$readquality\n";}
	else {$output="$readname\t16\t$refname\t$site\t42\t$rlen"."M\t*\t0\t0\t$sfs\t$readquality\n";}
	}
    else{
	$output="$readname\t0\t$refname\t$site\t42\t$rlen"."M\t*\t0\t0\t$sfs\t$readquality\n";	
    }
    $results{$output}=$site; 
}



open (file6,">./reads_1.sam")||die "cannot open reads_1.sam_file: $!";
print file6 "\@HD\tVN:1.0\tSO:sorted\n\@SQ\tSN:$refname\tLN:$reflen\n\@PG\tID:bowtie2\tPN:bowtie2\tVN:2.1.0\n";

foreach my $key(sort {$results{$a}<=>$results{$b}} keys%results){
    print file6 "$key";
}



sub single_read_search{

    $query=$_[0];
    my $qi=0;
    my $mmi=0;
    @q=split("",$query);
    $len=length($query);
    my $fi;
    my $li;
    my $cch;
    if($q[$qi] eq "N"){
	return 0;
    }
    my $left=$left{$q[$qi]};
    my $right=$right{$q[$qi]};
#####################main iteration
    my $mi=1;
    while ($qi<$len-1){
	$qi++;
	$cch=$q[$qi];
	if($cch eq "N"){
	    if ($left==$right){
		$cch=substr($lc,$left,1);
		return 0 if $cch eq "\$";
	    }
	    else{ return 0;}
	}
	$wlb=(int(($left-1)/100));
	$wrb=(int($right/100));
	$wordleft=substr($lc,$wlb*100+1,$left-$wlb*100-1);
	$wordright=substr($lc,$wrb*100+1,$right-$wrb*100);   
	$tofi= $wordleft =~ s/$cch/$cch/g || 0;
	$toli= $wordright =~ s/$cch/$cch/g || 0;
	$tofi++;
	$fi=$tally{$cch}[$wlb]+$tofi;
	$li=$tally{$cch}[$wrb]+$toli;    
	if ($fi>$li){ 
	    if ($left==$right){
		if($mmi==1){ return 0;}
		$cch=substr($lc,$left,1);
		if($cch eq "\$"){return 0;}
		$toli= $wordright=~s/$cch/$cch/g ||0;
		$li=$tally{$cch}[$wrb]+$toli;
		$fi=$li;
		$mmi++;
	    }
	    else{return 0;}
	}
	$left=$left{$cch}+$fi-1;
	$right=$left{$cch}+$li-1;
	$mi++;
	}
#######################match output
    if ($left == $right){
	$st=$S[$left];
  } 
    return $st;
}