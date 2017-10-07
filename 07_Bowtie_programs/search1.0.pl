#!/usr/bin/perl-w

open FILE, '../jiankuihe-exome-1.fq' || die "Can not open"; #打开文件
open FILE1,'../lambda_virus_lc.txt' || die "Can not open1";
open FILE2,'../lambda_virus_sa.txt' || die "Can not open2";
open FILE3,'../lambda_virus_tally.txt' || die "Can not open3";

$lc="";
while($lcline=<file1>){
    chomp $lcline;
    $lc.=$lcline;
}
print $lc;

$sai=0;
<file2>;
while($saline=<file2>){
    chomp $saline;
    $sa[$sai]=$saline;
    $sai++;
}
print $sa[0];

$tallyi=0;
while($tallyline=<file3>){
    chomp $tallyline;
    @a=split(/"\t"/,$tallyline);
    $tally{"A"}[$tallyi]=$a[0];
    $tally{"C"}[$tallyi]=$a[1];
    $tally{"G"}[$tallyi]=$a[2];
    $tally{"T"}[$tallyi]=$a[3];
    $tallyi++;
}
print $tally{"A"}[0];

$left{"A"}=1;
$right{"A"}=$tally{"A"}[$tallyi];
$left{"C"}=$right{"A"}+1;
$right{"C"}=$tally{"A"}[$tallyi]+$tally{"C"}[$tallyi];
$left{"G"}=$right{"C"}+1;
$right{"G"}=$tally{"A"}[$tallyi]+$tally{"C"}[$tallyi]+$tally{"G"}[$tallyi];
$left{"T"}=$right{"G"}+1;
$right{"T"}=$tally{"A"}[$tallyi]+$tally{"T"}[$tallyi]+$tally{"C"}[$tallyi]+$tally{"G"}[$tallyi];

$read = "AGCGCAGTGTCACTGCGCGCCTGTGCACTCTGTGGTGCTGCGGCCAGAATGCGGCGGGCCGTTTTCACGGTCATACCGGG";

@rd = split("", $read);

$site = &single_read_search($rd);

print "site".$site."site";

sub single_read_search{

    $query=$_[0];
    $qi=0;
    $mmi=0;
    @q=split("",$query);
    $len=length($query);
    $fi;
    $li;
    $cch;
    if($q[$qi] eq "N"){
    return 0;
    }
    $left=$left{$q[$qi]};
    $right=$right{$q[$qi]};
#####################main iteration
    $mi=1;
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