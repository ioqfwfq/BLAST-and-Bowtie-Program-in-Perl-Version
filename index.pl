#!/usr/bin/perl-w
#2017.03.09 this version is used for building BWT library of lambda_virus. copyright @Chris ZHU 

use warnings;

open FILE, 'lambda_virus.fa' || die "Can not open FASTA"; #打开文件

open OUT1,'>lambda.F';
open OUT2,'>lambda.L';
open OUT3,'>lambda.tly';

$line=<FILE>;
$line=<FILE>; chomp $line;
$content = "";
while ($line){
    $content=$content.$line;
    $line=<FILE>; chomp $line;
}

$content = $content.'$';
$len = length($content);
$i = 0;
for($i = 0; $i < $len; $i++){
    $up=substr($content,0,$i);
    $down=substr($content,$i);
    $seq{$down.$up}=$i;
}

$count_a=0;
$count_c=0;
$count_g=0;
$count_t=0;
foreach(sort keys %seq){
    $lastbase = substr($_,-1,1);
    print OUT1 "$lastbase";
    print OUT2 "$seq{$_}\n";
    $count_a++ if $lastbase eq 'A';
    $count_c++ if $lastbase eq 'C';
    $count_g++ if $lastbase eq 'G';
    $count_t++ if $lastbase eq 'T';
    print OUT3 "$count_a\t$count_c\t$count_g\t$count_t\n";
}

