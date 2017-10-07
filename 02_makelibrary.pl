#!/usr/bin/perl-w

use strict;

open FILE, '???.fa' || die "Can not open hg19.fa";     #open genome.fa  
local $/ = ">";     #设置读取结束标记为>, 默认为\n  
my $i=1;        #为输出文件编号  
while(<FILE>){  
        if(/^>$/){  
                next;  
        }  
        s/>//g; #去掉序列结束位置的>号
        s/(\w+)/>\1/; #为分割出来的序列第一行添加>号  
        my $name=$1; #以fasta序列第一行的一个单词作为文件名
        
        open OUTPUT, ">", '???'.'/'.$name; #write your words
        
        s/[\n\r]/ /g; #delate empty lines
        tr/>chr1234567890/           /; #delate names
        s/\s+//g;  #delate spaces
        $_ = uc($_);
        
        print OUTPUT $_; 
        close OUTPUT;  
        $i++;  
}  
close FILE;
