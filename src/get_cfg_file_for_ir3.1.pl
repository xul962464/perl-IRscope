#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use FindBin qw/$Bin/;
my $BEGIN_TIME=time();

my ($indir,$outfile,@infile,);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s{1,}"=>\@infile,
				"o:s"=>\$outfile,
				) or &USAGE;
&USAGE unless ($indir or @infile);
######################################################################################################

#获取文件

$outfile ||="irscope.v3.1.cfg";
open OUT ,"> $outfile" or die"$!";
for my $file(@infile){
	chomp (my $path = `readlink -f $file`);
	chomp (my $ir = `perl $Bin/yelvti_sc.pl -i $file  | grep "^IR" |cut -d ":" -f 2 |paste - - -d ","`);
	print OUT "$path\t$ir\n";
}

print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
sub USAGE {         
	my $usage=<<"USAGE";				
Program: 
Version: 2.1.0
Contact: xul<xul\@genepioneer.cn> <1275875706\@qq.com>
Description:
Usage:
  Options:
	-i	<infile>	input gbk_files
	-o	<output>	out.dir outfile default irscope.v3.1.cfg
	-h				Help
USAGE
	print $usage;
	exit;
}
