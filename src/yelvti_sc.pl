#!/usr/bin/perl -w
use FindBin qw($Bin $Script);
use lib $FindBin::RealBin;
use strict;
use File::Basename qw(basename); 
use Getopt::Long;
my $BEGIN_TIME=time();
use CPSTAT::CPseq;
no warnings 'experimental::smartmatch';

#CPseq::run();

#######################################################################################

my @infile;
my $outdir;
my $is_to_normal;
my $is_compute_gc;
my $is_find_region;
my $seed;
my $min_len;

GetOptions(
				"help|?" =>\&USAGE,
				"i:s{1,}"=>\@infile,
				"s:s"=>\$seed,
				"m:s"=>\$min_len,
				"o:s"=>\$outdir,
				) or &USAGE;
&USAGE unless (@infile);

#######################################################################################
# ------------------------------------------------------------------
# Main Body                 
# ------------------------------------------------------------------
$seed ||=5100;
$min_len ||=100;
my @file_type =("fasta","fa","fsa","gb","gbk");

for my $infile(@infile){
	my $sample = CPseq::seq($infile);  #get name and seq;
	
	my $name = $sample->[0];
	my $seq = $sample->[1];
	print "------$name------\n";

	my $region = CPseq::region($seq,$seed,$min_len);
	if(!$region){
		print "no repeat larger than $min_len bp\n";
		next;
	}	
	print "$region\n\n";    #get old region
}



#######################################################################################
print "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n"; 
#######################################################################################

sub USAGE {         
	my $usage=<<"USAGE";				

  Options:
	-i	<infile>	
					input gbk/fasta files

	-h	help			
					Help show this information and exit

USAGE
	print $usage;
	exit;
}
