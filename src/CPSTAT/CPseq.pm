#!/usr/bin/perl

package CPseq;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(run region seq to_normal_seq gc);

use warnings;
use strict;
use File::Basename qw(basename); 
no warnings 'experimental::smartmatch';

=head1 ****************************************************

=head1 **********Contact: jiangcb@genepioneer.com**********

=head1 ****************************************************

=head1 NAME AND INTRODUCTION
	CPseq:
		Find Chloroplast Region
	include:
		&seq
		&region
		&to_normal_seq
		&gc

=head2 usage:
	load model:
		use CPSTAT::CPseq
	
=head3 &seq:
		Get the sequence from the file.

		my $sample = CPseq::seq($file);
			
		then $sample contains a reference to two elements
		$sample -> [0]:is file basename
		$sample -> [1]:is seq

		the file can be gbk(.gb||.gbk) file or fasta(.fa||.fsa||.fasta) file

=head3 &region:
		Find four region in seq(LSC SSC IR).
		
		my $region = CPseq::seq($seq);

		$region:
		like this:
			LSC:1-86120
			IRb:86121-112511
			SSC:112512-131531
			IRa:131532-157922

		or this(If there is a region that is not continuous):
			LSC:157867-157922,1-86064
			IRb:86065-112455
			SSC:112456-131475
			IRa:131476-157866

=head3 &to_normal_seq
		Convert a sequence to a normal sequence(LSC,IRb,SSC,IRa) and all areas are continuous.
		
		my $normal_seq = CPseq::seq($seq);
		
		$normal_seq is the seq that is continuous

=head3 &gc
		Stactistic every region GC content

		my $GC = CPseq::gc($seq)
		
		$GC:
			full length Nucleotide Statistics
			length: 153482 bp
				    Size      %
			A:      48114   31.35%
			C:      28407   18.51%
			G:      27400   17.85%
			T:      49561   32.29%
			GC:     55807   36.36%
			All:    153482  100.00%

			LSC Nucleotide Statistics
			length: 83281 bp
				    Size      %
			A:      26693   32.05%
			C:      14614   17.55%
			G:      13805   16.58%
				...
=cut

sub gc{
	my $seq = shift;
	my $new_seq = to_normal_seq($seq);
	my $region = region($new_seq);
	my @regions = $region =~ /(\d+)/mg;
	my $content;

	my $lsc=[$regions[0],$regions[1],'l','a','a%','t','t%','g','g%','c','c%','LSC'];
	my $IRb=[$regions[2],$regions[3],'l','a','a%','t','t%','g','g%','c','c%','IRb'];
	my $ssc=[$regions[4],$regions[5],'l','a','a%','t','t%','g','g%','c','c%','SSC'];
	my $IRa=[$regions[6],$regions[7],'l','a','a%','t','t%','g','g%','c','c%','IRa'];

	my $all_length=[$regions[0],$regions[7],'l','a','a%','t','t%','g','g%','c','c%','full length'];
	
	compute($_,$new_seq) for (($all_length,$lsc,$IRb,$ssc,$IRa));
	

	for (($all_length,$lsc,$ssc,$IRb,$IRa)){
		$content .= out($_);
	}
	return $content;
}

sub compute{  #compute bp
	my $region = shift;
	my $seq = shift;
	my ($nu1,$nu2);

	$region->[2]=($nu1=substr($seq,$region->[0]-1,$region->[1]-$region->[0]+1)=~tr/GCTA/GCTA/);
	$region->[3]=($nu1=substr($seq,$region->[0]-1,$region->[1]-$region->[0]+1)=~tr/A/A/);
	$region->[4]=($nu1=substr($seq,$region->[0]-1,$region->[1]-$region->[0]+1)=~tr/A/A/)/($nu2=substr($seq,$region->[0]-1,$region->[1]-$region->[0]+1)=~tr/GCTA/GCTA/)*100;
	$region->[5]=($nu1=substr($seq,$region->[0]-1,$region->[1]-$region->[0]+1)=~tr/T/T/);
	$region->[6]=($nu1=substr($seq,$region->[0]-1,$region->[1]-$region->[0]+1)=~tr/T/T/)/($nu2=substr($seq,$region->[0]-1,$region->[1]-$region->[0]+1)=~tr/GCTA/GCTA/)*100;
	$region->[7]=($nu1=substr($seq,$region->[0]-1,$region->[1]-$region->[0]+1)=~tr/G/G/);	
	$region->[8]=($nu1=substr($seq,$region->[0]-1,$region->[1]-$region->[0]+1)=~tr/G/G/)/($nu2=substr($seq,$region->[0]-1,$region->[1]-$region->[0]+1)=~tr/GCTA/GCTA/)*100;
	$region->[9]=($nu1=substr($seq,$region->[0]-1,$region->[1]-$region->[0]+1)=~tr/C/C/);	
	$region->[10]=($nu1=substr($seq,$region->[0]-1,$region->[1]-$region->[0]+1)=~tr/C/C/)/($nu2=substr($seq,$region->[0]-1,$region->[1]-$region->[0]+1)=~tr/GCTA/GCTA/)*100;
}

sub out{     #print
	my $region = shift;

	my $result =sprintf "$region->[11] Nucleotide Statistics\nlength: $region->[2] bp\n\tSize\t  %%\n".
	"A:\t$region->[3]\t%.2f%%\n".
	"C:\t$region->[9]\t%.2f%%\n".
	"G:\t$region->[7]\t%.2f%%\n".
	"T:\t$region->[5]\t%.2f%%\n".
	"GC:\t%d\t%.2f%%\n".
	"All:\t$region->[2]\t%.2f%%\n\n",
	$region->[4],$region->[10],$region->[8],$region->[6],$region->[7]+$region->[9],$region->[10]+$region->[8],$region->[4]+$region->[10]+$region->[8]+$region->[6];	

	return $result; 
}

sub to_normal_seq{
	my $seq = shift;
	my $region = region($seq);
	my @region = split"\n",$region;
	
	my $new_seq;

	for(@region){
		my @pos = /(\d+)/g;
		for my $i(0..(@pos/2-1)){
			$new_seq .= extract_seq($seq,$pos[2*$i],$pos[2*$i+1]);
		}
	}
	return $new_seq;
}

sub extract_seq{
	my $seq = shift;
	my $s = shift;
	my $e = shift;
	return substr($seq,$s-1,$e-$s+1);
}

sub region{
	my $seq = shift;
	my $len_of_seed = shift;
	my $min_len = shift;
	$min_len ||= 100;
	my $re_seq = reverse $seq;
	$re_seq =~ tr/ATGC/TACG/;
#	print "$len_of_seed  $min_len\n";
	my $region1 = find_candidate_region($seq,$len_of_seed,$min_len);
	return 0 if(!$region1);
	my $region2 = find_candidate_region($re_seq,$len_of_seed,$min_len);
	my $lsc;
	my $ssc;
	my $ira;
	my $irb;
	my $result = Judge($region1,$region2,$seq);
	
	if($result =~ /,/){
		my @sort = $result=~/(\d+)/g;
		
		if($sort[0] - $sort[-1] > $sort[2] - $sort[1]){
				$lsc = "LSC:".($sort[-1]+1)."-".($sort[0]-1);
				$irb = "IRb:$sort[0]-$sort[1]";
				$ssc = "SSC:".($sort[1]+1)."-".($sort[2]-1);
				$ira = "IRa:$sort[2]-$sort[3],$sort[4]-$sort[5]";
			}else{
				$lsc = "LSC:".($sort[1]+1)."-".($sort[2]-1);
				$irb = "IRb:$sort[2]-$sort[3],$sort[4]-$sort[5]";
				$ssc = "SSC:".($sort[-1]+1)."-".($sort[0]-1);
				$ira = "IRa:$sort[0]-$sort[1]";
			}

	}else{
		my @sort = sort{$a <=> $b} $result=~/(\d+)/g;

		if($sort[-1] == length $seq){
			if($sort[2] - $sort[1] < $sort[0]){
				$lsc = "LSC:1-".($sort[0]-1);
				$irb = "IRb:$sort[0]-$sort[1]";
				$ssc = "SSC:".($sort[1]+1)."-".($sort[2]-1);
				$ira = "IRa:$sort[2]-$sort[3]";
			}else{
				$lsc = "LSC:".($sort[1]+1)."-".($sort[2]-1);
				$irb = "IRb:$sort[2]-$sort[3]";
				$ssc = "SSC:1-".($sort[0]-1);
				$ira = "IRa:$sort[0]-$sort[1]"; 
			}
		}else{
			if($sort[0]-1 == 0){
				if(($sort[2] - $sort[1]) < ($sort[0] + length($seq) - $sort[-1])){
					$lsc = "LSC:".($sort[-1]+1)."-".length($seq);
					$irb = "IRb:$sort[0]-$sort[1]";
					$ssc = "SSC:".($sort[1]+1)."-".($sort[2]-1);
					$ira = "IRa:$sort[2]-$sort[3]";
				}else{
					$lsc = "LSC:".($sort[1]+1)."-".($sort[2]-1);
					$irb = "IRa:$sort[2]-$sort[3]";
					$ssc = "SSC:".($sort[-1]+1)."-".length($seq);
					$ira = "IRb:$sort[0]-$sort[1]";
				}
			}else{
				if(($sort[2] - $sort[1]) < ($sort[0] + length($seq) - $sort[-1])){
					$lsc = "LSC:".($sort[-1]+1)."-".length($seq).",1-".($sort[0]-1);
					$irb = "IRb:$sort[0]-$sort[1]";
					$ssc = "SSC:".($sort[1]+1)."-".($sort[2]-1);
					$ira = "IRa:$sort[2]-$sort[3]";
				}else{
					$lsc = "LSC:".($sort[1]+1)."-".($sort[2]-1);
					$irb = "IRa:$sort[2]-$sort[3]";
					$ssc = "SSC:".($sort[-1]+1)."-".length($seq).",1-".($sort[0]-1);
					$ira = "IRb:$sort[0]-$sort[1]";
				}
			}
			
		}
	}
	return join("\n",($lsc,$irb,$ssc,$ira));
}

sub find_candidate_region{
	my $seq = shift;
	my $len_of_seed = shift;
	my $min_len = shift;
	$len_of_seed ||=5100;	
	my $re_seq = reverse $seq;
	$re_seq =~ tr/ATGC/TACG/;
	my $seq_long = length($seq);
	my $pos;
	my $pos2;
	my $seed;
	my $start;
	my $end;
	my $start2;
	my $end2;
	my $long;
	my $len;
	my $len2;
	my $flag = 0;

RESTAT:	

	my $nu_of_seed = int((length $seq)/$len_of_seed);

	for my $n(1..$nu_of_seed){
		$seed = substr($re_seq,$len_of_seed*($n-1),$len_of_seed);    #find seed
		$start = index($seq,$seed);	
		
		if($start < 0 and $n == $nu_of_seed){
			$len_of_seed = $len_of_seed/2;
			if($len_of_seed < $min_len){
				warn "no repeat larger than $min_len bp\n";
				return 0;	
			}
			goto RESTAT;
		}

		next if($start < 0);
		
		while( $start > 0){
			$long .= $seed;				#merge seed to long
			$seed = substr($re_seq,$len_of_seed*(++$n-1),$len_of_seed);
			$start = index($seq,$seed);
		}
	
		$pos  = index($re_seq,$long);    #find long pos in re_seq 
		
		$pos2 = $pos;

		$len = length($long);		#length of long;			
		$len2 = length($long);

		$start  = index($seq,$long);	#find long pos in seq 

		while(index($seq,$long)>=0 and $pos >=0){	#to left
			$start = $start - 1 ;
			$long = substr($re_seq,--$pos,$len) ;  
		}

		if(-1 == $pos){								#judge tail seq
			$long = substr($re_seq,0,$len);
			my $all_long = $long;
			$seq_long = $seq_long - 1;
			$long = substr($re_seq,$seq_long,length($seq) - $seq_long) . $all_long ;
			goto BBB if (index($seq,$long)<0);

			while(index($seq,$long)>=0){
				$start = $start - 1 ;
				$seq_long = $seq_long - 1;
				$long = substr($re_seq,$seq_long,length($seq) - $seq_long) . $all_long ; 
				$flag = 1;
			}

			$start = $start + 2;
			$start2 = $seq_long + 2;
		}else{
			BBB:
			$start = $start + 2;
			$start2 = $pos + 2;
		}
		
		$long = substr($re_seq,$pos2,$len2);		
	
		while(index($seq,$long)>=0){
			$end = index($seq,$long) + $len2;
			$long = substr($re_seq,$pos2,++$len2) ;		#to right
		}
		
		if($end == length($seq)){						
			$long = substr($re_seq,$pos2);
			my $all_long = $long;
			my $i = 1 ;
			$long = $all_long . substr($re_seq,0,$i);
			while(index($seq,$long)>=0){
				$i++;
				$long = $all_long . substr($re_seq,0,$i);
				$end = index($seq,$long) + $len2;
				$end2 = $i;
				$flag = 2;
			}
			
		}else{
			$end2 = $pos2 + $len2 - 1 ;
		}	
		
		my $tmp_e2 = $end2;

		$end2  = length($seq) - $start2 +1;
		$start2 = length($seq) - $tmp_e2 +1;

		if($flag == 0){
			return "$start-$end\n$start2-$end2";
		}else{
			return "$start-$end\n$start2-".length($seq).",1-$end2";
		}
	}
}

sub Judge{
	my $result1 = shift;
	my $result2 = shift;
	#print "$result1\n\n";	
	#print "$result2\n";
	my $seq = shift;
	
	my $len = length $seq;

	my $r1 = (split"\n",$result1)[0];
	my $r2 = (split"\n",$result2)[0];

	my $r22 = (split"\n",$result2)[1];
	
	my $s1 = (split"-",$r1)[0];
	my $s2 = (split"-",$r2)[0];

	my $e1 =  (split"-",$r1)[1]; 
	my $e2 =  (split"-",$r2)[1]; 

	my $size1 = $e1 - $s1;
	my $size2 = $e2 - $s2;
	
	if ($size1 >= $size2) {
		return $result1;
	}else{
		my $tmp_e2 = $e2;
		$e2 = $len - $s2 +1;
		$s2 = $len - $tmp_e2 +1;
		
		my ($s21,$e21,$s22,$e22) = $r22 =~/(\d+)/g;

		return "$s2-$e2\n".($len-$e22+1)."-".($len-$s22+1).",".($len-$e21+1)."-".($len-$s21+1);
	}
}

sub gbk2fa{
	my $infile = shift;
	open IN ,$infile or die "Could not open $infile $!";

	$/="ORIGIN";
	<IN>;
	my $seq=<IN>;
	close IN;
	$/="\n";
	$seq=~ s/\/\/|\d+|\s+|\n//g;#将序列中的空白，前面的数字，换行符，最后的//全部换成空--及seq是一条完整的序列
	$seq="\U$seq";
}

sub fa2fa{
	my $infile = shift;
	open IN,$infile or die"$!";
	my $seq;

	while(<IN>){
		chomp;
		next if(/^\s+$/);
		s/\s+//g;
		if(/>(.*)/){
			1;
		}else{
			$seq .= "\U$_";
		}
	}close IN;
	return $seq;
}

sub seq{
	my $file = shift;
	my $filename = basename $file;
	my ($type) = $filename =~ /\.(\w+?)$/;
	my $seq_name = $`;
	my @fasta = ("fasta","fa","fsa");
	my @gbk = ("gb","gbk");
	if($type ~~ @fasta){
		[$seq_name,fa2fa($file)];
	}elsif($type ~~ @gbk){
		[$seq_name,gbk2fa($file)];
	}else{
		die"Please chech your file name\n";
	}
}

sub run{
print"
*************************************************************
***************Hello World! Hello Gene Pioneer***************
*************************************************************\n\n";
} 

1;
