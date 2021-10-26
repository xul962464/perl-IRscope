#!/usr/bin/perl -w
use strict;
use SVG;
use Getopt::Long;
use File::Basename qw(basename dirname);

#############################################################################################################
#############################################################################################################
#############################################################################################################
#需要改的参数
#定义颜色
my %color=(	
		#定义每个区域的颜色
			"LSC" => "#D1EEEE",
			"IRb" => "#EEAD0E",
			"SSC" => "#B4EEB4",
			"IRa" => "#EEAD0E",

		#定义基因的颜色，未定义的话会默认使用红色，% 表示假基因
			"rps19" => "#228B22",
			"ycf1"=> "#1E90FF",
			"ndhF" => "#8B4500",
			"trnH" => "#8B8B83",
			"trnN" => "#8B7765",
			"rpl21" => "#00F5FF",
			"trnI" => "#8B8682",
			"rpl22" => "#836FFF",
			"rps15" => "#9BCD9B",
			"psbA" => "#9A32CD",
			"trnV" => "#adb7f5",
			"trnR" => "#6E8B3D",
			"rpl32" => "#D02090",
			"rps7" => "#CD3333",
			"rpl23" => "#9370DB",
			"rpl2" => "#FF4500",
			"trnL" => "#FF6EB4",
			"ccsA" => "#4F94CD",
			"trnK" => "#4F4F4F",
			"%rps19" => "#228B22",
			"%ycf1"=> "#1E90FF",  
			"%ndhF" => "#8B4500", 
			"ndhH" => "#EE9A00",
			"%trnH" => "#8B8B83", 
			"%trnN" => "#8B7765", 
			"%rpl22" => "#836FFF",
			"%rps15" => "#9BCD9B",
			"%psbA" => "#9A32CD", 
			"%trnV" => "#adb7f5", 
			"%trnR" => "#6E8B3D", 
			"%rpl32" => "#D02090",
			"%rps7" => "#CD3333", 
			"%rpl23" => "#9370DB",
			"%rpl2" => "#FF4500", 
			"%trnL" => "#FF6EB4", 
			"%ccsA" => "#4F94CD", 
			"%trnK" => "#4F4F4F", 
			"cbbx" => "#228B22",
			"rrns"=> "#1E90FF",  
			"rrn5" => "#8B4500", 
			"ycf19" => "#8B8B83", 
			"rpl21" => "#8B7765", 
			"ycf37" => "#836FFF",
			"trnA" => "#9BCD9B",
			"rrnl" => "#9A32CD", 
			"rpl16" => "#2E8B57", 
			"rps3" => "#DAA520", 
);


my $BEGIN_TIME=time();

my ($indir,$outfile,$infile,$left,$is_print_gene_len,$pseudo);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$infile,
				"l:s"=>\$left,
				"p!"=>\$is_print_gene_len,	#if($is_print_gene_len)
				"pseudo!"=>\$pseudo,
				"o:s"=>\$outfile,
				) or &USAGE;
&USAGE unless ($infile);
######################################################################################################


#从文件中获取文件及其对应的重复区
open IN,"$infile" or die;

my @gbfiles;

while(<IN>){
	chomp;
	next if(/^\s*$|^#/);
	my	$gb = (split)[0];
	my  $ir = (split)[1];
	
	if($ir =~ /^(\d+)-(\d+),(\d+)-(\d+)$){
		if($1 < $2 and $2 < $3 and $3 < $4){
			1;
		}else{
			print "#" x 80 ,"\n";
			print "this file :$gb IR pos is reverse. please check it;";
			die;
		}
	}else{
		print "#" x 80 ,"\n";
		print "this file :$gb IR region($ir) maybe wrong. please check it;\n";
		die;
	}
	push @gbfiles,[$gb,$ir];
}




#判断最长物种名称,确定左边留多少空白
my $longest_name = 0;

for my $info(@gbfiles){
	my $file = $info->[0];
	my $name = get_sample_name($file);
	my $name_len = length $name;
	$longest_name = $name_len > $longest_name ? $name_len:$longest_name;
}


#左端距离，左边距
$left ||=120;
my $zx = $left + $longest_name*3;


#画布大小
my $svg = SVG->new(width=>1300+$zx,height=>180+120*$#gbfiles);  #画布

#输出文件
$outfile ||="irscope";
my $pre = basename $outfile;

if($outfile and $outfile !~ /\.svg$/){
	$outfile .=".svg";
}
open OUT,">$outfile" or die"$!";

#作图
my $file_nu = 0;
my $y = 0;
$|= 1;
for my $info(@gbfiles){

	my $file = $info->[0];
	my $ir = $info->[1];

	#文件名称 及数量
	print "#########".($file_nu+1)."-$file########\n";	
	#获取样品名
	my $sample_name = get_sample_name($file);	
	
	#获取序列
	my $seq = get_seq($file);		
	my $seq_len = length $seq;


	#获取区域
	my $start_lsc;
	my $end_lsc ;
	my $start_irb;
	my $end_irb;
	my $start_ssc;
	my $end_ssc;
	my $start_ira;
	my $end_ira;
	my $T = length $seq;	#环形序列，长度相当于周期

	if($ir =~ /(\d+)-(\d+),(\d+)-(\d+)/){	#文件中给的重复区位置
		$start_lsc = 1;
		$end_lsc = $1 - 1;
		$start_irb = $1;
		$end_irb = $2;
		$start_ssc = $2 + 1;
		$end_ssc = $3 - 1;
		$start_ira = $3;
		$end_ira = $4;
	}else{	#如果没给位置，会尝试自己计算
		print"##WARN:not find IR in your config file,try find it use seq\n";
		my @region = find_region($seq);
		die unless(@region);
		$start_lsc = 1;
		$end_lsc = $region[0] - 1;
		$start_irb = $region[0];
		$end_irb = $region[1];
		$start_ssc = $region[1]+1;
		$end_ssc = $region[2] - 1;
		$start_ira = $region[2];
		$end_ira = $region[3];
	}

	#LSC增加了可能存在末尾的序列
	my @region_len = ($end_lsc-$start_lsc+1+$T-$end_ira,$end_irb-$start_irb+1,$end_ssc-$start_ssc+1,$end_ira-$start_ira+1);
	print"-" x 40;
	print"\n";
	print"LSC:",($start_lsc-$T+$end_ira-1),"-$end_lsc\nIRb:$start_irb-$end_irb\nSSC:$start_ssc-$end_ssc\nIRa:$start_ira-$end_ira\n" if($T - $end_ira != 0);
	print"LSC:",($start_lsc),"-$end_lsc\nIRb:$start_irb-$end_irb\nSSC:$start_ssc-$end_ssc\nIRa:$start_ira-$end_ira\n" if($T - $end_ira == 0);
	print"-" x 40;
	print"\n";

	#获取基因位置信息
	my %pos_for_gene = get_gene_pos($file,$T);	#T 为序列长度，一个周期

	#每个样品的间隔
	$y = 120*($file_nu++);
	my @region_name=qw/LSC IRb SSC IRa LSC/;

	#样品名称
	$svg->text(x=>$zx+130,y=>100+$y,'font-size'=>15,'font-family'=>'times new roman','font-style'=>'italic','-cdata'=>$sample_name,"text-anchor"=>"end");


	#区域的长度
	my $len_ir = 250;
	my $len_lsc = 150;
	my $len_ssc	= 300;

	#区域垂直线段
	my($px,$py) = ($zx+150+$len_lsc,170+$y);
	$svg->path("d"=>"M $px 50 L $px $py",style=>"stroke","stroke"=>"#9D9D9D","stroke-width"=>0.1,fill=>"none");
	($px,$py) = ($zx+150+$len_lsc+$len_ir,170+$y);
	$svg->path("d"=>"M $px 50 L $px $py",style=>"stroke","stroke"=>"#9D9D9D","stroke-width"=>0.1,fill=>"none");
	($px,$py) = ($zx+150+$len_lsc+$len_ir+$len_ssc,170+$y);
	$svg->path("d"=>"M $px 50 L $px $py",style=>"stroke","stroke"=>"#9D9D9D","stroke-width"=>0.1,fill=>"none");
	($px,$py) = ($zx+150+$len_lsc+$len_ir+$len_ssc+$len_ir,170+$y);
	$svg->path("d"=>"M $px 50 L $px $py",style=>"stroke","stroke"=>"#9D9D9D","stroke-width"=>0.1,fill=>"none");

	#JLB附近的基因
	get_gene_of_junction($end_lsc,\%pos_for_gene,$zx+150+$len_lsc,80+$y);
	get_gene_of_junction($end_irb,\%pos_for_gene,$zx+150+$len_lsc+$len_ir,80+$y);
	get_gene_of_junction($end_ssc,\%pos_for_gene,$zx+150+$len_lsc+$len_ir+$len_ssc,80+$y);
	get_gene_of_junction($end_ira,\%pos_for_gene,$zx+150+$len_lsc+$len_ir+$len_ssc+$len_ir,80+$y,$end_ira);
	
	#序列长度(基因组长度)
	$svg->text(x=>$zx+130,y=>120+$y,'font-size'=>15,'font-family'=>'times new roman','-cdata'=>$seq_len." bp","text-anchor"=>"end");

	#五个矩形方框
	$svg->rect(x=>$zx+150,y=>80+$y,width=>$len_lsc,height=>20,"stroke-opacity"=>"1","stroke"=>"black",fill=>$color{$region_name[0]});	#LSC
	$svg->rect(x=>$zx+150+$len_lsc,y=>80+$y,width=>$len_ir,height=>20,"stroke-opacity"=>"1","stroke"=>"black",fill=>$color{$region_name[1]});	#IRb
	$svg->rect(x=>$zx+150+$len_lsc+$len_ir,y=>80+$y,width=>$len_ssc,height=>20,"stroke-opacity"=>"1","stroke"=>"black",fill=>$color{$region_name[2]});	#SSC
	$svg->rect(x=>$zx+150+$len_lsc+$len_ir+$len_ssc,y=>80+$y,width=>$len_ir,height=>20,"stroke-opacity"=>"1","stroke"=>"black",fill=>$color{$region_name[3]});	#IRa
	$svg->rect(x=>$zx+150+$len_lsc+$len_ir+$len_ssc+$len_ir,y=>80+$y,width=>$len_lsc,height=>20,"stroke-opacity"=>"1","stroke"=>"black",fill=>$color{$region_name[4]});	#LSC

	#区域名称
	$svg->text(x=>$zx+170,y=>97.5+$y,'font-size'=>15,'font-family'=>'times new roman','-cdata'=>$region_name[0],"text-anchor"=>"middle");
	$svg->text(x=>$zx+170+$len_lsc,y=>97.5+$y,'font-size'=>15,'font-family'=>'times new roman','-cdata'=>$region_name[1],"text-anchor"=>"middle");
	$svg->text(x=>$zx+170+$len_lsc+$len_ir,y=>97.5+$y,'font-size'=>15,'font-family'=>'times new roman','-cdata'=>$region_name[2],"text-anchor"=>"middle");
	$svg->text(x=>$zx+170+$len_lsc+$len_ir+$len_ssc,y=>97.5+$y,'font-size'=>15,'font-family'=>'times new roman','-cdata'=>$region_name[3],"text-anchor"=>"middle");
	$svg->text(x=>$zx+170+$len_lsc+$len_ir+$len_ssc+$len_ir,y=>97.5+$y,'font-size'=>15,'font-family'=>'times new roman','-cdata'=>$region_name[4],"text-anchor"=>"middle");

	#区域长度
	$svg->text(x=>$zx+150+$len_lsc-10,y=>97.5+$y,'font-size'=>15,'font-family'=>'times new roman','-cdata'=>$region_len[0]." bp","text-anchor"=>"end");
	$svg->text(x=>$zx+150+$len_lsc+$len_ir-10,y=>97.5+$y,'font-size'=>15,'font-family'=>'times new roman','-cdata'=>$region_len[1]." bp","text-anchor"=>"end");
	$svg->text(x=>$zx+150+$len_lsc+$len_ir+$len_ssc-10,y=>97.5+$y,'font-size'=>15,'font-family'=>'times new roman','-cdata'=>$region_len[2]." bp","text-anchor"=>"end");
	$svg->text(x=>$zx+150+$len_lsc+$len_ir+$len_ssc+$len_ir-10,y=>97.5+$y,'font-size'=>15,'font-family'=>'times new roman','-cdata'=>$region_len[3]." bp","text-anchor"=>"end");

	#线段上方文字
	$svg->text(x=>$zx+150+$len_lsc,y=>35,'font-size'=>15,'font-family'=>'times new roman','-cdata'=>"JLB","text-anchor"=>"middle");
	$svg->text(x=>$zx+150+$len_lsc+$len_ir,y=>35,'font-size'=>15,'font-family'=>'times new roman','-cdata'=>"JSB","text-anchor"=>"middle");
	$svg->text(x=>$zx+150+$len_lsc+$len_ir+$len_ssc,y=>35,'font-size'=>15,'font-family'=>'times new roman','-cdata'=>"JSA","text-anchor"=>"middle");
	$svg->text(x=>$zx+150+$len_lsc+$len_ir+$len_ssc+$len_ir,y=>35,'font-size'=>15,'font-family'=>'times new roman','-cdata'=>"JLA","text-anchor"=>"middle");

}

print  OUT $svg->xmlify;

#转化图片，需要用到convert命令

print"converting...\n";
my $dirname1 = dirname $outfile;
my $basename1 = basename $outfile;
chomp(my $pwd = `pwd`);
chdir "$dirname1";
`convert -density 300 $basename1 $pre.png`;
`convert -density 300 $pre.png $pre.tif`;
`svg2xxx -t pdf $basename1`;
chdir "$pwd";


sub get_gene_of_junction{
	my $junction_pos = shift;
	my $gene_pos = shift;
	my $x = shift;
	my $y = shift;

	print"\nJunction_pos:$junction_pos\n";
	#有基因跨区域
	if(defined $gene_pos->{$junction_pos} and defined $gene_pos->{$junction_pos+1} and defined $gene_pos->{$junction_pos-1}){
		print"    Gene over junction\n";
		my $s = $gene_pos->{$junction_pos}->[2];
		my $e = $gene_pos->{$junction_pos}->[3];

		print"\tOver: $gene_pos->{$junction_pos}->[0] $gene_pos->{$junction_pos}->[2] $gene_pos->{$junction_pos}->[3] $gene_pos->{$junction_pos}->[1] ".($e - $s + 1)."\n";

		if($e - $junction_pos >= $junction_pos - $s){	#基因靠近右边（左边跨区域较短）
			#交界左边的基因
			print"\tRight\n";
			my $pos1;
			unless(1){	#如果跨区域超过100 bp则不进行寻找下一个基因
				1;
			}else{	#向左寻找
				$pos1 = $gene_pos->{$junction_pos}->[2]-1;
				while(!defined $gene_pos->{$pos1} or !defined $gene_pos->{$pos1-1}){
					$pos1--;
				}

				print"\tleft: $gene_pos->{$pos1}->[0] $gene_pos->{$pos1}->[2] $gene_pos->{$pos1}->[3] $gene_pos->{$pos1}->[1] ".($gene_pos->{$pos1}->[3]-$gene_pos->{$pos1}->[2]+1)."\n";
				if($gene_pos->{$pos1}->[3] > $junction_pos){	#有基因重叠
					print"##gene is overlap\n";
					my $s = $gene_pos->{$pos1}->[2];
					my $e = $gene_pos->{$pos1}->[3];
					my $len2 = get_gene_len($e - $s);
					my $over_lap2 = get_gap_len($e - $junction_pos);
					my $x_for_gene2 = $x - $len2 + $over_lap2;

					#左边基因一律放到下方
#					my $y_for_gene2 = $gene_pos->{$pos1}->[1] ? $y + 20 : $y - 15;
					my $y_for_gene2 = $y + 20;

					#基因矩形，加上方向
#					$svg->rect(x=>$x_for_gene2,y=>$y_for_gene2,width=>$len2,height=>15,"stroke-opacity"=>"1","stroke"=>$color{$gene_pos->{$pos1}->[0]}||="red"),fill=>$color{$gene_pos->{$pos1}->[0]}||="red");
					if($gene_pos->{$pos1}->[1]){#负链
						$svg->rect(x=>$x_for_gene2+5,y=>$y_for_gene2,width=>$len2-5,height=>15,"stroke-opacity"=>"1","stroke"=>$color{$gene_pos->{$pos1}->[0]}||="red",fill=>$color{$gene_pos->{$pos1}->[0]}||="red");
						my $p0 = $x_for_gene2;#顶点x 
						my $p1 = $x_for_gene2+5;#底边x
						my $p2 = $y_for_gene2 + 7.5;#顶点y
						my $p3 = $y_for_gene2 + 15;#底y
						
						$svg->path("d"=>"M $p0 $p2 L $p1 $y_for_gene2 L $p1 $p3",style=>"stroke","stroke"=>$color{$gene_pos->{$pos1}->[0]}||="red","stroke-width"=>0.5,fill=>$color{$gene_pos->{$pos1}->[0]}||="red");  #箭头
					}else{#正链
						$svg->rect(x=>$x_for_gene2,y=>$y_for_gene2,width=>$len2-5,height=>15,"stroke-opacity"=>"1","stroke"=>$color{$gene_pos->{$pos1}->[0]}||="red",fill=>$color{$gene_pos->{$pos1}->[0]}||="red");
						my $p0 = $x_for_gene2+$len2;      
						my $p1 = $x_for_gene2+$len2-5;    
						my $p2 = $y_for_gene2 + 7.5;
						my $p3 = $y_for_gene2 + 15; 
						$svg->path("d"=>"M $p0 $p2 L $p1 $y_for_gene2 L $p1 $p3",style=>"stroke","stroke"=>$color{$gene_pos->{$pos1}->[0]}||="red","stroke-width"=>0.5,fill=>$color{$gene_pos->{$pos1}->[0]}||="red");  #箭头
					}

					

					my ($line2_x11,$line2_x21) = ($x_for_gene2,$x - 2);
					my ($line2_x12,$line2_x22) = ($x+2,$x_for_gene2 + $len2);
#					my $line2_y = $gene_pos->{$pos1}->[1] ? $y + 20 + 15 + 10:$y - 15 -10;
					my $line2_y = $y + 20 + 15 + 10;

					#跨区域的两条平行线
					$svg->path("d"=>"M $line2_x11 $line2_y L $line2_x21 $line2_y",style=>"stroke","stroke"=>"black","stroke-width"=>0.1,fill=>"none");	
					$svg->path("d"=>"M $line2_x12 $line2_y L $line2_x22 $line2_y",style=>"stroke","stroke"=>"black","stroke-width"=>0.1,fill=>"none");

					my ($line2_vx1,$line2_vx2,$line2_vx3,$line2_vx4) = ($line2_x11,$line2_x21,$line2_x12,$line2_x22);
					my ($line2_vy1,$line2_vy2,$line2_vy3,$line2_vy4,$line2_vy5,$line2_vy6,$line2_vy7,$line2_vy8) = ($line2_y+5,$line2_y-5,$line2_y+5,$line2_y-5,$line2_y+5,$line2_y-5,$line2_y+5,$line2_y-5);

					#垂直的小短线
					$svg->path("d"=>"M $line2_vx1 $line2_vy1 L $line2_vx1 $line2_vy2",style=>"stroke","stroke"=>"black","stroke-width"=>0.1,fill=>"none");	
					$svg->path("d"=>"M $line2_vx2 $line2_vy3 L $line2_vx2 $line2_vy4",style=>"stroke","stroke"=>"black","stroke-width"=>0.1,fill=>"none");
					$svg->path("d"=>"M $line2_vx3 $line2_vy5 L $line2_vx3 $line2_vy6",style=>"stroke","stroke"=>"black","stroke-width"=>0.1,fill=>"none");
					$svg->path("d"=>"M $line2_vx4 $line2_vy7 L $line2_vx4 $line2_vy8",style=>"stroke","stroke"=>"black","stroke-width"=>0.1,fill=>"none");

					#基因名字
					my $name_x1 = ($line2_x11+$line2_x22)/2;
#					my $name_y1 = $gene_pos->{$pos1}->[1] ? $y + 20 + 15 - 2 : $y - 2;
					my $name_y1 = $y + 20 + 15 - 2;
					$svg->text(x=>$name_x1,y=>$name_y1,'font-size'=>13,'font-family'=>'times new roman','-cdata'=>"$gene_pos->{$pos1}->[0]","text-anchor"=>"middle",fill=>"white",'font-style'=>'italic');

					#长度
#					my $len_y = $gene_pos->{$pos1}->[1] ? $line2_y + 5 + 10 : $line2_y - 5 ;
					my $len_y = $line2_y + 5 + 10;
					$svg->text(x=>($line2_vx1 + $line2_vx2)/2,y=>$len_y,'font-size'=>13,'font-family'=>'times new roman','-cdata'=>($junction_pos-$gene_pos->{$pos1}->[2]+1)." bp","text-anchor"=>"middle",fill=>"black");
					$svg->text(x=>$line2_vx4 + 2,y=>$len_y,'font-size'=>13,'font-family'=>'times new roman','-cdata'=>($gene_pos->{$pos1}->[3] - $junction_pos)." bp","text-anchor"=>"start",fill=>"black");
				
				}else{
					my $len1 = get_gene_len(($gene_pos->{$pos1}->[3] - $gene_pos->{$pos1}->[2]));
					my $gap1 = get_gap_len($junction_pos - $gene_pos->{$pos1}->[3]);
					#$gene_pos->{$pos1}->[1] = !$gene_pos->{$junction_pos}->[1]?0:1;

					my $x_for_gene1 = $x - $gap1 - $len1;
#					my $y_for_gene1 = $gene_pos->{$pos1}->[1] ? $y + 20 : $y - 15;
					my $y_for_gene1 = $y + 20;

					#基因矩形
#					$svg->rect(x=>$x_for_gene1,y=>$y_for_gene1,width=>$len1,height=>15,"stroke-opacity"=>"1","stroke"=>"black",fill=>$color{$gene_pos->{$pos1}->[0]}||="red");
					#基因矩形，加上方向
					if($gene_pos->{$pos1}->[1]){#负链
						$svg->rect(x=>$x_for_gene1+5,y=>$y_for_gene1,width=>$len1-5,height=>15,"stroke-opacity"=>"1","stroke"=>$color{$gene_pos->{$pos1}->[0]}||="red",fill=>$color{$gene_pos->{$pos1}->[0]}||="red");
						my $p0 = $x_for_gene1;#顶点x 
						my $p1 = $x_for_gene1+5;#底边x
						my $p2 = $y_for_gene1 + 7.5;#顶点y
						my $p3 = $y_for_gene1 + 15;#底y
						
						$svg->path("d"=>"M $p0 $p2 L $p1 $y_for_gene1 L $p1 $p3",style=>"stroke","stroke"=>$color{$gene_pos->{$pos1}->[0]}||="red","stroke-width"=>0.5,fill=>$color{$gene_pos->{$pos1}->[0]}||="red");  #箭头
					}else{#正链
						$svg->rect(x=>$x_for_gene1,y=>$y_for_gene1,width=>$len1-5,height=>15,"stroke-opacity"=>"1","stroke"=>$color{$gene_pos->{$pos1}->[0]}||="red",fill=>$color{$gene_pos->{$pos1}->[0]}||="red");
						my $p0 = $x_for_gene1+$len1;      
						my $p1 = $x_for_gene1+$len1-5;    
						my $p2 = $y_for_gene1 + 7.5;
						my $p3 = $y_for_gene1 + 15; 
						$svg->path("d"=>"M $p0 $p2 L $p1 $y_for_gene1 L $p1 $p3",style=>"stroke","stroke"=>$color{$gene_pos->{$pos1}->[0]}||="red","stroke-width"=>0.5,fill=>$color{$gene_pos->{$pos1}->[0]}||="red");  #箭头
					}
					#基因名字
					my $name_x1;
					my $name_y1;
					my $bp_y1;
					
					if($len1 < 45){
						$name_x1 = ($x_for_gene1 + $x_for_gene1 + $len1)/2;
#						$name_y1 = 	$gene_pos->{$pos1}->[1] ? $y + 20 + 15 + 15 : $y - 15 - 2;
						$name_y1 = $y + 20 + 15 + 15;
#						$bp_y1 = $gene_pos->{$pos1}->[1] ? $y + 20 + 15 - 2 : $y - 2;
						$bp_y1 = $y + 20 + 15 - 2;

						$svg->text(x=>$name_x1,y=>$name_y1,'font-size'=>13,'font-family'=>'times new roman','-cdata'=>"$gene_pos->{$pos1}->[0]","text-anchor"=>"end",fill=>"black",'font-style'=>'italic');
						$svg->text(x=>$x_for_gene1-2,y=>$bp_y1,'font-size'=>13,'font-family'=>'times new roman','-cdata'=>($gene_pos->{$pos1}->[3] - $gene_pos->{$pos1}->[2]+1)." bp","text-anchor"=>"end",fill=>"black",'font-style'=>'italic') if($is_print_gene_len);
					}else{
						$name_x1 = ($x_for_gene1 + $x_for_gene1 + $len1)/2;
#						$name_y1 = 	$gene_pos->{$pos1}->[1] ? $y + 20 + 15 - 2 : $y - 2;
						$name_y1 = $y + 20 + 15 - 2;
#						$bp_y1 = $gene_pos->{$pos1}->[1] ? $y + 20 + 15 + 14 : $y - 15 - 1;
						$bp_y1 = $y + 20 + 15 + 14;

						$svg->text(x=>$name_x1,y=>$name_y1,'font-size'=>13,'font-family'=>'times new roman','-cdata'=>"$gene_pos->{$pos1}->[0]","text-anchor"=>"middle",fill=>"white",'font-style'=>'italic');
						$svg->text(x=>$name_x1,y=>$bp_y1,'font-size'=>13,'font-family'=>'times new roman','-cdata'=>($gene_pos->{$pos1}->[3] - $gene_pos->{$pos1}->[2]+1)." bp","text-anchor"=>"middle",fill=>"black");
					}
					
					{	#gap距离 1
						use Math::Trig;
						my $x0 = $x - $gap1/2 - 25;	#圆心
#						my $y0 = $gene_pos->{$pos1}->[1] ? $y + 20 + 15 : $y - 15;
						my $y0 = $y + 20 + 15;

						my $r = $x - $x0 - $gap1/2;	#半径

						my $x1 = $x0 + $r * cos(0);			#外圆弧的起点和中点
						my $y1 = $y0 + $r * sin(0);

						if(1 or $gene_pos->{$pos1}->[1]){#确定在下方
							my $x2 = $x0 + $r * cos(pi*0.5);
							my $y2 = $y0 + $r * sin(pi*0.5);

							$svg->path("d"=>"M $x1  $y1 A $r $r 0 0 1 $x2  $y2",style=>"stroke", "stroke"=>"black","stroke-width"=>0.5, fill=>"none");	#圆弧

							my ($x21,$y21) = ($x2+3,$y2-3);	#箭头
							my ($x22,$y22) = ($x2+3,$y2+3);

							$svg->path("d"=>"M $x2  $y2 L $x21 $y21",style=>"stroke","stroke"=>"black","stroke-width"=>0.5,fill=>"none");
							$svg->path("d"=>"M $x2  $y2 L $x22 $y22",style=>"stroke","stroke"=>"black","stroke-width"=>0.5,fill=>"none");
							$svg->text(x=>$x0-2,y=>$y0+$r+5,'font-size'=>13,'font-family'=>'times new roman','-cdata'=>$junction_pos - $gene_pos->{$pos1}->[3]." bp","text-anchor"=>"end",fill=>"black"); #gap
						}else{
							my $x2 = $x0 + $r * cos(-pi*0.5);
							my $y2 = $y0 + $r * sin(-pi*0.5);
							$svg->path("d"=>"M $x1  $y1 A $r $r 0 0 0 $x2  $y2",style=>"stroke", "stroke"=>"black","stroke-width"=>0.5, fill=>"none");	#圆弧

							my ($x21,$y21) = ($x2+3,$y2-3);	#箭头
							my ($x22,$y22) = ($x2+3,$y2+3);

							$svg->path("d"=>"M $x2  $y2 L $x21 $y21",style=>"stroke","stroke"=>"black","stroke-width"=>0.5,fill=>"none");
							$svg->path("d"=>"M $x2  $y2 L $x22 $y22",style=>"stroke","stroke"=>"black","stroke-width"=>0.5,fill=>"none");
							$svg->text(x=>$x0-2,y=>$y0-$r+5,'font-size'=>13,'font-family'=>'times new roman','-cdata'=>$junction_pos - $gene_pos->{$pos1}->[3]." bp","text-anchor"=>"end",fill=>"black"); #gap
						}
					}
				}
			}

			#画跨区域的偏右方的基因，放到上方
			my $len1 = get_gene_len($e - $s) > 20 ? get_gene_len($e - $s) : 45;
			my $over_lap1 = get_gap_len($junction_pos - $s);
			my $x_for_gene1 = $x - $over_lap1;
#			my $y_for_gene1 = !$gene_pos->{$pos1}->[1] ? $y + 20 : $y - 15;
			my $y_for_gene1 =  $y - 15;
			
			#基因矩形

#			$svg->rect(x=>$x_for_gene1,y=>$y_for_gene1,width=>$len1,height=>15,"stroke-opacity"=>"1","stroke"=>"black",fill=>$color{$gene_pos->{$junction_pos}->[0]}||="red");	
			if($gene_pos->{$junction_pos}->[1]){#负链
				$svg->rect(x=>$x_for_gene1+5,y=>$y_for_gene1,width=>$len1-5,height=>15,"stroke-opacity"=>"1","stroke"=>$color{$gene_pos->{$junction_pos}->[0]}||="red",fill=>$color{$gene_pos->{$junction_pos}->[0]}||="red");
				my $p0 = $x_for_gene1;#顶点x 
				my $p1 = $x_for_gene1+5;#底边x
				my $p2 = $y_for_gene1 + 7.5;#顶点y
				my $p3 = $y_for_gene1 + 15;#底y
				
				$svg->path("d"=>"M $p0 $p2 L $p1 $y_for_gene1 L $p1 $p3",style=>"stroke","stroke"=>$color{$gene_pos->{$junction_pos}->[0]}||="red","stroke-width"=>0.5,fill=>$color{$gene_pos->{$junction_pos}->[0]}||="red");  #箭头
			}else{#正链
				$svg->rect(x=>$x_for_gene1,y=>$y_for_gene1,width=>$len1-5,height=>15,"stroke-opacity"=>"1","stroke"=>$color{$gene_pos->{$junction_pos}->[0]}||="red",fill=>$color{$gene_pos->{$junction_pos}->[0]}||="red");
				my $p0 = $x_for_gene1+$len1;      
				my $p1 = $x_for_gene1+$len1-5;    
				my $p2 = $y_for_gene1 + 7.5;
				my $p3 = $y_for_gene1 + 15; 
				$svg->path("d"=>"M $p0 $p2 L $p1 $y_for_gene1 L $p1 $p3",style=>"stroke","stroke"=>$color{$gene_pos->{$junction_pos}->[0]}||="red","stroke-width"=>0.5,fill=>$color{$gene_pos->{$junction_pos}->[0]}||="red");  #箭头
			}

			my ($line1_x11,$line1_x21) = ($x_for_gene1,$x - 2);
			my ($line1_x12,$line1_x22) = ($x+2,$x_for_gene1 + $len1);
#			my $line1_y = !$gene_pos->{$pos1}->[1] ? $y + 20 + 15 + 10:$y - 15 -10;
			my $line1_y = $y - 15 -10;

			#跨区域的两条平行线
			$svg->path("d"=>"M $line1_x11 $line1_y L $line1_x21 $line1_y",style=>"stroke","stroke"=>"black","stroke-width"=>0.1,fill=>"none");	
			$svg->path("d"=>"M $line1_x12 $line1_y L $line1_x22 $line1_y",style=>"stroke","stroke"=>"black","stroke-width"=>0.1,fill=>"none");

			my ($line1_vx1,$line1_vx2,$line1_vx3,$line1_vx4) = ($line1_x11,$line1_x21,$line1_x12,$line1_x22);
			my ($line1_vy1,$line1_vy2,$line1_vy3,$line1_vy4,$line1_vy5,$line1_vy6,$line1_vy7,$line1_vy8) = ($line1_y+5,$line1_y-5,$line1_y+5,$line1_y-5,$line1_y+5,$line1_y-5,$line1_y+5,$line1_y-5);
			
			#垂直的小短线
			$svg->path("d"=>"M $line1_vx1 $line1_vy1 L $line1_vx1 $line1_vy2",style=>"stroke","stroke"=>"black","stroke-width"=>0.1,fill=>"none");	
			$svg->path("d"=>"M $line1_vx2 $line1_vy3 L $line1_vx2 $line1_vy4",style=>"stroke","stroke"=>"black","stroke-width"=>0.1,fill=>"none");
			$svg->path("d"=>"M $line1_vx3 $line1_vy5 L $line1_vx3 $line1_vy6",style=>"stroke","stroke"=>"black","stroke-width"=>0.1,fill=>"none");
			$svg->path("d"=>"M $line1_vx4 $line1_vy7 L $line1_vx4 $line1_vy8",style=>"stroke","stroke"=>"black","stroke-width"=>0.1,fill=>"none");

			#基因名字
			my $name_x1 = ($line1_x11+$line1_x22)/2;
#			my $name_y1 = !$gene_pos->{$pos1}->[1] ? $y + 20 + 15 - 2 : $y - 2;
			my $name_y1 = $y - 2;
			$svg->text(x=>$name_x1,y=>$name_y1,'font-size'=>13,'font-family'=>'times new roman','-cdata'=>"$gene_pos->{$junction_pos}->[0]","text-anchor"=>"middle",fill=>"white",'font-style'=>'italic');

			#长度
#			my $len_y = !$gene_pos->{$pos1}->[1] ? $line1_y + 5 + 10 : $line1_y - 5 ;
			my $len_y =  $line1_y - 5 ;
			$svg->text(x=>$line1_vx1-2,y=>$len_y,'font-size'=>13,'font-family'=>'times new roman','-cdata'=>($junction_pos-$gene_pos->{$junction_pos}->[2]+1)." bp","text-anchor"=>"end",fill=>"black");
			$svg->text(x=>($line1_vx3+$line1_vx4)/2,y=>$len_y,'font-size'=>13,'font-family'=>'times new roman','-cdata'=>($gene_pos->{$junction_pos}->[3] - $junction_pos)." bp","text-anchor"=>"middle",fill=>"black");

		}else{
			###############################
			##基因跨区域，且靠近左方
			###############################

			print"\tLeft\n";
			my $len2 = get_gene_len($e - $s);
			my $over_lap2 = get_gap_len($e - $junction_pos);
			my $x_for_gene2 = $x - $len2 + $over_lap2;
#左边的基因放到下面
#			my $y_for_gene2 = $gene_pos->{$junction_pos}->[1] ? $y + 20 : $y - 15;
			my $y_for_gene2 = $y + 20;

#			$svg->rect(x=>$x_for_gene2,y=>$y_for_gene2,width=>$len2,height=>15,"stroke-opacity"=>"1","stroke"=>"black",fill=>$color{$gene_pos->{$junction_pos}->[0]}||="red");
			if($gene_pos->{$junction_pos}->[1]){#负链
				$svg->rect(x=>$x_for_gene2+5,y=>$y_for_gene2,width=>$len2-5,height=>15,"stroke-opacity"=>"1","stroke"=>$color{$gene_pos->{$junction_pos}->[0]}||="red",fill=>$color{$gene_pos->{$junction_pos}->[0]}||="red");
				my $p0 = $x_for_gene2;#顶点x 
				my $p1 = $x_for_gene2+5;#底边x
				my $p2 = $y_for_gene2 + 7.5;#顶点y
				my $p3 = $y_for_gene2 + 15;#底y
				
				$svg->path("d"=>"M $p0 $p2 L $p1 $y_for_gene2 L $p1 $p3",style=>"stroke","stroke"=>$color{$gene_pos->{$junction_pos}->[0]}||="red","stroke-width"=>0.5,fill=>$color{$gene_pos->{$junction_pos}->[0]}||="red");  #箭头
			}else{#正链
				$svg->rect(x=>$x_for_gene2,y=>$y_for_gene2,width=>$len2-5,height=>15,"stroke-opacity"=>"1","stroke"=>$color{$gene_pos->{$junction_pos}->[0]}||="red",fill=>$color{$gene_pos->{$junction_pos}->[0]}||="red");
				my $p0 = $x_for_gene2+$len2;      
				my $p1 = $x_for_gene2+$len2-5;    
				my $p2 = $y_for_gene2 + 7.5;
				my $p3 = $y_for_gene2 + 15; 
				$svg->path("d"=>"M $p0 $p2 L $p1 $y_for_gene2 L $p1 $p3",style=>"stroke","stroke"=>$color{$gene_pos->{$junction_pos}->[0]}||="red","stroke-width"=>0.5,fill=>$color{$gene_pos->{$junction_pos}->[0]}||="red");  #箭头
			}
			my ($line2_x11,$line2_x21) = ($x_for_gene2,$x - 2);
			my ($line2_x12,$line2_x22) = ($x+2,$x_for_gene2 + $len2);
#			my $line2_y = $gene_pos->{$junction_pos}->[1] ? $y + 20 + 15 + 10:$y - 15 -10;
			my $line2_y = $y + 20 + 15 + 10;

			#跨区域的两条平行线
			$svg->path("d"=>"M $line2_x11 $line2_y L $line2_x21 $line2_y",style=>"stroke","stroke"=>"black","stroke-width"=>0.1,fill=>"none");	
			$svg->path("d"=>"M $line2_x12 $line2_y L $line2_x22 $line2_y",style=>"stroke","stroke"=>"black","stroke-width"=>0.1,fill=>"none");

			my ($line2_vx1,$line2_vx2,$line2_vx3,$line2_vx4) = ($line2_x11,$line2_x21,$line2_x12,$line2_x22);
			my ($line2_vy1,$line2_vy2,$line2_vy3,$line2_vy4,$line2_vy5,$line2_vy6,$line2_vy7,$line2_vy8) = ($line2_y+5,$line2_y-5,$line2_y+5,$line2_y-5,$line2_y+5,$line2_y-5,$line2_y+5,$line2_y-5);

			#垂直的小短线
			$svg->path("d"=>"M $line2_vx1 $line2_vy1 L $line2_vx1 $line2_vy2",style=>"stroke","stroke"=>"black","stroke-width"=>0.1,fill=>"none");	
			$svg->path("d"=>"M $line2_vx2 $line2_vy3 L $line2_vx2 $line2_vy4",style=>"stroke","stroke"=>"black","stroke-width"=>0.1,fill=>"none");
			$svg->path("d"=>"M $line2_vx3 $line2_vy5 L $line2_vx3 $line2_vy6",style=>"stroke","stroke"=>"black","stroke-width"=>0.1,fill=>"none");
			$svg->path("d"=>"M $line2_vx4 $line2_vy7 L $line2_vx4 $line2_vy8",style=>"stroke","stroke"=>"black","stroke-width"=>0.1,fill=>"none");

			#基因名字
			my $name_x1 = ($line2_x11+$line2_x21)/2;
#			my $name_y1 = $gene_pos->{$junction_pos}->[1] ? $y + 20 + 15 - 2 : $y - 2;
			my $name_y1 = $y + 20 + 15 - 2;

			$svg->text(x=>$name_x1,y=>$name_y1,'font-size'=>13,'font-family'=>'times new roman','-cdata'=>"$gene_pos->{$junction_pos}->[0]","text-anchor"=>"middle",fill=>"white",'font-style'=>'italic');

			#长度
#			my $len_y = $gene_pos->{$junction_pos}->[1] ? $line2_y + 5 + 10 : $line2_y - 5 ;
			my $len_y = $line2_y + 5 + 10;

			$svg->text(x=>($line2_vx1 + $line2_vx2)/2,y=>$len_y,'font-size'=>13,'font-family'=>'times new roman','-cdata'=>($junction_pos-$gene_pos->{$junction_pos}->[2]+1)." bp","text-anchor"=>"middle",fill=>"black");
			$svg->text(x=>$line2_vx4 + 2,y=>$len_y,'font-size'=>13,'font-family'=>'times new roman','-cdata'=>($gene_pos->{$junction_pos}->[3] - $junction_pos)." bp","text-anchor"=>"start",fill=>"black");

			unless(1){	#
				1
			}else{	#向右寻找
				my $pos2 = $gene_pos->{$junction_pos}->[3]+1;
				while(!defined $gene_pos->{$pos2} or !defined $gene_pos->{$pos2+1}){
					$pos2++;
				}
				print"\tleft: $gene_pos->{$pos2}->[0] $gene_pos->{$pos2}->[2] $gene_pos->{$pos2}->[3] $gene_pos->{$pos2}->[1] ".($gene_pos->{$pos2}->[3]-$gene_pos->{$pos2}->[2]+1)."\n";

				if($gene_pos->{$pos2}->[2] <= $junction_pos){
					print"    !Some gene is special\n";
					my $s = $gene_pos->{$pos2}->[2];
					my $e = $gene_pos->{$pos2}->[3];

					my $len1 = get_gene_len($e - $s) > 20 ? get_gene_len($e - $s) : 45;
					my $over_lap1 = get_gap_len($junction_pos - $s);
					my $x_for_gene1 = $x - $over_lap1;
#右边的放到上面
#					my $y_for_gene1 = !$gene_pos->{$junction_pos}->[1] ? $y + 20 : $y - 15;
					my $y_for_gene1 = $y - 15;

					#基因矩形
#					$svg->rect(x=>$x_for_gene1,y=>$y_for_gene1,width=>$len1,height=>15,"stroke-opacity"=>"1","stroke"=>"black",fill=>$color{$gene_pos->{$pos2}->[0]}||="red");	
					if($gene_pos->{$pos2}->[1]){#负链
						$svg->rect(x=>$x_for_gene1+5,y=>$y_for_gene1,width=>$len1-5,height=>15,"stroke-opacity"=>"1","stroke"=>$color{$gene_pos->{$pos2}->[0]}||="red",fill=>$color{$gene_pos->{$pos2}->[0]}||="red");
						my $p0 = $x_for_gene1;#顶点x 
						my $p1 = $x_for_gene1+5;#底边x
						my $p2 = $y_for_gene1 + 7.5;#顶点y
						my $p3 = $y_for_gene1 + 15;#底y
						
						$svg->path("d"=>"M $p0 $p2 L $p1 $y_for_gene1 L $p1 $p3",style=>"stroke","stroke"=>$color{$gene_pos->{$pos2}->[0]}||="red","stroke-width"=>0.5,fill=>$color{$gene_pos->{$pos2}->[0]}||="red");  #箭头
					}else{#正链
						$svg->rect(x=>$x_for_gene1,y=>$y_for_gene1,width=>$len1-5,height=>15,"stroke-opacity"=>"1","stroke"=>$color{$gene_pos->{$pos2}->[0]}||="red",fill=>$color{$gene_pos->{$pos2}->[0]}||="red");
						my $p0 = $x_for_gene1+$len1;      
						my $p1 = $x_for_gene1+$len1-5;    
						my $p2 = $y_for_gene1 + 7.5;
						my $p3 = $y_for_gene1 + 15; 
						$svg->path("d"=>"M $p0 $p2 L $p1 $y_for_gene1 L $p1 $p3",style=>"stroke","stroke"=>$color{$gene_pos->{$pos2}->[0]}||="red","stroke-width"=>0.5,fill=>$color{$gene_pos->{$pos2}->[0]}||="red");  #箭头
					}
					my ($line1_x11,$line1_x21) = ($x_for_gene1,$x - 2);
					my ($line1_x12,$line1_x22) = ($x+2,$x_for_gene1 + $len1);

#					my $line1_y = !$gene_pos->{$junction_pos}->[1] ? $y + 20 + 15 + 10:$y - 15 -10;
					my $line1_y = $y - 15 -10;

					#跨区域的两条平行线
					$svg->path("d"=>"M $line1_x11 $line1_y L $line1_x21 $line1_y",style=>"stroke","stroke"=>"black","stroke-width"=>0.1,fill=>"none");	
					$svg->path("d"=>"M $line1_x12 $line1_y L $line1_x22 $line1_y",style=>"stroke","stroke"=>"black","stroke-width"=>0.1,fill=>"none");

					my ($line1_vx1,$line1_vx2,$line1_vx3,$line1_vx4) = ($line1_x11,$line1_x21,$line1_x12,$line1_x22);
					my ($line1_vy1,$line1_vy2,$line1_vy3,$line1_vy4,$line1_vy5,$line1_vy6,$line1_vy7,$line1_vy8) = ($line1_y+5,$line1_y-5,$line1_y+5,$line1_y-5,$line1_y+5,$line1_y-5,$line1_y+5,$line1_y-5);
					
					#垂直的小短线
					$svg->path("d"=>"M $line1_vx1 $line1_vy1 L $line1_vx1 $line1_vy2",style=>"stroke","stroke"=>"black","stroke-width"=>0.1,fill=>"none");	
					$svg->path("d"=>"M $line1_vx2 $line1_vy3 L $line1_vx2 $line1_vy4",style=>"stroke","stroke"=>"black","stroke-width"=>0.1,fill=>"none");
					$svg->path("d"=>"M $line1_vx3 $line1_vy5 L $line1_vx3 $line1_vy6",style=>"stroke","stroke"=>"black","stroke-width"=>0.1,fill=>"none");
					$svg->path("d"=>"M $line1_vx4 $line1_vy7 L $line1_vx4 $line1_vy8",style=>"stroke","stroke"=>"black","stroke-width"=>0.1,fill=>"none");

					#基因名字
					my $name_x1 = ($line1_x12+$line1_x22)/2;
#					my $name_y1 = !$gene_pos->{$junction_pos}->[1] ? $y + 20 + 15 - 2 : $y - 2;
					my $name_y1 =  $y - 2;
					$svg->text(x=>$name_x1,y=>$name_y1,'font-size'=>13,'font-family'=>'times new roman','-cdata'=>"$gene_pos->{$pos2}->[0]","text-anchor"=>"middle",fill=>"white",'font-style'=>'italic');

					#长度
#					my $len_y = !$gene_pos->{$junction_pos}->[1] ? $line1_y + 5 + 10 : $line1_y - 5 ;
					my $len_y = $line1_y - 5 ;

					$svg->text(x=>$line1_vx1-2,y=>$len_y,'font-size'=>13,'font-family'=>'times new roman','-cdata'=>($junction_pos-$gene_pos->{$pos2}->[2]+1)." bp","text-anchor"=>"end",fill=>"black");
					$svg->text(x=>($line1_vx3+$line1_vx4)/2,y=>$len_y,'font-size'=>13,'font-family'=>'times new roman','-cdata'=>($gene_pos->{$pos2}->[3] - $junction_pos)." bp","text-anchor"=>"middle",fill=>"black");
				}else{
					my $len2 = get_gene_len(($gene_pos->{$pos2}->[3] - $gene_pos->{$pos2}->[2]));
					my $gap2 = get_gap_len($gene_pos->{$pos2}->[2] - $junction_pos);

#					$gene_pos->{$pos2}->[1] = $gene_pos->{$junction_pos}->[1] ? 0:1;

					my $x_for_gene2 = $x + $gap2;
#					my $y_for_gene2 = $gene_pos->{$pos2}->[1]?$y + 20 : $y - 15;
					my $y_for_gene2 = $y - 15;
					#基因矩形
#					$svg->rect(x=>$x_for_gene2,y=>$y_for_gene2,width=>$len2,height=>15,"stroke-opacity"=>"1","stroke"=>"black",fill=>$color{$gene_pos->{$pos2}->[0]}||="red");
					if($gene_pos->{$pos2}->[1]){#负链
						$svg->rect(x=>$x_for_gene2+5,y=>$y_for_gene2,width=>$len2-5,height=>15,"stroke-opacity"=>"1","stroke"=>$color{$gene_pos->{$pos2}->[0]}||="red",fill=>$color{$gene_pos->{$pos2}->[0]}||="red");
						my $p0 = $x_for_gene2;#顶点x 
						my $p1 = $x_for_gene2+5;#底边x
						my $p2 = $y_for_gene2 + 7.5;#顶点y
						my $p3 = $y_for_gene2 + 15;#底y
						
						$svg->path("d"=>"M $p0 $p2 L $p1 $y_for_gene2 L $p1 $p3",style=>"stroke","stroke"=>$color{$gene_pos->{$pos2}->[0]}||="red","stroke-width"=>0.5,fill=>$color{$gene_pos->{$pos2}->[0]}||="red");  #箭头
					}else{#正链
						$svg->rect(x=>$x_for_gene2,y=>$y_for_gene2,width=>$len2-5,height=>15,"stroke-opacity"=>"1","stroke"=>$color{$gene_pos->{$pos2}->[0]}||="red",fill=>$color{$gene_pos->{$pos2}->[0]}||="red");
						my $p0 = $x_for_gene2+$len2;      
						my $p1 = $x_for_gene2+$len2-5;    
						my $p2 = $y_for_gene2 + 7.5;
						my $p3 = $y_for_gene2 + 15; 
						$svg->path("d"=>"M $p0 $p2 L $p1 $y_for_gene2 L $p1 $p3",style=>"stroke","stroke"=>$color{$gene_pos->{$pos2}->[0]}||="red","stroke-width"=>0.5,fill=>$color{$gene_pos->{$pos2}->[0]}||="red");  #箭头
					}

					#基因名字
					my $name_x2;
					my $name_y2;
					my $bp_y2;
					if($len2 < 45){
						$name_x2 = ($x_for_gene2 + $x_for_gene2 + $len2)/2;
#						$name_y2 = 	$gene_pos->{$pos2}->[1] ? $y + 20 + 15 + 15 : $y - 15 - 2;
						$name_y2 = $y - 15 - 2;
#						$bp_y2 = $gene_pos->{$pos2}->[1] ? $y + 20 + 15 - 2 : $y - 2;
						$bp_y2 = $y - 2;

						$svg->text(x=>$name_x2,y=>$name_y2,'font-size'=>13,'font-family'=>'times new roman','-cdata'=>"$gene_pos->{$pos2}->[0]","text-anchor"=>"start",fill=>"black",'font-style'=>'italic');
						$svg->text(x=>$x_for_gene2 + 2 + $len2,y=>$bp_y2,'font-size'=>13,'font-family'=>'times new roman','-cdata'=>($gene_pos->{$pos2}->[3] - $gene_pos->{$pos2}->[2]+1)." bp","text-anchor"=>"start",fill=>"black",'font-style'=>'italic') if($is_print_gene_len);
					}else{
						$name_x2 = ($x_for_gene2 + $x_for_gene2 + $len2)/2;
#						$name_y2 = 	$gene_pos->{$pos2}->[1] ? $y + 20 + 15 - 2 : $y - 2;
						$name_y2 = $y - 2;
#						$bp_y2 = $gene_pos->{$pos2}->[1] ? $y + 20 + 15 + 14 : $y - 15 - 1;
						$bp_y2 = $y - 15 - 1;

						$svg->text(x=>$name_x2,y=>$name_y2,'font-size'=>13,'font-family'=>'times new roman','-cdata'=>"$gene_pos->{$pos2}->[0]","text-anchor"=>"middle",fill=>"white",'font-style'=>'italic');
						$svg->text(x=>$name_x2,y=>$bp_y2,'font-size'=>13,'font-family'=>'times new roman','-cdata'=>($gene_pos->{$pos2}->[3] - $gene_pos->{$pos2}->[2]+1)." bp","text-anchor"=>"middle",fill=>"black");
					}

					{	#gap距离 2
						use Math::Trig;
						my $x0 = $x + $gap2/2 + 25;	#圆心
#						my $y0 = $gene_pos->{$pos2}->[1] ? $y + 20 + 15 : $y - 15;
						my $y0 = $y - 15;
						my $r = $x0 - $x - $gap2/2;	#半径

						my $x1 = $x0 + $r * cos(pi);			#外圆弧的起点和中点
						my $y1 = $y0 + $r * sin(pi);

						if(0 and $gene_pos->{$pos2}->[1]){
							my $x2 = $x0 + $r * cos(pi*0.5);
							my $y2 = $y0 + $r * sin(pi*0.5);

							$svg->path("d"=>"M $x1  $y1 A $r $r 0 0 0 $x2  $y2",style=>"stroke", "stroke"=>"black","stroke-width"=>0.5, fill=>"none");	#圆弧

							my ($x21,$y21) = ($x2-3,$y2-3);	#箭头
							my ($x22,$y22) = ($x2-3,$y2+3);

							$svg->path("d"=>"M $x2  $y2 L $x21 $y21",style=>"stroke","stroke"=>"black","stroke-width"=>0.5,fill=>"none");
							$svg->path("d"=>"M $x2  $y2 L $x22 $y22",style=>"stroke","stroke"=>"black","stroke-width"=>0.5,fill=>"none");

							$svg->text(x=>$x0+2,y=>$y0+$r+5,'font-size'=>13,'font-family'=>'times new roman','-cdata'=>($gene_pos->{$pos2}->[2] - $junction_pos - 1)." bp","text-anchor"=>"start",fill=>"black"); #gap
						}else{
							my $x2 = $x0 + $r * cos(-pi*0.5);
							my $y2 = $y0 + $r * sin(-pi*0.5);
							$svg->path("d"=>"M $x1  $y1 A $r $r 0 0 1 $x2  $y2",style=>"stroke", "stroke"=>"black","stroke-width"=>0.5, fill=>"none");	#圆弧

							my ($x21,$y21) = ($x2-3,$y2-3);	#箭头
							my ($x22,$y22) = ($x2-3,$y2+3);

							$svg->path("d"=>"M $x2  $y2 L $x21 $y21",style=>"stroke","stroke"=>"black","stroke-width"=>0.5,fill=>"none");
							$svg->path("d"=>"M $x2  $y2 L $x22 $y22",style=>"stroke","stroke"=>"black","stroke-width"=>0.5,fill=>"none");
							$svg->text(x=>$x0+2,y=>$y0-$r+5,'font-size'=>13,'font-family'=>'times new roman','-cdata'=>($gene_pos->{$pos2}->[2] - $junction_pos - 1)." bp","text-anchor"=>"start",fill=>"black"); #gap
						}
					}
				}
			}
		}

	}else{
		print"    No gene in junction_pos:\n";

		#向左寻找
		my $pos1 = $junction_pos;
		while(!defined $gene_pos->{$pos1} or !defined $gene_pos->{$pos1-1}){
			$pos1--;
		}
		print"\tleft : $gene_pos->{$pos1}->[0] $gene_pos->{$pos1}->[2] $gene_pos->{$pos1}->[3] $gene_pos->{$pos1}->[1] ".($gene_pos->{$pos1}->[3] - $gene_pos->{$pos1}->[2] + 1)."\n";
		my $len1 = get_gene_len(($gene_pos->{$pos1}->[3] - $gene_pos->{$pos1}->[2]));
		my $gap1 = get_gap_len($junction_pos - $gene_pos->{$pos1}->[3]);
		
		my $x_for_gene1 = $x - $gap1 - $len1;
#左边放到下方
#		my $y_for_gene1 = $gene_pos->{$pos1}->[1] ? $y + 20 : $y - 15;
		my $y_for_gene1 = $y + 20;

		#基因矩形
#		$svg->rect(x=>$x_for_gene1,y=>$y_for_gene1,width=>$len1,height=>15,"stroke-opacity"=>"1","stroke"=>"black",fill=>$color{$gene_pos->{$pos1}->[0]}||="red");
		if($gene_pos->{$pos1}->[1]){#负链
			$svg->rect(x=>$x_for_gene1+5,y=>$y_for_gene1,width=>$len1-5,height=>15,"stroke-opacity"=>"1","stroke"=>$color{$gene_pos->{$pos1}->[0]}||="red",fill=>$color{$gene_pos->{$pos1}->[0]}||="red");
			my $p0 = $x_for_gene1;#顶点x 
			my $p1 = $x_for_gene1+5;#底边x
			my $p2 = $y_for_gene1 + 7.5;#顶点y
			my $p3 = $y_for_gene1 + 15;#底y
			
			$svg->path("d"=>"M $p0 $p2 L $p1 $y_for_gene1 L $p1 $p3",style=>"stroke","stroke"=>$color{$gene_pos->{$pos1}->[0]}||="red","stroke-width"=>0.5,fill=>$color{$gene_pos->{$pos1}->[0]}||="red");  #箭头
		}else{#正链
			$svg->rect(x=>$x_for_gene1,y=>$y_for_gene1,width=>$len1-5,height=>15,"stroke-opacity"=>"1","stroke"=>$color{$gene_pos->{$pos1}->[0]}||="red",fill=>$color{$gene_pos->{$pos1}->[0]}||="red");
			my $p0 = $x_for_gene1+$len1;      
			my $p1 = $x_for_gene1+$len1-5;    
			my $p2 = $y_for_gene1 + 7.5;
			my $p3 = $y_for_gene1 + 15; 
			$svg->path("d"=>"M $p0 $p2 L $p1 $y_for_gene1 L $p1 $p3",style=>"stroke","stroke"=>$color{$gene_pos->{$pos1}->[0]}||="red","stroke-width"=>0.5,fill=>$color{$gene_pos->{$pos1}->[0]}||="red");  #箭头
		}

		#基因名字
		my $name_x1;
		my $name_y1;
		my $bp_y1;
		if($len1 < 45){
			$name_x1 = ($x_for_gene1 + $x_for_gene1 + $len1)/2;
#			$name_y1 = 	$gene_pos->{$pos1}->[1] ? $y + 20 + 15 + 15 : $y - 15 - 2;
			$name_y1 = 	$y + 20 + 15 + 15;
#			$bp_y1 = $gene_pos->{$pos1}->[1] ? $y + 20 + 15 - 2 : $y - 2;
			$bp_y1 = $y + 20 + 15 - 2;

			$svg->text(x=>$name_x1,y=>$name_y1,'font-size'=>13,'font-family'=>'times new roman','-cdata'=>"$gene_pos->{$pos1}->[0]","text-anchor"=>"end",fill=>"black",'font-style'=>'italic');
			$svg->text(x=>$x_for_gene1-2,y=>$bp_y1,'font-size'=>13,'font-family'=>'times new roman','-cdata'=>($gene_pos->{$pos1}->[3] - $gene_pos->{$pos1}->[2]+1)." bp","text-anchor"=>"end",fill=>"black",'font-style'=>'italic') if($is_print_gene_len);
		}else{
			$name_x1 = ($x_for_gene1 + $x_for_gene1 + $len1)/2;
#			$name_y1 = 	$gene_pos->{$pos1}->[1] ? $y + 20 + 15 - 2 : $y - 2;
			$name_y1 = 	 $y + 20 + 15 - 2;
#			$bp_y1 = $gene_pos->{$pos1}->[1] ? $y + 20 + 15 + 14 : $y - 15 - 1;
			$bp_y1 =  $y + 20 + 15 + 14;
			$svg->text(x=>$name_x1,y=>$name_y1,'font-size'=>13,'font-family'=>'times new roman','-cdata'=>"$gene_pos->{$pos1}->[0]","text-anchor"=>"middle",fill=>"white",'font-style'=>'italic');
			$svg->text(x=>$name_x1,y=>$bp_y1,'font-size'=>13,'font-family'=>'times new roman','-cdata'=>($gene_pos->{$pos1}->[3] - $gene_pos->{$pos1}->[2]+1)." bp","text-anchor"=>"middle",fill=>"black");
		}
		
		{	#gap距离 1
			use Math::Trig;
			my $x0 = $x - $gap1/2 - 25;	#圆心
#			my $y0 = $gene_pos->{$pos1}->[1] ? $y + 20 + 15 : $y - 15;
			my $y0 = $y + 20 + 15;
			my $r = $x - $x0 - $gap1/2;	#半径

			my $x1 = $x0 + $r * cos(0);			#外圆弧的起点和中点
			my $y1 = $y0 + $r * sin(0);

			if(1 or $gene_pos->{$pos1}->[1]){
				my $x2 = $x0 + $r * cos(pi*0.5);
				my $y2 = $y0 + $r * sin(pi*0.5);

				$svg->path("d"=>"M $x1  $y1 A $r $r 0 0 1 $x2  $y2",style=>"stroke", "stroke"=>"black","stroke-width"=>0.5, fill=>"none");	#圆弧

				my ($x21,$y21) = ($x2+3,$y2-3);	#箭头
				my ($x22,$y22) = ($x2+3,$y2+3);

				$svg->path("d"=>"M $x2  $y2 L $x21 $y21",style=>"stroke","stroke"=>"black","stroke-width"=>0.5,fill=>"none");
				$svg->path("d"=>"M $x2  $y2 L $x22 $y22",style=>"stroke","stroke"=>"black","stroke-width"=>0.5,fill=>"none");
				$svg->text(x=>$x0-2,y=>$y0+$r+5,'font-size'=>13,'font-family'=>'times new roman','-cdata'=>$junction_pos - $gene_pos->{$pos1}->[3]." bp","text-anchor"=>"end",fill=>"black"); #gap
			}else{
				my $x2 = $x0 + $r * cos(-pi*0.5);
				my $y2 = $y0 + $r * sin(-pi*0.5);
				$svg->path("d"=>"M $x1  $y1 A $r $r 0 0 0 $x2  $y2",style=>"stroke", "stroke"=>"black","stroke-width"=>0.5, fill=>"none");	#圆弧

				my ($x21,$y21) = ($x2+3,$y2-3);	#箭头
				my ($x22,$y22) = ($x2+3,$y2+3);

				$svg->path("d"=>"M $x2  $y2 L $x21 $y21",style=>"stroke","stroke"=>"black","stroke-width"=>0.5,fill=>"none");
				$svg->path("d"=>"M $x2  $y2 L $x22 $y22",style=>"stroke","stroke"=>"black","stroke-width"=>0.5,fill=>"none");
				$svg->text(x=>$x0-2,y=>$y0-$r+5,'font-size'=>13,'font-family'=>'times new roman','-cdata'=>$junction_pos - $gene_pos->{$pos1}->[3]." bp","text-anchor"=>"end",fill=>"black"); #gap
			}
		}


		#向右寻找

		my $pos2 = $junction_pos;
		#$pos2 = $pos2 - $end_ira if($junction_pos == $end_ira);  #直接在基因中添加了周期，此处无需判断
		#$junction_pos = 0 if($junction_pos == $end_ira);

		while(!defined $gene_pos->{$pos2} or !defined $gene_pos->{$pos2+1}){
			$pos2++;
		}
		print"\trigth: $gene_pos->{$pos2}->[0] $gene_pos->{$pos2}->[2] $gene_pos->{$pos2}->[3] $gene_pos->{$pos2}->[1] ".($gene_pos->{$pos2}->[3] - $gene_pos->{$pos2}->[2] +1)."\n";


#		$gene_pos->{$pos2}->[1] = $gene_pos->{$pos1}->[1]?0:1;
		my $len2 = get_gene_len(($gene_pos->{$pos2}->[3] - $gene_pos->{$pos2}->[2]));
		my $gap2 = get_gap_len($gene_pos->{$pos2}->[2] - $junction_pos);

		my $x_for_gene2 = $x + $gap2;

#		my $y_for_gene2 = $gene_pos->{$pos2}->[1] ? $y + 20 : $y - 15;
		my $y_for_gene2 =  $y - 15;
		
		#基因矩形
#		$svg->rect(x=>$x_for_gene2,y=>$y_for_gene2,width=>$len2,height=>15,"stroke-opacity"=>"1","stroke"=>"black",fill=>$color{$gene_pos->{$pos2}->[0]}||="red");
		if($gene_pos->{$pos2}->[1]){#负链
			$svg->rect(x=>$x_for_gene2+5,y=>$y_for_gene2,width=>$len2-5,height=>15,"stroke-opacity"=>"1","stroke"=>$color{$gene_pos->{$pos2}->[0]}||="red",fill=>$color{$gene_pos->{$pos2}->[0]}||="red");
			my $p0 = $x_for_gene2;#顶点x 
			my $p1 = $x_for_gene2+5;#底边x
			my $p2 = $y_for_gene2 + 7.5;#顶点y
			my $p3 = $y_for_gene2 + 15;#底y
			
			$svg->path("d"=>"M $p0 $p2 L $p1 $y_for_gene2 L $p1 $p3",style=>"stroke","stroke"=>$color{$gene_pos->{$pos2}->[0]}||="red","stroke-width"=>0.5,fill=>$color{$gene_pos->{$pos2}->[0]}||="red");  #箭头
		}else{#正链
			$svg->rect(x=>$x_for_gene2,y=>$y_for_gene2,width=>$len2-5,height=>15,"stroke-opacity"=>"1","stroke"=>$color{$gene_pos->{$pos2}->[0]}||="red",fill=>$color{$gene_pos->{$pos2}->[0]}||="red");
			my $p0 = $x_for_gene2+$len2;      
			my $p1 = $x_for_gene2+$len2-5;    
			my $p2 = $y_for_gene2 + 7.5;
			my $p3 = $y_for_gene2 + 15; 
			$svg->path("d"=>"M $p0 $p2 L $p1 $y_for_gene2 L $p1 $p3",style=>"stroke","stroke"=>$color{$gene_pos->{$pos2}->[0]}||="red","stroke-width"=>0.5,fill=>$color{$gene_pos->{$pos2}->[0]}||="red");  #箭头
		}
		#基因名字
		my $name_x2;
		my $name_y2;
		my $bp_y2;
		if($len2 < 45){
			$name_x2 = ($x_for_gene2 + $x_for_gene2 + $len2)/2;
#			$name_y2 = 	$gene_pos->{$pos2}->[1] ? $y + 20 + 15 + 15 : $y - 15 - 2;
			$name_y2 =  $y - 15 - 2;
#			$bp_y2 = $gene_pos->{$pos2}->[1] ? $y + 20 + 15 - 2 : $y - 2;
			$bp_y2 =  $y - 2;
			$svg->text(x=>$name_x2,y=>$name_y2,'font-size'=>13,'font-family'=>'times new roman','-cdata'=>"$gene_pos->{$pos2}->[0]","text-anchor"=>"start",fill=>"black",'font-style'=>'italic');
			$svg->text(x=>$x_for_gene2 + 2 + $len2,y=>$bp_y2,'font-size'=>13,'font-family'=>'times new roman','-cdata'=>($gene_pos->{$pos2}->[3] - $gene_pos->{$pos2}->[2]+1)." bp","text-anchor"=>"start",fill=>"black",'font-style'=>'italic') if($is_print_gene_len);
		}else{
			$name_x2 = ($x_for_gene2 + $x_for_gene2 + $len2)/2;
#			$name_y2 = 	$gene_pos->{$pos2}->[1] ? $y + 20 + 15 - 2 : $y - 2;
			$name_y2 =  $y - 2;
#			$bp_y2 = $gene_pos->{$pos2}->[1] ? $y + 20 + 15 + 14 : $y - 15 - 1;
			$bp_y2 = $y - 15 - 1;
			$svg->text(x=>$name_x2,y=>$name_y2,'font-size'=>13,'font-family'=>'times new roman','-cdata'=>"$gene_pos->{$pos2}->[0]","text-anchor"=>"middle",fill=>"white",'font-style'=>'italic');
			$svg->text(x=>$name_x2,y=>$bp_y2,'font-size'=>13,'font-family'=>'times new roman','-cdata'=>($gene_pos->{$pos2}->[3] - $gene_pos->{$pos2}->[2]+1)." bp","text-anchor"=>"middle",fill=>"black");
		}

		{	#gap距离 2
			use Math::Trig;
			my $x0 = $x + $gap2/2 + 25;	#圆心
#			my $y0 = $gene_pos->{$pos2}->[1] ? $y + 20 + 15 : $y - 15;
			my $y0 =  $y - 15;
			my $r = $x0 - $x - $gap2/2;	#半径

			my $x1 = $x0 + $r * cos(pi);			#外圆弧的起点和中点
			my $y1 = $y0 + $r * sin(pi);

			if(0 and $gene_pos->{$pos2}->[1]){
				my $x2 = $x0 + $r * cos(pi*0.5);
				my $y2 = $y0 + $r * sin(pi*0.5);

				$svg->path("d"=>"M $x1  $y1 A $r $r 0 0 0 $x2  $y2",style=>"stroke", "stroke"=>"black","stroke-width"=>0.5, fill=>"none");	#圆弧

				my ($x21,$y21) = ($x2-3,$y2-3);	#箭头
				my ($x22,$y22) = ($x2-3,$y2+3);

				$svg->path("d"=>"M $x2  $y2 L $x21 $y21",style=>"stroke","stroke"=>"black","stroke-width"=>0.5,fill=>"none");
				$svg->path("d"=>"M $x2  $y2 L $x22 $y22",style=>"stroke","stroke"=>"black","stroke-width"=>0.5,fill=>"none");

				$svg->text(x=>$x0+2,y=>$y0+$r+5,'font-size'=>13,'font-family'=>'times new roman','-cdata'=>($gene_pos->{$pos2}->[2] - $junction_pos - 1)." bp","text-anchor"=>"start",fill=>"black"); #gap
			}else{
				my $x2 = $x0 + $r * cos(-pi*0.5);
				my $y2 = $y0 + $r * sin(-pi*0.5);
				$svg->path("d"=>"M $x1  $y1 A $r $r 0 0 1 $x2  $y2",style=>"stroke", "stroke"=>"black","stroke-width"=>0.5, fill=>"none");	#圆弧

				my ($x21,$y21) = ($x2-3,$y2-3);	#箭头
				my ($x22,$y22) = ($x2-3,$y2+3);

				$svg->path("d"=>"M $x2  $y2 L $x21 $y21",style=>"stroke","stroke"=>"black","stroke-width"=>0.5,fill=>"none");
				$svg->path("d"=>"M $x2  $y2 L $x22 $y22",style=>"stroke","stroke"=>"black","stroke-width"=>0.5,fill=>"none");
				$svg->text(x=>$x0+2,y=>$y0-$r+5,'font-size'=>13,'font-family'=>'times new roman','-cdata'=>($gene_pos->{$pos2}->[2] - $junction_pos - 1)." bp","text-anchor"=>"start",fill=>"black"); #gap
			}
		}
	}

}

sub get_gene_len{
	my $len = shift;
	if($len<10){
		5
	}elsif($len<20){
		7
	}elsif($len<30){
		9	
	}elsif($len<40){
		11	
	}elsif($len<50){
		20	
	}elsif($len<100){
		45
	}elsif($len<200){
		50
	}elsif($len<300){
		60
	}elsif($len<400){
		61
	}elsif($len<500){
		62
	}elsif($len<600){
		64
	}elsif($len<700){
		66
	}elsif($len<800){
		68
	}elsif($len<900){
		70
	}elsif($len<1000){
		72
	}elsif($len<2000){
		76
	}elsif($len<3000){
		80
	}elsif($len<4000){
		85
	}elsif($len<5000){
		90
	}else{
		100
	}
}


sub get_gap_len{
	my $len = shift;

	if($len == 0 ){
		0
	}elsif($len<5){
		2
	}elsif($len<10){
		5
	}elsif($len<30){
		9	
	}elsif($len<40){
		11	
	}elsif($len<50){
		13	
	}elsif($len<100){
		15
	}elsif($len<200){
		17
	}elsif($len<300){
		19
	}elsif($len<400){
		21
	}elsif($len<500){
		21
	}elsif($len<600){
		22
	}elsif($len<700){
		23
	}elsif($len<800){
		24
	}elsif($len<900){
		25
	}elsif($len<1000){
		25
	}elsif($len<2000){
		25
	}elsif($len<3000){
		27
	}elsif($len<4000){
		29
	}elsif($len<5000){
		29
	}else{
		30
	}
}

sub get_gene_pos{
	my $infile = shift;
	my $T = shift;
	die"no ira pos" unless($T);
	my @all = get_ref_gbk_gene_info($infile);

	my %for_gene;

	for my $gene_info(@all){
		#print "$gene_info\n";
		my $gene_name = (split":",$gene_info)[0];
		my $pos_info = (split":",$gene_info)[1];
		my @pos = $pos_info =~ /(\d+)/g;
		my ($gene_pos_s,$gene_pos_e) = (sort{$a <=> $b} @pos )[0,-1];
		my $gene_tpye = $pos_info =~ /complement/ ?1:0;
		if($gene_pos_e - $gene_pos_s <= 10_000 and $gene_name !~ /orf/i){
			#print"$gene_name:$gene_pos_s,$gene_pos_e,$gene_tpye\n";
			for($gene_pos_s..$gene_pos_e){
				$for_gene{$_} = [$gene_name,$gene_tpye,$gene_pos_s,$gene_pos_e];
				$for_gene{$_+$T} = [$gene_name,$gene_tpye,$gene_pos_s+$T,$gene_pos_e+$T];
			}
		}elsif($gene_name !~ /rps12/i and $gene_name !~ /orf/i and $pos[2] == 1){
			for($pos[0]..$pos[3]+$T){
				$for_gene{$_} = [$gene_name,$gene_tpye,$pos[0],$pos[3]+$T];
			}
		}
	}
	return %for_gene;
}

sub get_ref_gbk_gene_info{
	my $file = shift;
	open IN ,$file or die"$!";
	$/ = "--------adf32adsd-----\n\n";
	my $info = <IN>;
	$info =~ s/\n {20}//gm;
	$info =~ s/ +/ /g;
	$info =~ s/, /,/g;
	$info =~ s/\n\s+/\n/g;
	my @infos = split/\n/,$info;
	my @sample_gene;
	for(@infos){
	   if(/^CDS|^tRNA|^rRNA/){
			
			my $pos_info = (split/ /,$_)[1];		
			my ($gene_name) = /\/gene=(\S+)/;

	        if(!$gene_name){
				if(/^rRNA/ and /product="23/){
					$gene_name = "rrn23";
				}elsif(/^rRNA/ and /product="16/){
					$gene_name = "rrn16";
				}elsif(/^rRNA/ and /product="4\.5/){
					$gene_name = "rrn4.5";
				}elsif(/^rRNA/ and /product="5/){
					$gene_name = "rrn5";
				}else{
					 warn "$_ \n is not exists gene_name\n";
					 next;
				}
	        }

			if($gene_name =~ /trn/i){
				$gene_name =~ s/T/U/g;
				$gene_name =~ s/trnU/trnT/g;
				$gene_name =~ s/[-_]...//g;
			}
			$gene_name =~ s/"//g;
			$gene_name =~ s/rrn(.*)S/rrn$1/;

	        push @sample_gene,$gene_name.":".$pos_info;
	   }elsif($pseudo and /gene/ and /pseudo/){
			my $pos_info = (split/ /,$_)[1];
			my ($gene_name) = /\/gene=(\S+) /;

			if($gene_name =~ /trn/i){
				$gene_name =~ s/T/U/g;
				$gene_name =~ s/trnU/trnT/g;
			}
			$gene_name =~ s/"//g;
			$gene_name =~ s/rrn(.*)S/rrn$1/;
			push @sample_gene,"%".$gene_name.":".$pos_info;
	   }
	}
	$/ = "\n";
	close IN;
	return @sample_gene;
}



sub get_seq{
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

sub get_sample_name{
	my $file = shift;

	open IN,$file or die"$!";
	$/="ORIGIN";
	my $info = <IN>;
	my ($sample_name) = $info =~ /ORGANISM\s+(.*?)\n/;
	close IN;
	$/="\n";
	return $sample_name;
}

sub find_region{
	my $len_of_seed = 1100;	
	my $seq = shift;
	my $re_seq = reverse $seq;
	$re_seq =~ tr/ATGC/TACG/;
	my $seq_long = length($seq);
	die "seq is too short!" if length($seq) < 100;
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
			if($len_of_seed <200){
				warn "no repeat larger than 200bp\n";
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
			return ($start,$end,$start2,$end2);
		}else{
			die"your gbk file is not format!\n";
		}
	}
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
	-i	<infile>	input cfg file  "file  ir" eg.   sample.gbk   1000-2000,4000-5000
	-l	left distance default 120
	-pseudo	draw pseudo gene
	-o	<output>	out.name(.svg)  default irscope.svg
	-h				Help
USAGE
	print $usage;
	exit;
}
