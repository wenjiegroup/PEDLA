#!/usr/bin/env perl
use strict;
use warnings;
use 5.010;

###<discription>this program counts the number of tags (adds up the total signal intensity) of all the bins of all chromosomes, so we just need to calculate once and can use it for later calculation of RPKM-like value.

###input of the program is in the directory "Raw_Data/phastCons46way"
###output of the program is in the corresponding sub-directories of the directory "Tag_count_of_all_bins_of_all_chrs/phastCons46way-step200win200"


my $minimal_step=200; ## resolution, which means the genome is divided into 200 bp intervals.
my $win=$minimal_step;   #size of window

my $dir_tag="Raw_Data/phastCons46way";   

my $dir_out="Tag_count_of_all_bins_of_all_chrs";
mkdir $dir_out unless -e $dir_out;
$dir_out="$dir_out/phastCons46way-step${minimal_step}win$win";
mkdir $dir_out unless -e $dir_out;


my $marker='phastCons46way';
$dir_out="$dir_out/$marker";
mkdir $dir_out unless -e $dir_out;

say "calculating tag count of $marker...";

foreach my $ii(1..23)
{
	my $chr="chr$ii";
	$chr="chrX" if $ii==23;
	
	# next unless $chr eq "chr1";
	# next if -s "$dir_out/$marker.$chr.TagCount" && "$dir_out/$marker.$chr.Total";
	
	say "\t$chr";
	
	open(OUT,">","$dir_out/$marker.$chr.TagCount");
	open(TOTAL,">","$dir_out/$marker.$chr.Total");
	
	say "\tinitialize $chr...";
	open(BIN,"<","./human-hg19-size.bed") or die($!);
	my @tag_of_bin;
	while(<BIN>)
	{
		my @temp=split;
		next unless $temp[0] eq $chr;
		for (my $i=0;$i<$temp[-1]/$minimal_step;$i++)
		{
			$tag_of_bin[$i]=0;
		}
	}
	my $file_marker="$chr.phastCons46way.wigFix";
	open(TAG,"<","$dir_tag/$file_marker") or die("$dir_tag/$file_marker do not exsit!! download and prepare the data first.$!");
	my $N=0;
	
	say "\tcounting tag...";
	my ($chrom,$start,$step);
	while(<TAG>)
	{
		my @temp=split;
		if(@temp>1)
		{
			die("$_$!") unless /^fixedStep/ && @temp==4;
			$chrom=(split(/=/,$temp[1]))[1];
			$start=(split(/=/,$temp[2]))[1];
			$step=(split(/=/,$temp[3]))[1];
			die("$chr ne $chrom$!") unless $chrom eq $chr;
		}
		elsif(@temp==1)
		{
			
			$N+=$temp[-1];
			
			my ($b1,$b2)=($start-$win/2,$start+$win/2);	#middle of bin 
			($b1,$b2)=($b1-$minimal_step/2,$b2-$minimal_step/2);	#start of bin
			if($b1%$minimal_step==0)
			{
				$b1=$b1/$minimal_step;
				$b1=0 if $b1<0;
			}
			else
			{
				if($b1>=0){
					$b1=(int($b1/$minimal_step)+1);
				}
				else{
					$b1=0;
				}
			}
			if($b2%$minimal_step==0)
			{
				$b2=$b2/$minimal_step;
			}
			else
			{
				$b2=(int($b2/$minimal_step)+0);
			}
			for(my $i=$b1;$i<=$b2;$i++){
				if( defined($tag_of_bin[$i]) )
				{
					$tag_of_bin[$i]+=$temp[-1];
				}
			}
			$start+=$step;
		}
		else
		{
			die(scalar(@temp));
		}
	}
	say TOTAL $N;
	
	say "\twriting...";
	say OUT "#$chr\t0\t$minimal_step";
	say OUT foreach @tag_of_bin;
}



