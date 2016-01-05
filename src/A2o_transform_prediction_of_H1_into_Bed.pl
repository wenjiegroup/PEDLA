#!/usr/bin/env perl
use strict;
use warnings;
use 5.010;

###<discription>this program gets the segmentations of the whole genome using PEDLA and tansforms them to the 'bed' format. The fourth column of '.bed' file indicates the state with 'E1' indicating enhancer and 'E2' indicating non-enhancer. 

my $num_state=2;
my $len=200;


my $dir_in="Result_of_prediction.H1";
my $dir_out=$dir_in;
mkdir $dir_out unless -e $dir_out;

opendir(DIR,$dir_in) or die($!);
my @dirs=sort grep { -d "$dir_in/$_" && !/^\./ } readdir(DIR);

foreach my $dir(@dirs)
{
	# next unless $dir=~/Bj_Rep/;
	
	say $dir;
	
	my $dir_out="$dir_out/$dir";
	mkdir $dir_out unless -e $dir_out;
	open(OUT,">","$dir_out/${dir}_${num_state}_segments.bed");
	foreach my $i(1..23)
	{
		my $chr="chr$i";
		$chr="chrX" if $i==23;
		open(IN1,"<","$dir_in/$dir/PEDLA_prediction_of_${dir}.$chr.txt") or die($!);
		
		# next unless $chr eq "chr1";
		
		my $pre_state='';
		my $count=0;
		my ($start,$end)=(0,0);
		while( defined(my $line1=<IN1> ) )   #gets the refined segmentations
		{
			my @temp1=split(/\s+/,$line1);
			die("a row of new state contains more than 1 state!\n") unless @temp1==1;
			if($temp1[0] ne $pre_state)
			{
				say OUT "$chr\t$start\t$end\tE$pre_state" if $end;
				$start=$end;
			}
			$end=($count+1)*$len;
			
			
			$count++;
			$pre_state=$temp1[0];
		}
		say OUT "$chr\t$start\t$end\tE$pre_state" if $end;
	}
	
}
