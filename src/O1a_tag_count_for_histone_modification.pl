#!/usr/bin/env perl
use strict;
use warnings;
use 5.010;

###<discription>this program counts the number of tags of all the bins of all chromosomes, so we just need to count once and can use it for later calculation of RPKM. 

###input of the program is in sub-directories of the directory "Raw_Data/Histone_Modification"
###output of the program is in the corresponding sub-directories of the directory "Tag_count_of_all_bins_of_all_chrs/Histone_Modification-step200win200"


my $minimal_step=200; ## resolution, which means the genome is divided into 200 bp intervals.
my $win=$minimal_step;   #size of window

my $dir_tag="Raw_Data/Histone_Modification";
opendir(DIR,$dir_tag) or die("open $dir_tag error!\n$!");
my %cell_with_marker;
$cell_with_marker{$_}=0 foreach sort grep { -d "$dir_tag/$_" && !/^\./} readdir(DIR);

my @cells=sort keys %cell_with_marker;

my $dir_out="Tag_count_of_all_bins_of_all_chrs";
mkdir $dir_out unless -e $dir_out;
$dir_out="$dir_out/Histone_Modification-step${minimal_step}win$win";
mkdir $dir_out unless -e $dir_out;

foreach my $cell(@cells)
{
	# next unless $cell=~/H1_cell_line/;
	
	my $dir_out="$dir_out/$cell";
	mkdir $dir_out unless -e $dir_out;


	opendir(DIR_TAG,"$dir_tag/$cell") or die("open $dir_tag/$cell error!\n$!");
	my @marker_tag=sort grep { !/^\./ } readdir(DIR_TAG);
	foreach my $file_marker(@marker_tag)
	{
		# next unless $file_marker=~/test/;
		
		my $marker=$cell.'-'.(split(/\./,$file_marker))[1].'-'.(split(/\./,$file_marker))[2];
		my $dir_out="$dir_out/$marker";
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
			open(TAG,"<","$dir_tag/$cell/$file_marker") or die("$!");
			my $N=0;
			
			say "\tcounting tag...";
			while(<TAG>)
			{
				my @temp=split;
				next unless $temp[0] eq $chr;
				
				$N++;
				
				my ($b1,$b2)=($temp[1]-$win/2,$temp[1]+$win/2);	#middle of bin 
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
						$tag_of_bin[$i]++;
					}
				}
			}
			say TOTAL $N;
			
			say "\twriting...";
			say OUT "#$chr\t0\t$minimal_step";
			say OUT foreach @tag_of_bin;
		}
		

	}
}