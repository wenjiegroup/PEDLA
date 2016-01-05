#!/usr/bin/env perl
use strict;
use warnings;
use 5.010;

###<discription>this program calculates the RPKM-like value of all bins

###input of the program is in the directory "Tag_count_of_all_bins_of_all_chrs/phastCons46way-step200win200"
###output of the program is in the corresponding sub-directories of the directory "RPKM/phastCons46way-step200win200"


my $minimal_step=200; ## resolution, which means the genome is divided into 200 bp intervals.
my $win=$minimal_step;   #size of window


my $dir_tag_count="Tag_count_of_all_bins_of_all_chrs/phastCons46way-step${minimal_step}win$win";   
opendir(DIR_TAG,"$dir_tag_count") or die("open $dir_tag_count error!\n$!");
my @marker_tag=sort grep { -d "$dir_tag_count/$_" && !/^\./ } readdir(DIR_TAG);

my $dir_out="RPKM";
mkdir $dir_out unless -e $dir_out;
$dir_out="$dir_out/phastCons46way-step${minimal_step}win$win";
mkdir $dir_out unless -e $dir_out;


foreach my $marker(@marker_tag)
{
	
	# next unless $marker=~/V_AFP1_Q6/;
	
	my $dir_out="$dir_out/$marker";
	mkdir $dir_out unless -e $dir_out;
	
	say "calculating RPKM of $marker...";
	
	my $N=0;
	foreach my $ii(1..23)
	{
		my $chr="chr$ii";
		$chr="chrX" if $ii==23;
		open(TOTAL,"<","$dir_tag_count/$marker/$marker.$chr.Total") or die("$!");
		$_=<TOTAL>;
		s/[\r\n]+$//;
		$N+=$_;
	}
	say "N=$N";
	
	
	
	foreach my $ii(1..23)
	{
		my $chr="chr$ii";
		$chr="chrX" if $ii==23;
		
		# next if -s "$dir_out/".'all'."-".$marker.".$chr.rpkm";
		say "\t$chr";
		open(OUT,">","$dir_out/".'all'."-".$marker.".$chr.rpkm");
		open(TAG,"<","$dir_tag_count/$marker/$marker.$chr.TagCount") or die("$!");
		<TAG>;
		while(<TAG>)
		{
			s/[\r\n]+$//;
			my $rpkm.=$_*(10**9)/$N/$win;
			say OUT $rpkm;
		}
		close(TAG);
	}

}




