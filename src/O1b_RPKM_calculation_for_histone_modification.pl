#!/usr/bin/env perl
use strict;
use warnings;
use 5.010;

###<discription>this program calculates the RPKM of all bins based on counts

###input of the program is in sub-directories of the directory "Tag_count_of_all_bins_of_all_chrs/Histone_Modification-step200win200"
###output of the program is in the corresponding sub-directories of the directory "RPKM/Histone_Modification-step200win200"


my $minimal_step=200; ## resolution, which means the genome is divided into 200 bp intervals.
my $win=$minimal_step;   #size of window


my $dir_tag_count="Tag_count_of_all_bins_of_all_chrs/Histone_Modification-step${minimal_step}win$win";   
opendir(DIR,$dir_tag_count) or die("open $dir_tag_count error!\n$!");
my %cell_with_marker;
$cell_with_marker{$_}=0 foreach sort grep { -d "$dir_tag_count/$_" && !/^\./} readdir(DIR);

my $dir_out="RPKM";
mkdir $dir_out unless -e $dir_out;
$dir_out="$dir_out/Histone_Modification-step${minimal_step}win$win";
mkdir $dir_out unless -e $dir_out;


##read histone marker
foreach my $cell(sort keys %cell_with_marker)
{
	# next unless $cell=~/^Imr90/i;
	say "$cell";
	
	my $dir_out="$dir_out/$cell";
	mkdir $dir_out unless -e $dir_out;
	
	opendir(DIR_TAG,"$dir_tag_count/$cell") or die("open $dir_tag_count/$cell error!\n$!");
	my @marker_tag=sort grep { -d "$dir_tag_count/$cell/$_" && !/^\./ } readdir(DIR_TAG);
	foreach my $marker(@marker_tag)
	{
		# next if $marker=~/Methyl/;
		# next unless $marker=~/test/;
		
		
		say "\tcalculating RPKM of $marker...";
		
		my $dir_out="$dir_out/$marker";
		mkdir $dir_out unless -e $dir_out;
		
		my $N=0;
		foreach my $ii(1..23)
		{
			my $chr="chr$ii";
			$chr="chrX" if $ii==23;
			open(TOTAL,"<","$dir_tag_count/$cell/$marker/$marker.$chr.Total") or die("$!");
			$_=<TOTAL>;
			s/[\r\n]+$//;
			die("reads less than expected $chr$!") if $_<=500;
			$N+=$_;
		}
		say "\tN=$N";
		
		
		
		foreach my $ii(1..23)
		{
			my $chr="chr$ii";
			$chr="chrX" if $ii==23;
			
			# next if -s "$dir_out/".$marker.".$chr.rpkm";
			say "\t\t$chr";
			open(OUT,">","$dir_out/".$marker.".$chr.rpkm");
			open(TAG,"<","$dir_tag_count/$cell/$marker/$marker.$chr.TagCount") or die("$!");
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
}




