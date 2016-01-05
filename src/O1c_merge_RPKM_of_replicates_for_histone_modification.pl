#!/usr/bin/env perl
use strict;
use warnings;
use 5.010;

###<discription>merge different replicates by averaging them.

###input of the program is in sub-directories of the directory "RPKM/Histone_Modification-step200win200"
###output of the program is in the corresponding sub-directories of the directory "RPKM/Histone_Modification_average-step200win200"


my $minimal_step=200; ## resolution, which means the genome is divided into 200 bp intervals.
my $win=$minimal_step;   #size of window


my $dir_tag_count="RPKM/Histone_Modification-step${minimal_step}win$win";   
opendir(DIR,$dir_tag_count) or die("open $dir_tag_count error!\n$!");
my %cell_with_marker;
$cell_with_marker{$_}=0 foreach sort grep { -d "$dir_tag_count/$_" && !/^\./} readdir(DIR);

my $dir_out="RPKM";
mkdir $dir_out unless -e $dir_out;
$dir_out="$dir_out/Histone_Modification_average-step${minimal_step}win$win";
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
	my %rep_of_mark;
	$rep_of_mark{ (split(/-Rep/,$_))[0] }++ foreach @marker_tag;
	
	foreach my $marker(sort keys %rep_of_mark)
	{
		# next if $marker=~/Methyl/;
		# next unless $marker=~/test/;
		
		my $n=$rep_of_mark{$marker};
		
		say "\tcalculating average RPKM of $marker of $n replicates...";
		
		my $dir_out="$dir_out/$marker";
		mkdir $dir_out unless -e $dir_out;
		
		
		foreach my $ii(1..23)
		{
			my $chr="chr$ii";
			$chr="chrX" if $ii==23;
			
			# next if -s "$dir_out/".$marker.".$chr.rpkm";
			say "\t\t$chr";
			open(OUT,">","$dir_out/".$marker.".$chr.rpkm");
			
			my @data;
			foreach my $j(1..$n)
			{
				open(TAG,"<","$dir_tag_count/$cell/$marker-Rep$j/$marker-Rep$j.$chr.rpkm") or die("$dir_tag_count/$cell/$marker-Rep$j/$marker-Rep$j.$chr.rpkm$!");
				my @temp;
				while(<TAG>)
				{
					s/[\r\n]+$//;
					push @temp,$_;
				}
				close(TAG);
				if(@data)
				{
					die("number of rep$j is inconsistent$!") unless @data==@temp;
					foreach my $ind(0..$#data)
					{
						$data[$ind]+=$temp[$ind];
					}
				}
				else
				{
					@data=@temp;
				}
			}
			
			say OUT ($_/$n) foreach @data;
			close(OUT);
		}

	}
}
