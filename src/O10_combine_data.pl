#!/usr/bin/env perl
use strict;
use warnings;
use 5.010;

###<discription>combine histone modification/phastCons46way data to form matrix as we described in our paper. histone modification includes H3K4me1|H3K4me3|H3K9me3|H3K27me3|H3K36me3|H3K27ac|input, so the total dimensionality of histone modification/phastCons46way is 8.

###input of the program is in sub-directories of the directory "RPKM/Histone_Modification_average-step200win200","RPKM/phastCons46way-step200win200"
###output of the program is in the corresponding sub-directories of the directory "RPKM/Data_matrix_combined-step200win200"


my $minimal_step=200; ## resolution, which means the genome is divided into 200 bp intervals.
my $win=$minimal_step;   #size of window

my $dir_tag_count="RPKM/Histone_Modification_average-step${minimal_step}win$win";   
opendir(DIR,$dir_tag_count) or die("open $dir_tag_count error!\n$!");
my %cell_with_marker;
$cell_with_marker{$_}=0 foreach sort grep { -d "$dir_tag_count/$_" && !/^\./} readdir(DIR);

my $dir_out="RPKM";
mkdir $dir_out unless -e $dir_out;
$dir_out="$dir_out/Data_matrix_combined-step${minimal_step}win$win";
mkdir $dir_out unless -e $dir_out;



foreach my $cell(sort keys %cell_with_marker)
{
	# next unless $cell=~/^Imr90/i;
	say "$cell";
	
	my $dir_out="$dir_out/$cell";
	mkdir $dir_out unless -e $dir_out;
	
	opendir(DIR_TAG,"$dir_tag_count/$cell") or die("open $dir_tag_count/$cell error!\n$!");
	my @marker_tag=sort by_my_way grep { -d "$dir_tag_count/$cell/$_" && !/^\./ && /-(H3K4me1|H3K4me3|H3K9me3|H3K27me3|H3K36me3|H3K27ac|input)/ } readdir(DIR_TAG);	#histone modification
	
	opendir(DIR_TAG,"$dir_tag_count/../phastCons46way-step${minimal_step}win$win") or die("open error!\n$!");
	my @phast_tag=sort  grep { -d "$dir_tag_count/../phastCons46way-step${minimal_step}win$win/$_" && !/^\./ } readdir(DIR_TAG);	#phastCons46way
	
	
	foreach my $ii(1..23)
	{
		my $chr="chr$ii";
		$chr="chrX" if $ii==23;
		
		# next if -e "$dir_out/$cell.$chr.rpkm";
		say "\t$chr";
		open(OUT,">","$dir_out/$cell.$chr.rpkm");
		
		my $header="$cell\t$chr\n";
		my @data;
		
		print "\t";
		
		foreach my $marker(@marker_tag)#histone modification
		{
			print "  $marker";
			
			$header.=(split(/-/,$marker))[-1]."\t";
			my @temp;
			open(IN,"<","$dir_tag_count/$cell/$marker/$marker.$chr.rpkm") or die($!);
			
			while(<IN>)
			{
				s/[\r\n]+$//;
				push @temp,sprintf("%.4g",$_);
			}
			close(IN);
			if(@data)
			{
				die("number of $marker is inconsistent$!") unless @data==@temp;
				foreach my $ind(0..$#data)
				{
					$data[$ind].="\t".$temp[$ind];
				}
			}
			else
			{
				@data=@temp;
			}
			
		}
		
		foreach my $marker(@phast_tag)#phastCons46way
		{
			print "  $marker";
			
			$header.=$marker."\t";
			my @temp;
			open(IN,"<","$dir_tag_count/../phastCons46way-step${minimal_step}win$win/$marker/all-$marker.$chr.rpkm") or die($!);
			
			while(<IN>)
			{
				s/[\r\n]+$//;
				push @temp,sprintf("%.4g",$_);
			}
			close(IN);
			if(@data)
			{
				die("number of $marker is inconsistent$!") unless @data==@temp;
				foreach my $ind(0..$#data)
				{
					$data[$ind].="\t".$temp[$ind];
				}
			}
			else
			{
				@data=@temp;
			}
			
		}
		
		
		print "\n";
		
		say "\twriting...";
		$header=~s/\t$//;
		say OUT $header;
		say OUT foreach @data;
	}
}




sub by_my_way
{
	my ($score_a,$score_b);
	
	if($a=~/-input$/)
	{
		$score_a=0;
	}
	elsif($a=~/-DNase_hypersensitivity$/)
	{
		$score_a=2;
	}
	elsif($a=~/-H\d\w+$/)
	{
		$score_a=1;
	}
	elsif($a=~/-mRNA_Seq$/)
	{
		$score_a=3;
	}
	elsif($a=~/-\w+$/)
	{
		$score_a=4;
	}
	else
	{
		die("$a!$!");
	}
	
	if($b=~/-input$/)
	{
		$score_b=0;
	}
	elsif($b=~/-DNase_hypersensitivity$/)
	{
		$score_b=2;
	}
	elsif($b=~/-H\d\w+$/)
	{
		$score_b=1;
	}
	elsif($b=~/-mRNA_Seq$/)
	{
		$score_b=3;
	}
	elsif($b=~/-\w+$/)
	{
		$score_b=4;
	}
	else
	{
		die("$b!$!");
	}
	
	$score_a <=> $score_b or $a cmp $b;
}