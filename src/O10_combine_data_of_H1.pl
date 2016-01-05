#!/usr/bin/env perl
use strict;
use warnings;
use 5.010;

###<discription>combine 9 categories of heterogeneous data to form matrix as we described in our paper. the total dimensionality of all the heterogeneous data is 1114.

###input of the program is in sub-directories of the directory "RPKM/Histone_Modification_average-step200win200","RPKM/phastCons46way-step200win200" and so on.
###output of the program is in the corresponding sub-directories of the directory "RPKM/Data_matrix_combined-step200win200.H1"


my $minimal_step=200; ## resolution, which means the genome is divided into 200 bp intervals.
my $win=$minimal_step;   #size of window

my $dir_tag_count="RPKM/Histone_Modification_average-step${minimal_step}win$win";   
opendir(DIR,$dir_tag_count) or die("open $dir_tag_count error!\n$!");
my %cell_with_marker;
$cell_with_marker{$_}=0 foreach sort grep { -d "$dir_tag_count/$_" && !/^\./} readdir(DIR);

my $dir_out="RPKM";
mkdir $dir_out unless -e $dir_out;
$dir_out="$dir_out/Data_matrix_combined-step${minimal_step}win$win.H1";
mkdir $dir_out unless -e $dir_out;



foreach my $cell(sort keys %cell_with_marker)
{
	next unless $cell=~/^H1_cell_line/i;
	say "$cell";
	
	my $dir_out="$dir_out/$cell";
	mkdir $dir_out unless -e $dir_out;
	
	opendir(DIR_TAG,"$dir_tag_count/$cell") or die("open $dir_tag_count/$cell error!\n$!");
	my @marker_tag=sort by_my_way grep { -d "$dir_tag_count/$cell/$_" && !/^\./ } readdir(DIR_TAG);	#histone/DHS/mRNA/TF
	
	opendir(DIR_TAG,"$dir_tag_count/../CpGIsland-step${minimal_step}win$win") or die("open error!\n$!");
	my @cpg_tag=sort  grep { -d "$dir_tag_count/../CpGIsland-step${minimal_step}win$win/$_" && !/^\./ } readdir(DIR_TAG);	#CpGIsland
	
	opendir(DIR_TAG,"$dir_tag_count/../Methyl_RRBS_average-step${minimal_step}win$win/$cell") or die("open error!\n$!");
	my @rrbs_tag=sort  grep { -d "$dir_tag_count/../Methyl_RRBS_average-step${minimal_step}win$win/$cell/$_" && !/^\./ } readdir(DIR_TAG);	#Methyl_RRBS_average
	
	opendir(DIR_TAG,"$dir_tag_count/../phastCons46way-step${minimal_step}win$win") or die("open error!\n$!");
	my @phast_tag=sort  grep { -d "$dir_tag_count/../phastCons46way-step${minimal_step}win$win/$_" && !/^\./ } readdir(DIR_TAG);	#phastCons46way
	
	opendir(DIR_TAG,"$dir_tag_count/../motif-step${minimal_step}win$win") or die("open error!\n$!");
	my @motif_tag=sort  grep { -d "$dir_tag_count/../motif-step${minimal_step}win$win/$_" && !/^\./ } readdir(DIR_TAG);	#Motif
	
	opendir(DIR_TAG,"$dir_tag_count/../Sequence_feature-step${minimal_step}win$win") or die("open error!\n$!");
	my @seq_tag=sort { length($a) <=> length($b) or $a cmp $b } grep { -d "$dir_tag_count/../Sequence_feature-step${minimal_step}win$win/$_" && !/^\./ } readdir(DIR_TAG);	#Sequence_feature
	
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
		
		foreach my $marker(@marker_tag)#histone/DHS/mRNA/TF
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
		foreach my $marker(@cpg_tag)#CpGIsland
		{
			print "  $marker";
			
			$header.=$marker."\t";
			my @temp;
			open(IN,"<","$dir_tag_count/../CpGIsland-step${minimal_step}win$win/$marker/all-$marker.$chr.rpkm") or die($!);
			
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
		foreach my $marker(@rrbs_tag)#Methyl_RRBS_average
		{
			print "  $marker";
			
			$header.=(split(/-/,$marker))[-1]."\t";
			my @temp;
			open(IN,"<","$dir_tag_count/../Methyl_RRBS_average-step${minimal_step}win$win/$cell/$marker/$marker.$chr.rpkm") or die($!);
			
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
		foreach my $marker(@motif_tag)#motif
		{
			print "  $marker";
			
			$header.=$marker."\t";
			my @temp;
			open(IN,"<","$dir_tag_count/../motif-step${minimal_step}win$win/$marker/all-$marker.$chr.rpkm") or die($!);
			
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
		foreach my $marker(@seq_tag)#Sequence_feature
		{
			print "  $marker";
			
			$header.=$marker."\t";
			my @temp;
			open(IN,"<","$dir_tag_count/../Sequence_feature-step${minimal_step}win$win/$marker/all-$marker.$chr.rpkm") or die($!);
			
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