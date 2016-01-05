#!/usr/bin/env perl
use strict;
use warnings;
use 5.010;

###<discription>this program extracts unlabeled features of the whole genome for PEDLA to predict. 

###input features of the program are in sub-directories of the directory "RPKM/Data_matrix_combined-step200win200"
###output of the program is in the corresponding sub-directories of the directory "Data_for_training_and_prediction/unlabelled_data-step200win200". In the output directory there will be about 23 files, where each corresponds to a chromosome. The 23 files prefixed by "x_" are the features.


my $filter_chr=undef;
if(@ARGV)
{
	$filter_chr=$ARGV[0];
}


my $num_state=2;	
my $num_negative_set=9;	# ratio 

my $flag='training_set';  #  training_set

my $bin_size=200;	## resolution, which means the genome is divided into 200 bp intervals.
my $win=$bin_size;


my $dir_tag_count="RPKM/Data_matrix_combined-step${bin_size}win$win";   
opendir(DIR,$dir_tag_count) or die("open $dir_tag_count error!\n$!");
my %cell_with_marker;
$cell_with_marker{$_}=0 foreach sort grep { -d "$dir_tag_count/$_" && !/^\./} readdir(DIR);


foreach my $cell(sort keys %cell_with_marker)
{
	# next unless $cell=~/H1_cell_line/;
	
	say "$cell";
	

	my $dir_in="$dir_tag_count/$cell";
	my $dir_out="Data_for_training_and_prediction";
	mkdir $dir_out unless -e $dir_out;
	$dir_out.="/unlabelled_data-step${bin_size}win$win";
	mkdir $dir_out unless -e $dir_out;
	$dir_out.="/${cell}";
	mkdir $dir_out unless -e $dir_out;


	foreach my $i(1..23)
	{
		# next unless $i==1 || $i==3 || $i==4 || $i==10 || $i==14 || $i==17;  
		
		my $chr="chr$i";
		$chr="chrX" if $i==23;
		
		if(defined($filter_chr))
		{
			next unless $chr eq $filter_chr;
		}
		
		say "  $chr";
		
		open(IN,"<","$dir_in/$cell.$chr.rpkm") or die("$dir_in/$cell.$chr.rpkm\n$!");
		open(OUT,">","$dir_out/x_$chr.matrix");
		
		my @temp=split(/\s+/,<IN>);
		die('chr name not match!!!\n') unless $temp[-1] eq $chr;
		$_=<IN>;
		@temp=split;
		# die("H1 shoud be 1114D!!!\n") unless @temp==1114;
		print OUT while <IN>;
	}

}