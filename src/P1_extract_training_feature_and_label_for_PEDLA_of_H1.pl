#!/usr/bin/env perl
use strict;
use warnings;
use 5.010;

###<discription>this program extracts 1114-dimensional features and labels of H1 cell for PEDLA according to the training set. Each enhancher or non-enhancher in the training set is treated as a variable-length sequence that will be used to train PEDLA.

###1114-dimensional input features of the program are in sub-directory "H1_cell_line" of the directory "RPKM/Data_matrix_combined-step200win200.H1" and input labels are in sub-directory "H1_cell_line" of the directory "Training_set/Training_set9"
###output of the program is in the sub-directory "H1_cell_line" of the directory "Data_for_training_and_prediction.H1/Data_for_training_and_prediction-ratio9step200win200". In the output directory there will be about 2*N files, where N is the total number of ehnhancers and non-enhanchers in the training set. The N files prefixed by "x_" are the features and the other N files prefixed by "y_" are corresponding labels.


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


my $dir_tag_count="RPKM/Data_matrix_combined-step${bin_size}win$win.H1";   
opendir(DIR,$dir_tag_count) or die("open $dir_tag_count error!\n$!");
my %cell_with_marker;
$cell_with_marker{$_}=0 foreach sort grep { -d "$dir_tag_count/$_" && !/^\./} readdir(DIR);


foreach my $cell(sort keys %cell_with_marker)
{
	next unless $cell=~/H1_cell_line/;
	
	say "$cell";
	
	my $dir_in="./Training_set/Training_set$num_negative_set/$cell";
	my $dir_in2="$dir_tag_count/$cell";
	my $dir_out="Data_for_training_and_prediction";
	mkdir $dir_out unless -e $dir_out;
	$dir_out.=".H1/$dir_out-ratio${num_negative_set}step${bin_size}win$win";
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
		
		open(IN2,"<","$dir_in2/$cell.$chr.rpkm") or die($!);
		
		
		my @temp=split(/\s+/,<IN2>);
		die('chr name not match!!!\n') unless $temp[-1] eq $chr;
		$_=<IN2>;
		s/[\n\r]+$//;
		my @data;
		push @data,$_ while <IN2>;
		
		open(IN,"<","$dir_in/${flag}_segments.bed") or die($!);
		
		# my @data;
		# my %state;
		my $flag=0;
		my $seq=1;
		my $pre=-1;
		while(<IN>)
		{
			last if !/^$chr\t/ && $flag==1;
			next unless /^$chr\t/;
			$flag=1 unless $flag;
			my @temp=split;
			$temp[-1]=~s/^E//;
			if($pre!=$temp[1])
			{
				open(OUT,">","$dir_out/y_${chr}_seq$seq.matrix");
				open(OUT2,">","$dir_out/x_${chr}_seq$seq.matrix");
				say "    $seq" if $seq%1000==1;
				$seq++;
			}
			for(my $i=int($temp[1]/$bin_size);$i<=int($temp[2]/$bin_size);$i++)
			{
				say OUT $temp[3];
				print OUT2 $data[$i];
				# push @data,$temp[3];
				# $state{$temp[3]}++;
			}
			$pre=$temp[2];
		}

		# my $i=1;
		# $state{$_}=$i++ foreach sort keys %state;
		# say OUT $state{$_} foreach @data;
		#say OUT ("0\t" x ($state{$_}-1))."1".("\t0" x ((keys %state)-$state{$_})) foreach @data;
	}
}