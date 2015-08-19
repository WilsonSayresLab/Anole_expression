#! /usr/bin/perl -w
use strict;
use warnings;
use Term::ANSIColor; 

#----------------------------------------------------------------------------------------
# What: This program calls the the kaks calculator which will calculate nonsynonymous (Ka) and
#   synonymous (Ks) substitution rates through model selection and model averaging.
#
# Need: axt files made from fasta files
#----------------------------------------------------------------------------------------
#
my $usage = "perl 10_kakscall.pl [path_kaks_calculator] [input_dir] [output_dir]\n";
die $usage unless @ARGV == 3;
#----------------------------------------------------------------------------------------
my $kaks = "$ARGV[0]";
my $input_dir = "$ARGV[1]";
my $output_dir = "$ARGV[2]";

#Access the input directory
my $file;
opendir(DIR, $input_dir) || die "can't opendir $input_dir: $!";
while (defined($file = readdir(DIR))){
	next if $file =~ /^\.\.?$/;		# skip . and ..
	next if $file =~ /^\./;
	
	#define the input file
	my $input_file = "$input_dir"."$file";		# input file

	print "$input_file\n";
	my $output_file = "$output_dir"."$file".".kaks";		# input file
	
	system(`$kaks -i $input_file -o $output_file -m NG`);
}
