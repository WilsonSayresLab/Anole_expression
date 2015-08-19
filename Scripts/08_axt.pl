#! /usr/bin/perl -w
use strict;
use warnings;
use Term::ANSIColor; 

#----------------------------------------------------------------------------------------
# What: This program will convert fasta files to axt files
#
# Need: multiple FASTA formatted sequence file
#----------------------------------------------------------------------------------------
#
my $usage = "perl axt.pl [path_to_fasta_to_AXT_converter] [input_dir]\n";
die $usage unless @ARGV == 2;
#----------------------------------------------------------------------------------------
#
my $AXT = "$ARGV[0]";
my $input_dir = "$ARGV[1]";
#my $output_dir = "$ARGV[2]";

#Access the input directory
my $file;
opendir(DIR, $input_dir) || die "can't opendir $input_dir: $!";
while (defined($file = readdir(DIR))){
	next if $file =~ /^\.\.?$/;		# skip . and ..
	next if $file =~ /^\./;
	
	#define the input file
	my $input_file = "$input_dir"."$file";		# input file
	
#	print "$input_file\n";
#		#define the output file
#	my @get_filename = split(/\./,$file);		# output file
#	my $filename = $get_filename[0];
#	my $output_file = "$output_dir"."$filename".".kaks.txt";
	
	system(`perl $AXT $input_file`);
}
