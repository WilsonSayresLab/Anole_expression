#! /usr/bin/perl -w
use strict;
use warnings;

#----------------------------------------------------------------------------------------
# What: This program will read through a directory that contains 
#		aligned multiple FASTA files and remove the NM identifier 
#		(FASTA header) and all the trailing information (end lines)
#		for each gene/protein.  
#
# Need: multiple FASTA formatted sequence file
#----------------------------------------------------------------------------------------
#
my $usage = "usage = removeheader_endlines_onDir.pl [input_dir] [output_dir]\n",
#die $usage unless @ARGV == 2;
#----------------------------------------------------------------------------------------
my $file;
my $input_dir = "$ARGV[0]";
my $output_dir = "$ARGV[1]";

opendir(DIR, $input_dir) || die "can't opendir $input_dir: $!";
while (defined($file = readdir(DIR))){
	next if $file =~ /^\.\.?$/;		# skip . and ..
	next if $file =~ /^\./;

	my $input_file = "$input_dir"."$file";		# input file
	my $output_file = "$output_dir"."$file";		#output file

	open(INPUT, $input_file) or die("Could not open $input_file");
	open(OUT, ">$output_file") or die("Could not open $output_file");

	my $seq; 
	foreach $seq (<INPUT>)  {   
		$seq=~ s/ENS\S+\d+_//;
		$seq =~ s/\s\d+\s+\S+\d+\:\d+\S\d+\S//;	
			#remove trailing information
	
		print OUT "$seq";    
	} 
 
	close(INPUT);
	close(OUT);
}
