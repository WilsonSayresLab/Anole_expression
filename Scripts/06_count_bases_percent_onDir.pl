#! /usr/bin/perl -w
use strict;
use warnings;
use Term::ANSIColor;

#-------------------------------------------------------------------------------
# What: This program will check a multiple FASTA file to see that each species
# 		has at least 50% of its nucleotide sequence. If not, it will remove that
#		  sequence.
#
# Need: - multiple FASTA formatted sequence file
#		
#-------------------------------------------------------------------------------
# 
my $usage = "usage = count_bases_percent_onDir.pl [input_dir] [output_dir]\n";
die $usage unless @ARGV == 2; 
#-------------------------------------------------------------------------------
my $input_dir = "$ARGV[0]";
my $output_dir = "$ARGV[1]";

#Acess the input directory
my $file;
my $counter = 1;
opendir(DIR, $input_dir) || die "can't opendir $input_dir: $!";
while (defined($file = readdir(DIR))){
	next if $file =~ /^\.\.?$/;		# skip . and ..
	next if $file =~ /^\./;
	
	my $input_file = "$input_dir"."$file";		# input file
	
	my @get_filename = split(/_/,$file);		# output file
	my $filename = $get_filename[0];
	my $output_file = "$output_dir"."$filename"."_percent.fasta";
	open(OUT, ">$output_file");
	
	my $outfile2 = "$output_dir"."$filename"."_excludedspecies.txt";
	open(OUT2, ">$outfile2") or die "Cannot open the output file2 - $outfile2 :$!";
	
	open (MULTI_FASTA, "<$input_file") or die "Cannot open $input_file:$!\n";
	my @allspecies = <MULTI_FASTA>;		# Define array
	close(MULTI_FASTA);

#put all sequences into their own files
my $all_length = $#allspecies;
my @names;		my $names;
my @seq;		my $seq;
my $c = 0;
my $d = 1;
while($c < $all_length){
	push(@names, $allspecies[$c]);
	push(@seq, $allspecies[$d]);
	$c = $c + 2;
	$d = $d + 2;
}
@allspecies="";	#clear the array

#Name of the reference species
#my $ref_species = $ARGV[1];

#system(`grep -A1 "$ref_species" $file > tmp.txt`);

#open (TMP, "<tmp.txt") or die "Cannot open tmp.txt:$!\n";
#my @tmp = <TMP>;		# Define array
#close(TMP);

#system(`rm tmp.txt`);

# check to see if the reference sequence is a multiple of three

my $counter = 0; #counter for @seq, list of sequences
my $seq_length = $#seq + 1; #length of @seq

my $indexnumber = 0; #counter for @sequence, list of base pairs for each sequence from @seq
my $basepairs = 0; #counter for number of base pairs in @sequence
my $totalcount = 0; #counter for the length of each sequence in @sequence

#Loop for each species sequence
for ($counter = 0; $counter < $seq_length; $counter++){
	my @sequence = split(//, $seq[$counter]);
	my $sequence_length = $#sequence + 1;
	
	$basepairs = 0;
	$totalcount = 0;
	
	#Loop for each base pair
	BASE_COUNTER:
	for ($indexnumber = 0; $indexnumber < $sequence_length; $indexnumber++){
		
		#Counts up by one if sequence is a base pair
		if ($sequence[$indexnumber] =~ /A|T|G|C|a|t|g|c|N/){
			$basepairs++;
		}
		
		$totalcount = $indexnumber;
	}
	
	# Prints sequence if at least 50% is retained 
	if ($basepairs / $totalcount >= 0.5){
		print OUT "$names[$counter]";
		print OUT "$seq[$counter]";
	}
	else {
		print OUT2 "$names[$counter]";
	}
}
}
