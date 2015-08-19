#! /usr/bin/perl -w
use strict;
use warnings;
use Term::ANSIColor; 

#-------------------------------------------------------------------------------
# What: This program will remove the internal stop codons (TAA, TAG, TGA)
#		and replace with gaps (---) from the nucleotide alignment. 
#		Typically, this would be done to calculate dN/dS
#
# Need: - multiple FASTA formatted sequence file
#
# Note: This assumes prior filtering has been done so that:
#			-the human reference sequence is a conserved open reading frame.
#			-each nucleotide sequence is divisible by three evenly.
#			
#-------------------------------------------------------------------------------
# 
my $usage = "usage = ReplaceStopCodonWithGaps_onDir.pl [input_dir] [output_dir]\n";
die $usage unless @ARGV == 2;
#--------------------------------------------------------------------------------
#
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
	my $output_file = "$output_dir"."$filename"."_noStopCodon.fasta";
	open(OUT, ">$output_file");
	
	open (MULTI_FASTA, "<$input_file") or die "Cannot open $input_file:$!\n";
	my @allspecies = <MULTI_FASTA>;		# Define array
	close(MULTI_FASTA);

# put all sequences into their own files
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

# make an output array to store the processed sequences
my @processed_seq;

my $count = 0;
# loop to search in triplets, identify stop codons, and replace stop codon with gaps
foreach $seq (@seq){

	my @check_ref = split(//,$seq);
	my $length = $#check_ref + 1;

	my $e = 0;
	my $f = 1;
	my $g = 2;
	
	while ($g < $length){
	
		# search in triplets 
		my $codon = "$check_ref[$e]$check_ref[$f]$check_ref[$g]";
		# identify stop codons
		if(($codon eq "TAA")|($codon eq "TAG")|($codon eq "TGA")|($codon eq "tAG")|($codon eq "tAA")|($codon eq "tGA")|($codon eq "TaG")|($codon eq "TaA")|($codon eq "TgA")|($codon eq "TAg")|($codon eq "TAa")|($codon eq "TGa")|($codon eq "taG")|($codon eq "taA")|($codon eq "tgA")|($codon eq "tAg")|($codon eq "tAa")|($codon eq "tGa")|($codon eq "Tga")|($codon eq "Taa")|($codon eq "Tag")|($codon eq "tga")|($codon eq "taa")|($codon eq "tag")){            
			#replace the stop codons with gaps			
			$check_ref[$e] = "-";
			$check_ref[$f] = "-";
			$check_ref[$g] = "-";

			$e = $e + 3;
			$f = $f + 3;
			$g = $g + 3;	
		
		}
		else{
			$e = $e + 3;
			$f = $f + 3;
			$g = $g + 3;	
		
			next;
		}
	}

	# merge the processed nucleotides back into one string per sequence
	my $merged_seq = "";
	foreach my $check_ref (@check_ref){
		$merged_seq = join("", $merged_seq,$check_ref);
	}
	
	# now make an array that adds all of the merged sequences together
	push(@processed_seq, $merged_seq);
	
	#print ("$count\n");
	$count++;
}		 
	
# print the new sequences without stop codons
my $replace = $#names + 1;
my $final = 0;
while ($final < $replace){
	print OUT "$names[$final]";
	print OUT "$processed_seq[$final]\n";
	$final++;
}
}
