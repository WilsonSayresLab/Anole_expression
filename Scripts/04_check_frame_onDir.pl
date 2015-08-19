#! /usr/bin/perl -w
use strict;
use warnings;
use Term::ANSIColor;

#-------------------------------------------------------------------------------
# What: This program will check a multiple FASTA file in a directory 
#		to see that the open reading frame is conserved in each species, 
#		relative to reference species, human. If not it will correct it.
#		This works only for insertions in non-reference species relative the reference.
#		This program assumes that the reference species is in frame, so any
#		non-triplet gaps in the reference will be removed from the alignments.
#		This code is intended only for use with coding nucleotide alignments.
#
# For example it wil change this:
#
#		hg19    ATT-TCATAG
#		gorGor1 ATTTTCATAG
#
# To this:
#
#		hg19    ATTTCATAG
#		gorGor1 ATTTCATAG
#
# But this:
#
#		hg 19   ATT---TCATAG
#		gorGor1 ATTCTTTCATAG
#
# will not be removed because it doen't change the frame of the reference (human, hg19)
# 
# Need: - multiple FASTA formatted sequence file
#		- check the reference species name
#-------------------------------------------------------------------------------
my $usage = "usage = perl check_frame_onDir.pl [input_dir] [output_dir] [ref_species]\n";
die $usage unless @ARGV == 3; 
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
	
	my @get_filename = split(/\./,$file);		# output file
	my $filename = $get_filename[0];
	my $output_file = "$output_dir"."$filename"."_checkFrame.fasta";
	open(OUT, ">$output_file");
	
	open (MULTI_FASTA, "<$input_file") or die "Cannot open $input_file:$!\n";
	my @allspecies = <MULTI_FASTA>;		# Define array
	close(MULTI_FASTA);
	
	#put all sequences into their own arrays
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

	@allspecies="";	#clear the array, to save memory	

	#Extract the name and sequence of the reference species
	my $ref_species = $ARGV[2];

	system(`grep -A1 "$ref_species" $input_file > tmp.txt`);

	open (TMP, "<tmp.txt") or die "Cannot open tmp.txt:$!\n";
	my @tmp = <TMP>;	my $tmp;	# Define array
	close(TMP);

	system(`rm tmp.txt`);
	
	#print "$tmp[1]\n";

	#check to see if the reference sequence is a multiple of three 
	my @check_ref = split(//,$tmp[1]);
	#print ("@check_ref\n");
	my $length = $#check_ref + 1;
	#print "$length\n";

	my $e = 0;
	my $f = 1;
	my $g = 2;
	# search in triplets
	while ($e < $length){
		if(($check_ref[$e] =~ /N|A|T|G|C|n|a|t|g|c/)&&($check_ref[$f] =~ /N|A|T|G|C|n|a|t|g|c/)&&($check_ref[$g] =~ /N|A|T|G|C|n|a|t|g|c/)){
		#the sequence is in frame, so do nothing
		$e = $e + 3;	$f = $f + 3;	$g = $g + 3;
		}
		elsif(($check_ref[$e] =~ /-/)&&($check_ref[$f] =~ /-/)&&($check_ref[$g] =~ /-/)){
		#the gap is in frame, so do nothing
		$e = $e + 3;	$f = $f + 3;	$g = $g + 3;
		}
		elsif(($check_ref[$e] =~ /-/)&&($check_ref[$f] =~ /N|A|T|G|C|n|a|t|g|c/)&&($check_ref[$g] =~ /N|A|T|G|C|n|a|t|g|c/)){	
		#Replace this column with a placeholder to be removed later
		foreach $seq (@seq){
			my @findseq = split(//, $seq);
			$findseq[$e] = "X";
			$seq = join('', @findseq)
			}
			$e = $e + 3;	$f = $f + 3;	$g = $g + 3;
		}
		elsif(($check_ref[$e] =~ /N|A|T|G|C|n|a|t|g|c/)&&($check_ref[$f] =~ /-/)&&($check_ref[$g] =~ /N|A|T|G|C|n|a|t|g|c/)){	
			#Replace this column with a placeholder to be removed later
			foreach $seq (@seq){
				my @findseq = split(//, $seq);
				$findseq[$f] = "X";
				$seq = join('', @findseq)
			}
			$e = $e + 3;	$f = $f + 3;	$g = $g + 3;
		}
		elsif(($check_ref[$e] =~ /N|A|T|G|C|n|a|t|g|c/)&&($check_ref[$f] =~ /N|A|T|G|C|n|a|t|g|c/)&&($check_ref[$g] =~ /-/)){	
			#Replace this column with a placeholder to be removed later
			foreach $seq (@seq){
				my @findseq = split(//, $seq);
				$findseq[$g] = "X";
				$seq = join('', @findseq)
			}
			$e = $e + 3;	$f = $f + 3;	$g = $g + 3;
		}
		elsif(($check_ref[$e] =~ /-/)&&($check_ref[$f] =~ /-/)&&($check_ref[$g] =~ /N|A|T|G|C|n|a|t|g|c/)){	
			#Replace these columns with a placeholder to be removed later
			foreach $seq (@seq){
				my @findseq = split(//, $seq);
				$findseq[$e] = "X";
				$findseq[$f] = "X";
				$seq = join('', @findseq)
			}
			$e = $e + 3;	$f = $f + 3;	$g = $g + 3;
		}
		elsif(($check_ref[$e] =~ /-/)&&($check_ref[$f] =~ /N|A|T|G|C|n|a|t|g|c/)&&($check_ref[$g] =~ /-/)){	
			#Replace these columns with a placeholder to be removed later
			foreach $seq (@seq){
				my @findseq = split(//, $seq);
				$findseq[$e] = "X";
				$findseq[$g] = "X";
				$seq = join('', @findseq)
			}
			$e = $e + 3;	$f = $f + 3;	$g = $g + 3;
		}
		elsif(($check_ref[$e] =~ /N|A|T|G|C|n|a|t|g|c/)&&($check_ref[$f] =~ /-/)&&($check_ref[$g] =~ /-/)){	
			#Replace these columns with a placeholder to be removed later
			foreach $seq (@seq){
				my @findseq = split(//, $seq);
				$findseq[$f] = "X";
				$findseq[$g] = "X";
				$seq = join('', @findseq)
			}
			$e = $e + 3;	$f = $f + 3;	$g = $g + 3;
		}
		elsif(($check_ref[$e] =~ /N|A|T|G|C|n|a|t|g|c/)&&($check_ref[$f] =~ /N|A|T|G|C|n|a|t|g|c/)){	
			#We are at the end of the sequence, and don't need to do anything
			$e = $e + 3;	$f = $f + 3;	$g = $g + 3;
		}
		elsif(($check_ref[$e] =~ /-/)&&($check_ref[$f] =~ /N|A|T|G|C|n|a|t|g|c/)){	
			#Replace this column with a placeholder to be removed later
			foreach $seq (@seq){
				my @findseq = split(//, $seq);
				$findseq[$e] = "X";
				$seq = join('', @findseq)
			}
			$e = $e + 3;	$f = $f + 3;	$g = $g + 3;
		}
		elsif(($check_ref[$e] =~ /N|A|T|G|C|n|a|t|g|c/)&&($check_ref[$f] =~ /-/)){	
			#Replace this column with a placeholder to be removed later
			foreach $seq (@seq){
				my @findseq = split(//, $seq);
				$findseq[$f] = "X";
				$seq = join('', @findseq)
			}
			$e = $e + 3;	$f = $f + 3;	$g = $g + 3;
		}
		elsif(($check_ref[$e] =~ /-/)&&($check_ref[$f] =~ /-/)){	
			#Replace these columns with a placeholder to be removed later
			foreach $seq (@seq){
				my @findseq = split(//, $seq);
				$findseq[$e] = "X";
				$findseq[$f] = "X";
				$seq = join('', @findseq)
			}
			$e = $e + 3;	$f = $f + 3;	$g = $g + 3;
		}
		elsif($check_ref[$e] =~ /N|A|T|G|C|n|a|t|g|c/){	
			#We are at the end of the sequence, and don't need to do anything
			$e = $e + 3;	$f = $f + 3;	$g = $g + 3;
		}
		elsif($check_ref[$e] =~ /-/){	
			#Replace these columns with a placeholder to be removed later
			foreach $seq (@seq){
				my @findseq = split(//, $seq);
				$findseq[$e] = "X";
				$seq = join('', @findseq)
			}
			$e = $e + 3;	$f = $f + 3;	$g = $g + 3;
		}
		else{
			#do nothing
			$e = $e + 3;	$f = $f + 3;	$g = $g + 3;
		}
	}

	# Remove all placeholders
	foreach $seq (@seq){
		my @filtered;
		my @findseq = split(//, $seq);
		foreach my $findseq (@findseq){
			if ($findseq ne "X"){
				push (@filtered, "$findseq");
			}
		}
		$seq = join('', @filtered);
	}

	#Print the species names and filtered sequences
	my $get_length = $#names + 1;
	my $final = 0;
	while($final < $get_length){
		#push (@allspecies, "$names[$final]");
		#push (@allspecies, "$seq[$final]");
		print OUT "$names[$final]";
		print OUT "$seq[$final]";
		$final++;
	}
}	
