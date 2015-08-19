#! /usr/bin/perl -w
use strict;
use warnings;
use Term::ANSIColor;

#-------------------------------------------------------------------------------
# What: This program will check a multiple FASTA file to see that the open
# 		reading frame is conserved in each species, relative to chicken. If not
#		it will correct it.
#
# Need: - multiple FASTA formatted sequence file
#		- check the reference species name
#-------------------------------------------------------------------------------
# 
my $usage = "usage = fix_frame.pl [input_dir] [output_dir]\n";
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
	my $output_file = "$output_dir"."$filename"."_fixFrame.fasta";
	open(OUT, ">$output_file");
	
	open (MULTI_FASTA, "<$input_file") or die "Cannot open $input_file:$!\n";
	my @allspecies = <MULTI_FASTA>;		# Define array
	close(MULTI_FASTA);

	#print "$allspecies[1]\n";
	
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

	my @filtered;

	foreach my $seq (@seq){
		@filtered = "";

		my @check_ref = split(//,$seq);
		my $length = $#check_ref + 1;

		my $e = 0;
		my $f = 1;
		my $g = 2;

		# search in triplets
		while ($e < $length){
			#print "$check_ref[$e]$check_ref[$f]$check_ref[$g]\n";
			if(($check_ref[$e] =~ /N|A|T|G|C|n|a|t|g|c/)&&($check_ref[$f] =~ /N|A|T|G|C|n|a|t|g|c/)&&($check_ref[$g] =~ /N|A|T|G|C|n|a|t|g|c/)){
				#the sequence is in frame, so do nothing
				push(@filtered, $check_ref[$e]);
				push(@filtered, $check_ref[$f]);
				push(@filtered, $check_ref[$g]);
				$e = $e + 3;	$f = $f + 3;	$g = $g + 3;
			}
			elsif(($check_ref[$e] =~ /-/)&&($check_ref[$f] =~ /-/)&&($check_ref[$g] =~ /-/)){
				#the gap is in frame, so do nothing
				push(@filtered, "-");
				push(@filtered, "-");
				push(@filtered, "-");
				$e = $e + 3;	$f = $f + 3;	$g = $g + 3;
			}
			elsif(($check_ref[$e] =~ /-/)&&($check_ref[$f] =~ /N|A|T|G|C|n|a|t|g|c/)&&($check_ref[$g] =~ /N|A|T|G|C|n|a|t|g|c/)){	
				#Replace this trio with all gaps
				push(@filtered, "-");
				push(@filtered, "-");
				push(@filtered, "-");
				$e = $e + 3;	$f = $f + 3;	$g = $g + 3;
			}
			elsif(($check_ref[$e] =~ /N|A|T|G|C|n|a|t|g|c/)&&($check_ref[$f] =~ /-/)&&($check_ref[$g] =~ /N|A|T|G|C|n|a|t|g|c/)){	
				#Replace this trio with all gaps
				push(@filtered, "-");
				push(@filtered, "-");
				push(@filtered, "-");
				$e = $e + 3;	$f = $f + 3;	$g = $g + 3;
			}
			elsif(($check_ref[$e] =~ /N|A|T|G|C|n|a|t|g|c/)&&($check_ref[$f] =~ /N|A|T|G|C|n|a|t|g|c/)&&($check_ref[$g] =~ /-/)){	
				#Replace this trio with all gaps
				push(@filtered, "-");
				push(@filtered, "-");
				push(@filtered, "-");
				$e = $e + 3;	$f = $f + 3;	$g = $g + 3;
			}
			elsif(($check_ref[$e] =~ /-/)&&($check_ref[$f] =~ /-/)&&($check_ref[$g] =~ /N|A|T|G|C|n|a|t|g|c/)){	
				#Replace this trio with all gaps
				push(@filtered, "-");
				push(@filtered, "-");
				push(@filtered, "-");
				$e = $e + 3;	$f = $f + 3;	$g = $g + 3;
			}
			elsif(($check_ref[$e] =~ /-/)&&($check_ref[$f] =~ /N|A|T|G|C|n|a|t|g|c/)&&($check_ref[$g] =~ /-/)){	
				#Replace this trio with all gaps
				push(@filtered, "-");
				push(@filtered, "-");
				push(@filtered, "-");
				$e = $e + 3;	$f = $f + 3;	$g = $g + 3;
			}
			elsif(($check_ref[$e] =~ /N|A|T|G|C|n|a|t|g|c/)&&($check_ref[$f] =~ /-/)&&($check_ref[$g] =~ /-/)){	
				#Replace this trio with all gaps
				push(@filtered, "-");
				push(@filtered, "-");
				push(@filtered, "-");
				$e = $e + 3;	$f = $f + 3;	$g = $g + 3;
			}
			elsif(($check_ref[$e] =~ /N|A|T|G|C|n|a|t|g|c/)&&($check_ref[$f] =~ /N|A|T|G|C|n|a|t|g|c/)){	
				#We are at the end of the sequence, and don't need to do anything
				push(@filtered, $check_ref[$e]);
				push(@filtered, $check_ref[$f]);
				$e = $e + 3;	$f = $f + 3;	$g = $g + 3;
			}
			elsif(($check_ref[$e] =~ /-/)&&($check_ref[$f] =~ /N|A|T|G|C|n|a|t|g|c/)){	
				#Replace this duo with all gaps
				push(@filtered, "-");
				push(@filtered, "-");
				$e = $e + 3;	$f = $f + 3;	$g = $g + 3;
			}
			elsif(($check_ref[$e] =~ /N|A|T|G|C|n|a|t|g|c/)&&($check_ref[$f] =~ /-/)){	
				#Replace this duo with all gaps
				push(@filtered, "-");
				push(@filtered, "-");
				$e = $e + 3;	$f = $f + 3;	$g = $g + 3;
			}
			elsif(($check_ref[$e] =~ /-/)&&($check_ref[$f] =~ /-/)){	
				#Replace this duo with all gaps
				push(@filtered, "-");
				push(@filtered, "-");
				$e = $e + 3;	$f = $f + 3;	$g = $g + 3;
			}
			elsif($check_ref[$e] =~ /N|A|T|G|C|n|a|t|g|c/){	
				#We are at the end of the sequence, and don't need to do anything
				push(@filtered, $check_ref[$e]);
				$e = $e + 3;	$f = $f + 3;	$g = $g + 3;
			}
			elsif($check_ref[$e] =~ /-/){	
				#Replace these columns with a placeholder to be removed later
				push(@filtered, "-");
				$e = $e + 3;	$f = $f + 3;	$g = $g + 3;
			}
			else{
				#do nothing
				$e = $e + 3;	$f = $f + 3;	$g = $g + 3;
			}
	
		}
		$seq = join('', @filtered)
	}

	#Print the species names and filtered sequences
	my $get_length = $#names + 1;
	my $final = 0;
	while($final < $get_length){
		#push (@allspecies, "$names[$final]");
		#push (@allspecies, "$seq[$final]");
		print OUT "$names[$final]";
		print OUT "$seq[$final]\n";
		$final++;
	}
}

