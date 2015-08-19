#!/usr/bin/perl -w

# This will take the aligned multiple FASTA file with multiple genes/proteins 
# and create individual FASTA alignment files for each gene/protein

$usage = "splitFastaFiles.pl [multiple genes/proteins in FASTA format]";
die $usage unless @ARGV == 1;

# This variable makes perl read lines that end with "\n\n>" instead of a newline \n
$/="\n\n>";

while (<>) { # foreach line in the input files
	if(/^\s*(\S+)/) { # grab the first word of text
		my @filename = split(/_hg19/,$1);
		my $file = $filename[0];
		open(F,">$file") || # open a file named that word
		warn "$1 write failed:$!\n";
		chomp; # strip off the > at the end
		$temp = ">". $_;
		print F "$temp";
	}
}
