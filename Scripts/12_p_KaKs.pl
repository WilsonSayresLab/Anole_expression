#! /usr/bin/perl -w
use strict;
use warnings;
use Term::ANSIColor;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);

#-------------------------------------------------------------------------------
# What: This program will read through a txt file of KaKs values and
#		compute the weighted averages and median & mean and then compute the 95% CI
#       Will also compute the empirical p-value for the median and mean
#       To read through the Ka or Ks value independtly change the array to 2 or 3 respectively
#
# Need: Define chromosome X and autosomes
#		Define number of replicates
#
#-------------------------------------------------------------------------------

my $usage = "usage = perl 12_p.pl [input_file_X] [input_file_A]\n"; # X-Chromosome, A-autosomes
die $usage unless @ARGV == 2;

#--------------------------------------------------------------------------------

# -------------------------------------
# Define how many replicates
# -------------------------------------

my $replicates = 10000;

# -------------------------------------
# Input files
# -------------------------------------

# putative X
my $input_file_X = "$ARGV[0]";
open(FILE, "<$input_file_X") or die "Cannot open $input_file_X :$!";
my @lineX = <FILE>;
close(FILE);
shift(@lineX); #remove the header line
my $lengthX = @lineX;	#total number of sites in the file
#print "lengthX = $lengthX\n"; #print the number of lines

# autosomes 1-6
my $input_file_A = "$ARGV[1]";
open(FILE, "<$input_file_A") or die "Cannot open $input_file_A :$!";
my @lineA = <FILE>;
close(FILE);
shift(@lineA); #remove the header line
my $lengthA = @lineA;	#total number of sites in the file
#print "lengthA = $lengthA\n"; #print the number of lines

# -------------------------------------
# Open arrays for Ka/Ks values for X and A
# -------------------------------------

my @X;
my @A;

# Push Ka/Ks values for chromosome X into the X array
my $c_X = 0; # chromosome X
foreach my $lineX (@lineX){
    my @arrayX = split(/\t/,$lineX);
    push(@X, $arrayX[4]);
    $c_X++;
}

# Push Ka/Ks values for chromosome A into the A array
my $c_A = 0; # chromosomes 1-6
foreach my $lineA (@lineA){
    my @arrayA = split(/\t/,$lineA);
    push(@A, $arrayA[4]);
    $c_A++;
}

# -------------------------------------
# Compute mean and median values
# -------------------------------------

my $sum_X = 0;
my $sum_A = 0;

my $length_X;
my $length_A;

@X = sort(@X);
@A = sort(@A);

foreach my $X (@X){
    if ($X !~ "NA"){
        $sum_X = $sum_X+ $X;
        $length_X++;
    }
}
foreach my $A (@A){
    if ($A !~ "NA"){
        $sum_A = $sum_A + $A;
        $length_A++;
    }
}
print "length_X = $length_X\n";
print "length_A = $length_A\n";

# compute the mean
my $mean_X = $sum_X/$length_X;
my $mean_A = $sum_A/$length_A;

# sanity check for means
print "X mean:\t$mean_X\n";
print "A mean:\t$mean_A\n";

my $median_X;
my $mid_X = int $length_X/2;
if(@X % 2) {
    $median_X = $X[$mid_X];
} else {
    $median_X = ($X[$mid_X-1] + $X[$mid_X])/2;
}

my $median_A;
my $mid_A = int $length_A/2;
if(@A % 2) {
    $median_A = $A[$mid_A];
} else {
    $median_A = ($A[$mid_A-1] + $A[$mid_A])/2;
}

# sanity check for medians
print "X median:\t$median_X\n";
print "A median:\t$median_A\n";

# -------------------------------------
# compute the bootstrap replicates for the mean
# -------------------------------------

# merge X and A
my @merge = @X;
push(@merge, @A);
my $length_merge = @merge;

my $mean_p = 0;
my $r_mn = 0;
while ($r_mn < $replicates){
    @merge = shuffle(@merge);
    my @b_X; 	#define the bootstrap arrays
    my @b_A;
    my $c = 0;
    while ($c < $lengthX){
        push(@b_X, $merge[$c]);
        $c++;
    }
    while($c < $length_merge){
        push(@b_A, $merge[$c]);
        $c++;
    }
    
    @b_X = sort(@b_X);
    @b_A = sort(@b_A);
    
    my $length_b_X;
    my $length_b_A;
    my $sum_b_X = 0;
    my $sum_b_A = 0;
    
    foreach my $b_X (@b_X){
        if ($b_X !~ "NA"){
            $sum_b_X = $sum_b_X+ $b_X;
            $length_b_X++;
        }
    }
    foreach my $b_A (@b_A){
        if ($b_A !~ "NA"){
            $sum_b_A = $sum_b_A + $b_A;
            $length_b_A++;
        }
    }
    
    # compute the mean
    my $mean_b_X = $sum_b_X/$length_b_X;
    my $mean_b_A = $sum_b_A/$length_b_A;
    
    if (($mean_b_X - $mean_b_A)>=($mean_X - $mean_A)){
        $mean_p++;
    }
    
    $r_mn++;
}

# -------------------------------------
# compute the p-value for the mean
# -------------------------------------

my $empirical_mean_p = $mean_p/$replicates;
print "mean_p: $empirical_mean_p\n";

# -------------------------------------
# compute the bootstrap replicates for the median
# -------------------------------------

my $median_p = 0;
my $r_md = 0;
while ($r_md < $replicates){
    @merge = shuffle(@merge);
    my @b_X;
    my @b_A;
    my $c = 0;
    while ($c < $lengthX){
        push(@b_X, $merge[$c]);
        $c++;
    }
    while($c < $length_merge){
        push(@b_A, $merge[$c]);
        $c++;
    }
    
    @b_X = sort(@b_X);
    @b_A = sort(@b_A);
    
    my $length_b_X;
    my $length_b_A;
    
    foreach my $b_X (@b_X){
        if ($b_X !~ "NA"){
            $length_b_X++;
        }
    }
    foreach my $b_A (@b_A){
        if ($b_A !~ "NA"){
            $length_b_A++;
        }
    }
    
    my $median_b_X;
    my $mid_X = int @b_X/2;
    if(@b_X % 2) {
        $median_b_X = $b_X[$mid_X];
    } else {
        $median_b_X = ($b_X[$mid_X-1] + $b_X[$mid_X])/2;
    }
    
    my $median_b_A;
    my $mid_A = int @b_A/2;
    if(@b_A % 2) {
        $median_b_A = $b_A[$mid_A];
    } else {
        $median_b_A = ($b_A[$mid_A-1] + $b_A[$mid_A])/2;
    }
    
    if (($median_b_X - $median_b_A)>=($median_X - $median_A)){
        $median_p++;
    }	
    
    $r_md++;
}

# -------------------------------------
# compute the p-value for the median 
# -------------------------------------

my $empirical_median_p = $median_p/$replicates;
print "median_p: $empirical_median_p\n";
