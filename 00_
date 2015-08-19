###############################################################
#-------------- OverView ----------------
#
# Codes/programs needed to run this pipeline
#
#---Publicly available programs-------------------------
#
#		parseFastaIntoAXT.pl 	https://code.google.com/p/kaks-calculator/downloads/detail?name=parseFastaIntoAXT.pl
#		KaKs_Calculator1.2		https://code.google.com/p/kaks-calculator/downloads/detail?name=KaKs_Calculator1.2.tar.gz&can=2&q=
#
#---In house programs-------------------------
#
#	/AnoleChickeSubstitution/Scripts
#		01_splitFastaFiles.pl
#		02_removeheader_endlines_onDir.pl
#		03_changeNames_onDir.pl
#		04_check_frame_onDir.pl
#		05_fix_frame_onDir.pl
#		06_count_bases_percent_onDir.pl
#		07_ReplaceStopCodonWithGaps_onDir.pl
#		08_axt.pl
#		09_parseFastaIntoAXT.pl
#		10_kakscall.pl
#		11_kaks.txt.sh
#
#############################################################

######## KaKs README #######

#--------------------------
# 1. Prepare directories 
#--------------------------

# make directories
cd 				# move to home directory or location you would like to run the analysis
cd Desktop
mkdir AnoleChickenSubstitution
cd AnoleChickenSubstitution
mkdir Scripts 				# locally written perl scripts
mkdir Input					# downloaded fasta files 
mkdir Output				# output data
cd Output					# create output folder for the following positions:
mkdir AAWZ02037698
mkdir AAWZ02040127
mkdir AAWZ02041607
mkdir chr1
mkdir chr2
mkdir chr3
mkdir chr4
mkdir chr5
mkdir chr6
mkdir GL343282
mkdir GL343338
mkdir GL343364
mkdir GL343417
mkdir GL343422
mkdir GL343423
mkdir GL343439
mkdir GL343462
mkdir GL343516
mkdir GL343525
mkdir GL343550
mkdir GL343588
mkdir GL343731
mkdir GL343913
mkdir GL343947
mkdir GL344042
mkdir GL344393
mkdir GL344496
mkdir GL344539
mkdir LGb

#--------------------------
# 2. prepare paths
#--------------------------

kaks=path_to_kaks_calculator 			# /Users/olneykimberly/Desktop/KaKs_Calculator1.2/bin/Mac/KaKs_Calculator

#--------------------------
# 3. Download alignments to chicken from UCSC genome browser
#--------------------------

# https://genome.ucsc.edu 
#-> click on Table Browser 

# clade: "Vertebrate"	genome: "Lizard"	assembly: "May 2010 (Broad AnoCar2.0/anoCar2)"
# group: "Genes and Gene Predictions" 	track: "Ensemble Genes"
# table: ensGene
# region: position "chr1:1-263,920,457" ## will need to be done for all of the positions 
# output format "CDS FASTA alignment from multiple alignment"
# output file : "Chr1_gg.fasta" ## will need to done for all of the positions 
# file type returned "plain text"
#-> click "get output"

# Formatting options "Show nucleotides"
# Species selection "chicken" # unselected other species so that only chicken is selected
#-> click "get output" # Download will begin

# Move downloaded fasta files to your working directory 
mv Downloads/chr1_gg.fasta AnoleChickenSubstitutions/Input/

# Postions:
chrUn_AAWZ02037698:1-14,962
chrUn_AAWZ02040127:1-7,362
chrUn_AAWZ02041607:1-5,474
chr1:1-263,920,457
chr2:1-199,619,895
chr3:1-204,416,410 
chr4:1-156,502,444
chr5:1-150,641,573
chr6:1-80,741,955
chrUn_GL343282:1-1779868
chrUn_GL343338:1-1,258,094
chrUn_GL343364:1-1,083,274
chrUn_GL343417:1-831,895
chrUn_GL343422:1-839,908
chrUn_GL343423:1-834,740
chrUn_GL343439:1-810,073
chrUn_GL343462:1-693,399
chrUn_GL343516:1-504,621
chrUn_GL343525:1-498,372
chrUn_GL343550:1-526,944
chrUn_GL343588:1-426,171
chrUn_GL343731:1-266,425
chrUn_GL343913:1-147,151
chrUn_GL343947:1-117,443
chrUn_GL344042:1-95,694
chrUn_GL344393:1-41,178
chrUn_GL344496:1-31,185
chrUn_GL344539:1-36,94
chrLGb:1-3,271,537

#--------------------------
# 4. Split fasta file 
#--------------------------

# The output fasta files downloaded in the pervious step contains multiple genes per chromosome
# will need to split files into one file (one alignment) per gene:

#cd # move to the output directory each time for the different positions
	# Output/Chr1/1  # each time you run the 01_splitFastaFiles.pl for the different positions to keep them well organized
		
	perl 01_splitFastaFiles.pl <multiple FASTA file>

# example:
# cd Output/Chr1/1		
# perl Scripts/01_splitFastaFiles.pl AnoleChickenSunstitutions/Input/Chr1_gg.fasta 
	
#--------------------------
# 5. Remove header information except for genome build name  # will remove the gene id following the species id
#--------------------------	
		
	perl 02_removeheader_endlines_onDir.pl [input_directory] [output_directory]

# example:
# perl Scripts/02_removeheader_endlines_onDir.pl Output/Chr1/1/ Output/Chr1/2/

#--------------------------
# 6. Make a common name for all files in directory
#--------------------------	

	perl 03_changeNames_onDir.pl [input_directory] [output_directory]

# example:
# perl Scripts/03_changeNames_onDir.pl Output/Chr1/2/ Output/Chr1/3/

#--------------------------
# 7. Check nucleotide frame relative to green_anole
#--------------------------	

# This program will check a multiple FASTA file to see that the open reading frame is conserved in each species, relative to reference species, chicken. 
# If not it will correct it.This works only for insertions in non-reference species relative the reference.This program assumes that the reference species is in frame, so any non-triplet gaps in the reference will be removed from the alignments.

	perl 04_check_fram_onDir.pl [input_directory] [output_directory] green_anole

# example:
# perl Scripts/04_check_fram_onDir.pl Output/Chr1/3/ Output/Chr1/4/ green_anole

#--------------------------
# 8. Check the frame of the nucleotide sequences relative to itself
#--------------------------	

# program will check a multiple FASTA file to see that the open reading frame of the nucleotide sequence is conserved in each species, relative to itself, not to a reference species

	perl 05_fix_frame_onDir.pl [input_directory] [output_directory]

# example:
# perl Scripts/05_fix_frame_onDir.pl Output/Chr1/4/ Output/Chr1/5/

#--------------------------
# 9. Check to see that each species has at least 50% of its nucleotide sequences
#--------------------------	

# This program will check a multiple FASTA file to see that each species has at least 50% of its nucleotide sequence. 
# If not, it will remove that sequence. The output files include the FASTA sequence file and excluded species file.

	perl 06_count_base_percent.pl [input_directory] [output_directory]

# example: 
# perl Scripts/06_count_base_percent.pl Output/Chr1/5/ Output/Chr1/6/

#--------------------------
# 10. Remove stop codons for downstream analysis
#--------------------------	

# This program will remove the internal and terminate stop codons (TAA, TAG, TGA) and replace with gaps (—) from the nucleotide alignment.  
# Typically, this would be done to calculate the substitution rate value, dN/dS.

	perl 07_ReplaceStopCodonWithGaps_onDir.pl [input_directory] [output_directory]

# example:
# perl Scripts/07_ReplaceStopCodonWithGaps_onDir.pl Output/Chr1/6/ Output/Chr1/7/

#--------------------------
# 11. Parse fasta files into AXT files 
#--------------------------	

	perl 08_axt.pl 09_parseFastaIntoAXT.pl [input_directory]

# example: 
# perl 08_axt.pl 09_parseFastaIntoAXT.pl path_to_input_directory/Chr1/7/ 

#--------------------------
# 12. Compute substitution rates nonsynonymous (Ka) and synonymous (Ks) 
#--------------------------	

# KaKs_Calculator is a program that calculates nonsynonymous (Ka) and synonymous (Ks) substitution rates through model selection and model averaging.

	perl 10_kakscall.pl path_to_KaKs_Calculator [input_directory] [output_directory]

# exmple:
# perl 10_kakscall.pl Desktop/KaKs_Calculator1.2/bin/Mac/KaKs_Calculator Output/Chr1/7/ Output/Chr1/8/

#--------------------------
# 13. Combine output data into one file for contigs, autosomes, and scaffolds 
#--------------------------

# will output text files for the contigs, autosomes, and scaffolds in their respective folder "position"/8 and into the output/KaKs/ folder
# will output .csv files for the autosomes, LGb, and HypX into output/KaKs/ 

	sh Scripts/11_kaks.txt.sh 

#############################################################
######## Analysis KaKs values README #######

#--------------------------
# 14. Remove nan 
#--------------------------

# open the autosomes.csv, LGb.csv, and HypX.csv files that was created in the pervious step. 
# remove rows that have "nan" for Ka, Ks, and KaKs columns. Keep rows that "NA" in these columns

#--------------------------
# 15. Replace 0 with NA
#--------------------------

# open the autosomes.csv, LGb.csv, and HypX.csv files 
# sort the Ka.Ks column in ascending order and select expand all. 
# this will pull the rows that have a 0 in the KaKs column to the top. 
# find the rows that have an "NA" in the Ka column and "0" in the KaKs column
# replace the 0 with NA in the KaKs column, this is done so that you don't get a false low value for the KaKs 

#--------------------------
# 16. Find the mean, median, and P-value for the putative X and the autosomes
#--------------------------

# will print the mean, median, and P-value for the putative X and the autosomes to the screen 

perl perl 12_p.pl [input_file_X] [input_file_A]\n 

# example:
# perl perl 12_p.pl Output/KaKs/HypX.txt Output/KaKs/autosomes.txt 

#--------------------------
# 17. Check median and mean in R 
#--------------------------

# open R and run the R script KaKsMedian.R
# this will calculate the median and mean for autosomes, contigs, and scaffolds
# will also create a boxplot comparing the autosomes, LGb, and putative X

#############################################################
	
