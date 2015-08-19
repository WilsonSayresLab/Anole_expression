#! /usr/bin/perl -w
use strict;
use warnings;

#----------------------------------------------------------------------------------------
# What: This program will read through a directory that contains 
#		aligned multiple FASTA files and change all the 
#		genome-build-name to common name.
#
# Need: multiple FASTA formatted sequence file
#----------------------------------------------------------------------------------------
#
my $usage = "usage = changeNames_onDir.pl [input_dir] [output_dir]\n",
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
	
	#change genome build name to common name
	my $seq;
	foreach $seq (<INPUT>)	{
		$seq =~ s/hg19/Human/;
		$seq =~ s/panTro4/Chimp/;
		$seq =~ s/gorGor3/Gorilla/;
		$seq =~ s/ponAbe2/Orangutan/;
		$seq =~ s/nomLeu3/Gibbon/;
		$seq =~ s/rheMac3/Rhesus/;
		$seq =~ s/macFas5/Crab_eating_macaque/;
		$seq =~ s/papHam1/Baboon/;
		$seq =~ s/chlSab1/Green_monkey/;
		$seq =~ s/calJac3/Marmoset/;
		$seq =~ s/saiBol1/Squirrel_monkey/;
 		$seq =~ s/otoGar3/Bushbaby/;
 		$seq =~ s/tupChi1/Chinese_tree_shrew/; 
 		$seq =~ s/speTri2/Squirrel/; 
 		$seq =~ s/jacJac1/Lesser_Egyptian_jerboa/;
 		$seq =~ s/micOch1/Prairie_vole/; 
 		$seq =~ s/criGri1/Chinese_hamster/;
		$seq =~ s/mesAur1/Golden_hamster/; 
		$seq =~ s/mm10/Mouse/; 
		$seq =~ s/rn5/Rat/;
		$seq =~ s/hetGla2/Naked_mole_rat/;
		$seq =~ s/cavPor3/Guinea_pig/; 
		$seq =~ s/chiLan1/Chinchilla/; 
		$seq =~ s/octDeg1/Brush_tailed_rat/;
		$seq =~ s/oryCun2/Rabbit/; 
		$seq =~ s/ochPri3/Pika/; 
		$seq =~ s/susScr3/Pig/; 
		$seq =~ s/vicPac2/Alpaca/; 
		$seq =~ s/camFer1/Bactrian_camel/;
		$seq =~ s/turTru2/Dolphin/; 
		$seq =~ s/orcOrc1/Killer_whale/; 
		$seq =~ s/panHod1/Tibetan_antelope/; 
		$seq =~ s/bosTau7/Cow/;
		$seq =~ s/oviAri3/Sheep/; 
		$seq =~ s/capHir1/Domestic_goat/; 
		$seq =~ s/equCab2/Horse/; 
		$seq =~ s/cerSim1/White_rhinoceros/; 
		$seq =~ s/felCat5/Cat/;
		$seq =~ s/canFam3/Dog/;
		$seq =~ s/musFur1/Ferret/; 
		$seq =~ s/ailMel1/Panda/; 
		$seq =~ s/odoRosDiv1/Pacific_walrus/; 
		$seq =~ s/lepWed1/Weddell_seal/; 
		$seq =~ s/pteAle1/Black_flying_fox/;
		$seq =~ s/pteVam1/Megabat/;
		$seq =~ s/myoDav1/Davids_myotis/;
		$seq =~ s/myoLuc2/Microbat/;
		$seq =~ s/eptFus1/Big_brown_bat/; 
		$seq =~ s/eriEur2/Hedgehog/;
		$seq =~ s/sorAra2/Shrew/; 
		$seq =~ s/conCri1/Star_nosed_mole/; 
		$seq =~ s/loxAfr3/Elephant/; 
		$seq =~ s/eleEdw1/Cape_elephant_shrew/; 
		$seq =~ s/triMan1/Manatee/; 
		$seq =~ s/chrAsi1/Cape_golden_mole/; 
		$seq =~ s/echTel2/Tenrec/; 
		$seq =~ s/oryAfe1/Aardvark/; 
		$seq =~ s/dasNov3/Armadillo/;
		$seq =~ s/monDom5/Opossum/;
		$seq =~ s/sarHar1/Tasmanian_devil/;
		$seq =~ s/macEug2/Wallaby/;
		$seq =~ s/ornAna1/Platypus/;
		$seq =~ s/falChe1/Saker_falcon/; 
		$seq =~ s/falPer1/Peregrine_falcon/; 
		$seq =~ s/ficAlb2/Collared_flycatcher/; 
		$seq =~ s/zonAlb1/White_throated_sparrow/; 
		$seq =~ s/geoFor1/Medium_ground_finch/; 
		$seq =~ s/taeGut2/Zebra_finch/; 
		$seq =~ s/pseHum1/Tibetan_ground_jay/;
		$seq =~ s/melUnd1/Budgerigar/; 
		$seq =~ s/amaVit1/Puerto_Rican_parrot/; 
		$seq =~ s/araMac1/Scarlet_macaw/; 
		$seq =~ s/colLiv1/Rock_pigeon/;
		$seq =~ s/anaPla1/Mallard_duck/; 
		$seq =~ s/galGal4/Chicken/; 
		$seq =~ s/galGal3/Chicken/; 
		$seq =~ s/melGal1/Turkey/; 
		$seq =~ s/allMis1/American_alligator/; 
		$seq =~ s/cheMyd1/Green_seaturtle/; 
		$seq =~ s/chrPic1/Painted_turtle/; 
		$seq =~ s/pelSin1/Chinese_softshell_turtle/; 
		$seq =~ s/apaSpi1/Spiny_softshell_turtle/;
		$seq =~ s/anoCar2/green_anole/;
		$seq =~ s/xenTro7/Frog_X_tropicalis/; 
		$seq =~ s/latCha1/Coelacanth/;
		$seq =~ s/tetNig2/Tetraodon/; 
		$seq =~ s/fr3/Fugu/;
		$seq =~ s/takFla1/Yellowbelly_pufferfish/; 
		$seq =~ s/oreNil2/Nile_tilapia/;
		$seq =~ s/neoBri1/Princess_of_Burundi/;
		$seq =~ s/hapBur1/Burtons_mouthbreeder/; 
		$seq =~ s/mayZeb1/Zebra_mbuna/; 
		$seq =~ s/punNye1/Pundamilia_nyererei/; 
		$seq =~ s/oryLat2/Medaka/; 
		$seq =~ s/xipMac1/Southern_platyfish/;
		$seq =~ s/gasAcu1/Stickleback/; 
		$seq =~ s/gadMor1/Atlantic_cod/; 
		$seq =~ s/danRer7/Zebrafish/;
		$seq =~ s/astMex1/Mexican_tetra/; 
		$seq =~ s/lepOcu1/Spotted_gar/;
		$seq =~ s/petMar2/Lamprey/; 
		$seq =~ s/Aapl1/anolis_apletophallus/; 
		$seq =~ s/Aaur1/grass_anole/; 
		$seq =~ s/Afre1/bridled_anole/; 
		
		print OUT "$seq";
	}
	
	close(INPUT);
	close(OUT);
}
	
