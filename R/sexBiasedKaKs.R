################################################################################
# Overview
#		This is a Read Me for calculating sex-biased Ka/Ks values from 
#			the anole X chromosome with no FPKM threshold
#
#		Required programs:	R
#					boot package
#					
#################################################################################

#-----------------
# 1. Read in csv file containing Ka/Ks ouput and subset 
#-----------------

kaks<-read.csv("/home/yitzhak/Dropbox/Squamates/Rupp_lab_notes/ReadMe/DifExp_sra/biasedXGenes.csv",header=TRUE)

# Subset genes by bias (Record number of total genes prior to omitting NAs):
female<-subset(kaks, SexBias=="female",select=KaKs)
	female<-na.omit(female)
male<-subset(kaks, SexBias=="male",select=KaKs)
	male<-na.omit(male)

#-----------------
# 2. Calculate average Ka/Ks for male- and female-biased genes
#-----------------

	library("boot")

# Female-Biased Mean

mean(female$KaKs)
mkaks <- function(x,y) {return(mean(sample(female$KaKs, 300, replace=TRUE)))}
meankaks <- boot(data=female$KaKs,statistic=mkaks, R=1000)
cikaks <- boot.ci(meankaks, conf=0.95, type="norm")
cikaks$norm

# Female-Biased Median

median(female$KaKs)
mkaks <- function(x,y) {return(median(sample(female$KaKs, 300, replace=TRUE)))}
mediankaks <- boot(data=female$KaKs,statistic=mkaks, R=1000)
cikaks <- boot.ci(mediankaks, conf=0.95, type="norm")
cikaks$norm

# Male-Biased Mean

mean(male$KaKs)
mkaks <- function(x,y) {return(mean(sample(male$KaKs, 300, replace=TRUE)))}
meankaks <- boot(data=male$KaKs,statistic=mkaks, R=1000)
cikaks <- boot.ci(meankaks, conf=0.95, type="norm")
cikaks$norm

# Male-Biased Median

median(male$KaKs)
mkaks <- function(x,y) {return(median(sample(male$KaKs, 300, replace=TRUE)))}
mediankaks <- boot(data=male$KaKs,statistic=mkaks, R=1000)
cikaks <- boot.ci(mediankaks, conf=0.95, type="norm")
cikaks$norm																																																																			
	
#-----------------
# 5. Calculate P values for X-linked genes vs. Microchromosomal genes
# (I copied and pasted autosomal and X-linked genes from the KaKs csv file to use as input and manually removed lines with NAs)
#-----------------

python /home/yitzhak/Dropbox/Toxicofera/ReadMe/Cython/PermutationScripts/KaKs_permutation.py /home/yitzhak/Dropbox/Squamates/Rupp_lab_notes/ReadMe/DifExp_sra/Data/femaleBiasedXGenes.csv /home/yitzhak/Dropbox/Squamates/Rupp_lab_notes/ReadMe/DifExp_sra/Data/MaleBiasedXGenes.csv