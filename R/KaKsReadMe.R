################################################################################
# Overview
#		This is a Read Me for calculating average Ka/Ks values from 
#			KaKs_Calculator output with an FPKM threshold of 2
#
#		Required programs:	R
#					boot package
#					removeLines.py
#					KaKsConcatenate.py
#					
#################################################################################

#-----------------
# 0. Remove genes with premature stop codons from Ka/Ks output csv and join the output with list of gene IKs, chromosomes, and other locus data
#-----------------

# Remove genes:
python /home/yitzhak/Dropbox/Toxicofera/ReadMe/PythonScripts/removeLines.py /media/yitzhak/Data1/Alignments/AnoleChickenSunstitutions/Pairwise_ap/internalStops.txt /media/yitzhak/Data1/Alignments/AnoleChickenSunstitutions/Pairwise_ap/KaKs.txt

# Join files:
  # Sequences with premature stops removed
python /home/yitzhak/Dropbox/AvesAlignments/KaKsConcatenate.py /media/yitzhak/Data1/Alignments/AnoleChickenSunstitutions/anoCarGenes.csv /media/yitzhak/Data1/Alignments/AnoleChickenSunstitutions/Pairwise_ap/KaKs.rmLines.txt

  # Sequences with premature stops retained (change the output file name in the script)
python /home/yitzhak/Dropbox/AvesAlignments/KaKsConcatenate.py /media/yitzhak/Data1/Alignments/AnoleChickenSunstitutions/anoCarGenes.csv /media/yitzhak/Data1/Alignments/AnoleChickenSunstitutions/Pairwise_ap/KaKs.txt

#-----------------
# 1. Read in csv file containing Ka/Ks ouput and subset 
#-----------------

# Sequences with premature stops removed
kaks<-read.csv("/media/yitzhak/Data1/Alignments/AnoleChickenSunstitutions/Pairwise_ap/KaKs.scaffolds.noStops.csv",header=TRUE)

# Sequences with premature stops retained
kaks<-read.csv("/media/yitzhak/Data1/Alignments/AnoleChickenSunstitutions/Pairwise_ap/KaKs.scaffolds.withStops.csv",header=TRUE)

# Subset anole genes by locus:
auto<-subset(kaks,Chromosome=="1"|Chromosome=="2"|Chromosome=="3"|Chromosome=="4"|Chromosome=="5"|Chromosome=="6",select=c(Ka,Ks,KaKs))
auto<-na.omit(auto)
chrx<-subset(kaks,Chromosome=="GL343282.1"|Chromosome=="GL343338.1"|Chromosome=="GL343417.1"|Chromosome=="GL343423.1"|Chromosome=="GL343550.1"|Chromosome=="GL343947.1"|Chromosome=="GL343913.1"|Chromosome=="GL343364.1"|Chromosome=="LGb",select=c(Ka,Ks,KaKs))
chrx<-na.omit(chrx)

#-----------------
# 2. Calculate mean Ka, Ks, and Ka/Ks 
#-----------------

library("boot")

# Ka

mean(auto$Ka)
mka <- function(x,y) {return(mean(sample(auto$Ka, 300, replace=TRUE)))}
meanka <- boot(data=auto$Ka,statistic=mka, R=1000)
cika <- boot.ci(meanka, conf=0.95, type="norm")
cika$norm

mean(chrx$Ka)
macxKa <- function(x,y) {return(mean(sample(chrx$Ka, 300, replace=TRUE)))}
meanacxKa <- boot(data=chrx$Ka,statistic=macxKa, R=1000)
ciacxKa <- boot.ci(meanacxKa, conf=0.95, type="norm")
ciacxKa$norm

# Ks

mean(auto$Ks)
mks <- function(x,y) {return(mean(sample(auto$Ks, 300, replace=TRUE)))}
meanks <- boot(data=auto$Ks,statistic=mks, R=1000)
ciks <- boot.ci(meanks, conf=0.95, type="norm")
ciks$norm

mean(chrx$Ks)
macxKs <- function(x,y) {return(mean(sample(chrx$Ks, 300, replace=TRUE)))}
meanacxKs <- boot(data=chrx$Ks,statistic=macxKs, R=1000)
ciacxKs <- boot.ci(meanacxKs, conf=0.95, type="norm")
ciacxKs$norm

# Ka/Ks

mean(auto$KaKs)
mkaks <- function(x,y) {return(mean(sample(auto$KaKs, 300, replace=TRUE)))}
meankaks <- boot(data=auto$KaKs,statistic=mkaks, R=1000)
cikaks <- boot.ci(meankaks, conf=0.95, type="norm")
cikaks$norm

mean(chrx$KaKs)
macxKaKs <- function(x,y) {return(mean(sample(chrx$KaKs, 300, replace=TRUE)))}
meanacxKaKs <- boot(data=chrx$KaKs,statistic=macxKaKs, R=1000)
ciacxKaKs <- boot.ci(meanacxKaKs, conf=0.95, type="norm")
ciacxKaKs$norm


#-----------------
# 3. Calculate median Ka, Ks, and Ka/Ks for anole when compaired to human and chicken
#-----------------

# Ka

median(auto$Ka)
mka <- function(x,y) {return(median(sample(auto$Ka, 300, replace=TRUE)))}
medianka <- boot(data=auto$Ka,statistic=mka, R=1000)
cika <- boot.ci(medianka, conf=0.95, type="norm")
cika$norm

median(chrx$Ka)
macxKa <- function(x,y) {return(median(sample(chrx$Ka, 300, replace=TRUE)))}
medianacxKa <- boot(data=chrx$Ka,statistic=macxKa, R=1000)
ciacxKa <- boot.ci(medianacxKa, conf=0.95, type="norm")
ciacxKa$norm

# Ks

median(auto$Ks)
mks <- function(x,y) {return(median(sample(auto$Ks, 300, replace=TRUE)))}
medianks <- boot(data=auto$Ks,statistic=mks, R=1000)
ciks <- boot.ci(medianks, conf=0.95, type="norm")
ciks$norm

median(chrx$Ks)
macxKs <- function(x,y) {return(median(sample(chrx$Ks, 300, replace=TRUE)))}
medianacxKs <- boot(data=chrx$Ks,statistic=macxKs, R=1000)
ciacxKs <- boot.ci(medianacxKs, conf=0.95, type="norm")
ciacxKs$norm

# Ka/Ks

median(auto$KaKs)
mkaks <- function(x,y) {return(median(sample(auto$KaKs, 300, replace=TRUE)))}
mediankaks <- boot(data=auto$KaKs,statistic=mkaks, R=1000)
cikaks <- boot.ci(mediankaks, conf=0.95, type="norm")
cikaks$norm

median(chrx$KaKs)
macxKaKs <- function(x,y) {return(median(sample(chrx$KaKs, 300, replace=TRUE)))}
medianacxKaKs <- boot(data=chrx$KaKs,statistic=macxKaKs, R=1000)
ciacxKaKs <- boot.ci(medianacxKaKs, conf=0.95, type="norm")
ciacxKaKs$norm

#-----------------
# 4. Calculate P values for X-linked genes vs. Autosomal genes
# (I copied and pasted autosomal and X-linked genes from the KaKs csv file to use as input and manually removed lines with NAs)
#-----------------

# Sequences with premature stops removed
python /home/yitzhak/Dropbox/AvesAlignments/Auto_Xsubset.py /media/yitzhak/Data1/Alignments/AnoleChickenSunstitutions/Pairwise_ap/KaKs.scaffolds.noStops.csv 

python /home/yitzhak/Dropbox/Toxicofera/ReadMe/PythonScripts/permutation.py /media/yitzhak/Data1/Alignments/AnoleChickenSunstitutions/Pairwise_ap/ChrXKaKs.noStops.csv /media/yitzhak/Data1/Alignments/AnoleChickenSunstitutions/Pairwise_ap/AutoKaKs.noStops.csv

# Sequences with premature stops retained
python /home/yitzhak/Dropbox/AvesAlignments/Auto_Xsubset.py /media/yitzhak/Data1/Alignments/AnoleChickenSunstitutions/Pairwise_ap/KaKs.scaffolds.withStops.csv

python /home/yitzhak/Dropbox/Toxicofera/ReadMe/PythonScripts/permutation.py /media/yitzhak/Data1/Alignments/AnoleChickenSunstitutions/Pairwise_ap/ChrXKaKs.withStops.csv /media/yitzhak/Data1/Alignments/AnoleChickenSunstitutions/Pairwise_ap/AutoKaKs.withStops.csv
