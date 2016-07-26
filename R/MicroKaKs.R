################################################################################
# Overview
#		This is a Read Me for calculating microchromosomal Ka/Ks values from 
#			KaKs_Calculator output with an FPKM threshold of 2
#
#		Required programs:	R
#					boot package
#					
#################################################################################

#-----------------
# 1. Read in csv file containing Ka/Ks ouput and subset 
#-----------------

# Sequences with premature stops removed
kaks<-read.csv("/media/yitzhak/Data1/Alignments/AnoleChickenSunstitutions/Pairwise_ap/KaKs.scaffolds.noStops.csv",header=TRUE)

# Subset anole genes by locus:
auto<-subset(kaks,Chromosome=="1"|Chromosome=="2"|Chromosome=="3"|Chromosome=="4"|Chromosome=="5"|Chromosome=="6",select=c(Ka,Ks,KaKs))
	auto<-na.omit(auto)
chrx<-subset(kaks,Chromosome=="GL343282.1"|Chromosome=="GL343338.1"|Chromosome=="GL343417.1"|Chromosome=="GL343423.1"|Chromosome=="GL343550.1"|Chromosome=="GL343947.1"|Chromosome=="GL343913.1"|Chromosome=="GL343364.1"|Chromosome=="LGb",select=c(Ka,Ks,KaKs))
	chrx<-na.omit(chrx)

# Microchromosomes
lga <- subset(kaks, Chromosome=="LGa", select=c(Ka,Ks,KaKs))
	lga <- na.omit(lga)
lgb <- subset(kaks, Chromosome=="LGb", select=c(Ka,Ks,KaKs))
	lgb <- na.omit(lgb)
lgc <- subset(kaks, Chromosome=="LGc", select=c(Ka,Ks,KaKs))
	lgc <- na.omit(lgc)
lgd <- subset(kaks, Chromosome=="LGd", select=c(Ka,Ks,KaKs))
	lgd <- na.omit(lgd)
lgf <- subset(kaks, Chromosome=="LGf", select=c(Ka,Ks,KaKs))
	lgf <- na.omit(lgf)
lgg <- subset(kaks, Chromosome=="LGg", select=c(Ka,Ks,KaKs))
	lgg <- na.omit(lgg)
lgh <- subset(kaks, Chromosome=="LGh", select=c(Ka,Ks,KaKs))
	lgh <- na.omit(lgh)
micro <- subset(kaks, Chromosome=="LGa"|Chromosome=="LGc"|Chromosome=="LGd"|Chromosome=="LGf"|Chromosome=="LGg"|Chromosome=="LGh", select=c(Ka,Ks,KaKs))
	micro <- na.omit(micro)

#-----------------
# 2. Calculate average Ka, Ks, and Ka/Ks for LGa
#-----------------

library("boot")

	### Mean ###

# Ka

mean(lga$Ka)
mka <- function(x,y) {return(mean(sample(lga$Ka, 300, replace=TRUE)))}
meanka <- boot(data=lga$Ka,statistic=mka, R=1000)
cika <- boot.ci(meanka, conf=0.95, type="norm")
cika$norm

# Ks

mean(lga$Ks)
mks <- function(x,y) {return(mean(sample(lga$Ks, 300, replace=TRUE)))}
meanks <- boot(data=lga$Ks,statistic=mks, R=1000)
ciks <- boot.ci(meanks, conf=0.95, type="norm")
ciks$norm

# Ka/Ks

mean(lga$KaKs)
mkaks <- function(x,y) {return(mean(sample(lga$KaKs, 300, replace=TRUE)))}
meankaks <- boot(data=lga$KaKs,statistic=mkaks, R=1000)
cikaks <- boot.ci(meankaks, conf=0.95, type="norm")
cikaks$norm

	### Median ###

# Ka

median(lga$Ka)
mka <- function(x,y) {return(median(sample(lga$Ka, 300, replace=TRUE)))}
medianka <- boot(data=lga$Ka,statistic=mka, R=1000)
cika <- boot.ci(medianka, conf=0.95, type="norm")
cika$norm

# Ks

median(lga$Ks)
mks <- function(x,y) {return(median(sample(lga$Ks, 300, replace=TRUE)))}
medianks <- boot(data=lga$Ks,statistic=mks, R=1000)
ciks <- boot.ci(medianks, conf=0.95, type="norm")
ciks$norm

# Ka/Ks

median(lga$KaKs)
mkaks <- function(x,y) {return(median(sample(lga$KaKs, 300, replace=TRUE)))}
mediankaks <- boot(data=lga$KaKs,statistic=mkaks, R=1000)
cikaks <- boot.ci(mediankaks, conf=0.95, type="norm")
cikaks$norm																																																																					

#-----------------
# 3. Calculate average Ka, Ks, and Ka/Ks for LGb
#-----------------

	### Mean ###

# Ka

mean(lgb$Ka)
mka <- function(x,y) {return(mean(sample(lgb$Ka, 300, replace=TRUE)))}
meanka <- boot(data=lgb$Ka,statistic=mka, R=1000)
cika <- boot.ci(meanka, conf=0.95, type="norm")
cika$norm

# Ks

mean(lgb$Ks)
mks <- function(x,y) {return(mean(sample(lgb$Ks, 300, replace=TRUE)))}
meanks <- boot(data=lgb$Ks,statistic=mks, R=1000)
ciks <- boot.ci(meanks, conf=0.95, type="norm")
ciks$norm

# Ka/Ks

mean(lgb$KaKs)
mkaks <- function(x,y) {return(mean(sample(lgb$KaKs, 300, replace=TRUE)))}
meankaks <- boot(data=lgb$KaKs,statistic=mkaks, R=1000)
cikaks <- boot.ci(meankaks, conf=0.95, type="norm")
cikaks$norm

	### Median ###

# Ka

median(lgb$Ka)
mka <- function(x,y) {return(median(sample(lgb$Ka, 300, replace=TRUE)))}
medianka <- boot(data=lgb$Ka,statistic=mka, R=1000)
cika <- boot.ci(medianka, conf=0.95, type="norm")
cika$norm

# Ks

median(lgb$Ks)
mks <- function(x,y) {return(median(sample(lgb$Ks, 300, replace=TRUE)))}
medianks <- boot(data=lgb$Ks,statistic=mks, R=1000)
ciks <- boot.ci(medianks, conf=0.95, type="norm")
ciks$norm

# Ka/Ks

median(lgb$KaKs)
mkaks <- function(x,y) {return(median(sample(lgb$KaKs, 300, replace=TRUE)))}
mediankaks <- boot(data=lgb$KaKs,statistic=mkaks, R=1000)
cikaks <- boot.ci(mediankaks, conf=0.95, type="norm")
cikaks$norm		

#-----------------
# 4. Calculate average Ka, Ks, and Ka/Ks for LGc
#-----------------

	### Mean ###

# Ka

mean(lgc$Ka)
mka <- function(x,y) {return(mean(sample(lgc$Ka, 300, replace=TRUE)))}
meanka <- boot(data=lgc$Ka,statistic=mka, R=1000)
cika <- boot.ci(meanka, conf=0.95, type="norm")
cika$norm

# Ks

mean(lgc$Ks)
mks <- function(x,y) {return(mean(sample(lgc$Ks, 300, replace=TRUE)))}
meanks <- boot(data=lgc$Ks,statistic=mks, R=1000)
ciks <- boot.ci(meanks, conf=0.95, type="norm")
ciks$norm

# Ka/Ks

mean(lgc$KaKs)
mkaks <- function(x,y) {return(mean(sample(lgc$KaKs, 300, replace=TRUE)))}
meankaks <- boot(data=lgc$KaKs,statistic=mkaks, R=1000)
cikaks <- boot.ci(meankaks, conf=0.95, type="norm")
cikaks$norm

	### Median ###

# Ka

median(lgc$Ka)
mka <- function(x,y) {return(median(sample(lgc$Ka, 300, replace=TRUE)))}
medianka <- boot(data=lgc$Ka,statistic=mka, R=1000)
cika <- boot.ci(medianka, conf=0.95, type="norm")
cika$norm

# Ks

median(lgc$Ks)
mks <- function(x,y) {return(median(sample(lgc$Ks, 300, replace=TRUE)))}
medianks <- boot(data=lgc$Ks,statistic=mks, R=1000)
ciks <- boot.ci(medianks, conf=0.95, type="norm")
ciks$norm

# Ka/Ks

median(lgc$KaKs)
mkaks <- function(x,y) {return(median(sample(lgc$KaKs, 300, replace=TRUE)))}
mediankaks <- boot(data=lgc$KaKs,statistic=mkaks, R=1000)
cikaks <- boot.ci(mediankaks, conf=0.95, type="norm")
cikaks$norm		

#-----------------
# 5. Calculate average Ka, Ks, and Ka/Ks for LGd
#-----------------
	
	### Mean ###

# Ka

mean(lgd$Ka)
mka <- function(x,y) {return(mean(sample(lgd$Ka, 300, replace=TRUE)))}
meanka <- boot(data=lgd$Ka,statistic=mka, R=1000)
cika <- boot.ci(meanka, conf=0.95, type="norm")
cika$norm

# Ks

mean(lgd$Ks)
mks <- function(x,y) {return(mean(sample(lgd$Ks, 300, replace=TRUE)))}
meanks <- boot(data=lgd$Ks,statistic=mks, R=1000)
ciks <- boot.ci(meanks, conf=0.95, type="norm")
ciks$norm

# Ka/Ks

mean(lgd$KaKs)
mkaks <- function(x,y) {return(mean(sample(lgd$KaKs, 300, replace=TRUE)))}
meankaks <- boot(data=lgd$KaKs,statistic=mkaks, R=1000)
cikaks <- boot.ci(meankaks, conf=0.95, type="norm")
cikaks$norm

	### Median ###

# Ka

median(lgd$Ka)
mka <- function(x,y) {return(median(sample(lgd$Ka, 300, replace=TRUE)))}
medianka <- boot(data=lgd$Ka,statistic=mka, R=1000)
cika <- boot.ci(medianka, conf=0.95, type="norm")
cika$norm

# Ks

median(lgd$Ks)
mks <- function(x,y) {return(median(sample(lgd$Ks, 300, replace=TRUE)))}
medianks <- boot(data=lgd$Ks,statistic=mks, R=1000)
ciks <- boot.ci(medianks, conf=0.95, type="norm")
ciks$norm

# Ka/Ks

median(lgd$KaKs)
mkaks <- function(x,y) {return(median(sample(lgd$KaKs, 300, replace=TRUE)))}
mediankaks <- boot(data=lgd$KaKs,statistic=mkaks, R=1000)
cikaks <- boot.ci(mediankaks, conf=0.95, type="norm")
cikaks$norm	
	
#-----------------
# 6. Calculate average Ka, Ks, and Ka/Ks for LGf
#-----------------
	
	### Mean ###

# Ka

mean(lgf$Ka)
mka <- function(x,y) {return(mean(sample(lgf$Ka, 300, replace=TRUE)))}
meanka <- boot(data=lgf$Ka,statistic=mka, R=1000)
cika <- boot.ci(meanka, conf=0.95, type="norm")
cika$norm

# Ks

mean(lgf$Ks)
mks <- function(x,y) {return(mean(sample(lgf$Ks, 300, replace=TRUE)))}
meanks <- boot(data=lgf$Ks,statistic=mks, R=1000)
ciks <- boot.ci(meanks, conf=0.95, type="norm")
ciks$norm

# Ka/Ks

mean(lgf$KaKs)
mkaks <- function(x,y) {return(mean(sample(lgf$KaKs, 300, replace=TRUE)))}
meankaks <- boot(data=lgf$KaKs,statistic=mkaks, R=1000)
cikaks <- boot.ci(meankaks, conf=0.95, type="norm")
cikaks$norm

	### Median ###

# Ka

median(lgf$Ka)
mka <- function(x,y) {return(median(sample(lgf$Ka, 300, replace=TRUE)))}
medianka <- boot(data=lgf$Ka,statistic=mka, R=1000)
cika <- boot.ci(medianka, conf=0.95, type="norm")
cika$norm

# Ks

median(lgf$Ks)
mks <- function(x,y) {return(median(sample(lgf$Ks, 300, replace=TRUE)))}
medianks <- boot(data=lgf$Ks,statistic=mks, R=1000)
ciks <- boot.ci(medianks, conf=0.95, type="norm")
ciks$norm

# Ka/Ks

median(lgf$KaKs)
mkaks <- function(x,y) {return(median(sample(lgf$KaKs, 300, replace=TRUE)))}
mediankaks <- boot(data=lgf$KaKs,statistic=mkaks, R=1000)
cikaks <- boot.ci(mediankaks, conf=0.95, type="norm")
cikaks$norm	
	
#-----------------
# 7. Calculate average Ka, Ks, and Ka/Ks for LGg
#-----------------
	
	### Mean ###

# Ka

mean(lgg$Ka)
mka <- function(x,y) {return(mean(sample(lgg$Ka, 300, replace=TRUE)))}
meanka <- boot(data=lgg$Ka,statistic=mka, R=1000)
cika <- boot.ci(meanka, conf=0.95, type="norm")
cika$norm

# Ks

mean(lgg$Ks)
mks <- function(x,y) {return(mean(sample(lgg$Ks, 300, replace=TRUE)))}
meanks <- boot(data=lgg$Ks,statistic=mks, R=1000)
ciks <- boot.ci(meanks, conf=0.95, type="norm")
ciks$norm

# Ka/Ks

mean(lgg$KaKs)
mkaks <- function(x,y) {return(mean(sample(lgg$KaKs, 300, replace=TRUE)))}
meankaks <- boot(data=lgg$KaKs,statistic=mkaks, R=1000)
cikaks <- boot.ci(meankaks, conf=0.95, type="norm")
cikaks$norm

	### Median ###

# Ka

median(lgg$Ka)
mka <- function(x,y) {return(median(sample(lgg$Ka, 300, replace=TRUE)))}
medianka <- boot(data=lgg$Ka,statistic=mka, R=1000)
cika <- boot.ci(medianka, conf=0.95, type="norm")
cika$norm

# Ks

median(lgg$Ks)
mks <- function(x,y) {return(median(sample(lgg$Ks, 300, replace=TRUE)))}
medianks <- boot(data=lgg$Ks,statistic=mks, R=1000)
ciks <- boot.ci(medianks, conf=0.95, type="norm")
ciks$norm

# Ka/Ks

median(lgg$KaKs)
mkaks <- function(x,y) {return(median(sample(lgg$KaKs, 300, replace=TRUE)))}
mediankaks <- boot(data=lgg$KaKs,statistic=mkaks, R=1000)
cikaks <- boot.ci(mediankaks, conf=0.95, type="norm")
cikaks$norm	
	
#-----------------
# 8. Calculate average Ka, Ks, and Ka/Ks for LGh
#-----------------
	
	### Mean ###

# Ka

mean(lgh$Ka)
mka <- function(x,y) {return(mean(sample(lgh$Ka, 300, replace=TRUE)))}
meanka <- boot(data=lgh$Ka,statistic=mka, R=1000)
cika <- boot.ci(meanka, conf=0.95, type="norm")
cika$norm

# Ks

mean(lgh$Ks)
mks <- function(x,y) {return(mean(sample(lgh$Ks, 300, replace=TRUE)))}
meanks <- boot(data=lgh$Ks,statistic=mks, R=1000)
ciks <- boot.ci(meanks, conf=0.95, type="norm")
ciks$norm

# Ka/Ks

mean(lgh$KaKs)
mkaks <- function(x,y) {return(mean(sample(lgh$KaKs, 300, replace=TRUE)))}
meankaks <- boot(data=lgh$KaKs,statistic=mkaks, R=1000)
cikaks <- boot.ci(meankaks, conf=0.95, type="norm")
cikaks$norm

	### Median ###

# Ka

median(lgh$Ka)
mka <- function(x,y) {return(median(sample(lgh$Ka, 300, replace=TRUE)))}
medianka <- boot(data=lgh$Ka,statistic=mka, R=1000)
cika <- boot.ci(medianka, conf=0.95, type="norm")
cika$norm

# Ks

median(lgh$Ks)
mks <- function(x,y) {return(median(sample(lgh$Ks, 300, replace=TRUE)))}
medianks <- boot(data=lgh$Ks,statistic=mks, R=1000)
ciks <- boot.ci(medianks, conf=0.95, type="norm")
ciks$norm

# Ka/Ks

median(lgh$KaKs)
mkaks <- function(x,y) {return(median(sample(lgh$KaKs, 300, replace=TRUE)))}
mediankaks <- boot(data=lgh$KaKs,statistic=mkaks, R=1000)
cikaks <- boot.ci(mediankaks, conf=0.95, type="norm")
cikaks$norm	
	
#-----------------
# 9. Calculate average Ka, Ks, and Ka/Ks for all microchromosomes
#-----------------
	
	### Mean ###

# Ka

mean(micro$Ka)
mka <- function(x,y) {return(mean(sample(micro$Ka, 300, replace=TRUE)))}
meanka <- boot(data=micro$Ka,statistic=mka, R=1000)
cika <- boot.ci(meanka, conf=0.95, type="norm")
cika$norm

# Ks

mean(micro$Ks)
mks <- function(x,y) {return(mean(sample(micro$Ks, 300, replace=TRUE)))}
meanks <- boot(data=micro$Ks,statistic=mks, R=1000)
ciks <- boot.ci(meanks, conf=0.95, type="norm")
ciks$norm

# Ka/Ks

mean(micro$KaKs)
mkaks <- function(x,y) {return(mean(sample(micro$KaKs, 300, replace=TRUE)))}
meankaks <- boot(data=micro$KaKs,statistic=mkaks, R=1000)
cikaks <- boot.ci(meankaks, conf=0.95, type="norm")
cikaks$norm

	### Median ###

# Ka

median(micro$Ka)
mka <- function(x,y) {return(median(sample(micro$Ka, 300, replace=TRUE)))}
medianka <- boot(data=micro$Ka,statistic=mka, R=1000)
cika <- boot.ci(medianka, conf=0.95, type="norm")
cika$norm

# Ks

median(micro$Ks)
mks <- function(x,y) {return(median(sample(micro$Ks, 300, replace=TRUE)))}
medianks <- boot(data=micro$Ks,statistic=mks, R=1000)
ciks <- boot.ci(medianks, conf=0.95, type="norm")
ciks$norm

# Ka/Ks

median(micro$KaKs)
mkaks <- function(x,y) {return(median(sample(micro$KaKs, 300, replace=TRUE)))}
mediankaks <- boot(data=micro$KaKs,statistic=mkaks, R=1000)
cikaks <- boot.ci(mediankaks, conf=0.95, type="norm")
cikaks$norm	
	
#-----------------
# 5. Calculate P values for X-linked genes vs. Microchromosomal genes
# (I copied and pasted autosomal and X-linked genes from the KaKs csv file to use as input and manually removed lines with NAs)
#-----------------

python /home/yitzhak/Dropbox/Toxicofera/ReadMe/Cython/PermutationScripts/KaKs_permutation.py /media/yitzhak/Data1/Alignments/AnoleChickenSunstitutions/Pairwise_ap/ChrXKaKs.noStops.csv /media/yitzhak/Data1/Alignments/AnoleChickenSunstitutions/Pairwise_ap/MicroKaKs.noStops.csv

