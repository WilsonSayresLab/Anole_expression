###############################################################
#-------------- OverView ----------------
#
#  Anole genome relative expression
#
#	Median and mean relative expression with 95% Confidence intervals
#		Automosomes 1-6
# 		Chromosome x
#		Hypothetical chromosome x
#		LGb
#
#############################################################

######## Read in genome information README #######

#--------------------------
# 1. Prepare directories and files
#--------------------------

# make directories
#cd 				# move to home directory or location you would like to run the analysis
#mkdir /Anole

#--------------------------
# 2. Open R libraries and set working directory
#--------------------------

# setwd("/Users/olneykimberly/Desktop/AnoleChickenSunstitutions/Output/00_Tables")
# set working directory 

library(gdata)
library(boot)

#--------------------------
# 3. Read in genome information for male and female gene expression for chr1-6 and the autosomes 1-6, the complete X chromsome, and the hypothetical X chromosome
#--------------------------

#Table of Male and Female gene expression
mfx<-read.xls("MF_GeneExp.xlsx",header=TRUE)

#Table of the complete X chromosome gene expression
completechrx<-read.xls("Complete_Chr_X.xlsx",header=TRUE)

#Table of the hypothetical chromsome X
PutativeX<-read.xls("PutativeX.xlsx",header=TRUE)

#--------------------------
# 4. Create subset of the relative expression for:
#		 autosomes 1-6
#		 Putative X
#		 LGb
#		 Genome
#--------------------------

# autosomes
chr1<- subset(mfx, scaffold=='1', select=c(rel_exp))
chr2<- subset(mfx, scaffold=='2', select=c(rel_exp))
chr3<- subset(mfx, scaffold=='3', select=c(rel_exp))
chr4<- subset(mfx, scaffold=='4', select=c(rel_exp))
chr5<- subset(mfx, scaffold=='5', select=c(rel_exp))
chr6<- subset(mfx, scaffold=='6', select=c(rel_exp))
auto<-c(chr1,chr2,chr3,chr4,chr5,chr6)

# X chromosome
chrx<- subset(completechrx, select=c(relative_expression))
PutativeChrX<- subset(PutativeX, select=c(rel.exp))
LGb<- subset(completechrx, Source=='Ensembl', select=c(relative_expression))

# Genome
Genome<- subset(mfx, select=c(rel_exp))
Genome<-na.omit(Genome$rel_exp)

#--------------------------
# 5. Find the medians with a 95% CI
#--------------------------

# automes 1-6
rx1 <-chr1$rel_exp
rx1<-na.omit(chr1$rel_exp)
median(rx1)
md1 <- function(x,y) {return(median(sample(rx1, 2341, replace=TRUE)))}
median1 <- boot(data=chr1$rel_exp,statistic=md1, R=10000)
ci1 <- boot.ci(median1, conf=0.95, type="norm")
ci1$norm

rx2 <- chr2$rel_exp
rx2<-na.omit(chr2$rel_exp)
median(rx2)
md2 <- function(x,y) {return(median(sample(rx2, 2335, replace=TRUE)))}
median2 <- boot(data=chr2$rel_exp,statistic=md2, R=10000)
ci2 <- boot.ci(median2, conf=0.95, type="norm")
ci2$norm

rx3 <- chr3$rel_exp
rx3<-na.omit(chr3$rel_exp)
median(rx3)
md3 <- function(x,y) {return(median(sample(rx3, 1586, replace=TRUE)))}
median3 <- boot(data=chr3$rel_exp,statistic=md3, R=10000)
ci3 <- boot.ci(median3, conf=0.95, type="norm")
ci3$norm

rx4 <- chr4$rel_exp
rx4<-na.omit(chr4$rel_exp)
median(rx4)
md4 <- function(x,y) {return(median(sample(rx4, 1567, replace=TRUE)))}
median4 <- boot(data=chr4$rel_exp,statistic=md4, R=10000)
ci4 <- boot.ci(median4, conf=0.95, type="norm")
ci4$norm

rx5 <- chr5$rel_exp
rx5<-na.omit(chr5$rel_exp)
median(rx5)
md5 <- function(x,y) {return(median(sample(rx5, 1409, replace=TRUE)))}
median5 <- boot(data=chr5$rel_exp,statistic=md5, R=10000)
ci5 <- boot.ci(median5, conf=0.95, type="norm")
ci5$norm

rx6 <- chr6$rel_exp
rx6<-na.omit(chr6$rel_exp)
median(rx6)
md6 <- function(x,y) {return(median(sample(rx6, 1026, replace=TRUE)))}
median6 <- boot(data=chr6$rel_exp,statistic=md6, R=10000)
ci6 <- boot.ci(median6, conf=0.95, type="norm")
ci6$norm

# autosomes together 1-6
rxauto<-c(rx1,rx2,rx3,rx4,rx5,rx6)
median(rxauto)
mdauto <- function(x,y) {return(median(sample(rxauto, 10264, replace=TRUE)))}
medianauto <- boot(data=rxauto,statistic=mdauto, R=10000)
ciauto <- boot.ci(medianauto, conf=0.95, type="norm")
ciauto$norm


# Putative X chromosome
PutativeChrX<-PutativeX$rel.exp
PutativeChrX<-na.omit(PutativeX$rel.exp)
median(PutativeChrX)
PutativeChrXmdgen<-function(x,y) {return(median(sample(PutativeChrX, 166, replace=TRUE)))}
PutativeChrXmediangen <- boot(data=PutativeChrX,statistic=PutativeChrXmdgen, R=10000)
PutativeChrXci <- boot.ci(PutativeChrXmediangen, conf=0.95, type="norm")
PutativeChrXci$norm

# LGb scaffold
rxLGb <-LGb$relative_expression
rxLGb<-na.omit(LGb$relative_expression)
median(rxLGb)
mdLGb <- function(x,y) {return(median(sample(rxLGb, 56, replace=TRUE)))}
medianLGb <- boot(data=LGb$relative_expression,statistic=mdLGb, R=10000)
ci1 <- boot.ci(medianLGb, conf=0.95, type="norm")
ci1$norm

# whole geome
median(Genome)
mdgen<-function(x,y) {return(median(sample(Genome, 21838, replace=TRUE)))}
mediangen <- boot(data=Genome,statistic=mdgen, R=10000)
cigen <- boot.ci(mediangen, conf=0.95, type="norm")
cigen$norm

#--------------------------
# 6. Find the means with a 95% CI
#--------------------------

# automes 1-6
mean(rx1)
mn1 <- function(x,y) {return(mean(sample(rx1, 2341, replace=TRUE)))}
mean1 <- boot(data=chr1$rel_exp,statistic=mn1, R=1000)
ci1 <- boot.ci(mean1, conf=0.95, type="norm")
ci1$norm

mean(rx2)
mn2 <- function(x,y) {return(mean(sample(rx2, 2335, replace=TRUE)))}
mean2 <- boot(data=chr2$rel_exp,statistic=mn2, R=1000)
ci2 <- boot.ci(mean2, conf=0.95, type="norm")
ci2$norm

mean(rx3)
mn3 <- function(x,y) {return(mean(sample(rx3, 1586, replace=TRUE)))}
mean3 <- boot(data=chr3$rel_exp,statistic=mn3, R=1000)
ci3 <- boot.ci(mean3, conf=0.95, type="norm")
ci3$norm

mean(rx4)
mn4 <- function(x,y) {return(mean(sample(rx4, 1567, replace=TRUE)))}
mean4 <- boot(data=chr4$rel_exp,statistic=mn4, R=1000)
ci4 <- boot.ci(mean4, conf=0.95, type="norm")
ci4$norm

mean(rx5)
mn5 <- function(x,y) {return(mean(sample(rx5, 1409, replace=TRUE)))}
mean5 <- boot(data=chr5$rel_exp,statistic=mn5, R=1000)
ci5 <- boot.ci(mean5, conf=0.95, type="norm")
ci5$norm

mean(rx6)
mn6 <- function(x,y) {return(mean(sample(rx6, 1026, replace=TRUE)))}
mean6 <- boot(data=chr6$rel_exp,statistic=mn6, R=1000)
ci6 <- boot.ci(mean6, conf=0.95, type="norm")
ci6$norm

# autosomes together 1-6
mean(rxauto)
mnauto <- function(x,y) {return(mean(sample(rxauto, 10264, replace=TRUE)))}
meanauto <- boot(data=rxauto,statistic=mnauto, R=1000)
ciauto <- boot.ci(meanauto, conf=0.95, type="norm")
ciauto$norm



# Putative X
mean(PutativeChrX)
PutativeChrXmn<-function(x,y) {return(median(sample(PutativeChrX, 166, replace=TRUE)))}
PutativeChrXmean <- boot(data=PutativeChrX,statistic=PutativeChrXmn, R=10000)
PutativeChrXci <- boot.ci(PutativeChrXmean, conf=0.95, type="norm")
PutativeChrXci$norm

# LGb scaffold
mean(rxLGb)
mnLGb <- function(x,y) {return(mean(sample(rxLGb, 56, replace=TRUE)))}
meanLGb <- boot(data=LGb$relative_expression,statistic=mnLGb, R=10000)
ci1 <- boot.ci(meanLGb, conf=0.95, type="norm")
ci1$norm

# whole geome
mean(Genome)
mngen<-function(x,y) {return(mean(sample(Genome, 21838, replace=TRUE)))}
meangen <- boot(data=Genome,statistic=mngen, R=10000)
cigen <- boot.ci(meangen, conf=0.95, type="norm")
cigen$norm
