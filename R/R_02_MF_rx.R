###############################################################
#-------------- OverView ----------------
#
#  Anole genome relative expression
#
#	Median relative expression with 95% Confidence intervals for:
#   male autosome
#   female autsomes
#   male and female autosomes
#
#   male putative X
#   female putative X
#   male and female putative X
#
#   male LGb
#   female LGb
#   male and female LGb
#
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

setwd("/Users/olneykimberly/Desktop/AnoleChickenSunstitutions/Output/00_Tables")

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
# 4. Male relative expression for autosomes
#--------------------------

#create subset for each chromosome
maleschr6<- subset(mfx, scaffold=='6', select=c(value_1))
maleschr5<- subset(mfx, scaffold=='5', select=c(value_1))
maleschr4<- subset(mfx, scaffold=='4', select=c(value_1))
maleschr3<- subset(mfx, scaffold=='3', select=c(value_1))
maleschr2<- subset(mfx, scaffold=='2', select=c(value_1))
maleschr1<- subset(mfx, scaffold=='1', select=c(value_1))

#create subset for the autosomes for male values combined
maleauto<-c(maleschr1,maleschr2,maleschr3,maleschr4,maleschr5,maleschr6)

mrx1 <-maleschr1$value_1
mrx1<-na.omit(maleschr1$value_1)
median(mrx1)
md1 <- function(x,y) {return(median(sample(mrx1, 2447, replace=TRUE)))}
median1 <- boot(data=maleschr1$value_1,statistic=md1, R=10000)
ci1 <- boot.ci(median1, conf=0.95, type="norm")
ci1$norm

mrx2 <-maleschr2$value_1
mrx2<-na.omit(maleschr2$value_1)
median(mrx2)
md2 <- function(x,y) {return(median(sample(mrx2, 2426, replace=TRUE)))}
median2 <- boot(data=maleschr2$value_1,statistic=md2, R=10000)
ci2 <- boot.ci(median2, conf=0.95, type="norm")
ci2$norm

mrx3 <-maleschr3$value_1
mrx3<-na.omit(maleschr3$value_1)
median(mrx3)
md3 <- function(x,y) {return(median(sample(mrx3, 1661, replace=TRUE)))}
median3 <- boot(data=maleschr3$value_1,statistic=md3, R=10000)
ci3 <- boot.ci(median3, conf=0.95, type="norm")
ci3$norm

mrx4 <-maleschr4$value_1
mrx4<-na.omit(maleschr4$value_1)
median(mrx4)
md4 <- function(x,y) {return(median(sample(mrx4, 1659, replace=TRUE)))}
median4 <- boot(data=maleschr4$value_1,statistic=md4, R=10000)
ci4 <- boot.ci(median4, conf=0.95, type="norm")
ci4$norm

mrx5 <-maleschr5$value_1
mrx5<-na.omit(maleschr5$value_1)
median(mrx5)
md5 <- function(x,y) {return(median(sample(mrx5, 1489, replace=TRUE)))}
median5 <- boot(data=maleschr5$value_1,statistic=md5, R=10000)
ci5 <- boot.ci(median5, conf=0.95, type="norm")
ci5$norm

mrx6 <-maleschr6$value_1
mrx6<-na.omit(maleschr6$value_1)
median(mrx6)
md6 <- function(x,y) {return(median(sample(mrx6, 1083, replace=TRUE)))}
median6 <- boot(data=maleschr6$value_1,statistic=md6, R=10000)
ci6 <- boot.ci(median6, conf=0.95, type="norm")
ci6$norm

mrxauto<-c(mrx1,mrx2,mrx3,mrx4,mrx5,mrx6)
median(mrxauto)
mdauto <- function(x,y) {return(median(sample(mrxauto, 10765, replace=TRUE)))}
medianauto <- boot(data=mrxauto,statistic=mdauto, R=10000)
ciauto <- boot.ci(medianauto, conf=0.95, type="norm")
ciauto$norm

#--------------------------
# 5. Female relative expression for autosomes
#--------------------------

#create subset for each chromosome
femaleschr6<- subset(mfx, scaffold=='6', select=c(value_2))
femaleschr5<- subset(mfx, scaffold=='5', select=c(value_2))
femaleschr4<- subset(mfx, scaffold=='4', select=c(value_2))
femaleschr3<- subset(mfx, scaffold=='3', select=c(value_2))
femaleschr2<- subset(mfx, scaffold=='2', select=c(value_2))
femaleschr1<- subset(mfx, scaffold=='1', select=c(value_2))

#combin the chromomes to get an automosome data set
femaleauto<-c(femaleschr1,femaleschr2,femaleschr3,femaleschr4,femaleschr5,femaleschr6)

fmrx1 <-femaleschr1$value_2
fmrx1<-na.omit(femaleschr1$value_2)
median(fmrx1)
fmd1 <- function(x,y) {return(median(sample(fmrx1, 2447, replace=TRUE)))}
median1 <- boot(data=femaleschr1$value_2,statistic=fmd1, R=10000)
ci1 <- boot.ci(median1, conf=0.95, type="norm")
ci1$norm

fmrx2 <-femaleschr2$value_2
fmrx2<-na.omit(femaleschr2$value_2)
median(fmrx2)
fmd2 <- function(x,y) {return(median(sample(fmrx2, 2426, replace=TRUE)))}
median2 <- boot(data=femaleschr2$value_2,statistic=fmd2, R=10000)
ci2 <- boot.ci(median2, conf=0.95, type="norm")
ci2$norm

fmrx3 <-femaleschr3$value_2
fmrx3<-na.omit(femaleschr3$value_2)
median(fmrx3)
fmd3 <- function(x,y) {return(median(sample(fmrx3, 1661, replace=TRUE)))}
median3 <- boot(data=femaleschr3$value_2,statistic=fmd3, R=10000)
ci3 <- boot.ci(median3, conf=0.95, type="norm")
ci3$norm

fmrx4 <-femaleschr4$value_2
fmrx4<-na.omit(femaleschr4$value_2)
median(fmrx4)
fmd4 <- function(x,y) {return(median(sample(fmrx4, 1659, replace=TRUE)))}
median4 <- boot(data=femaleschr4$value_2,statistic=fmd4, R=10000)
ci4 <- boot.ci(median4, conf=0.95, type="norm")
ci4$norm

fmrx5 <-femaleschr5$value_2
fmrx5<-na.omit(femaleschr5$value_2)
median(fmrx5)
fmd5 <- function(x,y) {return(median(sample(fmrx5, 1489, replace=TRUE)))}
median5 <- boot(data=femaleschr5$value_2,statistic=fmd5, R=10000)
ci5 <- boot.ci(median5, conf=0.95, type="norm")
ci5$norm

fmrx6 <-femaleschr6$value_2
fmrx6<-na.omit(femaleschr6$value_2)
median(fmrx6)
fmd6 <- function(x,y) {return(median(sample(fmrx6, 1083, replace=TRUE)))}
median6 <- boot(data=femaleschr6$value_2,statistic=fmd6, R=10000)
ci6 <- boot.ci(median6, conf=0.95, type="norm")
ci6$norm

fmrxauto<-c(fmrx1,fmrx2,fmrx3,fmrx4,fmrx5,fmrx6)
median(fmrxauto)
fmdauto <- function(x,y) {return(median(sample(fmrxauto, 10765, replace=TRUE)))}
medianfmrxauto <- boot(data=fmrxauto,statistic=fmdauto, R=10000)
ciauto <- boot.ci(medianfmrxauto, conf=0.95, type="norm")
ciauto$norm

#--------------------------
# 6. Meidan relative expression of for male and female autosomes 1-6
#-------------------------

MFrxauto<-c(mfx$rel_exp)
MFrxauto<-na.omit(mfx$rel_exp)
median(MFrxauto)
MFrxautomd <-function(x,y) {return(median(sample(MFrxauto, 21838, replace=TRUE)))}
MFrxautomedian <- boot(data=MFrxauto,statistic=MFrxautomd, R=10000)
MFrxautoci <- boot.ci(MFrxautomedian, conf=0.95, type="norm")
MFrxautoci$norm

#--------------------------
# 7. Meidan relative expression of LGb
#-------------------------

#create subset
LGb<- subset(completechrx, Source=='Ensembl', select=c(relative_expression))

#median relative expression with 95% CI
rxLGb <-LGb$relative_expression
rxLGb<-na.omit(LGb$relative_expression)
median(rxLGb)
mdLGb <- function(x,y) {return(median(sample(rxLGb, 56, replace=TRUE)))}
medianLGb <- boot(data=LGb$relative_expression,statistic=mdLGb, R=10000)
ciLGb <- boot.ci(medianLGb, conf=0.95, type="norm")
ciLGb$norm

# Meidan relative expression of LGb for males
#create subset
LGbmales<- subset(completechrx, Source=='Ensembl', select=c(Male_Expression))

#median relative expression with 95% CI
rxLGbmales <-LGbmales$Male_Expression
rxLGbmales<-na.omit(LGbmales$Male_Expression)
median(rxLGbmales)
mdLGbmales <- function(x,y) {return(median(sample(rxLGbmales, 166, replace=TRUE)))}
medianLGbmales <- boot(data=LGbmales$Male_Expression,statistic=mdLGbmales, R=10000)
ciLGbmales <- boot.ci(medianLGbmales, conf=0.95, type="norm")
ciLGbmales$norm

# Meidan relative expression of LGb for females
#create subset
LGbfemales<- subset(completechrx, Source=='Ensembl', select=c(Female_Expression))

#median relative expression with 95% CI
rxLGbfemales <-LGbfemales$Female_Expression
rxLGbfemales<-na.omit(LGbfemales$Female_Expression)
median(rxLGbfemales)
mdLGbfemales <- function(x,y) {return(median(sample(rxLGbfemales, 166, replace=TRUE)))}
medianLGbfemales <- boot(data=LGbfemales$Female_Expression,statistic=mdLGbfemales, R=10000)
ciLGbfemales <- boot.ci(medianLGbfemales, conf=0.95, type="norm")
ciLGbfemales$norm

#--------------------------
# 8. Meidan relative expression of Hypothetical X
#-------------------------

putativechrx<-PutativeX$rel.exp
putativechrx<-na.omit(PutativeX$rel.exp)
median(putativechrx)
mdgen<-function(x,y) {return(median(sample(putativechrx, 166, replace=TRUE)))}
mediangen <- boot(data=putativechrx,statistic=mdgen, R=10000)
cigen <- boot.ci(mediangen, conf=0.95, type="norm")
cigen$norm

# Meidan relative expression of Hypothetical X males
putativechrxmales<-PutativeX$value_1
putativechrxmales<-na.omit(PutativeX$value_1)
median(putativechrxmales)
mdgenm <-function(x,y) {return(median(sample(putativechrxmales, 166, replace=TRUE)))}
mediangenm <- boot(data=putativechrxmales,statistic=mdgenm, R=10000)
cigenmales <- boot.ci(mediangenm, conf=0.95, type="norm")
cigenmales$norm

# Meidan relative expression of Hypothetical X females
putativechrxfemales<-PutativeX$value_2
putativechrxfemales<-na.omit(PutativeX$value_2)
median(putativechrxfemales)
mdgenf <-function(x,y) {return(median(sample(putativechrxfemales, 166, replace=TRUE)))}
mediangenf <- boot(data=putativechrxfemales,statistic=mdgenf, R=10000)
cigenfemales <- boot.ci(mediangenf, conf=0.95, type="norm")
cigenfemales$norm


#--------------------------
# 9. Meidan relative expression LGbrx/rx_male_autosomes and LGbrx/rx_female_autosomes and for hypothetical X HypXrx/rx_male_autosomes and HypXrx/rx_female_autosomes
#--------------------------

#LGb
median(rxLGbmales)/median(mrxauto)
median(rxLGbfemales)/median(fmrxauto)
median(rxLGb)/median(rxauto)

#Hypothetical X
median(putativechrxmales)/median(mrxauto)
median(putativechrxfemales)/median(fmrxauto)
median(putativechrx)/median(rxauto)
