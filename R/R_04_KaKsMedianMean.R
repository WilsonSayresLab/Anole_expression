###############################################################
#-------------- OverView ----------------
#
#  KaKs Median and Mean
#     nonsynonymous (Ka) and synonymous (Ks) substitution rates
#     
#############################################################

######## Read in KaKs information README #######

#--------------------------
# 1. Open R libraries and set working directory
#--------------------------

setwd("/Users/olneykimberly/Desktop/AnoleChickenSunstitutions/Output/KaKs_Out")

library(gdata)
library(boot)

#--------------------------
# 3. Read in KaKs information for:
#		 autosomes 1-6
#		 LGb
#        putative X
#--------------------------

# Contigs
AAWZ02037698<-read.csv("AAWZ02037698.csv", fill = TRUE, header=TRUE)
AAWZ02040127<-read.csv("AAWZ02040127.csv", fill = TRUE, header=TRUE)
AAWZ02041607<-read.csv("AAWZ02041607.csv", fill = TRUE, header=TRUE)

# autosomes
Chr1<-read.csv("Chr1.csv", fill = TRUE, header=TRUE)
Chr2<-read.csv("Chr2.csv", fill = TRUE, header=TRUE)
Chr3<-read.csv("Chr3.csv", fill = TRUE, header=TRUE)
Chr4<-read.csv("Chr4.csv", fill = TRUE, header=TRUE)
Chr5<-read.csv("Chr5.csv", fill = TRUE, header=TRUE)
Chr6<-read.csv("Chr6.csv", fill = TRUE, header=TRUE)
A<-read.csv("autosomes.csv", fill =TRUE, header=TRUE) # autosomes 1-6

# scaffolds
GL343282<-read.csv("GL343282.csv", fill = TRUE, header=TRUE)
GL343338<-read.csv("GL343338.csv", fill = TRUE, header=TRUE)
GL343364<-read.csv("GL343364.csv", fill = TRUE, header=TRUE)
GL343417<-read.csv("GL343417.csv", fill = TRUE, header=TRUE)
GL343422<-read.csv("GL343422.csv", fill = TRUE, header=TRUE)
GL343423<-read.csv("GL343423.csv", fill = TRUE, header=TRUE)
GL343439<-read.csv("GL343439.csv", fill = TRUE, header=TRUE)
GL343462<-read.csv("GL343462.csv", fill = TRUE, header=TRUE)
GL343516<-read.csv("GL343516.csv", fill = TRUE, header=TRUE)
GL343525<-read.csv("GL343525.csv", fill = TRUE, header=TRUE)
GL343550<-read.csv("GL343550.csv", fill = TRUE, header=TRUE)
GL343588<-read.csv("GL343588.csv", fill = TRUE, header=TRUE)
GL343731<-read.csv("GL343731.csv", fill = TRUE, header=TRUE)
GL343947<-read.csv("GL343947.csv", fill = TRUE, header=TRUE)
GL344042<-read.csv("GL344042.csv", fill = TRUE, header=TRUE)
GL344393<-read.csv("GL344393.csv", fill = TRUE, header=TRUE)
GL344496<-read.csv("GL344496.csv", fill = TRUE, header=TRUE)
GL344539<-read.csv("GL344539.csv", fill = TRUE, header=TRUE)
LGb<-read.csv("LGb.csv", fill = TRUE, header=TRUE)

px<-read.csv("HypX.csv", fill =TRUE, header=TRUE) # putative X


#--------------------------
# 5. Find the medians for Ka values with a 95% CI
#    Find the mean for Ka values with a 95% CI for autosomes, LGb, and putative X
#--------------------------

# contigs
AAWZ02037698Ka <-AAWZ02037698$Ka
AAWZ02037698Ka<-na.omit(AAWZ02037698$Ka)
median(AAWZ02037698Ka) # only one value no need to bootstrap
#AAWZ02037698md <- function(x,y) {return(median(sample(AAWZ02037698Ka, 2, replace=TRUE)))}
#AAWZ02037698median <- boot(data=AAWZ02037698$Ka,statistic=AAWZ02037698md, R=10000)
#AAWZ02037698ci <- boot.ci(AAWZ02037698median, conf=0.95, type="norm")
#AAWZ02037698ci$norm

AAWZ02040127Ka <-AAWZ02040127$Ka
AAWZ02040127Ka<-na.omit(AAWZ02040127$Ka)
median(AAWZ02040127Ka)
AAWZ02040127md <- function(x,y) {return(median(sample(AAWZ02040127Ka, 1, replace=TRUE)))}
AAWZ02040127median <- boot(data=AAWZ02040127$Ka,statistic=AAWZ02040127md, R=10000)
AAWZ02040127ci <- boot.ci(AAWZ02040127median, conf=0.95, type="norm")
AAWZ02040127ci$norm

# Has no value
AAWZ02041607Ka <-AAWZ02041607$Ka
AAWZ02041607Ka<-na.omit(AAWZ02041607$Ka)
median(AAWZ02041607Ka)
#AAWZ02041607md <- function(x,y) {return(median(sample(AAWZ02041607Ka, 1, replace=TRUE)))}
#AAWZ02041607median <- boot(data=AAWZ02041607$Ka,statistic=AAWZ02041607md, R=10000)
#AAWZ02041607ci <- boot.ci(AAWZ02041607median, conf=0.95, type="norm")
#AAWZ02041607ci$norm

Kacontigs<-c(AAWZ02037698Ka,AAWZ02040127Ka) #AAWZ02041607Ka
median(Kacontigs)
mdcontigs <- function(x,y) {return(median(sample(Kacontigs, 2, replace=TRUE)))} #only 2 observed values
mediancontigs <- boot(data=Kacontigs,statistic=mdcontigs, R=10000)
cicontigs <- boot.ci(mediancontigs, conf=0.95, type="norm")
cicontigs$norm

# autosomes 1-6
Chr1Ka<-Chr1$Ka
Chr1Ka<-na.omit(Chr1$Ka)
median(Chr1Ka)
Chr1md <- function(x,y) {return(median(sample(Chr1Ka, 603, replace=TRUE)))}
Chr1median <- boot(data=Chr1$Ka,statistic=Chr1md, R=10000)
Chr1ci <- boot.ci(Chr1median, conf=0.95, type="norm")
Chr1ci$norm

# doesn't have a header
Chr2Ka<-Chr2$Ka
Chr2Ka<-na.omit(Chr2$Ka)
median(Chr2Ka)
Chr2md <- function(x,y) {return(median(sample(Chr2Ka, 679, replace=TRUE)))}
Chr2median <- boot(data=Chr2$Ka,statistic=Chr2md, R=10000)
Chr2ci <- boot.ci(Chr2median, conf=0.95, type="norm")
Chr2ci$norm

Chr3Ka<-Chr3$Ka
Chr3Ka<-na.omit(Chr3$Ka)
median(Chr3Ka)
Chr3md <- function(x,y) {return(median(sample(Chr3Ka, 81, replace=TRUE)))}
Chr3median <- boot(data=Chr3$Ka,statistic=Chr3md, R=10000)
Chr3ci <- boot.ci(Chr3median, conf=0.95, type="norm")
Chr3ci$norm

Chr4Ka<-Chr4$Ka
Chr4Ka<-na.omit(Chr4$Ka)
median(Chr4Ka)
Chr4md <- function(x,y) {return(median(sample(Chr4Ka, 146, replace=TRUE)))}
Chr4median <- boot(data=Chr4$Ka,statistic=Chr4md, R=10000)
Chr4ci <- boot.ci(Chr4median, conf=0.95, type="norm")
Chr4ci$norm

Chr5Ka<-Chr5$Ka
Chr5Ka<-na.omit(Chr5$Ka)
median(Chr5Ka)
Chr5md <- function(x,y) {return(median(sample(Chr5Ka, 99, replace=TRUE)))}
Chr5median <- boot(data=Chr5$Ka,statistic=Chr5md, R=10000)
Chr5ci <- boot.ci(Chr5median, conf=0.95, type="norm")
Chr5ci$norm

Chr6Ka<-Chr6$Ka
Chr6Ka<-na.omit(Chr6$Ka)
median(Chr6Ka)
Chr6md <- function(x,y) {return(median(sample(Chr6Ka, 99, replace=TRUE)))}
Chr6median <- boot(data=Chr6$Ka,statistic=Chr6md, R=10000)
Chr6ci <- boot.ci(Chr6median, conf=0.95, type="norm")
Chr6ci$norm

# all autsomes 1-6
# Median
AKa<-na.omit(A$Ka)
median(AKa)
AKamd<- function(x,y) {return(median(sample(AKa, 1940, replace=TRUE)))}
AKamedian<- boot(data=AKa,statistic=AKamd, R=10000)
AKaci <- boot.ci(AKamedian, conf=0.95, type="norm")
AKaci$norm

# Mean
AKa<-na.omit(A$Ka)
mean(AKa)
AKamd<- function(x,y) {return(mean(sample(AKa, 1940, replace=TRUE)))}
AKamean<- boot(data=AKa,statistic=AKamd, R=10000)
AKaci <- boot.ci(AKamean, conf=0.95, type="norm")
AKaci$norm

# scaffolds
GL343282Ka<-GL343282$Ka
GL343282Ka<-na.omit(GL343282$Ka)
median(GL343282Ka)
GL343282md <- function(x,y) {return(median(sample(GL343282Ka, 9, replace=TRUE)))}
GL343282median <- boot(data=GL343282$Ka,statistic=GL343282md, R=10000)
GL343282ci <- boot.ci(GL343282median, conf=0.95, type="norm")
GL343282ci$norm

GL343338Ka<-GL343338$Ka
GL343338Ka<-na.omit(GL343338$Ka)
median(GL343338Ka)
GL343338md <- function(x,y) {return(median(sample(GL343338Ka, 11, replace=TRUE)))}
GL343338median <- boot(data=GL343338$Ka,statistic=GL343338md, R=10000)
GL343338ci <- boot.ci(GL343338median, conf=0.95, type="norm")
GL343338ci$norm

GL343364Ka<-GL343364$Ka
GL343364Ka<-na.omit(GL343364$Ka)
median(GL343364Ka)
GL343364md <- function(x,y) {return(median(sample(GL343364Ka, 8, replace=TRUE)))}
GL343364median <- boot(data=GL343364$Ka,statistic=GL343364md, R=10000)
GL343364ci <- boot.ci(GL343364median, conf=0.95, type="norm")
GL343364ci$norm

GL343417Ka<-GL343417$Ka
GL343417Ka<-na.omit(GL343417$Ka)
median(GL343417Ka)
GL343417md <- function(x,y) {return(median(sample(GL343417Ka, 3, replace=TRUE)))}
GL343417median <- boot(data=GL343417$Ka,statistic=GL343417md, R=10000)
GL343417ci <- boot.ci(GL343417median, conf=0.95, type="norm")
GL343417ci$norm

GL343422Ka<-GL343422$Ka
GL343422Ka<-na.omit(GL343422$Ka)
median(GL343422Ka)
GL343422md <- function(x,y) {return(median(sample(GL343422Ka, 3, replace=TRUE)))}
GL343422median <- boot(data=GL343422$Ka,statistic=GL343422md, R=10000)
GL343422ci <- boot.ci(GL343422median, conf=0.95, type="norm")
GL343422ci$norm

GL343423Ka<-GL343423$Ka
GL343423Ka<-na.omit(GL343423$Ka)
median(GL343423Ka)
GL343423md <- function(x,y) {return(median(sample(GL343423Ka, 4, replace=TRUE)))}
GL343423median <- boot(data=GL343423$Ka,statistic=GL343423md, R=10000)
GL343423ci <- boot.ci(GL343423median, conf=0.95, type="norm")
GL343423ci$norm

GL343439Ka<-GL343439$Ka
GL343439Ka<-na.omit(GL343439$Ka)
median(GL343439Ka)
GL343439md <- function(x,y) {return(median(sample(GL343439Ka, 1, replace=TRUE)))}
GL343439median <- boot(data=GL343439$Ka,statistic=GL343439md, R=10000)
GL343439ci <- boot.ci(GL343439median, conf=0.95, type="norm")
GL343439ci$norm

GL343462Ka<-GL343462$Ka
GL343462Ka<-na.omit(GL343462$Ka)
median(GL343462Ka)
GL343462md <- function(x,y) {return(median(sample(GL343462Ka, 1, replace=TRUE)))}
GL343462median <- boot(data=GL343462$Ka,statistic=GL343462md, R=10000)
GL343462ci <- boot.ci(GL343462median, conf=0.95, type="norm")
GL343462ci$norm

GL343516Ka<-GL343516$Ka
GL343516Ka<-na.omit(GL343516$Ka)
median(GL343516Ka)
GL343516md <- function(x,y) {return(median(sample(GL343516Ka, 4, replace=TRUE)))}
GL343516median <- boot(data=GL343516$Ka,statistic=GL343516md, R=10000)
GL343516ci <- boot.ci(GL343516median, conf=0.95, type="norm")
GL343516ci$norm

GL343525Ka<-GL343525$Ka
GL343525Ka<-na.omit(GL343525$Ka)
median(GL343525Ka)
GL343525md <- function(x,y) {return(median(sample(GL343525Ka, 1, replace=TRUE)))}
GL343525median <- boot(data=GL343525$Ka,statistic=GL343525md, R=10000)
GL343525ci <- boot.ci(GL343525median, conf=0.95, type="norm")
GL343525ci$norm

GL343550Ka<-GL343550$Ka
GL343550Ka<-na.omit(GL343550$Ka)
median(GL343550Ka)
GL343550md <- function(x,y) {return(median(sample(GL343550Ka, 7, replace=TRUE)))}
GL343550median <- boot(data=GL343550$Ka,statistic=GL343550md, R=10000)
GL343550ci <- boot.ci(GL343550median, conf=0.95, type="norm")
GL343550ci$norm

# Has no value
#GL343588Ka<-GL343588$Ka
#GL343588Ka<-na.omit(GL343588$Ka)
#median(GL343588Ka)
#GL343588md <- function(x,y) {return(median(sample(GL343588Ka, 1, replace=TRUE)))}
#GL343588median <- boot(data=GL343588$Ka,statistic=GL343588md, R=10000)
#GL343588ci <- boot.ci(GL343588median, conf=0.95, type="norm")
#GL343588ci$norm

# Has no value
#GL343731Ka<-GL343731$Ka
#GL343731Ka<-na.omit(GL343731$Ka)
#median(GL343731Ka)
#GL343731md <- function(x,y) {return(median(sample(GL343731Ka, 1, replace=TRUE)))}
#GL343731median <- boot(data=GL343731$Ka,statistic=GL343731md, R=10000)
#GL343731ci <- boot.ci(GL343731median, conf=0.95, type="norm")
#GL343731ci$norm

GL343947Ka<-GL343947$Ka
GL343947Ka<-na.omit(GL343947$Ka)
median(GL343947Ka)
GL343947md <- function(x,y) {return(median(sample(GL343947Ka, 6, replace=TRUE)))}
GL343947median <- boot(data=GL343947$Ka,statistic=GL343947md, R=10000)
GL343947ci <- boot.ci(GL343947median, conf=0.95, type="norm")
GL343947ci$norm

GL344042Ka<-GL344042$Ka
GL344042Ka<-na.omit(GL344042$Ka)
median(GL344042Ka)
GL344042md <- function(x,y) {return(median(sample(GL344042Ka, 1, replace=TRUE)))}
GL344042median <- boot(data=GL344042$Ka,statistic=GL344042md, R=10000)
GL344042ci <- boot.ci(GL344042median, conf=0.95, type="norm")
GL344042ci$norm

# Has no value
#GL344393Ka<-GL344393$Ka
#GL344393Ka<-na.omit(GL344393$Ka)
#median(GL344393Ka)
#GL344393md <- function(x,y) {return(median(sample(GL344393Ka, 1, replace=TRUE)))}
#GL344393median <- boot(data=GL344393$Ka,statistic=GL344393md, R=10000)
#GL344393ci <- boot.ci(GL344393median, conf=0.95, type="norm")
#GL344393ci$norm

# Has no value
#GL344496Ka<-GL344496$Ka
#GL344496Ka<-na.omit(GL344496$Ka)
#median(GL344496Ka)
#GL344496md <- function(x,y) {return(median(sample(GL344496Ka, 1, replace=TRUE)))}
#GL344496median <- boot(data=GL344496$Ka,statistic=GL344496md, R=10000)
#GL344496ci <- boot.ci(GL344496median, conf=0.95, type="norm")
#GL344496ci$norm

GL344539Ka<-GL344539$Ka
GL344539Ka<-na.omit(GL344539$Ka)
median(GL344539Ka)
GL344539md <- function(x,y) {return(median(sample(GL344539Ka, 1, replace=TRUE)))}
GL344539median <- boot(data=GL344539$Ka,statistic=GL344539md, R=10000)
GL344539ci <- boot.ci(GL344539median, conf=0.95, type="norm")
GL344539ci$norm

# Medain
LGbKa<-LGb$Ka
LGbKa<-na.omit(LGb$Ka)
median(LGbKa)
LGbmd <- function(x,y) {return(median(sample(LGbKa, 7, replace=TRUE)))}
LGbmedian <- boot(data=LGb$Ka,statistic=LGbmd, R=10000)
LGbci <- boot.ci(LGbmedian, conf=0.95, type="norm")
LGbci$norm

# Mean
LGbKa<-LGb$Ka
LGbKa<-na.omit(LGb$Ka)
mean(LGbKa)
LGbmd <- function(x,y) {return(mean(sample(LGbKa, 7, replace=TRUE)))}
LGbmean <- boot(data=LGb$Ka,statistic=LGbmd, R=10000)
LGbci <- boot.ci(LGbmean, conf=0.95, type="norm")
LGbci$norm

# putative X
# Median
pxKa<-na.omit(px$Ka)
median(pxKa)
pxKamd <- function(x,y) {return(median(sample(pxKa, 73, replace=TRUE)))}
pxKamedian<- boot(data=pxKa,statistic=pxKamd, R=10000)
pxKaci <- boot.ci(pxKamedian, conf=0.95, type="norm")
pxKaci$norm

# Mean
pxKa<-na.omit(px$Ka)
mean(pxKa)
pxKamd <- function(x,y) {return(mean(sample(pxKa, 73, replace=TRUE)))}
pxKamean<- boot(data=pxKa,statistic=pxKamd, R=10000)
pxKaci <- boot.ci(pxKamean, conf=0.95, type="norm")
pxKaci$norm

#--------------------------
# 5. Find the medians for Ks values with a 95% CI
#    Find the mean for Ks values with a 95% CI for autosomes, LGb, and putative X
#--------------------------

# contigs
# no value not included
AAWZ02037698Ks <-AAWZ02037698$Ks
AAWZ02037698Ks<-na.omit(AAWZ02037698$Ks)
median(AAWZ02037698Ks)
#AAWZ02037698md <- function(x,y) {return(median(sample(AAWZ02037698Ks, 0, replace=TRUE)))}
#AAWZ02037698median <- boot(data=AAWZ02037698$Ks,statistic=AAWZ02037698md, R=10000)
#AAWZ02037698ci <- boot.ci(AAWZ02037698median, conf=0.95, type="norm")
#AAWZ02037698ci$norm

AAWZ02040127Ks <-AAWZ02040127$Ks
AAWZ02040127Ks<-na.omit(AAWZ02040127$Ks)
median(AAWZ02040127Ks)
AAWZ02040127md <- function(x,y) {return(median(sample(AAWZ02040127Ks, 1, replace=TRUE)))}
AAWZ02040127median <- boot(data=AAWZ02040127$Ks,statistic=AAWZ02040127md, R=10000)
AAWZ02040127ci <- boot.ci(AAWZ02040127median, conf=0.95, type="norm")
AAWZ02040127ci$norm

# Has no value
AAWZ02041607Ks <-AAWZ02041607$Ks
AAWZ02041607Ks<-na.omit(AAWZ02041607$Ks)
median(AAWZ02041607Ks)
#AAWZ02041607md <- function(x,y) {return(median(sample(AAWZ02041607Ks, 1, replace=TRUE)))}
#AAWZ02041607median <- boot(data=AAWZ02041607$Ks,statistic=AAWZ02041607md, R=10000)
#AAWZ02041607ci <- boot.ci(AAWZ02041607median, conf=0.95, type="norm")
#AAWZ02041607ci$norm

Kscontigs<-c(AAWZ02037698Ks,AAWZ02040127Ks) #AAWZ02041607Ks
median(Kscontigs)
mdcontigs <- function(x,y) {return(median(sample(Kscontigs, 1, replace=TRUE)))}
mediancontigs <- boot(data=Kscontigs,statistic=mdcontigs, R=10000)
cicontigs <- boot.ci(mediancontigs, conf=0.95, type="norm")
cicontigs$norm

# autosomes 1-6
Chr1Ks<-Chr1$Ks
Chr1Ks<-na.omit(Chr1$Ks)
median(Chr1Ks)
Chr1md <- function(x,y) {return(median(sample(Chr1Ks, 590, replace=TRUE)))}
Chr1median <- boot(data=Chr1$Ks,statistic=Chr1md, R=10000)
Chr1ci <- boot.ci(Chr1median, conf=0.95, type="norm")
Chr1ci$norm

# doesn't have a header
Chr2Ks<-Chr2$Ks
Chr2Ks<-na.omit(Chr2$Ks)
median(Chr2Ks)
Chr2md <- function(x,y) {return(median(sample(Chr2Ks, 679, replace=TRUE)))}
Chr2median <- boot(data=Chr2$Ks,statistic=Chr2md, R=10000)
Chr2ci <- boot.ci(Chr2median, conf=0.95, type="norm")
Chr2ci$norm

Chr3Ks<-Chr3$Ks
Chr3Ks<-na.omit(Chr3$Ks)
median(Chr3Ks)
Chr3md <- function(x,y) {return(median(sample(Chr3Ks, 78, replace=TRUE)))}
Chr3median <- boot(data=Chr3$Ks,statistic=Chr3md, R=10000)
Chr3ci <- boot.ci(Chr3median, conf=0.95, type="norm")
Chr3ci$norm

Chr4Ks<-Chr4$Ks
Chr4Ks<-na.omit(Chr4$Ks)
median(Chr4Ks)
Chr4md <- function(x,y) {return(median(sample(Chr4Ks, 140, replace=TRUE)))}
Chr4median <- boot(data=Chr4$Ks,statistic=Chr4md, R=10000)
Chr4ci <- boot.ci(Chr4median, conf=0.95, type="norm")
Chr4ci$norm

Chr5Ks<-Chr5$Ks
Chr5Ks<-na.omit(Chr5$Ks)
median(Chr5Ks)
Chr5md <- function(x,y) {return(median(sample(Chr5Ks, 92, replace=TRUE)))}
Chr5median <- boot(data=Chr5$Ks,statistic=Chr5md, R=10000)
Chr5ci <- boot.ci(Chr5median, conf=0.95, type="norm")
Chr5ci$norm

Chr6Ks<-Chr6$Ks
Chr6Ks<-na.omit(Chr6$Ks)
median(Chr6Ks)
Chr6md <- function(x,y) {return(median(sample(Chr6Ks, 92, replace=TRUE)))}
Chr6median <- boot(data=Chr6$Ks,statistic=Chr6md, R=10000)
Chr6ci <- boot.ci(Chr6median, conf=0.95, type="norm")
Chr6ci$norm

# all autosomes 1-6
# Median
AKs<-na.omit(A$Ks)
median(AKs)
AKsmd<- function(x,y) {return(median(sample(AKs, 1853, replace=TRUE)))}
AKsmedian<- boot(data=AKs,statistic=AKsmd, R=10000)
AKsci <- boot.ci(AKsmedian, conf=0.95, type="norm")
AKsci$norm

# Mean
AKs<-na.omit(A$Ks)
mean(AKs)
AKsmd<- function(x,y) {return(mean(sample(AKs, 1853, replace=TRUE)))}
AKsmean<- boot(data=AKs,statistic=AKsmd, R=10000)
AKsci <- boot.ci(AKsmean, conf=0.95, type="norm")
AKsci$norm

# scaffolds
GL343282Ks<-GL343282$Ks
GL343282Ks<-na.omit(GL343282$Ks)
median(GL343282Ks)
GL343282md <- function(x,y) {return(median(sample(GL343282Ks, 9, replace=TRUE)))}
GL343282median <- boot(data=GL343282$Ks,statistic=GL343282md, R=10000)
GL343282ci <- boot.ci(GL343282median, conf=0.95, type="norm")
GL343282ci$norm

GL343338Ks<-GL343338$Ks
GL343338Ks<-na.omit(GL343338$Ks)
median(GL343338Ks)
GL343338md <- function(x,y) {return(median(sample(GL343338Ks, 12, replace=TRUE)))}
GL343338median <- boot(data=GL343338$Ks,statistic=GL343338md, R=10000)
GL343338ci <- boot.ci(GL343338median, conf=0.95, type="norm")
GL343338ci$norm

GL343364Ks<-GL343364$Ks
GL343364Ks<-na.omit(GL343364$Ks)
median(GL343364Ks)
GL343364md <- function(x,y) {return(median(sample(GL343364Ks, 7, replace=TRUE)))}
GL343364median <- boot(data=GL343364$Ks,statistic=GL343364md, R=10000)
GL343364ci <- boot.ci(GL343364median, conf=0.95, type="norm")
GL343364ci$norm

GL343417Ks<-GL343417$Ks
GL343417Ks<-na.omit(GL343417$Ks)
median(GL343417Ks)
GL343417md <- function(x,y) {return(median(sample(GL343417Ks, 3, replace=TRUE)))}
GL343417median <- boot(data=GL343417$Ks,statistic=GL343417md, R=10000)
GL343417ci <- boot.ci(GL343417median, conf=0.95, type="norm")
GL343417ci$norm

GL343422Ks<-GL343422$Ks
GL343422Ks<-na.omit(GL343422$Ks)
median(GL343422Ks)
GL343422md <- function(x,y) {return(median(sample(GL343422Ks, 3, replace=TRUE)))}
GL343422median <- boot(data=GL343422$Ks,statistic=GL343422md, R=10000)
GL343422ci <- boot.ci(GL343422median, conf=0.95, type="norm")
GL343422ci$norm

GL343423Ks<-GL343423$Ks
GL343423Ks<-na.omit(GL343423$Ks)
median(GL343423Ks)
GL343423md <- function(x,y) {return(median(sample(GL343423Ks, 4, replace=TRUE)))}
GL343423median <- boot(data=GL343423$Ks,statistic=GL343423md, R=10000)
GL343423ci <- boot.ci(GL343423median, conf=0.95, type="norm")
GL343423ci$norm

GL343439Ks<-GL343439$Ks
GL343439Ks<-na.omit(GL343439$Ks)
median(GL343439Ks)
GL343439md <- function(x,y) {return(median(sample(GL343439Ks, 4, replace=TRUE)))}
GL343439median <- boot(data=GL343439$Ks,statistic=GL343439md, R=10000)
GL343439ci <- boot.ci(GL343439median, conf=0.95, type="norm")
GL343439ci$norm

GL343462Ks<-GL343462$Ks
GL343462Ks<-na.omit(GL343462$Ks)
median(GL343462Ks)
GL343462md <- function(x,y) {return(median(sample(GL343462Ks, 1, replace=TRUE)))}
GL343462median <- boot(data=GL343462$Ks,statistic=GL343462md, R=10000)
GL343462ci <- boot.ci(GL343462median, conf=0.95, type="norm")
GL343462ci$norm

GL343516Ks<-GL343516$Ks
GL343516Ks<-na.omit(GL343516$Ks)
median(GL343516Ks)
GL343516md <- function(x,y) {return(median(sample(GL343516Ks, 2, replace=TRUE)))}
GL343516median <- boot(data=GL343516$Ks,statistic=GL343516md, R=10000)
GL343516ci <- boot.ci(GL343516median, conf=0.95, type="norm")
GL343516ci$norm

GL343525Ks<-GL343525$Ks
GL343525Ks<-na.omit(GL343525$Ks)
median(GL343525Ks)
GL343525md <- function(x,y) {return(median(sample(GL343525Ks, 1, replace=TRUE)))}
GL343525median <- boot(data=GL343525$Ks,statistic=GL343525md, R=10000)
GL343525ci <- boot.ci(GL343525median, conf=0.95, type="norm")
GL343525ci$norm

GL343550Ks<-GL343550$Ks
GL343550Ks<-na.omit(GL343550$Ks)
median(GL343550Ks)
GL343550md <- function(x,y) {return(median(sample(GL343550Ks, 7, replace=TRUE)))}
GL343550median <- boot(data=GL343550$Ks,statistic=GL343550md, R=10000)
GL343550ci <- boot.ci(GL343550median, conf=0.95, type="norm")
GL343550ci$norm

# Has no value
#GL343588Ks<-GL343588$Ks
#GL343588Ks<-na.omit(GL343588$Ks)
#median(GL343588Ks)
#GL343588md <- function(x,y) {return(median(sample(GL343588Ks, 1, replace=TRUE)))}
#GL343588median <- boot(data=GL343588$Ks,statistic=GL343588md, R=10000)
#GL343588ci <- boot.ci(GL343588median, conf=0.95, type="norm")
#GL343588ci$norm

# Has no value
#GL343731Ks<-GL343731$Ks
#GL343731Ks<-na.omit(GL343731$Ks)
#median(GL343731Ks)
#GL343731md <- function(x,y) {return(median(sample(GL343731Ks, 1, replace=TRUE)))}
#GL343731median <- boot(data=GL343731$Ks,statistic=GL343731md, R=10000)
#GL343731ci <- boot.ci(GL343731median, conf=0.95, type="norm")
#GL343731ci$norm

GL343947Ks<-GL343947$Ks
GL343947Ks<-na.omit(GL343947$Ks)
median(GL343947Ks)
GL343947md <- function(x,y) {return(median(sample(GL343947Ks, 4, replace=TRUE)))}
GL343947median <- boot(data=GL343947$Ks,statistic=GL343947md, R=10000)
GL343947ci <- boot.ci(GL343947median, conf=0.95, type="norm")
GL343947ci$norm

GL344042Ks<-GL344042$Ks
GL344042Ks<-na.omit(GL344042$Ks)
median(GL344042Ks)
GL344042md <- function(x,y) {return(median(sample(GL344042Ks, 1, replace=TRUE)))}
GL344042median <- boot(data=GL344042$Ks,statistic=GL344042md, R=10000)
GL344042ci <- boot.ci(GL344042median, conf=0.95, type="norm")
GL344042ci$norm

# Has no value
#GL344393Ks<-GL344393$Ks
#GL344393Ks<-na.omit(GL344393$Ks)
#median(GL344393Ks)
#GL344393md <- function(x,y) {return(median(sample(GL344393Ks, 1, replace=TRUE)))}
#GL344393median <- boot(data=GL344393$Ks,statistic=GL344393md, R=10000)
#GL344393ci <- boot.ci(GL344393median, conf=0.95, type="norm")
#GL344393ci$norm

# Has no value
#GL344496Ks<-GL344496$Ks
#GL344496Ks<-na.omit(GL344496$Ks)
#median(GL344496Ks)
#GL344496md <- function(x,y) {return(median(sample(GL344496Ks, 1, replace=TRUE)))}
#GL344496median <- boot(data=GL344496$Ks,statistic=GL344496md, R=10000)
#GL344496ci <- boot.ci(GL344496median, conf=0.95, type="norm")
#GL344496ci$norm

GL344539Ks<-GL344539$Ks
GL344539Ks<-na.omit(GL344539$Ks)
median(GL344539Ks)
GL344539md <- function(x,y) {return(median(sample(GL344539Ks, 1, replace=TRUE)))}
GL344539median <- boot(data=GL344539$Ks,statistic=GL344539md, R=10000)
GL344539ci <- boot.ci(GL344539median, conf=0.95, type="norm")
GL344539ci$norm

# Median
LGbKs<-LGb$Ks
LGbKs<-na.omit(LGb$Ks)
median(LGbKs)
LGbmd <- function(x,y) {return(median(sample(LGbKs, 7, replace=TRUE)))}
LGbmedian <- boot(data=LGb$Ks,statistic=LGbmd, R=10000)
LGbci <- boot.ci(LGbmedian, conf=0.95, type="norm")
LGbci$norm

# Mean
LGbKs<-LGb$Ks
LGbKs<-na.omit(LGb$Ks)
mean(LGbKs)
LGbmd <- function(x,y) {return(mean(sample(LGbKs, 7, replace=TRUE)))}
LGbmean <- boot(data=LGb$Ks,statistic=LGbmd, R=10000)
LGbci <- boot.ci(LGbmean, conf=0.95, type="norm")
LGbci$norm

# putative X
# Median
pxKs<-na.omit(px$Ks)
median(pxKs)
pxKsmd <- function(x,y) {return(median(sample(pxKs, 69, replace=TRUE)))}
pxKsmedian<- boot(data=pxKs,statistic=pxKsmd, R=10000)
pxKsci <- boot.ci(pxKsmedian, conf=0.95, type="norm")
pxKsci$norm

# Mean
pxKs<-na.omit(px$Ks)
mean(pxKs)
pxKsmd <- function(x,y) {return(mean(sample(pxKs, 69, replace=TRUE)))}
pxKsmean<- boot(data=pxKs,statistic=pxKsmd, R=10000)
pxKsci <- boot.ci(pxKsmean, conf=0.95, type="norm")
pxKsci$norm

#--------------------------
# 5. Find the medians for Ka.Ks values with a 95% CI
#    Find the mean for Ka.Ks values with a 95% CI for autosomes, LGb, and putative X
#--------------------------

#
# Has no value
#AAWZ02037698KaKs <-AAWZ02037698$Ka.Ks
#AAWZ02037698KaKs<-na.omit(AAWZ02037698$Ka.Ks)
#median(AAWZ02037698KaKs)
#AAWZ02037698md <- function(x,y) {return(median(sample(AAWZ02037698KaKs, 2, replace=TRUE)))}
#AAWZ02037698median <- boot(data=AAWZ02037698$Ka.Ks,statistic=AAWZ02037698md, R=10000)
#AAWZ02037698ci <- boot.ci(AAWZ02037698median, conf=0.95, type="norm")
#AAWZ02037698ci$norm

AAWZ02040127KaKs <-AAWZ02040127$Ka.Ks
AAWZ02040127KaKs<-na.omit(AAWZ02040127$Ka.Ks)
median(AAWZ02040127KaKs)
AAWZ02040127md <- function(x,y) {return(median(sample(AAWZ02040127KaKs, 1, replace=TRUE)))}
AAWZ02040127median <- boot(data=AAWZ02040127$Ka.Ks,statistic=AAWZ02040127md, R=10000)
AAWZ02040127ci <- boot.ci(AAWZ02040127median, conf=0.95, type="norm")
AAWZ02040127ci$norm

# Has no value
#AAWZ02041607KaKs <-AAWZ02041607$Ka.Ks
#AAWZ02041607KaKs<-na.omit(AAWZ02041607$Ka.Ks)
#median(AAWZ02041607KaKs)
#AAWZ02041607md <- function(x,y) {return(median(sample(AAWZ02041607KaKs, 1, replace=TRUE)))}
#AAWZ02041607median <- boot(data=AAWZ02041607$Ka.Ks,statistic=AAWZ02041607md, R=10000)
#AAWZ02041607ci <- boot.ci(AAWZ02041607median, conf=0.95, type="norm")
#AAWZ02041607ci$norm

KaKscontigs<-c(AAWZ02040127KaKs) #AAWZ02041607Ks #AAWZ02037698KaKs,
median(KaKscontigs)
mdcontigs <- function(x,y) {return(median(sample(KaKscontigs, 1, replace=TRUE)))}
mediancontigs <- boot(data=KaKscontigs,statistic=mdcontigs, R=10000)
cicontigs <- boot.ci(mediancontigs, conf=0.95, type="norm")
cicontigs$norm

# autosomes 1-6
Chr1KaKs<-Chr1$Ka.Ks
Chr1KaKs<-na.omit(Chr1$Ka.Ks)
median(Chr1KaKs)
Chr1md <- function(x,y) {return(median(sample(Chr1KaKs, 590, replace=TRUE)))}
Chr1median <- boot(data=Chr1$Ka.Ks,statistic=Chr1md, R=10000)
Chr1ci <- boot.ci(Chr1median, conf=0.95, type="norm")
Chr1ci$norm

# doesn't have a header
Chr2KaKs<-Chr2$Ka.Ks
Chr2KaKs<-na.omit(Chr2$Ka.Ks)
median(Chr2KaKs)
Chr2md <- function(x,y) {return(median(sample(Chr2KaKs, 628, replace=TRUE)))}
Chr2median <- boot(data=Chr2$Ka.Ks,statistic=Chr2md, R=10000)
Chr2ci <- boot.ci(Chr2median, conf=0.95, type="norm")
Chr2ci$norm

Chr3KaKs<-Chr3$Ka.Ks
Chr3KaKs<-na.omit(Chr3$Ka.Ks)
median(Chr3KaKs)
Chr3md <- function(x,y) {return(median(sample(Chr3KaKs, 78, replace=TRUE)))}
Chr3median <- boot(data=Chr3$Ka.Ks,statistic=Chr3md, R=10000)
Chr3ci <- boot.ci(Chr3median, conf=0.95, type="norm")
Chr3ci$norm

Chr4KaKs<-Chr4$Ka.Ks
Chr4KaKs<-na.omit(Chr4$Ka.Ks)
median(Chr4KaKs)
Chr4md <- function(x,y) {return(median(sample(Chr4KaKs, 140, replace=TRUE)))}
Chr4median <- boot(data=Chr4$Ka.Ks,statistic=Chr4md, R=10000)
Chr4ci <- boot.ci(Chr4median, conf=0.95, type="norm")
Chr4ci$norm

Chr5KaKs<-Chr5$Ka.Ks
Chr5KaKs<-na.omit(Chr5$Ka.Ks)
median(Chr5KaKs)
Chr5md <- function(x,y) {return(median(sample(Chr5KaKs, 92, replace=TRUE)))}
Chr5median <- boot(data=Chr5$Ka.Ks,statistic=Chr5md, R=10000)
Chr5ci <- boot.ci(Chr5median, conf=0.95, type="norm")
Chr5ci$norm

Chr6KaKs<-Chr6$Ka.Ks
Chr6KaKs<-na.omit(Chr6$Ka.Ks)
median(Chr6KaKs)
Chr6md <- function(x,y) {return(median(sample(Chr6KaKs, 92, replace=TRUE)))}
Chr6median <- boot(data=Chr6$Ka.Ks,statistic=Chr6md, R=10000)
Chr6ci <- boot.ci(Chr6median, conf=0.95, type="norm")
Chr6ci$norm

# all autosomes 1-6
# Medain
AKa.Ks<-na.omit(A$Ka.Ks)
median(AKa.Ks)
AKa.Ksmd<- function(x,y) {return(median(sample(AKa.Ks, 1845, replace=TRUE)))}
AKa.Ksmedian<- boot(data=AKa.Ks,statistic=AKa.Ksmd, R=10000)
AKa.Ksci <- boot.ci(AKa.Ksmedian, conf=0.95, type="norm")
AKa.Ksci$norm

# Mean
AKa.Ks<-na.omit(A$Ka.Ks)
mean(AKa.Ks)
AKa.Ksmd<- function(x,y) {return(mean(sample(AKa.Ks, 1845, replace=TRUE)))}
AKa.Ksmean<- boot(data=AKa.Ks,statistic=AKa.Ksmd, R=10000)
AKa.Ksci <- boot.ci(AKa.Ksmean, conf=0.95, type="norm")
AKa.Ksci$norm

# scaffolds
GL343282KaKs<-GL343282$Ka.Ks
GL343282KaKs<-na.omit(GL343282$Ka.Ks)
median(GL343282KaKs)
GL343282md <- function(x,y) {return(median(sample(GL343282KaKs, 9, replace=TRUE)))}
GL343282median <- boot(data=GL343282$Ka.Ks,statistic=GL343282md, R=10000)
GL343282ci <- boot.ci(GL343282median, conf=0.95, type="norm")
GL343282ci$norm

GL343338KaKs<-GL343338$Ka.Ks
GL343338KaKs<-na.omit(GL343338$Ka.Ks)
median(GL343338KaKs)
GL343338md <- function(x,y) {return(median(sample(GL343338KaKs, 12, replace=TRUE)))}
GL343338median <- boot(data=GL343338$Ka.Ks,statistic=GL343338md, R=10000)
GL343338ci <- boot.ci(GL343338median, conf=0.95, type="norm")
GL343338ci$norm

GL343364KaKs<-GL343364$Ka.Ks
GL343364KaKs<-na.omit(GL343364$Ka.Ks)
median(GL343364KaKs)
GL343364md <- function(x,y) {return(median(sample(GL343364KaKs, 7, replace=TRUE)))}
GL343364median <- boot(data=GL343364$Ka.Ks,statistic=GL343364md, R=10000)
GL343364ci <- boot.ci(GL343364median, conf=0.95, type="norm")
GL343364ci$norm

GL343417KaKs<-GL343417$Ka.Ks
GL343417KaKs<-na.omit(GL343417$Ka.Ks)
median(GL343417KaKs)
GL343417md <- function(x,y) {return(median(sample(GL343417KaKs, 3, replace=TRUE)))}
GL343417median <- boot(data=GL343417$Ka.Ks,statistic=GL343417md, R=10000)
GL343417ci <- boot.ci(GL343417median, conf=0.95, type="norm")
GL343417ci$norm

GL343422KaKs<-GL343422$Ka.Ks
GL343422KaKs<-na.omit(GL343422$Ka.Ks)
median(GL343422KaKs)
GL343422md <- function(x,y) {return(median(sample(GL343422KaKs, 3, replace=TRUE)))}
GL343422median <- boot(data=GL343422$Ka.Ks,statistic=GL343422md, R=10000)
GL343422ci <- boot.ci(GL343422median, conf=0.95, type="norm")
GL343422ci$norm

GL343423KaKs<-GL343423$Ka.Ks
GL343423KaKs<-na.omit(GL343423$Ka.Ks)
median(GL343423KaKs)
GL343423md <- function(x,y) {return(median(sample(GL343423KaKs, 4, replace=TRUE)))}
GL343423median <- boot(data=GL343423$Ka.Ks,statistic=GL343423md, R=10000)
GL343423ci <- boot.ci(GL343423median, conf=0.95, type="norm")
GL343423ci$norm

GL343439KaKs<-GL343439$Ka.Ks
GL343439KaKs<-na.omit(GL343439$Ka.Ks)
median(GL343439KaKs)
GL343439md <- function(x,y) {return(median(sample(GL343439KaKs, 4, replace=TRUE)))}
GL343439median <- boot(data=GL343439$Ka.Ks,statistic=GL343439md, R=10000)
GL343439ci <- boot.ci(GL343439median, conf=0.95, type="norm")
GL343439ci$norm

GL343462KaKs<-GL343462$Ka.Ks
GL343462KaKs<-na.omit(GL343462$Ka.Ks)
median(GL343462KaKs)
GL343462md <- function(x,y) {return(median(sample(GL343462KaKs, 1, replace=TRUE)))}
GL343462median <- boot(data=GL343462$Ka.Ks,statistic=GL343462md, R=10000)
GL343462ci <- boot.ci(GL343462median, conf=0.95, type="norm")
GL343462ci$norm

GL343516KaKs<-GL343516$Ka.Ks
GL343516KaKs<-na.omit(GL343516$Ka.Ks)
median(GL343516KaKs)
GL343516md <- function(x,y) {return(median(sample(GL343516KaKs, 2, replace=TRUE)))}
GL343516median <- boot(data=GL343516$Ka.Ks,statistic=GL343516md, R=10000)
GL343516ci <- boot.ci(GL343516median, conf=0.95, type="norm")
GL343516ci$norm

GL343525KaKs<-GL343525$Ka.Ks
GL343525KaKs<-na.omit(GL343525$Ka.Ks)
median(GL343525KaKs)
GL343525md <- function(x,y) {return(median(sample(GL343525KaKs, 1, replace=TRUE)))}
GL343525median <- boot(data=GL343525$Ka.Ks,statistic=GL343525md, R=10000)
GL343525ci <- boot.ci(GL343525median, conf=0.95, type="norm")
GL343525ci$norm

GL343550KaKs<-GL343550$Ka.Ks
GL343550KaKs<-na.omit(GL343550$Ka.Ks)
median(GL343550KaKs)
GL343550md <- function(x,y) {return(median(sample(GL343550KaKs, 7, replace=TRUE)))}
GL343550median <- boot(data=GL343550$Ka.Ks,statistic=GL343550md, R=10000)
GL343550ci <- boot.ci(GL343550median, conf=0.95, type="norm")
GL343550ci$norm

# Has no value
#GL343588KaKs<-GL343588$Ka.Ks
#GL343588KaKs<-na.omit(GL343588$Ka.Ks)
#median(GL343588KaKs)
#GL343588md <- function(x,y) {return(median(sample(GL343588KaKs, 1, replace=TRUE)))}
#GL343588median <- boot(data=GL343588$Ka.Ks,statistic=GL343588md, R=10000)
#GL343588ci <- boot.ci(GL343588median, conf=0.95, type="norm")
#GL343588ci$norm

# Has no value
#GL343731KaKs<-GL343731$Ka.Ks
#GL343731KaKs<-na.omit(GL343731$Ka.Ks)
#median(GL343731KaKs)
#GL343731md <- function(x,y) {return(median(sample(GL343731KaKs, 1, replace=TRUE)))}
#GL343731median <- boot(data=GL343731$Ka.Ks,statistic=GL343731md, R=10000)
#GL343731ci <- boot.ci(GL343731median, conf=0.95, type="norm")
#GL343731ci$norm

GL343947KaKs<-GL343947$Ka.Ks
GL343947KaKs<-na.omit(GL343947$Ka.Ks)
median(GL343947KaKs)
GL343947md <- function(x,y) {return(median(sample(GL343947KaKs, 4, replace=TRUE)))}
GL343947median <- boot(data=GL343947$Ka.Ks,statistic=GL343947md, R=10000)
GL343947ci <- boot.ci(GL343947median, conf=0.95, type="norm")
GL343947ci$norm

GL344042KaKs<-GL344042$Ka.Ks
GL344042KaKs<-na.omit(GL344042$Ka.Ks)
median(GL344042KaKs)
GL344042md <- function(x,y) {return(median(sample(GL344042KaKs, 1, replace=TRUE)))}
GL344042median <- boot(data=GL344042$Ka.Ks,statistic=GL344042md, R=10000)
GL344042ci <- boot.ci(GL344042median, conf=0.95, type="norm")
GL344042ci$norm

# Has no value
#GL344393KaKs<-GL344393$Ka.Ks
#GL344393KaKs<-na.omit(GL344393$Ka.Ks)
#median(GL344393KaKs)
#GL344393md <- function(x,y) {return(median(sample(GL344393KaKs, 1, replace=TRUE)))}
#GL344393median <- boot(data=GL344393$Ka.Ks,statistic=GL344393md, R=10000)
#GL344393ci <- boot.ci(GL344393median, conf=0.95, type="norm")
#GL344393ci$norm

# Has no value
#GL344496KaKs<-GL344496$Ka.Ks
#GL344496KaKs<-na.omit(GL344496$Ka.Ks)
#median(GL344496KaKs)
#GL344496md <- function(x,y) {return(median(sample(GL344496KaKs, 1, replace=TRUE)))}
#GL344496median <- boot(data=GL344496$Ka.Ks,statistic=GL344496md, R=10000)
#GL344496ci <- boot.ci(GL344496median, conf=0.95, type="norm")
#GL344496ci$norm

GL344539KaKs<-GL344539$Ka.Ks
GL344539KaKs<-na.omit(GL344539$Ka.Ks)
median(GL344539KaKs)
GL344539md <- function(x,y) {return(median(sample(GL344539KaKs, 1, replace=TRUE)))}
GL344539median <- boot(data=GL344539$Ka.Ks,statistic=GL344539md, R=10000)
GL344539ci <- boot.ci(GL344539median, conf=0.95, type="norm")
GL344539ci$norm

# Medain
LGbKaKs<-LGb$Ka.Ks
LGbKaKs<-na.omit(LGb$Ka.Ks)
median(LGbKaKs)
LGbmd <- function(x,y) {return(median(sample(LGbKaKs, 7, replace=TRUE)))}
LGbmedian <- boot(data=LGb$Ka.Ks,statistic=LGbmd, R=10000)
LGbci <- boot.ci(LGbmedian, conf=0.95, type="norm")
LGbci$norm

# Mean
mean(LGbKaKs)
LGbmd <- function(x,y) {return(mean(sample(LGbKaKs, 7, replace=TRUE)))}
LGbmean <- boot(data=LGb$Ka.Ks,statistic=LGbmd, R=10000)
LGbci <- boot.ci(LGbmean, conf=0.95, type="norm")
LGbci$norm

# putative X
# Medain
pxKa.Ks<-na.omit(px$Ka.Ks)
median(pxKa.Ks)
pxKa.Ksmd <- function(x,y) {return(median(sample(pxKa.Ks, 68, replace=TRUE)))}
pxKa.Ksmedian<- boot(data=pxKa.Ks,statistic=pxKa.Ksmd, R=10000)
pxKa.Ksci <- boot.ci(pxKa.Ksmedian, conf=0.95, type="norm")
pxKa.Ksci$norm

# Mean
mean(pxKa.Ks)
pxKa.Ksmd <- function(x,y) {return(mean(sample(pxKa.Ks, 68, replace=TRUE)))}
pxKa.Ksmean<- boot(data=pxKa.Ks,statistic=pxKa.Ksmd, R=10000)
pxKa.Ksci <- boot.ci(pxKa.Ksmean, conf=0.95, type="norm")
pxKa.Ksci$norm

#--------------------------
# 6. Median Ka.Ks boxplot
#--------------------------

boxplot(AKa.Ks,pxKa.Ks,LGbKaKs,
col=c("firebrick3","green4","blue3"),
names=c("Autosomes KaKs","Putative X KaKs","LGb KaKs"),
ylab = "KaKs",
main="Autosomes, Putative X, LGb KaKs values")
