###############################################################
#-------------- OverView ----------------
#
#  Male to female median and CI calculations
#	for anole genome relative expression based on mapping
#	to the chicken genome.
#
#	Required Programs:	R
#				Boot
#############################################################

#--------------------------
# 1. Import data into R and subset by locus
#--------------------------

# Autosomes
mf<-read.csv("/home/yitzhak/Documents/AnoleDifExp/RelExp/anoleRelExp_sra_fpkm2.csv",header=TRUE)

chr1<-subset(mf,locus=="1",select=c(relativeExpression))
chr2<-subset(mf,locus=="2",select=c(relativeExpression))
chr3<-subset(mf,locus=="3",select=c(relativeExpression))
chr4<-subset(mf,locus=="4",select=c(relativeExpression))
chr5<-subset(mf,locus=="5",select=c(relativeExpression))
chr6<-subset(mf,locus=="6",select=c(relativeExpression))
auto<-subset(mf,locus=="1"|locus=="2"|locus=="3"|locus=="4"|locus=="5"|locus=="6",select=relativeExpression)

# X-Linked Genes (NAs are omitted after susetting to preserve genes which may be missing data in other columns)
chrx<-read.csv("/home/yitzhak/Documents/AnoleDifExp/RelExp/XLinkedGenes_fpkm2.csv",header=TRUE)

propx<-na.omit(chrx$RelativeExpression)

lgb<-subset(chrx,locus=="LGb",select=RelativeExpression)
lgb<-na.omit(lgb)

scaf<-subset(chrx,locus!="LGb")
gg15<-subset(scaf,GG_Chr=="15",select=RelativeExpression)
gg15<-na.omit(gg15)
unmapped<-subset(chrx,GG_Chr=="unmapped",select=RelativeExpression)
unmapped<-na.omit(unmapped)

other<-subset(chrx,GG_Chr=="1"|GG_Chr=="2"|GG_Chr=="19",select=RelativeExpression)

#--------------------------
# 2. Calculate Median Relative Expression and bootstrap p-values
#--------------------------

library("boot")

# Autosomes
median(chr1$relativeExpression)
m1 <- function(x,y) {return(median(sample(chr1$relativeExpression, 300, replace=TRUE)))}
median1 <- boot(data=chr1$relativeExpression,statistic=m1, R=1000)
ci1 <- boot.ci(median1, conf=0.95, type="norm")
ci1$norm

median(chr2$relativeExpression)
m2 <- function(x,y) {return(median(sample(chr2$relativeExpression, 300, replace=TRUE)))}
median2 <- boot(data=chr2$relativeExpression,statistic=m2, R=1000)
ci2 <- boot.ci(median2, conf=0.95, type="norm")
ci2$norm

median(chr3$relativeExpression)
m3 <- function(x,y) {return(median(sample(chr3$relativeExpression, 300, replace=TRUE)))}
median3 <- boot(data=chr3$relativeExpression,statistic=m3, R=1000)
ci3 <- boot.ci(median3, conf=0.95, type="norm")
ci3$norm

median(chr4$relativeExpression)
m4 <- function(x,y) {return(median(sample(chr4$relativeExpression, 300, replace=TRUE)))}
median4 <- boot(data=chr4$relativeExpression,statistic=m4, R=1000)
ci4 <- boot.ci(median4, conf=0.95, type="norm")
ci4$norm

median(chr5$relativeExpression)
m5 <- function(x,y) {return(median(sample(chr5$relativeExpression, 300, replace=TRUE)))}
median5 <- boot(data=chr5$relativeExpression,statistic=m5, R=1000)
ci5 <- boot.ci(median5, conf=0.95, type="norm")
ci5$norm

median(chr6$relativeExpression)
m6 <- function(x,y) {return(median(sample(chr6$relativeExpression, 300, replace=TRUE)))}
median6 <- boot(data=chr6$relativeExpression,statistic=m6, R=1000)
ci6 <- boot.ci(median6, conf=0.95, type="norm")
ci6$norm

median(auto$relativeExpression)
ma<- function(x,y) {return(median(sample(auto$relativeExpression, 300, replace=TRUE)))}
mediana <- boot(data=auto$relativeExpression,statistic=ma, R=1000)
cia <- boot.ci(mediana, conf=0.95, type="norm")
cia$norm

# X-linked genes
median(lgb$RelativeExpression)
ml <- function(x,y) {return(median(sample(lgb$RelativeExpression, 300, replace=TRUE)))}
medianl <- boot(data=lgb$RelativeExpression,statistic=ml, R=1000)
cil <- boot.ci(medianl, conf=0.95, type="norm")
cil$norm

median(gg15$RelativeExpression)
mg <- function(x,y) {return(median(sample(gg15$RelativeExpression, 300, replace=TRUE)))}
mediang <- boot(data=gg15$RelativeExpression,statistic=mg, R=1000)
cig <- boot.ci(mediang, conf=0.95, type="norm")
cig$norm

median(other$RelativeExpression)
mo <- function(x,y) {return(median(sample(other$RelativeExpression, 300, replace=TRUE)))}
mediano <- boot(data=other$RelativeExpression,statistic=mo, R=1000)
cio <- boot.ci(mediano, conf=0.95, type="norm")
cio$norm

median(unmapped$RelativeExpression)
mu <- function(x,y) {return(median(sample(unmapped$RelativeExpression, 300, replace=TRUE)))}
medianu <- boot(data=unmapped$RelativeExpression,statistic=mu, R=1000)
ciu <- boot.ci(medianu, conf=0.95, type="norm")
ciu$norm

median(propx)
mp <- function(x,y) {return(median(sample(propx, 300, replace=TRUE)))}
medianp <- boot(data=propx,statistic=mp, R=1000)
cip <- boot.ci(medianp, conf=0.95, type="norm")
cip$norm
