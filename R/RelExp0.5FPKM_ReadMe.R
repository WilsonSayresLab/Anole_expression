###############################################################
#-------------- OverView ----------------
#
#  Male to female median and CI calculations
#	for anole genome relative expression with minimum FPKM of 0.5.
#
#	Required Programs:	R
#				              Boot
#############################################################

#--------------------------
# 1. Remove genes with an FPKM below 0.5 for either male or female expression
#--------------------------

# Removed 5503 genes
# Keep all extra data for further sorting

#--------------------------
# 2. Import data into R and subset by locus
#--------------------------

# Autosomes
mf<-read.csv("MF_GeneExp0.5FPKM.csv",header=TRUE)

chr1<-subset(mf,scaffold=="1",select=c(rel_exp))
chr2<-subset(mf,scaffold=="2",select=c(rel_exp))
chr3<-subset(mf,scaffold=="3",select=c(rel_exp))
chr4<-subset(mf,scaffold=="4",select=c(rel_exp))
chr5<-subset(mf,scaffold=="5",select=c(rel_exp))
chr6<-subset(mf,scaffold=="6",select=c(rel_exp))
auto<-c(chr1,chr2,chr3,chr4,chr5,chr6)

# X-Linked Genes (NAs are omitted after susetting to preserve genes which may be missing data in other columns)
chrx<-read.csv("PropX0.5FPKM.csv",header=TRUE)

propx<-na.omit(chrx$RelativeExpression)

lgb<-subset(chrx,Scaffold=="LGb",select=RelativeExpression)
lgb<-na.omit(lgb)

scaf<-subset(chrx,Scaffold!="LGb")
gg15<-subset(scaf,GG_Chr=="15",select=RelativeExpression)
gg15<-na.omit(gg15)
unmapped<-subset(scaf,GG_Chr=="nomatch",select=RelativeExpression)
unmapped<-na.omit(unmapped)

gg1<-subset(scaf,GG_Chr=="1",select=RelativeExpression)
gg2<-subset(scaf,GG_Chr=="2",select=RelativeExpression)
gg19<-subset(scaf,GG_Chr=="19",select=RelativeExpression)
other<-c(gg1,gg2,gg19)

#--------------------------
# 3. Calculate Median Relative Expression and bootstrap p-values
#--------------------------

library("boot")

# Autosomes
median(chr1$rel_exp)
m1 <- function(x,y) {return(median(sample(chr1$rel_exp, 300, replace=TRUE)))}
median1 <- boot(data=chr1$rel_exp,statistic=m1, R=1000)
ci1 <- boot.ci(median1, conf=0.95, type="norm")
ci1$norm

median(chr2$rel_exp)
m2 <- function(x,y) {return(median(sample(chr2$rel_exp, 300, replace=TRUE)))}
median2 <- boot(data=chr2$rel_exp,statistic=m2, R=1000)
ci2 <- boot.ci(median2, conf=0.95, type="norm")
ci2$norm

median(chr3$rel_exp)
m3 <- function(x,y) {return(median(sample(chr3$rel_exp, 300, replace=TRUE)))}
median3 <- boot(data=chr3$rel_exp,statistic=m3, R=1000)
ci3 <- boot.ci(median3, conf=0.95, type="norm")
ci3$norm

median(chr4$rel_exp)
m4 <- function(x,y) {return(median(sample(chr4$rel_exp, 300, replace=TRUE)))}
median4 <- boot(data=chr4$rel_exp,statistic=m4, R=1000)
ci4 <- boot.ci(median4, conf=0.95, type="norm")
ci4$norm

median(chr5$rel_exp)
m5 <- function(x,y) {return(median(sample(chr5$rel_exp, 300, replace=TRUE)))}
median5 <- boot(data=chr5$rel_exp,statistic=m5, R=1000)
ci5 <- boot.ci(median5, conf=0.95, type="norm")
ci5$norm

median(chr6$rel_exp)
m6 <- function(x,y) {return(median(sample(chr6$rel_exp, 300, replace=TRUE)))}
median6 <- boot(data=chr6$rel_exp,statistic=m6, R=1000)
ci6 <- boot.ci(median6, conf=0.95, type="norm")
ci6$norm

median(auto$rel_exp)
ma<- function(x,y) {return(median(sample(auto$rel_exp, 300, replace=TRUE)))}
mediana <- boot(data=auto$rel_exp,statistic=ma, R=1000)
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
