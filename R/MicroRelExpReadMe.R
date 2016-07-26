###############################################################
#-------------- OverView ----------------
#
#  Male to female median and CI calculations
#	for anole genome relative expression with a minimum 
#	FPKM of 2.
#
#	Required Programs:	R
#				Boot
#############################################################

#--------------------------
# 1. Import data into R and subset by locus 
# (NAs are omitted after susetting to preserve genes which may be missing data in other columns)
#--------------------------

mf<-read.csv("/home/yitzhak/Documents/AnoleDifExp/RelExp/anoleRelExp_sra_fpkm2.csv",header=TRUE)

lga<-subset(mf,locus=="LGa",select=relativeExpression)
lga<-na.omit(lga)
lgb<-subset(mf,locus=="LGb",select=relativeExpression)
lgb<-na.omit(lgb)
lgc<-subset(mf,locus=="LGc",select=relativeExpression)
lgc<-na.omit(lgc)
lgd<-subset(mf,locus=="LGd",select=relativeExpression)
lgd<-na.omit(lgd)
lgf<-subset(mf,locus=="LGf",select=relativeExpression)
lgf<-na.omit(lgf)
lgg<-subset(mf,locus=="LGg",select=relativeExpression)
lgg<-na.omit(lgg)
lgh<-subset(mf,locus=="LGh",select=relativeExpression)
lgh<-na.omit(lgh)
micro<-subset(mf,locus=="LGa"|locus=="LGc"|locus=="LGf"|locus=="LGg"|locus=="LGh",select=relativeExpression)
micro<-na.omit(micro)

auto<-subset(mf,locus=="1"|locus=="2"|locus=="3"|locus=="4"|locus=="5"|locus=="6",select=relativeExpression)
auto<-na.omit(auto)

propx<-subset(mf,locus=="GL343282.1"|locus=="GL343338.1"|locus=="GL343417.1"|locus=="GL343423.1"|locus=="GL343550.1"|locus=="GL343947.1"|locus=="GL343913.1"|locus=="GL343364.1"|locus=="LGb",select=relativeExpression)
propx<-na.omit(propx)

#--------------------------
# 2. Calculate Median Relative Expression and bootstrap p-values
#--------------------------

library("boot")

median(auto$relativeExpression)
ma <- function(x,y) {return(median(sample(auto$relativeExpression, 300, replace=TRUE)))}
mediana <- boot(data=auto$relativeExpression,statistic=ma, R=1000)
cia <- boot.ci(mediana, conf=0.95, type="norm")
cia$norm

median(lga$relativeExpression)
mlga <- function(x,y) {return(median(sample(lga$relativeExpression, 300, replace=TRUE)))}
medianlga <- boot(data=lga$relativeExpression,statistic=mlga, R=1000)
cilga <- boot.ci(medianlga, conf=0.95, type="norm")
cilga$norm

median(lgb$relativeExpression)
mlgb <- function(x,y) {return(median(sample(lgb$relativeExpression, 300, replace=TRUE)))}
medianlgb <- boot(data=lgb$relativeExpression,statistic=mlgb, R=1000)
cilgb <- boot.ci(medianlgb, conf=0.95, type="norm")
cilgb$norm

median(lgc$relativeExpression)
mlgc <- function(x,y) {return(median(sample(lgc$relativeExpression, 300, replace=TRUE)))}
medianlgc <- boot(data=lgc$relativeExpression,statistic=mlgc, R=1000)
cilgc <- boot.ci(medianlgc, conf=0.95, type="norm")
cilgc$norm

median(lgd$relativeExpression)
mlgd <- function(x,y) {return(median(sample(lgd$relativeExpression, 300, replace=TRUE)))}
medianlgd <- boot(data=lgd$relativeExpression,statistic=mlgd, R=1000)
cilgd <- boot.ci(medianlgd, conf=0.95, type="norm")
cilgd$norm

median(lgf$relativeExpression)
mlgf <- function(x,y) {return(median(sample(lgf$relativeExpression, 300, replace=TRUE)))}
medianlgf <- boot(data=lgf$relativeExpression,statistic=mlgf, R=1000)
cilgf <- boot.ci(medianlgf, conf=0.95, type="norm")
cilgf$norm

median(lgg$relativeExpression)
mlgg <- function(x,y) {return(median(sample(lgg$relativeExpression, 300, replace=TRUE)))}
medianlgg <- boot(data=lgg$relativeExpression,statistic=mlgg, R=1000)
cilgg <- boot.ci(medianlgg, conf=0.95, type="norm")
cilgg$norm

median(lgh$relativeExpression)
mlgh <- function(x,y) {return(median(sample(lgh$relativeExpression, 300, replace=TRUE)))}
medianlgh <- boot(data=lgh$relativeExpression,statistic=mlgh, R=1000)
cilgh <- boot.ci(medianlgh, conf=0.95, type="norm")
cilgh$norm

median(micro$relativeExpression)
mmicro <- function(x,y) {return(median(sample(micro$relativeExpression, 300, replace=TRUE)))}
medianmicro <- boot(data=micro$relativeExpression,statistic=mmicro, R=1000)
cimicro <- boot.ci(medianmicro, conf=0.95, type="norm")
cimicro$norm

median(propx$relativeExpression)
mp <- function(x,y) {return(median(sample(propx$relativeExpression, 300, replace=TRUE)))}
medianp <- boot(data=propx$relativeExpression,statistic=mp, R=1000)
cip <- boot.ci(medianp, conf=0.95, type="norm")
cip$norm

#--------------------------
# 3. Calculate Mean Relative Expression and bootstrap p-values
#--------------------------

mean(auto$relativeExpression)
ma <- function(x,y) {return(mean(sample(auto$relativeExpression, 300, replace=TRUE)))}
meana <- boot(data=auto$relativeExpression,statistic=ma, R=1000)
cia <- boot.ci(meana, conf=0.95, type="norm")
cia$norm

mean(lga$relativeExpression)
mlga <- function(x,y) {return(mean(sample(lga$relativeExpression, 300, replace=TRUE)))}
meanlga <- boot(data=lga$relativeExpression,statistic=mlga, R=1000)
cilga <- boot.ci(meanlga, conf=0.95, type="norm")
cilga$norm

mean(lgb$relativeExpression)
mlgb <- function(x,y) {return(mean(sample(lgb$relativeExpression, 300, replace=TRUE)))}
meanlgb <- boot(data=lgb$relativeExpression,statistic=mlgb, R=1000)
cilgb <- boot.ci(meanlgb, conf=0.95, type="norm")
cilgb$norm

mean(lgc$relativeExpression)
mlgc <- function(x,y) {return(mean(sample(lgc$relativeExpression, 300, replace=TRUE)))}
meanlgc <- boot(data=lgc$relativeExpression,statistic=mlgc, R=1000)
cilgc <- boot.ci(meanlgc, conf=0.95, type="norm")
cilgc$norm

mean(lgd$relativeExpression)
mlgd <- function(x,y) {return(mean(sample(lgd$relativeExpression, 300, replace=TRUE)))}
meanlgd <- boot(data=lgd$relativeExpression,statistic=mlgd, R=1000)
cilgd <- boot.ci(meanlgd, conf=0.95, type="norm")
cilgd$norm

mean(lgf$relativeExpression)
mlgf <- function(x,y) {return(mean(sample(lgf$relativeExpression, 300, replace=TRUE)))}
meanlgf <- boot(data=lgf$relativeExpression,statistic=mlgf, R=1000)
cilgf <- boot.ci(meanlgf, conf=0.95, type="norm")
cilgf$norm

mean(lgg$relativeExpression)
mlgg <- function(x,y) {return(mean(sample(lgg$relativeExpression, 300, replace=TRUE)))}
meanlgg <- boot(data=lgg$relativeExpression,statistic=mlgg, R=1000)
cilgg <- boot.ci(meanlgg, conf=0.95, type="norm")
cilgg$norm

mean(lgh$relativeExpression)
mlgh <- function(x,y) {return(mean(sample(lgh$relativeExpression, 300, replace=TRUE)))}
meanlgh <- boot(data=lgh$relativeExpression,statistic=mlgh, R=1000)
cilgh <- boot.ci(meanlgh, conf=0.95, type="norm")
cilgh$norm

mean(micro$relativeExpression)
mmicro <- function(x,y) {return(mean(sample(micro$relativeExpression, 300, replace=TRUE)))}
meanmicro <- boot(data=micro$relativeExpression,statistic=mmicro, R=1000)
cimicro <- boot.ci(meanmicro, conf=0.95, type="norm")
cimicro$norm

mean(propx$relativeExpression)
mp <- function(x,y) {return(mean(sample(propx$relativeExpression, 300, replace=TRUE)))}
meanp <- boot(data=propx$relativeExpression,statistic=mp, R=1000)
cip <- boot.ci(meanp, conf=0.95, type="norm")
cip$norm
