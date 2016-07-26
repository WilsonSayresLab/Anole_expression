###############################################################
#-------------- OverView ----------------
#
#  Male to female comparative regression plot and values for 
#	anole genome relative expression with FPKM threshold of 2.
#
#	Required Programs:	R
#
#############################################################

#--------------------------
# 1. Read in data for regression analysis
#--------------------------
### Male to Female ###
# Whole genome:
	mf <- read.csv("/home/yitzhak/Documents/AnoleDifExp/RelExp/anoleRelExp_sra_fpkm2.csv",header=TRUE)

# Subset LGb and scaffolds:
	lgb<-subset(mf,locus=="LGb",select=c(value2,value1))
	lgb<-na.omit(lgb)

	scaf<-subset(mf,locus=="GL343282.1"|locus=="GL343338.1"|locus=="GL343417.1"|locus=="GL343423.1"|locus=="GL343550.1"|locus=="GL343947.1"|locus=="GL343913.1"|locus=="GL343364.1",select=c(value2,value1))
	scaf<-na.omit(scaf)

	chrx<-subset(mf,locus=="GL343282.1"|locus=="GL343338.1"|locus=="GL343417.1"|locus=="GL343423.1"|locus=="GL343550.1"|locus=="GL343947.1"|locus=="GL343913.1"|locus=="GL343364.1"|locus=="LGb",select=c(value2,value1))
	chrx<-na.omit(chrx)


### Individual ###
indv<-read.csv("/home/yitzhak/Documents/AnoleDifExp/RelExp/indvRelExp_fpkm2.csv",header=TRUE)

# Subset by comparison:
D15D26<-subset(indv,sample1=="d15"&sample2=="d26",select=c(value1,value2))
D15D47<-subset(indv,sample1=="d15"&sample2=="d47",select=c(value1,value2))
D15D61<-subset(indv,sample1=="d15"&sample2=="d61",select=c(value1,value2))
D26D47<-subset(indv,sample1=="d26"&sample2=="d47",select=c(value1,value2))
D26D61<-subset(indv,sample1=="d26"&sample2=="d61",select=c(value1,value2))
D47D61<-subset(indv,sample1=="d47"&sample2=="d61",select=c(value1,value2))

#--------------------------
# 3. Obtain regression values:
#--------------------------

	### Male to Female ###

lm_mf <- lm(formula=value1~value2, data=mf, na.action=na.omit, singular.ok=TRUE)	
summary(lm_mf)
# r2 = 0.8233 
# p-value: < 2.2e-16

	### Individual ###
lm_d15d26 <- lm(formula=value2~value1, data=D15D26, na.action=na.omit, singular.ok=TRUE)  
summary(lm_d15d26)
# r2 = 0.9528
# p-value: < 2.2e-16

lm_d15d47 <- lm(formula=value2~value1, data=D15D47, na.action=na.omit, singular.ok=TRUE)	
summary(lm_d15d47)
# r2 = 0.6291
# p-value: < 2.2e-16

lm_d15d61 <- lm(formula=value2~value1, data=D15D61, na.action=na.omit, singular.ok=TRUE)	
summary(lm_d15d61)
# r2 = 0.7616
# p-value: < 2.2e-16

lm_d26d47 <- lm(formula=value2~value1, data=D26D47, na.action=na.omit, singular.ok=TRUE)	
summary(lm_d26d47)
# r2 = 0.7124
# p-value: < 2.2e-16

lm_d26d61 <- lm(formula=value2~value1, data=D26D61, na.action=na.omit, singular.ok=TRUE)	
summary(lm_d26d61)
# r2 = 0.802
# p-value: < 2.2e-16


lm_d47d61 <- lm(formula=value2~value1, data=D47D61, na.action=na.omit, singular.ok=TRUE)	
summary(lm_d47d61)
# r2 = 0.7836
# p-value: < 2.2e-16

#--------------------------
# 4. Make scattelplot for male to female expression, including LGb genes in blue and other X-linked genes in green:
#--------------------------

	plot(mf$value1~mf$value2,log="xy",xlab="Male Expression (FPKM)",ylab="Female Expression (FPKM)",
       xaxt="n",yaxt="n",pch=20,col="gray50")
  	axis(1,at=c("2","10","100","1000"),labels=c("2","10","100","1000"))
  	axis(2,at=c("2","10","100","1000"),labels=c("2","10","100","1000"))
		abline(a=0,b=1,col="red")
		points(lgb$value1~lgb$value2,col="blue",pch=18)
		points(scaf$value1~scaf$value2,col="green3",pch=20)
		text(11,1500,labels="R Squared = 0.8233")
    		text(11,700,labels="P < 2.2e-16")


#--------------------------
# 5. Make scattelplot for individual expression:
#--------------------------

par(mfrow=c(3,2))

plot(D15D47$value2~D15D47$value1,log="xy",main=("A. Male 1 to Female 1"),xlab="Male 1 Expression (FPKM)",ylab="Female 1 Expression (FPKM)",
     xaxt="n",yaxt="n",pch=20,col="gray50")
    axis(1,at=c("2","10","100","1000"),labels=c("2","10","100","1000"))
    axis(2,at=c("2","10","100","1000"),labels=c("2","10","100","1000"))
    abline(a=0,b=1,col="red")
    text(13,4000,labels="R Squared = 0.6291")

plot(D15D61$value2~D15D61$value1,log="xy",main=("B. Male 1 to Female 2"),xlab="Male 1 Expression (FPKM)",ylab="Female 2 Expression (FPKM)",
     xaxt="n",yaxt="n",pch=20,col="gray50")
    axis(1,at=c("2","10","100","1000"),labels=c("2","10","100","1000"))
    axis(2,at=c("2","10","100","1000"),labels=c("2","10","100","1000"))
    abline(a=0,b=1,col="red")
    text(13,3000,labels="R Squared = 0.7616")

plot(D26D47$value2~D26D47$value1,log="xy",main=("C. Male 2 to Female 1"),xlab="Male 2 Expression (FPKM)",ylab="Female 1 Expression (FPKM)",
     xaxt="n",yaxt="n",pch=20,col="gray50")
    axis(1,at=c("2","10","100","1000"),labels=c("2","10","100","1000"))
    axis(2,at=c("2","10","100","1000"),labels=c("2","10","100","1000"))
    abline(a=0,b=1,col="red")
    text(13,4000,labels="R Squared = 0.7124")

plot(D26D61$value2~D26D61$value1,log="xy",main=("D. Male 2 to Female 2"),xlab="Male 2 Expression (FPKM)",ylab="Female 2 Expression (FPKM)",
     xaxt="n",yaxt="n",pch=20,col="gray50")
    axis(1,at=c("2","10","100","1000"),labels=c("2","10","100","1000"))
    axis(2,at=c("2","10","100","1000"),labels=c("2","10","100","1000"))
    abline(a=0,b=1,col="red")
    text(13,3000,labels="R Squared = 0.8020")

plot(D15D26$value2~D15D26$value1,log="xy",main=("E. Male 1 to Male 2"),xlab="Male 1 Expression (FPKM)",ylab="Male 1 Expression (FPKM)",
     xaxt="n",yaxt="n",pch=20,col="gray50")
    axis(1,at=c("2","10","100","1000"),labels=c("2","10","100","1000"))
    axis(2,at=c("2","10","100","1000"),labels=c("2","10","100","1000"))
    abline(a=0,b=1,col="red")
    text(13,4000,labels="R Squared = 0.9528")

plot(D47D61$value2~D47D61$value1,log="xy",main=("F. Female 1 to Female 2"),xlab="Female 2 Expression (FPKM)",ylab="Female 2 Expression (FPKM)",
     xaxt="n",yaxt="n",pch=20,col="gray50")
    axis(1,at=c("2","10","100","1000"),labels=c("2","10","100","1000"))
    axis(2,at=c("2","10","100","1000"),labels=c("2","10","100","1000"))
    abline(a=0,b=1,col="red")
    text(13,2500,labels="R Squared = 0.7836")
