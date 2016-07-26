##########################################################################################################################
# Overview
#
# These are instructins for making boxplots of relative expression levels in on the proposed anole X chromosome  when a 
#	threshold of 2 FPKM is applied.
# 
#  Required programs		R
##########################################################################################################################

#--------------------
# 1. Import data for X chromosome and subset:
#--------------------

# Import whole gneome data
	mf<-read.csv("/home/yitzhak/Documents/AnoleDifExp/RelExp/anoleRelExp_sra_fpkm2.csv",header=TRUE)

# Subset X-linked genes
	chrx<-subset(mf,locus=="GL343282.1"|locus=="GL343338.1"|locus=="GL343417.1"|locus=="GL343423.1"|locus=="GL343550.1"|locus=="GL343947.1"|locus=="GL343913.1"|locus=="GL343364.1"|locus=="LGb",select=c(locus,relativeExpression))
	chrx<-na.omit(chrx)

	scaf<-chrx$locus
  scaf<-factor(scaf)
	relex<-chrx$relativeExpression
relex<-factor(relex)
	rx<-data.frame(scaf,relex)

#--------------------
# 2. Create boxplot:
#--------------------

par(mar=c(5.1,5.1,4.1,2.1))
	boxplot(rx$relex~rx$scaf,ylab="Relative Expression (M/F)",ylim=c(0,2),xlab="Scaffold",cex.lab=1.5,
        	col=c("green4","green4","green4","green4","green4","green4","green4","green4","blue1"))
		segments(0,1,10,1,col="red2",lwd=2,lty=2)
		segments(0,0.5,10,0.5,lty=2)

#--------------------
# 3. Add number of expressed genes
#--------------------

# Subset by scaffold, record number of genes after subsetting

GL343550.1<-subset(rx,scaf=="GL343550.1")
GL343550.1<-na.omit(GL343550.1)

GL343364.1<-subset(rx,scaf=="GL343364.1")
GL343364.1<-na.omit(GL343364.1)

GL343423.1<-subset(rx,scaf=="GL343423.1")
GL343423.1<-na.omit(GL343423.1)

GL343338.1<-subset(rx,scaf=="GL343338.1")
GL343338.1<-na.omit(GL343338.1)

GL343282.1<-subset(rx,scaf=="GL343282.1")
GL343282.1<-na.omit(GL343282.1)

GL343417.1<-subset(rx,scaf=="GL343417.1")
GL343417.1<-na.omit(GL343417.1)

GL343913.1<-subset(rx,scaf=="GL343913.1")
GL343913.1<-na.omit(GL343913.1)

GL343947.1<-subset(rx,scaf=="GL343947.1")
GL343947.1<-na.omit(GL343947.1)

LGb<-subset(rx,scaf=="LGb")
LGb<-na.omit(LGb)

# Add the following values to the plot in GIMP:

          expressed(total)
GL343282.1 33(87)
GL343338.1 33(49)
GL343364.1 19(21)
GL343417.1 20(42)
GL343423.1 15(18)
GL343550.1 11(11)
GL343913.1 6(8)
GL343947.1 6(12)
LGb 59(80)
