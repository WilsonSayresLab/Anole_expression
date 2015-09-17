####Plot proposed X vs. LGb in Mb####

#----------------
#1. Open Supplementary_Table1.xlsx in Excel, divide ac_locus_start by 1,000,000 an add following values to the resulting column so that genes will apear in the order they appear on the anole X chromosome:
#----------------

GL343282.1	0
GL343364.1	1.8
LGb		2.9
GL343550.1	6.2
GL343423.1	6.8
GL343913.1	7.7
GL343947.1	7.9
GL343338.1	8.1
GL343417.1	9.4

#-----------------
#2. Read file into R (load gdata):
#-----------------
	#Set working directory to file location or specify full path to file.
	
	chrx<-read.xls("Supplementary_Table1.xlsx",header=TRUE)

	lgb<-subset(chrx,source=="Alfoldi, et al.",select=c(ac_chrx_start,rel.exp))
	lgb<-na.omit(lgb)
	scaffolds<-subset(chrx,source!="Alfoldi, et al.",select=c(ac_chrx_start,rel.exp))
	locuslgb<-as.numeric(as.character(lgb$ac_chrx_start))
	locusscaf<-as.numeric(as.character(scaffolds$ac_chrx_start))

#------------------
#3. Plot scaffold points in green, LGb points as blue diamonds, scaffold median in green, LGb median in blue, and autosomal median as a dotted red line. Plot the sliding window median values in black:
#------------------

plot(locusscaf, scaffolds$rel.exp, pch=20, col='green4', main="Proposed X-Chromosome", xlab="Position on Propsoed Anole X Chromosome (Mb)", ylab="Relative Expression (M/F)", ylim=c(0,2))
	points(locuslgb, lgb$rel.exp, pch=18,col='blue3')
	segments(-2,0.8715495,14,0.8715495,lwd=2,col="green4")
	segments(-2,0.9722797,14,0.9722797,lwd=2,col="blue3")
	segments(-2,1.014611,14,1.014611,lwd=2,lty=3,col="firebrick4")
