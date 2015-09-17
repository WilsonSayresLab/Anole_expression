##########################################################################################################################
# Overview
#
#
# These are instructins for making histograms of expression levels in male, female, and both male and female green anoles.
# 
#	Required programs		R
##########################################################################################################################

#--------------------
# 1. Import data for male and female comparisons. Each one contains expression values for the both catgory, so importing a third file is unnecessary.
#--------------------

	male<-read.csv("/home/yitzhak/Documents/Squamates/Output_Data/MF_Both/Male_Both_GeneExp.csv",header=TRUE)
	female<-read.csv("/home/yitzhak/Documents/Squamates/Output_Data/MF_Both/Female_Both_GeneExp.csv",header=TRUE)

#--------------------
# 2. Specify desired expression columns. Na.omit is performed in this step to only exlude rows with no expression values for the sample of interest.
#--------------------

	maleexp<-male$value_2
	maleexp<-na.omit(maleexp)

	femaleexp<-female$value_2
	femaleexp<-na.omit(femaleexp)

	bothexp<-male$value_1
	bothexp<-na.omit(bothexp)

#--------------------
# 3. Create a 3 pane image containing a histogram for each comparison
#--------------------

par(mfrow=c(3,1))
	hist(maleexp,breaks=seq(0,3800),xlim=c(0,50),ylim=c(0,8000),main="Male Expression Levels",xlab="Expression Level (FPKM)")
	hist(femaleexp,breaks=seq(0,4500),xlim=c(0,50),ylim=c(0,8000),main="Female Expression Levels",xlab="Expression Level (FPKM)")
	hist(bothexp,breaks=seq(0,3800),xlim=c(0,50),ylim=c(0,8000),main="Combined Male and Female Expression Levels",xlab="Expression Level (FPKM)")
