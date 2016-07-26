###############################################################
#-------------- OverView ----------------
#
#  Male to female median and CI calculations
#	for anole genome relative expression.
#
#	Required Programs:	R
#				Boot
#############################################################

#--------------------------
# 1. Run TopHat and Cuffdiff on Ocotillo cluster using tophat.sh and cuffdiff.sh scripts, as well as the mf_samples file for cuffdiff.
#--------------------------

#--------------------------
# 2. Download output from Ocotillo and format gene_exp.diff for use in R.
#--------------------------

# join gene_exp.diff file with orthologs list (I saved both as csv files and sorted them ascending on their ASU gene IDs to make sure they were in the same order):

	join --header --check-order -t "," -1 1 -2 2 "/home/yitzhak/Documents/AnoleDifExp/ASU_Acar_v2.2.1_orthologs_201408_v6.csv" "/home/yitzhak/Documents/AnoleDifExp/RelExp/gene_exp.diff" > "/home/yitzhak/Documents/AnoleDifExp/RelExp/anoleRelExp_sra_noMin.csv"

# Delete unecessary columns and remove underscores from header ro (R will replace them with periods anyway).
# Use the text to columns tool to split the locus column into a column for chromosome, locus start, and locus end (I titled them locus, start, and end).
# Determine witch column contains expression results for males (it will be the value# column which corresponds to the sample column which says male)
# Divide the male expression column by the female expression column, copy the column, and paste the data back using "paste only numbers". 
# This will remove the formulas from the cells.
# Lastly, search for "#DIV/0!" and replace with "NA" so na.omit will remove the values in R.

#--------------------------
# 3. Create files for each FPKM threshold (1-4)
#--------------------------

# Manually delete rows with an FPKM lower than the respective cutoff for either male or female expression and save as a new file

#--------------------------
# 4. Import data into R (repeat the following script for each data set separately)
#--------------------------

# No Threshold
mf<-read.csv("/home/yitzhak/Documents/AnoleDifExp/RelExp/anoleRelExp_sra_noMin.csv",header=TRUE)

  # Save copy of all X-linked gene data:
  chrx<-subset(mf,locus=="GL343282.1"|locus=="GL343338.1"|locus=="GL343417.1"|locus=="GL343423.1"|locus=="GL343550.1"|locus=="GL343947.1"|locus=="GL343913.1"|locus=="GL343364.1"|locus=="LGb")
  write.table(chrx,"/home/yitzhak/Documents/AnoleDifExp/RelExp/propX_relExp.csv")

# FPKM 1
mf<-read.csv("/home/yitzhak/Documents/AnoleDifExp/RelExp/anoleRelExp_sra_fpkm1.csv",header=TRUE)

# FPKM 2
mf<-read.csv("/home/yitzhak/Documents/AnoleDifExp/RelExp/anoleRelExp_sra_fpkm2.csv",header=TRUE)

# FPKM 3
mf<-read.csv("/home/yitzhak/Documents/AnoleDifExp/RelExp/anoleRelExp_sra_fpkm3.csv",header=TRUE)

# FPKM 4
mf<-read.csv("/home/yitzhak/Documents/AnoleDifExp/RelExp/anoleRelExp_sra_fpkm4.csv",header=TRUE)

#--------------------------
# 5. Subset by locus
#--------------------------

# Obtain observed genes for whole genome:
rx<-na.omit(mf$relativeExpression)

# Autosomes
auto<-subset(mf,locus=="1"|locus=="2"|locus=="3"|locus=="4"|locus=="5"|locus=="6",select=relativeExpression)
auto<-na.omit(auto)

# X-Linked Genes (NAs are omitted after susetting to preserve genes which may be missing data in other columns)

lgb<-subset(mf,locus=="LGb",select=relativeExpression)
lgb<-na.omit(lgb)

scaf<-subset(mf,locus=="GL343282.1"|locus=="GL343338.1"|locus=="GL343417.1"|locus=="GL343423.1"|locus=="GL343550.1"|locus=="GL343947.1"|locus=="GL343913.1"|locus=="GL343364.1",select=relativeExpression)
scaf<-na.omit(scaf)

propx<-subset(mf,locus=="GL343282.1"|locus=="GL343338.1"|locus=="GL343417.1"|locus=="GL343423.1"|locus=="GL343550.1"|locus=="GL343947.1"|locus=="GL343913.1"|locus=="GL343364.1"|locus=="LGb",select=relativeExpression)
propx<-na.omit(propx)

#--------------------------
# 6. Calculate Median Relative Expression and bootstrap p-values
#--------------------------

library("boot")

# Autosomes
median(auto$relativeExpression)
ma<- function(x,y) {return(median(sample(auto$relativeExpression, 300, replace=TRUE)))}
mediana <- boot(data=auto$relativeExpression,statistic=ma, R=1000)
cia <- boot.ci(mediana, conf=0.95, type="norm")
cia$norm

# X-linked genes
median(lgb$relativeExpression)
ml <- function(x,y) {return(median(sample(lgb$relativeExpression, 300, replace=TRUE)))}
medianl <- boot(data=lgb$relativeExpression,statistic=ml, R=1000)
cil <- boot.ci(medianl, conf=0.95, type="norm")
cil$norm

median(scaf$relativeExpression)
ms <- function(x,y) {return(median(sample(scaf$relativeExpression, 300, replace=TRUE)))}
medians <- boot(data=scaf$relativeExpression,statistic=ms, R=1000)
cis <- boot.ci(medians, conf=0.95, type="norm")
cis$norm

median(propx$relativeExpression)
mp <- function(x,y) {return(median(sample(propx$relativeExpression, 300, replace=TRUE)))}
medianp <- boot(data=propx$relativeExpression,statistic=mp, R=1000)
cip <- boot.ci(medianp, conf=0.95, type="norm")
cip$norm

#--------------------------
# 7. Calculate Mean Relative Expression and bootstrap p-values
#--------------------------

# Autosomes
mean(auto$relativeExpression)
ma<- function(x,y) {return(mean(sample(auto$relativeExpression, 300, replace=TRUE)))}
meana <- boot(data=auto$relativeExpression,statistic=ma, R=1000)
cia <- boot.ci(meana, conf=0.95, type="norm")
cia$norm

# X-linked genes
mean(lgb$relativeExpression)
ml <- function(x,y) {return(mean(sample(lgb$relativeExpression, 300, replace=TRUE)))}
meanl <- boot(data=lgb$relativeExpression,statistic=ml, R=1000)
cil <- boot.ci(meanl, conf=0.95, type="norm")
cil$norm

mean(scaf$relativeExpression)
ms <- function(x,y) {return(mean(sample(scaf$relativeExpression, 300, replace=TRUE)))}
means <- boot(data=scaf$relativeExpression,statistic=ms, R=1000)
cis <- boot.ci(means, conf=0.95, type="norm")
cis$norm

mean(propx$relativeExpression)
mp <- function(x,y) {return(mean(sample(propx$relativeExpression, 300, replace=TRUE)))}
meanp <- boot(data=propx$relativeExpression,statistic=mp, R=1000)
cip <- boot.ci(meanp, conf=0.95, type="norm")
cip$norm
