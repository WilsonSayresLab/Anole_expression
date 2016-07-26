###############################################################
#-------------- OverView ----------------
#
#  Anole genome Wilcox test comparing autosomal and X-linked 
#	  FPKM for males and females.
#
#	Required Programs:	R
#       boot
#       sm
#
#############################################################

#--------------------------
# 1. Import data into R: 
#--------------------------

mf<-read.csv("/home/yitzhak/Documents/AnoleDifExp/RelExp/anoleRelExp_sra_fpkm2.csv",header=TRUE)
indv <- read.csv("/home/yitzhak/Documents/AnoleDifExp/RelExp/indvRelExp_fpkm2.csv",header=TRUE)

  # Optionally remove 4 genes with expression over 400 FPKM (for whatever reason it wouldn't work with a one line command):
  mf <- subset(mf, ensemblID!="ENSACAG00000003060")
  mf <- subset(mf, ensemblID!="ENSACAG00000003756")
  mf <- subset(mf, ensemblID!="ENSACAG00000007175")
  mf <- subset(mf, ensemblID!="ENSACAG00000025324")

  indv <- subset(indv, ensembl.75_ID!="ENSACAG00000003060")
  indv <- subset(indv, ensembl.75_ID!="ENSACAG00000003756")
  indv <- subset(indv, ensembl.75_ID!="ENSACAG00000007175")
  indv <- subset(indv, ensembl.75_ID!="ENSACAG00000025324")

#--------------------------
# 2. Subset by individual and locus 
#--------------------------

# Autosomes
auto<-subset(mf,locus=="1"|locus=="2"|locus=="3"|locus=="4"|locus=="5"|locus=="6",select=c(value1,value2))
auto<-na.omit(auto)

# X chromosome:
chrx<-subset(mf,locus=="GL343282.1"|locus=="GL343338.1"|locus=="GL343417.1"|locus=="GL343423.1"|
               locus=="GL343550.1"|locus=="GL343947.1"|locus=="GL343913.1"|locus=="GL343364.1"|
               locus=="LGb",select=c(ensemblID,value1,value2))

# Subset into male and female expression before removing na values to retain as many expression values as possible:
  malex<-chrx$value2
    malex<-na.omit(malex)
  femalex<-chrx$value1
    femalex<-na.omit(femalex)

# Individual
# Extract all FPKM values for each individual and remove duplicates to ensure that each unique
# expression value is included.
d15 <- subset(indv,sample1=="d15",select=c(locus,value1))
d15 <- unique(d15)
d26 <- subset(indv,sample1=="d26",select=c(locus,value1))
d26 <- unique(d26)
d47 <- subset(indv,sample2=="d47",select=c(locus,value2))
d47 <- unique(d47)
d61 <- subset(indv,sample2=="d61",select=c(locus,value2))
d61 <- unique(d61)

d15auto <- subset(d15,locus=="1"|locus=="2"|locus=="3"|locus=="4"|locus=="5"|locus=="6",select=value1)
d15x <-  subset(d15,locus=="GL343282.1"|locus=="GL343338.1"|locus=="GL343417.1"|locus=="GL343423.1"|
                  locus=="GL343550.1"|locus=="GL343947.1"|locus=="GL343913.1"|locus=="GL343364.1"|locus=="LGb",
                select=value1)
d26auto <- subset(d26,locus=="1"|locus=="2"|locus=="3"|locus=="4"|locus=="5"|locus=="6",select=value1)
d26x <-  subset(d26,locus=="GL343282.1"|locus=="GL343338.1"|locus=="GL343417.1"|locus=="GL343423.1"|
                  locus=="GL343550.1"|locus=="GL343947.1"|locus=="GL343913.1"|locus=="GL343364.1"|locus=="LGb",
                select=value1)

d47auto <- subset(d47,locus=="1"|locus=="2"|locus=="3"|locus=="4"|locus=="5"|locus=="6",select=value2)
d47x <-  subset(d47,locus=="GL343282.1"|locus=="GL343338.1"|locus=="GL343417.1"|locus=="GL343423.1"|
                  locus=="GL343550.1"|locus=="GL343947.1"|locus=="GL343913.1"|locus=="GL343364.1"|locus=="LGb",
                select=value2)
d61auto <- subset(d61,locus=="1"|locus=="2"|locus=="3"|locus=="4"|locus=="5"|locus=="6",select=value2)
d61x <-  subset(d61,locus=="GL343282.1"|locus=="GL343338.1"|locus=="GL343417.1"|locus=="GL343423.1"|
                  locus=="GL343550.1"|locus=="GL343947.1"|locus=="GL343913.1"|locus=="GL343364.1"|locus=="LGb",
                select=value2)


#--------------------------
# 3. Perform Wilcox tests  
#--------------------------

# Male:

  wilcox.test(auto$value2,malex)

#   Wilcoxon rank sum test with continuity correction
# data:  auto$value2 and malex
# W = 755520, p-value = 0.001134
# alternative hypothesis: true location shift is not equal to 0

# Female:
  
  wilcox.test(auto$value1,femalex)

# Wilcoxon rank sum test with continuity correction
# data:  auto$value1 and femalex
# W = 701300, p-value = 0.2001
# alternative hypothesis: true location shift is not equal to 0

# d15
wilcox.test(d15auto$value1,d15x$value1)

# d26
wilcox.test(d26auto$value1,d26x$value1)

# d47
wilcox.test(d47auto$value2,d47x$value2)

# d61
wilcox.test(d61auto$value2,d61x$value2)

#--------------------------
# 4. Find median FPKM of X-linked and autosmal genes for males and females: 
#--------------------------

library("boot")

# Male:
median(malex)
mdmx <- function(x,y) {return(median(sample(malex, 300, replace=TRUE)))}
medianmx <- boot(data=malex,statistic=mdmx, R=1000)
cimx <- boot.ci(medianmx, conf=0.95, type="norm")
cimx$norm
# 8.60206 (7.151706, 9.669533)
median(auto$value2)
mdma <- function(x,y) {return(median(sample(auto$value2, 300, replace=TRUE)))}
medianma <- boot(data=auto$value2,statistic=mdma, R=1000)
cima <- boot.ci(medianma, conf=0.95, type="norm")
cima$norm
# 11.7743 (8.066553, 11.83179)
median(malex)/median(auto$value2)

# Female:
median(femalex)
mdfx <- function(x,y) {return(median(sample(femalex, 300, replace=TRUE)))}
medianfx <- boot(data=femalex,statistic=mdfx, R=1000)
cifx <- boot.ci(medianfx, conf=0.95, type="norm")
cifx$norm
# 9.88219 (8.646827, 10.99483)
median(auto$value1)
mdfa <- function(x,y) {return(median(sample(auto$value1, 300, replace=TRUE)))}
medianfa <- boot(data=auto$value1,statistic=mdfa, R=1000)
cifa <- boot.ci(medianfa, conf=0.95, type="norm")
cifa$norm
# 12.0345 (11.35353, 15.15441)
median(femalex)/median(auto$value1)

#--------------------------
# 5. Find mean FPKM of X-linked and autosmal genes for males and females: 
#--------------------------

# Male:
mean(malex)
mdmx <- function(x,y) {return(mean(sample(malex, 300, replace=TRUE)))}
meanmx <- boot(data=malex,statistic=mdmx, R=1000)
cimx <- boot.ci(meanmx, conf=0.95, type="norm")
cimx$norm

mean(auto$value2)
mdma <- function(x,y) {return(mean(sample(auto$value2, 300, replace=TRUE)))}
meanma <- boot(data=auto$value2,statistic=mdma, R=1000)
cima <- boot.ci(meanma, conf=0.95, type="norm")
cima$norm

mean(malex)/mean(auto$value2)

# Female:
mean(femalex)
mdfx <- function(x,y) {return(mean(sample(femalex, 300, replace=TRUE)))}
meanfx <- boot(data=femalex,statistic=mdfx, R=1000)
cifx <- boot.ci(meanfx, conf=0.95, type="norm")
cifx$norm

mean(auto$value1)
mdfa <- function(x,y) {return(mean(sample(auto$value1, 300, replace=TRUE)))}
meanfa <- boot(data=auto$value1,statistic=mdfa, R=1000)
cifa <- boot.ci(meanfa, conf=0.95, type="norm")
cifa$norm

mean(femalex)/mean(auto$value1)

#--------------------------
# 6. Plot FPKM of X-linked and autosmal genes for males and females: 
#--------------------------

  plot(auto$value1~auto$value1,log="xy",xlab="Male Expression (FPKM)",ylab="Female Expression (FPKM)",
       xaxt="n",yaxt="n",pch=20,col="gray30")
    axis(1,at=c("0.01","1","100"),labels=c("0.01","1","100"))
    axis(2,at=c("0.01","1","100"),labels=c("0.01","1","100"))
    points(femalex~malex,col="red2",pch=21)
    abline(log10(3.733225),0,lwd=2,lty=3,col="black")
    abline(log10(5.1327),0,lwd=3,lty=2,col="red2")
    abline(v=log10(3.744375),lwd=2,lty=3,col="black")
    abline(v=log10(4.44215),lwd=3,lty=2,col="red2")
    legend("topleft",legend=c("Autosomal Expression","X Chromosome Expression"),col=c("gray30","red2"),pch=c(20,21),
           bty="n")

#--------------------------
# 7. Create FPKM histograms of X-linked and autosmal genes for males and females: 
#--------------------------

par(mfrow=c(1,2))
  
  hist(auto$value2,col="gray40",breaks=seq(0,10800,0.5),xlim=c(2,30),ylim=c(0,350),
     main="Male Autosomal vs. X Chromosome Expression",xlab="FPKM")
  hist(malex,col="blue",add=TRUE,breaks=seq(0,10800,0.5))
    legend("topright",legend=c("Autosomal Expression","X Chromosome Expression"),col=c("gray30","blue"),
           lty=c(2,5),lwd=10,bty="n")
  
  hist(auto$value1,col="gray40",breaks=seq(0,10800,0.5),xlim=c(2,30),ylim=c(0,350),
     main="Female Autosomal vs. X Chromosome Expression",xlab="FPKM")
  hist(femalex,col="blue",add=TRUE,breaks=seq(0,10800,0.5))
    legend("topright",legend=c("Autosomal Expression","X Chromosome Expression"),col=c("gray30","blue"),
           lty=c(1,1),lwd=10,bty="n")


#--------------------------
# 8. Find median FPKM of X-linked and autosmal genes for individuals
#--------------------------

median(d15x$value1)
mdmx <- function(x,y) {return(median(sample(d15x$value1, 300, replace=TRUE)))}
medianmx <- boot(data=d15x$value1,statistic=mdmx, R=1000)
cimx <- boot.ci(medianmx, conf=0.95, type="norm")
cimx$norm

median(d15auto$value1)
mdma <- function(x,y) {return(median(sample(d15auto$value1, 300, replace=TRUE)))}
medianma <- boot(data=d15auto$value1,statistic=mdma, R=1000)
cima <- boot.ci(medianma, conf=0.95, type="norm")
cima$norm

median(d15x$value1)/median(d15auto$value1)

median(d26x$value1)
mdmx <- function(x,y) {return(median(sample(d26x$value1, 300, replace=TRUE)))}
medianmx <- boot(data=d26x$value1,statistic=mdmx, R=1000)
cimx <- boot.ci(medianmx, conf=0.95, type="norm")
cimx$norm

median(d26auto$value1)
mdma <- function(x,y) {return(median(sample(d26auto$value1, 300, replace=TRUE)))}
medianma <- boot(data=d26auto$value1,statistic=mdma, R=1000)
cima <- boot.ci(medianma, conf=0.95, type="norm")
cima$norm

median(d26x$value1)/median(d26auto$value1)

median(d47x$value2)
mdmx <- function(x,y) {return(median(sample(d47x$value2, 300, replace=TRUE)))}
medianmx <- boot(data=d47x$value2,statistic=mdmx, R=1000)
cimx <- boot.ci(medianmx, conf=0.95, type="norm")
cimx$norm

median(d47auto$value2)
mdma <- function(x,y) {return(median(sample(d47auto$value2, 300, replace=TRUE)))}
medianma <- boot(data=d47auto$value2,statistic=mdma, R=1000)
cima <- boot.ci(medianma, conf=0.95, type="norm")
cima$norm

median(d47x$value2)/median(d47auto$value2)

median(d61x$value2)
mdmx <- function(x,y) {return(median(sample(d61x$value2, 300, replace=TRUE)))}
medianmx <- boot(data=d61x$value2,statistic=mdmx, R=1000)
cimx <- boot.ci(medianmx, conf=0.95, type="norm")
cimx$norm

median(d61auto$value2)
mdma <- function(x,y) {return(median(sample(d61auto$value2, 300, replace=TRUE)))}
medianma <- boot(data=d61auto$value2,statistic=mdma, R=1000)
cima <- boot.ci(medianma, conf=0.95, type="norm")
cima$norm

median(d61x$value2)/median(d61auto$value2)

#--------------------------
# 9. Find mean FPKM of X-linked and autosmal genes for individuals
#--------------------------

mean(d15x$value1)
mdmx <- function(x,y) {return(mean(sample(d15x$value1, 300, replace=TRUE)))}
meanmx <- boot(data=d15x$value1,statistic=mdmx, R=1000)
cimx <- boot.ci(meanmx, conf=0.95, type="norm")
cimx$norm

mean(d15auto$value1)
mdma <- function(x,y) {return(mean(sample(d15auto$value1, 300, replace=TRUE)))}
meanma <- boot(data=d15auto$value1,statistic=mdma, R=1000)
cima <- boot.ci(meanma, conf=0.95, type="norm")
cima$norm

mean(d15x$value1)/mean(d15auto$value1)

mean(d26x$value1)
mdmx <- function(x,y) {return(mean(sample(d26x$value1, 300, replace=TRUE)))}
meanmx <- boot(data=d26x$value1,statistic=mdmx, R=1000)
cimx <- boot.ci(meanmx, conf=0.95, type="norm")
cimx$norm

mean(d26auto$value1)
mdma <- function(x,y) {return(mean(sample(d26auto$value1, 300, replace=TRUE)))}
meanma <- boot(data=d26auto$value1,statistic=mdma, R=1000)
cima <- boot.ci(meanma, conf=0.95, type="norm")
cima$norm

mean(d26x$value1)/mean(d26auto$value1)

mean(d47x$value2)
mdmx <- function(x,y) {return(mean(sample(d47x$value2, 300, replace=TRUE)))}
meanmx <- boot(data=d47x$value2,statistic=mdmx, R=1000)
cimx <- boot.ci(meanmx, conf=0.95, type="norm")
cimx$norm

mean(d47auto$value2)
mdma <- function(x,y) {return(mean(sample(d47auto$value2, 300, replace=TRUE)))}
meanma <- boot(data=d47auto$value2,statistic=mdma, R=1000)
cima <- boot.ci(meanma, conf=0.95, type="norm")
cima$norm

mean(d47x$value2)/mean(d47auto$value2)

mean(d61x$value2)
mdmx <- function(x,y) {return(mean(sample(d61x$value2, 300, replace=TRUE)))}
meanmx <- boot(data=d61x$value2,statistic=mdmx, R=1000)
cimx <- boot.ci(meanmx, conf=0.95, type="norm")
cimx$norm

mean(d61auto$value2)
mdma <- function(x,y) {return(mean(sample(d61auto$value2, 300, replace=TRUE)))}
meanma <- boot(data=d61auto$value2,statistic=mdma, R=1000)
cima <- boot.ci(meanma, conf=0.95, type="norm")
cima$norm

mean(d61x$value2)/mean(d61auto$value2)

#--------------------------
# 10. Create density plots for of autosomal and X chromosome expression for each comparison
#--------------------------

  library("sm")
# First combine into a single data frame for each comparison

automale <- data.frame(auto$value2, "autosome")
colnames(automale) <- c("expression", "locus")
xmale <- data.frame(malex, "X-linked")
colnames(xmale) <- c("expression", "locus")
male <- rbind(automale,xmale)

autofemale <- data.frame(auto$value1, "autosome")
colnames(autofemale) <- c("expression", "locus")
xfemale <- data.frame(femalex, "X-linked")
colnames(xfemale) <- c("expression", "locus")
female <- rbind(autofemale,xfemale)

autod15 <- data.frame(d15auto$value1, "autosome")
colnames(autod15) <- c("expression", "locus")
xd15 <- data.frame(d15x, "X-linked")
colnames(xd15) <- c("expression", "locus")
d15 <- rbind(autod15,xd15)

autod26 <- data.frame(d26auto$value2, "autosome")
colnames(autod26) <- c("expression", "locus")
xd26 <- data.frame(d26x, "X-linked")
colnames(xd26) <- c("expression", "locus")
d26 <- rbind(autod26,xd26)

autod47 <- data.frame(d47auto$value1, "autosome")
colnames(autod47) <- c("expression", "locus")
xd47 <- data.frame(d47x, "X-linked")
colnames(xd47) <- c("expression", "locus")
d47 <- rbind(autod47,xd47)

autod61 <- data.frame(d61auto$value2, "autosome")
colnames(autod61) <- c("expression", "locus")
xd61 <- data.frame(d61x, "X-linked")
colnames(xd61) <- c("expression", "locus")
d61 <- rbind(autod61,xd61)

# Create a density plot for each comparison in one pane

par(mfrow=c(3,2))

sm.density.compare(male$expression, male$locus, col=c("black","blue"),lty=c(1,1),xlim=c(2,150),xlab="FPKM")
  title(main="A. Male Autosomal vs. X Chromosome Expression Density")
  segments(11.8,0,11.8,0.03, lty=3)
  segments(8.6,0,8.6,0.03, lty=3, col="blue")
  segments(24.1,0,24.1,0.03, lty=2)
  segments(27.2,0,27.4,0.03, lty=2, col="blue")
  axis(1,at=c("2","10","25","50","100","150"),labels=c("2","10","25","50","100","150"))
  legend(60,0.021, legend=c("Autosomal Expression","X Chromosome Expression","Autosomal Mean Expression",
        "X Chromosome Mean Expression","Autosomal Median Expression","X Chromosome Median Expression"),
         col=c("black","blue","black","blue","black","blue"),lty=c(1,1,2,2,3,3),lwd=1,bty="n")

sm.density.compare(female$expression, female$locus, col=c("black","blue"),lty=c(1,1),xlim=c(2,150),xlab="FPKM")
  title(main="B. Female Autosomal vs. X Chromosome Expression Density")
  segments(12.,0,12,0.03, lty=3)
  segments(9.9,0,9.0,0.03, lty=3, col="blue")
  segments(25.4,0,25.4,0.03, lty=2)
  segments(32.7,0,32.7,0.03, lty=2, col="blue")
  axis(1,at=c("2","10","25","50","100","150"),labels=c("2","10","25","50","100","150"))
  legend(60,0.018,legend=c("Autosomal Expression","X Chromosome Expression","Autosomal Mean Expression",
                           "X Chromosome Mean Expression","Autosomal Median Expression","X Chromosome Median Expression"),
       col=c("black","blue","black","blue","black","blue"),lty=c(1,1,2,2,3,3),lwd=1,bty="n")

sm.density.compare(d15$expression, d15$locus, col=c("black","blue"),lty=c(1,1),xlim=c(2,150),xlab="FPKM")
title(main="C. Male 1 Autosomal vs. X Chromosome Expression Density")
  segments(11.6,0,11.6,0.03, lty=3)
      segments(9.1,0,9.1,0.03, lty=3, col="blue")
      segments(24.6,0,24.6,0.03, lty=2)
      segments(29.5,0,29.5,0.03, lty=2, col="blue")
      axis(1,at=c("2","10","25","50","100","150"),labels=c("2","10","25","50","100","150"))
      legend(60,0.019, legend=c("Autosomal Expression","X Chromosome Expression","Autosomal Mean Expression",
      "X Chromosome Mean Expression","Autosomal Median Expression","X Chromosome Median Expression"),
      col=c("black","blue","black","blue","black","blue"),lty=c(1,1,2,2,3,3),lwd=1,bty="n")
      
sm.density.compare(d26$expression, d26$locus, col=c("black","blue"),lty=c(1,1),xlim=c(2,150),xlab="FPKM")
  title(main="D. Male 2 Autosomal vs. X Chromosome Expression Density")
  segments(11.9,0,11.9,0.03, lty=3)
  segments(8.8,0,8.8,0.03, lty=3, col="blue")
  segments(24.6,0,24.6,0.03, lty=2)
  segments(25.9,0,25.9,0.03, lty=2, col="blue")
  axis(1,at=c("2","10","25","50","100","150"),labels=c("2","10","25","50","100","150"))
  legend(60,0.021, legend=c("Autosomal Expression","X Chromosome Expression","Autosomal Mean Expression",
                           "X Chromosome Mean Expression","Autosomal Median Expression","X Chromosome Median Expression"),
       col=c("black","blue","black","blue","black","blue"),lty=c(1,1,2,2,3,3),lwd=1,bty="n")

sm.density.compare(d47$expression, d47$locus, col=c("black","blue"),lty=c(1,1),xlim=c(2,150),xlab="FPKM")
  title(main="E. Female 1 Autosomal vs. X Chromosome Expression Density")
  segments(12.1,0,12.1,0.03, lty=3)
  segments(10.9,0,10.9,0.03, lty=3, col="blue")
  segments(26.9,0,26.9,0.03, lty=2)
  segments(38.7,0,38.7,0.03, lty=2, col="blue")
  axis(1,at=c("2","10","25","50","100","150"),labels=c("2","10","25","50","100","150"))
  legend(60,0.0145, legend=c("Autosomal Expression","X Chromosome Expression","Autosomal Mean Expression",
                           "X Chromosome Mean Expression","Autosomal Median Expression","X Chromosome Median Expression"),
       col=c("black","blue","black","blue","black","blue"),lty=c(1,1,2,2,3,3),lwd=1,bty="n")

sm.density.compare(d61$expression, d61$locus, col=c("black","blue"),lty=c(1,1),xlim=c(2,150),xlab="FPKM")
  title(main="F. Female 2 Autosomal vs. X Chromosome Expression Density")
  segments(12.4,0,12.4,0.03, lty=3)
  segments(9.7,0,9.7,0.03, lty=3, col="blue")
  segments(25.5,0,25.5,0.03, lty=2)
  segments(27.8,0,27.8,0.03, lty=2, col="blue")
  axis(1,at=c("2","10","25","50","100","150"),labels=c("2","10","25","50","100","150"))
  legend(60,0.021, legend=c("Autosomal Expression","X Chromosome Expression","Autosomal Mean Expression",
                           "X Chromosome Mean Expression","Autosomal Median Expression","X Chromosome Median Expression"),
       col=c("black","blue","black","blue","black","blue"),lty=c(1,1,2,2,3,3),lwd=1,bty="n")

# Create a density plot for males and females in one pane

svg("/home/yitzhak/Dropbox/Squamates/Rupp_lab_notes/Plots/MaleFemaleDensity.svg")

par(mfrow=c(2,1))

sm.density.compare(male$expression, male$locus, col=c("black","blue"),lty=c(1,1),xlim=c(2,150),xlab="FPKM")
title(main="A. Male Autosomal vs. X Chromosome Expression Density")
segments(11.8,0,11.8,0.03, lty=3)
segments(8.6,0,8.6,0.03, lty=3, col="blue")
segments(24.1,0,24.1,0.03, lty=2)
segments(27.2,0,27.4,0.03, lty=2, col="blue")
axis(1,at=c("10","25","50","100","150"),labels=c("10","25","50","100","150"))
legend(60,0.021, legend=c("Autosomal Expression","X Chromosome Expression","Autosomal Mean Expression",
                          "X Chromosome Mean Expression","Autosomal Median Expression","X Chromosome Median Expression"),
       col=c("black","blue","black","blue","black","blue"),lty=c(1,1,2,2,3,3),lwd=1,bty="n")

sm.density.compare(female$expression, female$locus, col=c("black","blue"),lty=c(1,1),xlim=c(2,150),xlab="FPKM")
title(main="B. Female Autosomal vs. X Chromosome Expression Density")
segments(12.,0,12,0.03, lty=3)
segments(9.9,0,9.0,0.03, lty=3, col="blue")
segments(25.4,0,25.4,0.03, lty=2)
segments(32.7,0,32.7,0.03, lty=2, col="blue")
axis(1,at=c("10","25","50","100","150"),labels=c("10","25","50","100","150"))
legend(60,0.018,legend=c("Autosomal Expression","X Chromosome Expression","Autosomal Mean Expression",
                         "X Chromosome Mean Expression","Autosomal Median Expression","X Chromosome Median Expression"),
       col=c("black","blue","black","blue","black","blue"),lty=c(1,1,2,2,3,3),lwd=1,bty="n")

dev.off()
