###############################################################
#-------------- OverView ----------------
#
#  Anole genome relative expression histograms
#
#   Male and female relative expression autosomes
#   Male and female relative expression LGb
#   Male and female relative expression hypothetical X
#
#############################################################

#--------------------------
# 1. Set working directory
#--------------------------

setwd("~/Desktop/AnoleChickenSubstitutions/")

#--------------------------
# 2. Run 01_MedianAutoLGbHypoX_rx.R
#--------------------------

# Already calucaluted the relative expression for males and females for:
#   autosomes
#   LGb
#   hypothetical X

#--------------------------
# 3. relative expression males and females for autosomes
#--------------------------

# Males and Females Autosome Relative Expression
# Males yellow Females purple
hist(mrxauto,
main= "Males and Females Relative Expression
Autosomes", #Title of your histogram
xlab="relative expression", # X axis label
ylab="Frequency", # Y axis label
col=rgb(1,1,0,0.7), # color of the rectangles
xlim=c(0,40), # X axis start and end value
ylim=c(0,4500),
las=1,
breaks=2500) # will increase or decrease the number of rectangles displayed on the histogram
par(new=TRUE)
hist(fmrxauto,
main= "Males and Females Relative Expression
Autosomes", #Title of your histogram
xlab="relative expression", # X axis label
ylab="Frequency", # Y axis label
col=rgb(0,0,1,0.2), # color of the rectangles
xlim=c(0,40), # X axis start and end value
ylim=c(0,4500),
las=1,
breaks=2500) # will increase or decrease the number of rectangles displayed on the histogram
legend("topright", c("Males", "Females"), fill=c("yellow", "purple"))

#--------------------------
# 4. relative expression males and females for LGb
#--------------------------

# Males and Females Relative Expression LGb
# Males yellow Females purple
hist(rxLGbmales,
main= "Males and Females Relative Expression
LGb", #Title of your histogram
xlab="relative expression", # X axis label
ylab="Frequency", # Y axis label
col=rgb(1,1,0,0.7), # color of the rectangles
xlim=c(0,200), # X axis start and end value
ylim=c(0,40),
las=1,
breaks=20) # will increase or decrease the number of rectangles displayed on the histogram
par(new=TRUE)
hist(rxLGbfemales,
main = "Males and Females Relative Expression
LGb",
xlab="relative expression", # X axis label
ylab="Frequency", # Y axis label
col=rgb(0,0,1,0.2), # color of the rectangles
xlim=c(0,200), # X axis start and end value
ylim=c(0,40),
las=1,
breaks=20) # will increase or decrease the number of rectangles displayed on the histogram
legend("topright", c("Males", "Females"), fill=c("yellow", "purple"))

#--------------------------
# 5. relative expression males and females for hypothetical X
#--------------------------

# Males and Females hyothetical X Relative Expression
# Males yellow Females purple
hist(putativechrxmales,
main= "
Putative X", #Title of your histogram
xlab="relative expression", # X axis label
ylab="Frequency", # Y axis label
col=rgb(1,1,0,0.7), # color of the rectangles
xlim=c(0,600), # X axis start and end value
ylim=c(0,200),
las=1,
breaks=20) # will increase or decrease the number of rectangles displayed on the histogram
par(new=TRUE)
hist(putativechrxfemales,
main= "Males and Females Relative Expression
Putative X", #Title of your histogram
xlab="relative expression", # X axis label
ylab="Frequency", # Y axis label
col=rgb(0,0,1,0.2), # color of the rectangles
xlim=c(0,600), # X axis start and end value
ylim=c(0,200),
las=1,
breaks=20) # will increase or decrease the number of rectangles displayed on the histogram
legend("topright", c("Males", "Females"), fill=c("yellow", "purple"))

