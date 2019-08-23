getwd()
setwd("D:/Introduction to R")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()
BiocManager::install("Biobase")
BiocManager::install("genefilter")
BiocManager::install("affy")
install.packages("devtools")
install.packages("RDocumentation")
devtools::install_github("datacamp/RDocumentation")
#BiocManager::install("rSFFreader")
#BiocManager::install("rSFFreader", repo='https://git.bioconductor.org/packages/rSFFreader', dependencies = TRUE, ref = 'master')
BiocManager::install("oligo")
BiocManager::install("GEOquery")
install.packages("normtest")
install.packages("reshape2")

library(RDocumentation)
library(dint)
library(Biobase)
library(genefilter)
library(affy)
#library("rSFFreader")
library(oligo)
library(GEOquery)
library(normtest)
library(limma)
library(reshape2)
library(ggplot2)

variables<-ls()
rm(variables)

inc <- function(x){
return(x<- x + 1)
}



#0. obtain and read in the sequence data
#sequence <- readSff("sequencedata","SRR048068.sff",package="rSFFreader")
getGEOSuppFiles("GSM206007")
list.files("GSM206607")
untar("D:/Introduction to R/GSM206607.gz", exdir ="?D:/Introduction to R/GSM206607")
#celfile<- list.files(celfiles <- list.files("GSM", full = TRUE))
untar("D:/Introduction to R/GSM206608.gz", exdir ="?D:/Introduction to R/GSM206608")
filename<- c("GSM206607", "GSM206608", "GSM206609", "GSM206610", "GSM206611", "GSM206612", "GSM206613", "GSM206614", "GSM206615", "GSM206616")
sourcename<- c("GSM206607.CEL", "GSM206608.CEL", "GSM206609.CEL", "GSM206610.CEL", "GSM206611.CEL", "GSM206612.CEL", "GSM206613.CEL", "GSM206614.CEL", "GSM206615.CEL", "GSM206616.CEL")
#GSM<- list("Geneseq" = NULL)
GSM<- data.frame(Gene = seq(1, length(exprs(read.celfiles(sourcename[1])))))
GSM<- sapply(seq(1,10), function(i) GSM[filename[i]]<- exprs(read.celfiles(sourcename[i])))
colnames(GSM)<- filename
write.csv(GSM, file = "GSM.csv", col.names = T)
dat<-ls()
save(dat, file = "GSM.RData")
GSMstat<- data.frame("stats" = names(summary(0)))
GSMstat<-sapply(seq(1,10), function(i)  GSMstat[Samples[i]]<-(as.vector(summary(GSM_df[Genes[i]]))))
colnames(GSMstat)<- Samples
write.csv(GSMstat, file = "GSMstat.csv", col.names = T)

#read in the binded GSM dataframe
GSM_df<- read.csv("GSM.csv")
Samples<- names(GSM_df)
GSMstat<-  read.csv("GSMstat.csv")


#1. Overview the raw data and preprocess 

#boxplot(GSM_df[2:11], main = "rawData Gene Expressions")
rawData<- as.matrix(GSM_df[2:11])
#GSMstat["TotalRaw"]<- as.vector(summary(rawData))
boxplot(rawData, main = "rawData Gene Expressions")
GSMstat["TotalRaw"]<- as.vector(summary(as.numeric(rawData)))
                                                                                                  + +
#normalization:  
#Z-transform(lightly skewed data): normData<- (rawData - mean(rawData))/sd(rawData)
d<- density(rawData)  
hist(rawData, freq = F)
lines(d)
#normData<- rma(GSM_df[2:11])
normData<- rawData
normData<- sapply(seq(1,10), function(i) normData[,i]<- (rawData[,i] - mean(rawData[,i]))/sd(rawData[,i]))
GSMstat["TotalNorm"]<- as.vector(summary(as.numeric(normData)))
d<- density(normData)  
hist(normData, freq = F)
lines(d)   
boxplot(normData, main = "normData Gene Expressions")  
kurtosis.norm.test(normData)
######################

        Kurtosis test for normality

data:  normData
T = 117.67, p-value < 2.2e-16
##################################################################
#Log-Normalization(large-scale skewed data) and MAnormalization
#Microarray data is often normalized within arrays to control for 
#systematic biases in dye coupling and hybridization efficiencies, 
#as well as other technical biases in the DNA probes 
#and the print tip used to spot the array.
#By minimizing these systematic variations, 
#true biological differences can be found. 
M=log _{2}(R/G)=log _{2}(R)-\log _{2}(G)} M=\log _{2}(R/G)=\log _{2}(R)-\log _{2}(G)
displaystyle A={\frac {1}{2}}\log _{2}(RG)={\frac {1}{2}}(\log _{2}(R)+\log _{2}(G))} A={\frac  12}\log _{2}(RG)={\frac  12}(\log _{2}(R)+\log _{2}(G))
M is, therefore, the binary logarithm of the intensity ratio (or difference between log intensities) and A is the average log intensity for a dot in the plot. 
MA plots are then used to visualize intensity-dependent ratio of raw microarray data (microarrays typically show a bias here, with higher A resulting in higher |M|,
i.e. the brighter the spot the more likely an observed difference between sample and control).
The MA plot uses M as the y-axis and A as the x-axis and gives a quick overview of the distribution of the data.
In many microarray gene expression experiments, an underlying assumption is that most of the genes would not see any change in their expression; therefore, the majority of the points on the y-axis (M) would be located at 0, since Log(1) is 0. 
If this is not the case, then a normalization method such as LOESS should be applied to the data before statistical analysis. (On the diagram below see the red line running below the zero mark before normalization, it should be straight.
Since it is not straight, the data should be normalized. After being normalized, the red line is straight on the zero line and shows as pink/black.)
##################################################################

logData<- log2(rawData)
boxplot(logData)
GSMstat["TotalLog"]<- summary(as.numeric(logData))
d<- density(logData)  
hist(logData, freq = F)
lines(d) 
kurtosis.norm.test(logData)
####################################
        Kurtosis test for normality


data:  logData
T = 5.2636, p-value < 2.2e-16
#########################################
M<- logData
diff<- logData
MAData<- logData
A<- logData
avg<- logData
MAD<- logData[1,]
scaleData<- logData
#look into differences between genes(supposed to be large)
logData<- rbind(logData,logData[1,])
diff<- sapply(seq(1,length(diff)/10), function(i) diff[i,]<- (logData[i+1,]-logData[i,]))
avg<- sapply(seq(1,length(diff)/10), function(i) avg[i,]<- 0.5*(logData[i+1,]-logData[i,]))
smoothScatter(avg, diff, col = 1, main="MA plot for 10 samples", xlab="A", ylab="M")
abline(h=c(-1,1), col="red")
logData<- log2(rawData)

#look into differences between genes(supposed to be small,close to zero)
logData<- rbind(logData,logData[1,])
M<- sapply(seq(1,10), function(i) M[i,]<- (logData[i+1,]-logData[i,]))
A<- sapply(seq(1,10), function(i) A[i,]<- 0.5*(logData[i+1,]+logData[i,]))
#MA Plot of random unrepeated sampled 100 genes(from all length(diff))
seed(123)
Ind<- sample(length(diff)/10, 100)
smoothScatter(avg[Ind,], diff[Ind,], col = 1, main="MA plot for 100 random genes", xlab="A", ylab="M")
abline(h=c(-1,1), col="red")
# according to the result, there is no bias of M thus does not need to do the LOESS normalization
# there are some significantly abnormal points thus need to be LOESS normalized
#LOESS normalization
sapply(seq(1,length(diff)/10), function(i) M[,i]<- M[,i])
logData<- log2(rawData)
l = loess(diff~avg[])
mc<- p

diff<- sapply(seq(1,length(diff)/10), function(i) diff[,i]<- (logData[,i+1]-logData[,i]))
MAData<- sapply(seq(1,10), function(i) MAData[,i]<- diff[,i]/avg[i])
GSMstat["TotalMA"]<- summary(as.numeric(MAData))
d<- density(MAData)  
hist(MAData, freq = F)
lines(d) 
kurtosis.norm.test(MAData)


# Scale-normalization(if data is extremely skewed) 
MAD<- sapply(seq(1,10), function(i) MAD[i]<- meadian(abs(sd(logData[i])))
scaleData<- sapply(seq(1,10), function(i) scaleData[,i]<- (scaleData[,i] - meadian(logData[,i]))/MAD[i])
# look at one gene: 19th, 
m<- diff[19,]
g<- logData[19,]

layout(matrix(c(1,2), 1, 2, byrow = TRUE))
plot(g, col=2, type = 'o',pch = 19, xaxt="n", xlab="datatype", ylab="log2(Expression)", main=paste("log-ratio:", round(m,3)[1]))
axis(1, labels=c("Log-norm"), at=1)
plot(m, pch = 19, col = 5, type = 'o')
axis(2, labels=c("Log-ratio"), at=2)

# heatmap of first 100 genes(10 samples)
png(file = "heatmap(norm).png", bg = "transparent")
heatmap(normData[1:100,])
dev.off()


# load GPL Annotation file
gpl2872 <- getGEO(filename='GPL2872_old_annotations.txt')
#for the purposes of exploring in gene and protein, we are only interested in 

subset <- gpl2872[, gpl2872$GO_ID %like% "GO:0005-"]
heatmap(exprs(subset[1:100,]))
f<- factor(as.character(subset$GO_ID))
design<- model.matrix(~f)
fit<- eBayes(lmFit(subset,design))
topTable(fit, coef=length(f))
selected  <- p.adjust(fit$p.value[,2]) <0.05
subsetSel <- subset [selected, ]
heatmap(exprs(subsetSel))


