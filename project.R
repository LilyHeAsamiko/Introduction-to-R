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
install.packages("rnaseqWrapper")
install.packages("seqinr")
install.packages("ape")
BiocManager::install(c("DESeq","topGO"))
install.packages("Peptides") 



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
library(seqinr)
library(ape) 
library(stringr)
library(Peptides) 

variables<-ls()
rm(variables)

inc <- function(x){
return(x<- x + 1)
}

"

#0. obtain and read in the sequence data
#sequence <- readSff("sequencedata","SRR048068.sff",package="rSFFreader")
getGEOSuppFiles("GSM206007")
list.files("GSM206607")
untar("D:/Introduction to R/GSM206607.gz", exdir ="?D:/Introduction to R/GSM206607")
#celfile<- list.files(celfiles <- list.files("GSM", full = TRUE))
untar("D:/Introduction to R/GSM206608.gz", exdir ="D:/Introduction to R/GSM206608")
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

#look into differences between genes(supposed to be large)
logData<- rbind(logData,logData[1,])
diff<- sapply(seq(1,length(diff)/10), function(i) diff[i,]<- (logData[i+1,]-logData[i,]))
avg<- sapply(seq(1,length(diff)/10), function(i) avg[i,]<- 0.5*(logData[i+1,]-logData[i,]))
smoothScatter(avg, diff, col = 1, main="MA plot for 10 samples", xlab="A", ylab="M")
abline(h=c(-1,1), col="red")
logData<- log2(rawData)

#look into differences between samples(supposed to be small,close to zero)
logData<- cbind(logData,logData[,1])
M<- sapply(seq(1,10), function(i) M[,i]<- (logData[,i+1]-logData[,i]))
A<- sapply(seq(1,10), function(i) A[,i]<- 0.5*(logData[,i+1]+logData[,i]))
#MA Plot of random unrepeated sampled 100 genes(from all length(diff))
set.seed(123)
Ind<- sample(seq(1,length(M)/10), 100)
smoothScatter(avg[Ind,], diff[Ind,], col = 1, main="MA plot for 100 random genes", xlab="A", ylab="M")
abline(h=c(-1,1), col="red")
# according to the result, there is no bias of M thus does not need to do the LOESS normalization
# there are some significantly abnormal points thus need to be LOESS normalized
#LOESS normalization
#MAData<- sapply(seq(1,10), function(i) M[Ind,i]<- (M[Ind,i] - loess(M[Ind,i]~A[Ind,i])))
Mc<- A[,1]
MAData<- sapply(seq(1,10), function(i) {l<- loess(M[Ind,i]~A[Ind,i])
Mc<- predict(l, A[,i])
Mc[is.na(Mc)]<- 0
MAData[,i]<- (M[,i] - Mc)})
boxplot(MAData, main = "MA_Normed Gene Expressions")
logData<- log2(rawData)
GSMstat["TotalMA"]<- summary(as.numeric(MAData))
write.csv(GSMstat, file = "GSMstat.csv", col.names = T)
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

#look into exact genes: 19 
id = 19
layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
plot(rawData[id,], type = 'b',pch = 19, xlab="Samples", ylab="rawData(Expression)", main=paste("rawData of Gene:", id))
plot(normData[id,], type = 'b',pch = 19, xlab="Samples", ylab="normData(Expression)", main=paste("normData of Gene:", id))
plot(logData[id,], type = 'b',pch = 19, xlab="Samples", ylab="logData(Expression)", main=paste("logData of Gene:", id))
plot(MAData[id,], type = 'b',pch = 19, xlab="Samples", ylab="MAData(Expression)", main=paste("MAData of Gene:", id))

# heatmap of first 100 genes(10 samples)
png(file = "heatmap(norm).png", bg = "transparent")
heatmap(normData[1:100,])
dev.off()

#gpl

#ID = Agilent feature number
#COL = Column
#ROW = Row
#NAME = Name
#SPOT_ID = Spot identifier
#CONTROL_TYPE = Control type
#REFSEQ = RefSeqAccession
#GB_ACC = GenBankAccession LINK_PRE:"http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Search&db=Nucleotide&term="
#GENE = Entrez Gene ID LINK_PRE:"http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=gene&cmd=Retrieve&dopt=Graphics&list_uids="
#GENE_SYMBOL = Gene Symbol
#GENE_NAME = Gene Name
#UNIGENE_ID = UnigeneID
#ENSEMBL_ID = EnsemblID
#TIGR_ID = TIGRID
#ACCESSION_STRING = Accession String
#CHROMOSOMAL_LOCATION = Chromosomal Location
#CYTOBAND = Cytoband
#DESCRIPTION = Description
#GO_ID = GoIDs
#SEQUENCE = 

#prepare data
gpl2872<- read.csv('GPL2872_old_annotations.csv', skip = 22, fill = TRUE, header = FALSE)
gpl2872matrix <- as.matrix(read.csv('GPL2872_old_annotations.csv', skip = 21, fill = TRUE, header = FALSE))
colnames(gpl2872)<- gpl2872matrix[1,]
#rownames(gpl2872)<- t(gpl2872[1])
write.csv(gpl2872, file = 'gpl2872.csv', col.names= TRUE)

#read in data
gpl2872<- read.csv('gpl2872.csv', header= TRUE, fill = TRUE)
gpl2872Var<- colnames(gpl2872) 

#gpl2872matrix<- sapply(seq(1,nrow(gpl2872)), function(i) {cbind(gpl2872matri
#x[i,ncol(gpl2872)], grep('A', gpl2872matrix[i,ncol(gpl2872matrix)]))
#cbind(gpl2872matrix[i,ncol(gpl2872matrix)],length(grep('A', gpl2872matrix[i,ncol(gpl2872matrix)])))})
posCount<- function(ACTG, sequence) {
A<- as.vector(gregexpr(ACTG, sequence)[[1]])
LenA<- length(A)
if ( A == -1) 
{return(c(0,0))}
else {
return(c(A,LenA))}
}

posA<- sapply(seq(1,nrow(gpl2872)), function(i) {posCount('A', gpl2872[i,ncol(gpl2872)])[-1]}) 
countA<- sapply(seq(1,nrow(gpl2872)), function(i) {posCount('A', gpl2872[i,ncol(gpl2872)])[length(posCount('A', gpl2872[i,ncol(gpl2872)]))]}) 
posC<- sapply(seq(1,nrow(gpl2872)), function(i) {posCount('C', gpl2872[i,ncol(gpl2872)])[-1]}) 
countC<- sapply(seq(1,nrow(gpl2872)), function(i) {posCount('C', gpl2872[i,ncol(gpl2872)])[length(posCount('C', gpl2872[i,ncol(gpl2872)]))]}) 
posT<- sapply(seq(1,nrow(gpl2872)), function(i) {posCount('T', gpl2872[i,ncol(gpl2872)])[-1]}) 
countT<- sapply(seq(1,nrow(gpl2872)), function(i) {posCount('T', gpl2872[i,ncol(gpl2872)])[length(posCount('T', gpl2872[i,ncol(gpl2872)]))]}) 
posG<- sapply(seq(1,nrow(gpl2872)), function(i) {posCount('G', gpl2872[i,ncol(gpl2872)])[-1]}) 
countG<- sapply(seq(1,nrow(gpl2872)), function(i) {posCount('G', gpl2872[i,ncol(gpl2872)])[length(posCount('G', gpl2872[i,ncol(gpl2872)]))]}) 


#return(reshape(Adf, direction="wide", sep=""))}}) 
#gpl2872["PosA"]<- as.data.frame(posA)
LengthA<- max(countA)
LengthC<- max(countC)
LengthT<- max(countT)
LengthG<- max(countG)
fillLine<- function(L, sequence){
if (length(sequence)< L) {sequence<- c(sequence, rep(NA, L-length(sequence)))}
return(sequence)
}
posATb<- t(sapply(seq(1,length(posA)), function(i) {fillLine(LengthA, as.vector(posA[i][[1]]))})) 
posCTb<- t(sapply(seq(1,length(posC)), function(i) {fillLine(LengthC, as.vector(posC[i][[1]]))})) 
posTTb<- t(sapply(seq(1,length(posT)), function(i) {fillLine(LengthT, as.vector(posT[i][[1]]))})) 
posGTb<- t(sapply(seq(1,length(posG)), function(i) {fillLine(LengthG, as.vector(posG[i][[1]]))})) 

posCountTb<- functon(df){
posA<- sapply(seq(1,nrow(df)), function(i) {posCount('A', df[i,ncol(df)])[-1]}) 
countA<- sapply(seq(1,nrow(df)), function(i) {pCA<- posCount('A', df[i,ncol(df)])
return(pCA[length(pCA)))]})
posC<- sapply(seq(1,nrow(df)), function(i) {posCount('C', df[i,ncol(df)])[-1]}) 
countC<- sapply(seq(1,nrow(df)), function(i) {pCC<- posCount('C', df[i,ncol(df)])
return(pCC[length(pCC)))]})return(pCC[length(pCC)))]})
posT<- sapply(seq(1,nrow(df)), function(i) {posCount('T', df[i,ncol(df)])[-1]}) 
countT<- sapply(seq(1,nrow(df)), function(i) {pCT<- posCount('T', df[i,ncol(df)])
return(pCT[length(pCT)))]})
posG<- sapply(seq(1,nrow(df)), function(i) {posCount('G', df[i,ncol(df)])[-1]}) 
countG<- sapply(seq(1,nrow(df)), function(i) {pCG<- posCount('G', df[i,ncol(df)])
return(pCG[length(pCG)))]})
posATb<- t(sapply(seq(1,length(posA)), function(i) {fillLine(LengthA, as.vector(posA[i][[1]]))})) 
posCTb<- t(sapply(seq(1,length(posC)), function(i) {fillLine(LengthC, as.vector(posC[i][[1]]))})) 
posTTb<- t(sapply(seq(1,length(posT)), function(i) {fillLine(LengthT, as.vector(posT[i][[1]]))})) 
posGTb<- t(sapply(seq(1,length(posG)), function(i) {fillLine(LengthG, as.vector(posG[i][[1]]))})) 
df["PosA"]<- posATb
df["CountA"]<- countA
df["PosC"]<- posCTb
df["CountC"]<- countC
df["PosT"]<- posTTb
df["CountT"]<- countT
df["PosG"]<- posGTb
df["CountG"]<- countG
return(df)
}

gpl2872["PosA"]<- posATb
gpl2872["CountA"]<- countA
gpl2872["PosC"]<- posCTb
gpl2872["CountC"]<- countC
gpl2872["PosT"]<- posTTb
gpl2872["CountT"]<- countT
gpl2872["PosG"]<- posGTb
gpl2872["CountG"]<- countG

gpl2872<- posCountTb(gpl2872)

#The GC content tells us what the ratio of G (guanines) and C (cystosines) residues compared to A
#(adenine) and T (thymidine) residues in a sequence is
#This is important because coding regions tend to be higher in GC.
#GC content also affects the "melting" temperature of DNA.
#GCcontent<- (CountG + CountC)/(CountA + CountT)

GC<- funcion(sequence){
pCA<- posCount('A', sequence)
posA<- pCA[-1]
countA<-pCA[Length[pCA]]
pCC<- posCount('C', sequence)
posC<- pCC[-1]
countC<-pCC[Length[pCC]]
pCT<- posCount('T', sequence)
posT<- pCT[-1]
countT<-pCT[Length[pCT]]
pCG<- posCount('G',sequence)
posG<- pCG[-1]
countG<-pCG[Length[pCG]]
return((countG + countC)/(countA + countT))
} 

slidingwindowGCplot <- function(windowsize,inputseq){
GCwindow <- seq(1, length(inputseq) - windowsize, by = windowsize)
n <- length(GCwindow)
Chunks <- numeric(n)       
for(i in 1: Chunks){
	chunk<- inputseq[GCwindow[i]:(inputseq + windowsize - 1)]
	chunGC<- GC(chunk)
	print(chunkGC)
	chunks[i]<- chunkGC
}
plot(GCwindow,chunks,type="b",xlab="Nucleotide start position",ylab="GC content",main=paste("GC Plot with windowsize ", windowsize))
}
slidingwindowGCplot(20,gpl2872[1,])

#Protein Sequence Statistics                  
#ALternative
#with FASTA

wget -c http://noncode.org/datadownload/ncrna_NONCODE[v3.0].fasta.tar.gz
tar -xzf *.tar.gz 
mv *.fasta ncrna_noncode_v3.fa
cat ncrna_noncode_v3.fa | grep "^>" | wc -l
411553


# this time output a fasta

writeLines(d[[1]], con = "test.fasta");


url = "http://noncode.org/datadownload/ncrna_NONCODE[v3.0].fasta.tar.gz"
method = "wget"
download.file(url, destfile = "411553.fasta.tar.gz", model = "w", headers = NULL)
untar("D:/Introduction to R/411553.fasta.tar.gz", exdir ="D:/Introduction to R/411553")
cox1<- read.fasta(file = "411553.fasta", seqtyp = "AA")


#Retrieve a sequence from GenBank
#AB003468 <- read.GenBank("AB003468", as.character = "TRUE")
#Save as FASTA format
#write.dna(AB003468, file ="AB003468.fasta", format = "fasta", append = FALSE, nbcol = 6, colsep = " ", colw = 10)

#install.packages("rentrez")
#library(rentrez) 
#entrez_search(db="nucleotide", term="human superoxide dismutase")
#cox1<- read.fasta(file = "cox1.fasta", seqtyp = "AA")
#get the amino acid composition for all four
aaComp(cox1) 
#returns the aliphatic index of the protein sequence (an indicator of thermostability of globular proteins)
aIndex(cox1)
# predicts the net charge of the protein(charge requires additional input parameters:
charge(cox1)
#Recall from the Proteins course that different amino acids will have different charges based on
their R group.
#The pH of the medium that the amino acids are in will affect this charge; hence the need to include
the pH value in the charge command as well as a pKscale.)
charge(seq, pH = 7, pKscale = "Lehninger")
#Notethat: Here we use a different pKscale – there are several to choose from, and when you are writing your own analysis you should consult the Peptides documentation to see which is more appropriate.
#specify the sequence 
charge(seq="FLPVLAG", pH=7, pKscale="EMBOSS")
#submit a single amino acid to see the charge on that AA


# load GPL Annotation file
#gpl2872 <- getGEO(filename='GPL2872_old_annotations.txt', GSEMatrix=TRUE)

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

#pos = regexpr('pattern', x) # Returns position of 1st match in a string
#pos = gregexpr('pattern', x) # Returns positions of every match in a string
#pos = grep('pattern', x)
#keep = substr(x, first, last)
#sub('pattern', replacement, input) # Changes only the 1st pattern match per string
#gsub('pattern', replacement, input) # Changes every occurrence of a pattern match


