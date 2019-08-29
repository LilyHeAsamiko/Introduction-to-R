install.packages("installr") # install 
setInternet2(TRUE) # only for R versions older than 3.3.0
installr::updateR(T) 
remove.packages("data.table")
install.packages("data.table", type = "source",
repos = "http://Rdatatable.github.io/data.table")
###################################################################

getwd()
setwd("D:/Introduction to R")
setwd("P:/Introduction-to-R-master/Introduction-to-R-master")

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
install.packages("utils")
BiocManager::install(c("DESeq","topGO"))
install.packages("Peptides") 
BiocManager::install("Biostrings", version = "3.8")
install.packages("gplots")
install.packages("timeSeries")
install.packages("MASS")
install.packages("rgl")
#install.packages("CrossValidate", repos = "http://silicovore.com/OOMPA")
BiocManager::install("CrossValidate")
install.packages("random")
install.packages("RColorBrewer")
BiocManager::install(c("edgeR", "limma", "Glimma", "org.Mm.eg.db"))
BiocManager::install('easyRNASeq')

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
library(gplots)
library(ggplot2)
library(seqinr)
library(ape) 
library(stringr)
library(Peptides) 
library(Biostrings)
library(timeSeries)
library(MASS)
library(rgl)
library(CrossValidate)
library(random)
library(edgeR)
library(limma)
library(Glimma)
library(org.Mm.eg.db)
library(RColorBrewer)
library(Mus.musculus)

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
utils::untar("D:/Introduction to R/GSM206607.gz", exdir ="D:/Introduction to R/GSM206607")
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
T = 117.67, p-value < 2.6e-16
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
T = 5.2636, p-value < 1.2e-16
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
print(Mc)
return(MAData[,i]<- (M[,i] - Mc))})
boxplot(MAData, main = "MA_Normed Gene Expressions")
logData<- log2(rawData)
GSMstat["TotalMA"]<- summary(as.numeric(MAData))
write.csv(GSMstat, file = "GSMstat.csv", col.names = T)
d<- density(MAData)  
hist(MAData, freq = F)
lines(d) 
kurtosis.norm.test(MAData)
#######################################################
   Kurtosis test for normality

data:  MAData
T = 13.391, p-value < 2.2e-16
#######################################################
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
Chunks <- length(GCwindow)
#Chunks <- numeric(n)
chunks<-seq(1,Chunks)        
for(i in seq(1,Chunks)){
	chunk<- inputseq[GCwindow[i]:(GCwindow[i] + windowsize)]
	chunkGC<- GC(chunk)
	print(chunkGC)
	chunks[i]<- chunkGC
}
plot(GCwindow,chunks,type="b",xlab="Nucleotide start position",ylab="GC content",main=paste("GC Plot with windowsize ", windowsize))
}
slidingwindowGCplot(20,gpl2872[1,])

#Protein Sequence Statistics                  
#ALternative
#with FASTA
#wget -c http://noncode.org/datadownload/ncrna_NONCODE[v3.0].fasta.tar.gz
#tar -xzf *.tar.gz 
#mv *.fasta ncrna_noncode_v3.fa
#cat ncrna_noncode_v3.fa | grep "^>" | wc -l
#411553


# this time output a fasta

#writeLines(d[[1]], con = "test.fasta");

#Alternatively: Retrieve a sequence from GenBank
AB003468 <- read.GenBank("AB003468", as.character = "TRUE")
#Save as FASTA format
#write.dna(AB003468, file ="AB003468.fasta", format = "fasta", append = FALSE, nbcol = 6, colsep = " ", colw = 10)
#install.packages("rentrez")
#library(rentrez) 
#entrez_search(db="nucleotide", term="human superoxide dismutase")
#cox1<- read.fasta(file = "cox1.fasta", seqtyp = "AA")



url = "http://noncode.org/datadownload/ncrna_NONCODE[v3.0].fasta.tar.gz"
method = "wget"
download.file(url, destfile = "411553.fasta.tar.gz", model = "w", headers = NULL)
untar("D:/Introduction to R/411553.fasta.tar.gz", exdir ="D:/Introduction to R/411553")
ncrna<- read.fasta(file = "411553/ncrna_NONCODE[v3.0].fasta")
#Basic Information:
str(ncrna)
#List of 411553 dna sequences
#Look into one sequence: 19
seqid<-19
nc<- ncrna[[seqid]] nc
nc
###############################
attr(,"name")
[1] "n19"
attr(,"Annot")
[1] ">n19 | AB015469 | snRNA | Arabidopsis thaliana (thale cress) | U2 | NONCODE v2.0 | NULL | NULL | -1.2887100 | -0.2345336"
attr(,"class")
[1] "SeqFastadna"
#################################################
#count number of nucleotides
count(nc,1)
#############
 a  c  g  t 
47 45 40 67 
#############
#count number of dinucleotides
count(nc,2)
#################################################
aa ac ag at ca cc cg ct ga gc gg gt ta tc tg tt 
 9 10 10 17 11  8  7 19  8 12 11  9 18 15 12 22
################################################## 
#GC content
GC(nc)
[1] 0.4271357
#with sliding window size of 20
slidingwindowGCplot(20,nc)
###########################
[1] 0.4761905
[1] 0.3333333
[1] 0.2380952
[1] 0.4761905
[1] 0.3333333
[1] 0.4285714
[1] 0.3333333
[1] 0.5238095
[1] 0.5714286
###########################
# Useage of hydrophobicity to mesure amino acid sequences(similar to EMBOSS)
Hyph<- hydrophobicity(nc)
slidingwindowHyphplot <- function(windowsize,inputseq){
Hyphwindow <- seq(1, length(inputseq) - windowsize, by = windowsize)
Chunks <- length(Hyphwindow)
chunks<-seq(1,Chunks)        
for(i in seq(1,Chunks)){
	chunk<- inputseq[Hyphwindow[i]:(Hyphwindow[i] + windowsize)]
	chunkHyph<- hydrophobicity(chunk)
	print(chunkHyph)
	chunks[i]<- chunkHyph
}
plot(Hyphwindow,chunks,type="b",xlab="Nucleotide start position",ylab="Hydrophobicity",main=paste("Hydrophobicity Plot with windowsize ", windowsize))
}
slidingwindowHyphplot(20,nc)


#get the amino acid composition for nc
aaComp(nc) 
################################
[[1]]
          Number Mole%
Tiny           1   100
Small          1   100
Aliphatic      1   100
Aromatic       0     0
NonPolar       1   100
Polar          0     0
Charged        0     0
Basic          0     0
Acidic         0     0

###################################
#returns the aliphatic index of the protein sequence (an indicator of thermostability of globular proteins)
aliph<- aIndex(nc)

# predicts the net charge of the protein(charge requires additional input parameters:
Q<- charge(nc)

#Recall from the Proteins course that different amino acids will have different charges based on
their R group.
#The pH of the medium that the amino acids are in will affect this charge; hence the need to include
the pH value in the charge command as well as a pKscale.)
Q1<- charge(nc, pH = 7, pKscale = "Lehninger")

#Notethat: Here we use a different pKscale – there are several to choose from, and when you are writing your own analysis you should consult the Peptides documentation to see which is more appropriate.
#specify the sequence 
#submit a single amino acid to see the charge on that AA
Q2<- charge(seq="FLPVLAG", pH=7, pKscale="EMBOSS")

slidingwindowaaCompplot <- function(windowsize,inputseq){
window <- seq(1, length(inputseq) - windowsize, by = windowsize)
Chunks <- length(window)
chunkaliphs<-seq(1,Chunks) 
chunkQs<-seq(1,Chunks)
chunkQ1s<-seq(1,Chunks)
chunkQ2s<-seq(1,Chunks)       
for(i in seq(1,Chunks)){
	chunk<- inputseq[window[i]:(window[i] + windowsize)]
	chunkaliph<- aIndex(chunk)
	chunkQ<- charge(chunk)
	chunkQ1<- charge(chunk, pH = 7, pKscale = "Lehninger")
	chunkQ2<- charge(chunk, pH=7, pKscale="EMBOSS")
	chunkaliphs[i]<- chunkaliph
	chunkQs[i]<- chunkQ
	chunkQ1s[i]<- chunkQ1
	chunkQ2s[i]<- chunkQ2
}
write.csv(t(c(chunkaliphs,chunkQs, chunkQ1s, chunkQ2s)),file = "aaCompTb.csv", col.names = TRUE)
print(data.frame("Aliphatic" = chunkaliphs,"Net Charge" = chunkQs, "Charge with Lehninger" = chunkQ1s,"Charge with EMBOSS" = chunkQ2s))
layout(matrix(c(1,2,3,4), 2,2, byrow = TRUE))
plot(window,chunkaliphs,type="b",xlab="Nucleotide start position",ylab="Aliphatic",main=paste("Aliphatic Plot with windowsize ", windowsize))
plot(window,chunkQs,type="b",xlab="Nucleotide start position",ylab="Charge",main=paste("Net Charge Plot with windowsize ", windowsize))
plot(window,chunkQ1s,type="b",xlab="Nucleotide start position",ylab="Charge",main=paste("ChargePlot with ph: 7 and pKscale: Lehninger and windowsize ", windowsize))
plot(window,chunkQ2s,type="b",xlab="Nucleotide start position",ylab="Charge",main=paste("ChargePlot with ph: 7 and pKscale: EMBOSS and windowsize ", windowsize))
}

slidingwindowaaCompplot(20, nc)
####################################
  Aliphatic   Net.Charge Charge.with.Lehninger Charge.with.EMBOSS
1       100 -0.002015701          -0.002015701        -0.02410542
2         0 -0.002015701          -0.002015701        -0.02410542
3         0 -0.063990417          -0.063990417        -0.05475885
4         0 -0.002015701          -0.002015701        -0.02410542
5         0 -0.063990417          -0.063990417        -0.05475885
6         0 -0.002015701          -0.002015701        -0.02410542
7       100 -0.002015701          -0.002015701        -0.02410542
8         0 -0.002015701          -0.002015701        -0.02410542
9         0 -0.002015701          -0.002015701        -0.02410542
######################################################################
#Sequence Alignment
#Sequence alignments are generally the first step in any comparative analysis
#takes two or more sequences that are suspected of being similar and "lines
#them up", one on top of the other, matching similar residues (bases) 
#and inserting gaps (dashes) where necessary.

#"mismatches" – residues or bases that do not match but are
#"orthologous" – meaning that they are known to be in the same position in both sequences despite not being identical.

#The columns represent orthologous residues same position in the protein, thus ostensibly the
same function.
#The rows represent the different species – with the exception of the last row which indicates the
homology between the species.
#This can also be represented as the "consensus" sequence – meaning the sequence of agreement,
or the sequence in common between all represented species.
#Alignments do not have to be between sequences from different species – some genes within the
same species will also share similarity.
#Pairwise alignments can be used to analyze two sequences directly.
#Pairwise alignments can also be an intermediate step to creating more complex "multiple
   sequence alignments", or MSA’s.

url = "https://prod-edxapp.edx-cdn.org/assets/courseware/v1/c309caedb1530152a48468b74c2f03b3/asset-v1:USMx+BIF003x+3T2016+type@asset+block/prok.fasta"
download.file(url, destfile = "seq.fasta", model = "w", headers = NULL)
prokaryotes <- read.fasta(file = "seq.fasta", seqtype = "DNA")
seq1<-as.character(prokaryotes[[1]])
seq1 = paste(seq1, collapse = "")
seq2<-as.character(prokaryotes[[2]])
seq2 = paste(seq2, collapse = "")
pair21<- pairwiseAlignment(pattern=seq2, subject=seq1)
summary(pair21) 
########################################################
pattern: atgcccaagacgca-atctgccgcaggctataag...gttcgaattcgaaacgatggacaagctc---taa
subject: atg---ataaccctgacctaccgcat-cgaaacg...attc--------agcgacgg-caagggcgcgtga
score: -3134.096

Global Single Subject Pairwise Alignments
Number of Alignments:  1

Scores:
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  -3134   -3134   -3134   -3134   -3134   -3134 

Number of matches:
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    753     753     753     753     753     753 

Top 10 Mismatch Counts:
   SubjectPosition Subject Pattern Count Probability
1                5       t       a     1           1
2                6       a       g     1           1
3                9       c       g     1           1
4               11       t       a     1           1
5               14       c       t     1           1
6               17       a       g     1           1
7               23       t       g     1           1
8               25       g       t     1           1
9               27       a       t     1           1
10              29       c       a     1           1
>
############################################################################################
pair21S = BStringSet(c(toString(subject(pair21)), toString(pattern(pair21))))
writeXStringSet(pair21S, "alignedseq21.txt", format="FASTA")
#DotPlot
#One of the oldest sequence comparison visuals is a dotplot.
#In its simplest form, a dotplot is generated by placing a dot at position (i,j) if character number i in
#the first sequence is the same as character number j in the second sequence.
#Dotplots can be used to visually see how similar sequences are – and to see if there are
duplication, inversions, repeats, etc.. between the two.

AB015469 | snRNA | Arabidopsis thaliana (thale cress) 

seqid2<- 16
nc2<- ncrna[[seqid2]]
#####################################
[1] "n16"
attr(,"Annot")
[1] ">n16 | AB014881 | Bsr RNA | Rattus norvegicus (Norway rat) | Bsr RNA | NONCODE v2.0 | NULL | NULL | -0.8771300 | -0.3362236"
attr(,"class")
[1] "SeqFastadna"
#####################################################################
nc<- as.character(nc)
nc2<- as.character(nc2)
#nc[is.na(nc)]<-0
#nc2[is.na(nc2)]<-0
layout(matrix(c(1,2,3), 1, 3, byrow = TRUE))
dotPlot(nc, nc2[1:length(nc)], main = "AB015469|(snRNA) Arabidopsis thaliana(thale cress)  vs AB014881|(Bsr RNA)Rattus norvegicus(Norway rat) Dotplot")
#to improve the signal to noise ratio: adapt wsize, wstep, nmatch
dotPlot(nc, nc2, wsize = 3, wstep = 3, nmatch = 3, main = "AB015469|(snRNA) Arabidopsis thaliana(thale cress) vs AB014881|(Bsr RNA)Rattus norvegicus(Norway rat) Dotplot\nwsize = 3, wstep = 3, nmatch = 3")
#try the first 100 residues of the amino acid sequences
dotPlot(nc[1:100], nc2[1:100], wsize = 3, wstep = 3, nmatch = 3,  main = "AB015469|(snRNA) Arabidopsis thaliana(thale cress) vs AB014881|(Bsr RNA)Rattus norvegicus(Norway rat) first 100 AA Dotplot\nwsize = 3, wstep = 3, nmatch = 3")

#Note that:
#type = "<type>" to tell pairwiseAlignment which type of alignment you are using, e.g.
type = "local".
#There are other options for pairwise alignment that let you control the substitution matrix to use
(submat = "<type>") as well as the gap opening and extension penalties (gapOpening = X,
gapExtension = Y).

From your DNA and Proteins classes you should remember the difference between local and
global alignments:

#In a local alignment, you are matching your query with a substring (fragment) of your subject
sequence. 
#In a global alignment you perform an end to end alignment between the two.
üYou may end up with a lot of gaps in global alignment if the sizes of query and subject are
dissimilar).
#Local alignments may also have gaps.
#There is in fact a third option – the overlap alignment.  Overlap alignments are used when you are
#trying to assemble overlapping sequences, e.g. from multiple sequencing runs in a genome
assembly.
#Overlap alignments focus on making the best alignments between the end of one sequence and
the end of another.

#Protein Strcture, Expression, and Interaction
# mass-spec (MS)
#mass-spec procedures such as MALDI-TOF produce "time of
#flight" data, where each protein fragment is charged and propelled towards a detector:
#Larger fragments go shorter distances; smaller fragments go further.
#By looking at all possible fragment sizes generated, software can assemble a list of all probable
#proteins in a sample of mixed proteins – as well as their quantity.
url<- "https://prod-edxapp.edx-cdn.org/assets/courseware/v1/9191814573c2eb0c37ca5ddec221bda7/asset-v1:USMx+BIF003x+3T2016+type@asset+block/peptidefrags.txt"
download.file(url,destfile = "peptidefrags.txt", model = "w", headers = NULL)
Peptides<- read.table("peptidefrags.txt", header=TRUE)
id<- names(Peptides)
#peptides<- as.vector(Peptides$id)
peptides<- as.numeric(Peptides[id][[1]])
#QC test
hist(peptides,breaks=400) 
kurtosis.norm.test(peptides)

#url<- "https://prod-edxapp.edx-cdn.org/assets/courseware/v1/c66ff508fc42acd0e0b79da26792a203/asset-v1:USMx+BIF003x+3T2016+type@asset+block/mascot.txt"
#download.file(url,destfile = "mascot.txt", model = "w", headers = NULL)
#Mascot<- read.table("mascot.txt", header=TRUE)


GSM206617<- read.csv('GSM206617.csv', header= TRUE, skip = 7, fill = TRUE)
GSM206617dt<- GSM206617[,2:100]
GSM206617Var<- colnames(GSM206617dt[1,])
GSM206617dt[,11]
GSM206618<- read.csv('GSM206618.csv', header= TRUE, skip = 7, fill = TRUE)
GSM206618dt<- GSM206618[,2:100]
GSM206618Var<- colnames(GSM206618dt[1,])
GSM206618dt[,11]
GSM206619<- read.csv('GSM206617.csv', header= TRUE, skip = 7, fill = TRUE)
GSM206619dt<- GSM206619[,2:100]
GSM206619Var<- colnames(GSM206619dt[1,])
GSM206619dt[,11]

transferline<- function(l){
if(is.na(l[98]) && is.na(l[99])){
	l[6:98]<- l[5:97]
	l[5]<- 0
	l[14:99]<- l[13:98]
	l[13]<- 0}else{
	l[10:99]<- l[9:98]
	l[9]<- 0}
return(l)
}
transferdt<- function(dt){
sapply(seq(1,length(dt)/99), function(i) dt[i,]<-transferline(dt[i,]))
return(dt)
}

transferdt(GSM206617dt)
transferdt(GSM206618dt)
transferdt(GSM206619dt)
GSM206617dt[is.na(GSM206617dt)]<- 0
GSM206618dt[is.na(GSM206618dt)]<- 0
GSM206619dt[is.na(GSM206619dt)]<- 0

GSMlabel<- list(GSM206617=GSM206617dt[,10], GSM206618=GSM206618dt[,10], GSM206619=GSM206619dt[,10])
venn(GSMlabel)
#f_s<- factor(GSMlabel)
#fs_l<-levels(f_s)
#nlevels(fs_l)
#summary(f_s)
GSMvar <- list(GSM206617=GSM206617dt[,11], GSM206618=GSM206618dt[,11], GSM206619=GSM206619dt[,11])
venn(GSMvar)
GSMvar1 <- list(GSM206617=GSM206617dt[,1], GSM206618=GSM206618dt[,1], GSM206619=GSM206619dt[,1])
venn(GSMvar1)
GSMvar2 <- list(GSM206617=GSM206617dt[,5], GSM206618=GSM206618dt[,5], GSM206619=GSM206619dt[,5])
venn(GSMvar2)

#LDA Analysis: comparison of proteins found between different sampling conditions(disease and control) 
GSMtrain<- data.frame("GSMLabel" = as.numeric(unlist(GSMlabel)), "GSMVar" = as.numeric(unlist(GSMvar)))
#GSMtrain[is.na(GSMtrain)]<- 0
#pred<- predict(ldam, GSMtest)
set.seed(123)
id<- randomSequence(1,10000)
#cannot allocate vector of size 33.3 Gb
#GSMtest<- c(GSMtrain,as.numeric(unlist(GSMvar2))) 
GSMtr_df<- as.data.frame(GSMtrain[id,])
pred<- predict(ldam, GSMtr_df)
x<- GSMtrain[id,1]
y<- pred
class<- as.vector(unlist(GSMtr_df[1]))
y_dt<-  as.vector(unlist(y[3]))
plotdata<-data.frame(class, x, y_dt)

#generates a distance matrix from the dataframe centroids using straight-line (euclidean)geometry – and only generates the "upper" part of the matrix.
centroids <- aggregate(cbind(x,y_dt)~class,plotdata,mean)
CentroidDistances <- dist(centroids, method = "euclidean", diag = TRUE, upper = FALSE, p = 2)
attr(CentroidDistances, "Labels") <- centroids$class
plot1 <- ggplot(plotdata,aes(x,y_dt,color=factor(class))) + geom_point(size=3)
plot1
ggplot(plotdata,aes(x,y_dt,color=factor(class))) + geom_point(size=3)+ geom_point(data=centroids,size=7)
ggplot(plotdata,aes(x,y_dt,color=factor(class))) + geom_point(size=3)+ geom_point(data=centroids,size=7) + geom_text(data=centroids, size=7, label=centroids$class, colour="black")
ggplot(plotdata,aes(x,y_dt,color=factor(class))) + geom_point(size=3) + geom_point(data=centroids,size=7) + geom_text(data=centroids, size=7, label=centroids$class, colour="black") + ggtitle("LDA of Conditions 1-")
ggplot(plotdata,aes(x,y_dt,color=factor(class))) + geom_point(size=3) + geom_point(data=centroids,size=7) + geom_text(data=centroids, size=7, label=centroids$class, colour="black") + ggtitle("LDA of Conditions 1-") + geom_text(aes(label=MS),hjust=0, vjust=0, colour="black")
plot1 <- ggplot(plotdata,aes(x,y_dt,color=factor(class))) + geom_point(size=3) + geom_point(data=centroids,size=7) + geom_text(data=centroids, size=7, label=centroids$class, colour="black") + ggtitle("LDA of Conditions 1-3") + geom_text(aes(label=factor(class)),hjust=0, vjust=0, colour="black")
ggsave(filename="plot1.pdf")

write.csv(as.matrix(CentroidDistances, file = "centroiddistances.csv"))

#Proteomics
#Microarray Analysis
A Microarray is easily one of the densest ways of collecting bioinformatics data, whether it is Protein
or DNA-based.
#One of the most common types of arrays is in an expression analysis, where the "spots" on the
"chip" represent different genes, and mRNA isolated from an experimental sample is "washed" over
this:
The RNA binds to its complementary DNA.
Probes on the RNA that can be used to detect binding levels, e.g. "spot intensity" or expression 
levels.

url = "https://prod-edxapp.edx-cdn.org/assets/courseware/v1/3d701221d05b85ae97b9cd4ce238add6/asset-v1:USMx+BIF003x+3T2016+type@asset+block/brain.fetalbrain.2color.data__1_.txt"
download.file(url, destfile = "Brain.fetalbrain.2color.dat.txt", model = "w", header = NULL)

affy.data <- ReadAffy()
#normalize the data
eset.mas5 <- mas5(affy.data)

# +1/-1, to represent fold change in expression.
#We can do that by using a log transform of the data using the log command that is standard in R.
#Moreover, if we use log base 2 transforms, we end up with a binary representation – perfect for our needs.
exprEset.nologs <- exprs(eset.mas5)
exprEset <- log(exprEset.nologs, 2)
write.table(exprEset, file="mouse4302_mas5_matrix.txt", quote=F, sep="\t")
#ceilfile = 
#Not all of our microarray spots will be "populated" with data in any experiment – so it is useful to do something called an "Absent/Present" check.
#"Absent/Present" check.
#generates a vector with the A/P values for the dataset through computing p-values
and Making P/M/A Calls on probe level data
data.mas5calls <- mas5calls(affy.data)
# converts this into an expression matrix containing A or P for each tissue and gene combination.
data.mas5calls.calls <- exprs(data.mas5calls)
Brain.fetalbrain.2color <- read.maimages("Brain.fetalbrain.2color.dat.txt", columns=list(G="brain.1",
R="fetal.brain.1", Gb="bg1", Rb="bg2"))
#normalize the color data 
Brain.fetalbrain.2color.loess <- normalizeWithinArrays(Brain.fetalbrain.2color, method="loess")
par(mfrow=c(1,2))
plotMA(Brain.fetalbrain.2color)
plotMA(Brain.fetalbrain.2color.loess)

#result of using multiple chips (e.g. each tissue was run twice, so brain.1 and brain.2 are
the results of two separate chip runs
#The value of 1 following the tissue columns names (e.g. 1, means) is the MARGIN for the
#command, which indicates that it is the rows within each column that the mean is calculated for
brain1_mean <- apply(exprEset[, c(1, 2)], 1, mean)
brain12_mean <- apply(exprEset[, c(3, 5)], 1, mean)
brain2_mean <- apply(exprEset[, c(7, 9)], 1, mean)
brain22_mean <- apply(exprEset[, c(11, 12)], 1, mean)

brain1_mean_1to2<- brain1_mean - brain12_mean
brain2_mean_1to2<- brain2_mean- brain22_mean
all_data <- cbind(exprEset, brain1_mean, brain12_mean, brain2_mean,  brain22_mean, brain1_mean_1to2,brain2_mean_1to2)
write.table(all_data, file = "Microarray_ALL.txt", quote=F, sep="\t")

#see if there are any differentially expressed genes
#T-test
dataset1<- exprEset[1, c(1, 2)]
dataset2<- exprEset[1, c(3, 5)]
t_test_brain_1 <- t.test(dataset1, dataset2, "two.sided")
-------------------------------------------------------------

        Welch Two Sample t-test

data:  dataset1 and dataset2
t = -0.71672, df = 1.5391, p-value = 0.5665
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.8109782  0.6327863
sample estimates:
mean of x mean of y 
 10.99189  11.08099 

---------------------------------------------------------------------
brain1_p_value_all<- apply(exprEset, 1, function(x){ t.test(x[1:2], x[c(3,5)])$p.value } )

brain2_p_value_all<- apply(exprEset, 1, function(x) { t.test(x[7:9], x[11:12]) $p.value } )

AP <- apply(data.mas5calls.calls, 1, paste, collapse="")
genesP = names(AP[AP != "AAAAAAAA"])
length(genesP)
exprEsetP <- exprEset[genesP,]
 
#adjust for the FDR, or false discovery rate(FDR correction is necessary when looking at many thousands of tests as a way of minimizing the
#false positives that inevitably are found when so many samples are analyzed, e.g. 5% of the time a
95% confidence interval is wrong)

brain1_p_valueP <- brain1_p_value_all[genesP]
brain2_p_valueP <- brain2_p_value_all[genesP]

#The p values are automatically adjusted to compensate for the false discovery rate 
brain1_fdr_p_valsP <- p.adjust(brain1_p_valueP, method="fdr")
brain2_fdr_p_valsP <- p.adjust(brain2_p_valueP, method="fdr")

brain1_fdr_p_valsP_sorted <- brain1_fdr_p_valsP[order(brain1_fdr_p_valsP)]
brain2_fdr_p_valsP_sorted <- brain2_fdr_p_valsP[order(brain2_fdr_p_valsP)]
brain1_DE_probesets <- names(brain1_p_valueP[brain1_p_valueP < 0.01])
brain2_DE_probesets <- names(brain2_p_valueP[brain2_p_valueP < 0.01])

brain1_DE_log2ratios <- all_data[brain1_DE_probesets, c("brain1_mean_1to2", "brain2_mean_1to2")]
brain2_DE_log2ratios <- all_data[brain2_DE_probesets, c("brain1_mean_1to2", "brain2_mean_1to2")]

write.table(brain1_DE_log2ratios, "brain1_DE_log2ratios.txt", sep="\t", quote=F)
write.table(brain2_DE_log2ratios, "brain2_DE_log2ratios.txt", sep="\t", quote=F)
#Any gene that has all A values is uninformative (the data is absent for all chips)

#Note that:a two-sided, or two-tailed T test – meaning we are looking for outliers on both ends of the spectrum (the outer 5% of both ends of the normal distribution)

#Method for correcting the false positive rate: Menjamini Hochberg Method:
#To choose how many genes to call significantly differentially expressed
#p<= x*a/m 
#where p is the larest p value you will call significant, alpha is the ideal FDR(ex. 0.05), x is the number of the gene you will call significant, m is the number for the hypothesis test(i.e. total gene numbers))
#



#layout(matrix(c(1,2)),1,2,byrow = TRUE)
par(mflow = c(1, 2))
xdata<- all_data[, "brain1_mean"]
ydata<- all_data[, "brain2_mean"]
plot(xdata, ydata, main = "Log2 expression in brain1 (n=2) vs brain2 (n=2)", xlab="brain 1", ylab="brain 2", col="blue", cex=0.5)
abline(0,1)
#MA
M<- ydata - xdata
A<- 0.5*(ydata + xdata)
plot(A, M, main = "MA in brain (n=2) vs adult brain (n=2)", xlab="brain 1", ylab="brain 2", pch = 19, col="red", cex=0.5)
abline(h =c( -1,1))

#a volcano plot is simply the p-values on the Y axis and the fold-difference on the X.
#a quick log base 10 transform of the raw p-values is conducted to compare; this
gives us an easier to read plot (albeit the Y axis is now a log axis)
expr_pvals<- cbind(exprEsetP, brain1_p_valueP, brain1_fdr_p_valsP, brain2_p_valueP, brain2_fdr_p_valsP)
log2_ratios<- expr_pvals[, 4] -  expr_pvals[, 2]
log2_ratios_fdr<- expr_pvals[, 5] -  expr_pvals[, 3]
p_values<- expr_pvals[, "brain1_p_valueP"]
p_values_fdr<- expr_pvals[, "brain1_fdr_p_valsP"]


par(mfrow=c(2,2))
plot(log2_ratios, p_values)
plot(log2_ratios, -log(p_values, 10))
plot(log2_ratios_fdr, p_values_fdr)
plot(log2_ratios_fdr, -log(p_values_fdr, 10))

#RNA sequencing
RNAseqdata <- read.table("GSE60450_Lactation-GenewiseCounts.txt", stringsAsFactors = FALSE, fill = TRUE)
sampleinfo <- read.table("SampleInfo.txt", fill = TRUE)

RNAcount <- RNAseqdata[,-(1:2)]
rownames(RNAcount) <- RNAseqdata[,1]
colnames(RNAcount) <- substr(RNAcount[1,],start=1,stop=10)#12
# filter the data to eliminate genes expressed at very low levels: Counts Per Million(CPM) values
myCPM<- RNAcount[-1,]
myCPM<- sapply(seq(2,nrow(RNAcount)), function(i) myCPM[i,]<- cpm(as.numeric(as.vector(RNAcount[i,]))))

#there's a large discrepancy between the largest and the smallest values

thresh<- myCPM > 0.5

#if use all variables of rawData: filledcols = colSds(GSMtrain) != 0.0 to clean empty column
f<- factor(GSMtrain$GSMLabel)
hf<- hist(as.numeric(f))
hf_br<- hf$breaks
hf_c<- hf$counts
hf_d<- hf$density
hf_mids<- hf$mids 
f_l<- levels(f)
nf_l<- nlevels(f) #33993
sapply(seq(1, length(hf$counts)), function(i) {
c<- (GSMtrain[,1] > hf_mids[i] | GSMtrain[,2] > hf_mids[i]) & (GSMtrain[,1] > hf_mids[i] |GSMtrain<hf_mids[i+1])
cols<- (1: hf_mids[i])[c]
rows<- cols[!is.na(cols)]
DT<- GSMtrain[rows,]
fn<- as.numeric(f)
if(is.null(f[!is.na(fn[rows])])) ldam<- lda(as.matrix(f[rows]<- ifelse(as.numeric(f[!is.na(fn[rows])]) <hf_mids[i], 1, 0)), data = DT, na.action ="na.omit", CV= TRUE)
else ldam<- lda(rep(0,length(DT)), data = DT, na.action ="na.omit", CV= TRUE)
pred<- predict(ldam,DT)})
#cv <- CrossValidate(ldam, t(GSMtrain), factor(GSMtrain$GSMLabel), 0.5, nLoop = 5)
#summary(cv)

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
///////////////////////////
#pos = regexpr('pattern', x) # Returns position of 1st match in a string
#pos = gregexpr('pattern', x) # Returns positions of every match in a string
#pos = grep('pattern', x)
#keep = substr(x, first, last)
#sub('pattern', replacement, input) # Changes only the 1st pattern match per string
#gsub('pattern', replacement, input) # Changes every occurrence of a pattern match


######
##low level access with mzR::openMSfile
##convenient access with MSnbase::readMSData
##on disk access(as opposed to in memory)
##
##
#library("msdata")
fls<- proteomics(full.names = TRUE)
basename(fls)
fl<- fls[2]
fl

library(mzR)  #Rcpp
rw<- openMSfile(fl)
rw

sp1<- spectra(rw, 1)
head(sp1)
sp1<- spectra(rw, 1:2)
length(sp1)

hd<- header(rw) #dataframe
head(hd)

suppressPackageStartupMessages(library("MSnbase"))
mse<- readMSData(fl, mode = "onDisk") #earth data function returns df containing information similar to BIA
mse[[1]]
fData(mse)
mse
