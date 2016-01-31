library(edgeR)
library(DEXSeq)


DIR.IN<-"dexseq_counts/"
DIR.OUT<-paste("dexseq_analysis/")

sampleTable<-as.data.frame(read.table("sample_info.txt",header=T))
sampleTable<-sampleTable[-5,] #samples 46 is missing
levels(sampleTable$Line)<-c("C","C","S","S")

sampleTable_line<-sampleTable[,-c(3,4)]
sampleTable_trt<-sampleTable[,-c(2,4)]
sampleTable$Line.Trt<-as.factor(paste(sampleTable$Line,sampleTable$Treatment,sep="."))
sampleTable$Line.Trt.Time<-as.factor(paste(sampleTable$Line,sampleTable$Treatment,sampleTable$Time,sep="."))
sampleTable_LineTrt<-sampleTable[,-c(2:4,6)]
sampleTable_LineTrtTime<-sampleTable[,-c(2:5)]

colnames(sampleTable_line)[2]<-"condition"
colnames(sampleTable_trt)[2]<-"condition"
colnames(sampleTable_LineTrt)[2]<-"condition"
colnames(sampleTable_LineTrtTime)[2]<-"condition"

ecs <- read.HTSeqCounts(countfiles = file.path(DIR.IN, sampleTable_LineTrtTime$Samples),
 design = sampleTable_LineTrtTime,
 flattenedfile = "dmel_flattened.gff" )

#Estimate factor size
library(devtools)
library(statmod)

ecs <- estimateSizeFactors( ecs )
sizeFactors(ecs)
head( counts(ecs), )
head( fData(ecs), )
design(ecs)

##Estimate Dispersion and make table
ecs <- estimateDispersions(ecs)
ecs <- fitDispersionFunction(ecs)
ecs <- testForDEU(ecs, nCores=4)
res1 <- DEUresultTable(ecs)

# make annotation files for groups and genes
genes <-data.frame(FB=res1$geneID)
egFLYBASE<-toTable(org.Dm.egFLYBASE)
egSYMBOL<-toTable(org.Dm.egSYMBOL)
m<-match(genes$FB, egFLYBASE$flybase_id)
genes$EG<-egFLYplotDEXSeq( ecs2, "FBgn0040734", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )BASE$gene_id[m]
m2 <- match(genes$EG, egSYMBOL$gene_id)
genes$Symbol<-egSYMBOL$symbol[m2]
res1$geneSymbol<-genes$Symbol

##To find out how many exons and genes were significant
table(res1$padjust < 0.05)
res1$geneSymbol[which(res1$padjust<0.05)]
res1$geneID[which(res1$padjust<0.05)]

#Plot the significant genes for Ctlvspar
pdf("DexSeqplotCG15065.pdf")
plotDEXSeq( ecs, "FBgn0040734", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
dev.off()

pdf("DexSeqplotTep2.pdf")
plotDEXSeq( ecs, "FBgn0041182", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
dev.off()

pdf("DexSeqplotCG43133.pdf")
plotDEXSeq( ecs, "FBgn0262607", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
dev.off()

pdf("DexSeqplotFok.pdf")
plotDEXSeq( ecs, "FBgn0263773", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
dev.off()

#Plot the significant genes for CvsS-CtalvsPar
 plotDEXSeq( ecs, "FBgn0003175", displayTranscripts=TRUE, legend=TRUE,cex.axis=1.2, cex=1.3, lwd=2 )


