library(edgeR)
library(DEXSeq)
#library("org.Dm.eg.db")

#Get the counts by running the following command
#"samtools view nsorted_sample4dedup.bam | python /usr/local/lib/R/site-library/DEXSeq/python_scripts/dexseq_count.py --paired=yes -s no ../../../references/dmel_flattened_ed.gff - sample4dexseq.counts"

wd<-getwd()

DIR.IN<-"~/RNAseq/dexseq_counts/"#"~/PhDpapers/RNAseq/dexseq_counts/mel"
DIR.OUT<-#paste("~/PhDpapers/RNAseq/dexseq_analysis")

sampleTable<-as.data.frame(read.table(paste(DIR.IN,"sample_info.txt",sep=""),header=T))
levels(sampleTable$Line)<-c("C","C","S","S")
sampleTable$libType<-rep(c("a","b","c"),16)


#Test variables one at the time to test their effect
sampleTable_trt<-sampleTable[,-c(2,4)]
colnames(sampleTable_trt)[2]<-"condition"

formulaFullModel= ~ sample + exon + libType:exon + condition:exon
formulaReducedModel = ~ sample + exon + libType:exon


dxd_trt = DEXSeqDataSetFromHTSeq(
countfiles=file.path(DIR.IN, sampleTable_trt$Samples),
sampleData=sampleTable_trt,
design= ~ sample + exon + condition:exon,
flattenedfile= "dmel_flattened.gff")

library(BiocParallel)
BPPARAM = MulticoreParam(workers=8)
dxd_trt = estimateSizeFactors( dxd_trt )
dxd_trt = estimateDispersions( dxd_trt,BPPARAM=BPPARAM)
dxd_trt = testForDEU( dxd, BPPARAM=BPPARAM)
dxd_trt = estimateExonFoldChanges(dxd_trt, BPPARAM=BPPARAM)
dxr1_trt = DEXSeqResults( dxd_trt )
table ( dxr1_trt$padj < 0.1 )

#FALSE  TRUE 
#70773     7
table ( tapply( dxr1_trt$padj < 0.1, dxr1_trt$groupID, any ) )

#FALSE  TRUE 
#11497     4 
#FBgn0040734 (CG15065),FBgn0041182 (Tep2),FBgn0262607 (CG43133),FBgn0263773 (fok)

plotMA( dxr1_trt, cex=0.8 )

#### Test Batch effect (Full vs Reduced model)

formulaFullModel = ~ sample + exon + libType:exon + condition:exon
formulaReducedModel = ~ sample + exon + libType:exon

dxd_trt = estimateDispersions( dxd_trt, formula = formulaFullModel, BPPARAM=BPPARAM)
dxd_trt = testForDEU( dxd_trt,reducedModel = formulaReducedModel,fullModel = formulaFullModel,BPPARAM=BPPARAM)
dxr2_trt = DEXSeqResults( dxd_trt )

#significant DEU after including batch
table( dxr2_trt$padj < 0.1 )

#FALSE  TRUE 
#56217     5 

dxr2_trt$groupID[which(dxr2_trt$padj<0.1)]
#"FBgn0041182" "FBgn0262607" "FBgn0263773" 

#### Test time effect alone 
sampleTable_time<-sampleTable[,-c(2,3)]
colnames(sampleTable_time)[2]<-"condition"

dxd_time = DEXSeqDataSetFromHTSeq(
countfiles=file.path(DIR.IN, sampleTable_time$Samples),
sampleData=sampleTable_time,
design= ~ sample + exon + condition:exon,
flattenedfile= "dmel_flattened.gff")

BPPARAM = MulticoreParam(workers=8)
dxd_time = estimateSizeFactors( dxd_time )
dxd_time = estimateDispersions( dxd_time,BPPARAM=BPPARAM)
dxd_time = testForDEU( dxd_time, BPPARAM=BPPARAM)
dxd_time = estimateExonFoldChanges(dxd_time, BPPARAM=BPPARAM)
dxr1_time = DEXSeqResults( dxd_time )
table ( dxr1_time$padj < 0.1 )

#FALSE  TRUE 
#69315  1465 

#length(unique(dxrtime$groupID[which(dxrtime$padj<0.1)]))
#798

#Include treatment, with libType: par,ctl

sampleTable_time$libType<-sampleTable$Treatment

formulaFullModel = ~ sample + exon + libType:exon + condition:exon
formulaReducedModel = ~ sample + exon + libType:exon

dxd = estimateDispersions( dxd, formula = formulaFullModel, BPPARAM=BPPARAM)
dxd = testForDEU( dxd,reducedModel = formulaReducedModel,fullModel = formulaFullModel,BPPARAM=BPPARAM)
dxr2 = DEXSeqResults( dxd )

table ( dxr2$padj < 0.1 )

#FALSE  TRUE 
#69379  1401 

length(unique(dxr1$groupID[which(dxr2$padj<0.1)]))
#761

#### Test time effect with parasitization, correcting for line ####
condition<-paste(sampleTable$Treatment,sampleTable$Time,sep='.')
sampleTable_trtxtime<-data.frame(Samples=sampleTable$Samples,condition,libType=sampleTable$Line)

dxd = DEXSeqDataSetFromHTSeq(
countfiles=file.path(DIR.IN, sampleTable_trt$Samples),
sampleData=sampleTable_trtxtime,
design= ~ sample + exon + condition:exon,
flattenedfile= "dmel_flattened.gff")

BPPARAM = MulticoreParam(workers=8)
dxd = estimateSizeFactors( dxd )
dxd = estimateDispersions( dxd,BPPARAM=BPPARAM)
dxd = testForDEU( dxd, BPPARAM=BPPARAM)
dxd = estimateExonFoldChanges(dxd, BPPARAM=BPPARAM)
dxr1 = DEXSeqResults( dxd )
table ( dxr1$padj < 0.1 )

#FALSE  TRUE 
#69315  1465 

length(unique(dxr1$groupID[which(dxr1$padj<0.1)]))
#798

formulaFullModel = ~ sample + exon + libType:exon + condition:exon + condition:Time
formulaReducedModel = ~ sample + exon + libType:exon + condition:exon

##Conclusion so far: time introduces most of the variation, which is not probably related to treatment


#### Test selection effect alone 
sampleTable_SC<-sampleTable[,-c(3,4)]
colnames(sampleTable_SC)[2]<-"condition"

formulaFullModel= ~ sample + exon + libType:exon + condition:exon
formulaReducedModel = ~ sample + exon + libType:exon

dxd_SC = DEXSeqDataSetFromHTSeq(
countfiles=file.path(DIR.IN, sampleTable_SC$Samples),
sampleData=sampleTable_SC,
design= ~ sample + exon + condition:exon,
flattenedfile= "dmel_flattened.gff")

dxd_SC = estimateSizeFactors( dxd_SC )
dxd_SC = estimateDispersions( dxd_SC,BPPARAM=BPPARAM)
dxd_SC = testForDEU( dxd_SC, BPPARAM=BPPARAM)
dxd_SC = estimateExonFoldChanges(dxd_SC, BPPARAM=BPPARAM)
dxr1_SC = DEXSeqResults( dxd_SC )
table ( dxr1_SC$padj < 0.1 )

#No significant genes

#Let's look at the batch effect
dxd_SC = estimateDispersions( dxd_SC, formula = formulaFullModel, BPPARAM=BPPARAM)
dxd_SC = testForDEU( dxd_SC,reducedModel = formulaReducedModel,fullModel = formulaFullModel,BPPARAM=BPPARAM)
dxr2_SC = DEXSeqResults( dxd_SC )

length(unique(dxr2_SC$groupID[which(dxr2_SC$padj<0.1)]))
2
unique(dxr2_SC$groupID[which(dxr2_SC$padj<0.1)])
#"FBgn0014469" (CG2060)  "FBgn0038247" (CG3389)

plotDEXSeq( dxr2_SC, "FBgn0014469", expression=FALSE, splicing=TRUE,
legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )

plotDEXSeq( dxr2_SC, "FBgn0038247", expression=FALSE, splicing=TRUE,
legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
