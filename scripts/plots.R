#rm(list = ls())

source('/home/laura/PhDpapers/RNAseq/scripts/script_Dmel_glmLRTfit.edgeR.R')
source('/home/laura/PhDpapers/RNAseq/scripts/script_Dspp_all_glmLRTfit.edgeR.R')
source('/home/laura/PhDpapers/RNAseq/scripts/script_Dspp_specific_glmLRTfit.edgeR.R')


###################################
########### MDS plots #############
###################################
## Spp 5 hours
sp.targets<-targets
sp.targets$group<-paste(sp.targets$line,sp.targets$treatment,sep=".") # sp.group
points <- c(15,17)
colors <- c("blue", "darkgreen", "red", "orange")
pdf("figures/MDS_allspp_5h.pdf",height=6,width=6)
plotMDS(sp.dge.5, col=colors[sp.targets.5$line], 
	pch=points[sp.targets.5$treatment],
	main= "MDS all species 5h",cex=1.25)
legend("topright", legend=paste(rep(levels(sp.targets.5$line),each=2),
			  	    levels(sp.targets.5$treatment),sep="."),
	col=colors[rep(as.factor(levels(sp.targets.5$line)),each=2)], pch=points, ncol=2)
#labels=sp.targets$sample)
dev.off()
## Spp 50 hours
pdf("figures/MDS_allspp_50h.pdf",height=6,width=6)
plotMDS(sp.dge.50, col=colors[sp.targets.50$line], 
	pch=points[sp.targets.5$treatment],
	main= "MDS all species 50h",cex=1.25)
legend("bottomright", legend=paste(rep(levels(sp.targets.5$line),each=2),
			  	    levels(sp.targets.5$treatment),sep="."),
	col=colors[rep(as.factor(levels(sp.targets.5$line)),each=2)], 
	cex=0.8,pch=points, ncol=2)
text(8,2.5,"batch3",col="orange")
dev.off()

## Mel 5 hours
points <- c(15,17)
colors <- c("blue", "darkgreen", "red", "orange")
pdf("figures/MDS_mel_5h.pdf",height=6,width=6)
plotMDS(mel.dge.5, col=colors[mel.targets.5$line], 
	pch=points[mel.targets.5$treatment],
	main= "MDS D. melanogaster 5h",cex=1.25)
legend("topright", legend=paste(rep(levels(mel.targets.5$line),each=2),
			  	    levels(mel.targets.5$treatment),sep="."),
	col=colors[rep(as.factor(levels(mel.targets.5$line)),each=2)], 
	pch=points, ncol=2, cex=0.5)
dev.off()
#labels=sp.targets$sample)

## Mel 50 hours
points <- c(15,17)
colors <- c("blue", "darkgreen", "red", "orange")
pdf("figures/MDS_mel_50h.pdf",height=6,width=6)
plotMDS(mel.dge.50, col=colors[mel.targets.50$line], 
	pch=points[mel.targets.50$treatment],
	main= "MDS D. melanogaster 50h",cex=1.25)
legend("topright", legend=paste(rep(levels(mel.targets.50$line),each=2),
			  	    levels(mel.targets.50$treatment),sep="."),
	col=colors[rep(as.factor(levels(mel.targets.50$line)),each=2)], pch=points, ncol=2, cex=0.8)
dev.off()
###################################
########### MA plots ##############
###################################
pdf("figures/MAplot_mel_5h.pdf", paper="a4", height=10)
plotSmear(mel.dge.5, de.tags=,main="Line (CS) effect Unparasitized")
dev.off()

pdf("figures/MAplots_50h.pdf", paper="a4", height=10)
par(mfrow=c(3,1))
plotSmear(lrt.tagwise.50.CvsS_ctl, de.tags=detags.tagwise.50.CvsS_ctl, main="Line (CS) effect Unparasitized")
plotSmear(lrt.tagwise.50.Par, de.tags=detags.tagwise.50.Par, main="Effect Parasitized")
plotSmear(lrt.tagwise.50.CSxTrt, de.tags=detags.tagwise.50.CSxTrt, main="Interaction Term: CS x Trt")
dev.off()

###################################
########### CPM plots #############
###################################
#### Species-specific
show.DEgenes.sp_specific<-function(species=sp,time=t){ 
	sp=substr(species,0,3)
	tdata=get(paste(sp,"cpm",time,sep="."))
	pTable=get(paste("LRT.pTable",sp,time,sep="."))
	cov=get(paste(sp,"targets",time,sep="."))
	o<-match(pTable[rev(order(pTable$logFC)),]$FB_sp,rownames(tdata))
	DF.ctl <- tdata[o,as.character(cov[which(cov$treatment=="ctl"),]$sample)]
	DF.par <- tdata[o,as.character(cov[which(cov$treatment=="par"),]$sample)]
	mean.ctl<-apply(DF.ctl,1,mean)
	mean.par<-apply(DF.par,1,mean)
	min.ctl<-apply(DF.ctl,1,min)
	min.par<-apply(DF.par,1,min)
	max.ctl<-apply(DF.ctl,1,max)
	max.par<-apply(DF.par,1,max)
	plot(1:nrow(DF.ctl), log2(mean.ctl + 1),
	type="p", pch=17, col="blue", lty=3, ylim=c(0,log2(max(tdata))),
	ylab="Log2(Counts Per Million)",
     	xlab="",
     	xaxt="n",
     	main=paste("D",species,sep="."))
	points(1:nrow(DF.par), log2(mean.par),
     	type="p", pch=19, col="red", lty=3)
	len = .07
	for (i in 1:nrow(DF.ctl)) {
     		arrows(i, log2(mean.ctl[i]+1),
            		i, log2(max.ctl[i]+1)+0.1,
            		angle=90, length=len, col="blue")
     		arrows(i, log2(mean.ctl[i]+1),
            		i, log2(min.ctl[i]+1)-0.1,
            		angle=90, length=len, col="blue")
     		arrows(i, log2(mean.par[i]+1),
            		i, log2(max.par[i]+1)+0.1,
            		angle=90, length=len, col="red")
	     	arrows(i, log2(mean.par[i]+1),
            		i, log2(min.par[i]+1)-0.1,
		        angle=90, length=len, col="red")
}
text(seq(1,nrow(DF.ctl)), par("usr")[3]-0.25,
	srt = 60, adj= 1, xpd = TRUE,
	labels = paste(pTable[which(pTable$FB_sp%in%rownames(DF.par)),]$CG), cex=1) 
	#labels=rownames(DF.par) for sp_specific FB
}

pdf("figures/CPMplots_LRTspspecific_5h.pdf", paper="a4")
par(mfrow=c(3,1), mar=c(5, 4, 4, 0.5))
show.DEgenes.sp_specific(species="simulans",time=5)
show.DEgenes.sp_specific(species="sechellia",time=5)
show.DEgenes.sp_specific(species="yakuba",time=5)
dev.off()

pdf("figures/CPMplots_LRTspspecific_50h.pdf", paper="a4")
par(mfrow=c(0,1), mar=c(5, 4, 4, 0.5))
show.DEgenes.sp_specific(species="simulans",time=50)	
#show.DEgenes.sp_specific(species="sechellia",time=50)
#show.DEgenes.sp_specific(species="yakuba",time=50)
dev.off()

### All species + Dmel
show.DEgenes<-function(FBid=NULL, tdata=NULL, cov=NULL, cg=NULL, q=NULL){#tdata=sp.dge.keep$counts,cov=sp.targets, q=NULL){ 
    gene<-which(rownames(tdata)==FBid) 
    y<-as.numeric(tdata[gene,])
    y<-log2(y+1) #for log2 count data
    interaction.plot(droplevels(cov$line),droplevels(cov$treatment),y,
    ylim=c(min(y),max(y)),
    ylab="Log2(Counts Per Million)",
    xlab="Line",
    legend=FALSE,
    col=c("blue","red"))
    points(rep(1:4, each=6)[cov$treatment=="par"],y[cov$treatment=="par"],pch=19,col="red")
    points(rep(1:4, each=6)[cov$treatment=="ctl"],y[cov$treatment=="ctl"],pch=17,col="blue") 
    #points(rep(1:4, each=5)[cov$treatment=="ctl"],y[cov$treatment=="ctl"],pch=17,col="blue")#comment this out for yak.50h
    #comment this out when completely excluding yak
#    points(rep(1:3, each=6)[cov$treatment=="ctl"],y[cov$treatment=="ctl"],pch=17,col="blue")
#    points(rep(1:3, each=6)[cov$treatment=="par"],y[cov$treatment=="par"],pch=19,col="red")    
    #label the replicates
#    text(rep(1:4,1),y[cov$treatment=="par" & cov$SBatch==1]+0.1,"b1")
#    text(rep(1:4,1),y[cov$treatment=="par" & cov$Batch==2]+0.1,"b2")
#    text(rep(1:4,1),y[cov$treatment=="par" & cov$Batch==3]+0.1,"b3")
#    text(rep(1:4,1),y[cov$treatment=="ctl" & cov$Batch==3]+0.1,"b3")
#    text(rep(1:4,1),y[cov$treatment=="ctl" & cov$Batch==2]+0.1,"b2")
#    text(rep(1:4,1),y[cov$treatment=="ctl" & cov$Batch==1]+0.1,"b1")
#    text(rep(1:4,1),y[cov$treatment=="ctl" & cov$Batch==1]+0.1,"b1")
#    abline(v=2.5, lty=2, col="grey")
    title(main=paste(cg,"(", FBid, ")", sep=" "), sub=paste("FDR=", signif(q , digits=4)))
}

pdf("figures/CPMplots_LRTmel_5h.pdf", paper="a4", height=10)
par(mfrow=c(4,3))
for (j in LRT.pTable.mel.Par.5[order(LRT.pTable.mel.Par.5$FDR),]$FB_mel){
 show.DEgenes(j,mel.cpm.5,mel.targets.5,
 LRT.pTable.mel.Par.5[which(LRT.pTable.mel.Par.5$FB_mel==j),]$CG,
 LRT.pTable.mel.Par.5[which(LRT.pTable.mel.Par.5$FB_mel==j),]$FDR)
 }
dev.off()


pdf("figures/CPMplots_LRTmel_50h.pdf", paper="a4", height=10)
par(mfrow=c(4,3))
for (j in LRT.pTable.mel.Par.50[order(LRT.pTable.mel.Par.50$FDR),]$FB_mel){
 show.DEgenes(j,mel.cpm.50,mel.targets.50,
 LRT.pTable.mel.Par.50[which(LRT.pTable.mel.Par.50$FB_mel==j),]$CG,
 LRT.pTable.mel.Par.50[which(LRT.pTable.mel.Par.50$FB_mel==j),]$FDR)
 }
dev.off()

#This only shows spALL
pdf("figures/CPMplots_LRTspAll_5h.pdf", paper="a4", height=10)
par(mfrow=c(4,3))
for (j in LRT.pTable.sp.Par.5[order(LRT.pTable.sp.Par.5$FDR),]$FB_mel){
 show.DEgenes(j,sp.cpm.5,sp.targets.5,
 LRT.pTable.sp.Par.5[which(LRT.pTable.sp.Par.5$FB_mel==j),]$CG,
 LRT.pTable.sp.Par.5[which(LRT.pTable.sp.Par.5$FB_mel==j),]$FDR)
 }
dev.off()

# Excluding yak at 50 hours
pdf("figures/CPMplots_LRTspAll_50h.pdf", paper="a4", height=10)
par(mfrow=c(3,3))
for (j in LRT.pTable.sp.Par.50[order(LRT.pTable.sp.Par.50$FDR),]$FB_mel){
 show.DEgenes(j,
 sp.cpm.50,
 sp.targets.50,
 LRT.pTable.sp.Par.50[which(LRT.pTable.sp.Par.50$FB_mel==j),]$CG,
 LRT.pTable.sp.Par.50[which(LRT.pTable.sp.Par.50$FB_mel==j),]$FDR)
 }
dev.off()


#This binds spAll and spClade
pdf("figures/CPMplots_LRTsp_5h.pdf", paper="a4", height=10)
LRT.pTable.sp.All.5<-
rbind(LRT.pTable.sp.Clade.5,
LRT.pTable.sp.Par.5[!(LRT.pTable.sp.Par.5$FBgn_mel%in%LRT.pTable.sp.Clade.5$FBgn_mel),])
par(mfrow=c(4,3))
for (j in LRT.pTable.sp.All.5[order(LRT.pTable.sp.All.5$FDR),]$FBgn){
 show.DEgenes(j,sp.cpm.5,sp.targets.5,
 LRT.pTable.sp.All.5[which(LRT.pTable.sp.All.5$FBgn_mel==j),]$CG,
 LRT.pTable.sp.All.5[which(LRT.pTable.sp.All.5$FBgn_mel==j),]$FDR)
 }
dev.off()

pdf("figures/CPMplots_LRTsp_50h.pdf", paper="a4", height=10)
LRT.pTable.sp.All.50<-
rbind(LRT.pTable.sp.Clade.50,
LRT.pTable.sp.Par.50[!(LRT.pTable.sp.Par.50$FBgn_mel%in%LRT.pTable.sp.Clade.50$FBgn_mel),])
par(mfrow=c(4,3))
for (j in LRT.pTable.sp.All.50[order(LRT.pTable.sp.All.50$FDR),]$FBgn){
 show.DEgenes(j,sp.cpm.50,sp.targets.50,
 LRT.pTable.sp.All.50[which(LRT.pTable.sp.All.50$FBgn_mel==j),]$CG,
 LRT.pTable.sp.All.50[which(LRT.pTable.sp.All.50$FBgn_mel==j),]$FDR)
 }
dev.off()


######## Adhoc Analyses #########

### D. mel CPM of CvsS in Wertheim
CvsS_wertheim<-read.table("other_studies/DE_66h_FBid_symbol.csv",header=T)
 
pdf("figures/CPMplots_adhocCvsS_50h.pdf", paper="a4", height=10)
par(mfrow=c(4,3))
for (j in CvsS_wertheim$FBID_KEY){
 show.DEgenes(j,mel.cpm.50,mel.targets.50,
 CvsS_wertheim[which(CvsS_wertheim$FBID_KEY==j),]$ANNOTATION_SYMBOL,
 1.000) }
dev.off()

pdf("figures/CPMplots_adhocCvsS_5h.pdf", paper="a4", height=10)
par(mfrow=c(4,3))
for (j in CvsS_wertheim$FBID_KEY){
 show.DEgenes(j,mel.cpm.5,mel.targets.5,
 CvsS_wertheim[which(CvsS_wertheim$FBID_KEY==j),]$ANNOTATION_SYMBOL,
 1.000) }
dev.off()

### D. yakuba CPM of D.mel DE

FBgn_mel_only.5 <-as.character(
		pTable.mel.Par.5[-which(pTable.mel.Par.5$FBgn_mel%in%
		pTable.sp.Par.5$FBgn_mel),]$FBgn)

FBgn_mel_sp.5 <-as.character(
		FBgn_mel_only.5[which(FBgn_mel_only.5%in%
		rownames(sp.cpm.5))])

FBgn_mel_only.50 <-as.character(
		pTable.mel.Par.50$FBgn_mel)
#		pTable.mel.Par.50[-which(pTable.mel.Par.50$FBgn_mel%in%
#		pTable.sp.Par.50$FBgn_mel),]$FBgn)


FBgn_mel_sp.50 <-as.character(
		FBgn_mel_only.50[which(FBgn_mel_only.50%in%
		rownames(sp.cpm.5))])

pdf("figures/CPMplots_adhocYakMel5h.pdf", paper="a4", height=10)
par(mfrow=c(4,3))
for (j
 in FBgn_mel_sp.5){
 show.DEgenes(j,sp.cpm.5,sp.targets.5,
 pTable.mel.Par.5[which(pTable.mel.Par.5$FBgn_mel==j),]$CG,
 1.000) }
dev.off()

pdf("figures/CPMplots_adhocYakMel50h.pdf", paper="a4", height=10)
par(mfrow=c(4,3))
for (j in FBgn_mel_sp.50){
 show.DEgenes(j,sp.cpm.50,sp.targets.50,
 pTable.mel.Par.50[which(pTable.mel.Par.50$FBgn_mel==j),]$CG,
 1.000) }
dev.off()


pdf("figures/CPMplots_melPar_50h.pdf", paper="a4", height=10)
par(mfrow=c(4,3))
for (j in pTable.mel.Par.50[order(pTable.mel.Par.50$FDR),]$FBgn){
 show.DEgenes(j,mel.cpm.50,mel.targets.50,
 pTable.mel.Par.50[which(pTable.mel.Par.50$FBgn_mel==j),]$CG,
 pTable.mel.Par.50[which(pTable.mel.Par.50$FBgn_mel==j),]$FDR)
 }
dev.off()


pdf("figures/CPMplots_Par_5h.pdf", paper="a4", height=10)
par(mfrow=c(4,3))
for(i in LRT.pTable.sp.Clade.5[order(LRT.pTable.sp.Clade.5$FDR),]$FBgn_mel){ 
  show.DEgenes(FBid=i,group="sp", time=50, model="Par")
  }
dev.off()


pdf("figures/CPMplots_Par_50h.pdf", paper="a4", height=10)
par(mfrow=c(4,3))
for(i in detags.tagwise.50.Par[order(pTable.tagwise.50.Par$FDR)]){
  show.DEgenes.50(FBid=i, q=FDR.tagwise.50.Par)
  }
dev.off()

pdf("figures/CPMplots_CvsSctl_50h.pdf", paper="a4", height=10)
par(mfrow=c(4,3))
for(i in detags.tagwise.50.CvsS_ctl[order(pTable.tagwise.50.CvsS_ctl$FDR)]){
  show.DEgenes.50(FBid=i, q=FDR.tagwise.50.CvsS_ctl)
  }
dev.off()

pdf("figures/CPMplots_CSxTrt_50h.pdf", paper="a4", height=10)
par(mfrow=c(4,3))
for(i in detags.tagwise.50.CSxTrt[order(pTable.tagwise.50.CSxTrt$FDR)]){
  show.DEgenes.50(FBid=i, q=FDR.tagwise.50.CSxTrt)
  }
dev.off()


###################################
########### Heatmaps ##############
###################################
### Spp 
#5h
library(gplots)
matrix.sp.5<-sp.cpm.5[which(rownames(sp.cpm.5)%in%
LRT.pTable.sp.All.5$FBgn_mel),]

#replace FB for CG and sample names
m<-match(rownames(matrix.sp.5),LRT.pTable.sp.All.5$FBgn_mel)
s<-match(colnames(matrix.sp.5),sp.targets.5$sample)
rownames(matrix.sp.5)<-LRT.pTable.sp.All.5[m,]$CG
colnames(matrix.sp.5)<-
paste(sp.targets.5[s,]$line,sp.targets.5[s,]$treatment,sep=".")
#Ward's minimum variance criterion minimizes the total within-cluster variance. To implement this method, at each step find the pair of clusters that leads to minimum increase in total within-cluster variance after merging. This increase is a weighted squared distance between cluster centers. At the initial step, all clusters are singletons 
hr<-hclust(as.dist(1-cor(t(as.matrix(matrix.sp.5)), method = "pearson")))
hc<-hclust(as.dist(1-cor(as.matrix(matrix.sp.5), method = "pearson")),method="ward.D2" )

png("figures/heatmap_LRTallspp_5.png")
heatmap.2(log(matrix.sp.5+1),col=topo.colors(75), 
	  scale="none",key=TRUE, symkey=FALSE, 
	  Rowv = as.dendrogram(hr), Colv = as.dendrogram(hc),
	  density.info="none", trace="none", cexRow=0.8)
dev.off()

#50h

matrix.sp.50<-sp.cpm.50[which(rownames(sp.cpm.50)%in%
LRT.pTable.sp.All.50$FBgn_mel),]

#replace FB for CG and sample names
m<-match(rownames(matrix.sp.50),LRT.pTable.sp.All.50$FBgn_mel)
s<-match(colnames(matrix.sp.50),sp.targets.50$sample)
rownames(matrix.sp.50)<-LRT.pTable.sp.All.50[m,]$CG
colnames(matrix.sp.50)<-
paste(sp.targets.50[s,]$line,sp.targets.50[s,]$treatment,sep=".")
#Ward's minimum variance criterion minimizes the total within-cluster variance. To implement this method, at each step find the pair of clusters that leads to minimum increase in total within-cluster variance after merging. This increase is a weighted squared distance between cluster centers. At the initial step, all clusters are singletons 
hr<-hclust(as.dist(1-cor(t(as.matrix(matrix.sp.50)), method = "pearson")))
hc<-hclust(as.dist(1-cor(as.matrix(matrix.sp.50), method = "pearson")),method="ward.D2" )

png("figures/heatmap_LRTallspp_50.png")
heatmap.2(log(matrix.sp.50+1),col=topo.colors(75), 
	  scale="none",key=TRUE, symkey=FALSE, 
	  Rowv = as.dendrogram(hr), Colv = as.dendrogram(hc),
	  density.info="none", trace="none", cexRow=0.8)
dev.off()
### Mel
#5 hours
matrix.mel.5<-mel.cpm.5[which(rownames(mel.cpm.5)%in%
LRT.pTable.mel.Par.5$FB_mel),]

m<-match(rownames(matrix.mel.5),LRT.pTable.mel.Par.5$FB_mel)
s<-match(colnames(matrix.mel.5),mel.targets.5$sample)
rownames(matrix.mel.5)[which(!is.na(LRT.pTable.mel.Par.5[m,]$CG))]<-
as.character(LRT.pTable.mel.Par.5[which(!is.na(LRT.pTable.mel.Par.5[m,]$CG)),]$CG) #replace for CG except when it is NA
colnames(matrix.mel.5)<-
paste(mel.targets.5[s,]$line,mel.targets.5[s,]$treatment,mel.targets.5[s,]$Batch,sep=".")
hr<-hclust(as.dist(1-cor(t(as.matrix(matrix.mel.5)), method = "pearson")))
hc<-hclust(as.dist(1-cor(as.matrix(matrix.mel.5), method = "pearson")),method="ward.D2" )
heatmap.2(log(matrix.mel.5+1),col=topo.colors(75), 
	  scale="none",key=TRUE, symkey=FALSE, 
	  Rowv = as.dendrogram(hr), Colv = as.dendrogram(hc),
	  cexRow = 0.5, density.info="none", trace="none")

#50 hours
matrix.mel.50<-mel.cpm.50[which(rownames(mel.cpm.50)%in%
LRT.pTable.mel.Par.50$FB_mel),]
rownames(matrix.mel.50)[which(!is.na(LRT.pTable.mel.Par.50[m,]$CG))]<-
as.character(LRT.pTable.mel.Par.50[which(!is.na(LRT.pTable.mel.Par.50[m,]$CG)),]$CG)
m<-match(rownames(matrix.mel.50),LRT.pTable.mel.Par.50$FB_mel)
s<-match(colnames(matrix.mel.50),mel.targets.50$sample)

colnames(matrix.mel.50)<-
paste(mel.targets.50[s,]$line,mel.targets.50[s,]$treatment,mel.targets.50[s,]$Batch,sep=".")
hr<-hclust(as.dist(1-cor(t(as.matrix(matrix.mel.50)), method = "pearson")))
hc<-hclust(as.dist(1-cor(as.matrix(matrix.mel.50), method = "pearson")),method="ward.D2" )
heatmap.2(log(matrix.mel.50+1),col=topo.colors(75), 
	  scale="none",key=TRUE, symkey=FALSE, 
	  Rowv = as.dendrogram(hr), Colv = as.dendrogram(hc),
	  cexRow = 0.5, density.info="none", trace="none")

### All spp and time points
#Take a list of all significant FBgn that have an ortholog in mel 
#(does it make sense to also include those that do not have orthologs in Dmel?)
allFB<-c(as.character(LRT.pTable.mel.Par.5[,1]),
	    as.character(LRT.pTable.mel.Par.50[,1]),
	    as.character(LRT.pTable.sp.Par.5[!is.na(LRT.pTable.sp.Par.5$FB_mel),]$FB_mel), #has orth in mel
	    as.character(LRT.pTable.sp.Par.50[!is.na(LRT.pTable.sp.Par.50$FB_mel),]$FB_mel), 
	    as.character(LRT.pTable.sim.5[!is.na(LRT.pTable.sim.5$FB_mel),]$FB_mel),
	    as.character(LRT.pTable.sim.50[!is.na(LRT.pTable.sim.50$FB_mel),]$FB_mel),
	    as.character(LRT.pTable.sec.5[!is.na(LRT.pTable.sec.5$FB_mel),]$FB_mel),
	    as.character(LRT.pTable.sec.50[!is.na(LRT.pTable.sec.50$FB_mel),]$FB_mel),
	    as.character(LRT.pTable.yak.5[!is.na(LRT.pTable.yak.5$FB_mel),]$FB_mel),
	    as.character(LRT.pTable.yak.50[!is.na(LRT.pTable.yak.50$FB_mel),]$FB_mel))


# Filter duplicates
FB<-unique(allFB)
# Initiate a matrix with FB as rows and line.time as columns
targets<-read.table("../targets.txt",header=T)
targets<-targets[-which(targets$line=="yak" & targets$time==50),] #remove yak at 50
targets$line.time<-as.factor(paste(targets$line,targets$time,sep="."))
cols = levels(targets$line.time)

FoldChange<-matrix(ncol = length(cols), nrow = length(FB))
colnames(FoldChange)<-cols
rownames(FoldChange)<-FB

for (i in rownames(FoldChange)) {
	for (j in colnames(FoldChange)){
		samples_ctl <- targets[which(targets$line.time==j & 	
				targets$treatment=="ctl"),]$sample
		samples_par <- targets[which(targets$line.time==j & 	
				targets$treatment=="par"),]$sample
		sp <- unlist(strsplit(j,"[.]"))[1]
		time<-unlist(strsplit(j,"[.]"))[2]
#		i_sp = 'NA'
		if (length(grep("mel",j))!=0){
			cpm <- get(paste("mel.cpm",time,sep="."))
			cpm_ctl <- as.numeric(cpm[i,as.character(samples_ctl)])
			cpm_par <- as.numeric(cpm[i,as.character(samples_par)])
			median_ctl <- median(cpm_ctl)
			median_par <- median(cpm_par)
			FoldChange[i,j] <- as.numeric(log2((median_par+0.001)/(median_ctl+0.001)))
			} else if(length(grep("mel",j))==0){ 
				cpm <- get(paste("sp.cpm",time,sep="."))
#				x<-which(i==orthologs$V1)
#				ortho_subset<-orthologs[x, ]
#				cpm <- get(paste(sp,"cpm",time,sep="."))
				# redefine i with the sp-specific FB. 
				# Some genes have 2 IDs in mel, I'm taking only the first value (?)
#				i_sp <-as.character(ortho_subset[grep(sp,orthologs[x, ]$V7),]$V6)[1]
#				if (length(which(rownames(cpm)==i_sp))!=0){
				if (length(which(rownames(cpm)==i))!=0){
					cpm_ctl <- as.numeric(cpm[i,as.character(samples_ctl)]) #i_sp
					cpm_par <- as.numeric(cpm[i,as.character(samples_par)])
					median_ctl <- median(cpm_ctl)
					median_par <- median(cpm_par)
					#avoid dividing by 0 -> this can have big effect in the magnitud of FC 
					FoldChange[i,j] <- as.numeric(log2((median_par+0.001)/(median_ctl+0.001)))
					}	 
				}
		}
}

#Remove rows that have NA -> genes do not exist or do not expressed for that line
FoldChange<-
FoldChange[which(!rownames(FoldChange)%in%names(which(is.na(rowSums(FoldChange))))),]

fc<-rownames(FoldChange)
for(r in c(1:length(fc))){
	if(!is.na(as.character(orthologs[orthologs$V1==fc[r],]$V2[1]))){
		fc[r]<-as.character(orthologs[orthologs$V1==fc[r],]$V2[1])
	}
}

rownames(FoldChange)<-fc

annotation<-read.table("/home/laura/PhDpapers/RNAseq/analysis/metanalysis/LRT_annotation.csv",header=T)
m<-match(rownames(FoldChange),annotation$Gene)
annorder<-annotation[m,]
ann.dat<-data.frame(humoral=annorder$Humoral,cellular=annorder$Cellular,
		   signalling=annorder$Signalling,receptor=annorder$Receptor,
		   melanization=annorder$Melanization,stress=annorder$Stress,
		   proteolysis=annorder$Proteolysis,morphogenesis=annorder$Morphogenesis, 
		   Orthology=annorder$Orthology)
#species order for heatmap
sp.order<-c("melC1.5","melC2.5","melS1.5","melS2.5","sim.5","sec.5","yak.5",
"melC1.50","melC2.50","melS1.50","melS2.50","sim.50","sec.50")

#Define clusters
library(vegan)
hr<-hclust(as.dist(1-cor(t(as.matrix(FoldChange)), method = "pearson")),method="ward.D2" )

#Include phenotypes
resistance<-data.frame(resistance=c(20,20,50,50,0,80,90,20,20,50,50,0,80,90))

library("Heatplus")
library(RColorBrewer)

svg("figures/HeatmapFC_allspp_times_annotation_glmLRT_pearson.svg")
plot(annHeatmap2(FoldChange[,sp.order],
	col = colorRampPalette(c("black","darkblue", "blue","gold","yellow", "lightgoldenrodyellow"), 
	space = "rgb")(12),
	breaks = -6:6,#niceBreaks(c(-3,3),12), #-6:6,
	scale= "none",#"row", 
	dendrogram = list(Row = list(dendro = as.dendrogram(hr)), Col = list(status="no")), #list(dendro = as.dendrogram(hr))),
	legend = 3, #3
	labels = list(Row = list(ncol = 10, cex=0.6), Col = list(cex=0.7)),
	ann = list(Row = list(data = ann.dat), Col = NULL), #list(data=resistance))
	cluster = list(Row = list(cuth = 2.5), col = brewer.pal(5, "Set2"))), # cuth gives the height at which the dendrogram should be cut to form clusters, and col specifies the colours for the clusters
	 widths=c(1,6,1), heights=c(0.8,6,1))
dev.off()

#https://www.biostars.org/p/14156/
###################################
######### Venn diagrams ###########
###################################

#Venn diagrams
library(VennDiagram)

#VennDiagrams

venn.plot.5 <- venn.diagram(list(LRT.pTable.sp.MelxPar.5$FBgn_mel,
				 LRT.pTable.sp.MelSimxPar.5$FBgn_mel,
				 LRT.pTable.sp.MelSimSecxPar.5$FBgn_mel,
				 LRT.pTable.sp.Par.5$FBgn_mel), NULL, 
				 fill=c("yellow","green","blue","red"), 
				 alpha=c(0.5,0.5,0.5,0.5), cex = 2, 
				 cat.fontface=4, 
			category.names=c("mel","mel+sim","mel+sim+sec","mel+sim+sec+yak"))


png("figures/Venn_LRTspp_5h.png")
grid.draw(venn.plot.5)
dev.off()

venn.plot.50 <- venn.diagram(list(LRT.pTable.sp.MelxPar.50$FBgn_mel,
				 LRT.pTable.sp.MelSimxPar.50$FBgn_mel,
				 LRT.pTable.sp.MelSimSecxPar.50$FBgn_mel,
				 LRT.pTable.sp.Par.50$FBgn_mel), NULL, 
				 fill=c("yellow","green","blue","red"), 
				 alpha=c(0.5,0.5,0.5,0.5), cex = 2, 
				 cat.fontface=4, 
			category.names=c("mel","mel+sim","mel+sim+sec","mel+sim+sec+yak"))


png("figures/Venn_LRTspp_5h.png")
grid.draw(venn.plot.50)
dev.off()

#UpSetR
# Create a matrix similar to FoldChanges with contrasts in the rows and genes in 
# the columns. Quantify presence (1) absence(0)

### All spp and time points
#Take a list of all significant FBgn with the ID ortholog in mel 
# and the original FBgn if there are not orthologs

allFB_sp<-c(as.character(LRT.pTable.mel.Par.5[,1]),
	    as.character(LRT.pTable.mel.Par.50[,1]),
	    as.character(LRT.pTable.sp.Par.5[!is.na(LRT.pTable.sp.All.5$FB_mel),]$FB_mel), #has orth in mel
	    as.character(LRT.pTable.sp.Par.5[is.na(LRT.pTable.sp.All.5$FB_mel),]$FB_sp),   #does not have orth in mel
	    as.character(LRT.pTable.sp.Par.50[!is.na(LRT.pTable.sp.All.50$FB_mel),]$FB_mel), 
	    as.character(LRT.pTable.sp.Par.50[is.na(LRT.pTable.sp.All.50$FB_mel),]$FB_sp),
	    as.character(LRT.pTable.sim.5[!is.na(LRT.pTable.sim.5$FB_mel),]$FB_mel),
	    as.character(LRT.pTable.sim.5[is.na(LRT.pTable.sim.5$FB_mel),]$FB_sp),
	    as.character(LRT.pTable.sim.50[!is.na(LRT.pTable.sim.50$FB_mel),]$FB_mel),
	    as.character(LRT.pTable.sim.50[is.na(LRT.pTable.sim.50$FB_mel),]$FB_sp),
	    as.character(LRT.pTable.sec.5[!is.na(LRT.pTable.sec.5$FB_mel),]$FB_mel),
	    as.character(LRT.pTable.sec.5[is.na(LRT.pTable.sec.5$FB_mel),]$FB_sp),
	    as.character(LRT.pTable.sec.50[!is.na(LRT.pTable.sec.50$FB_mel),]$FB_mel),
	    as.character(LRT.pTable.sec.50[is.na(LRT.pTable.sec.50$FB_mel),]$FB_sp),
	    as.character(LRT.pTable.yak.5[!is.na(LRT.pTable.yak.5$FB_mel),]$FB_mel),
	    as.character(LRT.pTable.yak.5[is.na(LRT.pTable.yak.5$FB_mel),]$FB_sp),
	    as.character(LRT.pTable.yak.50[!is.na(LRT.pTable.yak.50$FB_mel),]$FB_mel),
	    as.character(LRT.pTable.yak.50[is.na(LRT.pTable.yak.50$FB_mel),]$FB_sp))

# Filter duplicates
FB_sp<-unique(allFB_sp[!is.na(allFB_sp)]) #unique(allFB_sp)#[which(!is.na(allFB_sp))])

#### 5hours ##########################
SpSpecific_mel_5<-LRT.pTable.mel.Par.5
SpSpecific_sim_5<-LRT.pTable.sim.5
SpSpecific_sec_5<-LRT.pTable.sec.5
SpSpecific_yak_5<-LRT.pTable.yak.5
sp_All_5<-LRT.pTable.sp.Par.5
allcontrasts_5<-c("SpSpecific_mel_5","SpSpecific_sim_5","SpSpecific_sec_5","SpSpecific_yak_5",
		"sp_All_5")


# Initiate a matrix with FB as columns and contrast as rows
ContrastSets_5<-matrix(nrow = length(FB_sp), ncol = length(allcontrasts_5))
rownames(ContrastSets_5)<-FB_sp
colnames(ContrastSets_5)<-allcontrasts_5

#Fill the matrix with 1 if the gene is significant in the given contrast
#and 0 otherwise. FIX: until now it only takes genes that have mel ID!!!!! 
for (i in rownames(ContrastSets_5)) {
	for (j in colnames(ContrastSets_5)){
		if(i%in%get(j)[,1]==TRUE | i%in%get(j)[,2]==TRUE){#if(i%in%get(j)[,2]==TRUE)
			ContrastSets_5[i,j]<-1
			#if(!is.na(get(j)$FB_mel[i])==FALSE){
			#	i<-get(j)$FB_mel[which(get(j)[,1]==i)]
			#	}	
		}else{ContrastSets_5[i,j]<-0}	
	}
}
library(UpSetR)

write.csv(ContrastSets_5,"ContrastSets_5.csv")
C5<-read.csv("ContrastSets_5.csv")

colnames(C5)<-c("X","Mel","Sim","Sec","Yak","All spp")

svg("figures/ContrastSets_5.svg",height=6,width=6)
upset(C5,nsets=5,
	sets.x.label="DEG",mainbar.y.label="Number of genes in category",
	name.size=8,
	sets.bar.color=c("blue","blue","orange","blue","blue"),
	point.size = 4, line.size = 1,
	main.bar.color= c("purple","purple","purple",
	"green","green","green","green","green","green","green"))
dev.off()

#### 50hours ##########################
SpSpecific_mel_50<-LRT.pTable.mel.Par.50
SpSpecific_sim_50<-LRT.pTable.sim.50
SpSpecific_sec_50<-LRT.pTable.sec.50
SpSpecific_yak_50<-LRT.pTable.yak.50
sp_All_50<-LRT.pTable.sp.All.50

allcontrasts_50<-c("SpSpecific_mel_50","SpSpecific_sim_50","SpSpecific_sec_50","sp_All_50")


# Initiate a matrix with FB as columns and contrast as rows
ContrastSets_50<-matrix(nrow = length(FB_sp), ncol = length(allcontrasts_50))
rownames(ContrastSets_50)<-FB_sp
colnames(ContrastSets_50)<-allcontrasts_50


#Fill the matrix with 1 if the gene is significant in the given contrast
#and 0 otherwise
for (i in rownames(ContrastSets_50)) {
	for (j in colnames(ContrastSets_50)){
		if(i%in%get(j)[,1]==TRUE | i%in%get(j)[,2]==TRUE){#if(i%in%get(j)[,2]==TRUE)
			ContrastSets_50[i,j]<-1
			#if(!is.na(get(j)$FB_mel[i])==FALSE){
			#	i<-get(j)$FB_mel[which(get(j)[,1]==i)]
			#	}	
		}else{ContrastSets_50[i,j]<-0}	
	}
}
write.csv(ContrastSets_50,"ContrastSets_50.csv")
C50<-read.csv("ContrastSets_50.csv")
colnames(C50)<-c("X","Mel","Sim","Sec","All spp")

svg("figures/ContrastSets_50.svg",height=6,width=6)
upset(C50,nsets=5,
	sets.x.label="DEG",mainbar.y.label="Genes shared by category",
	name.size=8,
	sets.bar.color=c("blue","orange","blue","blue"),
	point.size = 4, line.size = 1,
	main.bar.color= c("purple","purple","purple","green"))
dev.off()
