
########### MDS plots #############
## Spp 5 hours
sp.targets$group<-sp.group
points <- c(15,17)
colors <- c("blue", "darkgreen", "red", "orange")
plotMDS(sp.dge.5, col=colors[sp.targets.5$line], 
	pch=points[sp.targets.5$treatment],
	main= "MDS all species 5h",cex=1.25)
legend("topright", legend=paste(rep(levels(sp.targets.5$line),each=2),
			  	    levels(sp.targets.5$treatment),sep="."),
	col=colors[rep(as.factor(levels(sp.targets.5$line)),each=2)], pch=points, ncol=2)
#labels=sp.targets$sample)

## Spp 50 hours
plotMDS(sp.dge.50, col=colors[sp.targets.50$line], 
	pch=points[sp.targets.5$treatment],
	main= "MDS all species 50h",cex=1.25)
legend("bottomright", legend=paste(rep(levels(sp.targets.5$line),each=2),
			  	    levels(sp.targets.5$treatment),sep="."),
	col=colors[rep(as.factor(levels(sp.targets.5$line)),each=2)], pch=points, ncol=2)

## Mel 5 hours
points <- c(15,17)
colors <- c("blue", "darkgreen", "red", "orange")
plotMDS(mel.dge.5, col=colors[mel.targets.5$line], 
	pch=points[mel.targets.5$treatment],
	main= "MDS D. melanogaster 5h",cex=1.25)
legend("topright", legend=paste(rep(levels(mel.targets.5$line),each=2),
			  	    levels(mel.targets.5$treatment),sep="."),
	col=colors[rep(as.factor(levels(mel.targets.5$line)),each=2)], pch=points, ncol=2, cex=0.8)
#labels=sp.targets$sample)

## Mel 50 hours
points <- c(15,17)
colors <- c("blue", "darkgreen", "red", "orange")
plotMDS(mel.dge.50, col=colors[mel.targets.50$line], 
	pch=points[mel.targets.50$treatment],
	main= "MDS D. melanogaster 50h",cex=1.25)
legend("topright", legend=paste(rep(levels(mel.targets.50$line),each=2),
			  	    levels(mel.targets.50$treatment),sep="."),
	col=colors[rep(as.factor(levels(mel.targets.50$line)),each=2)], pch=points, ncol=2, cex=0.8)

########### MA plots #############
pdf("figures/MAplots_5h.pdf", paper="a4", height=10)
par(mfrow=c(3,1))
plotSmear(lrt.tagwise.5.CvsS_ctl, de.tags=detags.tagwise.5.CvsS_ctl, main="Line (CS) effect Unparasitized")
plotSmear(lrt.tagwise.5.Par, de.tags=detags.tagwise.5.Par, main="Effect Parasitized")
plotSmear(lrt.tagwise.5.CSxTrt, de.tags=detags.tagwise.5.CSxTrt, main="Interaction Term: CS x Trt")
dev.off()

pdf("figures/MAplots_50h.pdf", paper="a4", height=10)
par(mfrow=c(3,1))
plotSmear(lrt.tagwise.50.CvsS_ctl, de.tags=detags.tagwise.50.CvsS_ctl, main="Line (CS) effect Unparasitized")
plotSmear(lrt.tagwise.50.Par, de.tags=detags.tagwise.50.Par, main="Effect Parasitized")
plotSmear(lrt.tagwise.50.CSxTrt, de.tags=detags.tagwise.50.CSxTrt, main="Interaction Term: CS x Trt")
dev.off()


########### CPM plots #############

show.DEgenes<-function(FBid=NULL,group=NULL,time=NULL,model=NULL){ #, time=NULL, sp=NULL,,q=NULL
    tdata=get(paste(group,"dge.keep",time,sep="."))$counts
    cov=get(paste(group,"targets",sep="."))
    q=pTable$FDR#get(paste("pTable",group,model,time,sep="."))$FDR
    gene<-which(rownames(tdata)==FBid) 
    y<-as.numeric(tdata[gene,])
    y<-log2(y+1) #for log2 count data
    interaction.plot(cov$line,cov$treatment,y,
    ylim=c(min(y),max(y)),
    ylab="Log2(Counts Per Million)",
    xlab="Line",
    legend=FALSE,
    col=c("blue","red"))
    points(rep(1:4, each=6)[cov$treatment=="par"],y[cov$treatment=="par"],pch=19,col="red")
    points(rep(1:4, each=6)[cov$treatment=="ctl"],y[cov$treatment=="ctl"],pch=17,col="blue")
    abline(v=2.5, lty=2, col="grey")
    title(main=paste(FBgn_CG[which(FBgn_CG==FBid),]$CG,"(", FBid, ")", sep=" "), sub=paste("FDR=", signif(q[gene], digits=4)))
     
}


pdf("figures/CPMplots_Par_5h.pdf", paper="a4", height=10)
par(mfrow=c(4,3))
for(i in pTable.sp.Par.5[order(pTable.sp.Par.5$FDR),]$FBgn_mel[1:12]){ 
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


#########SPP###########
#plot CPM for spp
# Edit plots to include both IDs
show.DEgenes<-function(FBid=NULL, tdata=sp.dge.keep.50$counts,cov=sp.cov.50, q=NULL){ 
    gene<-which(rownames(tdata)==FBid) 
    y<-as.numeric(tdata[gene,])
    y<-log2(y+1) #for log2 count data
    interaction.plot(cov$line,cov$treatment,y,
    ylim=c(min(y),max(y)),
    ylab="Log2(Counts Per Million)",
    xlab="Line",
    legend=FALSE,
    col=c("blue","red"))
    points(rep(1:4, each=6)[cov$treatment=="par"],y[cov$treatment=="par"],pch=19,col="red")
    points(rep(1:4, each=6)[cov$treatment=="ctl"],y[cov$treatment=="ctl"],pch=17,col="blue")
    abline(v=2.5, lty=2, col="grey")
    title(main=paste(FBgn_CG[which(FBgn_CG==FBid),]$CG,"(", FBid, ")", sep=" "), sub=paste("FDR=", signif(q , digits=4)))
}


pdf("figures/CPMplots_sppPar_5h.pdf", paper="a4", height=10)
par(mfrow=c(4,3))
for(i in pTable.sp.Par.50[order(pTable.sp.Par.50$FDR),]$FBgn_mel[1:9]){ 
  show.DEgenes(FBid=i, q=pTable.sp.Par.50[which(pTable.sp.Par.50==i),]$FDR)
  }
dev.off()

########### Heatmaps #############
### Spp
library(gplots)
matrix.sp.5<-sp.cpm.5[which(rownames(sp.cpm.5)%in%
pTable.spp.5$FBgn_mel),]

#replace FB for CG and sample names
m<-match(rownames(matrix.sp.5),pTable.spp.5$FBgn)
s<-match(colnames(matrix.sp.5),sp.targets.5$sample)
rownames(matrix.sp.5)<-pTable.spp.5[m,]$CG
colnames(matrix.sp.5)<-
paste(sp.targets.5[s,]$line,sp.targets.5[s,]$treatment,sep=".")
logcpm <- cpm(y, prior.count=2, log=TRUE)
heatmap.2(log(matrix.sp.5+1),col=topo.colors(75), 
	  scale="none",key=TRUE, symkey=FALSE, 
	  density.info="none", trace="none")

### Mel
#5 hours
matrix.mel.5<-mel.cpm.5[which(rownames(mel.cpm.5)%in%
pTable.mel.Par.5$FBgn_mel),]

m<-match(rownames(matrix.mel.5),pTable.mel.Par.5$FBgn_mel)
s<-match(colnames(matrix.mel.5),mel.targets.5$sample)
rownames(matrix.mel.5)<-pTable.mel.Par.5[m,]$CG
colnames(matrix.mel.5)<-
paste(mel.targets.5[s,]$line,mel.targets.5[s,]$treatment,sep=".")
logcpm <- cpm(y, prior.count=2, log=TRUE)
heatmap.2(log(matrix.mel.5+1),col=topo.colors(75), 
	  scale="none",key=TRUE, symkey=FALSE, 
	  density.info="none", trace="none")

#05 hours
matrix.mel.50<-mel.cpm.50[which(rownames(mel.cpm.50)%in%
pTable.mel.Par.50$FBgn_mel),]

m<-match(rownames(matrix.mel.50),pTable.mel.Par.50$FBgn_mel)
s<-match(colnames(matrix.mel.50),mel.targets.50$sample)
rownames(matrix.mel.50)<-pTable.mel.Par.50[m,]$CG
colnames(matrix.mel.50)<-
paste(mel.targets.50[s,]$line,mel.targets.50[s,]$treatment,sep=".")
heatmap.2(log(matrix.mel.50+1),col=topo.colors(75), 
	  scale="none",key=TRUE, symkey=FALSE, 
	  density.info="none", trace="none")
