library(edgeR)

# set the settings for input and output
wd<-getwd()

DIR.IN<-"~/PhDpapers/RNAseq/counts/"
DIR.OUT<-"~/PhDpapers/RNAseq/analysis/spp_all/"
DIR.gff<-"~/PhDpapers/RNAseq/gff/"

# get information about samples
targets<-read.table("/home/laura/PhDpapers/RNAseq/targets.txt",header=T)
targets<-targets[-c(48),] #yak48_par50h
targets<-with(targets,targets[order(line,treatment),])

###################################
######    Get the Counts     ######
######  Create One matrix    ######
####  for all sp in each time #####
###################################

time<-c("5","50")
species = c("melC1","sim","sec", "yak")
#sp.targets.5<-droplevels(targets[which(targets$line%in%species & targets$time==5),])
#sp.targets.50<-droplevels(targets[which(targets$line%in%species & targets$time==50),])

#Read Flybase Ortholog file. Header lines removed to avoid issues with reading the table
orthologs<-read.table(paste(DIR.gff,"gene_orthologs_fb_2015_02.tsv",sep=""))

#For each species, open its directory, get the counts from the file and 
#bind it to the matrix, then change the rowname to its corresponding 
#orthologous in D. melanogaster
for(s in species){
	for(t in time){
		sp_dir=paste(DIR.IN,s,"/",sep="")
		samples<-targets[which(targets$line==s & targets$time==t),]$sample
		sp.counts<-vector()
		for(n in samples){
		assign(n,read.table(paste(sp_dir,n,".counts",sep="")))
		strip=nrow(get(n))-5
		sp.counts<-cbind(sp.counts,get(n)$V2[1:strip])
		}
		colnames(sp.counts)<-samples
		rownames(sp.counts)<-get(n)$V1[1:strip]
		if(any(grep("mel",s))==FALSE){
#			m<-match(rownames(sp.counts),orthologs$V1)
#			rownames(sp.counts)<-orthologs[m,]$V1
#			assign(paste(s,t,sep="."),sp.counts)
#			} else if(any(grep("mel",s))==FALSE){
			m<-match(rownames(sp.counts),orthologs$V6)
			rownames(sp.counts)<-orthologs[m,]$V1 	#replace Id for mel orth id (V1 or V2)
			}
		assign(paste(s,t,sep="."),sp.counts)
	}
}

# we have to subset the intersections and remove NA
intersect.5<-intersect(intersect(rownames(sec.5),rownames(melC1.5)),intersect(rownames(sim.5),rownames(yak.5)))
na<-which(!is.na(intersect.5))
intersect.5<-intersect.5[na]
intersect.50<-intersect(intersect(rownames(sec.50),rownames(melC1.50)),intersect(rownames(sim.50),rownames(yak.50)))
na<-which(!is.na(intersect.50))
intersect.50<-intersect.50[na]

#put all counts together
mC1.5<-match(intersect.5,rownames(melC1.5))
msec.5<-match(intersect.5,rownames(sec.5))
msim.5<-match(intersect.5,rownames(sim.5))
myak.5<-match(intersect.5,rownames(yak.5))

mC1.50<-match(intersect.50,rownames(melC1.50))
msec.50<-match(intersect.50,rownames(sec.50))
msim.50<-match(intersect.50,rownames(sim.50))
myak.50<-match(intersect.50,rownames(yak.50))

sp.counts.5<-cbind(melC1.5[mC1.5,],sec.5[msec.5,],sim.5[msim.5,],yak.5[myak.5,])
sp.counts.50<-cbind(melC1.50[mC1.50,],sec.50[msec.50,],sim.50[msim.50,],yak.50[myak.50,])

#need to order the columns according to the targets
sp.counts.5<-sp.counts.5[,as.character(targets[which(targets$line%in%species & targets$time==5),]$sample)] 
sp.counts.50<-sp.counts.50[,as.character(targets[which(targets$line%in%species & targets$time==50),]$sample)]

## collate all info
#sp.dge.5<-DGEList(counts=sp.counts.5, group=paste(sp.targets.5$line, sp.targets.5$treatment, sep="."), genes=rownames(sp.counts.5))
#sp.dge.50<-DGEList(counts=sp.counts.50, group=paste(sp.targets.50$line, sp.targets.50$treatment, sep="."), genes=rownames(sp.counts.50))

#filter out very low counts
#keep.5<-rowSums(cpm(sp.dge.5)>5) >= 3
#keep.50<-rowSums(cpm(sp.dge.50)>5) >= 3

#sp.dge.keep.5<-sp.dge.5[keep.5,]
#sp.dge.keep.5$samples$lib.size <- colSums(sp.dge.keep.5$counts) 

#sp.dge.keep.50<-sp.dge.50[keep.50,]
#sp.dge.keep.50$samples$lib.size <- colSums(sp.dge.keep.50$counts) 

# normalize using TTM
#sp.dge.keep.5<-calcNormFactors(sp.dge.keep.5)
#sp.dge.keep.50<-calcNormFactors(sp.dge.keep.50)
############################################################################
# start statistical analysis of differential expression
#using contrasts to compare differences in expression for all species
# 1) Par = Par vs Ctl for all species
# 2) SpxTrt = the interaction term: Differences in response to Par between spp. I made each possible contrast manually to explore the data.

for (t in time){
	#Get desgin matrix and contrasts
	sp.counts<-get(paste("sp.counts",t,sep="."))
	sp.targets<-droplevels(targets[which(targets$line%in%species & targets$time==t),])
	#times have different replicates (missing one yakuba)
	if(t==5){
		sp.targets$Batch<- factor(rep(1:3, dim(sp.targets)[1]/3))	
		sp.group<-paste(sp.targets$line, sp.targets$treatment, sep=".")
		sp.design<-model.matrix(~0+sp.group+Batch,data=sp.targets)
		colnames(sp.design)<-c("mel.ctl", "mel.par", "sec.ctl", "sec.par", 
			"sim.ctl", "sim.par", "yak.ctl", "yak.par","Batch2", "Batch3")
		sp.contrasts<-makeContrasts(Par=(mel.par+sec.par+sim.par+yak.par)/4-(mel.ctl+sec.ctl+sim.ctl+yak.ctl)/4, 
			MelxPar=(mel.par-mel.ctl), SecxPar=(sec.par-sec.ctl), SimxPar=(sim.par-sim.ctl), 
			YakxPar=(yak.par-yak.ctl), ResSp=(mel.par+sim.par+yak.par)/3-(mel.ctl+sim.ctl+yak.ctl)/3,
			MelSimxPar=(mel.par+sim.par)/2-(mel.ctl+sim.ctl)/2,
			MelSimSecxPar=(mel.par+sec.par+sim.par)/3-(mel.ctl+sec.ctl+sim.ctl)/3,
			levels=sp.design)
	}else if (t==50){
		#exclude yakuba for 50 because its biological replicates are to high and one is missing
		sp.targets$Batch <- factor(c(rep(1:3,7),c(1,3)))
		#sp.targets<-sp.targets[-which(sp.targets$line%in%"yak"),]
		sp.counts<-sp.counts[,which(sp.targets$sample%in%colnames(sp.counts))]
		sp.group<-paste(sp.targets$line, sp.targets$treatment, sep=".")
		sp.design<-model.matrix(~0+sp.group+Batch,data=sp.targets)
		colnames(sp.design)<-c("mel.ctl", "mel.par", "sec.ctl", "sec.par", 
			"sim.ctl", "sim.par", "yak.ctl", "yak.par","Batch2", "Batch3")
		sp.contrasts<-makeContrasts(Par=(mel.par+sec.par+sim.par+yak.par)/4-(mel.ctl+sec.ctl+sim.ctl+yak.ctl)/4, 
			MelxPar=(mel.par-mel.ctl), SecxPar=(sec.par-sec.ctl), SimxPar=(sim.par-sim.ctl), 
			YakxPar=(yak.par-yak.ctl), ResSp=(mel.par+sim.par+yak.par)/3-(mel.ctl+sim.ctl+yak.ctl)/3,
			MelSimxPar=(mel.par+sim.par)/2-(mel.ctl+sim.ctl)/2,
			MelSimSecxPar=(mel.par+sec.par+sim.par)/3-(mel.ctl+sec.ctl+sim.ctl)/3,
			levels=sp.design)
		#colnames(sp.design)<-c("mel.ctl", "mel.par", "sec.ctl", "sec.par", 
		#	"sim.ctl", "sim.par", "Batch2", "Batch3")
		#sp.contrasts<-makeContrasts(Par=(mel.par+sec.par+sim.par)/3-(mel.ctl+sec.ctl+sim.ctl)/3, 
		#	MelxPar=(mel.par-mel.ctl), SecxPar=(sec.par-sec.ctl), SimxPar=(sim.par-sim.ctl), 
		#	MelSimxPar=(mel.par+sim.par)/2-(mel.ctl+sim.ctl)/2,
		#	levels=sp.design)
		
	}
	#m<-sp.targets[which(sp.targets$line=="melC1"),]$sample
	#mel.mean<-apply(sp.counts.5[,which(colnames(sp.counts.5)%in%m)],1,function(x) mean(x))
	#mel.sd<-apply(sp.counts.5[,which(colnames(sp.counts.5)%in%m)],1,function(x) sd(x))
	#se<-sp.targets[which(sp.targets$line=="sec"),]$sample
	#sec.sd<-apply(sp.counts.5[,which(colnames(sp.counts.5)%in%se)],1,function(x) sd(x))
	#si<-sp.targets[which(sp.targets$line=="sim"),]$sample
	#sim.sd<-apply(sp.counts.5[,which(colnames(sp.counts.5)%in%si)],1,function(x) sd(x))
	#y<-sp.targets[which(sp.targets$line=="yak"),]$sample
	#yak.sd<-apply(sp.counts.5[,which(colnames(sp.counts.5)%in%y)],1,function(x) sd(x))

	## collate all info
	sp.dge<-DGEList(counts=sp.counts, group=sp.group, genes=rownames(sp.counts))

	# filter out the tags that are not expressed at >5 CPM in at least 3 libraries
	keep<-rowSums(cpm(sp.dge)>5) >= 3
	sp.dge.keep<-sp.dge[keep,]
	sp.dge.keep$samples$lib.size <- colSums(sp.dge.keep$counts)

	# normalize using TTM
	sp.dge.keep<-calcNormFactors(sp.dge.keep)
	assign(paste("sp.dge",t,sep="."),sp.dge.keep)

	# sp.targets$group<-sp.group
	# points <- c(15,16,17,20)
	# colors <- rep(c("blue", "darkgreen", "red", "yellow"), 2)
	# color <- c("blue","yellow")
	# plotMDS(sp.dge.keep, col=colors[sp.targets$line], pch=points[sp.targets$line],main= "MDS all species 50h",
	# labels=sp.targets$sample,main= "MDS all species 50h")
	# plotMDS(sp.dge.keep, col=colors[sp.targets$line], pch=points[sp.targets$line],
	# labels=sp.targets$sample,main= "MDS all species 5h") 
	# plot this with other D. mel lines to check for dpecificity of C1
	# legend("topright", legend=levels(sp.targets$line),pch=points, col=colors, ncol=2)
	# This indicates that the differences between groups are larger than those within groups
	# Check outlier genes/samples: genes with highest mean and sd per sample
	# One gene is an outlier for mel, sim sec: FBgn0264695. FBgn0032350 and FBgn0000639 for yak
	# finding the 10 or so genes with the lowest dispersions and look at their counts to see if there is anything 		# amiss dge$tagwise.dispersions

	#Estimate dispersion

	#sp.dge.keep<-get(paste("sp.dge.keep",t,sep="."))
	sp.disp <- estimateDisp(sp.dge.keep, sp.design, robust=TRUE)
	#  plotBCV(sp.disp)
	sp.fit <- glmQLFit(sp.disp, sp.design, robust=TRUE)
	# plotQLDisp(sp.fit)
	sp.qlf <- glmQLFTest(sp.fit, contrast=sp.contrasts)		
	sp.FDR <- p.adjust(sp.qlf$table$PValue, method="BH")
	sp.qlf[which(sp.FDR<0.1),]

	#get the distribution of the dispersion: hist(sp.disp$tagwise.dispersion)
	#it shows that most are between 0-1 and few above 2
	#genes.filtered<-c("FBgn0027527","FBgn0000559","FBgn0264695")
	#genes.filtered<-names(sim.sd[sim.sd>10000])
	#sp.counts.filtered<-sp.counts[-which(rownames(sp.counts)%in%genes.filtered),]
	#sp.dge<-DGEList(counts=sp.counts.filtered, group=sp.group, genes=rownames(sp.counts.filtered))
	#sp.Dglm<-estimateGLMCommonDisp(sp.dge.keep,sp.design, verbose=T)
	#sp.Dglm.trended<-estimateGLMTrendedDisp(sp.Dglm, sp.design, verbose=T)
	#sp.Dglm.tagwise<-estimateGLMTagwiseDisp(sp.Dglm.trended, sp.design, prior.df=15) #,prior.df=20
	#sp.fit.trended <- glmFit(sp.Dglm.trended,sp.design)
	#sp.fit.tagwise <- glmFit(sp.Dglm.tagwise,sp.design)
	#sp.lrt <- glmLRT(sp.fit.tagwise,contrast=sp.contrasts)

	for (i in colnames(sp.contrasts)){
		sp.qlf <- glmQLFTest(sp.fit, contrast=sp.contrasts[,i])	
		sp.FDR <- p.adjust(sp.qlf$table$PValue, method="BH")
		#sp.lrt.tagwise<-glmLRT(sp.fit.tagwise, contrast=sp.contrasts[,i])
		summary(sp.de <- decideTestsDGE(sp.qlf,adjust.method="BH", p=0.1))
		sp.detags <- rownames(sp.dge.keep)[as.logical(sp.de)]
		sp.logFC<-sp.qlf$table$logFC[as.logical(sp.de)]
		m<-match(sp.detags,orthologs$V1)
		cg<-orthologs[m,]$V2
		pTable<-data.frame(FBgn_mel=sp.dge.keep$genes[as.logical(sp.de),],
			CG=cg,logFC=sp.logFC,FDR=sp.FDR[as.logical(sp.de)])
		pTable.i.t<-assign(paste("pTable.sp",i,t,sep="."),pTable)
	}
		assign(paste("sp.cpm",t,sep="."),cpm(sp.counts))
		assign(paste("sp.targets",t,sep="."),sp.targets)
}


#Get FBid and CG conversion
#m<-match(rownames(sp.counts.5),orthologs$V1)
#FBgn_CG<-data.frame(FBgn=orthologs[m,]$V1,CG=orthologs[m,]$V2)

#plot CPM 
# Edit plots to include both IDs
show.DEgenes<-function(FBid=NULL, tdata=NULL, cov=NULL, cg=NULL, q=NULL){#tdata=sp.dge.keep$counts,cov=sp.targets, q=NULL){ 
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
    title(main=paste(cg,"(", FBid, ")", sep=" "), sub=paste("FDR=", signif(q , digits=4)))
}

#Combine all contrasts in one table
pTable.spp.5<-pTable.sp.Par.5
pTable.spp.5<-rbind(pTable.spp.5,
		pTable.sp.ResSp.5[-which(pTable.sp.ResSp.5$FBgn_mel%in%pTable.spp.5$FBgn_mel),])
pTable.spp.5<-rbind(pTable.spp.5,
		pTable.sp.MelSimSecxPar.5[-which(pTable.sp.MelSimSecxPar.5$FBgn_mel%in%pTable.spp.5$FBgn_mel),])
pTable.spp.5<-rbind(pTable.spp.5,
		pTable.sp.MelSimxPar.5[-which(pTable.sp.MelSimxPar.5$FBgn_mel%in%pTable.spp.5$FBgn_mel),])



pdf("figures/CPMplots_sppPar_5h.pdf", paper="a4", height=10)
par(mfrow=c(4,3))
for (j in pTable.spp.5[order(pTable.spp.5$FDR),]$FBgn){
 show.DEgenes(j,cpm.5,targets.5,
 pTable.spp.5[which(pTable.spp.5$FBgn_mel==j),]$CG,
 pTable.spp.5[which(pTable.spp.5$FBgn_mel==j),]$FDR)
 }
dev.off()




#https://stat.ethz.ch/pipermail/bioconductor/2014-April/058822.html
#https://support.bioconductor.org/p/69374/

