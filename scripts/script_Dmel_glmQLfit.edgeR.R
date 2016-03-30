#rm(list = ls())
library(edgeR)
library(GO.db)
library(org.Dm.eg.db)


# set the settings for input and output
wd<-getwd()

DIR.IN<-"~/PhDpapers/RNAseq/counts/"
DIR.OUT<-"~/PhDpapers/RNAseq/analysis/mel/"
DIR.gff<-"~/PhDpapers/RNAseq/gff/"

# get information about samples
targets<-read.table("/home/laura/PhDpapers/RNAseq/targets.txt",header=T)
targets<-with(targets,targets[order(line,treatment),])

###################################
######    Get the Counts     ######
######  Create One matrix    ######
####  for all lines in each time ##
###################################
mel = c("melC1","melC2", "melS1", "melS2")
time<-c("5","50")

for(s in mel){
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
		assign(paste(s,t,sep="."),sp.counts)

		}
	}

#Put all counts together
mel.counts.5<-cbind(melC1.5,melC2.5,melS1.5,melS2.5)
mel.counts.50<-cbind(melC1.50,melC2.50,melS1.50,melS2.50)

mel.genes <-data.frame(FB=rownames(mel.counts.5))

# Read flybase file to convert FBgn to CG
ID_converter<-read.table(paste(DIR.gff,"gene_orthologs_fb_2015_02.tsv",sep=""))
m<-match(rownames(mel.counts.5),ID_converter$V1)
FBgn_CG<-data.frame(FBgn=ID_converter[m,]$V1,CG=ID_converter[m,]$V2)


for (t in time){
	#Get desgin matrix and contrasts
	#assign(paste("sp.Dglm",t,sep="."),
	mel.targets<-droplevels(targets[which(targets$line%in%mel & targets$time==t),])
	mel.targets$SC<-substr(mel.targets$line,4,4) #get the line (C or S)
	mel.targets$Batch<- factor(rep(1:3, dim(mel.targets)[1]/3)) 
	mel.group<-paste(mel.targets$SC, mel.targets$treatment, sep=".")
	mel.design<-model.matrix(~0+mel.group+Batch, data=mel.targets)
	colnames(mel.design)<-c("C.ctl", "C.par", "S.ctl", "S.par", "Batch2", "Batch3")
	mel.contrasts<-makeContrasts(CvsS.ctl=S.ctl-C.ctl, Par=(C.par+S.par)/2-(C.ctl+S.ctl)/2, 
			CSxTrt=(C.par-C.ctl)-(S.par-S.ctl), levels=mel.design)

	## collate all info
	mel.counts=get(paste("mel.counts",t,sep="."))
	mel.dge<-DGEList(counts=mel.counts, group=mel.group,genes=mel.genes)

	# filter out the tags that are not expressed at >5 CPM in at least 3 libraries
	keep<-rowSums(cpm(mel.dge)>5) >= 3
	mel.dge.keep<-mel.dge[keep,]
	mel.dge.keep$samples$lib.size <- colSums(mel.dge.keep$counts) 

	# normalize using TTM
	mel.dge.keep<-calcNormFactors(mel.dge.keep)
	assign(paste("mel.dge",t,sep="."),mel.dge.keep)
	#assign(paste("mel.dge.keep",t,sep="."),mel.dge.keep)

	#Estimate dispersion
	#mel.dge.keep<-get(paste("sp.dge.keep",t,sep="."))
	mel.disp <- estimateDisp(mel.dge.keep, mel.design, robust=TRUE)
	#  plotBCV(sp.disp)
	mel.fit <- glmQLFit(mel.disp, mel.design, robust=TRUE)
	# plotQLDisp(sp.fit)
	mel.qlf <- glmQLFTest(mel.fit, contrast=mel.contrasts)		
	mel.FDR <- p.adjust(mel.qlf$table$PValue, method="BH")
	mel.qlf[which(mel.FDR<0.05),]



	#mel.Dglm.trended <- estimateGLMTrendedDisp(mel.Dglm,mel.design)
	#mel.Dglm.tagwise <- estimateGLMTagwiseDisp(mel.Dglm.trended,mel.design)
	#mel.fit.trended <- glmFit(mel.Dglm.trended,mel.design)
	#mel.fit.tagwise <- glmFit(mel.Dglm.tagwise,mel.design)
	#mel.lrt.tagwise <- glmLRT(mel.fit.tagwise, contrast=mel.contrasts)	
	for (i in colnames(mel.contrasts)){
		mel.qlf <- glmQLFTest(mel.fit, contrast=mel.contrasts[,i])	
		mel.FDR <- p.adjust(mel.qlf$table$PValue, method="BH")
		#sp.lrt.tagwise<-glmLRT(sp.fit.tagwise, contrast=sp.contrasts[,i])
		summary(mel.de <- decideTestsDGE(mel.qlf,adjust.method="BH", p=0.05))
		mel.detags <- rownames(mel.dge.keep)[as.logical(mel.de)]
		mel.logFC<-mel.qlf$table$logFC[as.logical(mel.de)]
		cg=FBgn_CG[mel.dge.keep$genes[as.logical(mel.de),],]$CG
		pTable<-data.frame(FBgn_mel=mel.dge.keep$genes[as.logical(mel.de),],
			CG=cg,logFC=mel.logFC,FDR=mel.FDR[as.logical(mel.de)])
		pTable.i.t<-assign(paste("QL.pTable.mel",i,t,sep="."),pTable)
	}
		assign(paste("mel.cpm",t,sep="."),cpm(mel.counts))
		assign(paste("mel.targets",t,sep="."),mel.targets)
}

###################
### GO analysis ###
###################
#FBids<-as.character(pTable.mel.Par.5$FBgn_mel)
#eids<-select(org.Dm.eg.db, FBids, "ENTREZID", "FLYBASE")[2]
#go<-goana(eids,species="Dm")
#topGO(go,n=50)
#keg<-kegga(eids,species="Dm")



		#mel.lrt.tagwise<-glmLRT(mel.fit.tagwise, contrast=mel.contrasts[,i])
		#summary(mel.de.tagwise <- decideTestsDGE(mel.lrt.tagwise, p=0.05))
		#mel.detags.tagwise <- rownames(mel.dge.keep)[as.logical(mel.de.tagwise)]
		#mel.FDR.tagwise<-p.adjust(mel.lrt.tagwise$table$PValue, method="BH")
		#mel.logFC.tagwise<-mel.lrt.tagwise$table$logFC
		#pTable<-data.frame(FBgn_mel=mel.dge.keep$genes[as.logical(mel.de.tagwise),],
		#CG=FBgn_CG[mel.dge.keep$genes[as.logical(mel.de.tagwise),],]$CG,
		#logFC=mel.logFC.tagwise[as.logical(mel.de.tagwise)],
		#FDR=mel.FDR.tagwise[as.logical(mel.de.tagwise)])
		#assign(paste("pTable.mel",i,t,sep="."),pTable)

### Tables

#write.table(pTable.tagwise.5.Par,"mel_5_par.txt")
#write.table(pTable.tagwise.50.Par,"mel_50_par.txt")
#write.table(pTable.tagwise.50.CvsS_ctl,"mel_50_csctl.txt")
#write.table(pTable.tagwise.50.CSxTrt,"mel_50_cspar.txt")

#print(xtable(pTable.tagwise3.Par[with(pTable.tagwise3.Par,order(pTable.tagwise3.Par$logFC,decreasing=TRUE)),-2]),include.rownames=FALSE)

#https://pythonhosted.org/bioservices/kegg_tutorial.html#look-for-pathway-by-genes-i-e-ids-or-usual-name
