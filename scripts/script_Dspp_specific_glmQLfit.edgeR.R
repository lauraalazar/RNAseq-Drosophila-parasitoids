#rm(list = ls())
#load required packages
library(edgeR)

# set the settings for input and output
wd<-getwd()

DIR.IN<-"~/PhDpapers/RNAseq/counts/"
DIR.OUT<-"~/storage/PhD/RNAseq/analysis/spp_specific/"
DIR.gff<-"~/PhDpapers/RNAseq/gff/"

# get information about samples
targets<-read.table("../targets.txt",header=T)
targets<-targets[-c(48),] #yak48


###################################
######    Get the Counts     ######
######  Create One matrix    ######
######  for each sp and time ######
###################################

time<-c("5","50")
species = c("sim", "sec", "melC1", "yak")


#For each species, open its directory, get the counts from the file and 
#bind it to the matrix
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
		assign(paste(s,t,sep="."),sp.counts)
	}

}


#########################################
###### Statistical analysys #############
###### species-specific     #############
#########################################
orthologs<-read.table(paste(DIR.gff,"gene_orthologs_fb_2015_02.tsv",sep=""))

for(s in species){
	for(t in time){
	st.counts<-paste(s,t,sep=".")
	st.targets<-targets[which(targets$line==s & targets$time==t),]
	st.group<-st.targets$treatment
	print(st.counts)
## collate all info and assign cpm before filtering low counts (for later steps and plotting)
	st.dge <-DGEList(counts=get(st.counts),group=st.group,genes=rownames(st.counts))
	assign(paste(s,"cpm",t,sep="."),cpm(st.dge$counts))

#filter out very low counts
	keep<-rowSums(cpm(st.dge)>5) >= 3
	st.dge.keep<-st.dge[keep,]
#	dge.keep$samples$lib.size <- colSums(dge.keep$counts) 

# normalize using TTM
	st.dge.keep<-calcNormFactors(st.dge.keep)
	plotMD(cpm(st.dge.keep,log=TRUE),column=1)
	abline(h=0, col="red", lty=2, lwd=2)
	colors=c("blue","green")
	plotMDS(st.dge.keep, col=colors[st.group])
# Estimate dispersion, create the design matrix (only compare ctlvspar) and calulate the glmfit 
	#design <- model.matrix(~group,data=st.dge$samples)
	st.design<-model.matrix(~0+st.group)
	colnames(st.design)<-levels(st.group)
	st.con <- makeContrasts(par - ctl, levels=st.design)
	st.dge.keep <- estimateDisp(st.dge.keep, st.design, robust=TRUE)
	plotBCV(st.dge.keep)
	st.fit <- glmQLFit(st.dge.keep, st.design, robust=TRUE)
	st.qlf <- glmQLFTest(st.fit, contrast=st.con)		
	st.FDR <- p.adjust(st.qlf$table$PValue, method="BH")
	summary(de <- decideTestsDGE(st.qlf,adjust.method="BH", p=0.1))
	st.logFC<-st.qlf$table$logFC[as.logical(de)]
#	dglm<-estimateGLMCommonDisp(st.dge, design, verbose=T)
#	dglm_tagwise<-estimateGLMTagwiseDisp(dglm, design, prior.df =30)
#	fit <- glmFit(dglm_tagwise,design)
#	lrt <- glmLRT(fit)
#	summary(de <- decideTestsDGE(lrt))

# extract FDR and fold changes and make table
#	FDR<-p.adjust(lrt$table$PValue, method="BH")
#	logFC<-lrt$table$logFC
#get the contrast matrix
	st.detags <- rownames(st.dge.keep)[as.logical(de)]
	if(any(grep("mel",s))==TRUE){
		m<-match(st.detags,orthologs$V1)
		cg<-orthologs[m,]$V2
		fb_mel = rownames(st.qlf$table[as.logical(de),])
		} else if(any(grep("mel",s))==FALSE){
		m<-match(st.detags,orthologs$V6)
		cg<-orthologs[m,]$V2
		fb_mel <- orthologs[m,]$V1
	}
	pTable<-data.frame(FBgn_sp=st.dge.keep$genes[as.logical(de),],
			CG=cg,logFC=st.logFC,FDR=st.FDR[as.logical(de)])
#	pTable<-data.frame(FB_sp=rownames(lrt$table[as.logical(de),]),FB_mel=fb_mel,CG=cg,logFC=logFC[as.logical(de)],FDR=FDR[as.logical(de)])
	sp.pTable<-assign(paste("QL.pTable",s,t,sep="."),pTable)
#	assign(paste(s,"cpm",t,sep="."),cpm(dge.keep$counts))
	assign(paste(s,"targets",t,sep="."),st.targets[which(targets$line==s & targets$time==t),])
	}
}



