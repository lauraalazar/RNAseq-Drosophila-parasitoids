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
	st.groups<-targets[which(targets$line==s & targets$time==t),]$treatment
	print(st.counts)
## collate all info and assign cpm before filtering low counts (for later steps and plotting)
	dge <-DGEList(counts=get(st.counts),group=st.groups)
	assign(paste(s,"cpm",t,sep="."),cpm(dge$counts))

#filter out very low counts
	keep<-rowSums(cpm(dge)>5) >= 3
	dge.keep<-dge[keep,]
	dge.keep$samples$lib.size <- colSums(dge.keep$counts) 

# normalize using TTM
	dge.keep<-calcNormFactors(dge.keep)
	st.dge<-assign(paste(s,t,"dge",sep="."),dge.keep)

# Estimate dispersion, create the design matrix (only compare ctlvspar) and calulate the glmfit 
	design <- model.matrix(~group,data=st.dge$samples)
	st.disp <- estimateDisp(st.dge, design)
#	dglm<-estimateGLMCommonDisp(st.dge, design, verbose=T)
#	dglm_tagwise<-estimateGLMTagwiseDisp(dglm, design, prior.df =30)
	fit <- glmFit(st.disp,design)#(dglm_tagwise,design)
	lrt <- glmLRT(fit)
	summary(de <- decideTestsDGE(lrt))

# extract FDR and fold changes and make table
	FDR<-p.adjust(lrt$table$PValue, method="BH")
	logFC<-lrt$table$logFC
#get the contrast matrix
	detags <- rownames(dge.keep)[as.logical(de)]
	if(any(grep("mel",s))==TRUE){
		m<-match(detags,orthologs$V1)
		cg<-orthologs[m,]$V2
		fb_mel = rownames(lrt$table[as.logical(de),])
		} else if(any(grep("mel",s))==FALSE){
		m<-match(detags,orthologs$V6)
		cg<-orthologs[m,]$V2
		fb_mel <- orthologs[m,]$V1
	}
	pTable<-data.frame(FB_sp=rownames(lrt$table[as.logical(de),]),FB_mel=fb_mel,CG=cg,logFC=logFC[as.logical(de)],FDR=FDR[as.logical(de)])
	sp.pTable<-assign(paste("LRT.pTable",s,t,sep="."),pTable)
#	assign(paste(s,"cpm",t,sep="."),cpm(dge.keep$counts))
	assign(paste(s,"targets",t,sep="."),targets[which(targets$line==s & targets$time==t),])
	}
}


#http://ww2.coastal.edu/kingw/statistics/R-tutorials/graphs.html
