rm(list = ls())

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

############################################################################
# start statistical analysis of differential expression
#using contrasts to compare differences in expression for all species
# 1) Par = Par vs Ctl for all species
# 2) SpxTrt = the interaction term: Differences in response to Par between spp. I made each possible contrast manually to explore the data.
############################################################################
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
		sp.contrasts<-makeContrasts(Par=(mel.par+sec.par+sim.par+yak.par)/4-(mel.ctl+sec.ctl+sim.ctl+yak.ctl)/4, #all species
			#Clade specific (all possible combinations)
			MelxPar=(mel.par-mel.ctl), SecxPar=(sec.par-sec.ctl), SimxPar=(sim.par-sim.ctl), YakxPar=(yak.par-yak.ctl), 
			MelSimxPar=(mel.par+sim.par)/2-(mel.ctl+sim.ctl)/2,
			MelSecxPar=(mel.par+sec.par)/2-(mel.ctl+sec.ctl)/2,
			MelYakxPar=(mel.par+yak.par)/2-(mel.ctl+yak.ctl)/2,
			SecSimxPar=(sec.par+sim.par)/2-(sec.ctl+sim.ctl)/2,	
			SecYacxPar=(sec.par+yak.par)/2-(sec.ctl+yak.ctl)/2,
			SimYakxPar=(yak.par+sim.par)/2-(yak.ctl+sim.ctl)/2,
			MelSimSecxPar=(mel.par+sec.par+sim.par)/3-(mel.ctl+sec.ctl+sim.ctl)/3,
			SimSecYakxPar=(sim.par+sec.par+yak.par)/3-(sim.ctl+sec.ctl+yak.ctl)/3,
			MelSecYakxPar=(mel.par+sec.par+yak.par)/3-(mel.ctl+sec.ctl+yak.ctl)/3,
			#Resistant species
			MelSimYakxPar=(mel.par+sim.par+yak.par)/3-(mel.ctl+sim.ctl+yak.ctl)/3,
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
			#Clade specific (all possible combinations)
			MelxPar=(mel.par-mel.ctl), SecxPar=(sec.par-sec.ctl), SimxPar=(sim.par-sim.ctl), YakxPar=(yak.par-yak.ctl), 
			MelSimxPar=(mel.par+sim.par)/2-(mel.ctl+sim.ctl)/2,
			MelSecxPar=(mel.par+sec.par)/2-(mel.ctl+sec.ctl)/2,
			MelYakxPar=(mel.par+yak.par)/2-(mel.ctl+yak.ctl)/2,
			SecSimxPar=(sec.par+sim.par)/2-(sec.ctl+sim.ctl)/2,	
			SecYacxPar=(sec.par+yak.par)/2-(sec.ctl+yak.ctl)/2,
			SimYakxPar=(yak.par+sim.par)/2-(yak.ctl+sim.ctl)/2,
			MelSimSecxPar=(mel.par+sec.par+sim.par)/3-(mel.ctl+sec.ctl+sim.ctl)/3,
			SimSecYakxPar=(sim.par+sec.par+yak.par)/3-(sim.ctl+sec.ctl+yak.ctl)/3,
			MelSecYakxPar=(mel.par+sec.par+yak.par)/3-(mel.ctl+sec.ctl+yak.ctl)/3,
			#Resistant species
			MelSimYakxPar=(mel.par+sim.par+yak.par)/3-(mel.ctl+sim.ctl+yak.ctl)/3,
			levels=sp.design)
	}

	## collate all info
	sp.dge<-DGEList(counts=sp.counts, group=sp.group, genes=rownames(sp.counts))

	# filter out the tags that are not expressed at >5 CPM in at least 3 libraries
	keep<-rowSums(cpm(sp.dge)>5) >= 3
	sp.dge.keep<-sp.dge[keep,]
	sp.dge.keep$samples$lib.size <- colSums(sp.dge.keep$counts)

	# normalize using TTM
	sp.dge.keep<-calcNormFactors(sp.dge.keep)
	assign(paste("sp.dge",t,sep="."),sp.dge.keep)


	#Estimate dispersion using glmLRT
	sp.disp <- estimateDisp(sp.dge.keep, sp.design)
	#  plotBCV(sp.disp)
	sp.fit <- glmFit(sp.disp, sp.design, robust=TRUE)
	sp.lrt <- glmLRT(sp.fit, contrast=sp.contrasts)		
	sp.FDR <- p.adjust(sp.lrt$table$PValue, method="BH")

	for (i in colnames(sp.contrasts)){
		sp.lrt <- glmLRT(sp.fit, contrast=sp.contrasts[,i])	
		sp.FDR <- p.adjust(sp.lrt$table$PValue, method="BH")
		summary(sp.de <- decideTestsDGE(sp.lrt,adjust.method="BH", p=0.1))
		sp.detags <- rownames(sp.dge.keep)[as.logical(sp.de)]
		sp.logFC<-sp.lrt$table$logFC[as.logical(sp.de)]
		m<-match(sp.detags,orthologs$V1)
		cg<-orthologs[m,]$V2
		pTable<-data.frame(FBgn_mel=sp.dge.keep$genes[as.logical(sp.de),],
			CG=cg,logFC=sp.logFC,FDR=sp.FDR[as.logical(sp.de)])
		# Keep only genes that have a minimum logFC of 2
		pTable<-pTable[which(abs(pTable$logFC)>1),]
		pTable.i.t<-assign(paste("LRT.pTable.sp",i,t,sep="."),pTable)
		
	}
		assign(paste("sp.cpm",t,sep="."),cpm(sp.counts))
		assign(paste("sp.targets",t,sep="."),sp.targets)
}


#Combine all clade specific contrasts in one table with unique ID and FC and FDR from the 1st contrast
LRT.pTable.sp.Clade.5 <- as.data.frame(setNames(replicate(4,numeric(0), simplify = F), 
			c("FBgn_mel","CG","logFC","FDR")))
LRT.pTable.sp.Clade.50 <- as.data.frame(setNames(replicate(4,numeric(0), simplify = F), 
			c("FBgn_mel","CG","logFC","FDR")))

for (i in colnames(sp.contrasts)){
	ptable.5<-(paste("LRT.pTable.sp",i,5,sep="."))
	ptable.50<-(paste("LRT.pTable.sp",i,50,sep="."))
	if(i!="Par"){
		LRT.pTable.sp.Clade.5<-rbind(LRT.pTable.sp.Clade.5,
					get(ptable.5)[which(!(get(ptable.5)$FBgn_mel%in%LRT.pTable.sp.Clade.5$FBgn_mel)),])
		LRT.pTable.sp.Clade.50<-rbind(LRT.pTable.sp.Clade.50,
					get(ptable.50)[which(!(get(ptable.50)$FBgn_mel%in%LRT.pTable.sp.Clade.50$FBgn_mel)),])
	}
}



#########################################
###### Statistical analysys #############
###### species one-by-one   #############
#########################################
# Here I start with D. melanogaster and
# add species one by one

for (t in time){
	sp.genes = NULL
	subgroup.list = NULL
	gene_list=NULL
	subgroup.counts = NULL
	subgroup.targets = NULL
	subgroup.intersect = NULL
	combination = NULL
	for(s in species){
		print(paste(s,t,sep="."))
		sp.targets<-droplevels(targets[which(targets$line%in%s & targets$time==t),])
		replicates<-c(1:length(which(sp.targets$treatment=="ctl")),1:length(which(sp.targets$treatment=="par")))
		sp.targets$Batch = replicates
		sp.counts<-get(paste(s,t,sep="."))
		s.t.genes<-paste(s,t,'genes',sep='.')
		assign(s.t.genes,rownames(sp.counts))
		gene_list<-append(gene_list,s.t.genes)
		subgroup.intersect<-Reduce(intersect,mget(gene_list))
		count.intersect<-sp.counts[subgroup.intersect,]	
		subgroup.counts<-cbind(subgroup.counts[subgroup.intersect,],count.intersect)
		subgroup.list<-append(subgroup.list,s)		
		combination<-paste(combination,s,sep="")
		subgroup.targets<-rbind(subgroup.targets,sp.targets)
		assign(combination,subgroup.counts)

		subgroup.group<-paste(subgroup.targets$line, subgroup.targets$treatment, sep=".")
		subgroup.design<-model.matrix(~0+subgroup.group+Batch,data=subgroup.targets)
		colnames(subgroup.design)<-gsub("subgroup.group","",colnames(subgroup.design))
		if(length(subgroup.list)==1){ 
			subgroup.contrast<-makeContrasts(Par=melC1.par-melC1.ctl,levels=subgroup.design)
		}else if(length(subgroup.list)==2){
			subgroup.contrast<-makeContrasts(Par=(melC1.par+sim.par)-(melC1.ctl+sim.ctl),levels=subgroup.design)
		}else if(length(subgroup.list)==3){
			subgroup.contrast<-makeContrasts(Par=(melC1.par+sim.par+sec.par)-(melC1.ctl+sim.ctl+sec.ctl),levels=subgroup.design)
		}else if(length(subgroup.list)==4){
			subgroup.contrast<-makeContrasts(Par=(melC1.par+sim.par+sec.par+yak.par)-(melC1.ctl+sim.ctl+sec.ctl+yak.ctl),levels=subgroup.design)
		}

		subgroup.dge <-DGEList(counts=subgroup.counts, group=subgroup.group, genes=rownames(subgroup.counts)) 
		keep<-rowSums(cpm(subgroup.dge)>5) >= 3
		subgroup.dge.keep<-subgroup.dge[keep,]
		subgroup.dge.keep<-calcNormFactors(subgroup.dge.keep)
		subgroup.dge.keep$samples$lib.size <- colSums(subgroup.dge.keep$counts) 
		subgroup.disp <- estimateDisp(subgroup.dge.keep, subgroup.design, robust=TRUE)
		subgroup.fit <- glmQLFit(subgroup.disp, subgroup.design, robust=TRUE)
		subgroup.qlf <- glmQLFTest(subgroup.fit, contrast=subgroup.contrast)	
		subgroup.FDR <- p.adjust(subgroup.qlf$table$PValue, method="BH")
		summary(subgroup.de <- decideTestsDGE(subgroup.qlf,adjust.method="BH", p=0.1))
		subgroup.detags <- rownames(subgroup.dge.keep)[as.logical(subgroup.de)]
		subgroup.logFC<-subgroup.qlf$table$logFC[as.logical(subgroup.de)]
		m<-match(subgroup.detags,orthologs$V1)
		cg<-orthologs[m,]$V2
		pTable<-data.frame(FBgn_mel=subgroup.dge.keep$genes[as.logical(subgroup.de),],
			CG=cg,logFC=subgroup.logFC,FDR=subgroup.FDR[as.logical(subgroup.de)])
		pTable.i.t<-assign(paste("pTable.subgroup",combination,t,sep="."),pTable)
	}
}
		

