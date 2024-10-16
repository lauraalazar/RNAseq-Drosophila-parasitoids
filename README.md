# RNAseq-Drosophila-parasitoids
RNAseq analysis of four closely related Drosophila species and two lines of D. melanogaster
experimentally selected for increased resistance after parasitization in two time points (5 and 50 hours). 
These samples are compared against non parasitized control samples. Each sample consisted of three biological replicates.
The samples were collected in the summer of 2011 and sequenced in 2012 in Illumina Hiseq2000, paired end 2x100.
The following steps were followed for each of 84 samples to obtain differentially expressed genes between control and parasitized. 
The results are presented in https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-3697-3 and the raw data can be downloaded from https://www.ebi.ac.uk/ena/browser/view/PRJEB15540

Quality Control
```
fastx_quality_stats -i dm_1_ACTTGA_L001_R1_001.fastq -Q33 -o dm1_R1.txt"
```

Alignment
```
gsnap --gunzip --nofails --novelsplicing=1 --format=sam -D {input.genome}  -d mel5.51 -t 24 dm_1_ACTTGA_L001_R1_001.fastq dm_1_ACTTGA_L001_R2_001.fastq > sample1.sam

samtools view -Sbh sample1.sam > sample.bam

samtools sort -n -@ 8 sample1.bam} sample1.sorted

samtools index sample1.sorted.bam
```
Remove duplicates

```
picard-tools MarkDuplicates INPUT=sample1.sorted.bam OUTPUT=sample1.sorted.dedp.bam METRICS_FILE=sample1.metrics \
		 MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=250  MAX_RECORDS_IN_RAM=50000  ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT \
		 REMOVE_DUPLICATES=True
```

Counts
```
samtools view sample1.sorted.dedup.bam | htseq-count - dmel5.51_converted.gff -s no -q > sample1.counts
```
The counts for all samples are in the folder "counts". These are used to find differential expression using edgeR and dexSeq in R. 

## Differential expression
Differential expression was calculated for different sets and normalization functions (LRT and QLfit) using the scripts in folder "scripts":
1) Species-specific (non-melanogaster sp): script_Dspp_specific_glmLRTfit.edgeR.R and script_Dspp_specific_glmQLfit.edgeR.R
2) All species together: script_Dspp_all_glmLRTfit.edgeR.R and script_Dspp_all_glmQLfit.edgeR.R
3) Within lines of D. melanogaster: script_Dmel_glmLRTfit.edgeR.R and script_Dmel_glmQLfit.edgeR.R, and script_Dmel.dexseq.R for differential exon usage

For figures see the script plots.R





