# RNAseq-Drosophila-parasitoids
RNAseq analysis of four closely related Drosophila species and two lines of D. melanogaster
experimentally selected for increased resistance after parasitization in two time points (5 and 50 hours). 
These samples are compared against non parasitized control samples. Each sample consisted of three biological replicates.
The samples were collected in the summer of 2011 and sequenced in 2012 in Illumina Hiseq2000, paired end 2x100.

RAW READS
...
QUALITY CONTROL
...
RE-ANNOTATING GENOMES
Following steps from:
http://bioinf.uni-greifswald.de/bioinf/wiki/pmwiki.php?n=IncorporatingRNAseq.GSNAP

Edit fastq files to include \1 or \2 to the headers:
cat dm_6_TAGCTT_L001_R2_001.fastq |  awk '{if (NR % 4 == 1) {split($1, arr, ":"); printf "%s_%s:%s:%s:%s:%s/%s\n", arr[1], arr[3], arr[4], arr[5], arr[6], arr[7],substr($2, 1, 1), $0} else if (NR % 4 == 3){print "+"} else {print $0} }' > dm_6_TAGCTT_L001_R2_001_ed.fastq

Build index:
gmap_build -d sec_r1.0 dsec-all-chromosome-r1.0.fasta -D sec_ref/

Align reads, convert sam to bam and sort first by 
gsnap --nofails --novelsplicing=1 --format=sam -D /gmap_db/sec_ref/ -d sec_r1.0 R1_ed.fastq R2_ed.fastq -t 32 | 
samtools view -buS -| samtools sort - | samtools sort -n dsec6.ss

filterBam --uniq --paired --in dsec6.ss.bam --out dsec6.f.bam
samtools sort dsec6.f.bam dsec6.fs
samtools view -H dsec6.fs.bam > header.txt

