An analysis of differential exon usage was done using the Bioconductor package DEXSeq 
(http://www.bioconductor.org/packages/release/bioc/vignettes/DEXSeq/inst/doc/DEXSeq.pdf).

1) Pre-processing:
-----------------
The analysis requires a pre-processing step of the annotation file using the script dexseq_prepare_annotation.py. 
I used the gtf file from ensemble Drosophila_melanogaster.BDGP5.74.gtf

$ python dexseq_prepare_annotation.py Drosophila_melanogaster.BDGP5.74.gtf dmel_flattened.gff

I needed to edit the output file to include "chr" in the first column:
$ awk '{$1="chr"$1}1' OFS='\t' dmel_flattened.gff > dmel_flattened_ed.gff

2) Counting:
-----------
Bam-files containing the mapped reads were used as input after removing duplicates and non-unique mapped reads

$ samtools view nsorted_sample_dedup.bam | python dexseq_count.py --paired=yes -s no dmel_flattened_ed.gff - sample_dexseq.counts

3) Dufferential Exon usage:
--------------------------
Steps are found in the script script_Dmel.dexseq.R
This scripts uses as input the count-files in the folder dexseq_counts and the file sample_info.txt


