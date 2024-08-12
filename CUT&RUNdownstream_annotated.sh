#!/bin/bash

# Arguments:
# 1 Name of the condition (assumes naming convention of Condition_replicate) for replicate 1 and 2 

# Requirements: 
# Required programs in path: MACS2, bedtools, HOMER
# Assumes 32 cores for threading, but adjust accordingly

# Call peaks with MACS2

macs2 callpeak -t $1_1.bowtie2.noBL.bam -c IgG_bowtie2.noBL.bam -f BAMPE --keep-dup all -n $1_1_p01 -p 1e-2

macs2 callpeak -t $1_2.bowtie2.noBL.bam -c IgG_bowtie2.noBL.bam -f BAMPE --keep-dup all -n $1_2_p01 -p 1e-2

# Overlap peaks to find consensus set

bedtools intersect -wa -a $1_1_p01_peaks.narrowPeak -b $1_2_p01_peaks.narrowPeak > $1overlap_Rep1width.bed

bedtools intersect -wa -a $1_2_p01_peaks.narrowPeak -b $1_1_p01_peaks.narrowPeak > $1overlap_Rep2width.bed

cat $1overlap_Rep1width.bed $1overlap_Rep2width.bed > $1cat.bed

bedtools sort -i $1cat.bed > $1sorted.bed

bedtools merge -i $1sorted.bed > $1_consensus.bed 

# Search for enriched motifs with HOMER, using the whole peak width 

findMotifsGenome.pl $1_consensus.bed hg38 $1consensus -p 32 -size given

# Merge replicates for visualizaiton and create normalized (signal per million reads) bedgraphs with MACS2

samtools merge M_$1.bam $1_1.bowtie2.noBL.bam $1_2.bowtie2.noBL.bam -@32 

samtools sort -@32 M_$1.bam > $1_merge.bam

macs2 callpeak -t $1_merge.bam -f BAMPE --keep-dup all --SPMR --bdg -n SPMR_$1 -p 1e-2

# Additional downstream analyses include signal intensity plots using ngsplots, annotating peaks to genes with GREAT (web server) and GO term enrichment with EnrichR. 






