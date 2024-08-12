#!/bin/bash

# Arguments:
# 1 Original file prefix: each sample contains 8 files (4 lanes read 1 + 4 lanes read 2). If different format, edit script accordingly. 
# 2 Desired output prefix: what the sample should be called downsteam (used to change file name) 

# Requirements: 
# Required programs in path: bbduk (bbmap), bowtie2, samtools, bedtools
# Edit script to add correct paths to bowtie2 genome assembly ("/pathtobowtie2index") and to CUT&RUN hg38 suspect list bed file("/pathtosuspectlist.bed")
# Assumes 32 cores for threading, but adjust accordingly. 

#Read trimming with bbduk
echo "trimming reads"

bash bbduk.sh in=$1_L001_R1_001.fastq.gz  in2=$1_L001_R2_001.fastq.gz out=$2.1_1.fastq out2=$2.1_2.fastq ref=adapters,artifacts literal=TATATATATATATATATATATATATATATATATATA,ATATATATATATATATATATATATATATATATATAT,GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG,CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

bash bbduk.sh in=$1_L002_R1_001.fastq.gz  in2=$1_L002_R2_001.fastq.gz out=$2.2_1.fastq out2=$2.2_2.fastq ref=adapters,artifacts literal=TATATATATATATATATATATATATATATATATATA,ATATATATATATATATATATATATATATATATATAT,GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG,CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

bash bbduk.sh in=$1_L003_R1_001.fastq.gz  in2=$1_L003_R2_001.fastq.gz out=$2.3_1.fastq out2=$2.3_2.fastq ref=adapters,artifacts literal=TATATATATATATATATATATATATATATATATATA,ATATATATATATATATATATATATATATATATATAT,GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG,CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

bash bbduk.sh in=$1_L004_R1_001.fastq.gz  in2=$1_L004_R2_001.fastq.gz out=$2.4_1.fastq out2=$2.4_2.fastq ref=adapters,artifacts literal=TATATATATATATATATATATATATATATATATATA,ATATATATATATATATATATATATATATATATATAT,GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG,CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

#Alignment to the hg38 human genome with bowtie2
echo "starting alignment to hg38"

bowtie2 "-p 32" "--local" "--very-sensitive-local" "--no-unal" "--no-mixed" "--no-discordant" "--phred33" "--dovetail" "-I 0" "-X 500" "-x" "/pathtobowtie2index" "-1" $2.1_1.fastq,$2.2_1.fastq,$2.3_1.fastq,$2.4_1.fastq "-2" $2.1_2.fastq,$2.2_2.fastq,$2.3_2.fastq,$2.4_2.fastq "-S" $2.bowtie2.sam

echo "creating bam file"

samtools view "-@ 32" "-b" $2.bowtie2.sam > $2.bowtie2.bam

#Samtools to fix improper mates, remove duplicates, and sort bam:

echo "fixing mates"

samtools fixmate "-@ 32" "-m" $2.bowtie2.bam $2.bowtie2.fixmate.bam

echo "sorting bam file"

samtools sort "-@ 32" $2.bowtie2.fixmate.bam "-o" $2.bowtie2.positionsort.bam

echo "marking up"

samtools markdup "-r" "-@ 32" $2.bowtie2.positionsort.bam $2.bowtie2.markdup.bam

echo "sorting marked up bam file"

samtools sort "-@ 32" $2.bowtie2.markdup.bam "-o" $2.bowtie2.final.bam

#Samtools to filter out regions in the CUT&RUN hg38 suspect list

echo "removing mitochondrial and blacklist reads"

bedtools intersect -abam $2.bowtie2.final.bam -b "/pathtosuspectlist.bed" -v | samtools sort -o $2.bowtie2.noBL.bam

echo "indexing bam file"

samtools index "-@ 32" $2.bowtie2.noBL.bam

echo "creating no blacklist bedgraph file"

bedtools genomecov "-bg" "-pc" "-ibam" $2.bowtie2.noBL.bam > $2.bowtie2.noBL.bedgraph

#Removing intermediate files

rm $2.bowtie2.sam
rm $2.bowtie2.bam
rm $2.bowtie2.fixmate.bam
rm $2.bowtie2.positionsort.bam
rm $2.bowtie2.markdup.bam
rm $2.bowtie2.final.bam

exit






