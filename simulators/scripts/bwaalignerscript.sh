#!/bin/bash
###PRESTEP - Varsim
###POSTSTEP - fullanalysis
# After using Varsim to create lanes of variants, this method aligns the variants using BWA mem algorithm (paired end manner). it also sorts the sam file by chromosome, converts it into a bam file and piles up all bamfiles in a region. It then indexes the piled up bam file
# MUST ENSURE THAT REFPATH IS ALREADY INDEXED WITH BWA

BWAPATH=/data/backup/metacaller/download/bwa-0.7.13/bwa
REFPATH=$3
SAMTOOLSPATH=/data/backup/metacaller/download/samtools-1.3/samtools

file=$1
SAMPLENAME="$(basename $file)" 
file2="$(sed s/read1/read2/g <<<$file)"
echo "aligning $file"
echo "aligning $file2"
echo "this is the sample name $SAMPLENAME"
${BWAPATH} mem -t 8 ${REFPATH} ${file} ${file2} > ${SAMPLENAME}.sam
${SAMTOOLSPATH} sort -o ${SAMPLENAME}.sorted.sam ${SAMPLENAME}.sam
${SAMTOOLSPATH} view -b -o ${SAMPLENAME}.sorted.bam ${SAMPLENAME}.sorted.sam
mv ${SAMPLENAME}.sorted.bam $2
