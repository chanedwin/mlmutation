#!/bin/bash
###PRESTEP - bwaaligner
###POSTSTEP - vcftrawl/pyconcatenator 
## This programme takes in bam files, and produces VCF files using 8 variant callers##
### INITIALISATION ###

ORIGINPUTPATH=$1
INPUTPATH=$1
REFPATH=$3
SAMTOOLSPATH=/data/backup/metacaller/download/samtools-1.3/samtools
BCFPATH=/data/backup/metacaller/download/bcftools/bcftools-1.3/bcftools
VCFUTIL=/data/backup/metacaller/download/bcftools/bcftools-1.3/vcfutils.pl


#${SAMTOOLSPATH} mpileup -f ${REFPATH} $1  > $1.raw.vcf 
${SAMTOOLSPATH} mpileup -uf ${REFPATH} $1 | ${BCFPATH} call -c -v > $1.raw.bcf 
${BCFPATH} view $1.raw.bcf > $2
