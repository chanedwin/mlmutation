#!/bin/bash
###PRESTEP - bwaaligner
###POSTSTEP - vcftrawl/pyconcatenator 
## This programme takes in bam files, and produces VCF files using 8 variant callers##
### INITIALISATION ###

ORIGINPUTPATH=$1
INPUTPATH=$1
REFPATH=$4

PINDPATH=/data/backup/metacaller/download/pindel/pindel
P2VCFPATH=/data/backup/metacaller/download/pindel/pindel2vcf

#generate pindel config file
echo "$1 250 $2" > $2.pindel

#run pindel file
${PINDPATH} -T 12 -f ${REFPATH} -i $2.pindel -o pindel.$2 -d 2 -m 1

#convert to vcf file
${P2VCFPATH} -P pindel.$2 -r ${REFPATH} -R hg19 -d 2009 -v $3
