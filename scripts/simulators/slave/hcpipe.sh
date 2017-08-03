#!/bin/bash
###PRESTEP - bwaaligner
###POSTSTEP - vcftrawl/pyconcatenator 
## This programme takes in bam files, and produces VCF files using 8 variant callers##
### INITIALISATION ###

ORIGINPUTPATH=$1
INPUTPATH=$1
REFPATH=$3
GATK=/data/backup/metacaller/download/GenomeAnalysisTK.jar

java -jar ${GATK} -R ${REFPATH} -nct 5 -T HaplotypeCaller -I $1 -o $2 -U ALLOW_SEQ_DICT_INCOMPATIBILITY -mbq 1 -stand_call_conf 0 -minReadsPerAlignStart 1 --minPruning 2

