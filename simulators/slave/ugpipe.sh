#!/bin/bash
###PRESTEP - bwaaligner
###POSTSTEP - vcftrawl/pyconcatenator 
## This programme takes in bam files, and produces VCF files using 8 variant callers##
### INITIALISATION ###

ORIGINPUTPATH=$1
INPUTPATH=$1
REFPATH=$3
GATK=/data/backup/metacaller/download/GenomeAnalysisTK.jar

java -jar ${GATK} -glm both -R ${REFPATH} -T UnifiedGenotyper -I $1 -o $2 -U ALLOW_SEQ_DICT_INCOMPATIBILITY -stand_call_conf 1 -mbq 1  -minIndelCnt 1  -minIndelFrac 1

