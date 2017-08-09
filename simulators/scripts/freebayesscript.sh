#!/bin/bash
###PRESTEP - bwaaligner
###POSTSTEP - vcftrawl/pyconcatenator 
## This programme takes in bam files, and produces VCF files using 8 variant callers##
### INITIALISATION ###

INPUTPATH=$1
REFPATH=$3
FREEBAYESPATH=/data/backup/metacaller/download/freebayes/bin/freebayes

${FREEBAYESPATH} --fasta-reference ${REFPATH} $1 > $2
