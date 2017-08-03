#!/bin/bash
### PRESTEP - Fullanalysis
### POSTSTEP - pyconcat
### this file processes the vcf files created by fullanalysis by compressing them and then indexing them. The mainstep comes after, which is to normalise using a variant
REFPATH=$2
BCFPATH=/data/backup/metacaller/download/bcftools/bcftools-1.3/bcftools
SAMTOOLSPATH=/data/backup/metacaller/download/samtools-1.3/samtools
PICARDPATH=/data/backup/metacaller/download/picard-tools-2.3.0/picard.jar
VTPATH=/data/backup/metacaller/download/vt/vt

bgzip -c $1 > $1.gz
tabix -f -p vcf $1.gz
${VTPATH} normalize $1.gz -r ${REFPATH} | ${VTPATH} uniq - -o $1.normalisedtrain.vcf
