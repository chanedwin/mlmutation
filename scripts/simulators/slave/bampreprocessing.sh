#/bin/bash
###PRESTEP - bwaaligner
###POSTSTEP - vcftrawl/pyconcatenator 
## This programme takes in bam files, and produces VCF files using 8 variant callers##
### INITIALISATION ###

INPUTPATH=$1
REFPATH=$2
PICARDPATH=/data/backup/metacaller/download/picard-tools-2.3.0/picard.jar

#VMC DEPENDENCIES
SAMTOOLSPATH=/data/backup/metacaller/download/samtools-1.3/samtools
GATK=/data/backup/metacaller/download/GenomeAnalysisTK.jar

java -jar ${PICARDPATH} AddOrReplaceReadGroups \
    I=$1\
    O=$1.addedgroups1.bam \
    RGID=4 \
    RGLB=lib1 \
    RGPL=illumina \
    RGPU=unit1 \
    RGSM=20
echo "done adding read groups"

java -Djava.io.tmpdir=/data/backup/metacaller/temp/java -jar ${PICARDPATH} ReorderSam \
    I= $1.addedgroups1.bam \
    O= $1.reordered2.bam \
    REFERENCE= ${REFPATH}\
    TMP_DIR=/data/backup/metacaller/temp/java

echo "reordering done"

${SAMTOOLSPATH} sort -o $1.sorted3.bam $1.reordered2.bam
echo "sorting done"

${SAMTOOLSPATH} index $1.sorted3.bam
echo "indexing done"

java -jar ${GATK} \
    -T RealignerTargetCreator \
    -R ${REFPATH} \
    -I $1.sorted3.bam \
    -o $1.target.intervals
echo "realignment first part done"

java -jar ${GATK} \
    -T IndelRealigner \
    -targetIntervals $1.target.intervals\
    -R ${REFPATH} \
    -I $1.sorted3.bam \
    -o $1.bam 
echo "realignment second part done"

${SAMTOOLSPATH} index $1.bam

echo "indexing done"
