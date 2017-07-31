#!/bin/bash
###PRESTEP -VARSIM
###POSTSTEP - testonlyTruthTable.sh
REFPATH=/data/reference/human/fasta/ucsc.hg19.fa
BCFPATH=/data/backup/metacaller/download/bcftools/bcftools-1.3/bcftools
SAMTOOLSPATH=/data/backup/metacaller/download/samtools-1.3/samtools
PICARDPATH=/data/backup/metacaller/download/picard-tools-2.3.0/picard.jar
VTPATH=/data/backup/metacaller/download/vt/vt
OUTPUTPATH=work/normtruth
GATK=/data/backup/metacaller/download/GenomeAnalysisTK.jar


sed -i '3i##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of Structural Variant">' $1
sed -i '4i##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of Structural Variant">' $1
cat $1 | grep \# > $1.headers
cat $1 | awk -F '[\t,]' '$4 != $5'| awk -F '[\t,]' '$4 != $6'| awk -F '[\t,]' '$4 != $7' | grep -v \# > $1.remodel
cat "$1.remodel" >> "$1.headers"
java -jar ${PICARDPATH} UpdateVcfSequenceDictionary I=$1.headers O=$1.temp SEQUENCE_DICTIONARY=${REFPATH}
awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' $1.temp > $1.vcf
java -jar ${PICARDPATH} SortVcf I=$1.vcf O=$1.sorted.vcf
java -jar ${GATK} \
    -T ValidateVariants \
    -R ${REFPATH} \
    -V $1.sorted.vcf
bgzip -c $1.sorted.vcf > $1.gz
tabix -f -p vcf $1.gz
${VTPATH} normalize $1.gz -n -r ${REFPATH} -o $1.normalisedtruth.vcf
cat $1.normalisedtruth.vcf | grep -v .\|.\: > $1.final.vcf
