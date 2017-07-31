#!/bin/bash
###PRESTEP - bwaaligner
###POSTSTEP - vcftrawl/pyconcatenator 
## This programme takes in bam files, and produces VCF files using 8 variant callers##
### INITIALISATION ###

ORIGINPUTPATH=$1
INPUTPATH=$1
REFPATH=$6
REFCHROMPATH=/data/reference/human/fasta/chrom/
OUTPUTFILE=/data/project/pdx/VMCoutput8.6.16/
BAMPATH=/data/project/pdx/VMCoutput/Reorderedbam.bam
PICARDPATH=/data/backup/metacaller/download/picard-tools-2.3.0/picard.jar

#VMC DEPENDENCIES
GATK=/data/backup/metacaller/download/GenomeAnalysisTK.jar
SAMTOOLSPATH=/data/backup/metacaller/download/samtools-1.3/samtools
FREEBAYESPATH=/data/backup/metacaller/download/freebayes/bin/freebayes
BCFPATH=/data/backup/metacaller/download/bcftools/bcftools-1.3/bcftools
VCFUTIL=/data/backup/metacaller/download/bcftools/bcftools-1.3/vcfutils.pl

BWAPATH=/data/backup/metacaller/download/bwa-0.7.13/bwa
SAMTOOLSPATH=/data/backup/metacaller/download/samtools-1.3/samtools
BREAKPOINTLIB=/data/backup/metacaller/download/breakpointlib/breakseq2_bplib_20150129.fna
THREADS=16
FILE_NAME=${INPUTPATH}
SAMPLE_NAME=$(basename $FILE_NAME)


##BREAKSEQ 2##

#check indexing of breakpoint lib by creating index and dictionary
if [ $1 = "ref" ]; then  
docker exec ubuntu-docker ${BWAPATH} index ${BREAKPOINTLIB}

docker exec ubuntu-docker ${SAMTOOLSPATH} faidx ${BREAKPOINTLIB}
 
docker exec ubuntu-docker java -jar ${PICARDPATH} CreateSequenceDictionary REFERENCE=/data/backup/metacaller/download/breakpointlib/breakseq2_bplib_20150129.fna OUTPUT=/data/backup/metacaller/download/breakpointlib/breakseq2_bplib_20150129.dict
fi 
#execute breakseq2 program

docker exec ubuntu-docker mkdir -p /data/backup/metacaller/temp/$2/work

FILE=$(cat $3)

docker cp $FILE ubuntu-docker:/data/backup/metacaller/temp/$2/

docker exec ubuntu-docker chmod +x /data/backup/metacaller/temp/$2/$4

docker exec ubuntu-docker ${SAMTOOLSPATH} index /data/backup/metacaller/temp/$2/$4

docker exec ubuntu-docker run_breakseq2.py --reference ${REFPATH} --sample ARTIFICIAL --bams /data/backup/metacaller/temp/$2/$4 --bwa ${BWAPATH} --samtools ${SAMTOOLSPATH} \
          --bplib ${BREAKPOINTLIB} --nthreads ${THREADS} --work /data/backup/metacaller/temp/$2/work

docker exec ubuntu-docker cp /data/backup/metacaller/temp/$2/work/breakseq.vcf.gz $5

#docker exec ubuntu-docker rm -r /data/backup/metacaller/temp/$2
