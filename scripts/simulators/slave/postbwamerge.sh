#!/bin/bash

SAMTOOLSPATH=/data/backup/metacaller/download/samtools-1.3/samtools
OUTPUT="\$(ls -d \-\1 $1\)"
OUTPUT2=`find $1 -name '*.bam' |xargs`
echo "merged files were"
echo "${OUTPUT2}"

${SAMTOOLSPATH} merge $1/finalmerge.bam ${OUTPUT2}
echo "merging done"
mkdir -p $2
mv -f $1/finalmerge.bam $2/
echo "movement done"
