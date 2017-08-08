#!/bin/bash
REF=$1

SNP_RATE=0.00008
INDEL_RATE=0.00002
NUM_INSERTS=0.00005
NUM_DELETES=0.00005
NUM_MISMATCH=0.004
NUM_READS=100000000
DELETE_NAME=.temp.file.fa

#mason_variator --snp-rate ${SNP_RATE} --small-indel-rate ${INDEL_RATE} -ir $1 -ov $2 -of ${DELETE_NAME}


mason_simulator --illumina-prob-insert ${NUM_INSERTS} --illumina-prob-deletion ${NUM_DELETES} --illumina-prob-mismatch ${NUM_MISMATCH} -ir $1 -n ${NUM_READS} -iv $2 -o $3 -or $4 --num-threads 8
