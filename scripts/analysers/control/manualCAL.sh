#!/bin/bash

BASE='/data/backup/metacaller/stage'                                                                                                                              pr
project='mason_version1.0'
PROJECTOUT='mason_version1.0'
inpath="${BASE}/data/${project}"
inputpath="${inpath}"
OUTPATH="${BASE}/output/${PROJECTOUT}"
OUTPUTPATH="${OUTPATH}"
referencepath="${BASE}/reference/hs37d5.fa"

BASEPATH="/data/backup/metacaller/stage/deeplearnvcf"
BASEPATH="${BASEPATH}/version7.0-fulldataset-test-mason"
CALPATH="${BASEPATH}/control/CALwrapper.sh"

${CALPATH} ${OUTPUTPATH}/ANN/ANNgenerateresults.py ${OUTPUTPATH}/ANN/myXdata.txt.npy ${OUTPUTPATH}/ANN/myydata.txt.npy ${OUTPUTPATH}/ANN/ ${OUTPUTPATH}/ANN/samplelist.p ${OUTPUTPATH}/ANN/truthdict.p ${OUTPUTPATH}/ANN/callerlengths.txt.npy ${OUTPUTPATH}/ANN/vcf_dictionary.txt
