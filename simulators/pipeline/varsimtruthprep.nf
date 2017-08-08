#!/usr/bin/env nextflow

project ='version2.0'
outpath ="/data/backup/metacaller/simulatedgenomes/varsim"
projectpath="$outpath/$project/"

process simulate {

input:

output :

script:
"""
rm -f $outpath/$project/allVCFtruth/*
cp $outpath/$project/samples/VARSIM1/out/simu.truth.vcf $outpath/$project/allVCFtruth/truth.vcf
${NORMALISERTRUTH} $outpath/$project/allVCFtruth/truth.vcf ${REFPATH}
"""
}
