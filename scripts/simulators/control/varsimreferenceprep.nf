#!/usr/bin/env nextflow

refpath ="/data/reference/human/GRCh37p5/fasta"
process prepare{

input:

output :

script:
"""
java -jar \${PICARDPATH} CreateSequenceDictionary R=$refpath/hs37d5.fa O=$refpath/hs37d5.dict
\${SAMTOOLSPATH} faidx $refpath/hs37d5.fa
"""

}
