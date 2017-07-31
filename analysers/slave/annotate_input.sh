#!/bin/bash

module load bio/annovar
ANNOVAR_LOCATION=/data/reference/human/annovar/


table_annovar.pl $1 ${ANNOVAR_LOCATION} -buildver hg19 -outfile $2 -remove -protocol refGene,cytoBand,snp138,clinvar_20150629,ljb26_all,caddgt20 -operation g,r,f,f,f,f -nastring . --otherinfo -vcfinput

