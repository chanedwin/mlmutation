#!/usr/bin/env nextflow


project ='version2.0'
outpath ="/data/backup/metacaller/simulatedgenomes/mason"
projectpath="$outpath/$project/"
left_alignment="$outpath/$project/samples/simulatedfqleft.fq"
right_alignment="$outpath/$project/samples/simulatedfqright.fq"
truth_path="$outpath/$project/samples/truth.vcf"
bam_output_file="$outpath/$project/output/aligned.bam"
sample_sam="$outpath/$project/output/sample.sam"
sortedsample_sam="$outpath/$project/output/sortedsample.sam"
bam_output_file_processed="$outpath/$project/output/aligned.bam.bam"
num =Channel.from(1)
num_1 =Channel.from(1)
num_2 =Channel.from(1)
num_3 =Channel.from(1)
num_4 =Channel.from(1)
num_5 =Channel.from(1)
num_6 =Channel.from(1)

process simulate {

input:
val x from num 

output :
file fq1 into fq11

script:
"""
mkdir -p $outpath/$project/samples/
mkdir -p $outpath/$project/output/
#${MASONPATH} \${REFPATH} $outpath/$project/samples/truth.vcf $left_alignment $right_alignment 
echo "message done" > fq1
"""
}

process alignment1{

input:
file fastq from fq11

output:
file signal into nextstage1

script:
"""
\${BWAPATH} mem -t 8 "\${REFPATH}" $left_alignment $right_alignment > $sample_sam;
\${SAMTOOLSPATH} sort -o $sortedsample_sam $sample_sam;
\${SAMTOOLSPATH} view -b -o $bam_output_file $sortedsample_sam;
echo "1" > signal
"""
}

process bringtogether{

input:
file signal from nextstage1

output:
file bam into bs2bam, pindelbam, fbbam,ugbam,hcbam,stbam

script:
"""
$POSTBWAPROCESS $bam_output_file \${REFPATH}
echo $bam_output_file_processed > bam
#${NORMALISERTRUTH} $truth_path \${REFPATH}
"""

}


process bs2train {

input:
file bam from bs2bam

script:
"""
"""
}

process pindeltrain {

input:
file bam from pindelbam


script:
"""
mkdir -p $outpath/$project/allVCFtrain/
temp=\$(cat $bam)
${PINDELPIPE} \$temp $project $outpath/$project/allVCFtrain/pindel.vcf \${REFPATH}
${NORMALISER} $outpath/$project/allVCFtrain/pindel.vcf \${REFPATH}
"""
}

process ugtrain {

input:
file bam from ugbam



script:
"""
mkdir -p $outpath/$project/allVCFtrain/
temp=\$(cat $bam)
${UGPIPE} \$temp $outpath/$project/allVCFtrain/ug.vcf \${REFPATH}
${NORMALISER} $outpath/$project/allVCFtrain/ug.vcf \${REFPATH}
"""
}

process hctrain {

input:
file bam from hcbam


script:
"""
mkdir -p $outpath/$project/allVCFtrain/
temp=\$(cat $bam)
${HCPIPE} \$temp $outpath/$project/allVCFtrain/hc.vcf \${REFPATH}
${NORMALISER} $outpath/$project/allVCFtrain/hc.vcf \${REFPATH}
"""
}

process fbtrain {

input:
file bam from fbbam

script:
"""
mkdir -p $outpath/$project/allVCFtrain/
temp=\$(cat $bam)
${FBPIPE} \$temp $outpath/$project/allVCFtrain/fb.vcf \${REFPATH}
${NORMALISER} $outpath/$project/allVCFtrain/fb.vcf \${REFPATH}
"""
}

process sttrain {

input:
file bam from stbam

script:
"""
mkdir -p $outpath/$project/allVCFtrain/
temp=\$(cat $bam)
${STPIPE} \$temp $outpath/$project/allVCFtrain/st.vcf \${REFPATH}
${NORMALISER} $outpath/$project/allVCFtrain/st.vcf \${REFPATH}
"""
}
