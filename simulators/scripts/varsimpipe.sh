#!/bin/bash
### NO PRESTEP
### POSTSTEP - truthpreprocessingVarSim.sh
#Varsim first creates a perturbed genome by adding known variants from a database

#SV -> (DGV)database of supporting variants
#SV -> insert_seq.txt
#SNPs
#The reads will be generated in the out directory. A script quickstart.sh is also provided which runs the above steps and generates reads for 1X coverage. Run the secondary analysis tools (alignment and variant calling) on those. The ground truth VCF file is also in the out directory called simu.truth.vcf
DATASIM=/data/backup/metacaller/download/varsim_run
REFERENCE=$1
OUTPATH=$2

mkdir -p $2

${DATASIM}/varsim.py --vc_in_vcf ${DATASIM}/All_20160601.vcf.gz --sv_insert_seq ${DATASIM}/insert_seq.txt --sv_dgv ${DATASIM}/GRCh37_hg19_supportingvariants_2013\-\07-\23.txt --reference ${REFERENCE} --id simu --read_length 100 --vc_num_snp 200000 --vc_num_ins 20000 --vc_num_del 20000 --vc_num_mnp 5000 --vc_num_complex 5000 --sv_num_ins 0 --sv_num_del 0 --sv_num_dup 0 --sv_num_inv 0 --sv_percent_novel 0.01 --vc_percent_novel 50 --mean_fragment_size 350 --sd_fragment_size 5 --vc_min_length_lim 0 --vc_max_length_lim 49 --sv_min_length_lim 50 --sv_max_length_lim 100000 --nlanes 1 --total_coverage 10 --simulator_executable ${DATASIM}/ART/art_bin_VanillaIceCream/art_illumina --out_dir ${OUTPATH}/out --log_dir ${OUTPATH}/log --work_dir ${OUTPATH}/work --simulator art
